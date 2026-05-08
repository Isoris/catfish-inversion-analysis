#!/usr/bin/env Rscript

# =============================================================================
# 04_detect_L1_localpca_zblocks.R
#
# Boundary-driven Level-1 envelope discovery on a single chromosome. Slides
# a small W x W upper-triangle cross-block at offset G along the diagonal
# of the sim_mat; each peak in -median(Z[cross]) is a candidate boundary.
# Validates each peak via the growing-W cross-block z test. Surviving peaks
# (STABLE_BLUE) partition the chromosome into L1 segments.
#
# Pipeline position:
#   03 precompute  -> 04 DETECT_L1  -> 05 plot_L1
#                                   -> 06 detect_L2  -> 07 plot_L2
#                                   -> 08b atlas_json
#
# ── Default NN scale: 80 ──────────────────────────────────────────────────
# Settled value (history line 19553). Chromosome-wide L1 detection needs to
# suppress short-range noise; nn80 is the right level. Override via --nn N.
#
# ── Inputs (canonical mode) ───────────────────────────────────────────────
#   --precomp_dir <dir>    auto-resolves <chr>.precomp.rds AND
#                          sim_mats/<chr>.sim_mat_nn{N}.rds from --chr + --nn
#   --chr         <label>  chromosome to operate on (e.g. C_gar_LG28)
#   --outdir      <dir>    output directory
#   [--nn         80]      NN smoothing level for the sim_mat
#
# ── Inputs (legacy / power-user mode) ─────────────────────────────────────
# Pass either or both to override auto-resolution:
#   --precomp     <path>   full path to <chr>.precomp.rds
#   --sim_mat     <path>   full path to <chr>.sim_mat_nn{N}.rds
#
# ── Boundary-scan parameters (canonical defaults from history 19553) ──────
#   --boundary_scan TRUE
#   --boundary_validator_mode grow      (alt: 'both' adds perpendicular ray)
#   --boundary_W              5
#   --boundary_offset         5
#   --boundary_min_dist       30        (sweep settled at 30, hist 19488-19510)
#
# ── Outputs (in --outdir) ─────────────────────────────────────────────────
#   <chr>.L1_envelopes.tsv     L1 segments (boundary-derived partition).
#                              One row per segment, status='ENVELOPE'.
#   <chr>.L1_boundaries.tsv    every detected peak with validation status:
#                                STABLE_BLUE = real long-range separation
#                                STABLE_RED  = rejected by grow validator
#                                DECAYS      = trigger fired but val score
#                                              collapses as W grows
#                                MARGINAL    = ambiguous
#                                EDGE        = too close to chr edge
#                              Filter --boundary_filter stable in step 05/07
#                              shows only STABLE_BLUE.
#   <chr>.L1_score_curve.tsv   per-window boundary score curve
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: STEP_D17_multipass_L1_only_v7.R)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ---- CLI --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

precomp_f      <- get_arg("--precomp")
sim_mat_f      <- get_arg("--sim_mat")
chr_label      <- get_arg("--chr", "chr")
outdir         <- get_arg("--outdir", ".")

# ---- Auto-resolution from --precomp_dir + --chr + --nn ----------------------
# Canonical run: user passes --precomp_dir <step03_outdir>/precomp + --chr +
# (optionally) --nn 80. Script resolves precomp + sim_mat paths automatically.
# Power users can still pass --precomp <path> and/or --sim_mat <path> directly
# to override. Default NN scale for L1 is 80 (settled value, history line 19553).
precomp_dir    <- get_arg("--precomp_dir")
nn_scale       <- as.integer(get_arg("--nn", "80"))

if (is.na(precomp_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, paste0(chr_label, ".precomp.rds"))
  if (file.exists(cand)) {
    precomp_f <- cand
    cat("[D17L1] resolved --precomp from --precomp_dir + --chr: ", cand, "\n")
  }
}
if (is.na(sim_mat_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, "sim_mats",
                    paste0(chr_label, ".sim_mat_nn", nn_scale, ".rds"))
  if (file.exists(cand)) {
    sim_mat_f <- cand
    cat("[D17L1] resolved --sim_mat from --precomp_dir + --chr + --nn=",
        nn_scale, ": ", cand, "\n", sep = "")
  }
}

# Diagonal boundary scan. Slides a small WxW upper-triangle cross-block
# at offset G along the diagonal. At each position i the boundary score
# is -median(Z) over the cross-block; peaks of this 1D signal mark
# positions where adjacent regions are unusually separated (i.e. the
# upper-triangle cross-block is unusually blue). The chromosome is then
# partitioned into L1 segments using the validated STABLE_BLUE peaks as
# cut points. This boundary-driven L1 is the only path now -- the
# legacy envelope-discovery scan (scan_pass_l1, NMS, trim_to_core,
# post-trim NMS) has been removed.
boundary_scan          <- !is.na(get_arg("--boundary_scan", NA_character_))
boundary_W             <- as.integer(get_arg("--boundary_W", "5"))
boundary_offset        <- as.integer(get_arg("--boundary_offset", "5"))
boundary_score_min     <- as.numeric(get_arg("--boundary_score_min", "2.0"))

# Which boundary statuses count as L1 cut points? Default just STABLE_BLUE.
l1_boundary_filter     <- tolower(get_arg("--l1_boundary_filter", "stable"))

# Minimum L1 segment size (in windows). Segments with nW <= this threshold
# are merged into the more-similar neighbor. Default 1 -- only catches
# exact 1-window adjacent-tie artifacts. Pass --l1_min_segment_nw 3 to
# also catch 2-3 window tight-cluster artifacts.
l1_min_segment_nw      <- as.integer(get_arg("--l1_min_segment_nw", "1"))

# v7-bp: --boundary_min_dist resolution priority:
#   1. --boundary_min_dist N         (absolute integer override, in windows)
#   2. --boundary_min_dist_kb K      (kb-based; resolves to windows via the
#                                     median window size in this chr)
#   3. --boundary_min_dist_pct P     (legacy percent-of-N; honored if given)
# Default: kb-based with K=150 kb. Reason: the right parameter for "min
# spacing between two REAL boundaries" is a biological quantity (smallest
# meaningful inversion size), measured in bp, not a fraction of the
# chromosome length. A 150 kb floor on inversion separation translates to
# 30 windows on a 5 kb-window chr -- matching the value you've been
# happy with for LG28.
# A floor (--boundary_min_dist_floor, default 5) keeps min_dist >= 5 even
# when windows are very large.
boundary_min_dist_arg     <- get_arg("--boundary_min_dist", NA_character_)
boundary_min_dist_kb_arg  <- get_arg("--boundary_min_dist_kb", NA_character_)
boundary_min_dist_pct_arg <- get_arg("--boundary_min_dist_pct", NA_character_)
boundary_min_dist_floor   <- as.integer(get_arg(
  "--boundary_min_dist_floor", "5"))
# boundary_min_dist itself is resolved later, after N and window sizes
# are known. Defaults: 150 kb (kb-based) when nothing else is given.
boundary_min_dist         <- NA_integer_

# Multi-W validation of each detected peak. From each peak position i,
# scan two 1D rays of diag-Z values perpendicular to the diagonal:
#   right ray: Z[i, i+d]   for d = 1, 2, ..., d_max
#   left  ray: Z[i-d, i]   for d = 1, 2, ..., d_max
# A REAL boundary stays mostly blue (Z < 0) on both rays — the matrix
# never "rebounds" to red as we look further from the diagonal. A FAKE
# boundary (a small blue cross inside an inversion) shows blue near the
# diagonal but jumps red abruptly as we exit the cross-shape into the
# surrounding red bulk.
boundary_validate      <- get_arg("--boundary_validate", "TRUE")
boundary_validate      <- !(toupper(boundary_validate) %in%
                            c("FALSE", "F", "0", "NO"))
# Max ray length (windows)
boundary_perp_d_max    <- as.integer(get_arg("--boundary_perp_d_max", "20"))
# A boundary is STABLE_BLUE if BOTH rays have:
#   (a) frac_blue >= boundary_perp_min_blue_frac      (mostly blue along ray)
#   (b) max_z     <= boundary_perp_max_red            (no big red excursion)
# It's DECAYS if either ray crosses red abruptly:
#   first_d_red  <= boundary_perp_first_d_red_max
boundary_perp_min_blue_frac    <- as.numeric(get_arg(
  "--boundary_perp_min_blue_frac", "0.70"))
boundary_perp_max_red          <- as.numeric(get_arg(
  "--boundary_perp_max_red", "0.50"))
boundary_perp_first_d_red_max  <- as.integer(get_arg(
  "--boundary_perp_first_d_red_max", "5"))
# z threshold defining "red" per pixel
boundary_perp_red_z            <- as.numeric(get_arg(
  "--boundary_perp_red_z", "0.50"))

# ---- Growing-W cross-block validator (v6) ----------------------------------
# Re-evaluates the upper-triangle cross-block at growing W and watches the
# median-Z of the cross-block.
#
#   At each peak i, for each W in --boundary_grow_W (default 5,20,50,100,150,200):
#     cross = sim_mat[(i-W+1):i, (i+G+1):(i+G+W)]
#     z = (cross - diag_mean) / diag_sd     (per-cell, by diagonal distance)
#     cross_z[W] = median(z)
#
# Decision (uses the LARGEST W in --boundary_grow_W for which the cross-block
# fits inside the chromosome):
#
#   REAL  : grow_max_z <=  boundary_grow_real_max     (default 0.10)
#                 -- never drifts positive at large W
#   FAKE  : grow_max_z >=  boundary_grow_fake_min     (default 0.20)
#                 -- drifts to slightly positive / positive at large W
#   MARGINAL : in between
#   EDGE  : largest fitting W is below boundary_grow_min_largest_W
#
# When --boundary_validator_mode=grow, validation_status is set from this
# decision (REAL -> STABLE_BLUE, FAKE -> DECAYS, MARGINAL -> MARGINAL,
# EDGE -> EDGE) so the overlay's existing filters keep working.
boundary_validator_mode        <- tolower(get_arg(
  "--boundary_validator_mode", "grow"))   # "grow" | "perp" | "both"

# v7-adaptive: the W series can be either an absolute comma-separated list
# (--boundary_grow_W) or, when not given, a percent-of-N list scaled to
# the chromosome length (--boundary_grow_W_pct). Defaults to a percent
# series spanning 0.1%, 0.5%, 1%, 2%, 4%, 7% of N -- ranges from W~5 on a
# 4000-window chromosome to W~300, which matches "5 to 300" intuition.
# Anything below W=5 is floored at 5 to keep the cross-block math
# meaningful. The actual absolute list used is logged on stderr.
boundary_grow_W_arg            <- get_arg("--boundary_grow_W", NA_character_)
boundary_grow_W_pct_str        <- get_arg("--boundary_grow_W_pct",
                                          "0.001,0.005,0.01,0.02,0.04,0.07")
boundary_grow_W_pct            <- as.numeric(strsplit(
  boundary_grow_W_pct_str, ",")[[1]])
boundary_grow_W_pct            <- sort(unique(boundary_grow_W_pct[
  is.finite(boundary_grow_W_pct) & boundary_grow_W_pct > 0]))
# Floor for the W series (whether from absolute or percent).
boundary_grow_W_floor          <- as.integer(get_arg(
  "--boundary_grow_W_floor", "5"))

# v7-adaptive: thresholds. Default is "adaptive" -- compute real_max and
# fake_min as quantiles of grow_max_z within this chromosome's detected
# peaks, with sanity floors so we don't call clearly-positive peaks REAL
# or clearly-negative peaks FAKE. Pass --boundary_grow_threshold_mode fixed
# to use the legacy hard cutoffs (--boundary_grow_real_max /
# --boundary_grow_fake_min).
boundary_grow_threshold_mode   <- tolower(get_arg(
  "--boundary_grow_threshold_mode", "adaptive"))   # "adaptive" | "fixed"
# Adaptive parameters (used only when threshold_mode = "adaptive")
boundary_grow_real_pct         <- as.numeric(get_arg(
  "--boundary_grow_real_pct", "0.50"))   # bottom 50% of grow_max_z = REAL
boundary_grow_fake_pct         <- as.numeric(get_arg(
  "--boundary_grow_fake_pct", "0.50"))   # above the median = FAKE
# Sanity floors: cap how lenient adaptive can get
boundary_grow_real_max_ceiling <- as.numeric(get_arg(
  "--boundary_grow_real_max_ceiling", "0.10"))
boundary_grow_fake_min_floor   <- as.numeric(get_arg(
  "--boundary_grow_fake_min_floor", "0.10"))
# Min sample size before we trust the adaptive quantile fit; fall back to
# fixed if fewer peaks have a finite grow_max_z than this.
boundary_grow_min_n_for_adaptive <- as.integer(get_arg(
  "--boundary_grow_min_n_for_adaptive", "10"))
# Fixed thresholds (used only when threshold_mode = "fixed")
boundary_grow_real_max         <- as.numeric(get_arg(
  "--boundary_grow_real_max", "0.10"))
boundary_grow_fake_min         <- as.numeric(get_arg(
  "--boundary_grow_fake_min", "0.20"))
boundary_grow_min_largest_W    <- as.integer(get_arg(
  "--boundary_grow_min_largest_W", "20"))

dry_run        <- !is.na(get_arg("--dry_run", NA_character_))

if (is.na(precomp_f) || !file.exists(precomp_f))
  stop("[D17L1] --precomp is required and must exist")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[D17L1] precomp:    ", precomp_f, "\n")
cat("[D17L1] sim_mat:    ", sim_mat_f %||% "(auto)", "\n")
cat("[D17L1] chr:        ", chr_label, "\n")
cat("[D17L1] outdir:     ", outdir, "\n")
cat("[D17L1] l1_boundary_filter=", l1_boundary_filter, "\n", sep = "")
cat("[D17L1] l1_min_segment_nw=", l1_min_segment_nw,
    " (segments smaller than this get merged into more-similar neighbor)\n",
    sep = "")
cat("[D17L1] boundary_scan=", boundary_scan, sep = "")
if (boundary_scan) {
  min_dist_desc <- if (!is.na(boundary_min_dist_arg)) {
    paste0(boundary_min_dist_arg, " windows (absolute override)")
  } else if (!is.na(boundary_min_dist_kb_arg)) {
    paste0(boundary_min_dist_kb_arg, " kb (window count resolved later)")
  } else if (!is.na(boundary_min_dist_pct_arg)) {
    paste0(boundary_min_dist_pct_arg, " of N (legacy percent; resolved later)")
  } else {
    paste0("150 kb (default; window count resolved later)")
  }
  cat("  W=", boundary_W,
      "  offset=", boundary_offset,
      "  score_min=", boundary_score_min,
      "  min_dist=", min_dist_desc,
      sep = "")
}
cat("\n")
if (boundary_scan && boundary_validate) {
  cat("[D17L1] boundary_validate=TRUE  (perp-ray mode)",
      "  d_max=", boundary_perp_d_max,
      "  min_blue_frac=", boundary_perp_min_blue_frac,
      "  max_red=", boundary_perp_max_red,
      "  first_d_red_max=", boundary_perp_first_d_red_max,
      "  red_z=", boundary_perp_red_z, "\n", sep = "")
}
if (boundary_scan) {
  cat("[D17L1] validator_mode=", boundary_validator_mode, sep = "")
  if (boundary_validator_mode %in% c("grow", "both")) {
    cat("\n[D17L1]   grow W series: ",
        if (!is.na(boundary_grow_W_arg)) paste0("absolute (", boundary_grow_W_arg, ")")
        else paste0("percent (", paste(boundary_grow_W_pct, collapse = ","),
                    " of N, floor=", boundary_grow_W_floor, ")"),
        sep = "")
    cat("\n[D17L1]   grow threshold mode: ", boundary_grow_threshold_mode, sep = "")
    if (boundary_grow_threshold_mode == "adaptive") {
      cat("\n[D17L1]     REAL side: real_max = min( quantile(grow_max_z, real_pct=",
          boundary_grow_real_pct, "), real_max_ceiling=",
          boundary_grow_real_max_ceiling, " )", sep = "")
      cat("\n[D17L1]     FAKE side: fake_min = max( quantile(grow_max_z, fake_pct=",
          boundary_grow_fake_pct, "), fake_min_floor=",
          boundary_grow_fake_min_floor, " )", sep = "")
    } else {
      cat("  (real_max=", boundary_grow_real_max,
          "  fake_min=", boundary_grow_fake_min, ")", sep = "")
    }
    cat("\n[D17L1]   min_largest_W=", boundary_grow_min_largest_W, sep = "")
  }
  cat("\n")
}

# ---- Load precomp -----------------------------------------------------------

cat("[D17L1] loading precomp\n")
pc_obj <- readRDS(precomp_f)
if (!is.null(pc_obj$dt)) {
  pc <- pc_obj
} else if (!is.null(pc_obj$pc)) {
  pc <- pc_obj$pc
} else {
  stop("[D17L1] cannot find pc$dt in precomp file")
}

dt_pc <- as.data.table(pc$dt)
n_windows_total <- nrow(dt_pc)

# ---- Normalize window coordinate columns -----------------------------------
if (all(c("start", "end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start
  window_end_bp   <- dt_pc$end
} else if (all(c("start_bp", "end_bp") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start_bp
  window_end_bp   <- dt_pc$end_bp
} else if (all(c("window_start", "window_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$window_start
  window_end_bp   <- dt_pc$window_end
} else if (all(c("bp_start", "bp_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$bp_start
  window_end_bp   <- dt_pc$bp_end
} else {
  stop("[D17L1] precomp dt missing usable coordinate columns. Found: ",
       paste(names(dt_pc), collapse = ", "))
}

cat("[D17L1] N windows: ", n_windows_total, "\n", sep = "")

# ---- Build boundary_grow_W (deferred until N is known) ---------------------
# Either honor the user's absolute --boundary_grow_W list, or scale the
# percent series --boundary_grow_W_pct against N. Always floored at
# --boundary_grow_W_floor (default 5) so the cross-block math stays valid.
if (!is.na(boundary_grow_W_arg)) {
  boundary_grow_W <- as.integer(strsplit(boundary_grow_W_arg, ",")[[1]])
  boundary_grow_W <- boundary_grow_W[is.finite(boundary_grow_W)]
} else {
  boundary_grow_W <- as.integer(round(boundary_grow_W_pct * n_windows_total))
}
boundary_grow_W <- pmax(boundary_grow_W, boundary_grow_W_floor)
boundary_grow_W <- sort(unique(boundary_grow_W))
if (boundary_validator_mode %in% c("grow", "both")) {
  cat("[D17L1] grow W resolved (N=", n_windows_total, "): ",
      paste(boundary_grow_W, collapse = ","), "\n", sep = "")
}

# ---- Resolve boundary_min_dist (deferred until N and window sizes known)
# Priority:
#   1. Absolute --boundary_min_dist N (windows, override)
#   2. --boundary_min_dist_kb K  (kb-based; uses median window size)
#   3. --boundary_min_dist_pct P (legacy percent-of-N)
#   4. Default: 150 kb (kb-based)
# All paths floored at --boundary_min_dist_floor (default 5).
window_size_bp_med <- as.numeric(median(window_end_bp - window_start_bp + 1L,
                                        na.rm = TRUE))
if (!is.finite(window_size_bp_med) || window_size_bp_med <= 0) {
  window_size_bp_med <- NA_real_
}

if (!is.na(boundary_min_dist_arg)) {
  boundary_min_dist <- as.integer(boundary_min_dist_arg)
  resolution_note <- sprintf("absolute override (%d windows)", boundary_min_dist)
} else if (!is.na(boundary_min_dist_kb_arg)) {
  kb <- as.numeric(boundary_min_dist_kb_arg)
  if (!is.finite(window_size_bp_med))
    stop("[D17L1] cannot resolve boundary_min_dist_kb: window size unknown")
  boundary_min_dist <- as.integer(round(kb * 1000 / window_size_bp_med))
  resolution_note <- sprintf("%g kb / %.0f bp/window = %d windows",
                             kb, window_size_bp_med, boundary_min_dist)
} else if (!is.na(boundary_min_dist_pct_arg)) {
  pct <- as.numeric(boundary_min_dist_pct_arg)
  boundary_min_dist <- as.integer(round(pct * n_windows_total))
  resolution_note <- sprintf("legacy %g of N=%d -> %d windows",
                             pct, n_windows_total, boundary_min_dist)
} else {
  # Default: 150 kb -> windows
  default_kb <- 150
  if (!is.finite(window_size_bp_med)) {
    boundary_min_dist <- 30L
    resolution_note <- "default 30 windows (window size unknown)"
  } else {
    boundary_min_dist <- as.integer(round(default_kb * 1000 / window_size_bp_med))
    resolution_note <- sprintf("default %d kb / %.0f bp/window = %d windows",
                               default_kb, window_size_bp_med, boundary_min_dist)
  }
}
boundary_min_dist <- max(boundary_min_dist, boundary_min_dist_floor)
if (boundary_scan) {
  cat("[D17L1] boundary_min_dist resolved: ", resolution_note,
      "  (final=", boundary_min_dist, " windows, floor=",
      boundary_min_dist_floor, ")\n", sep = "")
}

# ---- Load sim_mat -----------------------------------------------------------

if (is.na(sim_mat_f)) {
  candidates_path <- c(
    file.path(dirname(precomp_f), "sim_mat_nn80.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn160.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn40.rds")
  )
  hit <- candidates_path[file.exists(candidates_path)]
  if (length(hit) == 0)
    stop("[D17L1] --sim_mat not given and could not auto-find sim_mat_nn*.rds")
  sim_mat_f <- hit[1]
  cat("[D17L1] auto-found sim_mat: ", sim_mat_f, "\n")
}

cat("[D17L1] loading sim_mat\n")
sm_obj <- readRDS(sim_mat_f)
if (is.matrix(sm_obj)) {
  sim_mat <- sm_obj
} else if (is.list(sm_obj) && !is.null(sm_obj$sim_mat) && is.matrix(sm_obj$sim_mat)) {
  sim_mat <- sm_obj$sim_mat
} else if (is.list(sm_obj) && length(sm_obj) == 1 && is.matrix(sm_obj[[1]])) {
  sim_mat <- sm_obj[[1]]
} else {
  stop("[D17L1] sim_mat object structure not recognized; class=", class(sm_obj)[1])
}
storage.mode(sim_mat) <- "double"

if (!isTRUE(nrow(sim_mat) == n_windows_total &&
            ncol(sim_mat) == n_windows_total)) {
  stop("[D17L1] sim_mat dim ", nrow(sim_mat), "x", ncol(sim_mat),
       " does not match precomp N windows ", n_windows_total)
}


# ---- Optional: diagonal boundary scan --------------------------------------
# Slides a WxW upper-triangle cross-block at offset G along the diagonal.
# At each window position i:
#   cross = Z[(i-W+1):i, (i+G+1):(i+G+W)]
#   boundary_score[i] = -median(cross)   (bluer cross = higher score)
# Then 1D peak detection: a peak is a local maximum within +/- min_dist
# windows that exceeds boundary_score_min (in z-units).

if (boundary_scan) {
  cat("\n[D17L1] === DIAGONAL BOUNDARY SCAN ===\n", sep = "")
  cat("[D17L1]   W=", boundary_W, "  offset=", boundary_offset,
      "  score_min=", boundary_score_min,
      "  min_dist=", boundary_min_dist, "\n", sep = "")

  N <- n_windows_total
  W <- boundary_W
  G <- boundary_offset

  # Step 1: build per-diagonal mean and sd for Z normalization.
  cat("[D17L1]   precomputing diagonal mean/sd for Z\n")
  diag_mean <- numeric(N)
  diag_sd   <- numeric(N)
  for (d in 0:(N - 1L)) {
    if (d == 0L) {
      vals <- diag(sim_mat)
    } else {
      ii <- seq.int(1L, N - d)
      vals <- sim_mat[cbind(ii, ii + d)]
    }
    vals <- vals[is.finite(vals)]
    if (length(vals) >= 5L) {
      diag_mean[d + 1L] <- mean(vals)
      sg <- sd(vals)
      diag_sd[d + 1L]   <- if (is.finite(sg) && sg > 1e-9) sg else NA_real_
    } else {
      diag_mean[d + 1L] <- NA_real_
      diag_sd[d + 1L]   <- NA_real_
    }
  }

  # Step 2: at each diagonal position i, score the cross-block.
  # Valid i: cross-block must fit. Need (i - W + 1) >= 1 and (i + G + W) <= N.
  i_lo <- W
  i_hi <- N - G - W
  if (i_hi < i_lo) {
    cat("[D17L1]   chromosome too short for this W/offset — skipping scan\n")
  } else {
    centers <- seq.int(i_lo, i_hi)
    boundary_score <- rep(NA_real_, N)
    cat("[D17L1]   scoring ", length(centers),
        " positions along the diagonal\n", sep = "")
    for (ctr in centers) {
      ii <- seq.int(ctr - W + 1L, ctr)
      jj <- seq.int(ctr + G + 1L, ctr + G + W)
      block <- sim_mat[ii, jj]
      d_mat <- outer(ii, jj, FUN = function(a, b) b - a)
      idx <- d_mat + 1L
      mu <- diag_mean[idx]
      sg <- diag_sd[idx]
      z  <- (block - mu) / sg
      z[!is.finite(z)] <- NA_real_
      vals <- as.numeric(z)
      vals <- vals[is.finite(vals)]
      if (length(vals) >= 5L) {
        boundary_score[ctr] <- -median(vals)
      }
    }
    score_finite <- boundary_score[is.finite(boundary_score)]
    if (length(score_finite) > 0L) {
      cat(sprintf(
        "[D17L1]   boundary_score range: [%.2f, %.2f]   median=%.2f   p95=%.2f\n",
        min(score_finite), max(score_finite),
        median(score_finite),
        as.numeric(quantile(score_finite, 0.95))))
    }

    # Step 3: 1D peak detection.
    # A peak: (a) score >= boundary_score_min,
    #         (b) score is the max within +/- min_dist windows.
    peaks <- integer(0)
    for (ctr in centers) {
      s <- boundary_score[ctr]
      if (!is.finite(s) || s < boundary_score_min) next
      lo <- max(1L, ctr - boundary_min_dist)
      hi <- min(N, ctr + boundary_min_dist)
      win <- boundary_score[lo:hi]
      win <- win[is.finite(win)]
      if (length(win) == 0L) next
      if (s >= max(win)) peaks <- c(peaks, ctr)
    }
    cat("[D17L1]   peaks found: ", length(peaks), "\n", sep = "")

    # ---- Multi-W validation of each peak ----------------------------------
    # For each peak position i, scan two 1D rays of diag-Z values
    # perpendicular to the diagonal:
    #   right ray: Z[i, i+d]   for d = 1, 2, ..., d_max
    #   left  ray: Z[i-d, i]   for d = 1, 2, ..., d_max
    # By matrix symmetry these are the same data row/col-swapped, but we
    # keep both as separate rays in case one runs off the chromosome edge.
    #
    # A REAL boundary stays mostly blue (Z < 0) on both rays — the matrix
    # never rebounds to red as we look further from the diagonal. A FAKE
    # boundary (a small blue cross inside an inversion) shows blue near
    # the diagonal but jumps red abruptly as we exit the cross-shape.
    perp_z <- function(r, c) {
      # Single-pixel diagonal-distance Z at sim_mat[r, c]
      if (r < 1L || r > n_windows_total) return(NA_real_)
      if (c < 1L || c > n_windows_total) return(NA_real_)
      d <- abs(c - r)
      mu <- diag_mean[d + 1L]
      sg <- diag_sd[d + 1L]
      sm <- sim_mat[r, c]
      if (!is.finite(sm) || !is.finite(mu) || !is.finite(sg) || sg < 1e-9)
        return(NA_real_)
      (sm - mu) / sg
    }

    # Build output table.
    if (length(peaks) > 0L) {
      bdt <- data.table(
        chr            = chr_label,
        boundary_idx   = sprintf("%s.L1_b%04d", chr_label,
                                 seq_along(peaks)),
        boundary_w     = peaks,
        boundary_bp    = window_start_bp[peaks],
        boundary_score = boundary_score[peaks],
        boundary_W     = W,
        boundary_offset = G
      )
      bdt <- bdt[order(boundary_w)]

      if (boundary_validate) {
        cat("[D17L1]   validating peaks via perpendicular rays (d_max=",
            boundary_perp_d_max, ", red_z=", boundary_perp_red_z, ")\n",
            sep = "")
        # Per-peak per-ray stats
        n_peaks <- nrow(bdt)
        right_frac_blue <- numeric(n_peaks)
        left_frac_blue  <- numeric(n_peaks)
        right_max_z     <- numeric(n_peaks)
        left_max_z      <- numeric(n_peaks)
        right_first_red <- integer(n_peaks)
        left_first_red  <- integer(n_peaks)
        right_n_finite  <- integer(n_peaks)
        left_n_finite   <- integer(n_peaks)
        right_curves    <- character(n_peaks)
        left_curves     <- character(n_peaks)

        for (k in seq_len(n_peaks)) {
          i <- bdt$boundary_w[k]
          d_seq <- seq_len(boundary_perp_d_max)
          # Right ray
          r_vals <- vapply(d_seq, function(d) perp_z(i, i + d), numeric(1))
          # Left ray
          l_vals <- vapply(d_seq, function(d) perp_z(i - d, i), numeric(1))

          # Right stats
          ok_r <- is.finite(r_vals)
          right_n_finite[k] <- sum(ok_r)
          if (any(ok_r)) {
            right_frac_blue[k] <- mean(r_vals[ok_r] < 0)
            right_max_z[k]     <- max(r_vals[ok_r])
            red_idx <- which(r_vals > boundary_perp_red_z & ok_r)
            right_first_red[k] <- if (length(red_idx) > 0L)
                                    as.integer(red_idx[1]) else NA_integer_
          } else {
            right_frac_blue[k] <- NA_real_
            right_max_z[k]     <- NA_real_
            right_first_red[k] <- NA_integer_
          }
          # Left stats
          ok_l <- is.finite(l_vals)
          left_n_finite[k] <- sum(ok_l)
          if (any(ok_l)) {
            left_frac_blue[k] <- mean(l_vals[ok_l] < 0)
            left_max_z[k]     <- max(l_vals[ok_l])
            red_idx <- which(l_vals > boundary_perp_red_z & ok_l)
            left_first_red[k] <- if (length(red_idx) > 0L)
                                   as.integer(red_idx[1]) else NA_integer_
          } else {
            left_frac_blue[k] <- NA_real_
            left_max_z[k]     <- NA_real_
            left_first_red[k] <- NA_integer_
          }
          # Compact curve strings (every other pixel for readability)
          show_idx <- if (boundary_perp_d_max <= 10L) seq_len(boundary_perp_d_max)
                      else c(1L, seq.int(2L, boundary_perp_d_max - 1L, by = 2L),
                             boundary_perp_d_max)
          right_curves[k] <- paste(
            ifelse(is.finite(r_vals[show_idx]),
                   sprintf("%+.1f", r_vals[show_idx]), "NA"),
            collapse = ",")
          left_curves[k] <- paste(
            ifelse(is.finite(l_vals[show_idx]),
                   sprintf("%+.1f", l_vals[show_idx]), "NA"),
            collapse = ",")
        }

        bdt[, right_frac_blue := right_frac_blue]
        bdt[, left_frac_blue  := left_frac_blue]
        bdt[, right_max_z     := right_max_z]
        bdt[, left_max_z      := left_max_z]
        bdt[, right_first_red := right_first_red]
        bdt[, left_first_red  := left_first_red]
        bdt[, right_n_finite  := right_n_finite]
        bdt[, left_n_finite   := left_n_finite]
        bdt[, right_ray_curve := right_curves]
        bdt[, left_ray_curve  := left_curves]

        # Classify (perpendicular-ray validator — kept as a diagnostic)
        classify_perp <- function(rfb, lfb, rmz, lmz, rfr, lfr) {
          # EDGE if either ray has too few finite samples
          if (!is.finite(rfb) || !is.finite(lfb)) return("EDGE")
          # DECAYS: either ray crosses red close to the diagonal
          first_red_too_close <-
            (is.finite(rfr) && rfr <= boundary_perp_first_d_red_max) ||
            (is.finite(lfr) && lfr <= boundary_perp_first_d_red_max)
          if (first_red_too_close) return("DECAYS")
          # STABLE_BLUE: both rays stay mostly blue + no big red excursion
          stable <- (rfb >= boundary_perp_min_blue_frac) &&
                    (lfb >= boundary_perp_min_blue_frac) &&
                    (is.finite(rmz) && rmz <= boundary_perp_max_red) &&
                    (is.finite(lmz) && lmz <= boundary_perp_max_red)
          if (stable) return("STABLE_BLUE")
          "MARGINAL"
        }
        bdt[, perp_status := mapply(classify_perp,
                                    right_frac_blue, left_frac_blue,
                                    right_max_z,     left_max_z,
                                    right_first_red, left_first_red)]
      }  # end if (boundary_validate)

      # ---- Growing-W cross-block validator (v6) ---------------------------
      # For each peak i and each W in boundary_grow_W, compute the median Z
      # of the upper-triangle cross-block:
      #   block = sim_mat[(i-W+1):i, (i+G+1):(i+G+W)]
      # using the existing diag_mean / diag_sd. Watch how median(Z) evolves
      # with W. A REAL boundary's cross-block stays blue (median Z << 0) or
      # drifts toward 0 but does NOT become positive, even at large W.
      # A FAKE boundary (interior cross of a single inversion) drifts to
      # slightly positive / positive at large W.
      if (boundary_validator_mode %in% c("grow", "both") &&
          length(boundary_grow_W) > 0L) {
        cat("[D17L1]   growing-W validator: W=",
            paste(boundary_grow_W, collapse = ","),
            "  threshold_mode=", boundary_grow_threshold_mode, "\n", sep = "")

        n_peaks <- nrow(bdt)
        # Allocate one column per W
        grow_mat <- matrix(NA_real_, nrow = n_peaks,
                           ncol = length(boundary_grow_W))
        colnames(grow_mat) <- paste0("cross_z_W", boundary_grow_W)

        for (k in seq_len(n_peaks)) {
          i <- bdt$boundary_w[k]
          for (wj in seq_along(boundary_grow_W)) {
            Wg <- boundary_grow_W[wj]
            # Cross-block must fit inside the chromosome.
            i_lo <- i - Wg + 1L
            i_hi <- i
            j_lo <- i + G + 1L
            j_hi <- i + G + Wg
            if (i_lo < 1L || j_hi > n_windows_total) next
            ii <- seq.int(i_lo, i_hi)
            jj <- seq.int(j_lo, j_hi)
            block <- sim_mat[ii, jj]
            d_mat <- outer(ii, jj, FUN = function(a, b) b - a)
            idx <- d_mat + 1L
            mu <- diag_mean[idx]
            sg <- diag_sd[idx]
            z  <- (block - mu) / sg
            z[!is.finite(z)] <- NA_real_
            vals <- as.numeric(z)
            vals <- vals[is.finite(vals)]
            if (length(vals) >= 5L) {
              grow_mat[k, wj] <- median(vals)
            }
          }
        }
        # Attach as columns
        for (wj in seq_along(boundary_grow_W)) {
          col <- colnames(grow_mat)[wj]
          bdt[, (col) := grow_mat[, wj]]
        }

        # Per-peak grow_max_z = max(median Z) over fitting W values.
        # Largest fitting W for that peak determines EDGE classification.
        grow_max_z         <- numeric(n_peaks)
        grow_largest_W     <- integer(n_peaks)
        grow_z_at_largest  <- numeric(n_peaks)
        for (k in seq_len(n_peaks)) {
          row_z <- grow_mat[k, ]
          ok <- which(is.finite(row_z))
          if (length(ok) == 0L) {
            grow_max_z[k]        <- NA_real_
            grow_largest_W[k]    <- NA_integer_
            grow_z_at_largest[k] <- NA_real_
          } else {
            grow_max_z[k]        <- max(row_z[ok])
            grow_largest_W[k]    <- as.integer(boundary_grow_W[max(ok)])
            grow_z_at_largest[k] <- row_z[max(ok)]
          }
        }
        bdt[, grow_max_z        := grow_max_z]
        bdt[, grow_largest_W    := grow_largest_W]
        bdt[, grow_z_at_largest := grow_z_at_largest]

        # ---- Resolve thresholds (adaptive vs fixed) -----------------------
        # Adaptive: real_max = quantile(grow_max_z, real_pct), capped at
        #           real_max_ceiling so we never call a clearly-positive
        #           peak REAL.
        #           fake_min = quantile(grow_max_z, fake_pct), floored at
        #           fake_min_floor so we never call a clearly-negative
        #           peak FAKE.
        # If the threshold mode is "fixed" or there aren't enough finite
        # samples, use the fixed parameter values directly.
        finite_gmz <- grow_max_z[is.finite(grow_max_z)]
        use_adaptive <- (boundary_grow_threshold_mode == "adaptive") &&
                        (length(finite_gmz) >=
                           boundary_grow_min_n_for_adaptive)
        if (use_adaptive) {
          q_real <- as.numeric(quantile(finite_gmz,
                                         boundary_grow_real_pct,
                                         na.rm = TRUE))
          q_fake <- as.numeric(quantile(finite_gmz,
                                         boundary_grow_fake_pct,
                                         na.rm = TRUE))
          resolved_real_max <- min(q_real, boundary_grow_real_max_ceiling)
          resolved_fake_min <- max(q_fake, boundary_grow_fake_min_floor)
          # Guarantee the gap real_max < fake_min stays open. If quantiles
          # invert (rare on tiny chromosomes), nudge fake_min up by 0.05.
          gap_collision <- resolved_fake_min <= resolved_real_max
          if (gap_collision) {
            resolved_fake_min <- resolved_real_max + 0.05
          }
          # Two separate lines, one per side of the band. Shows which of
          # (quantile, guard) won the min/max -- useful when calibrating
          # whether the data or the guard is driving the threshold.
          cat(sprintf(
            "[D17L1]   adaptive thresholds (n=%d peaks):\n",
            length(finite_gmz)))
          cat(sprintf(
            "[D17L1]     REAL ceiling: real_max = min(q%.2f(gmz)=%+.3f, ceiling=%+.3f) = %+.3f  [%s won]  -- peaks with gmz <= real_max are REAL (cross stays bluish at all W)\n",
            boundary_grow_real_pct, q_real,
            boundary_grow_real_max_ceiling, resolved_real_max,
            if (q_real <= boundary_grow_real_max_ceiling) "quantile" else "ceiling"))
          cat(sprintf(
            "[D17L1]     FAKE floor:   fake_min = max(q%.2f(gmz)=%+.3f, floor=%+.3f)   = %+.3f  [%s won]  -- peaks with gmz >= fake_min are FAKE (cross drifts red at large W)\n",
            boundary_grow_fake_pct, q_fake,
            boundary_grow_fake_min_floor,
            if (gap_collision) resolved_real_max + 0.05 else max(q_fake, boundary_grow_fake_min_floor),
            if (q_fake >= boundary_grow_fake_min_floor) "quantile" else "floor"))
          if (gap_collision) {
            cat("[D17L1]     gap-collision adjust: fake_min nudged to real_max+0.05 to keep band open\n")
          }
        } else {
          if (boundary_grow_threshold_mode == "adaptive" &&
              length(finite_gmz) < boundary_grow_min_n_for_adaptive) {
            cat("[D17L1]   adaptive thresholds requested but only ",
                length(finite_gmz), " peaks have finite grow_max_z (need ",
                boundary_grow_min_n_for_adaptive,
                "). Falling back to fixed.\n", sep = "")
          }
          resolved_real_max <- boundary_grow_real_max
          resolved_fake_min <- boundary_grow_fake_min
          cat(sprintf(
            "[D17L1]   fixed thresholds: real_max=%+.3f  fake_min=%+.3f  -- peaks below real_max are REAL, above fake_min are FAKE, in-between are MARGINAL\n",
            resolved_real_max, resolved_fake_min))
        }
        bdt[, resolved_real_max := resolved_real_max]
        bdt[, resolved_fake_min := resolved_fake_min]

        # Classify using resolved thresholds
        classify_grow <- function(gmax, glargest) {
          if (!is.finite(gmax) || !is.finite(glargest)) return("EDGE")
          if (glargest < boundary_grow_min_largest_W)  return("EDGE")
          if (gmax <= resolved_real_max)               return("REAL")
          if (gmax >= resolved_fake_min)               return("FAKE")
          "MARGINAL"
        }
        bdt[, grow_status := mapply(classify_grow,
                                    grow_max_z, grow_largest_W)]
      }

      # ---- Final validation_status decision -------------------------------
      # The overlay reads validation_status with the legacy vocabulary
      # STABLE_BLUE / DECAYS / MARGINAL / EDGE. Map the chosen validator's
      # output into that vocabulary so the overlay's filters keep working.
      if (boundary_validator_mode == "perp") {
        if ("perp_status" %in% names(bdt)) {
          bdt[, validation_status := perp_status]
        } else {
          bdt[, validation_status := "STABLE_BLUE"]
        }
      } else if (boundary_validator_mode == "grow") {
        if ("grow_status" %in% names(bdt)) {
          bdt[, validation_status := fifelse(
            grow_status == "REAL",     "STABLE_BLUE",
            fifelse(grow_status == "FAKE",     "DECAYS",
            fifelse(grow_status == "MARGINAL", "MARGINAL",
                                                "EDGE")))]
        } else {
          bdt[, validation_status := "EDGE"]
        }
      } else if (boundary_validator_mode == "both") {
        # Conservative AND: a peak is STABLE_BLUE only if BOTH validators
        # call it real. If either calls it fake (DECAYS / FAKE), it's DECAYS.
        # Otherwise MARGINAL or EDGE based on what the perp said.
        if ("perp_status" %in% names(bdt) && "grow_status" %in% names(bdt)) {
          combine <- function(ps, gs) {
            if (is.na(ps) || is.na(gs)) return("EDGE")
            if (ps == "EDGE" || gs == "EDGE") return("EDGE")
            if (gs == "FAKE" || ps == "DECAYS") return("DECAYS")
            if (ps == "STABLE_BLUE" && gs == "REAL") return("STABLE_BLUE")
            "MARGINAL"
          }
          bdt[, validation_status := mapply(combine, perp_status, grow_status)]
        } else if ("grow_status" %in% names(bdt)) {
          bdt[, validation_status := fifelse(
            grow_status == "REAL",     "STABLE_BLUE",
            fifelse(grow_status == "FAKE",     "DECAYS",
            fifelse(grow_status == "MARGINAL", "MARGINAL",
                                                "EDGE")))]
        } else if ("perp_status" %in% names(bdt)) {
          bdt[, validation_status := perp_status]
        } else {
          bdt[, validation_status := "STABLE_BLUE"]
        }
      } else {
        # Unknown mode — be safe and pass everything through as STABLE_BLUE
        bdt[, validation_status := "STABLE_BLUE"]
      }

      if (boundary_validate || boundary_validator_mode %in% c("grow", "both")) {

        n_stable <- sum(bdt$validation_status == "STABLE_BLUE")
        n_decay  <- sum(bdt$validation_status == "DECAYS")
        n_marg   <- sum(bdt$validation_status == "MARGINAL")
        n_edge   <- sum(bdt$validation_status == "EDGE")
        cat("[D17L1]   validation summary: STABLE_BLUE=", n_stable,
            "  DECAYS=", n_decay,
            "  MARGINAL=", n_marg,
            "  EDGE=", n_edge,
            "  (mode=", boundary_validator_mode, ")\n", sep = "")
        if ("right_frac_blue" %in% names(bdt)) {
          cat(sprintf(
            "[D17L1]   right ray: median frac_blue=%.2f  max_z=%.2f\n",
            median(bdt$right_frac_blue, na.rm = TRUE),
            median(bdt$right_max_z, na.rm = TRUE)))
          cat(sprintf(
            "[D17L1]   left  ray: median frac_blue=%.2f  max_z=%.2f\n",
            median(bdt$left_frac_blue, na.rm = TRUE),
            median(bdt$left_max_z, na.rm = TRUE)))
        }
        if ("grow_max_z" %in% names(bdt)) {
          cat(sprintf(
            "[D17L1]   grow_max_z: min=%+.2f  median=%+.2f  max=%+.2f  (n_finite=%d/%d)\n",
            min(bdt$grow_max_z, na.rm = TRUE),
            median(bdt$grow_max_z, na.rm = TRUE),
            max(bdt$grow_max_z, na.rm = TRUE),
            sum(is.finite(bdt$grow_max_z)), nrow(bdt)))
          for (Wg in boundary_grow_W) {
            col <- paste0("cross_z_W", Wg)
            if (col %in% names(bdt)) {
              v <- bdt[[col]]
              if (any(is.finite(v))) {
                cat(sprintf(
                  "[D17L1]     W=%-3d  median(cross_z)=%+.2f  n_finite=%d\n",
                  Wg, median(v, na.rm = TRUE), sum(is.finite(v))))
              }
            }
          }
        }
      }

      cat("[D17L1]   top 10 boundaries by trigger score:\n")
      top <- bdt[order(-boundary_score)][1:min(10L, .N)]
      has_status <- "validation_status" %in% names(top)
      has_perp   <- "right_frac_blue"   %in% names(top)
      has_grow   <- "grow_max_z"        %in% names(top)
      for (k in seq_len(nrow(top))) {
        line <- sprintf("[D17L1]     %s  w=%d  %.2f Mb  trig=%.2f",
                        top$boundary_idx[k], top$boundary_w[k],
                        top$boundary_bp[k] / 1e6, top$boundary_score[k])
        if (has_status) {
          line <- sprintf("%s  status=%-11s",
                          line, top$validation_status[k])
        }
        if (has_grow) {
          line <- sprintf("%s  grow_max_z=%+.2f",
                          line, top$grow_max_z[k])
        }
        if (has_perp) {
          line <- sprintf("%s  rfb=%.2f lfb=%.2f rmz=%+.2f lmz=%+.2f",
                          line,
                          top$right_frac_blue[k], top$left_frac_blue[k],
                          top$right_max_z[k],     top$left_max_z[k])
        }
        cat(line, "\n", sep = "")
      }
    } else {
      bdt <- data.table(
        chr = character(0), boundary_idx = character(0),
        boundary_w = integer(0), boundary_bp = integer(0),
        boundary_score = numeric(0),
        boundary_W = integer(0), boundary_offset = integer(0)
      )
    }

    if (!dry_run) {
      out_b <- file.path(outdir, paste0(chr_label, ".L1_boundaries.tsv"))
      fwrite(bdt, out_b, sep = "\t")
      cat("[D17L1]   boundaries written: ", out_b, "\n", sep = "")

      # Also write the full per-window score curve, for diagnostics.
      curve_dt <- data.table(
        chr = chr_label,
        window_idx = seq_len(N),
        bp = window_start_bp,
        boundary_score = boundary_score
      )
      out_c <- file.path(outdir,
                         paste0(chr_label, ".L1_score_curve.tsv"))
      fwrite(curve_dt, out_c, sep = "\t")
      cat("[D17L1]   score curve written: ", out_c, "\n", sep = "")
    }
  }
}

# ---- L1-from-boundaries: derive segments and write envelope catalogue ------
# Partitions the chromosome into L1 segments using the validated
# STABLE_BLUE boundaries as cut points:
#   segment 1   : [1,         b_1]
#   segment k   : [b_{k-1}+1, b_k]    for k = 2 .. K
#   segment K+1 : [b_K+1,     N]
# where b_1 < ... < b_K are the boundary_w positions (in window units)
# of all surviving STABLE_BLUE peaks. If --l1_boundary_filter=non_decay
# or =all, MARGINAL/EDGE/DECAYS are also used as cuts.
if (!dry_run && exists("bdt")) {
  cut_keep_set <- switch(l1_boundary_filter,
    "stable"    = c("STABLE_BLUE"),
    "non_decay" = c("STABLE_BLUE", "MARGINAL", "EDGE"),
    "all"       = c("STABLE_BLUE", "DECAYS", "MARGINAL", "EDGE"),
    c("STABLE_BLUE")
  )
  if ("validation_status" %in% names(bdt)) {
    cut_dt <- bdt[validation_status %in% cut_keep_set]
  } else {
    cut_dt <- bdt   # no validation columns -> use all detected peaks
  }
  cuts <- sort(unique(as.integer(cut_dt$boundary_w)))
  cuts <- cuts[is.finite(cuts) & cuts >= 1L & cuts <= n_windows_total]

  if (length(cuts) == 0L) {
    cat("[D17L1] no surviving boundaries -- writing empty envelope catalogue.\n",
        sep = "")
    out_main <- file.path(outdir, paste0(chr_label, ".L1_envelopes.tsv"))
    fwrite(data.table(
      chr = character(0), candidate_id = character(0),
      start_w = integer(0), end_w = integer(0),
      start_bp = integer(0), end_bp = integer(0),
      n_windows = integer(0), scale_W = integer(0),
      mean_sim = numeric(0), density_p70 = numeric(0),
      density_adaptive = numeric(0), threshold = numeric(0),
      trim_start_w_original = integer(0), trim_end_w_original = integer(0),
      trim_left_dropped = integer(0), trim_right_dropped = integer(0),
      trim_total_dropped = integer(0), status = character(0)
    ), out_main, sep = "\t")
    cat("[D17L1] empty envelope catalogue written: ", out_main, "\n", sep = "")
  } else {
    # Build segment endpoints. Segment k: rows [start_k .. end_k] in window
    # space. Segment 1 starts at 1. Subsequent segments start one past the
    # previous cut. The last segment ends at N.
    seg_starts <- c(1L, cuts + 1L)
    seg_ends   <- c(cuts, n_windows_total)
    # Drop empty/degenerate segments (start > end).
    keep_seg <- seg_ends >= seg_starts
    seg_starts <- seg_starts[keep_seg]
    seg_ends   <- seg_ends[keep_seg]
    n_seg <- length(seg_starts)

    # ---- Merge tiny L1 segments into more-similar neighbor ---------------
    # When two close peaks pass the local-max filter, the segment between
    # them can be very small (1-3 windows). Such segments are usually
    # detector artifacts: the cross-block scan fired twice on adjacent
    # diagonal positions for the same architectural transition, and the
    # space "between" them is a sliver, not a real architectural region.
    #
    # For each tiny segment (nW <= l1_min_segment_nw), compute mean
    # similarity to its left and right neighbors using the full block
    # mean(sim_mat[tiny_windows, neighbor_windows]). Merge into whichever
    # side is more similar. This is more principled than dropping a peak
    # by score because it uses the heatmap data to decide where the
    # boundary really is.
    n_merged <- 0L
    repeat {
      seg_nw <- seg_ends - seg_starts + 1L
      tiny_idx <- which(seg_nw <= l1_min_segment_nw)
      if (length(tiny_idx) == 0L) break
      k <- tiny_idx[1]                       # process leftmost first
      tiny_s <- seg_starts[k]
      tiny_e <- seg_ends[k]
      tiny_n <- tiny_e - tiny_s + 1L
      has_left  <- k > 1L
      has_right <- k < length(seg_starts)
      sim_left  <- if (has_left) {
        Ls <- seg_starts[k - 1L]; Le <- seg_ends[k - 1L]
        v <- as.numeric(sim_mat[tiny_s:tiny_e, Ls:Le])
        v <- v[is.finite(v)]
        if (length(v) > 0L) mean(v) else NA_real_
      } else NA_real_
      sim_right <- if (has_right) {
        Rs <- seg_starts[k + 1L]; Re <- seg_ends[k + 1L]
        v <- as.numeric(sim_mat[tiny_s:tiny_e, Rs:Re])
        v <- v[is.finite(v)]
        if (length(v) > 0L) mean(v) else NA_real_
      } else NA_real_

      # Decide which side to merge into. Forced when only one side exists.
      merge_left <- if (!has_right) TRUE
                    else if (!has_left) FALSE
                    else if (!is.finite(sim_left) && !is.finite(sim_right)) TRUE
                    else if (!is.finite(sim_right)) TRUE
                    else if (!is.finite(sim_left))  FALSE
                    else (sim_left >= sim_right)

      if (merge_left) {
        seg_ends[k - 1L] <- seg_ends[k]      # extend left to absorb
        seg_starts <- seg_starts[-k]
        seg_ends   <- seg_ends[-k]
        cat(sprintf("[D17L1] merged %d-window L1 at w=%d-%d into LEFT neighbor (sim_left=%.3f vs sim_right=%.3f)\n",
                    tiny_n, tiny_s, tiny_e, sim_left, sim_right))
      } else {
        seg_starts[k + 1L] <- seg_starts[k]  # extend right to absorb
        seg_starts <- seg_starts[-k]
        seg_ends   <- seg_ends[-k]
        cat(sprintf("[D17L1] merged %d-window L1 at w=%d-%d into RIGHT neighbor (sim_left=%.3f vs sim_right=%.3f)\n",
                    tiny_n, tiny_s, tiny_e, sim_left, sim_right))
      }
      n_merged <- n_merged + 1L
    }
    if (n_merged > 0L) {
      cat("[D17L1] total tiny L1 segments merged (nW <= ",
          l1_min_segment_nw, "): ", n_merged, "\n", sep = "")
    }
    n_seg <- length(seg_starts)
    # ----------------------------------------------------------------------

    seg_dt <- data.table(
      chr          = chr_label,
      candidate_id = sprintf("%s.L1_%04d", chr_label, seq_len(n_seg)),
      start_w      = seg_starts,
      end_w        = seg_ends,
      start_bp     = window_start_bp[seg_starts],
      end_bp       = window_end_bp[seg_ends],
      n_windows    = seg_ends - seg_starts + 1L,
      scale_W      = seg_ends - seg_starts + 1L,
      status       = "ENVELOPE"
    )
    # Compute mean_sim per segment from the actual sim_mat block.
    seg_dt[, mean_sim := vapply(seq_len(.N), function(k) {
      s <- start_w[k]; e <- end_w[k]
      if (e <= s) return(NA_real_)
      block <- sim_mat[s:e, s:e]
      mean(block, na.rm = TRUE)
    }, numeric(1))]
    seg_dt[, density_p70 := vapply(seq_len(.N), function(k) {
      s <- start_w[k]; e <- end_w[k]
      if (e <= s) return(NA_real_)
      block <- sim_mat[s:e, s:e]
      mean(block >= 0.70, na.rm = TRUE)
    }, numeric(1))]
    # Schema parity columns (kept for downstream compatibility; the trim/
    # sees the same column set).
    seg_dt[, density_adaptive      := NA_real_]
    seg_dt[, threshold             := NA_real_]
    seg_dt[, trim_start_w_original := NA_integer_]
    seg_dt[, trim_end_w_original   := NA_integer_]
    seg_dt[, trim_left_dropped     := NA_integer_]
    seg_dt[, trim_right_dropped    := NA_integer_]
    seg_dt[, trim_total_dropped    := NA_integer_]

    # Reorder to a stable column order.
    setcolorder(seg_dt, c(
      "chr", "candidate_id", "start_w", "end_w",
      "start_bp", "end_bp", "n_windows", "scale_W",
      "mean_sim", "density_p70", "density_adaptive",
      "threshold",
      "trim_start_w_original", "trim_end_w_original",
      "trim_left_dropped", "trim_right_dropped", "trim_total_dropped",
      "status"))

    out_main <- file.path(outdir, paste0(chr_label, ".L1_envelopes.tsv"))
    fwrite(seg_dt, out_main, sep = "\t")
    cat("[D17L1] L1 segments derived from boundaries (filter='",
        l1_boundary_filter, "', n_cuts=", length(cuts),
        ", n_segments=", n_seg, ") written: ", out_main, "\n", sep = "")
    cat("[D17L1] L1 segments:\n", sep = "")
    for (k in seq_len(n_seg)) {
      cat(sprintf("[D17L1]   %s  %.2f-%.2f Mb  nW=%d  mean_sim=%.3f\n",
                  seg_dt$candidate_id[k],
                  seg_dt$start_bp[k] / 1e6, seg_dt$end_bp[k] / 1e6,
                  seg_dt$n_windows[k], seg_dt$mean_sim[k]))
    }
  }
}

cat("[D17L1] done\n")
