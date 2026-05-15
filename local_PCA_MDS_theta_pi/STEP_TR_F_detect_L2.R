#!/usr/bin/env Rscript

# =============================================================================
# 06_detect_L2_localpca_zblocks.R
#
# Level-2 sub-boundary discovery. Reads the L1 envelope catalogue from
# step 04 and runs the same boundary-scan logic INDEPENDENTLY INSIDE EACH
# L1 SEGMENT, producing a finer L2 partition of every L1 region.
#
# Pipeline position:
#   04 detect_L1  -> 06 DETECT_L2  -> 07 plot_L2  -> 08b atlas_json
#
# ── Default NN scale: 40 ──────────────────────────────────────────────────
# Finer than L1's nn80 because we're zooming inside an L1 segment - smoother
# sim_mats wash out sub-structure here. Override via --nn N.
#
# ── Two key v8 features ───────────────────────────────────────────────────
#
# WARD-ADAPTIVE THRESHOLDING (--ward_adaptive yes, default)
#   For each L1 segment, computes Ward.D2-on-windows clustering and
#   intensity statistics (mean_sim, blue_red_ratio). Derives two priors,
#   chromosome-rank-normalized:
#       internal_homogeneity     high => looks like one inversion inside
#       architectural_complexity high => real architectural transitions
#       strictness = internal_homogeneity - architectural_complexity (in [-1,+1])
#   Each segment then gets its OWN (real_max, fake_min) by scaling pooled
#   baselines: STRICT segments tighten (most peaks reject as FAKE); LOOSE
#   segments widen (more peaks pass as REAL). Self-calibrating per chrom,
#   no hardcoded class table or magic mean_sim cutoff.
#   Pass --ward_adaptive no for pooled-only thresholds (legacy behaviour).
#
# QUADRANT VALIDATOR (--quadrant_validator yes, default in canonical run)
#   2D test: a candidate L2 boundary inside an L1 segment is consistent
#   with a true sub-block iff the four sim_mat quadrants it divides
#   behave as a true block. Demote-on-fail rescues weak peaks via grow_z.
#
# ── Inputs (canonical mode) ───────────────────────────────────────────────
#   --precomp_dir <dir>    auto-resolves precomp + sim_mat_nn{N}
#   --L1_dir      <dir>    auto-resolves <chr>.L1_envelopes.tsv
#   --chr         <label>
#   --outdir      <dir>
#   [--nn         40]
#
# ── Inputs (legacy) ───────────────────────────────────────────────────────
#   --precomp     <path>   --sim_mat   <path>
#   --l1_catalogue <path>  full path to <chr>.L1_envelopes.tsv
#
# ── Boundary-scan + validator parameters (canonical from history 19585) ───
#   --boundary_scan TRUE
#   --boundary_validator_mode grow
#   --quadrant_validator yes
#   --weak_demote_score      0
#   --quad_rescue_max_grow_z 1.5
#   --quad_demote_on_fail    yes
#   --quad_demote_drift_floor -1.0
#
# ── Outputs (in --outdir) ─────────────────────────────────────────────────
#   <chr>.L2_envelopes.tsv             L2 sub-segments
#   <chr>.L2_boundaries.tsv            all detected peaks + validation
#   <chr>.L2_segment_stats.tsv         per-L1-segment Ward stats + strictness
#   <chr>.L2_quadrant_validator.tsv    per-peak quadrant test (when enabled)
#   <chr>.L2_quadrant_audit.tsv        extra audit info
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: STEP_D17_multipass_L2_v8.R)
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
l1_cat_f       <- get_arg("--l1_catalogue")
chr_label      <- get_arg("--chr", "chr")
outdir         <- get_arg("--outdir", ".")

# ---- Auto-resolution: --precomp_dir + --chr + --nn (default 40 for L2) ------
# Plus --L1_dir for the L1 catalogue. L2 uses nn40 by default (history
# 19556 onwards) -- finer resolution than L1's nn80 since we're zooming
# inside an L1 segment.
precomp_dir    <- get_arg("--precomp_dir")
nn_scale       <- as.integer(get_arg("--nn", "40"))
L1_dir         <- get_arg("--L1_dir")

if (is.na(precomp_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, paste0(chr_label, ".precomp.rds"))
  if (file.exists(cand)) {
    precomp_f <- cand
    cat("[D17L2] resolved --precomp: ", cand, "\n")
  }
}
if (is.na(sim_mat_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, "sim_mats",
                    paste0(chr_label, ".sim_mat_nn", nn_scale, ".rds"))
  if (file.exists(cand)) {
    sim_mat_f <- cand
    cat("[D17L2] resolved --sim_mat (nn=", nn_scale, "): ", cand, "\n", sep = "")
  }
}
if (is.na(l1_cat_f) && !is.na(L1_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(L1_dir, paste0(chr_label, ".L1_envelopes.tsv"))
  if (file.exists(cand)) {
    l1_cat_f <- cand
    cat("[D17L2] resolved --l1_catalogue: ", cand, "\n")
  }
}

if (is.na(precomp_f) || !file.exists(precomp_f))
  stop("[D17L2] --precomp is required and must exist (or --precomp_dir + --chr)")
if (is.na(sim_mat_f) || !file.exists(sim_mat_f))
  stop("[D17L2] --sim_mat is required and must exist (or --precomp_dir + --chr + --nn)")
if (is.na(l1_cat_f) || !file.exists(l1_cat_f))
  stop("[D17L2] --l1_catalogue is required and must exist (or --L1_dir + --chr)")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Boundary scan parameters
boundary_W             <- as.integer(get_arg("--boundary_W", "5"))
boundary_offset        <- as.integer(get_arg("--boundary_offset", "5"))
# Shift the catalogue boundary_w from kernel anchor (last row of upper
# trigger arm) to gap-center (= visible center of the blue patch). The
# kernel computes scores using rows [ctr-W+1..ctr] x cols [ctr+G+1..ctr+G+W]
# and records boundary_w_local = ctr (the kernel anchor). Empirically the
# *visible* architectural edge sits at ctr + ceiling((G+1)/2). Setting
# this flag to "yes" (default) makes the OUTPUT catalogue (boundary_w,
# boundary_bp in the boundaries TSV; cut positions in the envelopes TSV)
# refer to the visible-center coordinate. Internal validators (perp ray,
# grow-W, quadrant audit) continue to use boundary_w_local, which stays
# at the kernel anchor, so internal logic is unaffected.
# Set to "no" to keep the legacy kernel-anchor convention in the output
# (off by ceiling((G+1)/2) = 3 windows from the visible center with the
# default G=5).
anchor_at_gap_center <- tolower(get_arg("--anchor_at_gap_center", "yes"))
anchor_at_gap_center <- anchor_at_gap_center %in% c("yes", "y", "true", "t", "1", "on")
boundary_w_shift <- if (anchor_at_gap_center)
  as.integer(ceiling((boundary_offset + 1L) / 2L)) else 0L
boundary_score_min     <- as.numeric(get_arg("--boundary_score_min", "0.5"))
boundary_min_dist      <- as.integer(get_arg("--boundary_min_dist", "6"))
# After scan + validation, drop any peak within --boundary_dedup_dist of a
# higher-scoring already-kept peak. This handles the case where the scan
# kernel fires twice on the same architectural edge at slightly different
# offsets. Set to 0 to disable.
boundary_dedup_dist    <- as.integer(get_arg("--boundary_dedup_dist", "15"))

# Perp-ray validator parameters
boundary_validate      <- get_arg("--boundary_validate", "TRUE")
boundary_validate      <- !(toupper(boundary_validate) %in%
                            c("FALSE", "F", "0", "NO"))
boundary_perp_d_max    <- as.integer(get_arg("--boundary_perp_d_max", "20"))
boundary_perp_min_blue_frac    <- as.numeric(get_arg(
  "--boundary_perp_min_blue_frac", "0.70"))
boundary_perp_max_red          <- as.numeric(get_arg(
  "--boundary_perp_max_red", "0.50"))
boundary_perp_first_d_red_max  <- as.integer(get_arg(
  "--boundary_perp_first_d_red_max", "5"))
boundary_perp_red_z            <- as.numeric(get_arg(
  "--boundary_perp_red_z", "0.50"))

# Growing-W validator parameters
boundary_validator_mode        <- tolower(get_arg(
  "--boundary_validator_mode", "grow"))
boundary_grow_W_arg            <- get_arg("--boundary_grow_W", NA_character_)
boundary_grow_W_pct_str        <- get_arg("--boundary_grow_W_pct",
                                          "0.005,0.02,0.05,0.10,0.20,0.40")
# NOTE: at L2 the percent-of-N series is WIDER than at L1 because segments
# are smaller. 40% of a 500-window segment = 200 windows -- still meaningful.
boundary_grow_W_pct            <- as.numeric(strsplit(
  boundary_grow_W_pct_str, ",")[[1]])
boundary_grow_W_pct            <- sort(unique(boundary_grow_W_pct[
  is.finite(boundary_grow_W_pct) & boundary_grow_W_pct > 0]))
boundary_grow_W_floor          <- as.integer(get_arg(
  "--boundary_grow_W_floor", "5"))

# Threshold parameters (pooled across all L2 peaks)
boundary_grow_threshold_mode   <- tolower(get_arg(
  "--boundary_grow_threshold_mode", "adaptive"))
boundary_grow_real_pct         <- as.numeric(get_arg(
  "--boundary_grow_real_pct", "0.30"))
boundary_grow_fake_pct         <- as.numeric(get_arg(
  "--boundary_grow_fake_pct", "0.50"))
boundary_grow_real_max_ceiling <- as.numeric(get_arg(
  "--boundary_grow_real_max_ceiling", "0.20"))
boundary_grow_fake_min_floor   <- as.numeric(get_arg(
  "--boundary_grow_fake_min_floor", "0.20"))
boundary_grow_min_n_for_adaptive <- as.integer(get_arg(
  "--boundary_grow_min_n_for_adaptive", "10"))
boundary_grow_real_max         <- as.numeric(get_arg(
  "--boundary_grow_real_max", "0.10"))
boundary_grow_fake_min         <- as.numeric(get_arg(
  "--boundary_grow_fake_min", "0.20"))
boundary_grow_min_largest_W    <- as.integer(get_arg(
  "--boundary_grow_min_largest_W", "10"))

# Ward.D2-adaptive per-segment threshold tuning (default: yes)
ward_adaptive <- tolower(get_arg("--ward_adaptive", "yes"))
ward_adaptive <- ward_adaptive %in% c("yes", "y", "true", "t", "1", "on")
ward_min_segment_n <- as.integer(get_arg("--ward_min_segment_n", "30"))
ward_max_n         <- as.integer(get_arg("--ward_max_n",         "300"))
ward_strictness_max <- as.numeric(get_arg("--ward_strictness_max", "2.0"))
ward_ceil_max       <- as.numeric(get_arg("--ward_ceil_max",       "0.5"))

# Quadrant validator (opt-in). Auditing pass that runs AFTER grow-W
# classification. For each peak, builds the local diagonal square at
# multiple W values, splits each into 2x2 sub-quadrants, and tests
# whether the quadrants are mutually homogeneous (low max-min range)
# AND not red-drifted (mean not strongly positive).
#   --quadrant_validator           default "no". Set "yes" to enable.
#   --quad_w_sweep                 comma-separated W values; default "10,20,30"
#   --quad_homog_thresh            range across q1..q4 below which W passes
#                                  default 1.5
#   --quad_drift_thresh            mean of q1..q4 above which W fails
#                                  default 0.3 (toward red)
#   --quad_min_pass                minimum number of W values that must
#                                  pass to declare HOMOG_PASS; default 2
#   --quad_min_room                minimum windows of room around the peak
#                                  required before running the audit. Peaks
#                                  too close to a segment edge get verdict
#                                  NO_ROOM and are not rescued or flagged.
#                                  default 20.
# Behaviour on classified peaks:
#   - Grow-W said FAKE/DECAYS, quadrant says HOMOG_PASS -> RESCUE to
#     STABLE_BLUE (peak gains validation).
#   - Grow-W said STABLE_BLUE, quadrant says HOMOG_FAIL -> kept as
#     STABLE_BLUE but quadrant_failed=TRUE flag set for downstream audit.
#   - Verdict NO_ROOM or NO_W_FIT (no usable square): grow-W decision
#     is left untouched (no rescue, no flag).
#   - Other combinations leave validation_status unchanged.
quadrant_validator <- tolower(get_arg("--quadrant_validator", "no"))
quadrant_validator <- quadrant_validator %in% c("yes", "y", "true", "t", "1", "on")
quad_w_sweep_str <- get_arg("--quad_w_sweep", "10,20,30")
quad_w_sweep <- as.integer(strsplit(quad_w_sweep_str, ",")[[1]])
quad_w_sweep <- sort(unique(quad_w_sweep[is.finite(quad_w_sweep) & quad_w_sweep > 0]))
quad_homog_thresh <- as.numeric(get_arg("--quad_homog_thresh", "1.5"))
# Drift default tightened from 0.3 to 0.0: requires the 4 quadrant means
# to AVERAGE neutral or negative (not red-leaning). 0.3 was empirically
# too permissive on flat/red regions where local pockets of zero z let
# false positives slip through. Set to 0.3 to restore old behaviour.
quad_drift_thresh <- as.numeric(get_arg("--quad_drift_thresh", "0.0"))
quad_min_pass     <- as.integer(get_arg("--quad_min_pass", "2"))
quad_min_room     <- as.integer(get_arg("--quad_min_room", "20"))
# Minimum window-distance between this peak and its nearest neighbor peak
# within the same L1 segment. If a peak has any other peak (regardless of
# validation status) closer than this, the audit declines to opine and
# returns verdict NO_NEIGHBOR_ROOM. Catches dense peak clusters (e.g.
# 5+ peaks within 30 windows of each other in a uniformly-coloured patch)
# where the surrounding context is dominated by other peaks rather than
# the architectural neighborhood, making the audit unreliable. Default 20.
quad_min_neighbor_dist <- as.integer(get_arg("--quad_min_neighbor_dist", "20"))
# Cap on grow_max_z above which quadrant rescue is refused. Even if the
# local quadrants pass HOMOG_PASS, if grow-W found grow_max_z > this cap
# the broader cross is too red-leaning to trust. Catches the case of a
# tiny pocket of negative z right at the diagonal inside an otherwise red
# region. Default +0.2. Set very high (e.g. 5) to disable.
quad_rescue_max_grow_z <- as.numeric(get_arg("--quad_rescue_max_grow_z", "0.2"))
# When TRUE (default), STABLE_BLUE peaks with quadrant verdict HOMOG_FAIL
# are DEMOTED to MARGINAL, removing them from the L2 envelope cuts. When
# FALSE (legacy), they keep STABLE_BLUE status with quadrant_failed=TRUE
# flag for inspection but still cut envelopes.
quad_demote_on_fail <- tolower(get_arg("--quad_demote_on_fail", "yes"))
quad_demote_on_fail <- quad_demote_on_fail %in% c("yes", "y", "true", "t", "1", "on")
# Floor on quad_min_drift below which STABLE_BLUE peaks are NEVER demoted
# even with HOMOG_FAIL. A peak with mean quadrant z below this is in
# unambiguously deep-blue territory; the homogeneity failure is between
# two flavours of "blue" (e.g. -1 vs -3) and is not biologically
# meaningful. Default -1.0. Set to -Inf to disable (always demote on
# HOMOG_FAIL when quad_demote_on_fail=yes).
quad_demote_drift_floor <- as.numeric(get_arg("--quad_demote_drift_floor", "-1.0"))

# Final pre-output strictness gate: STABLE_BLUE peaks where BOTH the
# kernel score AND the grow_max_z signal are weak get demoted to MARGINAL.
# Catches peaks the kernel found but where neither piece of evidence is
# convincing on its own, AND where the quadrant validator couldn't help
# (e.g. NO_NEIGHBOR_ROOM peaks in clusters).
#   --weak_demote_score    threshold below which boundary_score is "weak"
#                          default 1.5
#   --weak_demote_gmz      threshold above which grow_max_z is "weak"
#                          (closer to zero / less blue). Default -0.5.
# A peak is demoted only when boundary_score < weak_demote_score AND
# grow_max_z > weak_demote_gmz (both must be weak). Set --weak_demote_score
# to 0 to disable.
weak_demote_score <- as.numeric(get_arg("--weak_demote_score", "1.5"))
weak_demote_gmz   <- as.numeric(get_arg("--weak_demote_gmz",   "-0.5"))

# Filter: which L1 segments to scan? Default all. Pass a regex against
# parent_l1_id to subset.
l1_filter_regex <- get_arg("--l1_filter_regex", NA_character_)

dry_run <- !is.na(get_arg("--dry_run", NA_character_))

# ---- Banner -----------------------------------------------------------------

cat("[D17L2] precomp:      ", precomp_f,    "\n", sep = "")
cat("[D17L2] sim_mat:      ", sim_mat_f,    "\n", sep = "")
cat("[D17L2] l1_catalogue: ", l1_cat_f,     "\n", sep = "")
cat("[D17L2] chr:          ", chr_label,    "\n", sep = "")
cat("[D17L2] outdir:       ", outdir,       "\n", sep = "")
cat("[D17L2] boundary scan:  W=", boundary_W,
    "  offset=", boundary_offset,
    "  score_min=", boundary_score_min,
    "  min_dist=", boundary_min_dist, "\n", sep = "")
cat("[D17L2] validator_mode=", boundary_validator_mode,
    "  threshold_mode=", boundary_grow_threshold_mode,
    "\n", sep = "")
cat("[D17L2] grow W series percent-of-segment-N: ",
    paste(boundary_grow_W_pct, collapse = ","),
    " (floor=", boundary_grow_W_floor, ")\n", sep = "")
cat("[D17L2] ward_adaptive=", ward_adaptive,
    "  min_segment_n=", ward_min_segment_n,
    "  max_n=", ward_max_n,
    "  strictness_max=", ward_strictness_max,
    "  ceil_max=", ward_ceil_max, "\n", sep = "")
cat("[D17L2] quadrant_validator=", quadrant_validator,
    "  W_sweep=", paste(quad_w_sweep, collapse = ","),
    "  homog_thresh=", quad_homog_thresh,
    "  drift_thresh=", quad_drift_thresh,
    "  min_pass=", quad_min_pass,
    "  min_room=", quad_min_room,
    "  min_neighbor_dist=", quad_min_neighbor_dist, "\n", sep = "")
cat("[D17L2] quadrant strictness: rescue_max_grow_z=", quad_rescue_max_grow_z,
    "  demote_on_fail=", quad_demote_on_fail,
    "  demote_drift_floor=", quad_demote_drift_floor, "\n", sep = "")
cat("[D17L2] weak-peak demote gate: score<", weak_demote_score,
    " AND gmz>", weak_demote_gmz, " -> demote STABLE -> MARGINAL\n", sep = "")
cat("[D17L2] anchor_at_gap_center=", anchor_at_gap_center,
    "  catalogue boundary_w shift = +", boundary_w_shift,
    " (kernel anchor -> ",
    if (anchor_at_gap_center) "visible center" else "kernel anchor (legacy)",
    ")\n", sep = "")

# ---- Load precomp -----------------------------------------------------------

cat("[D17L2] loading precomp\n")
pc <- readRDS(precomp_f)
dt_pc <- if (!is.null(pc$dt)) {
  as.data.table(pc$dt)
} else if (!is.null(pc$pc) && !is.null(pc$pc$dt)) {
  as.data.table(pc$pc$dt)
} else {
  stop("[D17L2] cannot find dt inside precomp")
}

if (!is.null(dt_pc$start)) {
  window_start_bp <- dt_pc$start
  window_end_bp   <- dt_pc$end
} else if (!is.null(dt_pc$start_bp)) {
  window_start_bp <- dt_pc$start_bp
  window_end_bp   <- dt_pc$end_bp
} else if (!is.null(dt_pc$window_start)) {
  window_start_bp <- dt_pc$window_start
  window_end_bp   <- dt_pc$window_end
} else if (!is.null(dt_pc$bp_start)) {
  window_start_bp <- dt_pc$bp_start
  window_end_bp   <- dt_pc$bp_end
} else {
  stop("[D17L2] precomp dt is missing window start/end columns")
}

n_windows_total <- length(window_start_bp)
cat("[D17L2] N windows: ", n_windows_total, "\n", sep = "")

# ---- Load sim_mat -----------------------------------------------------------

cat("[D17L2] loading sim_mat\n")
sm_obj <- readRDS(sim_mat_f)
sim_mat <- if (is.matrix(sm_obj)) {
  sm_obj
} else if (!is.null(sm_obj$sim_mat)) {
  sm_obj$sim_mat
} else if (is.list(sm_obj) && length(sm_obj) >= 1L) {
  sm_obj[[1]]
} else {
  stop("[D17L2] sim_mat structure not recognized")
}
if (nrow(sim_mat) != n_windows_total)
  stop(sprintf("[D17L2] sim_mat has %d rows but precomp has %d windows",
               nrow(sim_mat), n_windows_total))

# ---- Load L1 catalogue ------------------------------------------------------

cat("[D17L2] loading L1 catalogue\n")
l1_dt <- fread(l1_cat_f)
required_cols <- c("candidate_id", "start_w", "end_w")
missing <- setdiff(required_cols, names(l1_dt))
if (length(missing) > 0L)
  stop("[D17L2] L1 catalogue is missing columns: ",
       paste(missing, collapse = ", "))
if (!is.na(l1_filter_regex)) {
  l1_dt <- l1_dt[grepl(l1_filter_regex, candidate_id)]
}
n_l1_segments <- nrow(l1_dt)
cat("[D17L2] L1 segments to scan: ", n_l1_segments, "\n", sep = "")
if (n_l1_segments == 0L) {
  cat("[D17L2] no L1 segments to scan -- nothing to do\n")
  quit(status = 0L)
}

# ---- Helper: per-distance Z stats on a sub-matrix ---------------------------
build_diag_stats <- function(sim_seg) {
  N <- nrow(sim_seg)
  diag_mean <- numeric(N)
  diag_sd   <- numeric(N)
  for (d in 0:(N - 1L)) {
    if (d == 0L) {
      vals <- diag(sim_seg)
    } else {
      ii <- seq.int(1L, N - d)
      vals <- sim_seg[cbind(ii, ii + d)]
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
  list(diag_mean = diag_mean, diag_sd = diag_sd)
}

# ---- Helper: run a boundary scan on one segment -----------------------------
# Returns a data.table with both LOCAL (boundary_w_local) and ABSOLUTE
# (boundary_w) window indices. No final classification yet -- pool first.
scan_segment <- function(sim_seg, parent_id, abs_offset) {
  N <- nrow(sim_seg)
  W <- boundary_W
  G <- boundary_offset

  ds <- build_diag_stats(sim_seg)
  diag_mean <- ds$diag_mean
  diag_sd   <- ds$diag_sd

  # Per-segment grow W series
  if (!is.na(boundary_grow_W_arg)) {
    grow_W <- as.integer(strsplit(boundary_grow_W_arg, ",")[[1]])
    grow_W <- grow_W[is.finite(grow_W)]
  } else {
    grow_W <- as.integer(round(boundary_grow_W_pct * N))
  }
  grow_W <- pmax(grow_W, boundary_grow_W_floor)
  grow_W <- sort(unique(grow_W))

  # Trigger scoring
  i_lo <- W
  i_hi <- N - G - W
  if (i_hi < i_lo) return(NULL)
  centers <- seq.int(i_lo, i_hi)
  boundary_score <- rep(NA_real_, N)
  for (ctr in centers) {
    ii <- seq.int(ctr - W + 1L, ctr)
    jj <- seq.int(ctr + G + 1L, ctr + G + W)
    block <- sim_seg[ii, jj]
    d_mat <- outer(ii, jj, FUN = function(a, b) b - a)
    idx <- d_mat + 1L
    mu <- diag_mean[idx]; sg <- diag_sd[idx]
    z <- (block - mu) / sg
    z[!is.finite(z)] <- NA_real_
    vals <- as.numeric(z); vals <- vals[is.finite(vals)]
    if (length(vals) >= 5L) boundary_score[ctr] <- -median(vals)
  }

  # 1D peak detection
  peaks <- integer(0)
  for (ctr in centers) {
    s <- boundary_score[ctr]
    if (!is.finite(s) || s < boundary_score_min) next
    lo <- max(1L, ctr - boundary_min_dist)
    hi <- min(N, ctr + boundary_min_dist)
    win <- boundary_score[lo:hi]; win <- win[is.finite(win)]
    if (length(win) == 0L) next
    if (s >= max(win)) peaks <- c(peaks, ctr)
  }
  if (length(peaks) == 0L) return(NULL)

  perp_z <- function(r, c) {
    if (r < 1L || r > N) return(NA_real_)
    if (c < 1L || c > N) return(NA_real_)
    d <- abs(c - r)
    mu <- diag_mean[d + 1L]; sg <- diag_sd[d + 1L]
    sm <- sim_seg[r, c]
    if (!is.finite(sm) || !is.finite(mu) || !is.finite(sg) || sg < 1e-9)
      return(NA_real_)
    (sm - mu) / sg
  }

  bdt <- data.table(
    chr              = chr_label,
    parent_l1_id     = parent_id,
    boundary_w_local = peaks,                              # kernel anchor (local)
    # boundary_w is the CATALOGUE coordinate. It receives the
    # boundary_w_shift (= +ceiling((G+1)/2) when anchor_at_gap_center=yes,
    # 0 otherwise) so the TSV refers to the visible center of the blue
    # patch instead of the kernel anchor. Internal validators use
    # boundary_w_local (above), which is unaffected.
    boundary_w       = peaks + abs_offset + boundary_w_shift,
    boundary_bp      = window_start_bp[
      pmin(length(window_start_bp),
           pmax(1L, peaks + abs_offset + boundary_w_shift))],
    boundary_score   = boundary_score[peaks],
    boundary_W       = W,
    boundary_offset  = G
  )
  bdt <- bdt[order(boundary_w_local)]

  # Perp-ray validation
  if (boundary_validate) {
    n_peaks <- nrow(bdt)
    rfb <- numeric(n_peaks); lfb <- numeric(n_peaks)
    rmz <- numeric(n_peaks); lmz <- numeric(n_peaks)
    rfr <- integer(n_peaks); lfr <- integer(n_peaks)
    rnf <- integer(n_peaks); lnf <- integer(n_peaks)
    d_seq <- seq_len(boundary_perp_d_max)
    for (k in seq_len(n_peaks)) {
      i <- bdt$boundary_w_local[k]
      r_vals <- vapply(d_seq, function(d) perp_z(i, i + d), numeric(1))
      l_vals <- vapply(d_seq, function(d) perp_z(i - d, i), numeric(1))
      ok_r <- is.finite(r_vals); rnf[k] <- sum(ok_r)
      if (any(ok_r)) {
        rfb[k] <- mean(r_vals[ok_r] < 0); rmz[k] <- max(r_vals[ok_r])
        ri <- which(r_vals > boundary_perp_red_z & ok_r)
        rfr[k] <- if (length(ri) > 0L) as.integer(ri[1]) else NA_integer_
      } else { rfb[k] <- NA_real_; rmz[k] <- NA_real_; rfr[k] <- NA_integer_ }
      ok_l <- is.finite(l_vals); lnf[k] <- sum(ok_l)
      if (any(ok_l)) {
        lfb[k] <- mean(l_vals[ok_l] < 0); lmz[k] <- max(l_vals[ok_l])
        li <- which(l_vals > boundary_perp_red_z & ok_l)
        lfr[k] <- if (length(li) > 0L) as.integer(li[1]) else NA_integer_
      } else { lfb[k] <- NA_real_; lmz[k] <- NA_real_; lfr[k] <- NA_integer_ }
    }
    bdt[, right_frac_blue := rfb]
    bdt[, left_frac_blue  := lfb]
    bdt[, right_max_z     := rmz]
    bdt[, left_max_z      := lmz]
    bdt[, right_first_red := rfr]
    bdt[, left_first_red  := lfr]
    bdt[, right_n_finite  := rnf]
    bdt[, left_n_finite   := lnf]
    perp_status <- rep("STABLE_BLUE", n_peaks)
    for (k in seq_len(n_peaks)) {
      r_ok <- is.finite(rfb[k]) && rfb[k] >= boundary_perp_min_blue_frac &&
              is.finite(rmz[k]) && rmz[k] <= boundary_perp_max_red
      l_ok <- is.finite(lfb[k]) && lfb[k] >= boundary_perp_min_blue_frac &&
              is.finite(lmz[k]) && lmz[k] <= boundary_perp_max_red
      if (!(r_ok && l_ok)) {
        early_red <-
          (is.finite(rfr[k]) && rfr[k] <= boundary_perp_first_d_red_max) ||
          (is.finite(lfr[k]) && lfr[k] <= boundary_perp_first_d_red_max)
        perp_status[k] <- if (early_red) "DECAYS" else "MARGINAL"
      }
    }
    bdt[, perp_status := perp_status]
  }

  # Growing-W cross-block z
  if (boundary_validator_mode %in% c("grow", "both")) {
    n_peaks <- nrow(bdt)
    grow_mat <- matrix(NA_real_, nrow = n_peaks, ncol = length(grow_W))
    colnames(grow_mat) <- sprintf("cross_z_W%d", grow_W)
    for (k in seq_len(n_peaks)) {
      i <- bdt$boundary_w_local[k]
      for (wj in seq_along(grow_W)) {
        Wgrow <- grow_W[wj]
        ii <- (i - Wgrow + 1L):i
        jj <- (i + G + 1L):(i + G + Wgrow)
        if (ii[1] < 1L) next
        if (jj[length(jj)] > N) next
        block <- sim_seg[ii, jj]
        d_mat <- outer(ii, jj, FUN = function(a, b) b - a)
        idx <- d_mat + 1L
        mu <- diag_mean[idx]; sg <- diag_sd[idx]
        z  <- (block - mu) / sg
        z[!is.finite(z)] <- NA_real_
        vals <- as.numeric(z); vals <- vals[is.finite(vals)]
        if (length(vals) >= 5L) grow_mat[k, wj] <- median(vals)
      }
    }
    for (wj in seq_along(grow_W)) {
      col <- colnames(grow_mat)[wj]
      bdt[, (col) := grow_mat[, wj]]
    }
    grow_max_z <- numeric(n_peaks)
    grow_largest_W <- integer(n_peaks)
    grow_z_at_largest <- numeric(n_peaks)
    for (k in seq_len(n_peaks)) {
      row_z <- grow_mat[k, ]
      ok <- which(is.finite(row_z))
      if (length(ok) == 0L) {
        grow_max_z[k]        <- NA_real_
        grow_largest_W[k]    <- NA_integer_
        grow_z_at_largest[k] <- NA_real_
      } else {
        grow_max_z[k]        <- max(row_z[ok])
        grow_largest_W[k]    <- as.integer(grow_W[max(ok)])
        grow_z_at_largest[k] <- row_z[max(ok)]
      }
    }
    bdt[, grow_max_z        := grow_max_z]
    bdt[, grow_largest_W    := grow_largest_W]
    bdt[, grow_z_at_largest := grow_z_at_largest]
  }

  bdt
}

# ---- Run scan on every L1 segment -------------------------------------------

# Helper: compute Ward.D2 + intensity stats for a single L1 segment.
# Returns a one-row data.table or NULL if the segment is too small.
# Stats:
#   N_seg          : number of windows
#   mean_sim       : mean of upper-triangle of sim_seg (off-diagonal)
#   blue_frac      : fraction of off-diag Z cells with z < -1
#   red_frac       : fraction of off-diag Z cells with z > +1
#   blue_red_ratio : blue_frac / (blue_frac + red_frac), in [0, 1]
#   K_best         : silhouette-best K in 2..5, or 1 if no valid K
#   dominant_frac  : fraction of windows in the largest cluster
#   n_runs         : rle alternation count of cluster IDs
#   max_run        : longest contiguous same-cluster run
# Helper: build a signed per-distance z-matrix for a segment. Same logic
# as the inline block inside compute_segment_stats, factored out so the
# quadrant validator can reuse it without recomputing diag_stats twice.
# Returns z_mat clipped to [-z_clip, +z_clip], NA outside finite cells.
build_segment_z_mat <- function(sim_seg, z_clip = 5) {
  N_seg <- nrow(sim_seg)
  ds <- build_diag_stats(sim_seg)
  z_mat <- matrix(NA_real_, N_seg, N_seg)
  for (d in 0:(N_seg - 1L)) {
    mu <- ds$diag_mean[d + 1L]; sg <- ds$diag_sd[d + 1L]
    if (!is.finite(mu) || !is.finite(sg)) next
    if (d == 0L) {
      z_mat[cbind(seq_len(N_seg), seq_len(N_seg))] <- (diag(sim_seg) - mu) / sg
    } else {
      ii <- seq.int(1L, N_seg - d)
      z_vals <- (sim_seg[cbind(ii, ii + d)] - mu) / sg
      z_mat[cbind(ii, ii + d)] <- z_vals
      z_mat[cbind(ii + d, ii)] <- z_vals
    }
  }
  z_mat[!is.finite(z_mat)] <- NA_real_
  z_mat[z_mat >  z_clip] <-  z_clip
  z_mat[z_mat < -z_clip] <- -z_clip
  z_mat
}

compute_segment_stats <- function(sim_seg, parent_id, seg_index) {
  N_seg <- nrow(sim_seg)

  # mean_sim from upper triangle (matches handoff convention)
  ut <- sim_seg[upper.tri(sim_seg, diag = FALSE)]
  ut_f <- ut[is.finite(ut)]
  mean_sim <- if (length(ut_f) > 0L) mean(ut_f) else NA_real_

  # Local Z (per-distance standardization) — signed
  z_mat <- build_segment_z_mat(sim_seg, z_clip = 5)
  zv <- z_mat[upper.tri(z_mat, diag = FALSE)]
  zv <- zv[is.finite(zv)]
  if (length(zv) > 0L) {
    blue_frac <- mean(zv < -1)
    red_frac  <- mean(zv > +1)
    blue_red_ratio <- if ((blue_frac + red_frac) > 0)
      blue_frac / (blue_frac + red_frac) else NA_real_
  } else {
    blue_frac <- NA_real_; red_frac <- NA_real_; blue_red_ratio <- NA_real_
  }

  # Ward.D2-on-windows: row-correlation distance, hclust, silhouette over K
  K_best <- NA_integer_; dominant_frac <- NA_real_
  n_runs <- NA_integer_; max_run <- NA_integer_
  if (N_seg >= ward_min_segment_n) {
    # Subsample large segments for speed; cluster IDs are then mapped back
    # to the full segment by nearest-neighbour in row index (preserves
    # spatial contiguity for rle).
    if (N_seg > ward_max_n) {
      idx <- sort(sample.int(N_seg, ward_max_n))
      sub <- sim_seg[idx, idx]
    } else {
      idx <- seq_len(N_seg); sub <- sim_seg
    }
    sub[!is.finite(sub)] <- 0
    cm <- tryCatch(cor(sub, use = "pairwise.complete.obs"),
                   error = function(e) NULL)
    if (!is.null(cm) && all(is.finite(cm))) {
      d <- as.dist(1 - cm)
      hc <- tryCatch(hclust(d, method = "ward.D2"),
                     error = function(e) NULL)
      if (!is.null(hc)) {
        K_max <- min(5L, length(hc$order) - 1L)
        sil_K <- integer(0); sil_mean <- numeric(0)
        if (requireNamespace("cluster", quietly = TRUE) && K_max >= 2L) {
          for (Kk in 2:K_max) {
            cl <- cutree(hc, k = Kk)
            sw <- tryCatch(cluster::silhouette(cl, d),
                           error = function(e) NULL)
            if (!is.null(sw)) {
              sil_mean <- c(sil_mean, mean(sw[, "sil_width"], na.rm = TRUE))
              sil_K    <- c(sil_K, Kk)
            }
          }
        }
        if (length(sil_mean) > 0L && max(sil_mean, na.rm = TRUE) > 0.05) {
          K_best <- as.integer(sil_K[which.max(sil_mean)])
        } else {
          K_best <- 1L
        }
        cl_sub <- if (K_best >= 2L) cutree(hc, k = K_best) else
                  rep(1L, length(hc$order))
        # Map sub-cluster IDs back to full segment by nearest sub-index
        if (length(idx) == N_seg) {
          cl_full <- cl_sub
        } else {
          # for each full window i, find nearest j in idx
          nn <- vapply(seq_len(N_seg), function(i) {
            which.min(abs(idx - i))
          }, integer(1))
          cl_full <- cl_sub[nn]
        }
        tab <- table(cl_full)
        dominant_frac <- as.numeric(max(tab) / sum(tab))
        rl <- rle(cl_full)
        n_runs  <- as.integer(length(rl$lengths))
        max_run <- as.integer(max(rl$lengths))
      }
    }
  }

  data.table(
    parent_l1_id  = parent_id,
    seg_index     = seg_index,
    N_seg         = N_seg,
    mean_sim      = mean_sim,
    blue_frac     = blue_frac,
    red_frac      = red_frac,
    blue_red_ratio = blue_red_ratio,
    K_best        = K_best,
    dominant_frac = dominant_frac,
    n_runs        = n_runs,
    max_run       = max_run
  )
}

# Helper: map strictness in [-1, +1] to (real_mult, fake_mult).
# Convention reminder:
#   real_max is the upper bound of grow_max_z that classifies REAL.
#     LOWER real_max => fewer peaks pass as REAL (STRICT).
#   fake_min is the lower bound of grow_max_z that classifies FAKE.
#     LOWER fake_min => more peaks classify FAKE (STRICT).
# So in STRICT mode (s > 0): real_mult < 1, fake_mult < 1.
# In LOOSE mode (s < 0):     real_mult > 1, fake_mult > 1.
# Strict scaling uses ward_ceil_max; loose scaling of real_max uses
# ward_strictness_max so loosening can be more aggressive than tightening
# (matches Quentin's "rather see more cuts than miss real ones" preference).
strictness_to_mults <- function(s) {
  if (!is.finite(s)) return(list(real_mult = 1, fake_mult = 1))
  s <- max(-1, min(1, s))
  if (s >= 0) {
    real_mult <- 1 - ward_ceil_max * s         # s=+1 -> x(1-ceil_max)
    fake_mult <- 1 - ward_ceil_max * s         # symmetric tightening
  } else {
    real_mult <- 1 - ward_strictness_max * s   # s=-1 -> x(1+strictness_max)
    fake_mult <- 1 - ward_ceil_max * s         # s=-1 -> x(1+ceil_max)
  }
  list(real_mult = max(0, real_mult),
       fake_mult = max(0, fake_mult))
}

# ---- Quadrant validator ------------------------------------------------------
# For a given peak at local position bw_local, build the upper-left local
# square [bw_local - W .. bw_local] x [bw_local .. bw_local + W] of the
# Z-matrix, split it into a 2x2 grid, and return per-quadrant mean z.
#
# Real boundary signature: at multiple W scales, the 4 quadrant means are
# close to each other (homogeneity) AND the overall mean does not drift
# strongly positive (no "trending toward red"). The per-W test passes if
# both:
#   max(q1..q4) - min(q1..q4) <= quad_homog_thresh   (homogeneity)
#   mean(q1..q4)              <= quad_drift_thresh   (no red drift)
# The peak is "homogeneous" if the per-W test passes at >= quad_min_pass
# of the W values in the sweep.
#
# Returns a list: $w_used (vector of W values that fit inside the segment),
# $homog_per_w, $drift_per_w, $pass_per_w (booleans), $n_pass, $verdict
# (one of HOMOG_PASS / HOMOG_FAIL / NO_ROOM / NO_NEIGHBOR_ROOM / NO_W_FIT / NA).
# NO_ROOM is set when the peak is closer than min_room to either edge of
# its parent segment -- the local square would be too small for the test
# to discriminate. NO_NEIGHBOR_ROOM is set when another peak in the same
# segment is closer than min_neighbor_dist windows -- in dense peak clusters
# the audit's surrounding context is dominated by the cluster itself, not
# the true architectural neighborhood, so the homogeneity test is unreliable.
# NO_W_FIT is set when no W in the sweep fits inside the segment around
# the peak (rare, only on extremely small segments).
quadrant_audit_peak <- function(z_mat, bw_local,
                                w_sweep,
                                homog_thresh,
                                drift_thresh,
                                min_pass,
                                min_room,
                                neighbors = integer(0),
                                min_neighbor_dist = 0L) {
  N_seg <- nrow(z_mat)
  out <- list(
    w_used = integer(0),
    homog_per_w = numeric(0),
    drift_per_w = numeric(0),
    pass_per_w  = logical(0),
    n_pass      = 0L,
    verdict     = "NA"
  )
  if (!is.finite(bw_local) || bw_local < 1L || bw_local > N_seg) return(out)

  # Neighbor-peak check: if another peak in the same segment is within
  # min_neighbor_dist windows, the local square at this peak overlaps
  # neighbouring peaks' surroundings and the homogeneity test is
  # unreliable. neighbors should NOT include bw_local itself.
  if (length(neighbors) > 0L && min_neighbor_dist > 0L) {
    nearest <- min(abs(as.integer(neighbors) - bw_local))
    if (is.finite(nearest) && nearest < min_neighbor_dist) {
      out$verdict <- "NO_NEIGHBOR_ROOM"
      return(out)
    }
  }

  # Room check: peak must have at least min_room windows on BOTH sides
  # within its parent segment, otherwise the audit can't run a meaningful
  # test (the local square is too small for sub-quadrants to discriminate
  # signal from noise). Caller distinguishes NO_ROOM from NO_W_FIT for
  # diagnostic purposes.
  room_left  <- bw_local - 1L
  room_right <- N_seg - bw_local
  if (min(room_left, room_right) < min_room) {
    out$verdict <- "NO_ROOM"
    return(out)
  }

  for (W in w_sweep) {
    # The local square is [bw_local - W .. bw_local] (rows) x
    # [bw_local + 1 .. bw_local + W] (cols), i.e. the upper-right
    # off-diagonal block of a cross at the anchor. Use upper-triangle
    # geometry so all cells are well-defined per-distance Z.
    r_lo <- bw_local - W; r_hi <- bw_local
    c_lo <- bw_local + 1L; c_hi <- bw_local + W
    # Skip W if the square doesn't fit inside the segment
    if (r_lo < 1L || c_hi > N_seg) next
    # Sub-z block (W rows x W cols, where applicable)
    sq <- z_mat[r_lo:r_hi, c_lo:c_hi, drop = FALSE]
    if (any(dim(sq) < 2L)) next  # need at least 2x2 to split into quadrants

    half_r <- floor(W / 2)
    half_c <- floor(W / 2)
    if (half_r < 1L || half_c < 1L) next

    # 2x2 sub-quadrants of the WxW square. Use the floor() split point;
    # any leftover row/col goes to the second half. Naming:
    # q1 = top-left, q2 = top-right, q3 = bottom-left, q4 = bottom-right.
    q1 <- sq[1:half_r, 1:half_c, drop = FALSE]
    q2 <- sq[1:half_r, (half_c + 1L):ncol(sq), drop = FALSE]
    q3 <- sq[(half_r + 1L):nrow(sq), 1:half_c, drop = FALSE]
    q4 <- sq[(half_r + 1L):nrow(sq), (half_c + 1L):ncol(sq), drop = FALSE]

    means <- c(
      q1 = mean(q1, na.rm = TRUE),
      q2 = mean(q2, na.rm = TRUE),
      q3 = mean(q3, na.rm = TRUE),
      q4 = mean(q4, na.rm = TRUE)
    )
    if (any(!is.finite(means))) next

    homog_W <- max(means) - min(means)
    drift_W <- mean(means)
    pass_W  <- (homog_W <= homog_thresh) && (drift_W <= drift_thresh)

    out$w_used     <- c(out$w_used, W)
    out$homog_per_w <- c(out$homog_per_w, homog_W)
    out$drift_per_w <- c(out$drift_per_w, drift_W)
    out$pass_per_w  <- c(out$pass_per_w, pass_W)
  }

  if (length(out$w_used) == 0L) {
    out$verdict <- "NO_W_FIT"
    return(out)
  }
  out$n_pass <- sum(out$pass_per_w)
  out$verdict <- if (out$n_pass >= min_pass) "HOMOG_PASS" else "HOMOG_FAIL"
  out
}

cat("\n[D17L2] === SCANNING L1 SEGMENTS ===\n")
all_bdt_list   <- vector("list", n_l1_segments)
seg_stats_list <- vector("list", n_l1_segments)
seg_n_list     <- vector("list", n_l1_segments)
for (k in seq_len(n_l1_segments)) {
  row <- l1_dt[k]
  s <- as.integer(row$start_w); e <- as.integer(row$end_w)
  pid <- as.character(row$candidate_id)
  N_seg <- e - s + 1L
  if (N_seg < 2L * boundary_W + boundary_offset + 5L) {
    cat(sprintf("[D17L2]   %s  [%d..%d]  N=%d  (too small to scan, skipping)\n",
                pid, s, e, N_seg))
    next
  }
  # Track N per segment regardless of ward_adaptive (used by per-segment
  # min_largest_W scaling downstream)
  seg_n_list[[k]] <- data.table(parent_l1_id = pid, N_seg = N_seg)
  cat(sprintf("[D17L2]   %s  [%d..%d]  N=%d  (%.2f-%.2f Mb)\n",
              pid, s, e, N_seg,
              window_start_bp[s] / 1e6, window_end_bp[e] / 1e6))
  sim_seg <- sim_mat[s:e, s:e]

  # Ward+intensity stats (only when adaptive is on)
  if (ward_adaptive) {
    st <- compute_segment_stats(sim_seg, pid, k)
    seg_stats_list[[k]] <- st
    cat(sprintf("[D17L2]     stats: mean_sim=%.3f  blue_red=%s  K=%s  dom=%s  n_runs=%s\n",
                st$mean_sim,
                ifelse(is.na(st$blue_red_ratio), "NA",
                       sprintf("%.3f", st$blue_red_ratio)),
                ifelse(is.na(st$K_best), "NA", as.character(st$K_best)),
                ifelse(is.na(st$dominant_frac), "NA",
                       sprintf("%.3f", st$dominant_frac)),
                ifelse(is.na(st$n_runs), "NA", as.character(st$n_runs))))
  }

  bdt_seg <- scan_segment(sim_seg, pid, abs_offset = s - 1L)
  if (is.null(bdt_seg) || nrow(bdt_seg) == 0L) {
    cat("[D17L2]     no peaks\n"); next
  }
  cat(sprintf("[D17L2]     peaks: %d\n", nrow(bdt_seg)))

  # Per-peak quadrant audit (only when --quadrant_validator yes)
  if (quadrant_validator) {
    z_mat_seg <- build_segment_z_mat(sim_seg, z_clip = 5)
    n_pk <- nrow(bdt_seg)
    quad_verdict   <- character(n_pk)
    quad_n_pass    <- integer(n_pk)
    quad_n_w_used  <- integer(n_pk)
    quad_min_homog <- numeric(n_pk)
    quad_min_drift <- numeric(n_pk)
    # All peaks' local positions in this segment, used for the
    # neighbor-room check inside quadrant_audit_peak.
    all_bw_local <- as.integer(bdt_seg$boundary_w_local)
    for (pk in seq_len(n_pk)) {
      bw_local <- as.integer(bdt_seg$boundary_w_local[pk])
      # Other peaks in this segment (excluding the current one)
      neighbor_pos <- all_bw_local[-pk]
      qa <- quadrant_audit_peak(z_mat_seg, bw_local,
                                w_sweep            = quad_w_sweep,
                                homog_thresh       = quad_homog_thresh,
                                drift_thresh       = quad_drift_thresh,
                                min_pass           = quad_min_pass,
                                min_room           = quad_min_room,
                                neighbors          = neighbor_pos,
                                min_neighbor_dist  = quad_min_neighbor_dist)
      quad_verdict[pk]   <- qa$verdict
      quad_n_pass[pk]    <- qa$n_pass
      quad_n_w_used[pk]  <- length(qa$w_used)
      quad_min_homog[pk] <- if (length(qa$homog_per_w) > 0L)
                             min(qa$homog_per_w) else NA_real_
      quad_min_drift[pk] <- if (length(qa$drift_per_w) > 0L)
                             min(qa$drift_per_w) else NA_real_
    }
    bdt_seg[, quad_verdict   := quad_verdict]
    bdt_seg[, quad_n_pass    := quad_n_pass]
    bdt_seg[, quad_n_w_used  := quad_n_w_used]
    bdt_seg[, quad_min_homog := quad_min_homog]
    bdt_seg[, quad_min_drift := quad_min_drift]
    qtab <- table(quad_verdict)
    cat("[D17L2]     quadrant audit: ",
        paste(names(qtab), qtab, sep = "=", collapse = "  "),
        "\n", sep = "")
  }

  all_bdt_list[[k]] <- bdt_seg
}
all_bdt_list <- all_bdt_list[!vapply(all_bdt_list, is.null, logical(1))]
seg_stats_list <- seg_stats_list[!vapply(seg_stats_list, is.null, logical(1))]
seg_stats <- if (length(seg_stats_list) > 0L)
  rbindlist(seg_stats_list, fill = TRUE) else data.table()
seg_n_list <- seg_n_list[!vapply(seg_n_list, is.null, logical(1))]
seg_n_table <- if (length(seg_n_list) > 0L)
  rbindlist(seg_n_list, fill = TRUE) else data.table()
if (length(all_bdt_list) == 0L) {
  cat("\n[D17L2] no L2 peaks found in any L1 segment -- writing empty outputs\n")
  bdt_all <- data.table()
} else {
  bdt_all <- rbindlist(all_bdt_list, fill = TRUE)
  bdt_all[, boundary_idx := sprintf("%s.L2_b%04d", chr_label,
                                    seq_len(.N))]
  cat(sprintf("\n[D17L2] total L2 peaks across all segments: %d\n",
              nrow(bdt_all)))
}

# ---- Pooled baseline thresholds + (optional) per-segment Ward-adaptive ------
# When ward_adaptive=yes, the pooled (baseline_real_max, baseline_fake_min)
# computed below is the BASELINE per-chromosome threshold. Each L1 segment
# then receives its own (real_max, fake_min) by multiplying the baseline
# with strictness-derived multipliers from rank-normalised priors.

if (nrow(bdt_all) > 0L && boundary_validator_mode %in% c("grow", "both")) {
  finite_gmz <- bdt_all$grow_max_z
  finite_gmz <- finite_gmz[is.finite(finite_gmz)]
  use_adaptive <- (boundary_grow_threshold_mode == "adaptive") &&
                  (length(finite_gmz) >= boundary_grow_min_n_for_adaptive)
  if (use_adaptive) {
    q_real <- as.numeric(quantile(finite_gmz, boundary_grow_real_pct))
    q_fake <- as.numeric(quantile(finite_gmz, boundary_grow_fake_pct))
    baseline_real_max <- min(q_real, boundary_grow_real_max_ceiling)
    baseline_fake_min <- max(q_fake, boundary_grow_fake_min_floor)
    if (baseline_fake_min <= baseline_real_max)
      baseline_fake_min <- baseline_real_max + 0.05
    cat(sprintf("[D17L2] POOLED baseline thresholds (n=%d): q%.2f=%+.3f -> real_max=%+.3f  q%.2f=%+.3f -> fake_min=%+.3f\n",
                length(finite_gmz),
                boundary_grow_real_pct, q_real, baseline_real_max,
                boundary_grow_fake_pct, q_fake, baseline_fake_min))
  } else {
    if (boundary_grow_threshold_mode == "adaptive")
      cat(sprintf("[D17L2] adaptive requested but only %d finite peaks (need %d). Using fixed.\n",
                  length(finite_gmz), boundary_grow_min_n_for_adaptive))
    baseline_real_max <- boundary_grow_real_max
    baseline_fake_min <- boundary_grow_fake_min
    cat(sprintf("[D17L2] fixed baseline thresholds: real_max=%+.3f  fake_min=%+.3f\n",
                baseline_real_max, baseline_fake_min))
  }

  # ---- Per-segment Ward-adaptive thresholds -------------------------------
  if (ward_adaptive && nrow(seg_stats) > 0L) {
    cat("\n[D17L2] === ADAPTIVE PER-SEGMENT TUNING ===\n")

    # Rank-normalise mean_sim across segments (chromosome-relative). Higher
    # rank => higher value. NA-safe: NAs stay NA.
    rank01 <- function(x) {
      out <- rep(NA_real_, length(x))
      f <- is.finite(x)
      if (sum(f) >= 2L) {
        r <- rank(x[f], ties.method = "average")
        out[f] <- (r - 1) / (sum(f) - 1)
      } else if (sum(f) == 1L) {
        out[f] <- 0.5
      }
      out
    }
    seg_stats[, mean_sim_rank := rank01(mean_sim)]

    # Replace NAs in core inputs with neutral 0.5 so a segment that failed
    # Ward (e.g. too small) gets neutral priors and falls back to baseline.
    seg_stats[, dominant_frac_safe := ifelse(is.finite(dominant_frac),
                                             dominant_frac, 0.5)]
    seg_stats[, blue_red_safe := ifelse(is.finite(blue_red_ratio),
                                        blue_red_ratio, 0.5)]
    seg_stats[, mean_sim_rank_safe := ifelse(is.finite(mean_sim_rank),
                                             mean_sim_rank, 0.5)]
    seg_stats[, n_runs_per_w := ifelse(is.finite(n_runs) & N_seg > 0,
                                       n_runs / N_seg, NA_real_)]

    # K_clean_score: Ward says clean nesting when K>=2 and runs ~= K
    # (each cluster spans one contiguous block). High n_runs/W penalises.
    seg_stats[, K_clean_score := {
      out <- rep(0.5, .N)
      for (i in seq_len(.N)) {
        Kk <- K_best[i]; rr <- n_runs_per_w[i]
        if (is.na(Kk) || is.na(rr)) { out[i] <- 0.5; next }
        if (Kk == 1L) { out[i] <- 0.0; next }
        # Best when n_runs_per_w is small. >5/W of fragmentation -> 0.
        out[i] <- max(0, min(1, 1 - 5 * rr))
      }
      out
    }]

    # scaled_n_runs: tent peaking at n_runs/W ≈ 0.02 (a few transitions per
    # 100 windows = real architecture). Falls off on both sides — too low
    # means uniform; too high means fragmented noise.
    seg_stats[, scaled_n_runs := {
      out <- rep(0.5, .N)
      for (i in seq_len(.N)) {
        r <- n_runs_per_w[i]
        if (is.na(r)) { out[i] <- 0.5; next }
        out[i] <- if (r <= 0.02) min(1, r / 0.02)
                  else max(0, 1 - (r - 0.02) / 0.20)
      }
      out
    }]

    # Two priors in [0, 1]
    seg_stats[, internal_homogeneity := pmin(1, pmax(0,
      (mean_sim_rank_safe + dominant_frac_safe + (1 - blue_red_safe)) / 3))]
    seg_stats[, architectural_complexity := pmin(1, pmax(0,
      (scaled_n_runs + K_clean_score + blue_red_safe) / 3))]

    seg_stats[, strictness := pmax(-1, pmin(1,
                                  internal_homogeneity - architectural_complexity))]
    mults <- lapply(seg_stats$strictness, strictness_to_mults)
    seg_stats[, real_mult := vapply(mults, `[[`, numeric(1), "real_mult")]
    seg_stats[, fake_mult := vapply(mults, `[[`, numeric(1), "fake_mult")]

    seg_stats[, seg_real_max := baseline_real_max * real_mult]
    seg_stats[, seg_fake_min := baseline_fake_min * fake_mult]
    # Sanity: keep gap >= 0.05 between real_max and fake_min after scaling
    seg_stats[seg_fake_min <= seg_real_max + 0.02,
              seg_fake_min := seg_real_max + 0.05]

    for (i in seq_len(nrow(seg_stats))) {
      s_i <- seg_stats$strictness[i]
      label <- if (s_i > 0.20) "TIGHTENED"
               else if (s_i < -0.20) "LOOSENED"
               else "NEUTRAL"
      cat(sprintf("[D17L2] %s  N=%d  mean_sim=%.3f(rk=%.2f)  blue_red=%.3f  K=%s dom=%.2f n_runs=%s\n",
                  seg_stats$parent_l1_id[i], seg_stats$N_seg[i],
                  seg_stats$mean_sim[i], seg_stats$mean_sim_rank_safe[i],
                  seg_stats$blue_red_safe[i],
                  ifelse(is.na(seg_stats$K_best[i]), "NA",
                         as.character(seg_stats$K_best[i])),
                  seg_stats$dominant_frac_safe[i],
                  ifelse(is.na(seg_stats$n_runs[i]), "NA",
                         as.character(seg_stats$n_runs[i]))))
      cat(sprintf("[D17L2]   internal_homog=%.2f  arch_complex=%.2f  strictness=%+.2f  [%s]\n",
                  seg_stats$internal_homogeneity[i],
                  seg_stats$architectural_complexity[i],
                  s_i, label))
      cat(sprintf("[D17L2]   real_max: %+.3f -> %+.3f (x%.2f)   fake_min: %+.3f -> %+.3f (x%.2f)\n",
                  baseline_real_max, seg_stats$seg_real_max[i],
                  seg_stats$real_mult[i],
                  baseline_fake_min, seg_stats$seg_fake_min[i],
                  seg_stats$fake_mult[i]))
    }

    # Map each peak to its segment's thresholds via merge on parent_l1_id
    thr_map <- seg_stats[, .(parent_l1_id, seg_real_max, seg_fake_min)]
    bdt_all <- merge(bdt_all, thr_map, by = "parent_l1_id",
                     all.x = TRUE, sort = FALSE)
    # Fall back to baseline for any peak whose segment has no stats
    bdt_all[is.na(seg_real_max), seg_real_max := baseline_real_max]
    bdt_all[is.na(seg_fake_min), seg_fake_min := baseline_fake_min]
    bdt_all[, resolved_real_max := seg_real_max]
    bdt_all[, resolved_fake_min := seg_fake_min]
    bdt_all[, c("seg_real_max","seg_fake_min") := NULL]
  } else {
    bdt_all[, resolved_real_max := baseline_real_max]
    bdt_all[, resolved_fake_min := baseline_fake_min]
  }

  # Per-segment min_largest_W: on small L1 segments the validator's max
  # achievable W is small, so a global floor of 10 auto-EDGEs every peak.
  # Scale the floor down on small segments (cap at the global floor for
  # large ones). Formula: max(5, min(global_floor, floor(N_seg * 0.10))).
  if (nrow(seg_n_table) > 0L) {
    seg_min_lW <- seg_n_table[, .(parent_l1_id,
                                  seg_min_largest_W = pmax(5L,
                                    pmin(boundary_grow_min_largest_W,
                                         as.integer(floor(N_seg * 0.10)))))]
    bdt_all <- merge(bdt_all, seg_min_lW, by = "parent_l1_id",
                     all.x = TRUE, sort = FALSE)
    bdt_all[is.na(seg_min_largest_W),
            seg_min_largest_W := boundary_grow_min_largest_W]
    cat("[D17L2] per-segment min_largest_W floors:\n")
    for (i in seq_len(nrow(seg_min_lW))) {
      cat(sprintf("[D17L2]   %s  min_largest_W=%d\n",
                  seg_min_lW$parent_l1_id[i],
                  seg_min_lW$seg_min_largest_W[i]))
    }
  } else {
    bdt_all[, seg_min_largest_W := boundary_grow_min_largest_W]
  }

  classify_grow_row <- function(gmax, glargest, rmax, fmin, lW_min) {
    if (!is.finite(gmax) || !is.finite(glargest)) return("EDGE")
    if (glargest < lW_min)                        return("EDGE")
    if (gmax <= rmax)                             return("REAL")
    if (gmax >= fmin)                             return("FAKE")
    "MARGINAL"
  }
  bdt_all[, grow_status := mapply(classify_grow_row,
                                  grow_max_z, grow_largest_W,
                                  resolved_real_max, resolved_fake_min,
                                  seg_min_largest_W)]
  remap <- c("REAL" = "STABLE_BLUE", "FAKE" = "DECAYS",
             "MARGINAL" = "MARGINAL", "EDGE" = "EDGE")
  bdt_all[, validation_status := unname(remap[grow_status])]
  bdt_all[, seg_min_largest_W := NULL]

  tab <- table(bdt_all$validation_status)
  cat("[D17L2] validation summary: ",
      paste(names(tab), tab, sep = "=", collapse = "  "), "\n", sep = "")

  # ---- Quadrant rescue + flag pass --------------------------------------
  # Apply quadrant verdicts (already computed per-peak during the scan loop)
  # to validation_status:
  #   - DECAYS / FAKE peak with HOMOG_PASS => RESCUE to STABLE_BLUE
  #     (sets quadrant_rescued = TRUE) -- subject to the rescue cap on
  #     grow_max_z: peaks with grow_max_z above the cap stay DECAYS even
  #     if quadrant says HOMOG_PASS (cross is too red-leaning to trust).
  #   - STABLE_BLUE peak with HOMOG_FAIL    =>
  #       quad_demote_on_fail=TRUE  (default) DEMOTE to MARGINAL
  #       quad_demote_on_fail=FALSE (legacy)  keep status, flag
  #       Either way: quadrant_failed = TRUE for downstream audit.
  # MARGINAL peaks with HOMOG_PASS are also rescued to STABLE_BLUE for
  # consistency (they were on the bubble; quadrant says they're real),
  # also subject to the rescue cap.
  # EDGE peaks are never rescued (validator failure means we have no
  # confidence in the kernel output for them).
  if (quadrant_validator && "quad_verdict" %in% names(bdt_all)) {
    bdt_all[, quadrant_rescued    := FALSE]
    bdt_all[, quadrant_failed     := FALSE]
    bdt_all[, quadrant_rescue_blocked_by_cap := FALSE]
    bdt_all[, quadrant_demoted    := FALSE]
    # Rescue candidates: DECAYS or MARGINAL that the quadrant rated PASS.
    rescue_cand_idx <- which(bdt_all$validation_status %in% c("DECAYS", "MARGINAL") &
                             bdt_all$quad_verdict == "HOMOG_PASS")
    # Apply rescue cap on grow_max_z
    if (length(rescue_cand_idx) > 0L && "grow_max_z" %in% names(bdt_all)) {
      gmz_vals <- bdt_all$grow_max_z[rescue_cand_idx]
      cap_pass <- is.finite(gmz_vals) & gmz_vals <= quad_rescue_max_grow_z
      rescue_idx        <- rescue_cand_idx[cap_pass]
      rescue_blocked_idx <- rescue_cand_idx[!cap_pass]
    } else {
      rescue_idx <- rescue_cand_idx
      rescue_blocked_idx <- integer(0)
    }
    flag_idx <- which(bdt_all$validation_status == "STABLE_BLUE" &
                      bdt_all$quad_verdict == "HOMOG_FAIL")
    if (length(rescue_idx) > 0L) {
      bdt_all[rescue_idx, quadrant_rescued := TRUE]
      bdt_all[rescue_idx, validation_status := "STABLE_BLUE"]
    }
    if (length(rescue_blocked_idx) > 0L) {
      bdt_all[rescue_blocked_idx, quadrant_rescue_blocked_by_cap := TRUE]
      # validation_status stays as DECAYS or MARGINAL
    }
    if (length(flag_idx) > 0L) {
      bdt_all[flag_idx, quadrant_failed := TRUE]
      if (quad_demote_on_fail) {
        # Drift-floor protection: don't demote peaks where the quadrant
        # mean drift is strongly negative (deep blue). Their HOMOG_FAIL
        # is between flavours of "blue" (e.g. -1 vs -3), which isn't
        # biologically meaningful.
        drift_vals <- bdt_all$quad_min_drift[flag_idx]
        # Treat NA drift as "no protection" (demote)
        protect <- is.finite(drift_vals) & drift_vals <= quad_demote_drift_floor
        demote_idx <- flag_idx[!protect]
        protected_idx <- flag_idx[protect]
        if (length(demote_idx) > 0L) {
          bdt_all[demote_idx, validation_status := "MARGINAL"]
          bdt_all[demote_idx, quadrant_demoted := TRUE]
        }
        # protected_idx keep STABLE_BLUE; quadrant_failed already TRUE
        n_demoted <- length(demote_idx)
        n_drift_protected <- length(protected_idx)
      } else {
        n_demoted <- 0L
        n_drift_protected <- 0L
      }
    } else {
      n_demoted <- 0L
      n_drift_protected <- 0L
    }
    cat(sprintf(
      "[D17L2] quadrant pass: rescued=%d  blocked_by_cap=%d  flagged_failed=%d  demoted=%d  drift_protected=%d\n",
      length(rescue_idx), length(rescue_blocked_idx),
      length(flag_idx),
      n_demoted, n_drift_protected))
    if (length(rescue_idx) > 0L || length(flag_idx) > 0L ||
        length(rescue_blocked_idx) > 0L) {
      tab2 <- table(bdt_all$validation_status)
      cat("[D17L2] post-quadrant validation summary: ",
          paste(names(tab2), tab2, sep = "=", collapse = "  "),
          "\n", sep = "")
    }
  }

  # ---- Weak-peak demote pass --------------------------------------------
  # Demote STABLE_BLUE peaks where BOTH boundary_score and grow_max_z are
  # weak. Catches kernel-ambiguous peaks the quadrant validator couldn't
  # adjudicate (NO_NEIGHBOR_ROOM / NO_ROOM in clusters). Set
  # weak_demote_score to 0 to disable this pass entirely.
  if (weak_demote_score > 0 && nrow(bdt_all) > 0L &&
      "boundary_score" %in% names(bdt_all) &&
      "grow_max_z" %in% names(bdt_all)) {
    if (!"weak_demoted" %in% names(bdt_all))
      bdt_all[, weak_demoted := FALSE]
    weak_idx <- which(bdt_all$validation_status == "STABLE_BLUE" &
                      is.finite(bdt_all$boundary_score) &
                      is.finite(bdt_all$grow_max_z) &
                      bdt_all$boundary_score < weak_demote_score &
                      bdt_all$grow_max_z   > weak_demote_gmz)
    if (length(weak_idx) > 0L) {
      bdt_all[weak_idx, validation_status := "MARGINAL"]
      bdt_all[weak_idx, weak_demoted := TRUE]
    }
    cat(sprintf("[D17L2] weak-peak demote pass: %d STABLE -> MARGINAL (score<%g AND gmz>%g)\n",
                length(weak_idx), weak_demote_score, weak_demote_gmz))
    if (length(weak_idx) > 0L) {
      tab3 <- table(bdt_all$validation_status)
      cat("[D17L2] post-weak-demote summary: ",
          paste(names(tab3), tab3, sep = "=", collapse = "  "),
          "\n", sep = "")
    }
  }

  # Per-L1-segment validation breakdown so tuning blind-spots are visible
  cat("[D17L2] per-segment validation breakdown:\n")
  per_seg <- bdt_all[, .(
    raw         = .N,
    STABLE_BLUE = sum(validation_status == "STABLE_BLUE"),
    MARGINAL    = sum(validation_status == "MARGINAL"),
    DECAYS      = sum(validation_status == "DECAYS"),
    EDGE        = sum(validation_status == "EDGE")
  ), by = parent_l1_id][order(parent_l1_id)]
  for (i in seq_len(nrow(per_seg))) {
    cat(sprintf("[D17L2]   %s  raw=%d  STABLE=%d  MARGINAL=%d  DECAYS=%d  EDGE=%d\n",
                per_seg$parent_l1_id[i], per_seg$raw[i],
                per_seg$STABLE_BLUE[i], per_seg$MARGINAL[i],
                per_seg$DECAYS[i], per_seg$EDGE[i]))
  }

  # ---- Dedup pass: drop near-duplicate STABLE_BLUE peaks ------------------
  # The scan kernel can fire twice on the same architectural edge at
  # slightly different offsets / smoothing scales. We keep the higher-
  # scoring peak per cluster of nearby peaks within boundary_dedup_dist
  # windows. Done per-L1-segment so segments are independent. Non-STABLE
  # peaks are left untouched (they're not used downstream anyway).
  if (boundary_dedup_dist > 0L &&
      "validation_status" %in% names(bdt_all) &&
      any(bdt_all$validation_status == "STABLE_BLUE")) {
    bdt_all[, dedup_dropped := FALSE]
    n_dropped_total <- 0L
    pids <- unique(bdt_all$parent_l1_id)
    for (pid in pids) {
      # Indices of STABLE peaks in this segment, sorted by score desc
      ids_seg_stable <- which(bdt_all$parent_l1_id == pid &
                              bdt_all$validation_status == "STABLE_BLUE")
      if (length(ids_seg_stable) <= 1L) next
      ord <- ids_seg_stable[order(-bdt_all$boundary_score[ids_seg_stable])]
      kept_w <- integer(0)
      dropped_here <- integer(0)
      for (idx in ord) {
        bw <- bdt_all$boundary_w[idx]
        if (length(kept_w) > 0L &&
            min(abs(kept_w - bw)) < boundary_dedup_dist) {
          dropped_here <- c(dropped_here, idx)
        } else {
          kept_w <- c(kept_w, bw)
        }
      }
      if (length(dropped_here) > 0L) {
        bdt_all[dropped_here, dedup_dropped := TRUE]
        bdt_all[dropped_here, validation_status := "DEDUP"]
        n_dropped_total <- n_dropped_total + length(dropped_here)
      }
    }
    cat(sprintf("[D17L2] dedup pass (within %d windows, per-segment, score-prioritised): %d STABLE peaks reclassified to DEDUP\n",
                boundary_dedup_dist, n_dropped_total))
  }
} else if (nrow(bdt_all) > 0L) {
  bdt_all[, validation_status := perp_status]
}

# ---- Write per-segment stats TSV (consumed by overlay/audit later) ----------
if (!dry_run && ward_adaptive && nrow(seg_stats) > 0L) {
  out_s <- file.path(outdir, paste0(chr_label, ".L2_segment_stats.tsv"))
  fwrite(seg_stats, out_s, sep = "\t")
  cat(sprintf("[D17L2] segment stats written: %s  (rows=%d)\n",
              out_s, nrow(seg_stats)))
}

# ---- Write quadrant audit sidecar ------------------------------------------
if (!dry_run && quadrant_validator && nrow(bdt_all) > 0L &&
    "quad_verdict" %in% names(bdt_all)) {
  out_q <- file.path(outdir, paste0(chr_label, ".L2_quadrant_audit.tsv"))
  cols <- intersect(c("chr","boundary_idx","parent_l1_id",
                      "boundary_w_local","boundary_w","boundary_bp",
                      "boundary_score","grow_max_z","grow_largest_W",
                      "validation_status",
                      "quad_verdict","quad_n_pass","quad_n_w_used",
                      "quad_min_homog","quad_min_drift",
                      "quadrant_rescued","quadrant_failed",
                      "quadrant_rescue_blocked_by_cap","quadrant_demoted",
                      "weak_demoted"),
                    names(bdt_all))
  fwrite(bdt_all[, ..cols], out_q, sep = "\t")
  cat(sprintf("[D17L2] quadrant audit written: %s  (rows=%d)\n",
              out_q, nrow(bdt_all)))
}

# ---- Write boundaries -------------------------------------------------------

out_b <- file.path(outdir, paste0(chr_label, ".L2_boundaries.tsv"))
if (!dry_run) {
  if (nrow(bdt_all) > 0L) {
    setcolorder(bdt_all, c(
      "chr", "boundary_idx", "parent_l1_id",
      "boundary_w_local", "boundary_w", "boundary_bp",
      "boundary_score", "boundary_W", "boundary_offset",
      intersect(c("right_frac_blue","left_frac_blue","right_max_z","left_max_z",
                  "right_first_red","left_first_red","right_n_finite",
                  "left_n_finite","perp_status"), names(bdt_all)),
      intersect(c(grep("^cross_z_W", names(bdt_all), value = TRUE),
                  "grow_max_z","grow_largest_W","grow_z_at_largest",
                  "grow_status","resolved_real_max","resolved_fake_min"),
                names(bdt_all)),
      "validation_status"
    ))
    fwrite(bdt_all, out_b, sep = "\t")
    cat("[D17L2] L2 boundaries written: ", out_b, "\n", sep = "")
  } else {
    fwrite(data.table(), out_b, sep = "\t")
    cat("[D17L2] L2 boundaries (empty) written: ", out_b, "\n", sep = "")
  }
}

# ---- Derive L2 segments per L1 parent ---------------------------------------

if (!dry_run) {
  cut_keep_set <- c("STABLE_BLUE")
  l2_segs <- vector("list", n_l1_segments)
  for (k in seq_len(n_l1_segments)) {
    row <- l1_dt[k]
    s <- as.integer(row$start_w); e <- as.integer(row$end_w)
    pid <- as.character(row$candidate_id)

    if (nrow(bdt_all) > 0L && "validation_status" %in% names(bdt_all)) {
      cuts_dt <- bdt_all[parent_l1_id == pid &
                         validation_status %in% cut_keep_set]
    } else {
      cuts_dt <- bdt_all[FALSE]
    }
    cuts <- sort(unique(as.integer(cuts_dt$boundary_w)))
    cuts <- cuts[is.finite(cuts) & cuts >= s & cuts < e]

    if (length(cuts) == 0L) {
      seg_starts <- s; seg_ends <- e
    } else {
      seg_starts <- c(s, cuts + 1L)
      seg_ends   <- c(cuts, e)
      keep <- seg_ends >= seg_starts
      seg_starts <- seg_starts[keep]; seg_ends <- seg_ends[keep]
    }
    n_sub <- length(seg_starts)
    sub_dt <- data.table(
      chr           = chr_label,
      parent_l1_id  = pid,
      candidate_id  = sprintf("%s.L2_%04d_%02d", chr_label, k, seq_len(n_sub)),
      start_w       = seg_starts,
      end_w         = seg_ends,
      start_bp      = window_start_bp[seg_starts],
      end_bp        = window_end_bp[seg_ends],
      n_windows     = seg_ends - seg_starts + 1L,
      scale_W       = seg_ends - seg_starts + 1L,
      status        = "L2_ENVELOPE"
    )
    sub_dt[, mean_sim := vapply(seq_len(.N), function(j) {
      ss <- start_w[j]; ee <- end_w[j]
      if (ee <= ss) return(NA_real_)
      mean(sim_mat[ss:ee, ss:ee], na.rm = TRUE)
    }, numeric(1))]
    sub_dt[, density_p70 := vapply(seq_len(.N), function(j) {
      ss <- start_w[j]; ee <- end_w[j]
      if (ee <= ss) return(NA_real_)
      mean(sim_mat[ss:ee, ss:ee] >= 0.70, na.rm = TRUE)
    }, numeric(1))]
    l2_segs[[k]] <- sub_dt
  }
  l2_dt <- rbindlist(l2_segs, fill = TRUE)
  out_e <- file.path(outdir, paste0(chr_label, ".L2_envelopes.tsv"))
  fwrite(l2_dt, out_e, sep = "\t")
  cat(sprintf("[D17L2] L2 envelopes written (%d segments across %d L1 parents): %s\n",
              nrow(l2_dt), n_l1_segments, out_e))

  cat("[D17L2] L2 segments per L1 parent:\n")
  per_parent <- l2_dt[, .N, by = parent_l1_id]
  for (i in seq_len(nrow(per_parent))) {
    cat(sprintf("[D17L2]   %s -> %d sub-segments\n",
                per_parent$parent_l1_id[i], per_parent$N[i]))
  }
}

cat("[D17L2] done\n")
