#!/usr/bin/env Rscript
# =============================================================================
# STEP_D17_boundary_detect_L1.R  (v1, 2026-05-12)
# =============================================================================
# Boundary-driven Level-1 envelope discovery on a single chromosome.
#
# PATH-AGNOSTIC. The algorithm operates entirely on the sim_mat and reads
# only $dt$start_bp / $dt$end_bp from the precomp. It works on:
#
#   - θπ:   STEP_TR_C output  (<chr>.precomp.rds         + sim_mats/<chr>.sim_mat_nn{N}.rds)
#   - GHSL: STEP_GH_C output  (<chr>.ghsl_precomp.rds    + sim_mats/<chr>.sim_mat_nn{N}.rds)
#   - any future per-window-PC1 path that produces the same artifacts.
#
# Algorithm (verbatim from the path-1 z-blocks D17 v7-adaptive):
#
#   1. Diagonal cross-block scan. At each window position i, the score is
#      -median(z) over the WxW upper-triangle cross-block at offset G,
#      where z is per-diagonal-distance Z-normalized similarity. Peaks
#      mark candidate boundaries.
#   2. Grow-W validator (default). Re-evaluates the cross-block at a
#      percent-of-N W ladder. REAL boundaries keep median(z) negative or
#      near zero even at large W; FAKE drifts positive. Adaptive
#      thresholds (with ceiling/floor guards) calibrate per chromosome.
#   3. Optional perpendicular-ray validator. Two 1D rays from each peak
#      perpendicular to the diagonal: REAL stays mostly blue; FAKE turns
#      red close to the diagonal.
#   4. STABLE_BLUE peaks partition the chromosome into L1 envelopes.
#      Tiny segments (≤ --l1_min_segment_nw windows) are merged into the
#      more-similar neighbor by block-mean sim.
#
# Path-routing knobs (NEW vs the original D17):
#
#   --precomp_suffix <str>  filename suffix to resolve <chr>.<suffix> in
#                           --precomp_dir. Default: "precomp.rds" (θπ).
#                           For GHSL: pass "ghsl_precomp.rds".
#   --label <str>           short tag (default "theta") prefixed to log
#                           lines and output filenames so outputs from
#                           different paths can coexist in one outdir.
#                           For GHSL: pass "ghsl".
#
# All other knobs are unchanged from the original D17 — see header in
# the path-1 04_detect_L1_localpca_zblocks.R for the full menu and the
# settled-defaults rationale.
#
# Outputs (in --outdir, with --label prefix):
#   <label>_<chr>.L1_envelopes.tsv      L1 segments (partition from boundaries)
#   <label>_<chr>.L1_boundaries.tsv     every detected peak + validation
#   <label>_<chr>.L1_score_curve.tsv    per-window boundary score curve
#
# Run examples:
#
#   # θπ path (coarse precomp from STEP_TR_C --mode full)
#   Rscript STEP_D17_boundary_detect_L1.R \
#     --precomp_dir $OUTROOT/precomp_coarse \
#     --chr C_gar_LG28 --nn 80 \
#     --outdir $OUTROOT/L1_theta \
#     --label theta
#
#   # GHSL path (coarse precomp from STEP_GH_C --mode full)
#   Rscript STEP_D17_boundary_detect_L1.R \
#     --precomp_dir $OUTROOT/ghsl_precomp_coarse \
#     --chr C_gar_LG28 --nn 80 \
#     --outdir $OUTROOT/L1_ghsl \
#     --label ghsl --precomp_suffix ghsl_precomp.rds
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# ---- CLI --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

# Path-routing
precomp_dir     <- get_arg("--precomp_dir")
precomp_suffix  <- get_arg("--precomp_suffix", "precomp.rds")
label           <- get_arg("--label", "theta")
precomp_f       <- get_arg("--precomp")
sim_mat_f       <- get_arg("--sim_mat")
chr_label       <- get_arg("--chr", "chr")
outdir          <- get_arg("--outdir", ".")
nn_scale        <- as.integer(get_arg("--nn", "80"))

# Auto-resolve precomp + sim_mat from --precomp_dir + --chr + --nn
if (is.na(precomp_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, paste0(chr_label, ".", precomp_suffix))
  if (file.exists(cand)) {
    precomp_f <- cand
    cat("[", label, ".D17L1] resolved --precomp: ", cand, "\n", sep = "")
  }
}
if (is.na(sim_mat_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, "sim_mats",
                    paste0(chr_label, ".sim_mat_nn", nn_scale, ".rds"))
  if (file.exists(cand)) {
    sim_mat_f <- cand
    cat("[", label, ".D17L1] resolved --sim_mat (nn=", nn_scale, "): ",
        cand, "\n", sep = "")
  }
}

# Detection parameters
boundary_scan          <- !is.na(get_arg("--boundary_scan", NA_character_))
# default-on: if --boundary_scan not given as a string-with-arg, treat the
# original "presence" semantics. For convenience, also default the
# behavior to TRUE since the script's purpose IS the boundary scan.
if (!boundary_scan) boundary_scan <- TRUE
boundary_W             <- as.integer(get_arg("--boundary_W", "5"))
boundary_offset        <- as.integer(get_arg("--boundary_offset", "5"))
boundary_score_min     <- as.numeric(get_arg("--boundary_score_min", "2.0"))

l1_boundary_filter     <- tolower(get_arg("--l1_boundary_filter", "stable"))
l1_min_segment_nw      <- as.integer(get_arg("--l1_min_segment_nw", "1"))

# min_dist resolution priority
boundary_min_dist_arg     <- get_arg("--boundary_min_dist", NA_character_)
boundary_min_dist_kb_arg  <- get_arg("--boundary_min_dist_kb", NA_character_)
boundary_min_dist_pct_arg <- get_arg("--boundary_min_dist_pct", NA_character_)
boundary_min_dist_floor   <- as.integer(get_arg("--boundary_min_dist_floor", "5"))
boundary_min_dist         <- NA_integer_

# Perp-ray validator
boundary_validate              <- get_arg("--boundary_validate", "TRUE")
boundary_validate              <- !(toupper(boundary_validate) %in% c("FALSE","F","0","NO"))
boundary_perp_d_max            <- as.integer(get_arg("--boundary_perp_d_max", "20"))
boundary_perp_min_blue_frac    <- as.numeric(get_arg("--boundary_perp_min_blue_frac", "0.70"))
boundary_perp_max_red          <- as.numeric(get_arg("--boundary_perp_max_red", "0.50"))
boundary_perp_first_d_red_max  <- as.integer(get_arg("--boundary_perp_first_d_red_max", "5"))
boundary_perp_red_z            <- as.numeric(get_arg("--boundary_perp_red_z", "0.50"))

# Grow-W validator
boundary_validator_mode        <- tolower(get_arg("--boundary_validator_mode", "grow"))
boundary_grow_W_arg            <- get_arg("--boundary_grow_W", NA_character_)
boundary_grow_W_pct_str        <- get_arg("--boundary_grow_W_pct", "0.001,0.005,0.01,0.02,0.04,0.07")
boundary_grow_W_pct            <- as.numeric(strsplit(boundary_grow_W_pct_str, ",")[[1]])
boundary_grow_W_pct            <- sort(unique(boundary_grow_W_pct[
  is.finite(boundary_grow_W_pct) & boundary_grow_W_pct > 0]))
boundary_grow_W_floor          <- as.integer(get_arg("--boundary_grow_W_floor", "5"))

boundary_grow_threshold_mode   <- tolower(get_arg("--boundary_grow_threshold_mode", "adaptive"))
boundary_grow_real_pct         <- as.numeric(get_arg("--boundary_grow_real_pct", "0.50"))
boundary_grow_fake_pct         <- as.numeric(get_arg("--boundary_grow_fake_pct", "0.50"))
boundary_grow_real_max_ceiling <- as.numeric(get_arg("--boundary_grow_real_max_ceiling", "0.10"))
boundary_grow_fake_min_floor   <- as.numeric(get_arg("--boundary_grow_fake_min_floor", "0.10"))
boundary_grow_min_n_for_adaptive <- as.integer(get_arg("--boundary_grow_min_n_for_adaptive", "10"))
boundary_grow_real_max         <- as.numeric(get_arg("--boundary_grow_real_max", "0.10"))
boundary_grow_fake_min         <- as.numeric(get_arg("--boundary_grow_fake_min", "0.20"))
boundary_grow_min_largest_W    <- as.integer(get_arg("--boundary_grow_min_largest_W", "20"))

dry_run        <- !is.na(get_arg("--dry_run", NA_character_))

if (is.na(precomp_f) || !file.exists(precomp_f))
  stop("[", label, ".D17L1] --precomp could not be resolved or does not exist")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

LP <- paste0("[", label, ".D17L1] ")
cat(LP, "precomp:    ", precomp_f, "\n", sep = "")
cat(LP, "sim_mat:    ", sim_mat_f %||% "(auto)", "\n", sep = "")
cat(LP, "chr:        ", chr_label, "\n", sep = "")
cat(LP, "label:      ", label, "\n", sep = "")
cat(LP, "outdir:     ", outdir, "\n", sep = "")
cat(LP, "filter=", l1_boundary_filter,
    " min_seg_nw=", l1_min_segment_nw, "\n", sep = "")

# ---- Load precomp -----------------------------------------------------------
cat(LP, "loading precomp\n", sep = "")
pc_obj <- readRDS(precomp_f)
if (!is.null(pc_obj$dt))      pc <- pc_obj
else if (!is.null(pc_obj$pc)) pc <- pc_obj$pc
else stop(LP, "cannot find $dt in precomp file")

dt_pc <- as.data.table(pc$dt)
n_windows_total <- nrow(dt_pc)

# Coordinate columns — accept any of the standard schemas
if (all(c("start_bp", "end_bp") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start_bp; window_end_bp <- dt_pc$end_bp
} else if (all(c("start", "end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start;    window_end_bp <- dt_pc$end
} else if (all(c("window_start", "window_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$window_start; window_end_bp <- dt_pc$window_end
} else if (all(c("bp_start", "bp_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$bp_start; window_end_bp <- dt_pc$bp_end
} else {
  stop(LP, "precomp dt missing usable coordinate columns. Found: ",
       paste(names(dt_pc), collapse = ", "))
}

cat(LP, "N windows: ", n_windows_total, "\n", sep = "")

# ---- Resolve boundary_grow_W (needs N) -------------------------------------
if (!is.na(boundary_grow_W_arg)) {
  boundary_grow_W <- as.integer(strsplit(boundary_grow_W_arg, ",")[[1]])
  boundary_grow_W <- boundary_grow_W[is.finite(boundary_grow_W)]
} else {
  boundary_grow_W <- as.integer(round(boundary_grow_W_pct * n_windows_total))
}
boundary_grow_W <- sort(unique(pmax(boundary_grow_W, boundary_grow_W_floor)))
if (boundary_validator_mode %in% c("grow", "both")) {
  cat(LP, "grow W resolved (N=", n_windows_total, "): ",
      paste(boundary_grow_W, collapse = ","), "\n", sep = "")
}

# ---- Resolve boundary_min_dist (needs window sizes) ------------------------
window_size_bp_med <- as.numeric(median(window_end_bp - window_start_bp + 1L, na.rm = TRUE))
if (!is.finite(window_size_bp_med) || window_size_bp_med <= 0) window_size_bp_med <- NA_real_
if (!is.na(boundary_min_dist_arg)) {
  boundary_min_dist <- as.integer(boundary_min_dist_arg)
  resolution_note <- sprintf("absolute override (%d windows)", boundary_min_dist)
} else if (!is.na(boundary_min_dist_kb_arg)) {
  kb <- as.numeric(boundary_min_dist_kb_arg)
  if (!is.finite(window_size_bp_med)) stop(LP, "cannot resolve _kb: window size unknown")
  boundary_min_dist <- as.integer(round(kb * 1000 / window_size_bp_med))
  resolution_note <- sprintf("%g kb / %.0f bp/win = %d windows", kb, window_size_bp_med, boundary_min_dist)
} else if (!is.na(boundary_min_dist_pct_arg)) {
  pct <- as.numeric(boundary_min_dist_pct_arg)
  boundary_min_dist <- as.integer(round(pct * n_windows_total))
  resolution_note <- sprintf("legacy %g of N=%d -> %d windows", pct, n_windows_total, boundary_min_dist)
} else {
  default_kb <- 150
  if (!is.finite(window_size_bp_med)) {
    boundary_min_dist <- 30L; resolution_note <- "default 30 windows (window size unknown)"
  } else {
    boundary_min_dist <- as.integer(round(default_kb * 1000 / window_size_bp_med))
    resolution_note <- sprintf("default %d kb / %.0f bp/win = %d windows",
                               default_kb, window_size_bp_med, boundary_min_dist)
  }
}
boundary_min_dist <- max(boundary_min_dist, boundary_min_dist_floor)
cat(LP, "min_dist: ", resolution_note, "  (final=", boundary_min_dist, ")\n", sep = "")

# ---- Load sim_mat -----------------------------------------------------------
if (is.na(sim_mat_f)) {
  candidates_path <- c(
    file.path(dirname(precomp_f), "sim_mat_nn80.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn160.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn40.rds")
  )
  hit <- candidates_path[file.exists(candidates_path)]
  if (length(hit) == 0L) stop(LP, "--sim_mat not given and no auto-find candidates exist")
  sim_mat_f <- hit[1]
  cat(LP, "auto-found sim_mat: ", sim_mat_f, "\n", sep = "")
}
cat(LP, "loading sim_mat\n", sep = "")
sm_obj <- readRDS(sim_mat_f)
if (is.matrix(sm_obj)) {
  sim_mat <- sm_obj
} else if (is.list(sm_obj) && !is.null(sm_obj$sim_mat) && is.matrix(sm_obj$sim_mat)) {
  sim_mat <- sm_obj$sim_mat
} else if (is.list(sm_obj) && length(sm_obj) == 1L && is.matrix(sm_obj[[1]])) {
  sim_mat <- sm_obj[[1]]
} else {
  stop(LP, "sim_mat object structure not recognized; class=", class(sm_obj)[1])
}
storage.mode(sim_mat) <- "double"
if (!isTRUE(nrow(sim_mat) == n_windows_total && ncol(sim_mat) == n_windows_total)) {
  stop(LP, "sim_mat dim ", nrow(sim_mat), "x", ncol(sim_mat),
       " does not match precomp N ", n_windows_total)
}

# =============================================================================
# Boundary scan + validation (algorithm is verbatim from path-1 D17 v7-adaptive)
# =============================================================================

if (boundary_scan) {
  cat("\n", LP, "=== DIAGONAL BOUNDARY SCAN ===\n", sep = "")
  cat(LP, "W=", boundary_W, " offset=", boundary_offset,
      " score_min=", boundary_score_min, " min_dist=", boundary_min_dist, "\n", sep = "")

  N <- n_windows_total; W <- boundary_W; G <- boundary_offset

  cat(LP, "precomputing diagonal mean/sd for Z\n", sep = "")
  diag_mean <- numeric(N); diag_sd <- numeric(N)
  for (d in 0:(N - 1L)) {
    if (d == 0L) {
      vals <- diag(sim_mat)
    } else {
      ii <- seq.int(1L, N - d); vals <- sim_mat[cbind(ii, ii + d)]
    }
    vals <- vals[is.finite(vals)]
    if (length(vals) >= 5L) {
      diag_mean[d + 1L] <- mean(vals)
      sg <- sd(vals); diag_sd[d + 1L] <- if (is.finite(sg) && sg > 1e-9) sg else NA_real_
    } else {
      diag_mean[d + 1L] <- NA_real_; diag_sd[d + 1L] <- NA_real_
    }
  }

  i_lo <- W; i_hi <- N - G - W
  if (i_hi < i_lo) {
    cat(LP, "chromosome too short for this W/offset — skipping scan\n", sep = "")
  } else {
    centers <- seq.int(i_lo, i_hi)
    boundary_score <- rep(NA_real_, N)
    cat(LP, "scoring ", length(centers), " diagonal positions\n", sep = "")
    for (ctr in centers) {
      ii <- seq.int(ctr - W + 1L, ctr); jj <- seq.int(ctr + G + 1L, ctr + G + W)
      block <- sim_mat[ii, jj]
      d_mat <- outer(ii, jj, FUN = function(a, b) b - a)
      idx <- d_mat + 1L
      z <- (block - diag_mean[idx]) / diag_sd[idx]
      z[!is.finite(z)] <- NA_real_
      vals <- as.numeric(z); vals <- vals[is.finite(vals)]
      if (length(vals) >= 5L) boundary_score[ctr] <- -median(vals)
    }
    score_finite <- boundary_score[is.finite(boundary_score)]
    if (length(score_finite) > 0L) {
      cat(sprintf("%sboundary_score range: [%.2f, %.2f]  median=%.2f  p95=%.2f\n",
                  LP, min(score_finite), max(score_finite),
                  median(score_finite), as.numeric(quantile(score_finite, 0.95))))
    }

    # 1D peak detection
    peaks <- integer(0)
    for (ctr in centers) {
      s <- boundary_score[ctr]
      if (!is.finite(s) || s < boundary_score_min) next
      lo <- max(1L, ctr - boundary_min_dist); hi <- min(N, ctr + boundary_min_dist)
      win <- boundary_score[lo:hi]; win <- win[is.finite(win)]
      if (length(win) == 0L) next
      if (s >= max(win)) peaks <- c(peaks, ctr)
    }
    cat(LP, "peaks found: ", length(peaks), "\n", sep = "")

    perp_z <- function(r, c) {
      if (r < 1L || r > n_windows_total) return(NA_real_)
      if (c < 1L || c > n_windows_total) return(NA_real_)
      d <- abs(c - r)
      mu <- diag_mean[d + 1L]; sg <- diag_sd[d + 1L]; sm <- sim_mat[r, c]
      if (!is.finite(sm) || !is.finite(mu) || !is.finite(sg) || sg < 1e-9) return(NA_real_)
      (sm - mu) / sg
    }

    if (length(peaks) > 0L) {
      bdt <- data.table(
        chr            = chr_label,
        boundary_idx   = sprintf("%s.%s.L1_b%04d", label, chr_label, seq_along(peaks)),
        boundary_w     = peaks,
        boundary_bp    = window_start_bp[peaks],
        boundary_score = boundary_score[peaks],
        boundary_W     = W,
        boundary_offset = G
      )
      bdt <- bdt[order(boundary_w)]

      if (boundary_validate) {
        cat(LP, "perp-ray validation (d_max=", boundary_perp_d_max,
            ", red_z=", boundary_perp_red_z, ")\n", sep = "")
        n_peaks <- nrow(bdt)
        right_frac_blue <- numeric(n_peaks); left_frac_blue <- numeric(n_peaks)
        right_max_z     <- numeric(n_peaks); left_max_z     <- numeric(n_peaks)
        right_first_red <- integer(n_peaks); left_first_red <- integer(n_peaks)
        right_n_finite  <- integer(n_peaks); left_n_finite  <- integer(n_peaks)

        for (k in seq_len(n_peaks)) {
          i <- bdt$boundary_w[k]
          d_seq <- seq_len(boundary_perp_d_max)
          r_vals <- vapply(d_seq, function(d) perp_z(i, i + d), numeric(1))
          l_vals <- vapply(d_seq, function(d) perp_z(i - d, i), numeric(1))
          ok_r <- is.finite(r_vals); right_n_finite[k] <- sum(ok_r)
          if (any(ok_r)) {
            right_frac_blue[k] <- mean(r_vals[ok_r] < 0)
            right_max_z[k] <- max(r_vals[ok_r])
            red_idx <- which(r_vals > boundary_perp_red_z & ok_r)
            right_first_red[k] <- if (length(red_idx) > 0L) as.integer(red_idx[1]) else NA_integer_
          } else {
            right_frac_blue[k] <- NA_real_; right_max_z[k] <- NA_real_; right_first_red[k] <- NA_integer_
          }
          ok_l <- is.finite(l_vals); left_n_finite[k] <- sum(ok_l)
          if (any(ok_l)) {
            left_frac_blue[k] <- mean(l_vals[ok_l] < 0)
            left_max_z[k] <- max(l_vals[ok_l])
            red_idx <- which(l_vals > boundary_perp_red_z & ok_l)
            left_first_red[k] <- if (length(red_idx) > 0L) as.integer(red_idx[1]) else NA_integer_
          } else {
            left_frac_blue[k] <- NA_real_; left_max_z[k] <- NA_real_; left_first_red[k] <- NA_integer_
          }
        }
        bdt[, `:=`(right_frac_blue = right_frac_blue, left_frac_blue = left_frac_blue,
                   right_max_z = right_max_z, left_max_z = left_max_z,
                   right_first_red = right_first_red, left_first_red = left_first_red,
                   right_n_finite = right_n_finite, left_n_finite = left_n_finite)]

        classify_perp <- function(rfb, lfb, rmz, lmz, rfr, lfr) {
          if (!is.finite(rfb) || !is.finite(lfb)) return("EDGE")
          first_red_too_close <-
            (is.finite(rfr) && rfr <= boundary_perp_first_d_red_max) ||
            (is.finite(lfr) && lfr <= boundary_perp_first_d_red_max)
          if (first_red_too_close) return("DECAYS")
          stable <- (rfb >= boundary_perp_min_blue_frac) &&
                    (lfb >= boundary_perp_min_blue_frac) &&
                    (is.finite(rmz) && rmz <= boundary_perp_max_red) &&
                    (is.finite(lmz) && lmz <= boundary_perp_max_red)
          if (stable) return("STABLE_BLUE")
          "MARGINAL"
        }
        bdt[, perp_status := mapply(classify_perp,
                                    right_frac_blue, left_frac_blue,
                                    right_max_z, left_max_z,
                                    right_first_red, left_first_red)]
      }

      if (boundary_validator_mode %in% c("grow", "both") && length(boundary_grow_W) > 0L) {
        cat(LP, "grow-W validator: W=", paste(boundary_grow_W, collapse = ","),
            " thr_mode=", boundary_grow_threshold_mode, "\n", sep = "")
        n_peaks <- nrow(bdt)
        grow_mat <- matrix(NA_real_, nrow = n_peaks, ncol = length(boundary_grow_W))
        colnames(grow_mat) <- paste0("cross_z_W", boundary_grow_W)
        for (k in seq_len(n_peaks)) {
          i <- bdt$boundary_w[k]
          for (wj in seq_along(boundary_grow_W)) {
            Wg <- boundary_grow_W[wj]
            if ((i - Wg + 1L) < 1L || (i + G + Wg) > n_windows_total) next
            ii <- seq.int(i - Wg + 1L, i); jj <- seq.int(i + G + 1L, i + G + Wg)
            block <- sim_mat[ii, jj]
            d_mat <- outer(ii, jj, FUN = function(a, b) b - a); idx <- d_mat + 1L
            z <- (block - diag_mean[idx]) / diag_sd[idx]
            z[!is.finite(z)] <- NA_real_
            vals <- as.numeric(z); vals <- vals[is.finite(vals)]
            if (length(vals) >= 5L) grow_mat[k, wj] <- median(vals)
          }
        }
        for (wj in seq_along(boundary_grow_W)) {
          col <- colnames(grow_mat)[wj]; bdt[, (col) := grow_mat[, wj]]
        }
        grow_max_z <- numeric(n_peaks); grow_largest_W <- integer(n_peaks); grow_z_at_largest <- numeric(n_peaks)
        for (k in seq_len(n_peaks)) {
          row_z <- grow_mat[k, ]; ok <- which(is.finite(row_z))
          if (length(ok) == 0L) {
            grow_max_z[k] <- NA_real_; grow_largest_W[k] <- NA_integer_; grow_z_at_largest[k] <- NA_real_
          } else {
            grow_max_z[k] <- max(row_z[ok])
            grow_largest_W[k] <- as.integer(boundary_grow_W[max(ok)])
            grow_z_at_largest[k] <- row_z[max(ok)]
          }
        }
        bdt[, `:=`(grow_max_z = grow_max_z, grow_largest_W = grow_largest_W,
                   grow_z_at_largest = grow_z_at_largest)]

        # Threshold resolution
        finite_gmz <- grow_max_z[is.finite(grow_max_z)]
        use_adaptive <- (boundary_grow_threshold_mode == "adaptive") &&
                        (length(finite_gmz) >= boundary_grow_min_n_for_adaptive)
        if (use_adaptive) {
          q_real <- as.numeric(quantile(finite_gmz, boundary_grow_real_pct, na.rm = TRUE))
          q_fake <- as.numeric(quantile(finite_gmz, boundary_grow_fake_pct, na.rm = TRUE))
          resolved_real_max <- min(q_real, boundary_grow_real_max_ceiling)
          resolved_fake_min <- max(q_fake, boundary_grow_fake_min_floor)
          if (resolved_fake_min <= resolved_real_max) resolved_fake_min <- resolved_real_max + 0.05
          cat(sprintf("%sadaptive thresholds (n=%d peaks): real_max=%+.3f fake_min=%+.3f\n",
                      LP, length(finite_gmz), resolved_real_max, resolved_fake_min))
        } else {
          if (boundary_grow_threshold_mode == "adaptive")
            cat(LP, "adaptive requested but only ", length(finite_gmz), " finite peaks — falling back to fixed\n", sep = "")
          resolved_real_max <- boundary_grow_real_max; resolved_fake_min <- boundary_grow_fake_min
          cat(sprintf("%sfixed thresholds: real_max=%+.3f fake_min=%+.3f\n",
                      LP, resolved_real_max, resolved_fake_min))
        }
        bdt[, resolved_real_max := resolved_real_max]
        bdt[, resolved_fake_min := resolved_fake_min]

        classify_grow <- function(gmax, glargest) {
          if (!is.finite(gmax) || !is.finite(glargest)) return("EDGE")
          if (glargest < boundary_grow_min_largest_W) return("EDGE")
          if (gmax <= resolved_real_max) return("REAL")
          if (gmax >= resolved_fake_min) return("FAKE")
          "MARGINAL"
        }
        bdt[, grow_status := mapply(classify_grow, grow_max_z, grow_largest_W)]
      }

      # Final validation_status
      if (boundary_validator_mode == "perp") {
        bdt[, validation_status := if ("perp_status" %in% names(bdt)) perp_status else "STABLE_BLUE"]
      } else if (boundary_validator_mode == "grow") {
        bdt[, validation_status := if ("grow_status" %in% names(bdt))
              fifelse(grow_status == "REAL",     "STABLE_BLUE",
              fifelse(grow_status == "FAKE",     "DECAYS",
              fifelse(grow_status == "MARGINAL", "MARGINAL", "EDGE")))
            else "EDGE"]
      } else if (boundary_validator_mode == "both") {
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
          bdt[, validation_status := fifelse(grow_status == "REAL", "STABLE_BLUE",
                fifelse(grow_status == "FAKE", "DECAYS",
                fifelse(grow_status == "MARGINAL", "MARGINAL", "EDGE")))]
        } else if ("perp_status" %in% names(bdt)) {
          bdt[, validation_status := perp_status]
        } else bdt[, validation_status := "STABLE_BLUE"]
      } else {
        bdt[, validation_status := "STABLE_BLUE"]
      }

      n_stable <- sum(bdt$validation_status == "STABLE_BLUE")
      n_decay  <- sum(bdt$validation_status == "DECAYS")
      n_marg   <- sum(bdt$validation_status == "MARGINAL")
      n_edge   <- sum(bdt$validation_status == "EDGE")
      cat(LP, "validation summary: STABLE_BLUE=", n_stable,
          " DECAYS=", n_decay, " MARGINAL=", n_marg, " EDGE=", n_edge,
          " (mode=", boundary_validator_mode, ")\n", sep = "")
    } else {
      bdt <- data.table(
        chr = character(0), boundary_idx = character(0),
        boundary_w = integer(0), boundary_bp = integer(0),
        boundary_score = numeric(0),
        boundary_W = integer(0), boundary_offset = integer(0)
      )
    }

    if (!dry_run) {
      out_b <- file.path(outdir, paste0(label, "_", chr_label, ".L1_boundaries.tsv"))
      fwrite(bdt, out_b, sep = "\t")
      cat(LP, "boundaries: ", out_b, "\n", sep = "")
      curve_dt <- data.table(chr = chr_label, window_idx = seq_len(N),
                             bp = window_start_bp, boundary_score = boundary_score)
      out_c <- file.path(outdir, paste0(label, "_", chr_label, ".L1_score_curve.tsv"))
      fwrite(curve_dt, out_c, sep = "\t")
      cat(LP, "score curve: ", out_c, "\n", sep = "")
    }
  }
}

# ---- Derive L1 envelopes from STABLE_BLUE boundaries -----------------------
if (!dry_run && exists("bdt")) {
  cut_keep_set <- switch(l1_boundary_filter,
    "stable"    = c("STABLE_BLUE"),
    "non_decay" = c("STABLE_BLUE", "MARGINAL", "EDGE"),
    "all"       = c("STABLE_BLUE", "DECAYS", "MARGINAL", "EDGE"),
    c("STABLE_BLUE")
  )
  cut_dt <- if ("validation_status" %in% names(bdt)) bdt[validation_status %in% cut_keep_set] else bdt
  cuts <- sort(unique(as.integer(cut_dt$boundary_w)))
  cuts <- cuts[is.finite(cuts) & cuts >= 1L & cuts <= n_windows_total]

  out_main <- file.path(outdir, paste0(label, "_", chr_label, ".L1_envelopes.tsv"))
  if (length(cuts) == 0L) {
    cat(LP, "no surviving boundaries — writing empty envelope catalogue\n", sep = "")
    fwrite(data.table(
      chr = character(0), candidate_id = character(0),
      start_w = integer(0), end_w = integer(0),
      start_bp = integer(0), end_bp = integer(0),
      n_windows = integer(0), scale_W = integer(0),
      mean_sim = numeric(0), density_p70 = numeric(0),
      status = character(0)
    ), out_main, sep = "\t")
  } else {
    seg_starts <- c(1L, cuts + 1L); seg_ends <- c(cuts, n_windows_total)
    keep_seg <- seg_ends >= seg_starts
    seg_starts <- seg_starts[keep_seg]; seg_ends <- seg_ends[keep_seg]

    # Tiny-segment merge into more-similar neighbor
    n_merged <- 0L
    repeat {
      seg_nw <- seg_ends - seg_starts + 1L
      tiny_idx <- which(seg_nw <= l1_min_segment_nw)
      if (length(tiny_idx) == 0L) break
      k <- tiny_idx[1]
      tiny_s <- seg_starts[k]; tiny_e <- seg_ends[k]; tiny_n <- tiny_e - tiny_s + 1L
      has_left  <- k > 1L; has_right <- k < length(seg_starts)
      sim_left  <- if (has_left) {
        Ls <- seg_starts[k - 1L]; Le <- seg_ends[k - 1L]
        v <- as.numeric(sim_mat[tiny_s:tiny_e, Ls:Le]); v <- v[is.finite(v)]
        if (length(v) > 0L) mean(v) else NA_real_
      } else NA_real_
      sim_right <- if (has_right) {
        Rs <- seg_starts[k + 1L]; Re <- seg_ends[k + 1L]
        v <- as.numeric(sim_mat[tiny_s:tiny_e, Rs:Re]); v <- v[is.finite(v)]
        if (length(v) > 0L) mean(v) else NA_real_
      } else NA_real_
      merge_left <- if (!has_right) TRUE
                    else if (!has_left) FALSE
                    else if (!is.finite(sim_left) && !is.finite(sim_right)) TRUE
                    else if (!is.finite(sim_right)) TRUE
                    else if (!is.finite(sim_left)) FALSE
                    else (sim_left >= sim_right)
      if (merge_left) {
        seg_ends[k - 1L] <- seg_ends[k]
        seg_starts <- seg_starts[-k]; seg_ends <- seg_ends[-k]
      } else {
        seg_starts[k + 1L] <- seg_starts[k]
        seg_starts <- seg_starts[-k]; seg_ends <- seg_ends[-k]
      }
      n_merged <- n_merged + 1L
    }
    if (n_merged > 0L) cat(LP, "merged ", n_merged, " tiny segments\n", sep = "")
    n_seg <- length(seg_starts)

    seg_dt <- data.table(
      chr          = chr_label,
      candidate_id = sprintf("%s.%s.L1_%04d", label, chr_label, seq_len(n_seg)),
      start_w      = seg_starts,
      end_w        = seg_ends,
      start_bp     = window_start_bp[seg_starts],
      end_bp       = window_end_bp[seg_ends],
      n_windows    = seg_ends - seg_starts + 1L,
      scale_W      = seg_ends - seg_starts + 1L,
      status       = "ENVELOPE"
    )
    seg_dt[, mean_sim := vapply(seq_len(.N), function(k) {
      s <- start_w[k]; e <- end_w[k]; if (e <= s) return(NA_real_)
      mean(sim_mat[s:e, s:e], na.rm = TRUE)
    }, numeric(1))]
    seg_dt[, density_p70 := vapply(seq_len(.N), function(k) {
      s <- start_w[k]; e <- end_w[k]; if (e <= s) return(NA_real_)
      mean(sim_mat[s:e, s:e] >= 0.70, na.rm = TRUE)
    }, numeric(1))]

    setcolorder(seg_dt, c("chr","candidate_id","start_w","end_w","start_bp","end_bp",
                          "n_windows","scale_W","mean_sim","density_p70","status"))
    fwrite(seg_dt, out_main, sep = "\t")
    cat(LP, "L1 segments (n=", n_seg, ", cuts=", length(cuts),
        ", filter='", l1_boundary_filter, "'): ", out_main, "\n", sep = "")
    for (k in seq_len(n_seg)) {
      cat(sprintf("%s  %s  %.2f-%.2f Mb  nW=%d  mean_sim=%.3f\n",
                  LP, seg_dt$candidate_id[k],
                  seg_dt$start_bp[k] / 1e6, seg_dt$end_bp[k] / 1e6,
                  seg_dt$n_windows[k], seg_dt$mean_sim[k]))
    }
  }
}

cat(LP, "done\n", sep = "")
