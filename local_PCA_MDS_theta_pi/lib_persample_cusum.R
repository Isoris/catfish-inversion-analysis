# =============================================================================
# lib_persample_cusum.R
# =============================================================================
# Phase 4 / shared_lib — pure CUSUM utility, observation-only.
#
# DESIGN PHILOSOPHY:
#   This library produces per-carrier observations, NOT cohort-level summaries.
#   We do not fit kernel densities, do not run distribution-shape tests, and
#   do not classify spread into named regimes. Those steps presuppose a model
#   for the distribution of per-carrier changepoints, and we don't have one —
#   the distribution is whatever biology produces (founder block, ragged tail,
#   bimodal, delta + uniform, etc).
#
#   The right place to combine observations into a breakpoint call is the
#   downstream consensus step (phase_4 / 4c, formerly phase_6, also formerly
#   `breakpoint_pipeline/`), which already gathers per-carrier evidence from
#   multiple streams (ancestral fragments, regime transitions, SV calls) and
#   takes the modal position across the pooled evidence — not per-stream
#   pre-aggregates.
#
#   See breakpoint_pipeline/docs/METHODOLOGY.md and the per-carrier mode
#   approach in 02_ancestral_fragments.R / 03_consensus_merge.R for the
#   established framework that this CUSUM module feeds into as one more
#   evidence stream alongside the ancestral-fragment scanner.
#
# INPUT:  a sample × window numeric matrix (rows = samples, cols = windows),
#         restricted to a candidate's window range. NAs tolerated per row.
# OUTPUT: a single data.table with one row per sample:
#           sample_id, candidate_id, cp_idx, cp_bp, strength, asymmetry,
#           left_mean, right_mean, n_used, informative
#
# Math: one-sided CUSUM with mean-removal (Lancaster MATH337 §1.3). For each
# row x[1..W]:
#   1. Restrict to non-NA entries.
#   2. Subtract row mean → centred series x_c.
#   3. cumsum(x_c) gives the running deviation.
#   4. Changepoint position = argmax_w |cumsum(x_c)[w]|.
#   5. Strength = max|cumsum(x_c)| / sd(x_c) — scale-invariant statistic.
#   6. Asymmetry sign = sign(cumsum(x_c)[w*]) — +1 if signal rises past the
#      changepoint, -1 if it falls.
#
# `informative` flag: scales as sqrt(n_used). Under pure noise, max|cumsum|/sd
# scales as sqrt(n) (random-walk variance). Empirically (calibrate2.R quantiles
# at 0.5 / 0.95 / 0.99 / 0.999 = 0.79 / 1.29 / 1.59 / 1.93), a threshold of
# 2.0 × sqrt(n) keeps roughly 1-in-1000 noise rows; real signal carriers score
# ~5× sqrt(n). This is empirical, not a model fit.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Strength threshold scaling: see header for derivation.
# Under pure noise, max|cumsum|/sd ~ sqrt(n). Real signal scores ~5× sqrt(n).
# A threshold of 2.0 × sqrt(n_used) admits <1-in-1000 noise rows.
.MIN_STRENGTH_PER_SQRT_N <- 2.0

# =============================================================================
# Per-sample CUSUM on one row
# =============================================================================
.cusum_one_row <- function(x, window_pos_bp) {
  # x: numeric vector length W, may contain NAs
  # window_pos_bp: integer vector length W giving the BP coord of each window
  ok <- is.finite(x)
  if (sum(ok) < 5L) {
    return(list(
      cp_idx        = NA_integer_,
      cp_bp         = NA_integer_,
      strength      = NA_real_,
      asymmetry     = NA_integer_,
      left_mean     = NA_real_,
      right_mean    = NA_real_,
      n_used        = sum(ok)
    ))
  }
  x_ok   <- x[ok]
  pos_ok <- window_pos_bp[ok]
  m      <- mean(x_ok)
  s      <- sd(x_ok)
  if (!is.finite(s) || s <= 0) {
    return(list(
      cp_idx     = NA_integer_, cp_bp = NA_integer_,
      strength   = NA_real_,    asymmetry = NA_integer_,
      left_mean  = m,            right_mean = m,
      n_used     = length(x_ok)
    ))
  }
  cs <- cumsum(x_ok - m)
  k  <- which.max(abs(cs))
  cp_idx <- which(ok)[k]                 # index back into the original W
  list(
    cp_idx     = as.integer(cp_idx),
    cp_bp      = as.integer(window_pos_bp[cp_idx]),
    strength   = abs(cs[k]) / s,
    asymmetry  = as.integer(sign(cs[k])),
    left_mean  = mean(x_ok[seq_len(k)]),
    right_mean = if (k < length(x_ok)) mean(x_ok[(k + 1L):length(x_ok)]) else NA_real_,
    n_used     = length(x_ok)
  )
}

# =============================================================================
# Public entry point
# =============================================================================
# matrix:         sample × window numeric, rownames = sample IDs
# window_pos_bp:  integer length ncol(matrix), the BP coordinate of each
#                 window (typically the midpoint or start_bp)
# candidate_id:   string, propagated to outputs for join keys
# Returns:        a single data.table with one row per sample.
# =============================================================================
persample_cusum <- function(matrix,
                            window_pos_bp,
                            candidate_id = NA_character_) {

  stopifnot(is.matrix(matrix))
  stopifnot(length(window_pos_bp) == ncol(matrix))
  if (is.null(rownames(matrix))) {
    stop("persample_cusum: matrix must have rownames (sample IDs)")
  }

  n_samp <- nrow(matrix)
  n_win  <- ncol(matrix)
  sample_ids <- rownames(matrix)

  # ── Per-sample CUSUM ────────────────────────────────────────────────────
  rows <- vector("list", n_samp)
  for (i in seq_len(n_samp)) {
    r <- .cusum_one_row(matrix[i, ], window_pos_bp)
    # Strength threshold scales with sqrt(n_used) — see .MIN_STRENGTH_PER_SQRT_N
    # comment for derivation. Falls back to FALSE if strength or n_used unknown.
    is_informative <- is.finite(r$strength) &&
                      is.finite(r$n_used) && r$n_used >= 5L &&
                      r$strength >= .MIN_STRENGTH_PER_SQRT_N * sqrt(r$n_used)
    rows[[i]] <- data.table(
      sample_id     = sample_ids[i],
      candidate_id  = candidate_id,
      cp_idx        = r$cp_idx,
      cp_bp         = r$cp_bp,
      strength      = r$strength,
      asymmetry     = r$asymmetry,
      left_mean     = r$left_mean,
      right_mean    = r$right_mean,
      n_used        = r$n_used,
      informative   = is_informative
    )
  }
  persample_dt <- rbindlist(rows)
  persample_dt
}
