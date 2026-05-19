#!/usr/bin/env Rscript
# =============================================================================
# test_lib_persample_cusum.R — observation-only tests
# =============================================================================
# Verifies the per-sample CUSUM observations only:
#   - cp_bp lands near the synthetic changepoint
#   - strength threshold separates signal from noise
#   - asymmetry sign is correct
#   - NA windows are tolerated row-by-row
# Cohort-level aggregation (modes, dispersion, multimodality) is the consensus
# step's job and is intentionally NOT tested here.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source("/home/claude/work/phase_4_resolution/shared_lib/lib_persample_cusum.R")

set.seed(42)

.failures <- list()
.successes <- 0L

ok   <- function(msg) { cat("[lib]   ok :", msg, "\n"); .successes <<- .successes + 1L }
fail <- function(msg) { cat("[lib] FAIL :", msg, "\n"); .failures[[length(.failures) + 1L]] <<- msg }
near <- function(a, b, tol) abs(a - b) <= tol

N_S    <- 226L
N_W    <- 100L
WIN_BP <- (seq_len(N_W) - 1L) * 10000L + 5000L
sample_ids <- sprintf("CGA%04d", seq_len(N_S))

# =============================================================================
# Test 1: SHARP INVERSION — 60 carriers shift at window 50, the rest noise
# =============================================================================
{
  M <- matrix(rnorm(N_S * N_W, 0, 0.1), nrow = N_S, ncol = N_W,
              dimnames = list(sample_ids, NULL))
  for (i in 1:60) M[i, 50:N_W] <- M[i, 50:N_W] + 1.0

  ps <- persample_cusum(M, WIN_BP, candidate_id = "TEST_SHARP")

  # Output shape
  if (nrow(ps) != N_S) fail(sprintf("Test 1: nrow=%d, expected %d", nrow(ps), N_S))
  else ok(sprintf("Test 1 sharp: persample table has %d rows", nrow(ps)))

  expected_cols <- c("sample_id", "candidate_id", "cp_idx", "cp_bp",
                     "strength", "asymmetry", "left_mean", "right_mean",
                     "n_used", "informative")
  miss <- setdiff(expected_cols, names(ps))
  if (length(miss) > 0) fail(sprintf("Test 1: missing cols %s", paste(miss, collapse=",")))
  else ok("Test 1 sharp: all 10 columns present")

  # Carriers should be informative, non-carriers should not
  inf_carriers      <- ps$informative[1:60]
  inf_noncarriers   <- ps$informative[61:226]
  if (sum(inf_carriers) < 55) {
    fail(sprintf("Test 1: only %d of 60 carriers informative", sum(inf_carriers)))
  } else {
    ok(sprintf("Test 1 sharp: %d/60 carriers informative", sum(inf_carriers)))
  }
  if (sum(inf_noncarriers) > 5) {
    fail(sprintf("Test 1: %d/166 non-carriers leaked through threshold (target ~0)",
                 sum(inf_noncarriers)))
  } else {
    ok(sprintf("Test 1 sharp: only %d/166 noise rows passed threshold", sum(inf_noncarriers)))
  }

  # Carriers' cp_bp should land near WIN_BP[50]
  expected_bp <- WIN_BP[50]
  carrier_cps <- ps$cp_bp[1:60]
  carrier_cps <- carrier_cps[is.finite(carrier_cps)]
  median_cp <- median(carrier_cps)
  if (!near(median_cp, expected_bp, 30000)) {
    fail(sprintf("Test 1: median carrier cp_bp = %d, expected ~%d (±30 kb)",
                 median_cp, expected_bp))
  } else {
    ok(sprintf("Test 1 sharp: median carrier cp_bp = %d (target ~%d)",
               median_cp, expected_bp))
  }
}

# =============================================================================
# Test 2: ASYMMETRY SIGN — half rise, half fall at the same position
# =============================================================================
{
  M <- matrix(rnorm(N_S * N_W, 0, 0.1), nrow = N_S, ncol = N_W,
              dimnames = list(sample_ids, NULL))
  for (i in 1:30)  M[i, 50:N_W] <- M[i, 50:N_W] + 1.0    # rises
  for (i in 31:60) M[i, 50:N_W] <- M[i, 50:N_W] - 1.0    # falls

  ps <- persample_cusum(M, WIN_BP, candidate_id = "TEST_ASYM")
  inf <- ps[informative == TRUE]

  n_pos <- sum(inf$asymmetry == 1L,  na.rm = TRUE)
  n_neg <- sum(inf$asymmetry == -1L, na.rm = TRUE)

  if (!near(n_pos, 30, 5) || !near(n_neg, 30, 5)) {
    fail(sprintf("Test 2: asymmetry counts pos=%d neg=%d (target ~30/30)", n_pos, n_neg))
  } else {
    ok(sprintf("Test 2 asym: pos=%d, neg=%d (correctly distinguishes rise vs fall)",
               n_pos, n_neg))
  }
}

# =============================================================================
# Test 3: PURE NOISE — almost no rows should pass the strength threshold
# =============================================================================
{
  M <- matrix(rnorm(N_S * N_W, 0, 0.1), nrow = N_S, ncol = N_W,
              dimnames = list(sample_ids, NULL))
  ps <- persample_cusum(M, WIN_BP, candidate_id = "TEST_NULL")
  n_inf <- sum(ps$informative)
  if (n_inf > 5) {
    fail(sprintf("Test 3: %d noise rows pass threshold (expected <5)", n_inf))
  } else {
    ok(sprintf("Test 3 null: %d/226 noise rows pass threshold (clean separation)", n_inf))
  }
}

# =============================================================================
# Test 4: PARTIAL NA — random NA windows on 30 of 60 carriers
# =============================================================================
{
  M <- matrix(rnorm(N_S * N_W, 0, 0.1), nrow = N_S, ncol = N_W,
              dimnames = list(sample_ids, NULL))
  for (i in 1:60) M[i, 50:N_W] <- M[i, 50:N_W] + 1.0
  for (i in 1:30) {
    na_pos <- sample(1:N_W, 30)
    M[i, na_pos] <- NA_real_
  }

  ps <- persample_cusum(M, WIN_BP, candidate_id = "TEST_NA")
  n_inf_carriers <- sum(ps$informative[1:60])
  if (n_inf_carriers < 50) {
    fail(sprintf("Test 4: only %d/60 carriers informative with NA tolerance",
                 n_inf_carriers))
  } else {
    ok(sprintf("Test 4 NA: %d/60 carriers informative (NA tolerance works)",
               n_inf_carriers))
  }

  # cp_bp should still land near WIN_BP[50]
  carrier_cps <- ps$cp_bp[1:60]
  carrier_cps <- carrier_cps[is.finite(carrier_cps)]
  if (length(carrier_cps) > 0) {
    med <- median(carrier_cps)
    if (!near(med, WIN_BP[50], 50000)) {
      fail(sprintf("Test 4 NA: median cp_bp %d not near %d", med, WIN_BP[50]))
    } else {
      ok(sprintf("Test 4 NA: median cp_bp %d (target ~%d)", med, WIN_BP[50]))
    }
  }
}

# =============================================================================
# Summary
# =============================================================================
cat("\n================================================\n")
if (length(.failures) == 0) {
  cat(sprintf("[lib] ALL %d CHECKS PASSED\n", .successes))
  quit(status = 0)
} else {
  cat(sprintf("[lib] %d passed, %d FAILED:\n", .successes, length(.failures)))
  for (f in .failures) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
