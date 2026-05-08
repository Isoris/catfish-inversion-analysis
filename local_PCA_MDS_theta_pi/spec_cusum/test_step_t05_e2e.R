#!/usr/bin/env Rscript
# =============================================================================
# test_step_t05_e2e.R — end-to-end test of STEP_T05_theta_cusum.R
# =============================================================================
# Builds a synthetic TR_B-shaped JSON with a known signal, runs the driver,
# verifies the per-sample TSV looks correct.
#
# Synthetic setup:
#   226 samples, 100 windows each 10 kb wide (chrom ~1 Mb)
#   60 carriers shift their θπ at window 50 (= ~500 kb)
#   166 non-carriers are pure noise
#   2 candidates: one covering windows 30-70 (should show signal), one
#   covering windows 80-99 (should NOT — no carriers in that range)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

set.seed(42)

# --- Configuration -----------------------------------------------------------
TEST_DIR <- "/tmp/t05_e2e"
dir.create(TEST_DIR, recursive = TRUE, showWarnings = FALSE)

N_SAMP <- 226L
N_WIN  <- 100L
WIN_BP_WIDTH <- 10000L
sample_ids   <- sprintf("CGA%04d", seq_len(N_SAMP))
win_starts   <- (seq_len(N_WIN) - 1L) * WIN_BP_WIDTH + 1L
win_ends     <- win_starts + WIN_BP_WIDTH - 1L

# --- Build the synthetic θπ matrix -------------------------------------------
M <- matrix(rnorm(N_SAMP * N_WIN, mean = 0.005, sd = 0.001),
            nrow = N_SAMP, ncol = N_WIN,
            dimnames = list(sample_ids, NULL))
# 60 carriers: θπ jumps from ~0.005 to ~0.012 starting at window 50
for (i in 1:60) {
  M[i, 50:N_WIN] <- M[i, 50:N_WIN] + 0.007
}

# --- Pack into TR_B JSON shape -----------------------------------------------
# values is column-major: window 1 [all samples], window 2 [all samples], ...
values_flat <- as.numeric(M)   # already column-major from R matrix layout

tpw <- list(
  schema_version = 1L,
  layer          = "theta_pi_per_window",
  chrom          = "TEST_CHR",
  scale          = "per_site",
  grid_mode      = "native",
  available_modes = c("per_site"),
  default_mode    = "per_site",
  n_samples      = N_SAMP,
  n_windows      = N_WIN,
  sample_ids     = sample_ids,
  values         = round(values_flat, 6),
  windows        = lapply(seq_len(N_WIN), function(wi) list(
    idx           = wi,
    start_bp      = win_starts[wi],
    end_bp        = win_ends[wi],
    n_sites_floor = 100L
  )),
  samples        = list()  # not used by T05; can be empty
)

J <- list(
  schema_version       = 2L,
  chrom                = "TEST_CHR",
  `_layers_present`    = c("theta_pi_per_window"),
  theta_pi_per_window  = tpw
)

json_path <- file.path(TEST_DIR, "TEST_CHR_phase2_theta.json")
write(toJSON(J, auto_unbox = TRUE, digits = 8, null = "null"), json_path)
cat("[t05-test] Wrote synthetic JSON:", json_path, "\n")

# --- Build candidate TSV -----------------------------------------------------
cand_dt <- data.table(
  candidate_id = c("TEST_inv_1", "TEST_no_signal"),
  chrom        = c("TEST_CHR",   "TEST_CHR"),
  start_bp     = c(win_starts[30], win_starts[80]),
  end_bp       = c(win_ends[70],   win_ends[99])
)
cand_path <- file.path(TEST_DIR, "candidates.tsv")
fwrite(cand_dt, cand_path, sep = "\t")
cat("[t05-test] Wrote candidate TSV:", cand_path, "\n")

# --- Run the driver ----------------------------------------------------------
out_dir <- file.path(TEST_DIR, "out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

t05_script <- "/home/claude/work/phase_4_resolution/4b_theta_resolution/STEP_T05_theta_cusum.R"
lib_path   <- "/home/claude/work/phase_4_resolution/shared_lib/lib_persample_cusum.R"

cat("[t05-test] Running driver...\n")
status <- system2("Rscript",
                  args = c(t05_script,
                           "--json", json_path,
                           "--candidates", cand_path,
                           "--out-dir", out_dir,
                           "--lib", lib_path),
                  stdout = "", stderr = "")
if (status != 0) stop("[t05-test] driver exited with status ", status)

# --- Verify outputs ----------------------------------------------------------
ps_path <- file.path(out_dir, "theta_cusum_per_sample.tsv.gz")
ss_path <- file.path(out_dir, "theta_cusum_summary.tsv")
stopifnot(file.exists(ps_path), file.exists(ss_path))

ps <- fread(ps_path)
ss <- fread(ss_path)

cat("\n=== per-sample TSV (head) ===\n")
print(head(ps[informative == TRUE], 6))
cat("\n=== summary TSV ===\n")
print(ss)

# --- Assertions --------------------------------------------------------------
.fail <- 0L
.pass <- 0L
chk <- function(cond, msg) {
  if (cond) { cat("  ok :", msg, "\n"); .pass <<- .pass + 1L }
  else      { cat("  FAIL:", msg, "\n"); .fail <<- .fail + 1L }
}

cat("\n=== assertions ===\n")
# Should have 226 rows for inv_1, 226 rows for no_signal = 452 total
chk(nrow(ps) == 2L * N_SAMP,
    sprintf("per-sample TSV has 2*%d = %d rows (got %d)",
            N_SAMP, 2L * N_SAMP, nrow(ps)))

# inv_1: ~60 carriers should be informative
inv1_inf <- ps[candidate_id == "TEST_inv_1" & informative == TRUE]
chk(nrow(inv1_inf) >= 50L && nrow(inv1_inf) <= 70L,
    sprintf("TEST_inv_1: %d informative (target ~60, range 50-70)",
            nrow(inv1_inf)))

# Their cp_bp should center near win_mids[50]
expected_bp <- as.integer((win_starts[50] + win_ends[50]) / 2)
if (nrow(inv1_inf) > 0L) {
  med <- median(inv1_inf$cp_bp)
  chk(abs(med - expected_bp) < 50000L,
      sprintf("TEST_inv_1: median cp_bp = %d (target ~%d, ±50 kb)",
              med, expected_bp))
}

# no_signal candidate: should have very few informative
nosig_inf <- ps[candidate_id == "TEST_no_signal" & informative == TRUE]
chk(nrow(nosig_inf) <= 10L,
    sprintf("TEST_no_signal: %d informative (target ~0, range 0-10)",
            nrow(nosig_inf)))

# Summary should have 2 rows
chk(nrow(ss) == 2L, sprintf("summary has 2 rows (got %d)", nrow(ss)))

# Check column completeness on per-sample TSV
expected_cols <- c("candidate_id", "chrom", "stream", "sample_id",
                   "candidate_start_bp", "candidate_end_bp",
                   "n_windows_in_range", "cp_idx", "cp_bp",
                   "dist_to_left_kb", "dist_to_right_kb",
                   "cp_side_inferred", "strength", "asymmetry",
                   "left_mean", "right_mean", "n_used", "informative")
miss <- setdiff(expected_cols, names(ps))
chk(length(miss) == 0L,
    sprintf("per-sample TSV has all expected columns (missing: %s)",
            if (length(miss) > 0) paste(miss, collapse = ", ") else "none"))

# Stream tag is "theta"
chk(all(ps$stream == "theta"), "stream column is uniformly 'theta'")

cat("\n", .pass, " passed,  ", .fail, " failed\n", sep = "")
quit(status = if (.fail > 0L) 1L else 0L)
