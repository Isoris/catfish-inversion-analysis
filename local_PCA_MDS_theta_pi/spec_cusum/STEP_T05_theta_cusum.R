#!/usr/bin/env Rscript
# =============================================================================
# STEP_T05_theta_cusum.R
# =============================================================================
# Phase 4 / 4b_theta_resolution — per-sample CUSUM on θπ matrices.
#
# Reads:   <chrom>_phase2_theta.json  (TR_B output; contains theta_pi_per_window)
#          --candidates <tsv>  (optional; columns: candidate_id chrom start_bp end_bp)
#                              If absent, runs whole-chromosome mode.
#
# Calls:   persample_cusum(matrix, window_pos_bp, candidate_id) from
#          shared_lib/lib_persample_cusum.R, returning per-sample observations.
#
# Writes:  <out_dir>/theta_cusum_per_sample.tsv.gz
#            One row per (sample × candidate). Columns mirror the grain of
#            02_ancestral_fragments.R per-sample output so streams can be
#            stacked for future consensus designs.
#
#          <out_dir>/theta_cusum_summary.tsv
#            One row per candidate. n_total, n_informative, plus the empirical
#            distribution of cp_bp (median, MAD, IQR, min, max). Deliberately
#            NOT the modal-position estimate — that's a consensus job, not
#            an observation job.
#
# Run modes
# ---------
#   per_candidate    Default. Slice the per-window matrix to the windows
#                    overlapping each candidate. CUSUM finds ONE changepoint
#                    per (sample × candidate) — the strongest one within that
#                    candidate's window range.
#   whole_chrom      Run CUSUM on the full chromosome per sample. Yields
#                    one cp_bp per sample = the location of that sample's
#                    strongest signal anywhere on the chromosome. Diagnostic /
#                    exploratory; not for production breakpoint refinement.
#
# Usage
# -----
#   Rscript STEP_T05_theta_cusum.R \
#       --json   /path/to/<chr>_phase2_theta.json \
#       --candidates /path/to/candidate_intervals.tsv \
#       --out-dir /path/to/output/dir \
#       [--mode per_candidate|whole_chrom]   (default: per_candidate)
#       [--lib /path/to/lib_persample_cusum.R]
#
# Walltime: ~10s per chromosome at 16,500 windows × ~10 candidates.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# Argument parsing
# =============================================================================
.args <- commandArgs(trailingOnly = TRUE)
opts <- list(
  json        = NA_character_,
  candidates  = NA_character_,
  out_dir     = NA_character_,
  mode        = "per_candidate",
  lib         = NA_character_
)
.ai <- 1L
while (.ai <= length(.args)) {
  k <- .args[.ai]
  v <- if (.ai < length(.args)) .args[.ai + 1L] else NA_character_
  if      (k == "--json")        { opts$json <- v;        .ai <- .ai + 2L }
  else if (k == "--candidates")  { opts$candidates <- v;  .ai <- .ai + 2L }
  else if (k == "--out-dir")     { opts$out_dir <- v;     .ai <- .ai + 2L }
  else if (k == "--mode")        { opts$mode <- v;        .ai <- .ai + 2L }
  else if (k == "--lib")         { opts$lib <- v;         .ai <- .ai + 2L }
  else { stop("Unknown arg: ", k) }
}

if (is.na(opts$json) || !file.exists(opts$json)) {
  stop("[T05] --json is required and must point to a TR_B output JSON")
}
if (is.na(opts$out_dir)) {
  stop("[T05] --out-dir is required")
}
if (!opts$mode %in% c("per_candidate", "whole_chrom")) {
  stop("[T05] --mode must be 'per_candidate' or 'whole_chrom', got: ", opts$mode)
}
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

# Locate the lib. Default: same directory as this script's `../shared_lib/`.
if (is.na(opts$lib)) {
  this_script <- normalizePath(sub("--file=", "",
                                   grep("^--file=", commandArgs(trailingOnly = FALSE),
                                        value = TRUE)[1]),
                               mustWork = FALSE)
  if (is.na(this_script) || !nzchar(this_script)) this_script <- "STEP_T05_theta_cusum.R"
  candidates_path <- c(
    file.path(dirname(this_script), "..", "shared_lib", "lib_persample_cusum.R"),
    file.path(dirname(this_script), "shared_lib", "lib_persample_cusum.R"),
    "phase_4_resolution/shared_lib/lib_persample_cusum.R",
    "shared_lib/lib_persample_cusum.R"
  )
  hit <- candidates_path[file.exists(candidates_path)]
  if (length(hit) == 0L) {
    stop("[T05] Could not locate lib_persample_cusum.R — pass --lib explicitly")
  }
  opts$lib <- hit[1]
}
source(opts$lib)
message("[T05] Library loaded from ", opts$lib)

# =============================================================================
# Read TR_B JSON and reconstruct the per-window matrix
# =============================================================================
message("[T05] Reading θπ JSON: ", opts$json)
J <- fromJSON(opts$json, simplifyVector = TRUE, simplifyDataFrame = FALSE)

tpw <- J$theta_pi_per_window
if (is.null(tpw)) stop("[T05] JSON has no theta_pi_per_window layer")

chrom    <- tpw$chrom %||% J$chrom %||% NA_character_
n_samp   <- as.integer(tpw$n_samples)
n_win    <- as.integer(tpw$n_windows)
sids_raw <- tpw$sample_ids
if (is.list(sids_raw)) sids_raw <- unlist(sids_raw, use.names = FALSE)
sample_ids <- as.character(sids_raw)
if (length(sample_ids) != n_samp) {
  stop(sprintf("[T05] sample_ids length %d != n_samples %d",
               length(sample_ids), n_samp))
}

# `values` is row-major over windows: window 1 [all samples], window 2 [...], ...
# Reconstruct a sample × window matrix: M[s, w] = values[(w-1)*n_samp + s]
vals <- as.numeric(unlist(tpw$values, use.names = FALSE))
if (length(vals) != n_samp * n_win) {
  stop(sprintf("[T05] values length %d != n_samp * n_win = %d * %d = %d",
               length(vals), n_samp, n_win, n_samp * n_win))
}
M <- matrix(vals, nrow = n_samp, ncol = n_win, byrow = FALSE,
            dimnames = list(sample_ids, NULL))
# byrow = FALSE means values fills the matrix column-major, which matches
# the JSON's "all samples for window 1, then window 2, ..." layout.

# Window coordinates: take the midpoint as the canonical position.
win_starts <- vapply(tpw$windows, function(w) as.integer(w$start_bp), integer(1))
win_ends   <- vapply(tpw$windows, function(w) as.integer(w$end_bp),   integer(1))
win_mids   <- as.integer((win_starts + win_ends) %/% 2L)

message(sprintf("[T05] Loaded matrix: %d samples × %d windows on %s",
                n_samp, n_win, chrom))
message(sprintf("[T05] Window range: %d - %d bp",
                min(win_starts), max(win_ends)))

# =============================================================================
# Determine the set of intervals to process
# =============================================================================
intervals <- NULL
if (opts$mode == "whole_chrom") {
  intervals <- data.table(
    candidate_id = paste0(chrom, "_whole"),
    chrom        = chrom,
    start_bp     = min(win_starts),
    end_bp       = max(win_ends)
  )
  message("[T05] Mode: whole_chrom — 1 interval covering the whole chromosome")
} else {
  if (is.na(opts$candidates) || !file.exists(opts$candidates)) {
    stop("[T05] --candidates is required in per_candidate mode")
  }
  intervals <- fread(opts$candidates)
  required_cols <- c("candidate_id", "chrom", "start_bp", "end_bp")
  miss <- setdiff(required_cols, names(intervals))
  if (length(miss) > 0L) {
    stop("[T05] Candidate TSV missing columns: ", paste(miss, collapse = ", "))
  }
  # Filter to this chromosome. Rename local var to avoid `chrom == chrom`
  # ambiguity inside data.table's [.
  this_chrom <- chrom
  intervals <- intervals[chrom == this_chrom]
  if (nrow(intervals) == 0L) {
    message("[T05] No candidates on chromosome ", this_chrom, " — exiting cleanly")
    quit(status = 0)
  }
  message(sprintf("[T05] Mode: per_candidate — %d candidate(s) on %s",
                  nrow(intervals), this_chrom))
}

# =============================================================================
# Run CUSUM per interval
# =============================================================================
all_persample <- vector("list", nrow(intervals))
all_summary   <- vector("list", nrow(intervals))

for (i in seq_len(nrow(intervals))) {
  cid     <- as.character(intervals$candidate_id[i])
  c_start <- as.integer(intervals$start_bp[i])
  c_end   <- as.integer(intervals$end_bp[i])

  # Window selection: any window overlapping the interval
  in_range <- which(win_ends >= c_start & win_starts <= c_end)
  if (length(in_range) < 5L) {
    message(sprintf("[T05] cid=%s SKIP — only %d overlapping windows (need >=5)",
                    cid, length(in_range)))
    next
  }

  M_sub        <- M[, in_range, drop = FALSE]
  win_mids_sub <- win_mids[in_range]

  ps <- persample_cusum(M_sub, win_mids_sub, candidate_id = cid)

  # Decorate: chrom, candidate boundaries, side-of-candidate inference.
  # cp_side_inferred: which candidate edge is the changepoint closer to?
  ps[, chrom := chrom]
  ps[, candidate_start_bp := c_start]
  ps[, candidate_end_bp   := c_end]
  ps[, dist_to_left_kb  := round(abs(cp_bp - c_start) / 1000, 2)]
  ps[, dist_to_right_kb := round(abs(cp_bp - c_end)   / 1000, 2)]
  ps[, cp_side_inferred := fifelse(
    is.na(cp_bp), NA_character_,
    fifelse(abs(cp_bp - c_start) <= abs(cp_bp - c_end), "left", "right")
  )]
  ps[, stream := "theta"]
  ps[, n_windows_in_range := length(in_range)]

  setcolorder(ps, c("candidate_id", "chrom", "stream", "sample_id",
                    "candidate_start_bp", "candidate_end_bp",
                    "n_windows_in_range",
                    "cp_idx", "cp_bp",
                    "dist_to_left_kb", "dist_to_right_kb",
                    "cp_side_inferred",
                    "strength", "asymmetry",
                    "left_mean", "right_mean",
                    "n_used", "informative"))
  all_persample[[i]] <- ps

  # Per-candidate summary: empirical distribution shape only (no parametric fits).
  inf <- ps[informative == TRUE & is.finite(cp_bp)]
  s_dt <- data.table(
    candidate_id       = cid,
    chrom              = chrom,
    stream             = "theta",
    candidate_start_bp = c_start,
    candidate_end_bp   = c_end,
    n_windows_in_range = length(in_range),
    n_total_samples    = nrow(ps),
    n_informative      = nrow(inf),
    cp_min_bp          = if (nrow(inf) > 0) as.integer(min(inf$cp_bp)) else NA_integer_,
    cp_q25_bp          = if (nrow(inf) > 0) as.integer(round(quantile(inf$cp_bp, 0.25))) else NA_integer_,
    cp_median_bp       = if (nrow(inf) > 0) as.integer(round(median(inf$cp_bp))) else NA_integer_,
    cp_q75_bp          = if (nrow(inf) > 0) as.integer(round(quantile(inf$cp_bp, 0.75))) else NA_integer_,
    cp_max_bp          = if (nrow(inf) > 0) as.integer(max(inf$cp_bp)) else NA_integer_,
    cp_iqr_kb          = if (nrow(inf) > 0)
      round(as.numeric(IQR(inf$cp_bp)) / 1000, 2) else NA_real_,
    cp_mad_kb          = if (nrow(inf) > 0)
      round(as.numeric(mad(inf$cp_bp)) / 1000, 2) else NA_real_,
    n_left             = if (nrow(inf) > 0) sum(inf$cp_side_inferred == "left",  na.rm = TRUE) else 0L,
    n_right            = if (nrow(inf) > 0) sum(inf$cp_side_inferred == "right", na.rm = TRUE) else 0L,
    n_asym_pos         = if (nrow(inf) > 0) sum(inf$asymmetry ==  1L, na.rm = TRUE) else 0L,
    n_asym_neg         = if (nrow(inf) > 0) sum(inf$asymmetry == -1L, na.rm = TRUE) else 0L
  )
  all_summary[[i]] <- s_dt

  message(sprintf("[T05] cid=%s: %d windows, %d/%d informative, cp_median=%s, cp_IQR=%s kb",
                  cid, length(in_range),
                  nrow(inf), nrow(ps),
                  if (nrow(inf) > 0) format(s_dt$cp_median_bp, big.mark = ",") else "NA",
                  if (nrow(inf) > 0) sprintf("%.1f", s_dt$cp_iqr_kb) else "NA"))
}

# =============================================================================
# Write outputs
# =============================================================================
all_persample <- Filter(Negate(is.null), all_persample)
all_summary   <- Filter(Negate(is.null), all_summary)

if (length(all_persample) == 0L) {
  message("[T05] No candidates produced output. Exiting.")
  quit(status = 0)
}

ps_out <- rbindlist(all_persample, fill = TRUE)
ss_out <- rbindlist(all_summary,   fill = TRUE)

ps_path <- file.path(opts$out_dir, "theta_cusum_per_sample.tsv.gz")
ss_path <- file.path(opts$out_dir, "theta_cusum_summary.tsv")
fwrite(ps_out, ps_path, sep = "\t", compress = "gzip")
fwrite(ss_out, ss_path, sep = "\t")

message("[T05] DONE. Wrote:")
message("  - ", ps_path, " (", nrow(ps_out), " rows)")
message("  - ", ss_path, " (", nrow(ss_out), " rows)")
