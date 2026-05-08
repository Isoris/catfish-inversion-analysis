#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_C_detect_L1.R   (drafts/v5)
# =============================================================================
# L1 envelope detection — wide-net cohort-level scan. Mirrors the role of
# local_PCA_MDS_z/04_detect_L1: identifies broad regions of elevated
# inversion-like signal that L2 will then refine.
#
# Logic (z-blocks parity)
# -----------------------
# Seed-and-grow against the rich precomp$dt produced by TR_B v5:
#
#   1. SEED windows = those flagged adaptive_seed (Beta-fit p<0.01 on
#      inv_likeness) OR with inv_likeness >= INV_LIKE_HARD_FLOOR (default
#      0.50) OR max_abs_z >= Z_HARD_FLOOR (default 2.5).
#   2. GROW each seed left/right while neighbours satisfy a lenient
#      condition: inv_likeness >= INV_LIKE_GROW_FLOOR (default 0.30)
#      AND morphology-coherent (flat_inv_score + spiky_inv_score - 0.5 *
#      fragmentation_score > 0). This prevents L1s from absorbing
#      incoherent jitter.
#   3. MERGE adjacent runs separated by <= MERGE_GAP windows (default 5).
#   4. FILTER runs shorter than MIN_L1_WIN (default 10).
#
# Outputs the same TSV shape as the simpler version (the simpler version
# is preserved as the fallback path when TR_B was run without morphology
# enrichment — script auto-detects).
#
# Reads:   <OUTROOT>/03_per_chrom/<chr>/precomp.rds
# Writes:  <OUTROOT>/03_per_chrom/<chr>/L1_envelopes.tsv
# Usage:   Rscript STEP_TR_C_detect_L1.R --chrom <CHR>
#                                        [--inv-like-seed 0.50]
#                                        [--inv-like-grow 0.30]
#                                        [--z-hard-floor 2.5]
#                                        [--min-l1-windows 10]
#                                        [--merge-gap 5]
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
CHROM                <- NULL
INV_LIKE_HARD_FLOOR  <- 0.50
INV_LIKE_GROW_FLOOR  <- 0.30
Z_HARD_FLOOR         <- 2.5
MIN_L1_WIN           <- 10L
MERGE_GAP            <- 5L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--chrom"          && i < length(args)) { CHROM               <- args[i + 1]; i <- i + 2L }
  else if (a == "--inv-like-seed"  && i < length(args)) { INV_LIKE_HARD_FLOOR <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--inv-like-grow"  && i < length(args)) { INV_LIKE_GROW_FLOOR <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--z-hard-floor"   && i < length(args)) { Z_HARD_FLOOR        <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--min-l1-windows" && i < length(args)) { MIN_L1_WIN          <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--merge-gap"      && i < length(args)) { MERGE_GAP           <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}

OUTROOT <- Sys.getenv("OUTROOT", unset = NA)
stopifnot(!is.na(OUTROOT))
per_chrom_dir <- file.path(OUTROOT, "03_per_chrom")
chroms <- if (!is.null(CHROM)) CHROM else list.files(per_chrom_dir, pattern = "^C_gar_LG[0-9]+$")

# Seed-grow: returns 1-indexed (start, end) pairs.
seed_grow <- function(dt, seed_flag, grow_flag, merge_gap, min_run) {
  n <- nrow(dt)
  if (!any(seed_flag, na.rm = TRUE)) return(data.table())
  # Find seed windows; for each, expand left/right while grow_flag holds.
  starts <- ends <- integer(0)
  visited <- logical(n)
  for (i in which(seed_flag)) {
    if (visited[i]) next
    lo <- i; hi <- i
    while (lo > 1L && grow_flag[lo - 1L]) lo <- lo - 1L
    while (hi < n  && grow_flag[hi + 1L]) hi <- hi + 1L
    visited[lo:hi] <- TRUE
    starts <- c(starts, lo); ends <- c(ends, hi)
  }
  if (length(starts) == 0L) return(data.table())
  envs <- data.table(win_start = starts, win_end = ends)
  setorder(envs, win_start)
  # Merge overlapping / close envelopes.
  merged <- envs[1]
  for (k in seq.int(2L, nrow(envs))) {
    if (envs$win_start[k] - merged$win_end[nrow(merged)] <= merge_gap) {
      merged[nrow(merged), win_end := pmax(win_end, envs$win_end[k])]
    } else {
      merged <- rbind(merged, envs[k])
    }
  }
  merged[, n_windows := win_end - win_start + 1L]
  merged[n_windows >= min_run]
}

# Fallback: pure |Z|-runs (used if precomp lacks the morphology enrichment).
fallback_z_runs <- function(z_vec, threshold, min_run, merge_gap) {
  flag <- !is.na(z_vec) & z_vec > threshold
  if (!any(flag)) return(data.table())
  runs <- rle(flag); ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1) + 1L)
  hits <- which(runs$values)
  envs <- data.table(win_start = starts[hits], win_end = ends[hits],
                     n_windows = runs$lengths[hits])
  envs <- envs[n_windows >= min_run]
  if (nrow(envs) == 0L) return(envs)
  setorder(envs, win_start)
  m <- envs[1]
  for (k in seq.int(2L, nrow(envs))) {
    if (envs$win_start[k] - m$win_end[nrow(m)] <= merge_gap) {
      m[nrow(m), win_end := envs$win_end[k]]
      m[nrow(m), n_windows := win_end - win_start + 1L]
    } else m <- rbind(m, envs[k])
  }
  m
}

for (chrom in chroms) {
  rds <- file.path(per_chrom_dir, chrom, "precomp.rds")
  if (!file.exists(rds)) { message("[TR_C] ", chrom, ": no precomp — skip"); next }
  precomp <- readRDS(rds)
  dt <- precomp$dt; n_win <- nrow(dt)
  z_vec <- dt$max_abs_z

  has_morph <- "inv_likeness" %in% names(dt)
  if (has_morph) {
    seed <- (!is.na(dt$adaptive_seed) & dt$adaptive_seed) |
            (is.finite(dt$inv_likeness) & dt$inv_likeness >= INV_LIKE_HARD_FLOOR) |
            (is.finite(z_vec) & z_vec >= Z_HARD_FLOOR)
    coh <- (dt$flat_inv_score %||% 0) + (dt$spiky_inv_score %||% 0) -
           0.5 * (dt$fragmentation_score %||% 0)
    grow <- (is.finite(dt$inv_likeness) & dt$inv_likeness >= INV_LIKE_GROW_FLOOR) &
            (is.finite(coh) & coh > 0)
    envs <- seed_grow(dt, seed, grow, MERGE_GAP, MIN_L1_WIN)
    mode <- "morphology"
  } else {
    envs <- fallback_z_runs(z_vec, Z_HARD_FLOOR, MIN_L1_WIN, MERGE_GAP)
    mode <- "z_runs_fallback"
  }

  if (nrow(envs) > 0) {
    envs[, `:=`(
      chrom     = chrom,
      start_bp  = dt$start_bp[win_start],
      end_bp    = dt$end_bp[win_end],
      span_kb   = round((dt$end_bp[win_end] - dt$start_bp[win_start]) / 1000, 1),
      peak_z    = vapply(seq_len(.N), function(k) max(z_vec[envs$win_start[k]:envs$win_end[k]],
                                                       na.rm = TRUE), numeric(1)),
      mean_z    = vapply(seq_len(.N), function(k) mean(z_vec[envs$win_start[k]:envs$win_end[k]],
                                                        na.rm = TRUE), numeric(1)),
      mean_inv_like = if (has_morph)
        vapply(seq_len(.N), function(k) mean(dt$inv_likeness[envs$win_start[k]:envs$win_end[k]],
                                              na.rm = TRUE), numeric(1)) else NA_real_,
      l1_id     = paste0(chrom, "_L1_", sprintf("%03d", seq_len(.N))),
      win_start_idx0 = win_start - 1L,
      win_end_idx0   = win_end   - 1L,
      detection_mode = mode
    )]
    setcolorder(envs, c("l1_id", "chrom", "win_start", "win_end",
                        "win_start_idx0", "win_end_idx0", "n_windows",
                        "start_bp", "end_bp", "span_kb",
                        "peak_z", "mean_z", "mean_inv_like", "detection_mode"))
  }
  fwrite(envs, file.path(per_chrom_dir, chrom, "L1_envelopes.tsv"), sep = "\t")
  message(sprintf("[TR_C] %s: %d L1 envelopes (mode=%s)", chrom, nrow(envs), mode))
}

message("[TR_C] DONE")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a