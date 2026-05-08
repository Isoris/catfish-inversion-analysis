#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_D_detect_L2.R   (drafts/v5)
# =============================================================================
# L2 envelope detection — refined cohort-level scan inside each L1.
# Mirrors local_PCA_MDS_z/06_detect_L2 in the z-blocks pipeline.
#
# Logic (z-blocks parity)
# -----------------------
# Within each L1's window range, identify L2s as contiguous runs of windows
# whose morphology indicates a coherent flat or spiky inversion-like block:
#
#   L2 window must satisfy:
#     max_abs_z          >= Z_L2 (default 2.5)
#     AND (flat_inv_score >= FLAT_FLOOR OR spiky_inv_score >= SPIKY_FLOOR)
#     AND fragmentation_score < FRAG_CEILING
#
# Then merge with MERGE_GAP and require min length MIN_L2_WIN. L1s with no
# qualifying L2 windows are filtered (a wide L1 with no peak is just
# elevated chromosomal jitter).
#
# Reads:   <OUTROOT>/03_per_chrom/<chr>/precomp.rds
#          <OUTROOT>/03_per_chrom/<chr>/L1_envelopes.tsv
# Writes:  <OUTROOT>/03_per_chrom/<chr>/L2_envelopes.tsv
#          <OUTROOT>/02_mds/candidate_regions.tsv.gz   (genome-wide L2 rollup)
# Usage:   Rscript STEP_TR_D_detect_L2.R --chrom <CHR>
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
CHROM <- NULL
Z_L2          <- 2.5
FLAT_FLOOR    <- 0.30
SPIKY_FLOOR   <- 0.30
FRAG_CEILING  <- 0.50
MIN_L2_WIN    <- 5L
L2_MERGE_GAP  <- 3L
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--chrom"          && i < length(args)) { CHROM        <- args[i + 1]; i <- i + 2L }
  else if (a == "--z-l2"           && i < length(args)) { Z_L2         <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--flat-floor"     && i < length(args)) { FLAT_FLOOR   <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--spiky-floor"    && i < length(args)) { SPIKY_FLOOR  <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--frag-ceiling"   && i < length(args)) { FRAG_CEILING <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--min-l2-windows" && i < length(args)) { MIN_L2_WIN   <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--merge-gap"      && i < length(args)) { L2_MERGE_GAP <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}

OUTROOT <- Sys.getenv("OUTROOT", unset = NA)
stopifnot(!is.na(OUTROOT))
per_chrom_dir <- file.path(OUTROOT, "03_per_chrom")
mds_dir       <- file.path(OUTROOT, "02_mds")
dir.create(mds_dir, recursive = TRUE, showWarnings = FALSE)
chroms <- if (!is.null(CHROM)) CHROM else list.files(per_chrom_dir, pattern = "^C_gar_LG[0-9]+$")

# Run-detection within a fixed window range, against a precomputed flag vec.
run_within <- function(flag, lo, hi, min_run, merge_gap) {
  if (lo > hi || lo < 1L || hi > length(flag)) return(data.table())
  fl <- flag[lo:hi]
  if (!any(fl, na.rm = TRUE)) return(data.table())
  runs <- rle(fl); ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1) + 1L)
  hits <- which(runs$values)
  envs <- data.table(win_start = starts[hits] + lo - 1L,
                     win_end   = ends[hits]   + lo - 1L,
                     n_windows = runs$lengths[hits])
  envs <- envs[n_windows >= min_run]
  if (nrow(envs) == 0L) return(envs)
  setorder(envs, win_start)
  m <- envs[1]
  for (k in seq.int(2L, nrow(envs))) {
    if (envs$win_start[k] - m$win_end[nrow(m)] <= merge_gap) {
      m[nrow(m), win_end   := envs$win_end[k]]
      m[nrow(m), n_windows := win_end - win_start + 1L]
    } else m <- rbind(m, envs[k])
  }
  m
}

genome_wide_l2 <- list()

for (chrom in chroms) {
  rds <- file.path(per_chrom_dir, chrom, "precomp.rds")
  l1f <- file.path(per_chrom_dir, chrom, "L1_envelopes.tsv")
  if (!file.exists(rds) || !file.exists(l1f)) { message("[TR_D] ", chrom, ": missing — skip"); next }
  precomp <- readRDS(rds); dt <- precomp$dt
  l1_dt <- fread(l1f)
  if (nrow(l1_dt) == 0L) {
    fwrite(data.table(), file.path(per_chrom_dir, chrom, "L2_envelopes.tsv"), sep = "\t")
    next
  }
  z_vec <- dt$max_abs_z
  has_morph <- "flat_inv_score" %in% names(dt)

  # Build the L2-eligibility flag vector once.
  if (has_morph) {
    flag <- (is.finite(z_vec) & z_vec >= Z_L2) &
            ((is.finite(dt$flat_inv_score)  & dt$flat_inv_score  >= FLAT_FLOOR) |
             (is.finite(dt$spiky_inv_score) & dt$spiky_inv_score >= SPIKY_FLOOR)) &
            (!is.finite(dt$fragmentation_score) | dt$fragmentation_score < FRAG_CEILING)
    mode <- "morphology"
  } else {
    flag <- is.finite(z_vec) & z_vec >= Z_L2
    mode <- "z_runs_fallback"
  }

  l2_rows <- list()
  for (k in seq_len(nrow(l1_dt))) {
    inner <- run_within(flag, l1_dt$win_start[k], l1_dt$win_end[k], MIN_L2_WIN, L2_MERGE_GAP)
    if (nrow(inner) == 0L) next
    inner[, l1_id := l1_dt$l1_id[k]]
    l2_rows[[length(l2_rows) + 1L]] <- inner
  }
  l2 <- if (length(l2_rows) > 0) rbindlist(l2_rows) else data.table()
  if (nrow(l2) > 0) {
    l2[, `:=`(
      chrom        = chrom,
      start_bp     = dt$start_bp[win_start],
      end_bp       = dt$end_bp[win_end],
      span_kb      = round((dt$end_bp[win_end] - dt$start_bp[win_start]) / 1000, 1),
      peak_z       = vapply(seq_len(.N), function(j) max(z_vec[l2$win_start[j]:l2$win_end[j]],
                                                          na.rm = TRUE), numeric(1)),
      mean_z       = vapply(seq_len(.N), function(j) mean(z_vec[l2$win_start[j]:l2$win_end[j]],
                                                           na.rm = TRUE), numeric(1)),
      mean_flat    = if (has_morph) vapply(seq_len(.N), function(j)
        mean(dt$flat_inv_score[l2$win_start[j]:l2$win_end[j]], na.rm = TRUE), numeric(1)) else NA_real_,
      mean_spiky   = if (has_morph) vapply(seq_len(.N), function(j)
        mean(dt$spiky_inv_score[l2$win_start[j]:l2$win_end[j]], na.rm = TRUE), numeric(1)) else NA_real_,
      mean_frag    = if (has_morph) vapply(seq_len(.N), function(j)
        mean(dt$fragmentation_score[l2$win_start[j]:l2$win_end[j]], na.rm = TRUE), numeric(1)) else NA_real_,
      l2_id        = paste0(chrom, "_L2_", sprintf("%03d", seq_len(.N))),
      candidate_id = paste0(chrom, "_C_", sprintf("%03d", seq_len(.N))),
      win_start_idx0 = win_start - 1L,
      win_end_idx0   = win_end   - 1L,
      detection_mode = mode
    )]
    setcolorder(l2, c("candidate_id", "l2_id", "l1_id", "chrom",
                      "win_start", "win_end", "win_start_idx0", "win_end_idx0",
                      "n_windows", "start_bp", "end_bp", "span_kb",
                      "peak_z", "mean_z", "mean_flat", "mean_spiky", "mean_frag",
                      "detection_mode"))
    genome_wide_l2[[chrom]] <- l2
  }
  fwrite(l2, file.path(per_chrom_dir, chrom, "L2_envelopes.tsv"), sep = "\t")
  message(sprintf("[TR_D] %s: %d L2 envelopes inside %d L1s (mode=%s)",
                  chrom, nrow(l2), nrow(l1_dt), mode))
}

if (length(genome_wide_l2) > 0) {
  fwrite(rbindlist(genome_wide_l2, fill = TRUE),
         file.path(mds_dir, "candidate_regions.tsv.gz"), sep = "\t", compress = "gzip")
  message("[TR_D] candidate_regions.tsv.gz: ",
          sum(vapply(genome_wide_l2, nrow, integer(1))), " L2s across ",
          length(genome_wide_l2), " chroms")
}

message("[TR_D] DONE")