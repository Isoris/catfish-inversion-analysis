#!/usr/bin/env Rscript

# =============================================================================
# STEP_ZO_F_mds_merge.R
#
# Coordination step that runs ONCE after all per-chrom ZO_E MDS tasks
# complete. Stitches per-focal-chrom MDS RDS files into the unified
# <outprefix>.mds.rds with the $per_chr structure ZO_G iterates over.
#
# Pipeline position:
#   ZO_E (per-chrom MDS)  ->  ZO_F (this script: merge)  ->  ZO_G (precomp)
#
# Inputs:
#   --tmpdir       <dir>   directory holding per-focal-chr .mds_perchr.rds
#   --outprefix    <path>  full path stem for outputs
#   [--mds_mode    chromosome]   default; "chunked_2x" still supported
#
# Outputs (at <outprefix>.*):
#   .mds.rds                            list($dt, $per_chr, $mds_mode)
#                                       — ZO_G reads $per_chr from here
#   .window_mds.tsv.gz                  per-window MDS coords (flat TSV)
#   .mds_mode_metadata.tsv              run config
#   .mds_background_<chr>.txt           per-focal-chr background indices
#                                       (only under legacy chunked_2x mode)
#
# Removed 2026-05-13:
#   - cluster_outliers_bp gap-bp candidate-region merge
#   - .candidate_regions.tsv.gz, .candidate_window_membership.tsv.gz outputs
#   Nothing downstream consumed them; the lostruct-style nominator over-
#   merged in dense-outlier regimes. L1/L2 read sim_mat directly.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# PARSE ARGS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
tmpdir     <- NULL
outprefix  <- NULL
MDS_MODE   <- "chromosome"

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--tmpdir" && i < length(args)) {
    tmpdir <- args[i + 1]; i <- i + 2L
  } else if (a == "--outprefix" && i < length(args)) {
    outprefix <- args[i + 1]; i <- i + 2L
  } else if (a == "--mds_mode" && i < length(args)) {
    MDS_MODE <- args[i + 1]; i <- i + 2L
  } else if (a == "--gap_bp" || a == "--min_windows") {
    # Legacy flags from the removed candidate-region clustering — silently
    # consume the value for CLI compatibility but ignore it.
    i <- i + 2L
  } else {
    i <- i + 1L
  }
}

if (is.null(tmpdir) || is.null(outprefix)) {
  stop("Usage: Rscript STEP10v2_stage2_merge_mds.R --tmpdir <dir> --outprefix <prefix> ...")
}

message("[STEP10v2-S2] ═══════ MDS merge ═══════")
message("[STEP10v2-S2] Reading from: ", tmpdir)

# =============================================================================
# LOAD PER-CHR RESULTS
# =============================================================================

rds_files <- sort(list.files(tmpdir, pattern = "\\.mds_perchr\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .mds_perchr.rds files found in: ", tmpdir)

all_chr_results <- list()
metadata_rows <- list()

for (f in rds_files) {
  obj <- readRDS(f)

  # Skip sentinel entries
  if (!is.null(obj$skip) && obj$skip) {
    message("[ZO_F] SKIP ", obj$chrom, ": ", obj$reason)
    next
  }

  chr <- obj$out_dt$chrom[1]
  all_chr_results[[chr]] <- list(
    out_dt = obj$out_dt,
    dmat   = obj$dmat,
    mds    = obj$mds
  )

  # Read metadata
  metaf <- file.path(tmpdir, paste0(chr, ".metadata.tsv"))
  if (file.exists(metaf)) {
    metadata_rows[[chr]] <- fread(metaf)
  }

  # Copy background IDs to final location (only present under legacy chunked_2x
  # mode; absent under the default chromosome mode).
  bgf <- file.path(tmpdir, paste0(chr, ".background_ids.txt"))
  if (file.exists(bgf)) {
    final_bgf <- paste0(outprefix, ".mds_background_", chr, ".txt")
    file.copy(bgf, final_bgf, overwrite = TRUE)
  }

  message("[ZO_F] ", chr, ": ", nrow(obj$out_dt), " windows")
}

if (length(all_chr_results) == 0) stop("No chromosomes produced results")
message("[ZO_F] Loaded ", length(all_chr_results), " chromosomes")

# =============================================================================
# MERGE (per-chrom MDS results -> unified $per_chr container)
# =============================================================================
# This script only stitches the 28 per-chrom .mds_perchr.rds files into one
# unified .mds.rds container (the $per_chr keyed list ZO_G iterates over)
# and assigns global_window_id. The lostruct-style gap-bp candidate-region
# clustering that used to live here has been removed (2026-05-13). It was
# never consumed by L1/L2 stripe detection (which reads sim_mat directly).
# =============================================================================

out_dt <- rbindlist(lapply(all_chr_results, function(x) x$out_dt), fill = TRUE)
message("[ZO_F] Merged: ", nrow(out_dt), " windows")

# ── Write outputs ────────────────────────────────────────────────────────────

meta_dt <- if (length(metadata_rows) > 0) rbindlist(metadata_rows, fill = TRUE) else data.table()

f1 <- paste0(outprefix, ".window_mds.tsv.gz")
f3 <- paste0(outprefix, ".mds.rds")
f4 <- paste0(outprefix, ".mds_mode_metadata.tsv")

dir.create(dirname(outprefix), recursive = TRUE, showWarnings = FALSE)

fwrite(out_dt, f1, sep = "\t")

# .mds.rds is the load-bearing artifact: ZO_G reads $per_chr to iterate.
saveRDS(list(
  dt = out_dt,
  per_chr = all_chr_results,
  mds_mode = MDS_MODE
), f3)

fwrite(meta_dt, f4, sep = "\t")

message("\n[DONE] ZO_F merge complete (mode=", MDS_MODE, ")")
message("  ", f1)
message("  ", f3, "  <- ZO_G reads this")
message("  ", f4)
