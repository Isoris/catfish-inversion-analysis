#!/usr/bin/env Rscript

# =============================================================================
# 02b_mds_merge.R
#
# STAGE 2 (merge) — coordination step that runs ONCE after all 02a array
# tasks complete. Stitches per-focal-chrom MDS RDS files into the unified
# <outprefix>.mds.rds with the $per_chr structure every downstream consumer
# (step 03 precompute, etc.) expects. Also runs candidate-region detection.
#
# Pipeline position:
#   02a mds_compute  -> 02b MDS_MERGE  -> 03 precompute  -> ...
#
# ── Inputs ────────────────────────────────────────────────────────────────
#   --tmpdir       <dir>   directory holding per-focal-chr .mds_perchr.rds
#   --outprefix    <path>  full path stem for outputs (e.g. .../inversion_localpca)
#   [--mds_mode    chunked_2x]
#   [--gap_bp      500000] gap-merge tolerance for candidate regions
#   [--min_windows 3]      min windows per candidate region
#
# ── Outputs (at <outprefix>.*) ────────────────────────────────────────────
#   .mds.rds                            main object ($dt, $per_chr,
#                                       $candidate_regions, $mds_mode)
#   .window_mds.tsv.gz                  per-window MDS coords
#   .candidate_regions.tsv.gz           coarse z-outlier regions
#   .candidate_window_membership.tsv.gz region <-> window mapping
#   .mds_mode_metadata.tsv              run config
#   .mds_background_<chr>.txt           per-focal-chr background indices
#                                       (chunked modes only)
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: STEP10v2_stage2_merge_mds.R v8.3-parallel)
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
MDS_MODE   <- "chunked_2x"
GAP_BP     <- 500000
MIN_WINDOWS <- 3L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--tmpdir" && i < length(args)) {
    tmpdir <- args[i + 1]; i <- i + 2L
  } else if (a == "--outprefix" && i < length(args)) {
    outprefix <- args[i + 1]; i <- i + 2L
  } else if (a == "--mds_mode" && i < length(args)) {
    MDS_MODE <- args[i + 1]; i <- i + 2L
  } else if (a == "--gap_bp" && i < length(args)) {
    GAP_BP <- as.numeric(args[i + 1]); i <- i + 2L
  } else if (a == "--min_windows" && i < length(args)) {
    MIN_WINDOWS <- as.integer(args[i + 1]); i <- i + 2L
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
# CLUSTER OUTLIERS (identical to monolithic STEP10v2)
# =============================================================================

cluster_outliers_bp <- function(dt, flag_col, gap_bp, min_windows) {
  idx <- which(dt[[flag_col]])
  if (length(idx) == 0) return(NULL)
  clusters <- list(); cur <- idx[1]
  if (length(idx) > 1) {
    for (ii in idx[-1]) {
      if ((dt$start_bp[ii] - dt$end_bp[max(cur)]) <= gap_bp) {
        cur <- c(cur, ii)
      } else {
        if (length(cur) >= min_windows) clusters[[length(clusters) + 1]] <- cur
        cur <- ii
      }
    }
  }
  if (length(cur) >= min_windows) clusters[[length(clusters) + 1]] <- cur
  clusters
}

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
    message("[STEP10v2-S2] SKIP ", obj$chrom, ": ", obj$reason)
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

  # Copy background IDs to final location
  bgf <- file.path(tmpdir, paste0(chr, ".background_ids.txt"))
  if (file.exists(bgf)) {
    final_bgf <- paste0(outprefix, ".mds_background_", chr, ".txt")
    file.copy(bgf, final_bgf, overwrite = TRUE)
  }

  n_outlier <- sum(obj$out_dt$MDS1_outlier, na.rm = TRUE)
  message("[STEP10v2-S2] ", chr, ": ", nrow(obj$out_dt), " windows, ",
          n_outlier, " MDS1 outliers")
}

if (length(all_chr_results) == 0) stop("No chromosomes produced results")
message("[STEP10v2-S2] Loaded ", length(all_chr_results), " chromosomes")

# =============================================================================
# MERGE + CANDIDATE REGIONS (identical logic to monolithic STEP10v2)
# =============================================================================

out_dt <- rbindlist(lapply(all_chr_results, function(x) x$out_dt), fill = TRUE)
message("[STEP10v2-S2] Merged: ", nrow(out_dt), " windows")

mds_cols <- grep("^MDS\\d+$", names(out_dt), value = TRUE)
cand_list <- list()
membership_list <- list()
cand_id <- 0L

for (ax in seq_along(mds_cols)) {
  flag_col <- paste0("MDS", ax, "_outlier")
  if (!(flag_col %in% names(out_dt))) next

  for (chr in unique(out_dt$chrom)) {
    sub <- out_dt[chrom == chr][order(start_bp)]
    clusters <- cluster_outliers_bp(sub, flag_col, GAP_BP, MIN_WINDOWS)
    if (!is.null(clusters)) {
      for (cl in clusters) {
        cand_id <- cand_id + 1L
        xx <- sub[cl]
        cand_list[[length(cand_list) + 1]] <- data.table(
          candidate_id = cand_id, mds_axis = ax, chrom = chr,
          start_bp = min(xx$start_bp), end_bp = max(xx$end_bp),
          center_bp = as.integer((min(xx$start_bp) + max(xx$end_bp)) / 2),
          n_windows = nrow(xx),
          first_global_window_id = min(xx$global_window_id),
          last_global_window_id = max(xx$global_window_id)
        )
        for (k in seq_len(nrow(xx))) {
          membership_list[[length(membership_list) + 1]] <- data.table(
            candidate_id = cand_id,
            global_window_id = xx$global_window_id[k],
            chrom = chr,
            start_bp = xx$start_bp[k],
            end_bp = xx$end_bp[k],
            mds_axis = ax
          )
        }
      }
    }
  }
}

cand_dt <- if (length(cand_list) > 0) rbindlist(cand_list) else {
  data.table(candidate_id = integer(), chrom = character())
}
membership_dt <- if (length(membership_list) > 0) rbindlist(membership_list) else {
  data.table(candidate_id = integer(), global_window_id = integer(), chrom = character())
}

message("[STEP10v2-S2] Candidate regions: ", nrow(cand_dt),
        " (", nrow(membership_dt), " member windows)")

# =============================================================================
# WRITE OUTPUTS (identical file names and structure to monolithic STEP10v2)
# =============================================================================

meta_dt <- if (length(metadata_rows) > 0) rbindlist(metadata_rows, fill = TRUE) else data.table()

f1 <- paste0(outprefix, ".window_mds.tsv.gz")
f2 <- paste0(outprefix, ".candidate_regions.tsv.gz")
f3 <- paste0(outprefix, ".mds.rds")
f4 <- paste0(outprefix, ".mds_mode_metadata.tsv")
f5 <- paste0(outprefix, ".candidate_window_membership.tsv.gz")

# Create output dir
dir.create(dirname(outprefix), recursive = TRUE, showWarnings = FALSE)

fwrite(out_dt, f1, sep = "\t")
fwrite(cand_dt, f2, sep = "\t")

# THE critical .mds.rds structure — must match what snakes expect
saveRDS(list(
  dt = out_dt,
  candidate_regions = cand_dt,
  per_chr = all_chr_results,
  mds_mode = MDS_MODE
), f3)

fwrite(meta_dt, f4, sep = "\t")
fwrite(membership_dt, f5, sep = "\t")

message("\n[DONE] STEP10v2 Stage 2 — MDS merge complete (mode=", MDS_MODE, ")")
message("  ", f1)
message("  ", f2)
message("  ", f5)
message("  ", f3, "  ← snakes read this")
message("  ", f4)
