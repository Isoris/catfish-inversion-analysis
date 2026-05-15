#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_H_classify_carriers.R
# =============================================================================
# Per-candidate carrier classification via k-means on local-PCA loadings.
# For each L2 candidate, take the per-sample mean PC1 over the candidate's
# window range, run k-means with k = 2..max_k, pick best k via approximate
# silhouette, label clusters by ascending centre as LOW_DIV / MID_DIV /
# HIGH_DIV (or DIV_TIER_k for k > 3). These cluster labels feed TR_I as the
# carrier partition for per-band CUSUM.
#
# Dual-scale aware (May 2026): in the dual-scale design, L2 envelopes come
# from the COARSE precomp (TR_F), but the k-means inside each envelope wants
# the DENSE per-window per-sample PC1 for finer carrier resolution. Default
# precomp dir is therefore $OUTROOT/precomp_dense if it exists, falling back
# to $OUTROOT/precomp for single-scale legacy runs. The L2 envelope's bp
# range is scale-agnostic (start_bp/end_bp), so windows_in_range() resolves
# it against whichever precomp's window grid is loaded.
#
# Reads:   <PRECOMP_DIR>/<chr>.precomp.rds          (precomp_dense if available)
#          <L2_DIR>/<chr>.L2_envelopes.tsv          (output of TR_F, COARSE scale)
# Writes:  <OUT_DIR>/<chr>.carrier_assignments.tsv
#            One row per (sample × candidate). Columns:
#              candidate_id  chrom  sample_id  band  k_chosen  silhouette
#              mean_pc1  mean_pc2  inside_l2_windows
#
# Defaults:
#   PRECOMP_DIR = $OUTROOT/precomp_dense   (or $OUTROOT/precomp if dense missing)
#   L2_DIR      = $OUTROOT/L2_detect       (coarse-scale envelopes from TR_F)
#   OUT_DIR     = $OUTROOT/carriers
#
# Usage:
#   Rscript STEP_TR_H_classify_carriers.R --chr <CHR>
#                                          [--precomp_dir <dir>]
#                                          [--l2_dir <dir>]
#                                          [--out_dir <dir>]
#                                          [--max-k 3]
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}

CHR          <- get_arg("--chr")
PRECOMP_DIR  <- get_arg("--precomp_dir")
L2_DIR       <- get_arg("--l2_dir")
OUT_DIR      <- get_arg("--out_dir")
MAX_K        <- as.integer(get_arg("--max-k", "3"))
# Precomp filename suffix — defaults to ".precomp.rds" (theta_pi) but can be
# pointed at GHSL outputs by passing ".ghsl_precomp.rds".
PRECOMP_SUFFIX <- get_arg("--precomp_suffix", ".precomp.rds")

OUTROOT <- Sys.getenv("OUTROOT", unset = NA)
if (is.na(PRECOMP_DIR)) {
  dense_dir <- file.path(OUTROOT, "precomp_dense")
  PRECOMP_DIR <- if (dir.exists(dense_dir)) dense_dir else file.path(OUTROOT, "precomp")
}
if (is.na(L2_DIR))      L2_DIR      <- file.path(OUTROOT, "L2_detect")
if (is.na(OUT_DIR))     OUT_DIR     <- file.path(OUTROOT, "carriers")
message("[TR_H] precomp dir: ", PRECOMP_DIR, " suffix=", PRECOMP_SUFFIX)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

precomp_pat <- paste0(gsub("\\.", "\\\\.", PRECOMP_SUFFIX, fixed = FALSE), "$")
chroms <- if (!is.na(CHR)) CHR else
  sub(precomp_pat, "",
      list.files(PRECOMP_DIR, pattern = precomp_pat))

approx_sil <- function(v, cl, k) {
  sil <- numeric(length(v))
  for (j in seq_along(v)) {
    own <- cl$cluster[j]
    own_members <- v[cl$cluster == own & seq_along(v) != j]
    a_j <- if (length(own_members) > 0) mean(abs(v[j] - own_members)) else 0
    od <- vapply(setdiff(seq_len(k), own),
                 function(o) {
                   m <- v[cl$cluster == o]
                   if (length(m) > 0) mean(abs(v[j] - m)) else NA_real_
                 }, numeric(1))
    b_j <- if (length(od) > 0) min(od, na.rm = TRUE) else 0
    sil[j] <- if (max(a_j, b_j) > 0) (b_j - a_j) / max(a_j, b_j) else 0
  }
  mean(sil)
}

label_clusters <- function(centres, k) {
  o <- order(centres)
  if      (k == 2L) c("LOW_DIV", "HIGH_DIV")[match(seq_len(k), o)]
  else if (k == 3L) c("LOW_DIV", "MID_DIV", "HIGH_DIV")[match(seq_len(k), o)]
  else paste0("DIV_TIER_", match(seq_len(k), o))
}

# Map an L2 row's bp range back to the precomp window-idx range.
windows_in_range <- function(dt, l2_row) {
  if ("win_start" %in% names(l2_row) && "win_end" %in% names(l2_row) &&
      is.finite(l2_row$win_start) && is.finite(l2_row$win_end)) {
    seq.int(as.integer(l2_row$win_start), as.integer(l2_row$win_end))
  } else if ("start_bp" %in% names(l2_row) && "end_bp" %in% names(l2_row)) {
    which(dt$end_bp >= l2_row$start_bp & dt$start_bp <= l2_row$end_bp)
  } else {
    integer(0)
  }
}

for (chrom in chroms) {
  rds <- file.path(PRECOMP_DIR, paste0(chrom, PRECOMP_SUFFIX))
  l2f <- file.path(L2_DIR,      paste0(chrom, ".L2_envelopes.tsv"))
  if (!file.exists(rds) || !file.exists(l2f)) {
    message("[TR_H] ", chrom, ": missing precomp or L2 — skip"); next
  }
  precomp <- readRDS(rds); dt <- precomp$dt
  l2_dt <- fread(l2f)
  if (nrow(l2_dt) == 0L) {
    fwrite(data.table(),
           file.path(OUT_DIR, paste0(chrom, ".carrier_assignments.tsv")), sep = "\t")
    next
  }

  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  pc2_cols <- grep("^PC_2_", names(dt), value = TRUE)
  if (length(pc1_cols) == 0L) {
    message("[TR_H] ", chrom, ": precomp has no PC_1_<sample> columns — skip"); next
  }
  pc1_mat <- as.matrix(dt[, ..pc1_cols])
  pc2_mat <- if (length(pc2_cols) > 0) as.matrix(dt[, ..pc2_cols]) else NULL
  sample_ids <- sub("^PC_1_", "", pc1_cols)

  # candidate_id may or may not be in z-blocks 06's output; fall back to row idx.
  if (!"candidate_id" %in% names(l2_dt)) {
    l2_dt[, candidate_id := paste0(chrom, "_C_", sprintf("%03d", seq_len(.N)))]
  }

  out <- list()
  for (k_row in seq_len(nrow(l2_dt))) {
    cid    <- as.character(l2_dt$candidate_id[k_row])
    in_win <- windows_in_range(dt, l2_dt[k_row])
    if (length(in_win) < 3L) next

    mean_pc1 <- colMeans(pc1_mat[in_win, , drop = FALSE], na.rm = TRUE)
    mean_pc2 <- if (!is.null(pc2_mat))
                  colMeans(pc2_mat[in_win, , drop = FALSE], na.rm = TRUE)
                else rep(NA_real_, length(mean_pc1))
    valid <- is.finite(mean_pc1)
    if (sum(valid) < 10L) next
    v <- mean_pc1[valid]; sids <- sample_ids[valid]

    best_k <- NA_integer_; best_sil <- -Inf; best_cl <- NULL
    for (k in seq.int(2L, min(MAX_K, length(unique(round(v, 6))) - 1L))) {
      if (k >= length(v)) next
      cl <- tryCatch(kmeans(v, centers = k, nstart = 25L, iter.max = 100L),
                     error = function(e) NULL)
      if (is.null(cl)) next
      s <- approx_sil(v, cl, k)
      if (s > best_sil) { best_sil <- s; best_k <- k; best_cl <- cl }
    }
    if (is.null(best_cl)) next
    labels <- label_clusters(best_cl$centers[, 1], best_k)

    out[[length(out) + 1L]] <- data.table(
      candidate_id      = cid,
      chrom             = chrom,
      sample_id         = sids,
      band              = labels[best_cl$cluster],
      k_chosen          = as.integer(best_k),
      silhouette        = round(best_sil, 4),
      mean_pc1          = round(v, 6),
      mean_pc2          = round(mean_pc2[valid], 6),
      inside_l2_windows = length(in_win)
    )
  }
  out_dt <- if (length(out) > 0) rbindlist(out) else data.table()
  fwrite(out_dt, file.path(OUT_DIR, paste0(chrom, ".carrier_assignments.tsv")), sep = "\t")
  message(sprintf("[TR_H] %s: classified %d candidates (%d sample-rows)",
                  chrom, length(unique(out_dt$candidate_id)), nrow(out_dt)))
}

message("[TR_H] DONE")