#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_E_classify_carriers.R   (drafts/v5)
# =============================================================================
# Per-candidate carrier classification via k-means on local-PCA loadings.
# Mirrors the z-blocks per-window k=3 band classification but applied at
# the L2-candidate level (not per-window).
#
# For each L2 candidate, take the focal-window-or-mean PC1 loading per
# sample, run k-means with k = 2..max_k, pick best k via approximate
# silhouette, label clusters by ascending centre as LOW/MID/HIGH (or
# DIV_TIER_1..k for k > 3). These cluster labels feed STEP_TR_F as the
# carrier partition for per-band CUSUM.
#
# Reads:   <OUTROOT>/03_per_chrom/<chr>/precomp.rds
#          <OUTROOT>/03_per_chrom/<chr>/L2_envelopes.tsv
# Writes:  <OUTROOT>/03_per_chrom/<chr>/carrier_assignments.tsv
#            One row per (sample × candidate). Columns:
#              candidate_id  chrom  sample_id  band  k_chosen  silhouette
#              mean_pc1  mean_pc2  inside_l2_windows
# Usage:   Rscript STEP_TR_E_classify_carriers.R --chrom <CHR>
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
CHROM <- NULL
MAX_K <- 3L
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--chrom" && i < length(args)) { CHROM <- args[i + 1]; i <- i + 2L }
  else if (a == "--max-k" && i < length(args)) { MAX_K <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}

OUTROOT <- Sys.getenv("OUTROOT", unset = NA)
stopifnot(!is.na(OUTROOT))
per_chrom_dir <- file.path(OUTROOT, "03_per_chrom")
chroms <- if (!is.null(CHROM)) CHROM else list.files(per_chrom_dir, pattern = "^C_gar_LG[0-9]+$")

# Approximate silhouette for 1D k-means; mirrors v4's helper.
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

label_clusters_by_centre <- function(centres, k) {
  o <- order(centres)
  if (k == 2L) c("LOW_DIV", "HIGH_DIV")[match(seq_len(k), o)]
  else if (k == 3L) c("LOW_DIV", "MID_DIV", "HIGH_DIV")[match(seq_len(k), o)]
  else paste0("DIV_TIER_", match(seq_len(k), o))
}

for (chrom in chroms) {
  rds <- file.path(per_chrom_dir, chrom, "precomp.rds")
  l2f <- file.path(per_chrom_dir, chrom, "L2_envelopes.tsv")
  if (!file.exists(rds) || !file.exists(l2f)) {
    message("[TR_E] ", chrom, ": missing precomp or L2 — skip"); next
  }
  precomp <- readRDS(rds); dt <- precomp$dt
  l2_dt <- fread(l2f)
  if (nrow(l2_dt) == 0L) {
    fwrite(data.table(), file.path(per_chrom_dir, chrom, "carrier_assignments.tsv"), sep = "\t")
    next
  }
  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  pc2_cols <- grep("^PC_2_", names(dt), value = TRUE)
  if (length(pc1_cols) == 0L) {
    message("[TR_E] ", chrom, ": precomp has no PC_1_<sample> columns — skip"); next
  }
  pc1_mat <- as.matrix(dt[, ..pc1_cols])
  pc2_mat <- if (length(pc2_cols) > 0) as.matrix(dt[, ..pc2_cols]) else NULL
  sample_ids <- sub("^PC_1_", "", pc1_cols)

  out <- list()
  for (k_row in seq_len(nrow(l2_dt))) {
    cid    <- l2_dt$candidate_id[k_row]
    w_lo   <- l2_dt$win_start[k_row]
    w_hi   <- l2_dt$win_end[k_row]
    in_win <- w_lo:w_hi
    if (length(in_win) < 3L) next

    # Mean PC1 over the candidate window range, per sample.
    mean_pc1 <- colMeans(pc1_mat[in_win, , drop = FALSE], na.rm = TRUE)
    mean_pc2 <- if (!is.null(pc2_mat)) colMeans(pc2_mat[in_win, , drop = FALSE], na.rm = TRUE)
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
    labels <- label_clusters_by_centre(best_cl$centers[, 1], best_k)
    band   <- labels[best_cl$cluster]

    out[[length(out) + 1L]] <- data.table(
      candidate_id      = cid,
      chrom             = chrom,
      sample_id         = sids,
      band              = band,
      k_chosen          = as.integer(best_k),
      silhouette        = round(best_sil, 4),
      mean_pc1          = round(v, 6),
      mean_pc2          = round(mean_pc2[valid], 6),
      inside_l2_windows = length(in_win)
    )
  }
  out_dt <- if (length(out) > 0) rbindlist(out) else data.table()
  fwrite(out_dt, file.path(per_chrom_dir, chrom, "carrier_assignments.tsv"), sep = "\t")
  message(sprintf("[TR_E] %s: classified %d candidates (%d sample-rows)",
                  chrom, length(unique(out_dt$candidate_id)), nrow(out_dt)))
}

message("[TR_E] DONE")