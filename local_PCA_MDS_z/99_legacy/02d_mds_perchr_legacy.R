#!/usr/bin/env Rscript

# =============================================================================
# STEP10_window_mds_outliers.R  (v7.4 — per-chromosome lostruct)
#
# Takes all STEP09 window PCA summaries, computes lostruct-style distances
# PER CHROMOSOME (not genome-wide), performs MDS per chromosome, then merges
# MDS coordinates and identifies outlier windows. Clusters nearby outliers
# into candidate inversion regions.
#
# v7.4 FIX: The original computed distances across all ~39K windows genome-wide
# in one O(n²) loop — ~770M pairs, hours of runtime. Inversions are within-
# chromosome, so cross-chromosome distances are meaningless. Now computes per-chr
# (~1-2K windows each, seconds per chr), then concatenates MDS results.
#
# Usage:
#   Rscript STEP10_window_mds_outliers.R \
#     <step09_rds_dir> <outprefix> [npc=2] [mds_dims=20] [z_thresh=3] \
#     [gap_bp=500000] [min_windows=3]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript STEP10_window_mds_outliers.R <step09_rds_dir> <outprefix> [npc=2] [mds_dims=20] [z_thresh=3] [gap_bp=500000] [min_windows=3]")
}

rds_dir     <- args[1]
outprefix   <- args[2]
npc         <- if (length(args) >= 3) as.integer(args[3])  else 4L  # v9.4: 2 → 4
mds_dims    <- if (length(args) >= 4) as.integer(args[4])  else 20L
z_thresh    <- if (length(args) >= 5) as.numeric(args[5])  else 3
gap_bp      <- if (length(args) >= 6) as.numeric(args[6])  else 500000
min_windows <- if (length(args) >= 7) as.integer(args[7])  else 3L

rds_files <- list.files(rds_dir, pattern = "\\.window_pca\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("No STEP09 .window_pca.rds files found in: ", rds_dir)

message("[INFO] Found ", length(rds_files), " STEP09 RDS files")

# ── Lostruct distance function (corrected) ─────────────────────────────────
# Frobenius distance between two rank-k covariance reconstructions:
#   ||S1 - S2||^2_F = sum(lam1^2) + sum(lam2^2) - 2 * sum_jk lam1_k * lam2_j * (u2_j . u1_k)^2
dist_sq_from_pcs <- function(values1, vectors1, values2, vectors2) {
  Xt <- crossprod(vectors2, vectors1)
  cross <- sum(outer(values2, values1) * Xt^2)
  sum(values1^2) + sum(values2^2) - 2 * cross
}

# ── Compute pairwise distances from a PCA data.frame (for ONE chromosome) ──
pc_dist_from_step09 <- function(pca_df, sample_names, npc, normalize = "L1") {
  values <- as.matrix(pca_df[, paste0("lam_", seq_len(npc)), drop = FALSE])

  vec_cols <- unlist(lapply(seq_len(npc), function(pc) paste0("PC_", pc, "_", sample_names)))
  vectors <- as.matrix(pca_df[, vec_cols, drop = FALSE])

  n_samples <- length(sample_names)

  if (normalize == "L1") {
    rs <- rowSums(abs(values))
    rs[rs == 0] <- NA
    values <- values / rs
  }

  n <- nrow(values)
  out <- matrix(NA_real_, n, n)

  emat <- function(u) matrix(u, nrow = n_samples, ncol = npc)

  for (i in seq_len(n)) {
    vi <- values[i, ]
    Ui <- emat(vectors[i, ])
    out[i, i] <- 0

    if (i < n) {
      for (j in (i + 1L):n) {
        vj <- values[j, ]
        Uj <- emat(vectors[j, ])

        if (anyNA(vi) || anyNA(vj) || anyNA(Ui) || anyNA(Uj)) {
          d2 <- NA_real_
        } else {
          d2 <- dist_sq_from_pcs(vi, Ui, vj, Uj)
          if (!is.finite(d2) || d2 < 0) d2 <- 0
        }
        out[i, j] <- sqrt(d2)
        out[j, i] <- out[i, j]
      }
    }
  }
  out
}

# ── Cluster outlier windows into candidate regions ─────────────────────────
cluster_outliers_bp <- function(dt, flag_col, gap_bp, min_windows) {
  idx <- which(dt[[flag_col]])
  if (length(idx) == 0) return(NULL)

  clusters <- list()
  cur <- idx[1]

  if (length(idx) > 1) {
    for (ii in idx[-1]) {
      prev_end   <- dt$end_bp[max(cur)]
      this_start <- dt$start_bp[ii]
      if ((this_start - prev_end) <= gap_bp) {
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

# ── Load all STEP09 outputs, grouped by chromosome ────────────────────────
chr_data <- list()
sample_names_ref <- NULL

for (f in sort(rds_files)) {
  obj <- readRDS(f)
  if (is.null(sample_names_ref)) {
    sample_names_ref <- obj$sample_names
  } else if (!identical(sample_names_ref, obj$sample_names)) {
    stop("Sample names differ across STEP09 files")
  }
  chr <- obj$chrom
  chr_data[[chr]] <- list(
    meta = as.data.table(obj$window_meta),
    pca  = as.data.table(obj$pca),
    chrom = chr
  )
}

total_windows <- sum(sapply(chr_data, function(x) nrow(x$meta)))
message("[INFO] Total windows across ", length(chr_data), " chromosomes: ", total_windows)

# ── Per-chromosome: distance → MDS ────────────────────────────────────────
global_id_offset <- 0L
all_chr_results <- list()

for (chr in names(chr_data)) {
  cd <- chr_data[[chr]]
  n_win <- nrow(cd$meta)

  if (n_win < 3) {
    message("[SKIP] ", chr, ": only ", n_win, " windows")
    global_id_offset <- global_id_offset + n_win
    next
  }

  message("[INFO] ", chr, ": ", n_win, " windows — computing lostruct distances...")
  t0 <- proc.time()

  # Assign global IDs
  cd$meta[, global_window_id := global_id_offset + seq_len(.N)]
  cd$pca[, global_window_id := global_id_offset + seq_len(.N)]

  # Merge meta + pca
  dt_chr <- merge(cd$meta, cd$pca, by = "global_window_id")

  # Compute pairwise distances (within-chromosome only)
  dmat <- pc_dist_from_step09(as.data.frame(dt_chr), sample_names_ref,
                               npc = npc, normalize = "L1")

  # Keep windows with all finite distances
  keep <- which(apply(dmat, 1, function(x) all(is.finite(x))))
  if (length(keep) < 3) {
    message("[SKIP] ", chr, ": only ", length(keep), " windows with finite distances")
    global_id_offset <- global_id_offset + n_win
    next
  }

  dt_keep <- dt_chr[keep]
  dmat_keep <- dmat[keep, keep, drop = FALSE]

  # MDS
  k_mds <- min(mds_dims, nrow(dmat_keep) - 1L)
  mds <- tryCatch(
    cmdscale(as.dist(dmat_keep), k = k_mds, eig = TRUE),
    error = function(e) { message("[WARN] MDS failed for ", chr, ": ", e$message); NULL }
  )

  if (is.null(mds)) {
    global_id_offset <- global_id_offset + n_win
    next
  }

  mds_dt <- as.data.table(mds$points)
  setnames(mds_dt, paste0("MDS", seq_len(ncol(mds_dt))))
  mds_dt[, global_window_id := dt_keep$global_window_id]

  # Z-scores and outlier flags (per-chromosome normalization)
  for (i in seq_len(ncol(mds$points))) {
    coln <- paste0("MDS", i)
    zcol <- paste0("MDS", i, "_z")
    vv   <- mds_dt[[coln]]
    mds_dt[[zcol]] <- (vv - mean(vv, na.rm = TRUE)) / sd(vv, na.rm = TRUE)
    mds_dt[[paste0("MDS", i, "_outlier")]] <- abs(mds_dt[[zcol]]) >= z_thresh
  }

  setkey(mds_dt, global_window_id)
  setkey(dt_keep, global_window_id)
  out_chr <- merge(dt_keep, mds_dt, by = "global_window_id")

  elapsed <- (proc.time() - t0)[3]
  message("[INFO] ", chr, ": done in ", round(elapsed, 1), "s (",
          length(keep), "/", n_win, " windows kept, k_mds=", k_mds, ")")

  all_chr_results[[chr]] <- list(out_dt = out_chr, dmat = dmat_keep, mds = mds)
  global_id_offset <- global_id_offset + n_win
}

# ── Merge all chromosomes ────────────────────────────────────────────────
if (length(all_chr_results) == 0) stop("No chromosomes produced results")

out_dt <- rbindlist(lapply(all_chr_results, function(x) x$out_dt), fill = TRUE)
message("[INFO] Merged: ", nrow(out_dt), " windows across ", length(all_chr_results), " chromosomes")

# ── Cluster outlier windows into candidate regions ─────────────────────────
# Use the per-chromosome MDS axes (already in out_dt)
mds_cols <- grep("^MDS\\d+$", names(out_dt), value = TRUE)
k_mds_max <- length(mds_cols)

cand_list <- list()
cand_id <- 1L

for (i in seq_len(k_mds_max)) {
  flag_col <- paste0("MDS", i, "_outlier")
  if (!(flag_col %in% names(out_dt))) next

  for (chr in unique(out_dt$chrom)) {
    sub <- out_dt[chrom == chr][order(start_bp)]
    if (!(flag_col %in% names(sub))) next
    if (all(is.na(sub[[flag_col]]))) next

    clusters <- cluster_outliers_bp(sub, flag_col, gap_bp = gap_bp,
                                     min_windows = min_windows)
    if (!is.null(clusters)) {
      for (cl in clusters) {
        xx <- sub[cl]
        cand_list[[length(cand_list) + 1]] <- data.table(
          candidate_id           = cand_id,
          mds_axis               = i,
          chrom                  = chr,
          start_bp               = min(xx$start_bp),
          end_bp                 = max(xx$end_bp),
          n_windows              = nrow(xx),
          first_global_window_id = min(xx$global_window_id),
          last_global_window_id  = max(xx$global_window_id)
        )
        cand_id <- cand_id + 1L
      }
    }
  }
}

cand_dt <- if (length(cand_list) > 0) rbindlist(cand_list) else data.table(
  candidate_id = integer(), mds_axis = integer(), chrom = character(),
  start_bp = numeric(), end_bp = numeric(), n_windows = integer(),
  first_global_window_id = integer(), last_global_window_id = integer()
)

message("[INFO] Candidate regions found: ", nrow(cand_dt))

# ── Write outputs ──────────────────────────────────────────────────────────
mds_out      <- paste0(outprefix, ".window_mds.tsv.gz")
cand_out     <- paste0(outprefix, ".candidate_regions.tsv.gz")
control_out  <- paste0(outprefix, ".candidate_regions_posmerge500kb.tsv.gz")
mds_rds_out  <- paste0(outprefix, ".mds.rds")

fwrite(out_dt, mds_out, sep = "\t")
fwrite(cand_dt, cand_out, sep = "\t")
fwrite(cand_dt, control_out, sep = "\t")  # explicit positional control copy
saveRDS(list(dt = out_dt, candidate_regions = cand_dt,
             per_chr = all_chr_results), mds_rds_out)

message("[DONE] Wrote:")
message("  ", mds_out)
message("  ", cand_out)
message("  ", control_out, "  ← positional merge 500 kb control")
message("  ", mds_rds_out)
