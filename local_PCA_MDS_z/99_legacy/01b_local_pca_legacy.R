#!/usr/bin/env Rscript

# =============================================================================
# STEP_A02_local_pca_windows_by_chr.R
#
# Read one chromosome's dosage file + sites file, split into fixed SNP windows,
# compute local PCA summaries per window: top-NPC eigenvectors of the
# individual covariance matrix + the FULL eigenvalue spectrum (for downstream
# scree plots / variance-explained diagnostics).
#
# v9.4 changes (vs. legacy v7.4):
#   - npc default raised from 4 → 4 (kept; was already 4)
#   - Full eigenvalue spectrum stored per window (scree_full = lam_1..lam_N).
#     Eigenvectors still capped at top npc to keep RDS size bounded.
#   - Eigenvalues stored as proportion-of-variance-explained (PVE) AND as
#     raw eigenvalues. Both cheap, both useful.
#
# Usage:
#   Rscript STEP_A02_local_pca_windows_by_chr.R \
#     <sites.tsv.gz> <dosage.tsv.gz> <outprefix> [npc=4] [winsize=100]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(RSpectra)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript STEP_A02_local_pca_windows_by_chr.R <sites.tsv.gz> <dosage.tsv.gz> <outprefix> [npc=4] [winsize=100]")
}

sites_file  <- args[1]
dosage_file <- args[2]
outprefix   <- args[3]
npc         <- if (length(args) >= 4) as.integer(args[4]) else 4L
winsize     <- if (length(args) >= 5) as.integer(args[5]) else 100L

# ── Read inputs ─────────────────────────────────────────────────────────────
sites <- fread(sites_file)
dos   <- fread(dosage_file)

req_sites <- c("marker", "chrom", "pos", "allele1", "allele2")
miss_sites <- setdiff(req_sites, names(sites))
if (length(miss_sites) > 0) stop("Missing required columns in sites file: ", paste(miss_sites, collapse = ", "))
if (!("marker" %in% names(dos))) stop("Dosage file must contain 'marker' column")

# Align sites and dosage by marker order
if (!identical(sites$marker, dos$marker)) {
  setkeyv(sites, "marker")
  setkeyv(dos, "marker")
  dos <- dos[sites$marker]
  if (!identical(sites$marker, dos$marker)) stop("Sites and dosage markers do not match")
}

sample_names <- setdiff(names(dos), "marker")
chroms <- unique(sites$chrom)
if (length(chroms) != 1) stop("STEP_A02 expects one chromosome per input file")
chrom <- chroms[1]

# ── Build numeric matrix (SNPs × samples) ──────────────────────────────────
X <- as.matrix(dos[, ..sample_names])
storage.mode(X) <- "double"

n_snps <- nrow(X)
n_samples <- ncol(X)
n_windows <- floor(n_snps / winsize)
if (n_windows < 1) stop("Not enough SNPs for one full window")

message("[A02] Chromosome:  ", chrom)
message("[A02] SNPs:        ", n_snps)
message("[A02] Samples:     ", n_samples)
message("[A02] Windows:     ", n_windows, "  (", winsize, " SNPs/window)")
message("[A02] npc kept:    ", npc, "  (eigenvectors)")
message("[A02] scree_full:  ", n_samples, " eigenvalues per window  (~",
        round(n_samples * 8 * n_windows / 1e6, 1), " MB scree storage)")

# ── Per-window PCA function ────────────────────────────────────────────────
# Input x:  SNPs × samples matrix
# Returns:
#   total          = sum of squared covariance entries (unchanged from v7.4)
#   eigval_top     = top npc eigenvalues (kept for back-compat with B01/lostruct)
#   eigvec_top     = top npc eigenvectors (n_samples × npc)
#   scree_full     = ALL eigenvalues (length n_samples), descending
#   pve_full       = scree_full / sum(scree_full), descending
cov_pca_one_window <- function(x, k) {
  rm <- rowMeans(x, na.rm = TRUE)
  xc <- x - rm
  covmat <- suppressWarnings(cov(xc, use = "pairwise.complete.obs"))

  ns <- ncol(x)
  if (anyNA(covmat) || nrow(covmat) < k) {
    return(list(
      total      = NA_real_,
      eigval_top = rep(NA_real_, k),
      eigvec_top = matrix(NA_real_, nrow = ns, ncol = k),
      scree_full = rep(NA_real_, ns),
      pve_full   = rep(NA_real_, ns)
    ))
  }

  total <- sum(covmat^2)

  # Always do full eigendecomposition: cheap on n×n where n = n_samples (~226).
  # This gives both the top-k for downstream lostruct AND the full scree.
  ee_full <- eigen(covmat, symmetric = TRUE)
  vals_full <- ee_full$values
  scree_sum <- sum(vals_full[is.finite(vals_full) & vals_full > 0])
  pve_full <- if (is.finite(scree_sum) && scree_sum > 0) {
    pmax(0, vals_full) / scree_sum
  } else {
    rep(NA_real_, length(vals_full))
  }

  vals_top <- vals_full[seq_len(min(k, length(vals_full)))]
  vecs_top <- ee_full$vectors[, seq_len(min(k, ncol(ee_full$vectors))), drop = FALSE]

  list(
    total      = total,
    eigval_top = vals_top,
    eigvec_top = vecs_top,
    scree_full = vals_full,
    pve_full   = pve_full
  )
}

# ── Process windows ────────────────────────────────────────────────────────
meta_list  <- vector("list", n_windows)
pca_rows   <- vector("list", n_windows)
scree_mat  <- matrix(NA_real_, nrow = n_windows, ncol = n_samples)
pve_mat    <- matrix(NA_real_, nrow = n_windows, ncol = n_samples)

for (w in seq_len(n_windows)) {
  i1 <- (w - 1L) * winsize + 1L
  i2 <- w * winsize
  xw <- X[i1:i2, , drop = FALSE]
  sw <- sites[i1:i2]
  pca <- cov_pca_one_window(xw, npc)

  meta_list[[w]] <- data.table(
    window_id       = w,
    chrom           = chrom,
    snp_start_index = i1,
    snp_end_index   = i2,
    n_snps          = nrow(xw),
    start_bp        = sw$pos[1],
    end_bp          = sw$pos[nrow(sw)],
    start_marker    = sw$marker[1],
    end_marker      = sw$marker[nrow(sw)]
  )

  rowvec <- c(total = pca$total, setNames(pca$eigval_top, paste0("lam_", seq_len(npc))))
  for (pc in seq_len(npc)) {
    rowvec <- c(rowvec, setNames(pca$eigvec_top[, pc], paste0("PC_", pc, "_", sample_names)))
  }
  pca_rows[[w]] <- rowvec
  scree_mat[w, ] <- pca$scree_full
  pve_mat[w, ]   <- pca$pve_full
}

# ── Write outputs ──────────────────────────────────────────────────────────
window_meta <- rbindlist(meta_list)
pca_mat     <- do.call(rbind, pca_rows)
pca_mat     <- as.data.frame(pca_mat)
pca_mat$window_id <- seq_len(n_windows)
pca_mat     <- pca_mat[, c("window_id", setdiff(names(pca_mat), "window_id"))]

meta_out    <- paste0(outprefix, ".window_meta.tsv.gz")
pca_rds_out <- paste0(outprefix, ".window_pca.rds")
pca_tsv_out <- paste0(outprefix, ".window_pca.tsv.gz")

fwrite(window_meta, meta_out, sep = "\t")
saveRDS(list(
  chrom        = chrom,
  sample_names = sample_names,
  npc          = npc,
  winsize      = winsize,
  window_meta  = window_meta,
  pca          = pca_mat,
  scree_full   = scree_mat,   # n_windows × n_samples — eigenvalues
  pve_full     = pve_mat      # n_windows × n_samples — proportion variance
), pca_rds_out)
fwrite(as.data.table(pca_mat), pca_tsv_out, sep = "\t")

message("[A02] DONE")
message("  ", meta_out)
message("  ", pca_rds_out)
message("  ", pca_tsv_out)
