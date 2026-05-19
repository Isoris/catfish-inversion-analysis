#!/usr/bin/env Rscript

# =============================================================================
# 01b_local_pca_compute.R
#
# STAGE 1 (compute) — sliding-window local PCA on the per-chromosome dosage
# matrix. Designed to be run as a SLURM array, ONE TASK PER CHROMOSOME.
# Per-task work is fully independent; global window IDs are deferred to the
# merge step (01c).
#
# Pipeline position (path 1 — local-PCA z-blocks):
#   01a beagle_to_dosage  -> 01b LOCAL_PCA_COMPUTE  -> 01c local_pca_merge
#                                                  -> 02a/02b MDS
#                                                  -> 03 precompute  -> ...
#
# ── Inputs ────────────────────────────────────────────────────────────────
#   --dosage_dir <dir>   directory holding per-chr files from step 01a:
#                          <chr>.dosage.tsv.gz  (site x sample dosage matrix)
#                          <chr>.sites.tsv.gz   (per-site allele/position)
#   --outdir     <dir>   output dir (harmonized 01_local_pca/ slot;
#                        v10: no <outdir>/tmp/ subdir)
#   --chrom      <name>  chromosome label (e.g. C_gar_LG28)
#   [--winsize   100]    SNPs per window
#   [--step      20]     SNP step between windows
#   [--npc       4]      number of top eigenvectors kept per window
#                        (full eigenvalue spectrum is also kept for scree)
#
# ── Outputs (in <outdir>/) ────────────────────────────────────────────────
# v10: no tmp/ subdir; stage-1 files coexist with stage-2 finals via the
# _tmp filename suffix (ZO_D consumes _tmp and writes finalized .rds).
#   <chrom>.window_pca_tmp.rds   per-window PCA (top-npc eigvecs + scree),
#                                window_id = NA (filled in by 01c / ZO_D)
#   <chrom>.chr_meta.tsv.gz      per-window metadata, local indices only
#   <chrom>.pca_table.tsv.gz     per-window PCA loadings, window_id = NA
#   <chrom>.summary.tsv          one-row chromosome summary
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: STEP_A03_dense_registry_stage1.R v9.4-parallel)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(RSpectra)
})

# =============================================================================
# PARSE ARGS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

dosage_dir <- NULL
outdir     <- NULL
CHROM      <- NULL
WINSIZE    <- 100L
STEP       <- NA_integer_
NPC        <- 4L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--dosage_dir" && i < length(args)) {
    dosage_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--outdir" && i < length(args)) {
    outdir <- args[i + 1]; i <- i + 2L
  } else if (a == "--chrom" && i < length(args)) {
    CHROM <- args[i + 1]; i <- i + 2L
  } else if (a == "--winsize" && i < length(args)) {
    WINSIZE <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--step" && i < length(args)) {
    STEP <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--npc" && i < length(args)) {
    NPC <- as.integer(args[i + 1]); i <- i + 2L
  } else {
    i <- i + 1L
  }
}

if (is.null(dosage_dir) || is.null(outdir) || is.null(CHROM)) {
  stop("Usage: Rscript STEP_A03_dense_registry_stage1.R --dosage_dir <dir> --outdir <dir> --chrom <chr> ",
       "[--winsize 100] [--step 20] [--npc 4]")
}

if (is.na(STEP)) STEP <- WINSIZE
WINDOW_SET_ID <- paste0("snp", WINSIZE, "_step", STEP)

# v10: write per-chrom files directly into <outdir> (the 01_local_pca slot).
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SKIP-IF-DONE
# =============================================================================

meta_out <- file.path(outdir, paste0(CHROM, ".chr_meta.tsv.gz"))
rds_out  <- file.path(outdir, paste0(CHROM, ".window_pca_tmp.rds"))

if (file.exists(meta_out) && file.exists(rds_out)) {
  message("[A03-S1] ", CHROM, ": outputs already exist, skipping")
  message("  ", meta_out)
  message("  ", rds_out)
  quit(status = 0)
}

# =============================================================================
# PCA HELPER (v9.4: full spectrum stored)
# =============================================================================

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

  # Full eigendecomp: covmat is n×n where n = n_samples (~226), so cheap.
  # Gives both top-k (for lostruct) AND full scree (for variance diagnostics).
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

# =============================================================================
# PROGRESS BAR HELPER
# =============================================================================

make_progress <- function(total, label = "A03-S1") {
  t0 <- proc.time()[3]
  last_print <- 0
  PRINT_INTERVAL <- 2  # seconds between progress updates

  function(current) {
    now <- proc.time()[3]
    if (now - last_print < PRINT_INTERVAL && current < total) return(invisible())
    last_print <<- now

    elapsed <- now - t0
    pct <- 100 * current / total
    if (current > 0) {
      eta_sec <- elapsed * (total - current) / current
      eta_str <- if (eta_sec >= 3600) {
        sprintf("%dh%02dm", as.integer(eta_sec / 3600), as.integer((eta_sec %% 3600) / 60))
      } else if (eta_sec >= 60) {
        sprintf("%dm%02ds", as.integer(eta_sec / 60), as.integer(eta_sec %% 60))
      } else {
        sprintf("%ds", as.integer(eta_sec))
      }
    } else {
      eta_str <- "?"
    }

    bar_width <- 30
    filled <- as.integer(bar_width * current / total)
    bar <- paste0("[", strrep("=", filled), strrep(" ", bar_width - filled), "]")

    message(sprintf("\r[%s] %s  %s  %d/%d  (%.1f%%)  ETA %s",
                    label, CHROM, bar, current, total, pct, eta_str))
  }
}

# =============================================================================
# LOAD DATA
# =============================================================================

message("[A03-S1] Dense window PCA — ", CHROM)
message("[A03-S1] Window: ", WINSIZE, " SNPs, Step: ", STEP, " SNPs, npc=", NPC)
message("[A03-S1] Set ID: ", WINDOW_SET_ID)

sf <- file.path(dosage_dir, paste0(CHROM, ".sites.tsv.gz"))
dos_file <- file.path(dosage_dir, paste0(CHROM, ".dosage.tsv.gz"))

if (!file.exists(sf))       stop("Missing sites file: ", sf)
if (!file.exists(dos_file)) stop("Missing dosage file: ", dos_file)

message("[A03-S1] Reading sites + dosage...")
t_load <- proc.time()[3]
sites <- fread(sf)
dos   <- fread(dos_file)
message("[A03-S1] Loaded in ", round(proc.time()[3] - t_load, 1), "s")

# Align
if (!identical(sites$marker, dos$marker)) {
  setkeyv(sites, "marker"); setkeyv(dos, "marker")
  dos <- dos[sites$marker]
}

sample_names <- setdiff(names(dos), "marker")
X <- as.matrix(dos[, ..sample_names])
storage.mode(X) <- "double"
n_snps <- nrow(X)
n_samples <- ncol(X)

starts <- seq(1L, n_snps - WINSIZE + 1L, by = STEP)
n_windows <- length(starts)

if (n_windows < 1) {
  message("[A03-S1] ", CHROM, " — ", n_snps, " SNPs, not enough for one window")
  fwrite(data.table(chrom = CHROM, n_snps = n_snps, n_windows = 0L),
         file.path(outdir, paste0(CHROM, ".summary.tsv")), sep = "\t")
  quit(status = 0)
}

message("[A03-S1] ", CHROM, ": ", n_snps, " SNPs → ", n_windows, " windows")
message("[A03-S1] scree_full storage: ", n_samples, " eigenvalues × ", n_windows,
        " windows ≈ ", round(n_samples * 8 * n_windows / 1e6, 1), " MB")

# =============================================================================
# PCA LOOP WITH PROGRESS (v9.4 — full spectrum)
# =============================================================================

meta_list <- vector("list", n_windows)
pca_rows  <- vector("list", n_windows)
scree_mat <- matrix(NA_real_, nrow = n_windows, ncol = n_samples)
pve_mat   <- matrix(NA_real_, nrow = n_windows, ncol = n_samples)

pb <- make_progress(n_windows)

for (wi in seq_len(n_windows)) {
  i1 <- starts[wi]
  i2 <- min(i1 + WINSIZE - 1L, n_snps)
  xw <- X[i1:i2, , drop = FALSE]
  sw <- sites[i1:i2]
  pca <- cov_pca_one_window(xw, NPC)

  # Local metadata — window_id = NA (assigned by Stage 2)
  meta_list[[wi]] <- data.table(
    window_id        = NA_integer_,
    window_set_id    = WINDOW_SET_ID,
    chrom            = CHROM,
    window_index_chr = wi - 1L,
    start_bp         = sw$pos[1],
    end_bp           = sw$pos[nrow(sw)],
    center_bp        = as.integer((sw$pos[1] + sw$pos[nrow(sw)]) / 2),
    start_snp_idx    = i1,
    end_snp_idx      = i2,
    n_snps           = nrow(xw),
    step_snps        = STEP,
    window_snps      = WINSIZE,
    focal_mode       = "focal_dense"
  )

  rowvec <- c(total = pca$total, setNames(pca$eigval_top, paste0("lam_", seq_len(NPC))))
  for (pc in seq_len(NPC)) {
    rowvec <- c(rowvec, setNames(pca$eigvec_top[, pc], paste0("PC_", pc, "_", sample_names)))
  }
  pca_rows[[wi]] <- rowvec
  scree_mat[wi, ] <- pca$scree_full
  pve_mat[wi, ]   <- pca$pve_full

  pb(wi)
}
message("")  # newline after progress bar

# =============================================================================
# ASSEMBLE + WRITE
# =============================================================================

chr_meta <- rbindlist(meta_list)
pca_mat  <- as.data.frame(do.call(rbind, pca_rows))
pca_mat$window_id <- NA_integer_
pca_mat <- pca_mat[, c("window_id", setdiff(names(pca_mat), "window_id"))]

# Per-chr meta
fwrite(chr_meta, meta_out, sep = "\t")

# Per-chr PCA (temp — Stage 2 will patch window_id and re-save to final location)
saveRDS(list(
  chrom        = CHROM,
  sample_names = sample_names,
  npc          = NPC,
  winsize      = WINSIZE,
  step         = STEP,
  window_meta  = chr_meta,
  pca          = pca_mat,
  scree_full   = scree_mat,
  pve_full     = pve_mat
), rds_out)

# PCA table
pca_tsv <- file.path(outdir, paste0(CHROM, ".pca_table.tsv.gz"))
fwrite(as.data.table(pca_mat), pca_tsv, sep = "\t")

# Summary
fwrite(data.table(
  chrom = CHROM, n_snps = n_snps, n_windows = n_windows, n_samples = n_samples,
  step_snps = STEP, window_snps = WINSIZE, npc = NPC
), file.path(outdir, paste0(CHROM, ".summary.tsv")), sep = "\t")

message("[A03-S1] DONE — ", CHROM, ": ", n_windows, " windows")
message("  ", meta_out)
message("  ", rds_out)
