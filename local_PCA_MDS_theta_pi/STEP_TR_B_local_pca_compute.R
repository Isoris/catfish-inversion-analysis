#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_B_local_pca_compute.R  (v6)
# =============================================================================
# Per-window heteroscedastic local PCA on the theta-pi matrix produced by
# STEP_TR_A. Anchor-flipped sign alignment is applied here so PC scores are
# directly comparable across windows.
#
# v6 change vs v5: the per-window SVD loop in local_pca_npc() is now
# parallelized with parallel::mclapply. Numerical output is bit-equivalent
# to v5 (same SVD per window; only the loop driver changed). Single-core
# fallback when N_CORES <= 1 or when running on Windows.
#
# Inputs (configured via 00_theta_config.sh):
#   $THETA_TSV_DIR/theta_native.<CHR>.<SCALE>.tsv.gz   (from STEP_TR_A)
#   $SAMPLE_LIST                                       cohort id list
#
# Output (under $OUT_LOCAL_PCA_DIR, default = $OUTROOT/01_local_pca):
#   <chr>.window_pca.rds   per-chrom local PCA bundle (schema unchanged):
#     $pcs            list of NPC matrices (n_samp x n_win), anchor-flipped
#     $lambda         (n_win x NPC) eigenvalues from each window's SVD
#     $theta_z_direct (length n_win) per-window |Z| from per-sample dev
#     $window_median  (length n_win) per-window cohort median theta_pi
#     $win_grid       data.table: chrom, window_idx, start_bp, end_bp
#     $anchor_idx     1-based anchor window index (used for the flip)
#     $unflipped_windows  1-based indices that couldn't be sign-aligned
#     $chrom, $n_windows, $n_samples, $npc, $pad, $sample_order
#
# Knobs:
#   --chrom <CHR>     single chrom (default: iterate CHROM_LIST from env)
#   --pad <int>       local-PCA neighbourhood half-width (default 1)
#   --npc <int>       number of PCs to keep (default 4; or $LOCAL_PCA_NPC)
#   --cores <int>     mclapply mc.cores (default $N_CORES or 8)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

PAD     <- 1L
NPC     <- 4L
N_CORES <- as.integer(Sys.getenv("N_CORES", unset = "8"))
if (!is.finite(N_CORES) || N_CORES < 1L) N_CORES <- 1L

env_npc <- suppressWarnings(as.integer(Sys.getenv("LOCAL_PCA_NPC", unset = "")))
if (length(env_npc) == 1L && is.finite(env_npc) && env_npc > 0L) NPC <- env_npc

args <- commandArgs(trailingOnly = TRUE)
CHROM <- NULL
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--chrom" && i < length(args)) { CHROM   <- args[i + 1];        i <- i + 2L }
  else if (a == "--pad"   && i < length(args)) { PAD     <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--npc"   && i < length(args)) { NPC     <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--cores" && i < length(args)) { N_CORES <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}
if (.Platform$OS.type == "windows" && N_CORES > 1L) {
  message("[TR_B PCA] Windows detected: forcing mc.cores=1 (mclapply is unix-only)")
  N_CORES <- 1L
}

THETA_TSV_DIR <- Sys.getenv("THETA_TSV_DIR", unset = NA)
PESTPG_SCALE  <- Sys.getenv("PESTPG_SCALE",  unset = "win10000.step2000")
OUTROOT       <- Sys.getenv("OUTROOT",       unset = NA)
SAMPLE_LIST   <- Sys.getenv("SAMPLE_LIST",   unset = NA)
OUT_LOCAL_PCA_DIR <- Sys.getenv("OUT_LOCAL_PCA_DIR",
                                unset = file.path(OUTROOT, "01_local_pca"))
stopifnot(!is.na(THETA_TSV_DIR), !is.na(OUTROOT), !is.na(SAMPLE_LIST))

if (is.null(CHROM)) {
  cl <- Sys.getenv("CHROM_LIST", unset = "")
  CHROM_LIST <- if (nchar(cl) > 0) strsplit(cl, "[ ,]+")[[1]] else
    sprintf("C_gar_LG%02d", 1:28)
} else {
  CHROM_LIST <- CHROM
}

dir.create(OUT_LOCAL_PCA_DIR, recursive = TRUE, showWarnings = FALSE)

sample_order <- readLines(SAMPLE_LIST)
sample_order <- sample_order[nchar(sample_order) > 0]
n_samp <- length(sample_order)
message("[TR_B PCA] cohort: ", n_samp, " samples, NPC=", NPC,
        ", PAD=", PAD, ", cores=", N_CORES)
message("[TR_B PCA] chroms: ", length(CHROM_LIST))

# =============================================================================
# Helpers
# =============================================================================

load_theta_matrix <- function(chrom) {
  tsv <- file.path(THETA_TSV_DIR, sprintf("theta_native.%s.%s.tsv.gz", chrom, PESTPG_SCALE))
  if (!file.exists(tsv)) return(NULL)
  long_dt <- fread(tsv)
  this_chrom <- chrom
  long_dt <- long_dt[chrom == this_chrom]
  if (nrow(long_dt) == 0) return(NULL)
  win_grid <- unique(long_dt[, .(window_idx, start_bp, end_bp)])
  setkey(win_grid, window_idx)
  n_win <- nrow(win_grid)

  samp_to_row <- setNames(seq_along(sample_order), sample_order)
  win_to_col  <- setNames(seq_along(win_grid$window_idx), as.character(win_grid$window_idx))

  theta_mat   <- matrix(NA_real_,    nrow = n_samp, ncol = n_win, dimnames = list(sample_order, NULL))
  n_sites_mat <- matrix(NA_integer_, nrow = n_samp, ncol = n_win, dimnames = list(sample_order, NULL))

  rows <- samp_to_row[long_dt$sample]
  cols <- win_to_col[as.character(long_dt$window_idx)]
  good <- !is.na(rows) & !is.na(cols)
  theta_mat[cbind(rows[good], cols[good])]   <- long_dt$theta_pi[good]
  n_sites_mat[cbind(rows[good], cols[good])] <- long_dt$n_sites[good]
  list(theta_mat = theta_mat, n_sites_mat = n_sites_mat, win_grid = win_grid)
}

# Per-window heteroscedastic local PCA, NPC PCs. v6: parallelized with
# mclapply. Each worker produces a per-window result list and the driver
# stitches them into the pcs[[k]] and lambda matrices. Numerical output is
# bit-equivalent to v5 (same per-window SVD, identical pad/weighting/centering).
local_pca_npc <- function(theta_mat, n_sites_mat, pad, npc, n_cores) {
  n_samp <- nrow(theta_mat); n_win <- ncol(theta_mat)
  ns_med <- median(n_sites_mat, na.rm = TRUE)
  if (!is.finite(ns_med) || ns_med <= 0) ns_med <- 1

  one_window <- function(wi) {
    lo <- max(1L, wi - pad); hi <- min(n_win, wi + pad)
    block <- theta_mat[, lo:hi, drop = FALSE]
    ok <- complete.cases(block)
    if (sum(ok) < max(20L, ncol(block) + 2L)) {
      return(list(ok = NULL, u = NULL, lam = rep(NA_real_, npc)))
    }
    block_ok <- block[ok, , drop = FALSE]
    nf <- n_sites_mat[ok, wi]
    nf[!is.finite(nf) | nf <= 0] <- 1L
    w <- sqrt(nf / ns_med)
    bw <- block_ok * w
    centred <- sweep(bw, 2, colMeans(bw), FUN = "-")
    nu_req <- min(npc, ncol(centred))
    sv <- tryCatch(svd(centred, nu = nu_req, nv = 0), error = function(e) NULL)
    if (is.null(sv)) return(list(ok = NULL, u = NULL, lam = rep(NA_real_, npc)))
    lam_full <- rep(NA_real_, npc)
    d2 <- sv$d^2
    lam_full[seq_len(min(npc, length(d2)))] <- d2[seq_len(min(npc, length(d2)))]
    list(ok = ok, u = sv$u, lam = lam_full)
  }

  if (n_cores > 1L) {
    results <- mclapply(seq_len(n_win), one_window,
                        mc.cores = n_cores, mc.preschedule = TRUE)
  } else {
    results <- lapply(seq_len(n_win), one_window)
  }

  # Stitch results
  pcs <- vector("list", npc)
  for (k in seq_len(npc))
    pcs[[k]] <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                       dimnames = list(rownames(theta_mat), NULL))
  lambda <- matrix(NA_real_, nrow = n_win, ncol = npc)
  for (wi in seq_len(n_win)) {
    r <- results[[wi]]
    # mclapply may return try-error objects on worker crash — treat as NA
    if (inherits(r, "try-error") || !is.list(r)) next
    if (!is.null(r$lam)) lambda[wi, ] <- r$lam
    if (is.null(r$u) || is.null(r$ok)) next
    nu_keep <- min(npc, ncol(r$u))
    for (k in seq_len(nu_keep)) pcs[[k]][r$ok, wi] <- as.numeric(r$u[, k])
  }
  list(pcs = pcs, lambda = lambda)
}

# Sign-flip every PC against an anchor window's PC. PC1's sign drives the
# flip decision, but each PC is flipped independently against its own anchor
# (eigenvectors are orthogonal so signs are independent).
anchor_flip <- function(pcs, anchor_signal) {
  npc <- length(pcs); n_win <- ncol(pcs[[1]])
  anchor <- which.max(anchor_signal)
  if (length(anchor) == 0L || !is.finite(anchor_signal[anchor])) {
    med_z <- median(anchor_signal, na.rm = TRUE); if (!is.finite(med_z)) med_z <- 0
    anchor <- which.min(abs(anchor_signal - med_z))[1]
  }
  anchor <- as.integer(anchor)
  aligned <- vector("list", npc)
  for (k in seq_len(npc)) aligned[[k]] <- pcs[[k]]
  unflipped <- integer(0)
  for (w in seq_len(n_win)) {
    if (w == anchor) next
    a1 <- pcs[[1]][, anchor]; w1 <- pcs[[1]][, w]
    ok <- is.finite(a1) & is.finite(w1)
    if (sum(ok) < 10L) { unflipped <- c(unflipped, w); next }
    for (k in seq_len(npc)) {
      ak <- pcs[[k]][, anchor]; wk <- pcs[[k]][, w]
      okk <- is.finite(ak) & is.finite(wk)
      if (sum(okk) < 10L) next
      rk <- cor(ak[okk], wk[okk])
      if (is.finite(rk) && rk < 0) aligned[[k]][, w] <- -pcs[[k]][, w]
    }
  }
  list(aligned = aligned, anchor_idx = anchor, unflipped = unflipped)
}

# =============================================================================
# Per-chrom PCA
# =============================================================================

for (chrom in CHROM_LIST) {
  t_chr <- proc.time()
  loaded <- load_theta_matrix(chrom)
  if (is.null(loaded)) { message("[TR_B PCA] ", chrom, ": missing TSV — skip"); next }
  theta_mat <- loaded$theta_mat; n_sites_mat <- loaded$n_sites_mat
  win_grid  <- loaded$win_grid;  n_win <- nrow(win_grid)

  # Per-window cohort median + per-sample-deviation |Z|. The anchor for sign
  # flips is whichever window has the largest |Z|.
  window_median <- apply(theta_mat, 2, median, na.rm = TRUE)
  dev_mat <- sweep(theta_mat, 2, window_median, FUN = "-")
  all_devs <- as.vector(dev_mat); all_devs <- all_devs[is.finite(all_devs)]
  chrom_dev_mad <- if (length(all_devs) > 100) mad(all_devs, constant = 1.4826) else NA_real_
  if (!is.finite(chrom_dev_mad) || chrom_dev_mad <= 0) chrom_dev_mad <- 1e-6
  theta_z_direct <- apply(dev_mat, 2, function(v) {
    v <- v[is.finite(v)]; if (length(v) < 10) return(NA_real_)
    max(abs(v / chrom_dev_mad))
  })

  t_pca <- proc.time()
  pca <- local_pca_npc(theta_mat, n_sites_mat, PAD, NPC, N_CORES)
  message(sprintf("[TR_B PCA] %s: local_pca_npc (%d windows, %d cores) %.1fs",
                  chrom, n_win, N_CORES, (proc.time() - t_pca)[3]))
  flipped <- anchor_flip(pca$pcs, theta_z_direct)

  win_grid_dt <- data.table(
    chrom      = chrom,
    window_idx = win_grid$window_idx,
    start_bp   = win_grid$start_bp,
    end_bp     = win_grid$end_bp
  )

  out <- list(
    pcs               = flipped$aligned,
    lambda            = pca$lambda,
    theta_z_direct    = theta_z_direct,
    window_median     = window_median,
    win_grid          = win_grid_dt,
    anchor_idx        = flipped$anchor_idx,
    unflipped_windows = flipped$unflipped,
    chrom             = chrom,
    n_windows         = as.integer(n_win),
    n_samples         = as.integer(n_samp),
    npc               = as.integer(NPC),
    pad               = as.integer(PAD),
    sample_order      = sample_order
  )
  out_file <- file.path(OUT_LOCAL_PCA_DIR, sprintf("%s.window_pca.rds", chrom))
  saveRDS(out, out_file)
  message(sprintf("[TR_B PCA] %s: %d windows, anchor=%d, unflipped=%d (%.1fs) -> %s",
                  chrom, n_win, flipped$anchor_idx, length(flipped$unflipped),
                  (proc.time() - t_chr)[3], out_file))
}

message("[TR_B PCA] DONE")
