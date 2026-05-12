#!/usr/bin/env Rscript
# =============================================================================
# STEP_GH_B_local_pca_compute.R  (v1, 2026-05-12)
# =============================================================================
# Per-window heteroscedastic local PCA on the GHSL divergence matrix
# produced by STEP_GH_A v2. Output schema mirrors θπ STEP_TR_B v6 so the
# downstream MDS / sketch / sparse-graph step (GH_C) can be a near-direct
# port of TR_C v6.
#
# Feature: rolling[[--scale]] from the ghsl_matrices.rds.
#   --scale s50   (coarse, default for --mode full)
#   --scale s10   (dense, recommended for --mode local)
#   --scale s100  (overview)
#   --scale raw   (raw 5 kb base from div_mat; noisy, NOT recommended for
#                  the local PCA layer — kept for diagnostic runs only)
#
# Heteroscedastic weight: for each sample, the SVD input is multiplied by
#   w_s = sqrt(n_total_s / median(n_total_*))
# where n_total_s is the aggregated denominator at the chosen scale for
# the focal window. This down-weights samples with sparse evidence in
# the current window (their divergence ratio is noisier).
#
# Output (one file per chrom):
#   <chr>.ghsl_window_pca.rds (schema identical in spirit to θπ TR_B's
#   window_pca.rds, with field names renamed for clarity):
#     $pcs               — list of NPC matrices [n_samp × n_win], anchor-flipped
#     $lambda            — [n_win × NPC] eigenvalues from each window's SVD
#     $div_z_direct      — [n_win] per-window |Z| from per-sample div deviation
#                          (analog of theta_z_direct)
#     $window_median     — [n_win] per-window cohort median divergence
#     $win_grid          — data.table: chrom, window_idx, start_bp, end_bp
#     $anchor_idx        — 1-based anchor window index
#     $unflipped_windows — 1-based indices that couldn't be sign-aligned
#     $chrom, $n_windows, $n_samples, $npc, $pad, $sample_order, $scale
#
# Knobs:
#   --chrom <CHR>          single chrom (default: iterate all *.ghsl_matrices.rds)
#   --matrices_dir <dir>   input directory (output of GH_A v2)
#   --outdir <dir>         output directory for *.ghsl_window_pca.rds
#   --scale <s10|s50|...>  which rolling scale to use (default: s50)
#   --pad <int>            local-PCA neighbourhood half-width (default 1)
#   --npc <int>            number of PCs to keep (default 4)
#   --cores <int>          mclapply mc.cores (default $N_CORES or 8)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# ── Defaults ────────────────────────────────────────────────────────────────
MATRICES_DIR <- NA_character_
OUTDIR       <- NA_character_
CHROM        <- NULL
SCALE        <- "s50"
PAD          <- 1L
NPC          <- 4L
N_CORES      <- as.integer(Sys.getenv("N_CORES", unset = "8"))
if (!is.finite(N_CORES) || N_CORES < 1L) N_CORES <- 1L

args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--matrices_dir" && i < length(args)) { MATRICES_DIR <- args[i + 1]; i <- i + 2L }
  else if (a == "--outdir"       && i < length(args)) { OUTDIR       <- args[i + 1]; i <- i + 2L }
  else if (a == "--chrom"        && i < length(args)) { CHROM        <- args[i + 1]; i <- i + 2L }
  else if (a == "--scale"        && i < length(args)) { SCALE        <- args[i + 1]; i <- i + 2L }
  else if (a == "--pad"          && i < length(args)) { PAD     <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--npc"          && i < length(args)) { NPC     <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--cores"        && i < length(args)) { N_CORES <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}
stopifnot(!is.na(MATRICES_DIR), !is.na(OUTDIR))
if (.Platform$OS.type == "windows" && N_CORES > 1L) N_CORES <- 1L
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[GH_B] matrices_dir=", MATRICES_DIR)
message("[GH_B] outdir=", OUTDIR)
message("[GH_B] scale=", SCALE, ", PAD=", PAD, ", NPC=", NPC, ", cores=", N_CORES)

# ── Chrom list ──────────────────────────────────────────────────────────────
mfiles <- sort(list.files(MATRICES_DIR, pattern = "\\.ghsl_matrices\\.rds$",
                          full.names = TRUE))
if (length(mfiles) == 0L) stop("[GH_B] no *.ghsl_matrices.rds in ", MATRICES_DIR)
chrom_from_file <- sub("\\.ghsl_matrices\\.rds$", "", basename(mfiles))
names(mfiles) <- chrom_from_file
chroms <- if (is.null(CHROM)) chrom_from_file else intersect(chrom_from_file, CHROM)
if (length(chroms) == 0L) stop("[GH_B] chrom ", CHROM, " not in ", MATRICES_DIR)
message("[GH_B] chroms: ", length(chroms))

# =============================================================================
# Helpers
# =============================================================================

# Pull the feature matrix and its denominator for the requested scale.
# Returns (feat_mat, ntot_mat) both [n_samp × n_win] aligned with sample_names.
get_feature_for_scale <- function(gm, scale_label) {
  if (scale_label == "raw") {
    list(feat = gm$div_mat, ntot = gm$n_total_mat)
  } else {
    feat <- gm$rolling[[scale_label]]
    ntot <- gm$rolling_n_total[[scale_label]]
    if (is.null(feat) || is.null(ntot))
      stop("[GH_B] scale '", scale_label, "' not in ghsl_matrices.rds (have: ",
           paste(names(gm$rolling), collapse = ","), ")")
    list(feat = feat, ntot = ntot)
  }
}

# Per-window heteroscedastic local PCA. Mirrors θπ TR_B v6 local_pca_npc
# with the GHSL-appropriate weight (sqrt(n_total / median)). Parallelized
# via mclapply.
local_pca_npc <- function(feat_mat, ntot_mat, pad, npc, n_cores) {
  n_samp <- nrow(feat_mat); n_win <- ncol(feat_mat)
  ntot_med <- median(ntot_mat[is.finite(ntot_mat) & ntot_mat > 0],
                     na.rm = TRUE)
  if (!is.finite(ntot_med) || ntot_med <= 0) ntot_med <- 1

  one_window <- function(wi) {
    lo <- max(1L, wi - pad); hi <- min(n_win, wi + pad)
    block <- feat_mat[, lo:hi, drop = FALSE]
    ok <- complete.cases(block)
    if (sum(ok) < max(20L, ncol(block) + 2L)) {
      return(list(ok = NULL, u = NULL, lam = rep(NA_real_, npc)))
    }
    block_ok <- block[ok, , drop = FALSE]
    # Hetero weight uses the FOCAL window's n_total per sample
    nf <- ntot_mat[ok, wi]
    nf[!is.finite(nf) | nf <= 0] <- 1L
    w <- sqrt(nf / ntot_med)
    bw <- block_ok * w
    centred <- sweep(bw, 2L, colMeans(bw), FUN = "-")
    nu_req <- min(npc, ncol(centred))
    sv <- tryCatch(svd(centred, nu = nu_req, nv = 0), error = function(e) NULL)
    if (is.null(sv)) return(list(ok = NULL, u = NULL, lam = rep(NA_real_, npc)))
    lam_full <- rep(NA_real_, npc)
    d2 <- sv$d^2
    lam_full[seq_len(min(npc, length(d2)))] <- d2[seq_len(min(npc, length(d2)))]
    list(ok = ok, u = sv$u, lam = lam_full)
  }

  results <- if (n_cores > 1L) {
    mclapply(seq_len(n_win), one_window, mc.cores = n_cores, mc.preschedule = TRUE)
  } else {
    lapply(seq_len(n_win), one_window)
  }

  pcs <- vector("list", npc)
  for (k in seq_len(npc))
    pcs[[k]] <- matrix(NA_real_, nrow = n_samp, ncol = n_win)
  lambda <- matrix(NA_real_, nrow = n_win, ncol = npc)
  for (wi in seq_len(n_win)) {
    r <- results[[wi]]
    if (inherits(r, "try-error") || !is.list(r)) next
    if (!is.null(r$lam)) lambda[wi, ] <- r$lam
    if (is.null(r$u) || is.null(r$ok)) next
    nu_keep <- min(npc, ncol(r$u))
    for (k in seq_len(nu_keep)) pcs[[k]][r$ok, wi] <- as.numeric(r$u[, k])
  }
  list(pcs = pcs, lambda = lambda)
}

# Anchor-flip — identical to θπ TR_B.
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

for (chr in chroms) {
  t_chr <- proc.time()
  gm <- readRDS(mfiles[[chr]])
  feat <- get_feature_for_scale(gm, SCALE)
  feat_mat <- feat$feat; ntot_mat <- feat$ntot
  sample_order <- gm$sample_names
  n_samp <- length(sample_order); n_win <- ncol(feat_mat)
  rownames(feat_mat) <- sample_order
  rownames(ntot_mat) <- sample_order

  # Cohort median + per-sample-deviation |Z| (anchor signal)
  window_median <- apply(feat_mat, 2L, median, na.rm = TRUE)
  dev_mat <- sweep(feat_mat, 2L, window_median, FUN = "-")
  all_devs <- as.vector(dev_mat); all_devs <- all_devs[is.finite(all_devs)]
  chrom_dev_mad <- if (length(all_devs) > 100L) mad(all_devs, constant = 1.4826) else NA_real_
  if (!is.finite(chrom_dev_mad) || chrom_dev_mad <= 0) chrom_dev_mad <- 1e-6
  div_z_direct <- apply(dev_mat, 2L, function(v) {
    v <- v[is.finite(v)]; if (length(v) < 10L) return(NA_real_)
    max(abs(v / chrom_dev_mad))
  })

  t_pca <- proc.time()
  pca <- local_pca_npc(feat_mat, ntot_mat, PAD, NPC, N_CORES)
  message(sprintf("[GH_B] %s: local PCA (%d × %d, scale=%s, cores=%d) %.1fs",
                  chr, n_samp, n_win, SCALE, N_CORES,
                  (proc.time() - t_pca)[3]))
  flipped <- anchor_flip(pca$pcs, div_z_direct)

  win_grid_dt <- data.table(
    chrom      = chr,
    window_idx = gm$window_info$window_idx,
    start_bp   = gm$window_info$start_bp,
    end_bp     = gm$window_info$end_bp
  )

  out <- list(
    pcs               = flipped$aligned,
    lambda            = pca$lambda,
    div_z_direct      = div_z_direct,
    window_median     = window_median,
    win_grid          = win_grid_dt,
    anchor_idx        = flipped$anchor_idx,
    unflipped_windows = flipped$unflipped,
    chrom             = chr,
    n_windows         = as.integer(n_win),
    n_samples         = as.integer(n_samp),
    npc               = as.integer(NPC),
    pad               = as.integer(PAD),
    sample_order      = sample_order,
    scale             = SCALE
  )
  out_file <- file.path(OUTDIR, paste0(chr, ".ghsl_window_pca.rds"))
  saveRDS(out, out_file)
  message(sprintf("[GH_B] %s: %d windows, anchor=%d, unflipped=%d (%.1fs) -> %s",
                  chr, n_win, flipped$anchor_idx, length(flipped$unflipped),
                  (proc.time() - t_chr)[3], out_file))
  rm(gm, feat, feat_mat, ntot_mat, pca, flipped); invisible(gc(verbose = FALSE))
}

message("[GH_B] DONE")
