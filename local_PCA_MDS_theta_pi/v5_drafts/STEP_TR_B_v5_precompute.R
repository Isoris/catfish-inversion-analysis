#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_B_v5_precompute.R   (drafts/v5)
# =============================================================================
# Per-chromosome precompute for the theta-pi local-PCA + MDS path.
# Mirrors the z-blocks 03_precompute output shape so the verbatim ports
# of 04_detect_L1 / 05_plot_L1 / 06_detect_L2 / 07_plot_L2 (alongside
# this script as STEP_TR_C / TR_C_plot / TR_D / TR_D_plot) can run on
# theta-pi precomp without changes.
#
# Inputs (from STEP_TR_A; configured in 00_theta_config.sh):
#   $THETA_TSV_DIR/theta_native.<CHR>.<SCALE>.tsv.gz   per-chrom long TSV
#   $SAMPLE_LIST                                       cohort id list
#
# Outputs (under $OUTROOT, default = results_inversions/local_PCA_MDS_theta_pi):
#   precomp/<chr>.precomp.rds                  z-blocks-shaped precomp:
#     $dt   per-window data.table:
#       chrom, window_idx, start_bp, end_bp, mid_bp,
#       theta_pi_median,  theta_z_direct,                 (theta-pi specific)
#       MDS1..MDSk, MDS1_z..MDSk_z, max_abs_z, max_z_axis,
#       lambda_1..lambda_NPC, lambda_ratio,
#       anchor_window_idx,
#       PC_1_<sample>..PC_NPC_<sample>                    (per-sample loadings)
#     $sim_mat                window×window similarity (full or banded)
#     $sim_mat_format         "upper_triangle_float32" | "banded_float32_pmN"
#     $sim_band_half          if banded
#     $mds_mat                n_win × k_mds
#     $bg_continuity_quantiles
#     $chrom, $n_windows, $n_samples, $sample_order, $unflipped_windows
#   precomp/sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
#   window_dt.tsv.gz          genome-wide per-window scalar table
#   precomp_summary.tsv       per-chrom QC stats
#
# 04_detect_L1 expects --precomp_dir <OUTROOT>/precomp; this layout matches.
#
# Knobs:
#   --chrom <CHR>                 single chrom (default: iterate CHROM_LIST)
#   --pad <int>                   local-PCA neighbourhood half-width (default 1)
#   --npc <int>                   number of PCs to keep (default 4)
#   --kmds <int>                  number of MDS axes to keep (default 5)
#   --sim-band-half <int>         banded sim_mat half-width (default 200)
#   --sim-n-full-threshold <int>  switch full↔banded (default 6000)
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

PAD                  <- 1L
NPC                  <- 4L     # sim_mat itself only needs PC1; v5 keeps NPC=4 so
                               # PC_3_<sample>/PC_4_<sample> land in the precomp dt
                               # for atlas / carrier classification. (v4 = 2.)
                               # Override precedence: default < $LOCAL_PCA_NPC < --npc.
K_MDS                <- 20L    # MDS axes to compute/store (matches z-blocks 02a MDS_DIMS).
                               # Only the first SEED_MDS_AXES (= 5, hardcoded below) feed
                               # `max_abs_z`; the higher axes are kept for atlas / downstream.
SEED_MDS_AXES        <- 5L     # # of MDS axes whose z-scores contribute to max_abs_z
                               # (matches z-blocks 03_precompute SEED_MDS_AXES).
SIM_BAND_HALF        <- 200L
SIM_N_FULL_THRESHOLD <- 6000L
NN_SIM_SCALES        <- c(20, 40, 80, 120, 160, 200, 240, 320)

env_npc <- suppressWarnings(as.integer(Sys.getenv("LOCAL_PCA_NPC", unset = "")))
if (length(env_npc) == 1L && is.finite(env_npc) && env_npc > 0L) NPC <- env_npc

args <- commandArgs(trailingOnly = TRUE)
CHROM <- NULL
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--chrom"                && i < length(args)) { CHROM <- args[i + 1]; i <- i + 2L }
  else if (a == "--pad"                  && i < length(args)) { PAD <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--npc"                  && i < length(args)) { NPC <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--kmds"                 && i < length(args)) { K_MDS <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--sim-band-half"        && i < length(args)) { SIM_BAND_HALF <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--sim-n-full-threshold" && i < length(args)) { SIM_N_FULL_THRESHOLD <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}

THETA_TSV_DIR <- Sys.getenv("THETA_TSV_DIR", unset = NA)
PESTPG_SCALE  <- Sys.getenv("PESTPG_SCALE",  unset = "win10000.step2000")
OUTROOT       <- Sys.getenv("OUTROOT",       unset = NA)
SAMPLE_LIST   <- Sys.getenv("SAMPLE_LIST",   unset = NA)
stopifnot(!is.na(THETA_TSV_DIR), !is.na(OUTROOT), !is.na(SAMPLE_LIST))

if (is.null(CHROM)) {
  cl <- Sys.getenv("CHROM_LIST", unset = "")
  CHROM_LIST <- if (nchar(cl) > 0) strsplit(cl, "[ ,]+")[[1]] else
    sprintf("C_gar_LG%02d", 1:28)
} else {
  CHROM_LIST <- CHROM
}

# z-blocks-shaped output layout: precomp/<chr>.precomp.rds at top level,
# sim_mats/<chr>.sim_mat_nn{N}.rds in a sibling subfolder.
precomp_dir <- file.path(OUTROOT, "precomp")
sim_dir     <- file.path(precomp_dir, "sim_mats")
dir.create(precomp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sim_dir,     recursive = TRUE, showWarnings = FALSE)

sample_order <- readLines(SAMPLE_LIST)
sample_order <- sample_order[nchar(sample_order) > 0]
n_samp <- length(sample_order)
message("[TR_B v5] cohort: ", n_samp, " samples, NPC=", NPC, ", K_MDS=", K_MDS)
message("[TR_B v5] chroms: ", length(CHROM_LIST))

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

# Per-window heteroscedastic local PCA, NPC PCs.
local_pca_npc <- function(theta_mat, n_sites_mat, pad, npc) {
  n_samp <- nrow(theta_mat); n_win <- ncol(theta_mat)
  pcs    <- vector("list", npc)
  for (k in seq_len(npc))
    pcs[[k]] <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                       dimnames = list(rownames(theta_mat), NULL))
  lambda <- matrix(NA_real_, nrow = n_win, ncol = npc)

  ns_med <- median(n_sites_mat, na.rm = TRUE)
  if (!is.finite(ns_med) || ns_med <= 0) ns_med <- 1

  for (wi in seq_len(n_win)) {
    lo <- max(1L, wi - pad); hi <- min(n_win, wi + pad)
    block <- theta_mat[, lo:hi, drop = FALSE]
    ok <- complete.cases(block)
    if (sum(ok) < max(20L, ncol(block) + 2L)) next
    block_ok <- block[ok, , drop = FALSE]
    nf <- n_sites_mat[ok, wi]
    nf[!is.finite(nf) | nf <= 0] <- 1L
    w <- sqrt(nf / ns_med)
    bw <- block_ok * w
    centred <- sweep(bw, 2, colMeans(bw), FUN = "-")
    nu_req <- min(npc, ncol(centred))
    sv <- tryCatch(svd(centred, nu = nu_req, nv = 0), error = function(e) NULL)
    if (is.null(sv)) next
    for (k in seq_len(min(nu_req, ncol(sv$u))))
      pcs[[k]][ok, wi] <- as.numeric(sv$u[, k])
    d2 <- sv$d^2
    for (k in seq_len(min(npc, length(d2))))
      lambda[wi, k] <- d2[k]
  }
  list(pcs = pcs, lambda = lambda)
}

# Banded or full upper-triangle sim_mat from |cor(pc1[,i], pc1[,j])|.
build_sim_mat <- function(pc1_mat, n_win, band_half, n_full_threshold) {
  if (n_win <= n_full_threshold) {
    M <- matrix(NA_real_, n_win, n_win)
    for (i in seq_len(n_win)) {
      pi_ <- pc1_mat[, i]; if (all(!is.finite(pi_))) next
      for (j in seq.int(i, n_win)) {
        pj_ <- pc1_mat[, j]; ok <- is.finite(pi_) & is.finite(pj_)
        if (sum(ok) < 10L) next
        r <- cor(pi_[ok], pj_[ok])
        if (is.finite(r)) M[i, j] <- abs(r)
      }
    }
    diag(M) <- 1
    list(sim_mat = M, sim_band = NULL,
         format = "upper_triangle_float32", band_half = NA_integer_)
  } else {
    nb <- 2L * band_half + 1L
    B  <- matrix(NA_real_, n_win, nb)
    for (i in seq_len(n_win)) {
      pi_ <- pc1_mat[, i]; if (all(!is.finite(pi_))) next
      jl <- max(1L, i - band_half); jh <- min(n_win, i + band_half)
      for (j in jl:jh) {
        pj_ <- pc1_mat[, j]; ok <- is.finite(pi_) & is.finite(pj_)
        if (sum(ok) < 10L) next
        r <- cor(pi_[ok], pj_[ok])
        if (is.finite(r)) B[i, j - i + band_half + 1L] <- abs(r)
      }
    }
    B[, band_half + 1L] <- 1
    list(sim_mat = NULL, sim_band = B,
         format = paste0("banded_float32_pm", band_half),
         band_half = as.integer(band_half))
  }
}

sim_full <- function(simbox, n_win) {
  if (simbox$format == "upper_triangle_float32") {
    M <- simbox$sim_mat; M[lower.tri(M)] <- t(M)[lower.tri(M)]; M
  } else {
    band_med <- median(simbox$sim_band, na.rm = TRUE)
    if (!is.finite(band_med)) band_med <- 0.0
    M <- matrix(band_med, n_win, n_win); bh <- simbox$band_half
    for (i in seq_len(n_win)) {
      jl <- max(1L, i - bh); jh <- min(n_win, i + bh)
      for (j in jl:jh) {
        v <- simbox$sim_band[i, j - i + bh + 1L]
        if (is.finite(v)) M[i, j] <- v
      }
    }
    diag(M) <- 1; M
  }
}

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
    # Use PC1 to decide whether to flip THIS window — but flip each PC
    # independently against its own anchor (PC1 anchor for PC1 column,
    # PC2 anchor for PC2, ...). Eigenvectors are orthogonal so signs are
    # independent.
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

make_nn_sim <- function(mds_mat, dmat, k_use) {
  n <- nrow(mds_mat)
  smoothed <- matrix(0, n, ncol(mds_mat))
  for (wi in seq_len(n)) {
    d <- dmat[wi, ]; d[wi] <- Inf
    nn_idx <- order(d)[seq_len(k_use)]
    smoothed[wi, ] <- colMeans(mds_mat[c(wi, nn_idx), , drop = FALSE], na.rm = TRUE)
  }
  d_nn <- as.matrix(dist(smoothed))
  fv <- d_nn[is.finite(d_nn)]
  dmax <- if (length(fv) > 0) quantile(fv, 0.95, na.rm = TRUE) else 1
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  sim <- 1 - pmin(d_nn / dmax, 1)
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  sim
}

# =============================================================================
# Per-chrom precompute
# =============================================================================
all_window_dt <- list()
all_summary   <- list()

for (chrom in CHROM_LIST) {
  t_chr <- proc.time()
  loaded <- load_theta_matrix(chrom)
  if (is.null(loaded)) { message("[TR_B v5] ", chrom, ": missing TSV — skip"); next }
  theta_mat <- loaded$theta_mat; n_sites_mat <- loaded$n_sites_mat
  win_grid  <- loaded$win_grid; n_win <- nrow(win_grid)

  # Theta-direct |Z| (per-sample dev from window cohort median, max across
  # samples). Useful as a SECOND signal for the atlas, separate from the
  # MDS-axis |Z| z-blocks 04 expects. Stored as theta_z_direct.
  window_median <- apply(theta_mat, 2, median, na.rm = TRUE)
  dev_mat <- sweep(theta_mat, 2, window_median, FUN = "-")
  all_devs <- as.vector(dev_mat); all_devs <- all_devs[is.finite(all_devs)]
  chrom_dev_mad <- if (length(all_devs) > 100) mad(all_devs, constant = 1.4826) else NA_real_
  if (!is.finite(chrom_dev_mad) || chrom_dev_mad <= 0) chrom_dev_mad <- 1e-6
  theta_z_direct <- apply(dev_mat, 2, function(v) {
    v <- v[is.finite(v)]; if (length(v) < 10) return(NA_real_)
    max(abs(v / chrom_dev_mad))
  })

  # Per-window NPC-generalized local PCA.
  pca <- local_pca_npc(theta_mat, n_sites_mat, PAD, NPC)
  lambda_ratio <- ifelse(is.finite(pca$lambda[, 2]) & pca$lambda[, 2] > 0,
                         pca$lambda[, 1] / pca$lambda[, 2], NA_real_)

  # Sim_mat (banded or full UT) from PC1 correlations.
  simbox <- build_sim_mat(pca$pcs[[1]], n_win, SIM_BAND_HALF, SIM_N_FULL_THRESHOLD)
  message("[TR_B v5] ", chrom, ": sim_mat format = ", simbox$format)

  # Anchor-flip sign alignment, all NPC PCs, anchor on theta_z_direct.
  flipped <- anchor_flip(pca$pcs, theta_z_direct)
  pcs_aligned <- flipped$aligned

  # Full sim matrix transient → MDS + NN smoothing baseline.
  sim_M <- sim_full(simbox, n_win)
  dmat  <- 1 - sim_M; diag(dmat) <- 0
  mds_fit <- tryCatch(cmdscale(as.dist(dmat), k = K_MDS), error = function(e) NULL)
  if (is.null(mds_fit) || nrow(mds_fit) != n_win) {
    mds_mat <- matrix(NA_real_, n_win, K_MDS)
  } else {
    mds_mat <- mds_fit
    if (ncol(mds_mat) < K_MDS) {
      mds_mat <- cbind(mds_mat, matrix(NA_real_, n_win, K_MDS - ncol(mds_mat)))
    }
  }

  # MDS-axis robust z-scores. We compute K_MDS columns for parity with
  # z-blocks (default 20 axes stored in the precomp); but only the first
  # SEED_MDS_AXES (= 5, matches z-blocks 03_precompute) feed max_abs_z —
  # which is what 04_detect_L1 / 06_detect_L2 actually consume.
  mds_z <- matrix(NA_real_, nrow = n_win, ncol = K_MDS)
  for (k in seq_len(K_MDS)) {
    v <- mds_mat[, k]
    if (sum(is.finite(v)) >= 10) {
      med <- median(v, na.rm = TRUE); md  <- mad(v, na.rm = TRUE)
      if (is.finite(md) && md > 1e-12) {
        mds_z[, k] <- (v - med) / md
      } else {
        sdv <- sd(v, na.rm = TRUE)
        if (is.finite(sdv) && sdv > 0) mds_z[, k] <- (v - mean(v, na.rm = TRUE)) / sdv
      }
    }
  }
  z_for_max <- mds_z[, seq_len(min(SEED_MDS_AXES, K_MDS)), drop = FALSE]
  max_abs_z <- apply(z_for_max, 1, function(r) {
    r <- r[is.finite(r)]; if (length(r) == 0) NA_real_ else max(abs(r))
  })
  max_z_axis <- apply(z_for_max, 1, function(r) {
    a <- which.max(abs(r)); if (length(a) == 0) NA_integer_ else as.integer(a)
  })

  # bg_continuity quantiles (off-diagonal sim distribution near baseline).
  bg_q <- if (n_win >= 10) {
    adj <- vapply(seq_len(n_win - 1), function(i) sim_M[i, i + 1L], numeric(1))
    quantile(adj, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
  } else rep(NA_real_, 6)

  # Save sim_mat_nn0 (=raw sim) + NN-smoothed sim_mats at each scale.
  saveRDS(sim_M, file.path(sim_dir, sprintf("%s.sim_mat_nn0.rds", chrom)))
  for (k in NN_SIM_SCALES) {
    k_use <- min(k, n_win - 1L)
    if (k_use < 2L) next
    t_nn <- proc.time()
    nn_sim <- make_nn_sim(mds_mat, dmat, k_use)
    saveRDS(nn_sim, file.path(sim_dir, sprintf("%s.sim_mat_nn%d.rds", chrom, k)))
    message(sprintf("[TR_B v5] %s: sim_mat_nn%d (%.1fs)",
                    chrom, k, (proc.time() - t_nn)[3]))
  }
  rm(sim_M); gc(verbose = FALSE)

  # Per-window data.table — z-blocks-shaped. Per-sample PC_k_<sample> columns
  # for k = 1..NPC. Atlas reads PC_1_*/PC_2_* (existing); NPC > 2 makes the
  # higher-PC columns available too.
  dt <- data.table(
    chrom            = chrom,
    window_idx       = win_grid$window_idx,
    start_bp         = win_grid$start_bp,
    end_bp           = win_grid$end_bp,
    mid_bp           = as.integer((win_grid$start_bp + win_grid$end_bp) / 2L),
    theta_pi_median  = round(window_median, 6),
    theta_z_direct   = round(theta_z_direct, 4),
    max_abs_z        = round(max_abs_z, 4),                          # MDS-axis (z-blocks)
    max_z_axis       = max_z_axis,
    lambda_ratio     = round(lambda_ratio, 4),
    anchor_window_idx = flipped$anchor_idx - 1L                       # 0-indexed
  )
  for (k in seq_len(NPC))     set(dt, j = paste0("lambda_", k), value = round(pca$lambda[, k], 6))
  for (k in seq_len(K_MDS))   set(dt, j = paste0("MDS",    k), value = round(mds_mat[, k], 6))
  for (k in seq_len(K_MDS))   set(dt, j = paste0("MDS",    k, "_z"), value = round(mds_z[, k], 4))
  for (k in seq_len(NPC)) {
    pck <- pcs_aligned[[k]]
    for (s_idx in seq_len(n_samp)) {
      sid <- sample_order[s_idx]
      set(dt, j = paste0("PC_", k, "_", sid), value = round(pck[s_idx, ], 6))
    }
  }

  precomp <- list(
    dt                       = dt,
    sim_mat                  = simbox$sim_mat,
    sim_band                 = simbox$sim_band,
    sim_mat_format           = simbox$format,
    sim_band_half            = simbox$band_half,
    mds_mat                  = mds_mat,
    bg_continuity_quantiles  = bg_q,
    chrom                    = chrom,
    n_windows                = as.integer(n_win),
    n_samples                = as.integer(n_samp),
    npc                      = as.integer(NPC),
    k_mds                    = as.integer(K_MDS),
    sample_order             = sample_order,
    unflipped_windows        = flipped$unflipped - 1L
  )
  saveRDS(precomp, file.path(precomp_dir, sprintf("%s.precomp.rds", chrom)))
  message("[TR_B v5] ", chrom, ": precomp.rds saved (", n_win, " windows, NPC=", NPC, ")")

  # Genome-wide rollup row + per-chrom QC summary.
  all_window_dt[[chrom]] <- dt[, .(chrom, window_idx, start_bp, end_bp, mid_bp,
                                    theta_pi_median, theta_z_direct,
                                    max_abs_z, max_z_axis, lambda_ratio,
                                    MDS1, MDS2)]
  elapsed <- round((proc.time() - t_chr)[3], 1)
  all_summary[[chrom]] <- data.table(
    chrom            = chrom,
    n_windows        = as.integer(n_win),
    n_samples        = as.integer(n_samp),
    npc              = as.integer(NPC),
    k_mds            = as.integer(K_MDS),
    sim_format       = simbox$format,
    n_unflipped      = length(flipped$unflipped),
    median_max_z     = round(median(max_abs_z, na.rm = TRUE), 3),
    q95_max_z        = round(quantile(max_abs_z, 0.95, na.rm = TRUE), 3),
    median_theta_z   = round(median(theta_z_direct, na.rm = TRUE), 3),
    bg_q50           = round(bg_q["50%"], 4),
    bg_q90           = round(bg_q["90%"], 4),
    bg_q95           = round(bg_q["95%"], 4),
    elapsed_sec      = elapsed
  )
  message("[TR_B v5] ", chrom, ": done in ", elapsed, "s")
}

# =============================================================================
# Genome-wide outputs
# =============================================================================
if (length(all_window_dt) > 0) {
  fwrite(rbindlist(all_window_dt, fill = TRUE),
         file.path(OUTROOT, "window_dt.tsv.gz"), sep = "\t", compress = "gzip")
  message("[TR_B v5] wrote ", file.path(OUTROOT, "window_dt.tsv.gz"))
}
if (length(all_summary) > 0) {
  fwrite(rbindlist(all_summary, fill = TRUE),
         file.path(OUTROOT, "precomp_summary.tsv"), sep = "\t")
  message("[TR_B v5] wrote ", file.path(OUTROOT, "precomp_summary.tsv"))
}

message("[TR_B v5] DONE")