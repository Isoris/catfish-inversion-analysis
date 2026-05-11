#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_C_mds_compute.R
# =============================================================================
# MDS + sim_mat + per-window features. Reads the per-chrom local-PCA bundle
# from STEP_TR_B and emits the precomp.rds + (optionally) NN-smoothed sim_mats
# that downstream L1/L2 detection (TR_D…G) consumes.
#
# Two modes — pick the right one for the input scale:
#
#   --mode full   (default; intended for COARSE scale, e.g. win50000.step10000)
#       * Builds banded sim_mat, then reconstructs the FULL N×N sim_mat in
#         memory.
#       * Runs full cmdscale MDS to k = K_MDS axes.
#       * Writes precomp.rds with mds_mat populated + full per-window MDS-z
#         columns + max_abs_z (consumed by TR_D detect_L1).
#       * Writes sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
#         (consumed by TR_D/E/F/G).
#       * Memory cost: ~N²·8B during sim_full + cmdscale. For LG28 at
#         win10000.step2000 (~16,500 windows) this is ~10–15 GB peak — DON'T
#         use full mode at the dense scale on a laptop.
#
#   --mode local  (intended for DENSE scale, e.g. win10000.step2000)
#       * Builds banded sim_mat only — no sim_full reconstruction.
#       * Skips cmdscale MDS (mds_mat / MDS{i} / MDS{i}_z columns are NA).
#       * Writes a "dense" precomp.rds with dt + sim_band + per-sample PC
#         loadings + theta_z_direct + lambda + morphology features.
#       * Does NOT write sim_mats/<chr>.sim_mat_nn{N}.rds (those require
#         full sim_mat reconstruction). TR_D/E/F/G cannot run on a local
#         precomp by design — it's intended as the input to TR_H carrier
#         classification, TR_I CUSUM refinement, and any future
#         boundary-refinement step that reads coarse L2 envelopes and
#         refines them locally using dense per-window PC scores.
#       * Memory cost: linear in N (banded only). LG28 dense fits in <2 GB.
#
# Inputs (configured via 00_theta_config.sh):
#   $OUT_LOCAL_PCA_DIR/<chr>.window_pca.rds   (from STEP_TR_B)
#
# Output (under $OUTROOT/<--out_subdir>, default subdir = "precomp"):
#   <chr>.precomp.rds          per-window dt + (in full mode) mds_mat + sim_mat
#   sim_mats/<chr>.sim_mat_nn{N}.rds   (full mode only)
#   ../window_dt.tsv.gz        genome-wide rollup (when iterating CHROM_LIST)
#   ../precomp_summary.tsv     per-chrom QC stats
#
# Knobs:
#   --chrom <CHR>                 single chrom (default: iterate CHROM_LIST)
#   --mode <full|local>           default "full"
#   --kmds <int>                  MDS axes (default 5; only used in full mode)
#   --sim-band-half <int>         banded sim_mat half-width in WINDOWS (default 200)
#   --sim-n-full-threshold <int>  switch full↔banded storage (default 6000)
#                                 In full mode, n_win > threshold ⇒ banded
#                                 storage but still full reconstruction.
#                                 In local mode, this knob is ignored — always
#                                 banded.
#   --out_subdir <name>           output subfolder under $OUTROOT (default "precomp")
#                                 Use different names for coarse/dense to keep
#                                 outputs isolated, e.g. "precomp_coarse" /
#                                 "precomp_dense".
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

K_MDS                <- 5L      # was 20 in original v5; SEED_MDS_AXES has always been 5
SEED_MDS_AXES        <- 5L      # # of MDS axes whose z-scores feed max_abs_z
SIM_BAND_HALF        <- 200L
SIM_N_FULL_THRESHOLD <- 6000L
NN_SIM_SCALES        <- c(20, 40, 80, 120, 160, 200, 240, 320)
MODE                 <- "full"
OUT_SUBDIR           <- "precomp"

args <- commandArgs(trailingOnly = TRUE)
CHROM <- NULL
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--chrom"                && i < length(args)) { CHROM <- args[i + 1]; i <- i + 2L }
  else if (a == "--mode"                 && i < length(args)) { MODE <- args[i + 1]; i <- i + 2L }
  else if (a == "--kmds"                 && i < length(args)) { K_MDS <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--sim-band-half"        && i < length(args)) { SIM_BAND_HALF <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--sim-n-full-threshold" && i < length(args)) { SIM_N_FULL_THRESHOLD <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--out_subdir"           && i < length(args)) { OUT_SUBDIR <- args[i + 1]; i <- i + 2L }
  else { i <- i + 1L }
}
stopifnot(MODE %in% c("full", "local"))

OUTROOT       <- Sys.getenv("OUTROOT", unset = NA)
OUT_LOCAL_PCA_DIR <- Sys.getenv("OUT_LOCAL_PCA_DIR",
                                unset = file.path(OUTROOT, "01_local_pca"))
stopifnot(!is.na(OUTROOT))

precomp_dir <- file.path(OUTROOT, OUT_SUBDIR)
sim_dir     <- file.path(precomp_dir, "sim_mats")
dir.create(precomp_dir, recursive = TRUE, showWarnings = FALSE)
if (MODE == "full") dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)

if (is.null(CHROM)) {
  cl <- Sys.getenv("CHROM_LIST", unset = "")
  CHROM_LIST <- if (nchar(cl) > 0) strsplit(cl, "[ ,]+")[[1]] else
    sprintf("C_gar_LG%02d", 1:28)
} else {
  CHROM_LIST <- CHROM
}

message(sprintf("[TR_C MDS] mode=%s, K_MDS=%d, SIM_BAND_HALF=%d, out=%s",
                MODE, K_MDS, SIM_BAND_HALF, precomp_dir))
message("[TR_C MDS] chroms: ", length(CHROM_LIST))

# =============================================================================
# Helpers
# =============================================================================

# Banded or full upper-triangle sim_mat from |cor(pc1[,i], pc1[,j])|.
# In local mode we always emit banded regardless of n_win.
build_sim_mat <- function(pc1_mat, n_win, band_half, n_full_threshold, force_banded = FALSE) {
  if (!force_banded && n_win <= n_full_threshold) {
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

# Reconstruct full N×N sim_mat from banded form (off-band cells = band median).
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
# Per-chrom MDS + features
# =============================================================================
all_window_dt <- list()
all_summary   <- list()

for (chrom in CHROM_LIST) {
  t_chr <- proc.time()
  pca_file <- file.path(OUT_LOCAL_PCA_DIR, sprintf("%s.window_pca.rds", chrom))
  if (!file.exists(pca_file)) {
    message("[TR_C MDS] ", chrom, ": missing PCA RDS (", pca_file, ") — skip")
    next
  }
  pca <- readRDS(pca_file)
  n_win   <- pca$n_windows
  n_samp  <- pca$n_samples
  NPC     <- pca$npc
  win_grid <- pca$win_grid
  pcs_aligned    <- pca$pcs            # already anchor-flipped by TR_B
  lambda         <- pca$lambda
  theta_z_direct <- pca$theta_z_direct
  window_median  <- pca$window_median
  sample_order   <- pca$sample_order
  unflipped      <- pca$unflipped_windows
  anchor_idx     <- pca$anchor_idx

  lambda_ratio <- ifelse(is.finite(lambda[, 2]) & lambda[, 2] > 0,
                         lambda[, 1] / lambda[, 2], NA_real_)

  # Build sim_mat. Full mode: full or banded depending on threshold.
  # Local mode: always banded (force_banded = TRUE).
  simbox <- build_sim_mat(pcs_aligned[[1]], n_win, SIM_BAND_HALF,
                          SIM_N_FULL_THRESHOLD, force_banded = (MODE == "local"))
  message("[TR_C MDS] ", chrom, ": sim_mat format = ", simbox$format)

  if (MODE == "full") {
    # Reconstruct full N×N matrix → cmdscale → MDS axes z-scores → NN sim_mats.
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

    bg_q <- if (n_win >= 10) {
      adj <- vapply(seq_len(n_win - 1), function(i) sim_M[i, i + 1L], numeric(1))
      quantile(adj, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
    } else rep(NA_real_, 6)

    saveRDS(sim_M, file.path(sim_dir, sprintf("%s.sim_mat_nn0.rds", chrom)))
    for (k in NN_SIM_SCALES) {
      k_use <- min(k, n_win - 1L)
      if (k_use < 2L) next
      t_nn <- proc.time()
      nn_sim <- make_nn_sim(mds_mat, dmat, k_use)
      saveRDS(nn_sim, file.path(sim_dir, sprintf("%s.sim_mat_nn%d.rds", chrom, k)))
      message(sprintf("[TR_C MDS] %s: sim_mat_nn%d (%.1fs)",
                      chrom, k, (proc.time() - t_nn)[3]))
    }
    rm(sim_M, dmat); gc(verbose = FALSE)
  } else {
    # Local mode: no full reconstruction, no MDS, no NN sim_mats.
    mds_mat   <- matrix(NA_real_, n_win, K_MDS)
    mds_z     <- matrix(NA_real_, n_win, K_MDS)
    max_abs_z <- rep(NA_real_, n_win)
    max_z_axis <- rep(NA_integer_, n_win)
    bg_q      <- rep(NA_real_, 6)
  }

  # Per-window data.table.
  dt <- data.table(
    chrom            = chrom,
    window_idx       = win_grid$window_idx,
    start_bp         = win_grid$start_bp,
    end_bp           = win_grid$end_bp,
    mid_bp           = as.integer((win_grid$start_bp + win_grid$end_bp) / 2L),
    theta_pi_median  = round(window_median, 6),
    theta_z_direct   = round(theta_z_direct, 4),
    max_abs_z        = round(max_abs_z, 4),
    max_z_axis       = max_z_axis,
    lambda_ratio     = round(lambda_ratio, 4),
    anchor_window_idx = anchor_idx - 1L
  )
  for (k in seq_len(NPC))    set(dt, j = paste0("lambda_", k), value = round(lambda[, k], 6))
  for (k in seq_len(K_MDS))  set(dt, j = paste0("MDS",    k),       value = round(mds_mat[, k], 6))
  for (k in seq_len(K_MDS))  set(dt, j = paste0("MDS",    k, "_z"), value = round(mds_z[, k],   4))
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
    unflipped_windows        = unflipped - 1L,
    mode                     = MODE
  )
  saveRDS(precomp, file.path(precomp_dir, sprintf("%s.precomp.rds", chrom)))
  message(sprintf("[TR_C MDS] %s: precomp.rds saved (%d windows, NPC=%d, K_MDS=%d, mode=%s)",
                  chrom, n_win, NPC, K_MDS, MODE))

  all_window_dt[[chrom]] <- dt[, .(chrom, window_idx, start_bp, end_bp, mid_bp,
                                    theta_pi_median, theta_z_direct,
                                    max_abs_z, max_z_axis, lambda_ratio,
                                    MDS1, MDS2)]
  elapsed <- round((proc.time() - t_chr)[3], 1)
  all_summary[[chrom]] <- data.table(
    chrom            = chrom,
    mode             = MODE,
    n_windows        = as.integer(n_win),
    n_samples        = as.integer(n_samp),
    npc              = as.integer(NPC),
    k_mds            = as.integer(K_MDS),
    sim_format       = simbox$format,
    n_unflipped      = length(unflipped),
    median_max_z     = round(median(max_abs_z, na.rm = TRUE), 3),
    q95_max_z        = round(quantile(max_abs_z, 0.95, na.rm = TRUE), 3),
    median_theta_z   = round(median(theta_z_direct, na.rm = TRUE), 3),
    bg_q50           = round(bg_q[1], 4),
    bg_q90           = round(bg_q[5], 4),
    bg_q95           = round(bg_q[6], 4),
    elapsed_sec      = elapsed
  )
  message("[TR_C MDS] ", chrom, ": done in ", elapsed, "s")
}

# Genome-wide rollups (one per --out_subdir).
if (length(all_window_dt) > 0) {
  fwrite(rbindlist(all_window_dt, fill = TRUE),
         file.path(precomp_dir, "window_dt.tsv.gz"), sep = "\t", compress = "gzip")
  message("[TR_C MDS] wrote ", file.path(precomp_dir, "window_dt.tsv.gz"))
}
if (length(all_summary) > 0) {
  fwrite(rbindlist(all_summary, fill = TRUE),
         file.path(precomp_dir, "precomp_summary.tsv"), sep = "\t")
  message("[TR_C MDS] wrote ", file.path(precomp_dir, "precomp_summary.tsv"))
}

message("[TR_C MDS] DONE")
