#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_B_v5_precompute.R   (drafts/v5)
# =============================================================================
# Replaces v4's all-in-one TR_B with a precompute step that mirrors the
# z-blocks pipeline shape (local_PCA_MDS_z/03_precompute_localpca_zblocks.R).
# Writes one precomp.rds per chromosome plus NN-smoothed sim_mats + a
# genome-wide window_dt rollup. NO envelope detection, NO atlas JSON —
# those are TR_C, TR_D, TR_G in this draft.
#
# Inputs (from STEP_TR_A; configured in 00_theta_config.sh):
#   $THETA_TSV_DIR/theta_native.<CHR>.<SCALE>.tsv.gz   per-chrom long TSV
#   $SAMPLE_LIST                                       226-id list
#
# Outputs (under $OUTROOT, default = results_inversions/local_PCA_MDS_theta_pi):
#   03_per_chrom/<chr>/precomp.rds       full precomp object:
#     $dt        per-window data.table with start_bp, end_bp, mid_bp,
#                window_idx, theta_pi_median, max_abs_z, top10_abs_z,
#                lambda_1, lambda_2, lambda_ratio, anchor_window_idx,
#                MDS1, MDS2, plus per-sample PC_1_<sample>, PC_2_<sample>
#     $sim_mat        full upper-tri OR banded ±SIM_BAND_HALF (see format)
#     $sim_mat_format "upper_triangle_float32" or "banded_float32_pmN"
#     $sim_band_half  if banded
#     $mds_mat        n_win x 2 (MDS1, MDS2)
#     $bg_continuity_quantiles
#     $chrom, $n_windows, $n_samples
#   03_per_chrom/<chr>/sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
#                                     NN-smoothed sim_mats (MDS-space k-NN)
#   window_dt.tsv.gz   genome-wide per-window scalar table (all chroms)
#   precomp_summary.tsv  per-chrom QC stats
#
# Usage:
#   source 00_theta_config.sh
#   Rscript STEP_TR_B_v5_precompute.R --chrom <CHR>           # single chrom
#   Rscript STEP_TR_B_v5_precompute.R                         # all chroms in CHROM_LIST
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

# -----------------------------------------------------------------------------
# Config knobs (overridable via --flag)
# -----------------------------------------------------------------------------
PAD                  <- 1L         # local-PCA neighbourhood half-width
SIM_BAND_HALF        <- 200L       # banded sim_mat half-width (windows)
SIM_N_FULL_THRESHOLD <- 6000L      # n_win > threshold => banded; else full UT
NN_SIM_SCALES        <- c(20, 40, 80, 120, 160, 200, 240, 320)

args <- commandArgs(trailingOnly = TRUE)
CHROM <- NULL
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--chrom"                && i < length(args)) { CHROM <- args[i + 1]; i <- i + 2L }
  else if (a == "--pad"                  && i < length(args)) { PAD <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--sim-band-half"        && i < length(args)) { SIM_BAND_HALF <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--sim-n-full-threshold" && i < length(args)) { SIM_N_FULL_THRESHOLD <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}

# Required env vars from 00_theta_config.sh
THETA_TSV_DIR <- Sys.getenv("THETA_TSV_DIR", unset = NA)
PESTPG_SCALE  <- Sys.getenv("PESTPG_SCALE",  unset = "win10000.step2000")
OUTROOT       <- Sys.getenv("OUTROOT",       unset = NA)
SAMPLE_LIST   <- Sys.getenv("SAMPLE_LIST",   unset = NA)
stopifnot(!is.na(THETA_TSV_DIR), !is.na(OUTROOT), !is.na(SAMPLE_LIST))

# Determine chrom set
if (is.null(CHROM)) {
  cl <- Sys.getenv("CHROM_LIST", unset = "")
  CHROM_LIST <- if (nchar(cl) > 0) strsplit(cl, "[ ,]+")[[1]] else
    sprintf("C_gar_LG%02d", 1:28)
} else {
  CHROM_LIST <- CHROM
}

per_chrom_dir <- file.path(OUTROOT, "03_per_chrom")
dir.create(per_chrom_dir, recursive = TRUE, showWarnings = FALSE)

sample_order <- readLines(SAMPLE_LIST)
sample_order <- sample_order[nchar(sample_order) > 0]
n_samp <- length(sample_order)
message("[TR_B v5] cohort: ", n_samp, " samples")
message("[TR_B v5] chroms: ", length(CHROM_LIST), " (", paste(head(CHROM_LIST, 3), collapse = ","),
        if (length(CHROM_LIST) > 3) ",..." else "", ")")

# =============================================================================
# Per-chrom helpers
# =============================================================================

# Build samples × windows matrix from the long TSV produced by TR_A.
load_theta_matrix <- function(chrom) {
  tsv <- file.path(THETA_TSV_DIR, sprintf("theta_native.%s.%s.tsv.gz", chrom, PESTPG_SCALE))
  if (!file.exists(tsv)) {
    message("[TR_B v5] ", chrom, ": missing TSV — skip (", tsv, ")")
    return(NULL)
  }
  long_dt <- fread(tsv)
  long_dt <- long_dt[chrom == ..chrom]
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
  fill <- cbind(rows[good], cols[good])
  theta_mat[fill]   <- long_dt$theta_pi[good]
  n_sites_mat[fill] <- long_dt$n_sites[good]
  list(theta_mat = theta_mat, n_sites_mat = n_sites_mat, win_grid = win_grid)
}

# Per-window heteroscedastic-weighted local PCA (carried over from v4).
local_pca_per_window <- function(theta_mat, n_sites_mat, pad) {
  n_samp <- nrow(theta_mat); n_win <- ncol(theta_mat)
  pc1 <- pc2 <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                       dimnames = list(rownames(theta_mat), NULL))
  l1 <- l2 <- lr <- rep(NA_real_, n_win)

  n_sites_chrom_median <- median(n_sites_mat, na.rm = TRUE)
  if (!is.finite(n_sites_chrom_median) || n_sites_chrom_median <= 0) n_sites_chrom_median <- 1

  for (wi in seq_len(n_win)) {
    lo <- max(1L, wi - pad); hi <- min(n_win, wi + pad)
    block <- theta_mat[, lo:hi, drop = FALSE]
    ok <- complete.cases(block)
    if (sum(ok) < max(20L, ncol(block) + 2L)) next
    block_ok <- block[ok, , drop = FALSE]
    nf <- n_sites_mat[ok, wi]
    nf[!is.finite(nf) | nf <= 0] <- 1L
    w <- sqrt(nf / n_sites_chrom_median)
    bw <- block_ok * w
    centred <- sweep(bw, 2, colMeans(bw), FUN = "-")
    sv <- tryCatch(svd(centred, nu = 2, nv = 0), error = function(e) NULL)
    if (is.null(sv)) next
    pc1[ok, wi] <- as.numeric(sv$u[, 1])
    if (ncol(sv$u) >= 2) pc2[ok, wi] <- as.numeric(sv$u[, 2])
    d2 <- sv$d^2
    if (length(d2) >= 1) l1[wi] <- d2[1]
    if (length(d2) >= 2 && d2[2] > 0) {
      l2[wi] <- d2[2]
      lr[wi] <- d2[1] / d2[2]
    }
  }
  list(pc1 = pc1, pc2 = pc2, lambda_1 = l1, lambda_2 = l2, lambda_ratio = lr)
}

# Banded OR full upper-triangle sim_mat from |cor(pc1[,i], pc1[,j])|.
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

# Reconstruct full sim matrix (transient; for cmdscale and dmat-derivable steps).
sim_full <- function(simbox, n_win) {
  if (simbox$format == "upper_triangle_float32") {
    M <- simbox$sim_mat
    M[lower.tri(M)] <- t(M)[lower.tri(M)]
    M
  } else {
    band_med <- median(simbox$sim_band, na.rm = TRUE)
    if (!is.finite(band_med)) band_med <- 0.0
    M <- matrix(band_med, n_win, n_win)
    bh <- simbox$band_half
    for (i in seq_len(n_win)) {
      jl <- max(1L, i - bh); jh <- min(n_win, i + bh)
      for (j in jl:jh) {
        v <- simbox$sim_band[i, j - i + bh + 1L]
        if (is.finite(v)) M[i, j] <- v
      }
    }
    diag(M) <- 1
    M
  }
}

# Anchor-flip sign alignment.
anchor_flip <- function(pc1_mat, pc2_mat, max_abs_z) {
  n_win <- ncol(pc1_mat)
  anchor <- which.max(max_abs_z)
  if (length(anchor) == 0L || !is.finite(max_abs_z[anchor])) {
    med_z <- median(max_abs_z, na.rm = TRUE)
    if (!is.finite(med_z)) med_z <- 0
    anchor <- which.min(abs(max_abs_z - med_z))[1]
  }
  anchor <- as.integer(anchor)
  p1 <- pc1_mat; p2 <- pc2_mat
  unflipped <- integer(0)
  a1 <- pc1_mat[, anchor]; a2 <- pc2_mat[, anchor]
  for (w in seq_len(n_win)) {
    if (w == anchor) next
    ok <- is.finite(a1) & is.finite(pc1_mat[, w])
    if (sum(ok) < 10L) { unflipped <- c(unflipped, w); next }
    r <- cor(a1[ok], pc1_mat[ok, w])
    if (is.finite(r) && r < 0) p1[, w] <- -pc1_mat[, w]
    ok2 <- is.finite(a2) & is.finite(pc2_mat[, w])
    if (sum(ok2) >= 10L) {
      r2 <- cor(a2[ok2], pc2_mat[ok2, w])
      if (is.finite(r2) && r2 < 0) p2[, w] <- -pc2_mat[, w]
    }
  }
  list(pc1_aligned = p1, pc2_aligned = p2, anchor_idx = anchor, unflipped = unflipped)
}

# NN-smoothed sim_mat via MDS-space k-NN averaging (mirrors z-blocks recipe).
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
  if (is.null(loaded)) next
  theta_mat <- loaded$theta_mat; n_sites_mat <- loaded$n_sites_mat
  win_grid  <- loaded$win_grid
  n_win     <- nrow(win_grid)

  # --- Per-window cohort metrics ---
  window_median <- apply(theta_mat, 2, median, na.rm = TRUE)
  dev_mat <- sweep(theta_mat, 2, window_median, FUN = "-")
  all_devs <- as.vector(dev_mat); all_devs <- all_devs[is.finite(all_devs)]
  chrom_dev_mad <- if (length(all_devs) > 100) mad(all_devs, constant = 1.4826) else NA_real_
  if (!is.finite(chrom_dev_mad) || chrom_dev_mad <= 0) chrom_dev_mad <- 1e-6

  max_abs_z <- apply(dev_mat, 2, function(v) {
    v <- v[is.finite(v)]; if (length(v) < 10) return(NA_real_)
    max(abs(v / chrom_dev_mad))
  })
  top10_abs_z <- apply(dev_mat, 2, function(v) {
    v <- abs(v[is.finite(v)] / chrom_dev_mad)
    if (length(v) < 10) return(NA_real_)
    k <- max(1L, ceiling(0.1 * length(v))); mean(sort(v, decreasing = TRUE)[seq_len(k)])
  })

  # --- Per-window local PCA (sample loadings + lambdas) ---
  pca <- local_pca_per_window(theta_mat, n_sites_mat, PAD)

  # --- Sim_mat (banded or full UT) ---
  simbox <- build_sim_mat(pca$pc1, n_win, SIM_BAND_HALF, SIM_N_FULL_THRESHOLD)
  message("[TR_B v5] ", chrom, ": sim_mat format = ", simbox$format)

  # --- Anchor-flip sign alignment ---
  flipped <- anchor_flip(pca$pc1, pca$pc2, max_abs_z)

  # --- Full sim matrix transient (for MDS + NN smoothing baseline) ---
  sim_M <- sim_full(simbox, n_win)
  dmat  <- 1 - sim_M
  diag(dmat) <- 0

  # --- 2D MDS ---
  mds_fit <- tryCatch(cmdscale(as.dist(dmat), k = 2L), error = function(e) NULL)
  if (is.null(mds_fit) || nrow(mds_fit) != n_win) {
    mds_mat <- matrix(NA_real_, n_win, 2)
  } else {
    mds_mat <- mds_fit
  }

  # --- bg_continuity quantiles (off-diagonal sim distribution near baseline) ---
  bg_q <- if (n_win >= 10) {
    adj_sims <- vapply(seq_len(n_win - 1), function(i) sim_M[i, i + 1L], numeric(1))
    quantile(adj_sims, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
  } else {
    rep(NA_real_, 6)
  }

  # --- NN-smoothed sim_mats ---
  per_chrom_outdir <- file.path(per_chrom_dir, chrom)
  dir.create(per_chrom_outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sim_M, file.path(per_chrom_outdir, "sim_mat_nn0.rds"))

  for (k in NN_SIM_SCALES) {
    k_use <- min(k, n_win - 1L)
    if (k_use < 2L) next
    t_nn <- proc.time()
    nn_sim <- make_nn_sim(mds_mat, dmat, k_use)
    saveRDS(nn_sim, file.path(per_chrom_outdir, sprintf("sim_mat_nn%d.rds", k)))
    message(sprintf("[TR_B v5] %s: sim_mat_nn%d (k_use=%d, %.1fs)",
                    chrom, k, k_use, (proc.time() - t_nn)[3]))
  }
  rm(sim_M); gc(verbose = FALSE)

  # --- Per-window data.table (mirrors z-blocks $dt fields) ---
  dt <- data.table(
    chrom            = chrom,
    window_idx       = win_grid$window_idx,
    start_bp         = win_grid$start_bp,
    end_bp           = win_grid$end_bp,
    mid_bp           = as.integer((win_grid$start_bp + win_grid$end_bp) / 2L),
    theta_pi_median  = round(window_median, 6),
    max_abs_z        = round(max_abs_z, 4),
    top10_abs_z      = round(top10_abs_z, 4),
    lambda_1         = round(pca$lambda_1, 6),
    lambda_2         = round(pca$lambda_2, 6),
    lambda_ratio     = round(pca$lambda_ratio, 4),
    MDS1             = round(mds_mat[, 1], 6),
    MDS2             = round(mds_mat[, 2], 6),
    anchor_window_idx = flipped$anchor_idx - 1L          # 0-indexed everywhere from v5
  )

  # ─── MORPHOLOGY FEATURES ─────────────────────────────────────────────────
  # Ported from local_PCA_MDS_z/03_precompute_localpca_zblocks.R (signal-
  # agnostic 1D shape descriptors of the |Z| score curve). Identical math,
  # different signal source — for theta_pi we feed `max_abs_z` here just as
  # z-blocks feeds its MDS-axis |Z|.
  z_vec <- max_abs_z
  z_vec[!is.finite(z_vec)] <- 0
  if (n_win >= 5 && any(z_vec > 0)) {
    # 1a. local jaggedness
    JAG <- 5L
    jagged <- vapply(seq_len(n_win), function(i) {
      lo <- max(2L, i - JAG); hi <- min(n_win, i + JAG)
      lz <- z_vec[lo:hi]
      if (length(lz) >= 3) mean(abs(diff(lz)), na.rm = TRUE) else NA_real_
    }, numeric(1))
    # 1b. elevated runs
    z_thresh <- max(quantile(z_vec, 0.75, na.rm = TRUE), 1.5)
    is_elev <- z_vec >= z_thresh
    rid <- integer(n_win); cur <- 0L
    for (i in seq_len(n_win)) {
      if (is_elev[i]) {
        if (i == 1L || !is_elev[i - 1L]) cur <- cur + 1L
        rid[i] <- cur
      }
    }
    run_len <- rep(0L, n_win); run_mean <- rep(NA_real_, n_win); run_sd <- rep(NA_real_, n_win)
    if (cur > 0) for (r in seq_len(cur)) {
      idx <- which(rid == r); rl <- length(idx)
      run_len[idx]  <- rl
      run_mean[idx] <- mean(z_vec[idx], na.rm = TRUE)
      run_sd[idx]   <- if (rl >= 2) sd(z_vec[idx], na.rm = TRUE) else 0
    }
    # 1c. peak prominence
    PROM <- 20L
    prom <- vapply(seq_len(n_win), function(i) {
      lo <- max(1L, i - PROM); hi <- min(n_win, i + PROM)
      bg <- quantile(z_vec[lo:hi], 0.25, na.rm = TRUE)
      z_vec[i] - bg
    }, numeric(1))
    # 1d. plateau flatness
    pf <- mapply(function(rl, rs, jg) {
      if (rl >= 3) rl * (1/(1 + (rs %||% 0))) * (1/(1 + (jg %||% 0))) else 0
    }, run_len, run_sd, jagged)
    pf_max <- max(pf, na.rm = TRUE); if (is.finite(pf_max) && pf_max > 0) pf <- pf / pf_max

    # 2. neighborhood support
    nb <- list()
    for (rad in c(5L, 10L, 20L)) {
      nb[[paste0("nbhood_support_", rad)]] <- vapply(seq_len(n_win), function(i) {
        lo <- max(1L, i - rad); hi <- min(n_win, i + rad)
        mean(is_elev[lo:hi], na.rm = TRUE)
      }, numeric(1))
    }
    # 3. sim_mat block morphology (uses banded or full M; we have sim_M already)
    SIM_HALF <- 10L
    sim_for_block <- if (exists("sim_M", inherits = FALSE)) sim_M else
                     sim_full(simbox, n_win)
    chr_sim_med <- median(sim_for_block[upper.tri(sim_for_block)], na.rm = TRUE)
    blk_compact <- blk_coher <- blk_frag <- sq_supp <- rep(NA_real_, n_win)
    for (i in seq_len(n_win)) {
      lo <- max(1L, i - SIM_HALF); hi <- min(n_win, i + SIM_HALF)
      if (hi - lo < 4L) next
      B <- sim_for_block[lo:hi, lo:hi, drop = FALSE]
      bv <- B[upper.tri(B)]
      if (length(bv) >= 3) blk_compact[i] <- mean(bv, na.rm = TRUE)
      dd <- abs(row(B) - col(B))
      nd <- B[dd <= 2 & dd > 0]; fd <- B[dd > max(2, (hi - lo) %/% 3)]
      if (length(nd) >= 2 && length(fd) >= 2 && mean(nd, na.rm = TRUE) > 0.01) {
        blk_coher[i] <- mean(fd, na.rm = TRUE) / mean(nd, na.rm = TRUE)
      }
      if (length(bv) >= 3) {
        bm <- mean(bv, na.rm = TRUE); bsd <- sd(bv, na.rm = TRUE)
        blk_frag[i] <- if (bm > 0.01) bsd / bm else 0
      }
      if (is.finite(chr_sim_med) && length(bv) >= 3) {
        sq_supp[i] <- mean(bv > chr_sim_med, na.rm = TRUE)
      }
    }
    # 4. composite scores (clamped to [0,1])
    clamp01 <- function(x) { x[!is.finite(x)] <- 0; pmin(1, pmax(0, x)) }
    na0     <- function(x) { x[!is.finite(x)] <- 0; x }

    flat_raw <- 0.22 * clamp01(log2(pmax(run_len, 1)) / 5) +
                0.16 * clamp01(na0(run_mean) / 4) +
                0.12 * clamp01(1 - na0(run_sd) / 2) +
                0.10 * clamp01(1 - na0(jagged) / 2) +
                0.16 * clamp01(na0(blk_compact)) +
                0.12 * clamp01(na0(blk_coher)) +
                0.12 * clamp01(1 - na0(blk_frag))
    flat_raw[run_len < 3] <- 0

    spiky_raw <- 0.28 * clamp01(na0(prom) / 3) +
                 0.17 * clamp01(na0(nb$nbhood_support_5)) +
                 0.17 * clamp01(na0(blk_compact)) +
                 0.12 * clamp01(1 - na0(blk_frag) * 0.5) +
                 0.12 * clamp01(ifelse(run_len >= 3 & run_len <= 12, 1,
                                       ifelse(run_len > 12, 0.5, 0.3))) +
                 0.14 * clamp01(1 - na0(pf) * 0.5)
    spiky_raw[z_vec < z_thresh * 0.5] <- 0

    sup5  <- na0(nb$nbhood_support_5);  sup20 <- na0(nb$nbhood_support_20)
    frag_raw <- 0.25 * clamp01(na0(jagged) / 2) +
                0.20 * clamp01(na0(blk_frag)) +
                0.15 * clamp01(1 - na0(blk_coher)) +
                0.20 * clamp01(pmax(0, sup5 - sup20)) +
                0.20 * clamp01(ifelse(run_len >= 1 & run_len <= 2, 1,
                                      ifelse(run_len == 3, 0.5, 0)))
    frag_raw[z_vec < 1.0] <- 0

    dt[, `:=`(
      local_jaggedness        = round(jagged, 4),
      local_run_len           = run_len,
      local_run_mean_z        = round(run_mean, 4),
      local_run_sd_z          = round(run_sd, 4),
      local_peak_prominence   = round(prom, 4),
      plateau_flatness        = round(pf, 4),
      nbhood_support_5        = round(nb$nbhood_support_5, 4),
      nbhood_support_10       = round(nb$nbhood_support_10, 4),
      nbhood_support_20       = round(nb$nbhood_support_20, 4),
      local_block_compactness = round(blk_compact, 4),
      local_block_coherence   = round(blk_coher, 4),
      local_block_fragmentation = round(blk_frag, 4),
      local_square_support    = round(sq_supp, 4),
      flat_inv_score          = round(flat_raw, 4),
      spiky_inv_score         = round(spiky_raw, 4),
      fragmentation_score     = round(frag_raw, 4)
    )]
    rm(sim_for_block); gc(verbose = FALSE)
  }

  # ─── θπ-flavored inv_likeness composite ──────────────────────────────────
  # 50% normalized max_abs_z + 30% sim_mat block compactness + 20% λ-ratio
  # elevation. Replaces z-blocks' dosage-PC1-trimodality formula. Higher
  # values mean a window looks like an inversion-anchored region under
  # cohort-θπ + local-structure criteria.
  z_norm <- pmin(1, pmax(0, max_abs_z / 6))         # |Z|=6 saturates
  bc <- if ("local_block_compactness" %in% names(dt)) dt$local_block_compactness else NA
  bc_n <- pmin(1, pmax(0, bc))
  lr <- pca$lambda_ratio
  lr_n <- pmin(1, pmax(0, (lr - 1) / 4))            # λ_ratio in [1,5] → [0,1]
  inv_like <- 0.50 * ifelse(is.finite(z_norm), z_norm, 0) +
              0.30 * ifelse(is.finite(bc_n),  bc_n,  0) +
              0.20 * ifelse(is.finite(lr_n),  lr_n,  0)
  dt[, inv_likeness := round(inv_like, 4)]

  # ─── Beta-adaptive p-values ──────────────────────────────────────────────
  # Fits Beta(α, β) to the chrom-wide inv_likeness distribution and scores
  # each window's right-tail p-value. Adaptive_seed = TRUE for p < 0.01.
  beta_pvals <- function(scores) {
    valid <- scores[is.finite(scores) & scores > 0.001 & scores < 0.999]
    if (length(valid) < 50) return(list(p = rep(NA_real_, length(scores)),
                                         a = NA_real_, b = NA_real_))
    fit <- tryCatch({
      if (requireNamespace("MASS", quietly = TRUE)) {
        MASS::fitdistr(valid, "beta", start = list(shape1 = 1, shape2 = 5))
      } else {
        m <- mean(valid); v <- var(valid)
        if (v >= m * (1 - m)) v <- m * (1 - m) * 0.99
        a <- m * (m * (1 - m) / v - 1); b <- (1 - m) * (m * (1 - m) / v - 1)
        list(estimate = c(shape1 = max(0.1, a), shape2 = max(0.1, b)))
      }
    }, error = function(e) list(estimate = c(shape1 = 1, shape2 = 5)))
    a <- fit$estimate["shape1"]; b <- fit$estimate["shape2"]
    list(p = 1 - pbeta(scores, a, b), a = a, b = b)
  }
  bp <- beta_pvals(dt$inv_likeness)
  dt[, `:=`(
    beta_pval     = round(bp$p, 6),
    adaptive_seed = is.finite(bp$p) & bp$p < 0.01,
    beta_alpha    = round(bp$a, 3),
    beta_beta     = round(bp$b, 3)
  )]
  message(sprintf("[TR_B v5] %s: inv_likeness median=%.3f q95=%.3f, %d adaptive seeds (p<0.01)",
                  chrom, median(dt$inv_likeness, na.rm = TRUE),
                  quantile(dt$inv_likeness, 0.95, na.rm = TRUE),
                  sum(dt$adaptive_seed, na.rm = TRUE)))
  # Per-sample PC1/PC2 columns (PC_1_<sample>, PC_2_<sample>) — matches z-blocks
  # so the atlas JSON exporter and downstream scripts can grep `^PC_[12]_`.
  for (s_idx in seq_len(n_samp)) {
    sid <- sample_order[s_idx]
    set(dt, j = paste0("PC_1_", sid), value = round(flipped$pc1_aligned[s_idx, ], 6))
    set(dt, j = paste0("PC_2_", sid), value = round(flipped$pc2_aligned[s_idx, ], 6))
  }

  # --- Save precomp.rds ---
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
    sample_order             = sample_order,
    unflipped_windows        = flipped$unflipped - 1L
  )
  saveRDS(precomp, file.path(per_chrom_outdir, "precomp.rds"))
  message("[TR_B v5] ", chrom, ": precomp.rds saved (", n_win, " windows)")

  # --- Genome-wide rollup row + per-chrom QC summary ---
  all_window_dt[[chrom]] <- dt[, .(chrom, window_idx, start_bp, end_bp, mid_bp,
                                    theta_pi_median, max_abs_z, top10_abs_z,
                                    lambda_ratio, MDS1, MDS2)]
  elapsed <- round((proc.time() - t_chr)[3], 1)
  all_summary[[chrom]] <- data.table(
    chrom         = chrom,
    n_windows     = as.integer(n_win),
    n_samples     = as.integer(n_samp),
    sim_format    = simbox$format,
    n_unflipped   = length(flipped$unflipped),
    median_max_z  = round(median(max_abs_z, na.rm = TRUE), 3),
    q95_max_z     = round(quantile(max_abs_z, 0.95, na.rm = TRUE), 3),
    bg_q50        = round(bg_q["50%"], 4),
    bg_q90        = round(bg_q["90%"], 4),
    bg_q95        = round(bg_q["95%"], 4),
    elapsed_sec   = elapsed
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