#!/usr/bin/env Rscript
# =============================================================================
# STEP_GH_C_mds_compute.R  (v1, 2026-05-12)
# =============================================================================
# MDS + sim_mat + per-window features + (dense) sketch & sparse k-NN graph,
# for the GHSL path (path 3). Direct mirror of θπ STEP_TR_C v6 — same BLAS
# sim_mat, same fast_cmdscale via RSpectra, same sketch + sparse_edges in
# --mode local, same cohort sketch basis (reused from θπ when cohort_id
# matches, giving free cross-path comparability in the atlas).
#
# Inputs:
#   $OUT_GH_LOCAL_PCA_DIR/<chr>.ghsl_window_pca.rds  (from STEP_GH_B)
#
# Outputs:
#   $OUTROOT/<--out_subdir>/<chr>.ghsl_precomp.rds         (per-chrom dt + sim)
#   $OUTROOT/<--out_subdir>/sim_mats/<chr>.sim_mat_nn{N}.rds (full mode only)
#   $OUTROOT/<--out_subdir>/window_dt.tsv.gz               (genome rollup)
#   $OUTROOT/<--out_subdir>/precomp_summary.tsv            (per-chrom QC)
#
# Local-mode (dense) additional outputs (ADR-0001 Layer C):
#   $SKETCH_DIR/sample_sketch_basis.rds       (shared with θπ — reuse if exists)
#   $SKETCH_DIR/<chr>.ghsl_sketch.rds         (n_win × d_sketch)
#   $SKETCH_DIR/<chr>.ghsl_sparse_edges.tsv.gz
#
# Note the "ghsl_" filename prefix on sketch/sparse_edges so they coexist
# with the θπ outputs (<chr>.sketch.rds, <chr>.sparse_edges.tsv.gz) under
# the same $SKETCH_DIR without collision.
#
# Knobs (identical to TR_C v6 except where flagged GHSL):
#   --chrom <CHR>                 single chrom
#   --mode <full|local>           default "full"
#   --kmds <int>                  default 5
#   --sim-band-half <int>         default 200
#   --sim-n-full-threshold <int>  default 6000
#   --out_subdir <name>           default "ghsl_precomp"           [GHSL]
#   --d_sketch <int>              default 32
#   --k_nn <int>                  default 50
#   --cohort_id <str>             default $COHORT_ID or "default_cohort"
#   --rng_seed <int>              default 20260512
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

K_MDS                <- 5L  # Top-5 MDS axes for sim_mat NN-smoothing space.
                             # sim_mat itself is correlation-based on per-window
                             # PC1 vectors (see sim_mat_blas), so K_MDS only
                             # affects the smoothing dimensionality. Yeah if u
                             # so homeless u can K=2 to try, but here we don't
                             # have time to test. K=5 is the safe default.
SEED_MDS_AXES        <- 5L
SIM_BAND_HALF        <- 200L
SIM_N_FULL_THRESHOLD <- 6000L
NN_SIM_SCALES        <- c(20, 40, 80, 120, 160, 200, 240, 320)
MODE                 <- "full"
OUT_SUBDIR           <- "ghsl_precomp"
D_SKETCH             <- 32L
K_NN                 <- 50L
COHORT_ID            <- Sys.getenv("COHORT_ID", unset = "default_cohort")
RNG_SEED             <- 20260512L

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
  else if (a == "--d_sketch"             && i < length(args)) { D_SKETCH <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--k_nn"                 && i < length(args)) { K_NN <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--cohort_id"            && i < length(args)) { COHORT_ID <- args[i + 1]; i <- i + 2L }
  else if (a == "--rng_seed"             && i < length(args)) { RNG_SEED <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}
stopifnot(MODE %in% c("full", "local"))

OUTROOT             <- Sys.getenv("OUTROOT", unset = NA)
OUT_GH_LOCAL_PCA_DIR <- Sys.getenv("OUT_GH_LOCAL_PCA_DIR",
                                   unset = file.path(OUTROOT, "01_ghsl_local_pca"))
SKETCH_DIR <- Sys.getenv("SKETCH_DIR",
                         unset = file.path(OUTROOT, "03_dense_sketch"))
stopifnot(!is.na(OUTROOT))

precomp_dir <- file.path(OUTROOT, OUT_SUBDIR)
sim_dir     <- file.path(precomp_dir, "sim_mats")
dir.create(precomp_dir, recursive = TRUE, showWarnings = FALSE)
if (MODE == "full")  dir.create(sim_dir,    recursive = TRUE, showWarnings = FALSE)
if (MODE == "local") dir.create(SKETCH_DIR, recursive = TRUE, showWarnings = FALSE)

if (is.null(CHROM)) {
  cl <- Sys.getenv("CHROM_LIST", unset = "")
  CHROM_LIST <- if (nchar(cl) > 0L) strsplit(cl, "[ ,]+")[[1]] else
    sub("\\.ghsl_window_pca\\.rds$", "",
        basename(list.files(OUT_GH_LOCAL_PCA_DIR,
                            pattern = "\\.ghsl_window_pca\\.rds$")))
} else {
  CHROM_LIST <- CHROM
}

HAVE_RSPECTRA <- requireNamespace("RSpectra", quietly = TRUE)
HAVE_HNSW     <- requireNamespace("RcppHNSW", quietly = TRUE)
if (MODE == "full" && !HAVE_RSPECTRA)
  warning("RSpectra not installed — falling back to cmdscale (slow).")
if (MODE == "local" && !HAVE_HNSW)
  warning("RcppHNSW not installed — sketch + sparse_edges will be skipped.")

message(sprintf("[GH_C] mode=%s K_MDS=%d SIM_BAND_HALF=%d out=%s",
                MODE, K_MDS, SIM_BAND_HALF, precomp_dir))
if (MODE == "local")
  message(sprintf("[GH_C] sketch: d=%d k_nn=%d cohort=%s seed=%d dir=%s",
                  D_SKETCH, K_NN, COHORT_ID, RNG_SEED, SKETCH_DIR))
message("[GH_C] chroms: ", length(CHROM_LIST))

# =============================================================================
# Helpers (identical math to TR_C v6 — see that file for derivations)
# =============================================================================

sim_mat_blas <- function(pc1_mat) {
  X <- pc1_mat
  col_means <- colMeans(X, na.rm = TRUE); col_means[!is.finite(col_means)] <- 0
  na_idx <- which(!is.finite(X), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) X[na_idx] <- col_means[na_idx[, 2L]]
  X <- sweep(X, 2L, col_means, FUN = "-")
  col_norm <- sqrt(colSums(X * X)); col_norm[col_norm < 1e-12] <- 1
  X <- sweep(X, 2L, col_norm, FUN = "/")
  M <- abs(crossprod(X)); diag(M) <- 1
  M
}

pack_banded <- function(M, band_half) {
  n <- nrow(M); nb <- 2L * band_half + 1L
  B <- matrix(NA_real_, n, nb)
  for (i in seq_len(n)) {
    jl <- max(1L, i - band_half); jh <- min(n, i + band_half)
    B[i, (jl - i + band_half + 1L):(jh - i + band_half + 1L)] <- M[i, jl:jh]
  }
  B
}

build_sim_mat <- function(pc1_mat, n_win, band_half, n_full_threshold,
                          force_banded = FALSE) {
  M <- sim_mat_blas(pc1_mat)
  if (force_banded || n_win > n_full_threshold) {
    B <- pack_banded(M, band_half); rm(M); gc(verbose = FALSE)
    list(sim_mat = NULL, sim_band = B,
         format = paste0("banded_float32_pm", band_half),
         band_half = as.integer(band_half))
  } else {
    list(sim_mat = M, sim_band = NULL,
         format = "upper_triangle_float32", band_half = NA_integer_)
  }
}

sim_full <- function(simbox, n_win) {
  if (simbox$format == "upper_triangle_float32") return(simbox$sim_mat)
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

fast_cmdscale <- function(dmat, k) {
  n <- nrow(dmat)
  if (!HAVE_RSPECTRA) {
    return(tryCatch(cmdscale(as.dist(dmat), k = k), error = function(e) NULL))
  }
  D2 <- dmat * dmat
  rm <- rowMeans(D2); gm <- mean(D2)
  B <- -0.5 * (D2 - outer(rm, rep(1, n)) - outer(rep(1, n), rm) + gm)
  B <- 0.5 * (B + t(B))
  eig <- tryCatch(RSpectra::eigs_sym(B, k = k, which = "LA"), error = function(e) NULL)
  if (is.null(eig)) return(NULL)
  vals <- pmax(eig$values, 0)
  sweep(eig$vectors, 2L, sqrt(vals), FUN = "*")
}

make_nn_sim <- function(mds_mat, dmat, k_use) {
  n <- nrow(mds_mat); smoothed <- matrix(0, n, ncol(mds_mat))
  for (wi in seq_len(n)) {
    d <- dmat[wi, ]; d[wi] <- Inf
    cutoff <- sort(d, partial = k_use)[k_use]
    nn_idx <- which(d <= cutoff); if (length(nn_idx) > k_use) nn_idx <- nn_idx[seq_len(k_use)]
    smoothed[wi, ] <- colMeans(mds_mat[c(wi, nn_idx), , drop = FALSE], na.rm = TRUE)
  }
  d_nn <- as.matrix(dist(smoothed)); fv <- d_nn[is.finite(d_nn)]
  dmax <- if (length(fv) > 0L) quantile(fv, 0.95, na.rm = TRUE) else 1
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  sim <- 1 - pmin(d_nn / dmax, 1); sim[!is.finite(sim)] <- 0
  diag(sim) <- 1; sim
}

# Layer C: sketch basis is COHORT-scoped, not path-scoped. If $SKETCH_DIR
# already has sample_sketch_basis.rds (likely from a prior θπ TR_C run),
# we reuse it — the projection is feature-agnostic and shared basis means
# θπ and GHSL sketches live in the same 32-d space (atlas can ask "which
# dense GHSL windows look like this dense θπ window").
load_or_create_basis <- function(sketch_dir, n_samp, sample_order, d_sketch,
                                 seed, cohort_id) {
  basis_path <- file.path(sketch_dir, "sample_sketch_basis.rds")
  if (file.exists(basis_path)) {
    b <- readRDS(basis_path)
    if (b$d_sketch == d_sketch && b$cohort_id == cohort_id &&
        b$seed == seed && nrow(b$basis) == n_samp &&
        identical(b$sample_order, sample_order)) {
      message("[GH_C] sketch basis: reusing ", basis_path,
              " (shared with θπ path)")
      return(b)
    } else {
      stop("[GH_C] sketch basis at ", basis_path,
           " has mismatched config (cohort_id/seed/d_sketch/sample_order). ",
           "Delete it or change --cohort_id / --rng_seed.")
    }
  }
  set.seed(seed)
  basis_R <- matrix(rnorm(n_samp * d_sketch, 0, 1 / sqrt(d_sketch)),
                    nrow = n_samp, ncol = d_sketch)
  b <- list(basis = basis_R, sample_order = sample_order,
            d_sketch = as.integer(d_sketch), seed = as.integer(seed),
            cohort_id = cohort_id, created_at = Sys.time(),
            toolkit_version = "inversion-popgen-toolkit (GH_C v1)")
  saveRDS(b, basis_path)
  message("[GH_C] sketch basis: created ", basis_path)
  b
}

build_window_sketch <- function(feat_mat, basis_R) {
  X <- feat_mat
  col_means <- colMeans(X, na.rm = TRUE); col_means[!is.finite(col_means)] <- 0
  na_idx <- which(!is.finite(X), arr.ind = TRUE)
  na_count <- if (nrow(na_idx) > 0L) tabulate(na_idx[, 2L], nbins = ncol(X))
              else integer(ncol(X))
  if (nrow(na_idx) > 0L) X[na_idx] <- col_means[na_idx[, 2L]]
  S <- crossprod(X, basis_R)
  rn <- sqrt(rowSums(S * S)); rn[rn < 1e-12] <- 1
  S <- sweep(S, 1L, rn, FUN = "/")
  list(sketch = S, na_count = na_count)
}

# See TR_C v6 docstring for the sign / anchor-flip note. Same caveat here.
build_sparse_edges <- function(sketch_mat, k_nn, chrom) {
  n_win <- nrow(sketch_mat); k_query <- min(k_nn + 1L, n_win)
  ann <- RcppHNSW::hnsw_knn(sketch_mat, k = k_query, distance = "cosine",
                            M = 16L, ef_construction = 200L, ef = 100L)
  idx <- ann$idx; dist <- ann$dist
  edges <- vector("list", n_win)
  for (wa in seq_len(n_win)) {
    nbrs <- idx[wa, ]; dists <- dist[wa, ]
    keep <- nbrs != wa
    nbrs <- nbrs[keep]; dists <- dists[keep]
    if (length(nbrs) > k_nn) { nbrs <- nbrs[seq_len(k_nn)]; dists <- dists[seq_len(k_nn)] }
    if (length(nbrs) == 0L) next
    sim_cos <- pmax(0, pmin(1, abs(1 - dists)))   # |cor|-semantics
    edges[[wa]] <- data.table(
      chrom = chrom, window_a = as.integer(wa - 1L),
      window_b = as.integer(nbrs - 1L), rank = seq_along(nbrs),
      sim_cosine = round(sim_cos, 6), hop_distance = 1L
    )
  }
  rbindlist(edges, use.names = TRUE)
}

compute_sketch_novelty <- function(edges_dt, n_win, excl = 5L) {
  nov <- rep(NA_real_, n_win)
  if (nrow(edges_dt) == 0L) return(nov)
  ed <- edges_dt[abs(window_a - window_b) > excl]
  if (nrow(ed) == 0L) return(nov)
  agg <- ed[, .(max_sim = max(sim_cosine, na.rm = TRUE)), by = window_a]
  nov[agg$window_a + 1L] <- 1 - agg$max_sim
  nov
}

# =============================================================================
# Per-chrom MDS + features
# =============================================================================
all_window_dt <- list(); all_summary <- list(); basis_cached <- NULL

for (chrom in CHROM_LIST) {
  t_chr <- proc.time()
  pca_file <- file.path(OUT_GH_LOCAL_PCA_DIR,
                        sprintf("%s.ghsl_window_pca.rds", chrom))
  if (!file.exists(pca_file)) {
    message("[GH_C] ", chrom, ": missing ", pca_file, " — skip"); next
  }
  pca <- readRDS(pca_file)
  n_win <- pca$n_windows; n_samp <- pca$n_samples; NPC <- pca$npc
  win_grid <- pca$win_grid
  pcs_aligned    <- pca$pcs
  lambda         <- pca$lambda
  div_z_direct   <- pca$div_z_direct
  window_median  <- pca$window_median
  sample_order   <- pca$sample_order
  unflipped      <- pca$unflipped_windows
  anchor_idx     <- pca$anchor_idx
  scale_used     <- pca$scale

  lambda_ratio <- ifelse(is.finite(lambda[, 2L]) & lambda[, 2L] > 0,
                         lambda[, 1L] / lambda[, 2L], NA_real_)

  t_sim <- proc.time()
  simbox <- build_sim_mat(pcs_aligned[[1L]], n_win, SIM_BAND_HALF,
                          SIM_N_FULL_THRESHOLD, force_banded = (MODE == "local"))
  message(sprintf("[GH_C] %s: sim_mat (%s) %.1fs",
                  chrom, simbox$format, (proc.time() - t_sim)[3]))

  if (MODE == "full") {
    sim_M <- sim_full(simbox, n_win)
    dmat <- 1 - sim_M; diag(dmat) <- 0
    t_mds <- proc.time()
    mds_fit <- fast_cmdscale(dmat, K_MDS)
    mds_mat <- if (is.null(mds_fit) || nrow(mds_fit) != n_win) {
      matrix(NA_real_, n_win, K_MDS)
    } else if (ncol(mds_fit) < K_MDS) {
      cbind(mds_fit, matrix(NA_real_, n_win, K_MDS - ncol(mds_fit)))
    } else mds_fit
    message(sprintf("[GH_C] %s: fast_cmdscale k=%d %.1fs",
                    chrom, K_MDS, (proc.time() - t_mds)[3]))

    mds_z <- matrix(NA_real_, n_win, K_MDS)
    for (k in seq_len(K_MDS)) {
      v <- mds_mat[, k]
      if (sum(is.finite(v)) >= 10L) {
        med <- median(v, na.rm = TRUE); md <- mad(v, na.rm = TRUE)
        if (is.finite(md) && md > 1e-12) mds_z[, k] <- (v - med) / md
        else {
          sdv <- sd(v, na.rm = TRUE)
          if (is.finite(sdv) && sdv > 0) mds_z[, k] <- (v - mean(v, na.rm = TRUE)) / sdv
        }
      }
    }
    z_for_max <- mds_z[, seq_len(min(SEED_MDS_AXES, K_MDS)), drop = FALSE]
    max_abs_z <- apply(z_for_max, 1L, function(r) {
      r <- r[is.finite(r)]; if (length(r) == 0L) NA_real_ else max(abs(r))
    })
    max_z_axis <- apply(z_for_max, 1L, function(r) {
      a <- which.max(abs(r)); if (length(a) == 0L) NA_integer_ else as.integer(a)
    })
    bg_q <- if (n_win >= 10L) {
      adj <- vapply(seq_len(n_win - 1L), function(i) sim_M[i, i + 1L], numeric(1))
      quantile(adj, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
    } else rep(NA_real_, 6L)

    saveRDS(sim_M, file.path(sim_dir, sprintf("%s.sim_mat_nn0.rds", chrom)))
    for (kk in NN_SIM_SCALES) {
      k_use <- min(kk, n_win - 1L); if (k_use < 2L) next
      t_nn <- proc.time()
      nn_sim <- make_nn_sim(mds_mat, dmat, k_use)
      saveRDS(nn_sim, file.path(sim_dir, sprintf("%s.sim_mat_nn%d.rds", chrom, kk)))
      message(sprintf("[GH_C] %s: sim_mat_nn%d (%.1fs)",
                      chrom, kk, (proc.time() - t_nn)[3]))
    }
    rm(sim_M, dmat); gc(verbose = FALSE)
    sketch_path <- NA_character_; sparse_edges_path <- NA_character_
    sketch_novelty <- rep(NA_real_, n_win)
  } else {
    mds_mat   <- matrix(NA_real_, n_win, K_MDS)
    mds_z     <- matrix(NA_real_, n_win, K_MDS)
    max_abs_z <- rep(NA_real_, n_win)
    max_z_axis <- rep(NA_integer_, n_win)
    bg_q       <- rep(NA_real_, 6L)
    sketch_novelty <- rep(NA_real_, n_win)
    sketch_path <- NA_character_; sparse_edges_path <- NA_character_

    if (HAVE_HNSW) {
      if (is.null(basis_cached)) {
        basis_cached <- load_or_create_basis(
          SKETCH_DIR, n_samp, sample_order, D_SKETCH, RNG_SEED, COHORT_ID)
      }
      t_sk <- proc.time()
      sk <- build_window_sketch(pcs_aligned[[1L]], basis_cached$basis)
      message(sprintf("[GH_C] %s: sketch (%d × %d) %.1fs",
                      chrom, n_win, D_SKETCH, (proc.time() - t_sk)[3]))

      sketch_obj <- list(
        chrom = chrom, sketch = sk$sketch,
        sketch_track = paste0("PC1_GHSL_", scale_used),
        window_idx = as.integer(win_grid$window_idx),
        d_sketch = as.integer(D_SKETCH),
        basis_id = paste0(COHORT_ID, "_seed", RNG_SEED, "_d", D_SKETCH),
        n_win = as.integer(n_win), na_count = sk$na_count
      )
      sketch_path <- file.path(SKETCH_DIR, sprintf("%s.ghsl_sketch.rds", chrom))
      saveRDS(sketch_obj, sketch_path)

      t_gr <- proc.time()
      edges_dt <- build_sparse_edges(sk$sketch, K_NN, chrom)
      sparse_edges_path <- file.path(SKETCH_DIR,
                                     sprintf("%s.ghsl_sparse_edges.tsv.gz", chrom))
      fwrite(edges_dt, sparse_edges_path, sep = "\t", compress = "gzip")
      message(sprintf("[GH_C] %s: sparse_edges (k=%d, %d rows) %.1fs",
                      chrom, K_NN, nrow(edges_dt), (proc.time() - t_gr)[3]))

      sketch_novelty <- compute_sketch_novelty(edges_dt, n_win, excl = 5L)
    } else {
      message("[GH_C] ", chrom, ": skipping sketch (RcppHNSW missing)")
    }
  }

  dt <- data.table(
    chrom            = chrom,
    window_idx       = win_grid$window_idx,
    start_bp         = win_grid$start_bp,
    end_bp           = win_grid$end_bp,
    mid_bp           = as.integer((win_grid$start_bp + win_grid$end_bp) / 2L),
    div_median       = round(window_median, 6),
    div_z_direct     = round(div_z_direct, 4),
    max_abs_z        = round(max_abs_z, 4),
    max_z_axis       = max_z_axis,
    lambda_ratio     = round(lambda_ratio, 4),
    anchor_window_idx = anchor_idx - 1L,
    sketch_novelty   = round(sketch_novelty, 4)
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
    dt = dt, sim_mat = simbox$sim_mat, sim_band = simbox$sim_band,
    sim_mat_format = simbox$format, sim_band_half = simbox$band_half,
    mds_mat = mds_mat, bg_continuity_quantiles = bg_q,
    chrom = chrom, n_windows = as.integer(n_win),
    n_samples = as.integer(n_samp), npc = as.integer(NPC),
    k_mds = as.integer(K_MDS), sample_order = sample_order,
    unflipped_windows = unflipped - 1L, mode = MODE,
    scale = scale_used,
    sketch_path = sketch_path, sparse_edges_path = sparse_edges_path
  )
  saveRDS(precomp, file.path(precomp_dir, sprintf("%s.ghsl_precomp.rds", chrom)))
  message(sprintf("[GH_C] %s: ghsl_precomp.rds saved (%d wins, scale=%s, mode=%s)",
                  chrom, n_win, scale_used, MODE))

  all_window_dt[[chrom]] <- dt[, .(chrom, window_idx, start_bp, end_bp, mid_bp,
                                    div_median, div_z_direct,
                                    max_abs_z, max_z_axis, lambda_ratio,
                                    MDS1, MDS2, sketch_novelty)]
  elapsed <- round((proc.time() - t_chr)[3], 1)
  all_summary[[chrom]] <- data.table(
    chrom = chrom, mode = MODE, scale = scale_used,
    n_windows = as.integer(n_win), n_samples = as.integer(n_samp),
    npc = as.integer(NPC), k_mds = as.integer(K_MDS),
    sim_format = simbox$format, n_unflipped = length(unflipped),
    median_max_z = round(median(max_abs_z, na.rm = TRUE), 3),
    q95_max_z = round(quantile(max_abs_z, 0.95, na.rm = TRUE), 3),
    median_div_z = round(median(div_z_direct, na.rm = TRUE), 3),
    bg_q50 = round(bg_q[1], 4), bg_q90 = round(bg_q[5], 4),
    bg_q95 = round(bg_q[6], 4),
    median_novelty = round(median(sketch_novelty, na.rm = TRUE), 4),
    elapsed_sec = elapsed
  )
  message("[GH_C] ", chrom, ": done in ", elapsed, "s")
}

if (length(all_window_dt) > 0L) {
  fwrite(rbindlist(all_window_dt, fill = TRUE),
         file.path(precomp_dir, "window_dt.tsv.gz"), sep = "\t", compress = "gzip")
}
if (length(all_summary) > 0L) {
  fwrite(rbindlist(all_summary, fill = TRUE),
         file.path(precomp_dir, "precomp_summary.tsv"), sep = "\t")
}
message("[GH_C] DONE")
