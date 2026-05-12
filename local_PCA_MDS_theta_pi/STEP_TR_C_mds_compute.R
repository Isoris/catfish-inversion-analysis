#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_C_mds_compute.R  (v6)
# =============================================================================
# MDS + sim_mat + per-window features + (dense) sketch & sparse k-NN graph.
#
# v6 changes vs v5:
#   1. build_sim_mat() replaced with BLAS-based tcrossprod path. The
#      nested-loop cor() over O(N^2) pairs becomes one matrix product.
#      Numerical caveat: NA cells in PC1 are imputed to the column mean
#      before centering, so high-missingness windows get attenuated
#      similarity rather than the v5 per-pair "intersection of finite"
#      behavior. On well-covered data (LG28 coarse, near-full coverage
#      after local-PCA mask) cor(v5_sim, v6_sim) >= 0.9999 is the
#      acceptance threshold (ADR-0001 §6.1).
#   2. cmdscale() replaced with fast_cmdscale() (RSpectra::eigs_sym on
#      the double-centered B matrix). Top-K_MDS axes, signs may flip
#      vs v5 — compare loadings via procrustes or |abs|.
#   3. NN-smoothed sim_mat construction uses a partial-sort top-k
#      neighbor query instead of full order(). Smoothed-MDS distance
#      matrix is computed once and reused across NN scales.
#   4. --mode local (dense scale) now also writes a sample-dim random
#      projection sketch and an HNSW-based sparse k-NN edge list. See
#      ADR-0001 §3.2–3.3, §4. These are additive outputs — banded
#      sim_mat is still written, so existing consumers continue working.
#      Sketch outputs are skipped if RcppHNSW is not installed (warning,
#      not error).
#
# Inputs (configured via 00_theta_config.sh):
#   $OUT_LOCAL_PCA_DIR/<chr>.window_pca.rds   (from STEP_TR_B v6)
#
# Output (under $OUTROOT/<--out_subdir>, default subdir = "precomp"):
#   <chr>.precomp.rds          per-window dt + (full mode) mds_mat + sim_mat
#   sim_mats/<chr>.sim_mat_nn{N}.rds   (full mode only)
#   ../window_dt.tsv.gz        genome-wide rollup (when iterating CHROM_LIST)
#   ../precomp_summary.tsv     per-chrom QC stats
#
# Additional output (local mode only, ADR-0001 Layer C):
#   $SKETCH_DIR/sample_sketch_basis.rds          (one per cohort, created once)
#   $SKETCH_DIR/<chr>.sketch.rds                 (n_win × d_sketch sketch)
#   $SKETCH_DIR/<chr>.sparse_edges.tsv.gz        (top-k_nn HNSW neighbors)
#
# Knobs:
#   --chrom <CHR>                 single chrom (default: iterate CHROM_LIST)
#   --mode <full|local>           default "full"
#   --kmds <int>                  MDS axes (default 5; only used in full mode)
#   --sim-band-half <int>         banded sim_mat half-width in WINDOWS (default 200)
#   --sim-n-full-threshold <int>  switch full↔banded storage (default 6000)
#   --out_subdir <name>           output subfolder under $OUTROOT (default "precomp")
#   --d_sketch <int>              sketch dimension (default 32; only used in local mode)
#   --k_nn <int>                  HNSW neighbors per window (default 50; local mode)
#   --sketch_track <str>          "PC1" (default) or "theta_raw"  — see ADR §3.2
#   --cohort_id <str>             cohort tag stored in basis file
#                                 (default $COHORT_ID env var or "default_cohort")
#   --rng_seed <int>              seed for sketch basis (default 20260512)
#
# Memory cost: same as v5 in full mode (N²·8B during sim_full + MDS), but
# wall-clock is seconds instead of hours at N=6k. Local mode is linear in N.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

K_MDS                <- 5L
SEED_MDS_AXES        <- 5L
SIM_BAND_HALF        <- 200L
SIM_N_FULL_THRESHOLD <- 6000L
NN_SIM_SCALES        <- c(20, 40, 80, 120, 160, 200, 240, 320)
MODE                 <- "full"
OUT_SUBDIR           <- "precomp"

# Layer C defaults (only used when MODE == "local")
D_SKETCH      <- 32L
K_NN          <- 50L
SKETCH_TRACK  <- "PC1"
COHORT_ID     <- Sys.getenv("COHORT_ID", unset = "default_cohort")
RNG_SEED      <- 20260512L

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
  else if (a == "--sketch_track"         && i < length(args)) { SKETCH_TRACK <- args[i + 1]; i <- i + 2L }
  else if (a == "--cohort_id"            && i < length(args)) { COHORT_ID <- args[i + 1]; i <- i + 2L }
  else if (a == "--rng_seed"             && i < length(args)) { RNG_SEED <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}
stopifnot(MODE %in% c("full", "local"))
stopifnot(SKETCH_TRACK %in% c("PC1", "theta_raw"))

OUTROOT       <- Sys.getenv("OUTROOT", unset = NA)
OUT_LOCAL_PCA_DIR <- Sys.getenv("OUT_LOCAL_PCA_DIR",
                                unset = file.path(OUTROOT, "01_local_pca"))
SKETCH_DIR    <- Sys.getenv("SKETCH_DIR",
                            unset = file.path(OUTROOT, "03_dense_sketch"))
stopifnot(!is.na(OUTROOT))

precomp_dir <- file.path(OUTROOT, OUT_SUBDIR)
sim_dir     <- file.path(precomp_dir, "sim_mats")
dir.create(precomp_dir, recursive = TRUE, showWarnings = FALSE)
if (MODE == "full")  dir.create(sim_dir,    recursive = TRUE, showWarnings = FALSE)
if (MODE == "local") dir.create(SKETCH_DIR, recursive = TRUE, showWarnings = FALSE)

if (is.null(CHROM)) {
  cl <- Sys.getenv("CHROM_LIST", unset = "")
  CHROM_LIST <- if (nchar(cl) > 0) strsplit(cl, "[ ,]+")[[1]] else
    sprintf("C_gar_LG%02d", 1:28)
} else {
  CHROM_LIST <- CHROM
}

# RSpectra and RcppHNSW are loaded lazily so the script still runs (in
# degraded mode) if one is missing. Full mode requires RSpectra; local
# mode's sketch step requires RcppHNSW.
HAVE_RSPECTRA <- requireNamespace("RSpectra", quietly = TRUE)
HAVE_HNSW     <- requireNamespace("RcppHNSW", quietly = TRUE)

if (MODE == "full" && !HAVE_RSPECTRA) {
  warning("RSpectra not installed — full mode will fall back to cmdscale ",
          "(slow). install.packages('RSpectra') to enable fast MDS.")
}
if (MODE == "local" && !HAVE_HNSW) {
  warning("RcppHNSW not installed — local mode will skip sketch + sparse ",
          "edges (banded sim_mat still written). ",
          "install.packages('RcppHNSW') to enable Layer C.")
}

message(sprintf("[TR_C MDS] mode=%s, K_MDS=%d, SIM_BAND_HALF=%d, out=%s",
                MODE, K_MDS, SIM_BAND_HALF, precomp_dir))
if (MODE == "local") {
  message(sprintf("[TR_C MDS] sketch: d=%d, k_nn=%d, track=%s, cohort=%s, seed=%d, dir=%s",
                  D_SKETCH, K_NN, SKETCH_TRACK, COHORT_ID, RNG_SEED, SKETCH_DIR))
}
message("[TR_C MDS] chroms: ", length(CHROM_LIST))

# =============================================================================
# Helpers — BLAS sim_mat, fast cmdscale, sketch, sparse k-NN
# =============================================================================

# Build the n_win × n_win |cor|-similarity matrix from a samples × windows
# PC1 matrix in one BLAS call. NAs are imputed to column means before
# centering, so high-missingness windows contribute attenuated similarity
# (matches the v5 behavior on well-covered windows; differs on heavy-NA
# pairs — see ADR-0001 §6.1 validation threshold).
#
# Memory: peak ~3 × n_win² × 8 bytes during crossprod + abs + intermediate.
# At n_win = 6,000 (coarse): peak ~900 MB, fine on 16 GB nodes.
# At n_win = 30,000 (dense): peak ~22 GB, requires --mem=64G.
# At n_win > ~45,000: consider blocked tcrossprod (not implemented in v6;
# raise an issue if you hit this scale).
sim_mat_blas <- function(pc1_mat) {
  X <- pc1_mat                                   # n_samp × n_win
  col_means <- colMeans(X, na.rm = TRUE)
  col_means[!is.finite(col_means)] <- 0
  na_idx <- which(!is.finite(X), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) X[na_idx] <- col_means[na_idx[, 2L]]
  X <- sweep(X, 2L, col_means, FUN = "-")        # column-center
  col_norm <- sqrt(colSums(X * X))
  col_norm[col_norm < 1e-12] <- 1
  X <- sweep(X, 2L, col_norm, FUN = "/")         # column L2-normalize
  M <- abs(crossprod(X))                         # n_win × n_win, |cor|
  diag(M) <- 1
  M
}

# Pack a full n × n sim into banded representation (n × (2*band_half+1)).
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
    B <- pack_banded(M, band_half)
    rm(M); gc(verbose = FALSE)
    list(sim_mat = NULL, sim_band = B,
         format = paste0("banded_float32_pm", band_half),
         band_half = as.integer(band_half))
  } else {
    list(sim_mat = M, sim_band = NULL,
         format = "upper_triangle_float32", band_half = NA_integer_)
  }
}

# Reconstruct full N×N sim_mat from banded form (off-band cells = band median).
# v5-compatible behavior; only used when sim_mat is banded but MDS step
# needs a full distance matrix.
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
  diag(M) <- 1
  M
}

# Top-k MDS via RSpectra::eigs_sym on the double-centered Gram matrix.
# Falls back to cmdscale if RSpectra is unavailable.
fast_cmdscale <- function(dmat, k) {
  n <- nrow(dmat)
  if (!HAVE_RSPECTRA) {
    fit <- tryCatch(cmdscale(as.dist(dmat), k = k), error = function(e) NULL)
    return(fit)
  }
  D2 <- dmat * dmat
  rm <- rowMeans(D2); gm <- mean(D2)
  # B = -0.5 (D2 - row_mean - col_mean + grand_mean)
  B <- -0.5 * (D2 - outer(rm, rep(1, n)) - outer(rep(1, n), rm) + gm)
  B <- 0.5 * (B + t(B))    # symmetrize against FP drift
  eig <- tryCatch(RSpectra::eigs_sym(B, k = k, which = "LA"),
                  error = function(e) NULL)
  if (is.null(eig)) return(NULL)
  vals <- pmax(eig$values, 0)
  sweep(eig$vectors, 2L, sqrt(vals), FUN = "*")
}

# NN-smoothed similarity. Uses partial-sort top-k neighbor selection on the
# (n_win × n_win) distance matrix; then computes pairwise dist of the
# smoothed MDS coords. The smoothed MDS dist matrix is reused across scales
# by the caller (passed in via cached_dmat) when present.
make_nn_sim <- function(mds_mat, dmat, k_use) {
  n <- nrow(mds_mat)
  smoothed <- matrix(0, n, ncol(mds_mat))
  for (wi in seq_len(n)) {
    d <- dmat[wi, ]; d[wi] <- Inf
    # partial sort: indices of the k smallest distances
    cutoff <- sort(d, partial = k_use)[k_use]
    nn_idx <- which(d <= cutoff)
    if (length(nn_idx) > k_use) nn_idx <- nn_idx[seq_len(k_use)]
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

# -----------------------------------------------------------------------------
# Layer C: sketch + sparse k-NN graph
# -----------------------------------------------------------------------------

# Load or create the cohort-scoped sample sketch basis (one file per cohort,
# shared across chroms). Persisted at $SKETCH_DIR/sample_sketch_basis.rds.
load_or_create_basis <- function(sketch_dir, n_samp, sample_order,
                                 d_sketch, seed, cohort_id) {
  basis_path <- file.path(sketch_dir, "sample_sketch_basis.rds")
  if (file.exists(basis_path)) {
    b <- readRDS(basis_path)
    if (b$d_sketch == d_sketch && b$cohort_id == cohort_id &&
        b$seed == seed && nrow(b$basis) == n_samp &&
        identical(b$sample_order, sample_order)) {
      message("[TR_C MDS] sketch basis: reusing ", basis_path)
      return(b)
    } else {
      stop("[TR_C MDS] sketch basis at ", basis_path,
           " has mismatched config (cohort_id/seed/d_sketch/sample_order). ",
           "Delete it or change --cohort_id / --rng_seed.")
    }
  }
  set.seed(seed)
  R <- matrix(rnorm(n_samp * d_sketch, mean = 0, sd = 1 / sqrt(d_sketch)),
              nrow = n_samp, ncol = d_sketch)
  b <- list(
    basis           = R,
    sample_order    = sample_order,
    d_sketch        = as.integer(d_sketch),
    seed            = as.integer(seed),
    cohort_id       = cohort_id,
    created_at      = Sys.time(),
    toolkit_version = "inversion-popgen-toolkit (TR_C v6)"
  )
  saveRDS(b, basis_path)
  message("[TR_C MDS] sketch basis: created ", basis_path,
          " (n_samp=", n_samp, ", d_sketch=", d_sketch, ")")
  b
}

# Build the n_win × d_sketch sketch from a per-window per-sample feature
# matrix (n_samp × n_win). NA → column mean before projection. Rows of the
# resulting sketch are L2-normalized so cosine sim == inner product.
build_window_sketch <- function(feat_mat, basis_R) {
  X <- feat_mat
  col_means <- colMeans(X, na.rm = TRUE)
  col_means[!is.finite(col_means)] <- 0
  na_idx <- which(!is.finite(X), arr.ind = TRUE)
  na_count <- if (nrow(na_idx) > 0L) {
    tabulate(na_idx[, 2L], nbins = ncol(X))
  } else {
    integer(ncol(X))
  }
  if (nrow(na_idx) > 0L) X[na_idx] <- col_means[na_idx[, 2L]]
  # Project: t(X) is n_win × n_samp, basis_R is n_samp × d_sketch
  S <- crossprod(X, basis_R)                # n_win × d_sketch
  rn <- sqrt(rowSums(S * S))
  rn[rn < 1e-12] <- 1
  S <- sweep(S, 1L, rn, FUN = "/")          # L2-normalize rows
  list(sketch = S, na_count = na_count)
}

# Build sparse k-NN edge list using RcppHNSW. Cosine metric (== inner product
# on L2-normalized vectors). Returns a data.table with columns
# (window_a, window_b, rank, sim_cosine).
#
# NOTE on sign / anchor-flip: TR_B anchor-flips PC1 so signs are globally
# consistent within a chromosome. Two windows in the same regime should
# produce sketch vectors with the same sign. Windows in pca$unflipped_windows
# may have an inconsistent sign and their HNSW neighbors may be wrong; check
# the unflipped count in precomp_summary.tsv. The downstream |cor|-semantic
# is recovered via sim_cosine = |1 - dist|, but that only helps once the
# neighbor was returned by HNSW in the first place.
build_sparse_edges <- function(sketch_mat, k_nn, chrom) {
  n_win <- nrow(sketch_mat)
  k_query <- min(k_nn + 1L, n_win)            # +1 to drop self
  ann <- RcppHNSW::hnsw_knn(sketch_mat, k = k_query, distance = "cosine",
                            M = 16L, ef_construction = 200L, ef = 100L)
  # ann$idx: n_win × k_query (1-based indices into sketch_mat rows)
  # ann$dist: same shape, cosine distance in [0, 2] for unit-norm vectors
  idx  <- ann$idx
  dist <- ann$dist
  edges <- vector("list", n_win)
  for (wa in seq_len(n_win)) {
    nbrs  <- idx[wa, ]
    dists <- dist[wa, ]
    keep  <- nbrs != wa
    nbrs  <- nbrs[keep]
    dists <- dists[keep]
    if (length(nbrs) > k_nn) {
      nbrs  <- nbrs[seq_len(k_nn)]
      dists <- dists[seq_len(k_nn)]
    }
    if (length(nbrs) == 0L) next
    # RcppHNSW cosine distance = 1 - cos(theta). On L2-normed vectors,
    # cos(theta) ∈ [-1, 1]. We want |cor|-semantics (matching the BLAS
    # sim_mat convention), so sim = |cos| = |1 - dist|. This treats
    # anchor-flipped-opposite PC1 vectors as similar (same regime, opposite
    # sign), which is what the inversion atlas wants.
    sim_cos <- pmax(0, pmin(1, abs(1 - dists)))
    edges[[wa]] <- data.table(
      chrom        = chrom,
      window_a     = as.integer(wa - 1L),       # 0-based to match window_idx
      window_b     = as.integer(nbrs - 1L),
      rank         = seq_along(nbrs),
      sim_cosine   = round(sim_cos, 6),
      hop_distance = 1L
    )
  }
  rbindlist(edges, use.names = TRUE)
}

# Sketch-novelty: per window, 1 - max(sim_cosine) over its top-k neighbors
# that lie outside an exclusion zone of ±excl windows. High novelty = the
# window's profile is unusual vs nearby windows. This is the dense-scale
# replacement for "MDS outlier ribbon" (ADR-0001 §5.2).
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
all_window_dt <- list()
all_summary   <- list()
basis_cached  <- NULL   # lazily created once per run in local mode

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

  # ── sim_mat: BLAS path (v6) ───────────────────────────────────────────
  t_sim <- proc.time()
  simbox <- build_sim_mat(pcs_aligned[[1]], n_win, SIM_BAND_HALF,
                          SIM_N_FULL_THRESHOLD, force_banded = (MODE == "local"))
  message(sprintf("[TR_C MDS] %s: sim_mat (%s) %.1fs",
                  chrom, simbox$format, (proc.time() - t_sim)[3]))

  # ── full mode: MDS + NN sims ──────────────────────────────────────────
  if (MODE == "full") {
    sim_M <- sim_full(simbox, n_win)
    dmat  <- 1 - sim_M; diag(dmat) <- 0

    t_mds <- proc.time()
    mds_fit <- fast_cmdscale(dmat, K_MDS)
    if (is.null(mds_fit) || nrow(mds_fit) != n_win) {
      mds_mat <- matrix(NA_real_, n_win, K_MDS)
    } else {
      mds_mat <- mds_fit
      if (ncol(mds_mat) < K_MDS) {
        mds_mat <- cbind(mds_mat, matrix(NA_real_, n_win, K_MDS - ncol(mds_mat)))
      }
    }
    message(sprintf("[TR_C MDS] %s: fast_cmdscale k=%d %.1fs",
                    chrom, K_MDS, (proc.time() - t_mds)[3]))

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

    sketch_path <- NA_character_
    sparse_edges_path <- NA_character_
    sketch_novelty <- rep(NA_real_, n_win)
  } else {
    # ── local mode: skip global MDS; build sketch + sparse k-NN graph ─
    mds_mat   <- matrix(NA_real_, n_win, K_MDS)
    mds_z     <- matrix(NA_real_, n_win, K_MDS)
    max_abs_z <- rep(NA_real_, n_win)
    max_z_axis <- rep(NA_integer_, n_win)
    bg_q      <- rep(NA_real_, 6)
    sketch_novelty <- rep(NA_real_, n_win)
    sketch_path <- NA_character_
    sparse_edges_path <- NA_character_

    if (HAVE_HNSW) {
      # Lazy basis creation (once per run, shared across chroms)
      if (is.null(basis_cached)) {
        basis_cached <- load_or_create_basis(
          SKETCH_DIR, n_samp, sample_order, D_SKETCH, RNG_SEED, COHORT_ID)
      }
      # Choose the per-window feature matrix for the sketch
      feat_mat <- if (SKETCH_TRACK == "PC1") {
        pcs_aligned[[1]]
      } else {
        # theta_raw track: would require the raw theta_pi matrix from TR_A.
        # Not currently passed through TR_B's RDS — placeholder for now.
        stop("[TR_C MDS] --sketch_track theta_raw not yet wired; ",
             "TR_B v6 does not carry the raw theta_pi matrix. ",
             "Use --sketch_track PC1 for v6.")
      }
      t_sk <- proc.time()
      sk <- build_window_sketch(feat_mat, basis_cached$basis)
      message(sprintf("[TR_C MDS] %s: sketch (%d × %d, track=%s) %.1fs",
                      chrom, n_win, D_SKETCH, SKETCH_TRACK,
                      (proc.time() - t_sk)[3]))

      sketch_obj <- list(
        chrom        = chrom,
        sketch       = sk$sketch,
        sketch_track = SKETCH_TRACK,
        window_idx   = as.integer(win_grid$window_idx),
        d_sketch     = as.integer(D_SKETCH),
        basis_id     = paste0(COHORT_ID, "_seed", RNG_SEED, "_d", D_SKETCH),
        n_win        = as.integer(n_win),
        na_count     = sk$na_count
      )
      sketch_path <- file.path(SKETCH_DIR, sprintf("%s.sketch.rds", chrom))
      saveRDS(sketch_obj, sketch_path)

      t_gr <- proc.time()
      edges_dt <- build_sparse_edges(sk$sketch, K_NN, chrom)
      sparse_edges_path <- file.path(SKETCH_DIR,
                                     sprintf("%s.sparse_edges.tsv.gz", chrom))
      fwrite(edges_dt, sparse_edges_path, sep = "\t", compress = "gzip")
      message(sprintf("[TR_C MDS] %s: sparse_edges (k=%d, %d rows) %.1fs",
                      chrom, K_NN, nrow(edges_dt), (proc.time() - t_gr)[3]))

      sketch_novelty <- compute_sketch_novelty(edges_dt, n_win, excl = 5L)
    } else {
      message("[TR_C MDS] ", chrom,
              ": skipping sketch + sparse_edges (RcppHNSW not installed)")
    }
  }

  # ── per-window data.table ─────────────────────────────────────────────
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
    anchor_window_idx = anchor_idx - 1L,
    sketch_novelty   = round(sketch_novelty, 4)    # NA in full mode / no HNSW
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
    mode                     = MODE,
    sketch_path              = sketch_path,
    sparse_edges_path        = sparse_edges_path
  )
  saveRDS(precomp, file.path(precomp_dir, sprintf("%s.precomp.rds", chrom)))
  message(sprintf("[TR_C MDS] %s: precomp.rds saved (%d windows, NPC=%d, K_MDS=%d, mode=%s)",
                  chrom, n_win, NPC, K_MDS, MODE))

  all_window_dt[[chrom]] <- dt[, .(chrom, window_idx, start_bp, end_bp, mid_bp,
                                    theta_pi_median, theta_z_direct,
                                    max_abs_z, max_z_axis, lambda_ratio,
                                    MDS1, MDS2, sketch_novelty)]
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
    median_novelty   = round(median(sketch_novelty, na.rm = TRUE), 4),
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
