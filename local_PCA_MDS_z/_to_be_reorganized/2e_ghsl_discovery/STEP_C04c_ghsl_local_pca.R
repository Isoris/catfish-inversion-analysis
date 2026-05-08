#!/usr/bin/env Rscript
# =============================================================================
# STEP_C04c_ghsl_local_pca.R
# =============================================================================
# Phase 2 / 2e_ghsl_discovery — Layer C local PCA companion to STEP_C04
# (heavy divergence engine) and STEP_C04b (light classifier).
#
# Purpose
# -------
# Page 3 ("local PCA GHSL") in the atlas mirrors page 1 (dosage local PCA |z|)
# and page 12 (local PCA θπ). It needs window-by-window local PCA outputs
# computed from the GHSL divergence matrix:
#   * pc1 / pc2 loadings per sample per window
#   * lambda_1 / lambda_2 / lambda_ratio (1D-ness indicator)
#   * z_profile (chromosome-baseline robust |Z| of per-sample population deviations)
#   * sim_mat (window × window similarity from |cor(pc1_i, pc1_j)|)
#   * mds_coords (cmdscale of 1 - sim_mat to MDS1, MDS2)
#   * sign-aligned loadings (anchor-window flipping for cursor-coherent display)
#   * L2 / L1 envelopes from contiguous high-|Z| runs
#
# This is the cluster-side precomp; turn-3 JSON exporter consumes the RDS
# this script writes.
#
# Architecture decisions (locked in phase2_ghsl_arrangement_v1.md and
# schema_v2_addendum_ghsl_page3.md)
# ------------------------------------------------------------------
#  1. Local PCA input is the RAW div_mat (5kb base scale), NOT the rolling-
#     smoothed div_roll matrices. Cross-sample averaging across 226 samples
#     reduces covariance noise by ~1/sqrt(226·pad) ≈ 0.04, well below
#     biological signal. Smoothed input would create artificial autocorrelation
#     between adjacent windows that isn't in the underlying biology.
#     A `--smoothing-scale s50` option falls back to rolling input if raw
#     proves too noisy on real LG28 data — empirically reversible.
#
#  2. Heteroscedastic weighting: each sample's contribution to a window's
#     local-PCA covariance is weighted by sqrt(n_phased_het / median).
#     Samples with sparse phased calls in this window contribute less.
#     Standard heteroscedastic-PCA practice.
#
#  3. Sign alignment: per-window PC1/PC2 are sign-ambiguous from SVD. After
#     local PCA, the highest-|Z| window is picked as the anchor; every
#     other window's PC1 is flipped if its correlation with the anchor's
#     PC1 is negative. Same for PC2 (against anchor's PC2, after PC1 align).
#     Both raw and aligned loadings are emitted so consumers can pick.
#
#  4. Sim_mat is dense (no banding, no quantization) at GHSL's 5-kb base.
#     LG28 = ~4,300 windows → sim_mat = 4,300² × 4 bytes = ~74 MB raw,
#     ~18 MB gzip. Within the 120 MB-per-chromosome JSON budget.
#     Stored as upper triangle to save half the bytes (turn-3 exporter
#     does the upper-triangle pack; this script keeps full symmetric in RDS).
#
#  5. MDS via cmdscale(1 - sim_mat, k = 2). Each point in MDS space = one
#     window's geometric fingerprint; the page-3 #ghMDSPanel scatters them.
#
#  6. Secondary envelopes (L2 / L1 from contiguous high-|Z| runs in z_profile)
#     ARE emitted, but as a SECONDARY layer. The primary GHSL candidate set
#     is STEP_C04b's PASS-status runs (calibrated, denominator-confound-
#     mitigated). Both ship to the page-3 atlas; the user sees both, can
#     compare where they agree vs disagree.
#
#     Asymmetry vs θπ: STEP_TR_B_classify_theta.R DOES emit primary L2/L1
#     envelopes from |Z| profile because no upstream production θπ
#     candidate detector exists (no θπ-equivalent of STEP_C04b is yet
#     calibrated). For θπ, the |Z|-derived envelopes ARE the candidate
#     set. For GHSL, they're a cross-check on top of the production set.
#     This isn't a contradiction — it reflects which streams have
#     calibrated upstream detectors, and which don't (yet).
#
# Inputs
# ------
#   <ghsl_v6_dir>/<CHROM>.ghsl_v6_matrices.rds  (from STEP_C04 v6 heavy)
#     Carries: div_mat, het_mat, n_sites_mat, n_phased_het_mat, rolling[],
#              rolling_het[], window_info, sample_names, chrom, params
#
# Output
# ------
#   <out_dir>/<CHROM>.ghsl_v6_localpca.rds
#     Carries: pc1_loadings_mat, pc2_loadings_mat (raw, sign-ambiguous)
#              pc1_loadings_aligned_mat, pc2_loadings_aligned_mat
#              lambda_1, lambda_2, lambda_ratio, z_profile, z_top10_mean
#              sim_mat (full symmetric float64 matrix)
#              mds1, mds2
#              anchor_window_idx
#              secondary_l2_envelopes, secondary_l1_envelopes (data.tables)
#              window_info, sample_names, chrom, params
#
#   PRIMARY envelopes (the canonical GHSL candidate set) are NOT in this
#   output — turn-3 exporter reads them from <chr>.ghsl_v6.annot.rds (the
#   production STEP_C04b PASS-runs), packs them into the page-3 JSON as
#   `ghsl_envelopes`. The SECONDARY ones from this script ship as
#   `ghsl_secondary_envelopes` — orthogonal cross-check.
#
# Wall time: ~5–10 min/chrom (~4,300 windows × 226 samples × pad=1 SVD).
#
# Usage
# -----
#   Rscript STEP_C04c_ghsl_local_pca.R \
#     <ghsl_v6_matrices_dir> <localpca_out_dir> \
#     [--chrom C_gar_LG28] \
#     [--pad 1] \
#     [--smoothing-scale none|s10|s20|s50|s100] \
#     [--anchor auto|<window_idx>] \
#     [--min-samples-per-window 50] \
#     [--z-threshold 2.5] \
#     [--min-l2-windows 5] \
#     [--merge-gap 3]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# Argument parsing
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(paste(
    "Usage: Rscript STEP_C04c_ghsl_local_pca.R <matrices_dir> <out_dir> [opts]",
    "  --chrom <chr>                 process single chromosome (default: all)",
    "  --pad 1                       local-PCA neighbourhood half-width (windows each side)",
    "  --smoothing-scale none        input matrix: 'none' for raw div_mat, or s10/s20/s50/s100",
    "  --anchor auto                 sign-alignment anchor window: auto = argmax(z_profile)",
    "  --min-samples-per-window 50   skip windows with fewer non-NA samples",
    "  --z-threshold 2.5             |Z| threshold for SECONDARY envelope seeding",
    "  --min-l2-windows 5            minimum contiguous-window run length for secondary L2",
    "  --merge-gap 3                 merge adjacent secondary L2 envelopes within this gap",
    "",
    "  Note: secondary envelopes are an orthogonal local-PCA cross-check.",
    "  Primary GHSL candidates come from STEP_C04b's PASS-status runs",
    "  (calibrated, denominator-confound-mitigated). The two layers ship",
    "  side-by-side in the page-3 atlas JSON; the user sees both.",
    sep = "\n"
  ))
}

MATRICES_DIR <- args[1]
OUT_DIR      <- args[2]

CHROM_FILTER <- NULL
PAD          <- 1L
SMOOTHING    <- "none"
ANCHOR_ARG   <- "auto"
MIN_SAMPLES  <- 50L
Z_THRESHOLD  <- 2.5
MIN_L2_WIN   <- 5L
MERGE_GAP    <- 3L

i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--chrom"                  && i < length(args)) { CHROM_FILTER <- args[i + 1]; i <- i + 2L }
  else if (a == "--pad"                    && i < length(args)) { PAD          <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--smoothing-scale"        && i < length(args)) { SMOOTHING    <- args[i + 1]; i <- i + 2L }
  else if (a == "--anchor"                 && i < length(args)) { ANCHOR_ARG   <- args[i + 1]; i <- i + 2L }
  else if (a == "--min-samples-per-window" && i < length(args)) { MIN_SAMPLES  <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--z-threshold"            && i < length(args)) { Z_THRESHOLD  <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--min-l2-windows"         && i < length(args)) { MIN_L2_WIN   <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--merge-gap"              && i < length(args)) { MERGE_GAP    <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}

stopifnot(SMOOTHING %in% c("none", "s10", "s20", "s30", "s40", "s50", "s100"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("================================================================")
message("[S3v6c] GHSL local PCA — turn-2 cluster-side precomp")
message("================================================================")
message("[S3v6c] Matrices dir: ", MATRICES_DIR)
message("[S3v6c] Output dir:   ", OUT_DIR)
message("[S3v6c] pad=",            PAD,
        " smoothing=",              SMOOTHING,
        " min_samples=",            MIN_SAMPLES,
        " anchor=",                 ANCHOR_ARG)
message("[S3v6c] secondary envelope params: z>", Z_THRESHOLD,
        " min_l2=", MIN_L2_WIN,
        " merge_gap=", MERGE_GAP)
message("[S3v6c] (primary GHSL candidates: STEP_C04b PASS-runs;",
        " these are local-PCA secondary cross-check)")

# =============================================================================
# Helpers
# =============================================================================

# ---- Per-window local PCA on a 226 × (2·pad+1) neighbourhood -----------------
# Returns a list with pc1, pc2, lambda_1, lambda_2, n_used.
# Uses heteroscedastic weighting: each sample's row in the centred matrix is
# scaled by sqrt(n_phased_het[s] / median(n_phased_het)) at the focal window.
local_pca_one_window <- function(div_block, weight_vec, min_samples = 50L) {
  # div_block: n_samples × (2·pad+1) numeric matrix; rows are samples,
  #            cols are windows in the focal neighbourhood
  # weight_vec: length n_samples vector of sqrt(n_phased_het / median)
  ok <- complete.cases(div_block) & is.finite(weight_vec) & weight_vec > 0
  if (sum(ok) < min_samples || sum(ok) < ncol(div_block) + 2L) {
    return(list(pc1 = NULL, pc2 = NULL, lambda_1 = NA_real_,
                lambda_2 = NA_real_, n_used = sum(ok)))
  }
  block_ok <- div_block[ok, , drop = FALSE]
  w        <- weight_vec[ok]

  # Apply heteroscedastic weighting (multiply each row by its weight) then
  # column-centre (so the PCA captures sample-vs-sample structure, not
  # mean-window structure).
  block_w  <- block_ok * w
  centred  <- sweep(block_w, 2, colMeans(block_w), FUN = "-")

  sv <- tryCatch(svd(centred, nu = 2, nv = 0), error = function(e) NULL)
  if (is.null(sv) || length(sv$d) < 1) {
    return(list(pc1 = NULL, pc2 = NULL, lambda_1 = NA_real_,
                lambda_2 = NA_real_, n_used = sum(ok)))
  }

  # Scatter eigenvectors back to full sample dimension (NA for rows we skipped)
  n_samp <- nrow(div_block)
  pc1_full <- rep(NA_real_, n_samp); pc1_full[ok] <- as.numeric(sv$u[, 1])
  pc2_full <- rep(NA_real_, n_samp)
  if (ncol(sv$u) >= 2) pc2_full[ok] <- as.numeric(sv$u[, 2])

  d2 <- sv$d ^ 2
  list(
    pc1      = pc1_full,
    pc2      = pc2_full,
    lambda_1 = if (length(d2) >= 1) d2[1] else NA_real_,
    lambda_2 = if (length(d2) >= 2) d2[2] else NA_real_,
    n_used   = sum(ok)
  )
}

# ---- Robust |Z| profile from per-window per-sample deviations ----------------
# For each window w, compute robust |z| of per-sample deviation from chromosome
# pattern. This is sign-invariant by construction (uses absolute values), so
# it does NOT need sim_mat machinery — a 1D scan, n_windows long.
compute_z_profile <- function(input_mat) {
  # input_mat: n_samples × n_windows numeric matrix
  n_win <- ncol(input_mat)
  win_median <- apply(input_mat, 2, function(v) median(v, na.rm = TRUE))
  dev_mat    <- sweep(input_mat, 2, win_median, FUN = "-")

  all_devs <- as.vector(dev_mat); all_devs <- all_devs[is.finite(all_devs)]
  scale_mad <- if (length(all_devs) > 100) mad(all_devs, constant = 1.4826) else NA_real_
  if (!is.finite(scale_mad) || scale_mad <= 0) scale_mad <- 1e-6

  # Per-window summaries
  z_max <- apply(dev_mat, 2, function(v) {
    v <- v[is.finite(v)]
    if (length(v) < MIN_SAMPLES) return(NA_real_)
    max(abs(v / scale_mad))
  })
  z_top10 <- apply(dev_mat, 2, function(v) {
    v <- abs(v[is.finite(v)] / scale_mad)
    if (length(v) < MIN_SAMPLES) return(NA_real_)
    k <- max(1L, ceiling(0.1 * length(v)))
    mean(sort(v, decreasing = TRUE)[seq_len(k)])
  })
  list(z_max = z_max, z_top10 = z_top10, scale_mad = scale_mad)
}

# ---- Sign alignment of per-window PC vectors to an anchor window -------------
# After SVD, each window's PC1 vector has arbitrary sign. To make the lines
# panel render coherently across windows, flip each window's PC1 if its
# correlation with the anchor's PC1 is negative. Then do PC2 against the
# anchor's PC2 (after PC1 alignment).
sign_align_loadings <- function(pc1_mat, pc2_mat, anchor_idx) {
  n_win <- ncol(pc1_mat)
  pc1_aligned <- pc1_mat
  pc2_aligned <- pc2_mat

  anchor_pc1 <- pc1_mat[, anchor_idx]
  anchor_pc2 <- pc2_mat[, anchor_idx]

  for (w in seq_len(n_win)) {
    if (w == anchor_idx) next
    # PC1 alignment
    v1 <- pc1_mat[, w]
    valid1 <- is.finite(v1) & is.finite(anchor_pc1)
    if (sum(valid1) >= MIN_SAMPLES) {
      r1 <- suppressWarnings(cor(v1[valid1], anchor_pc1[valid1]))
      if (is.finite(r1) && r1 < 0) pc1_aligned[, w] <- -v1
    }
    # PC2 alignment (against anchor's PC2)
    v2 <- pc2_mat[, w]
    valid2 <- is.finite(v2) & is.finite(anchor_pc2)
    if (sum(valid2) >= MIN_SAMPLES) {
      r2 <- suppressWarnings(cor(v2[valid2], anchor_pc2[valid2]))
      if (is.finite(r2) && r2 < 0) pc2_aligned[, w] <- -v2
    }
  }

  list(pc1_aligned = pc1_aligned, pc2_aligned = pc2_aligned)
}

# ---- Sim_mat: window × window similarity from |cor(pc1_i, pc1_j)| ------------
# Sign-invariant by construction (absolute value of correlation). Uses raw
# pc1_loadings (not aligned), since |cor| doesn't depend on sign.
compute_sim_mat <- function(pc1_mat) {
  n_win <- ncol(pc1_mat)
  # cor() column-wise, NA-aware
  sim <- suppressWarnings(abs(cor(pc1_mat, use = "pairwise.complete.obs")))
  # cor() returns NA for columns with no overlap — replace with 0 for safety
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  storage.mode(sim) <- "double"
  sim
}

# ---- MDS via cmdscale on (1 - sim) -------------------------------------------
compute_mds <- function(sim_mat) {
  d <- 1 - sim_mat
  # Numerical safety: enforce non-negative, symmetric, zero diagonal
  d[d < 0] <- 0
  d <- (d + t(d)) / 2
  diag(d) <- 0
  mds <- tryCatch(cmdscale(d, k = 2, eig = FALSE), error = function(e) NULL)
  if (is.null(mds) || ncol(mds) < 2) {
    return(list(mds1 = rep(NA_real_, nrow(d)), mds2 = rep(NA_real_, nrow(d))))
  }
  list(mds1 = mds[, 1], mds2 = mds[, 2])
}

# ---- Secondary envelope helpers ---------------------------------------------
# STEP_C04c's local-PCA-derived envelopes are SECONDARY: the production
# STEP_C04b classifier (calibrated, denominator-confound-mitigated) is the
# primary GHSL candidate detector. We compute these alongside as an
# orthogonal cross-check — when the two detectors agree, that's strong
# convergent evidence; when they disagree, that's diagnostic information
# (local-PCA sees a peak v6 doesn't flag = potential discovery, or potential
# noise needing review). The atlas can render both layers' strips on the
# |Z| panel, distinguished by color.
#
# This mirrors what STEP_TR_B does for θπ — except θπ has no upstream
# production detector (no equivalent of STEP_C04b exists for nucleotide
# diversity), so its envelopes ARE the candidate set. GHSL is the inverse
# situation: production detector exists upstream, so local-PCA envelopes
# are secondary.

call_envelopes <- function(z_vec, threshold, min_run, merge_gap) {
  flag <- !is.na(z_vec) & z_vec > threshold
  if (!any(flag)) return(data.table())
  runs <- rle(flag)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1) + 1L)
  hits <- which(runs$values)
  envs <- data.table(
    win_start = starts[hits],
    win_end   = ends[hits],
    n_windows = runs$lengths[hits]
  )
  envs <- envs[n_windows >= min_run]
  if (nrow(envs) == 0) return(envs)

  setorder(envs, win_start)
  merged <- envs[1]
  for (k in seq.int(2L, nrow(envs))) {
    if (envs$win_start[k] - merged$win_end[nrow(merged)] <= merge_gap) {
      merged[nrow(merged), win_end   := envs$win_end[k]]
      merged[nrow(merged), n_windows := win_end - win_start + 1L]
    } else {
      merged <- rbind(merged, envs[k])
    }
  }
  merged
}

merge_to_l1 <- function(l2_dt, l1_merge_gap) {
  if (nrow(l2_dt) == 0) return(data.table())
  if (nrow(l2_dt) == 1) {
    return(l2_dt[, .(win_start, win_end, n_l2 = 1L)])
  }
  setorder(l2_dt, win_start)
  merged <- l2_dt[1, .(win_start, win_end, n_l2 = 1L)]
  for (k in seq.int(2L, nrow(l2_dt))) {
    if (l2_dt$win_start[k] - merged$win_end[nrow(merged)] <= l1_merge_gap) {
      merged[nrow(merged), win_end := l2_dt$win_end[k]]
      merged[nrow(merged), n_l2 := n_l2 + 1L]
    } else {
      merged <- rbind(merged, l2_dt[k, .(win_start, win_end, n_l2 = 1L)])
    }
  }
  merged
}

# =============================================================================
# Main loop — one chromosome at a time
# =============================================================================

rds_files <- sort(list.files(MATRICES_DIR, pattern = "\\.ghsl_v6_matrices\\.rds$",
                             full.names = TRUE))
if (length(rds_files) == 0) {
  stop("[S3v6c] FATAL: No .ghsl_v6_matrices.rds files in: ", MATRICES_DIR)
}
message("[S3v6c] Found ", length(rds_files), " matrix files")

if (!is.null(CHROM_FILTER)) {
  rds_files <- rds_files[grepl(paste0("/", CHROM_FILTER, "\\."), rds_files)]
  if (length(rds_files) == 0) {
    stop("[S3v6c] FATAL: --chrom ", CHROM_FILTER, " not found among matrix files")
  }
}

for (rds_f in rds_files) {
  t0 <- proc.time()
  message("\n----------------------------------------------------------------")
  message("[S3v6c] Loading: ", basename(rds_f))

  m <- readRDS(rds_f)
  chr     <- m$chrom
  win_inf <- m$window_info
  snames  <- m$sample_names
  n_samp  <- length(snames)

  # Pick input matrix per --smoothing-scale
  input_mat <- if (SMOOTHING == "none") {
    m$div_mat
  } else {
    if (is.null(m$rolling[[SMOOTHING]])) {
      stop("[S3v6c] FATAL: rolling scale ", SMOOTHING,
           " not in matrices RDS for ", chr,
           "\n  Available: ", paste(names(m$rolling), collapse = ", "))
    }
    m$rolling[[SMOOTHING]]
  }
  n_win <- ncol(input_mat)

  # Heteroscedastic weights from n_phased_het_mat (per-window per-sample
  # phased-het count). Weight = sqrt(n_phased_het / median across all cells).
  # The median across all cells is the natural "typical" denominator —
  # samples and windows with ≥ median get full weight, those below scale down.
  nph <- m$n_phased_het_mat
  med_nph <- median(nph[is.finite(nph) & nph > 0], na.rm = TRUE)
  if (!is.finite(med_nph) || med_nph <= 0) med_nph <- 1
  weight_mat <- sqrt(pmax(nph, 0) / med_nph)
  weight_mat[!is.finite(weight_mat)] <- 0

  message("[S3v6c] ", chr, ": ", n_samp, " samples × ", n_win, " windows")
  message("[S3v6c]   smoothing=", SMOOTHING, "  median(n_phased_het)=", med_nph)

  # ── Per-window local PCA ──
  message("[S3v6c] Stage 1: per-window local PCA (pad=", PAD, ")...")
  t1 <- proc.time()

  pc1_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                    dimnames = list(snames, NULL))
  pc2_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                    dimnames = list(snames, NULL))
  lambda_1 <- rep(NA_real_, n_win)
  lambda_2 <- rep(NA_real_, n_win)
  n_used   <- rep(0L,        n_win)

  for (w in seq_len(n_win)) {
    lo <- max(1L, w - PAD); hi <- min(n_win, w + PAD)
    block <- input_mat[, lo:hi, drop = FALSE]
    # Use weight from the focal window only — heteroscedasticity is
    # primarily about the focal sample's phased-het coverage at this position
    weights <- weight_mat[, w]
    res <- local_pca_one_window(block, weights, min_samples = MIN_SAMPLES)
    if (!is.null(res$pc1)) pc1_mat[, w] <- res$pc1
    if (!is.null(res$pc2)) pc2_mat[, w] <- res$pc2
    lambda_1[w] <- res$lambda_1
    lambda_2[w] <- res$lambda_2
    n_used[w]   <- res$n_used

    if (w %% 1000 == 0) {
      message("[S3v6c]   window ", w, " / ", n_win,
              " (", round(100 * w / n_win, 1), "%)")
    }
  }

  lambda_ratio <- ifelse(is.finite(lambda_2) & lambda_2 > 0,
                         lambda_1 / lambda_2, NA_real_)

  message("[S3v6c]   local PCA in ", round((proc.time() - t1)[3], 1), "s")
  message("[S3v6c]   λ₁/λ₂ ratio: median=",
          round(median(lambda_ratio, na.rm = TRUE), 2),
          " p95=", round(quantile(lambda_ratio, 0.95, na.rm = TRUE), 2))

  # ── Robust |Z| profile (1D, sign-invariant) ──
  message("[S3v6c] Stage 2: robust |Z| profile...")
  z_res <- compute_z_profile(input_mat)
  z_profile    <- z_res$z_max
  z_top10_mean <- z_res$z_top10
  message("[S3v6c]   |Z|: median=", round(median(z_profile, na.rm = TRUE), 2),
          " p95=", round(quantile(z_profile, 0.95, na.rm = TRUE), 2),
          " p99=", round(quantile(z_profile, 0.99, na.rm = TRUE), 2))

  # ── Sign alignment of loadings ──
  message("[S3v6c] Stage 3: sign-aligning loadings to anchor window...")
  if (ANCHOR_ARG == "auto") {
    # argmax(z_profile), tie-broken by lowest index, ignoring NAs
    z_for_anchor <- z_profile
    z_for_anchor[!is.finite(z_for_anchor)] <- -Inf
    anchor_idx <- which.max(z_for_anchor)
  } else {
    anchor_idx <- as.integer(ANCHOR_ARG)
    if (is.na(anchor_idx) || anchor_idx < 1L || anchor_idx > n_win) {
      stop("[S3v6c] FATAL: --anchor ", ANCHOR_ARG, " out of range 1..", n_win)
    }
  }
  message("[S3v6c]   anchor window = ", anchor_idx,
          "  (z_profile=", round(z_profile[anchor_idx], 3), ")")

  aligned <- sign_align_loadings(pc1_mat, pc2_mat, anchor_idx)
  pc1_aligned_mat <- aligned$pc1_aligned
  pc2_aligned_mat <- aligned$pc2_aligned

  # ── Sim_mat ──
  message("[S3v6c] Stage 4: window×window sim_mat (|cor(pc1)| dense)...")
  t4 <- proc.time()
  sim_mat <- compute_sim_mat(pc1_mat)
  size_mb <- round(object.size(sim_mat) / 1024 / 1024, 1)
  message("[S3v6c]   sim_mat: ", n_win, "×", n_win, " (", size_mb,
          " MB in memory)  in ", round((proc.time() - t4)[3], 1), "s")

  # ── MDS coords ──
  message("[S3v6c] Stage 5: cmdscale → MDS1, MDS2...")
  mds <- compute_mds(sim_mat)
  message("[S3v6c]   MDS1 range = [",
          round(min(mds$mds1, na.rm = TRUE), 3), ", ",
          round(max(mds$mds1, na.rm = TRUE), 3), "]   MDS2 range = [",
          round(min(mds$mds2, na.rm = TRUE), 3), ", ",
          round(max(mds$mds2, na.rm = TRUE), 3), "]")

  # ── Stage 6: secondary envelopes from local-PCA z_profile ──
  # These are SECONDARY to the production STEP_C04b PASS-status runs.
  # Keep them: when they agree with the v6 classifier's calls, that's
  # convergent evidence; when they disagree, that's diagnostic. Cheap
  # to compute, near-zero cost to ship even if junk on real data
  # (atlas just doesn't render them).
  message("[S3v6c] Stage 6: secondary envelopes (z>", Z_THRESHOLD,
          ", min=", MIN_L2_WIN, ", gap=", MERGE_GAP, ")...")
  sec_l2 <- call_envelopes(z_profile, Z_THRESHOLD, MIN_L2_WIN, MERGE_GAP)
  if (nrow(sec_l2) > 0) {
    sec_l2[, `:=`(
      start_bp = win_inf$start_bp[win_start],
      end_bp   = win_inf$end_bp[win_end],
      span_kb  = round((win_inf$end_bp[win_end] -
                        win_inf$start_bp[win_start]) / 1000, 1),
      peak_z   = vapply(seq_len(.N),
                        function(k) max(z_profile[win_start[k]:win_end[k]],
                                        na.rm = TRUE),
                        numeric(1)),
      mean_z   = vapply(seq_len(.N),
                        function(k) mean(z_profile[win_start[k]:win_end[k]],
                                         na.rm = TRUE),
                        numeric(1)),
      sec_id   = paste0(chr, "_ghsl_sec_", sprintf("%03d", seq_len(.N)))
    )]
  }
  l1_merge_gap <- 3L * MERGE_GAP
  sec_l1 <- merge_to_l1(sec_l2, l1_merge_gap)
  if (nrow(sec_l1) > 0) {
    sec_l1[, `:=`(
      start_bp = win_inf$start_bp[win_start],
      end_bp   = win_inf$end_bp[win_end],
      span_kb  = round((win_inf$end_bp[win_end] -
                        win_inf$start_bp[win_start]) / 1000, 1),
      sec_id   = paste0(chr, "_ghsl_sec_L1_", sprintf("%03d", seq_len(.N)))
    )]
  }
  message("[S3v6c]   secondary L2: ", nrow(sec_l2), "  secondary L1: ", nrow(sec_l1))
  if (nrow(sec_l2) > 0) {
    message("[S3v6c]   sec L2 span range: ",
            round(min(sec_l2$span_kb), 0), "–",
            round(max(sec_l2$span_kb), 0), " kb")
  }
  message("[S3v6c]   (primary GHSL candidates come from STEP_C04b PASS-runs;",
          " these are an orthogonal local-PCA cross-check)")

  # ── Save ──
  out_rds <- file.path(OUT_DIR, paste0(chr, ".ghsl_v6_localpca.rds"))
  saveRDS(list(
    chrom                    = chr,
    sample_names             = snames,
    window_info              = win_inf,
    n_samples                = n_samp,
    n_windows                = n_win,
    pc1_loadings_mat         = pc1_mat,
    pc2_loadings_mat         = pc2_mat,
    pc1_loadings_aligned_mat = pc1_aligned_mat,
    pc2_loadings_aligned_mat = pc2_aligned_mat,
    lambda_1                 = lambda_1,
    lambda_2                 = lambda_2,
    lambda_ratio             = lambda_ratio,
    n_used_per_window        = n_used,
    z_profile                = z_profile,
    z_top10_mean             = z_top10_mean,
    z_scale_mad              = z_res$scale_mad,
    sim_mat                  = sim_mat,
    mds1                     = mds$mds1,
    mds2                     = mds$mds2,
    anchor_window_idx        = anchor_idx,
    secondary_l2_envelopes   = sec_l2,
    secondary_l1_envelopes   = sec_l1,
    params                   = list(
      pad                    = PAD,
      smoothing_scale        = SMOOTHING,
      anchor_strategy        = ANCHOR_ARG,
      min_samples_per_window = MIN_SAMPLES,
      z_threshold            = Z_THRESHOLD,
      min_l2_windows         = MIN_L2_WIN,
      merge_gap              = MERGE_GAP,
      l1_merge_gap           = l1_merge_gap,
      generated_at           = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      generator              = "STEP_C04c_ghsl_local_pca.R",
      primary_envelopes      = "STEP_C04b annot RDS (ghsl_v6_status == 'PASS') — turn-3 sources",
      secondary_envelopes    = "this script — local-PCA-derived |Z| threshold, orthogonal cross-check"
    )
  ), out_rds)

  fsize <- round(file.info(out_rds)$size / 1024 / 1024, 1)
  message("[S3v6c]   saved: ", out_rds, " (", fsize, " MB on disk)")
  message("[S3v6c] ", chr, " DONE in ", round((proc.time() - t0)[3], 1), "s")

  rm(m, input_mat, weight_mat, pc1_mat, pc2_mat,
     pc1_aligned_mat, pc2_aligned_mat, sim_mat)
  invisible(gc(verbose = FALSE, full = FALSE))
}

message("\n================================================================")
message("[DONE] Turn 2: GHSL local PCA precomp complete")
message("================================================================")
message("  Output: ", OUT_DIR, "/<chr>.ghsl_v6_localpca.rds")
message("  Next:   turn 3 — JSON exporter v3 reads these RDS files and")
message("          emits the ghsl_local_pca layer in <chr>_phase2_ghsl.json")
