#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_B_classify_theta.R
# =============================================================================
# Phase 2 / 2f_theta_discovery — light classifier and JSON emitter.
#
# Consumes the long-format TSV produced by STEP_TR_A_compute_theta_matrices.R
# and produces:
#   1. A samples × windows θπ matrix in memory.
#   2. Per-window population metrics: median, MAD, IQR, |Z| from per-sample
#      deviations against the chromosome baseline.
#   3. Per-window local PCA over a (2·pad + 1) window neighbourhood: PC1 / PC2
#      loadings per sample, λ₁/λ₂ ratio, and a one-dimensionality flag.
#   4. L2 envelopes from contiguous runs of |Z| > threshold (configurable,
#      default 2.5) with minimum length and merge gap parameters.
#   5. Optional per-interval karyotype-class assignment via k-means on
#      interval-mean θπ when --intervals is supplied (k = 2..max_k with
#      silhouette selection, mirroring the GHSL v6b Part C output).
#   6. Page-12 atlas JSON layer set:
#        - theta_pi_per_window      (samples × windows θπ + n_sites)
#        - theta_pi_local_pca       (PC1/PC2 loadings, λ ratios, |Z| profile)
#        - theta_pi_envelopes       (L2 + L1 envelope coordinates)
#        - tracks                   (theta_pi_median, theta_pi_z, lambda_ratio)
#
# Architectural notes
# -------------------
# This pipeline mirrors the dosage scrubber's L1/L2 envelope detection but
# operates on a 1D |Z| profile rather than an n × n similarity matrix. The
# dosage scrubber computes a window×window sim_mat to recover sign-invariant
# similarity from sign-ambiguous PC eigenvectors; the θπ engine doesn't
# need this step because per-sample θπ values are themselves sign-stable
# (higher diversity = larger value, no eigenvector flip ambiguity).
# Envelopes are therefore detected from contiguous high-|Z| runs in 1D —
# the same algorithmic spirit, no n² intermediate, no browser memory cap.
#
# Usage
# -----
#   source 00_theta_config.sh
#   Rscript STEP_TR_B_classify_theta.R --chrom <CHROM> [--out-json <path>]
#                                                      [--intervals <bed.gz>]
#                                                      [--z-threshold 2.5]
#                                                      [--min-l2-windows 5]
#                                                      [--merge-gap 3]
#                                                      [--pad 1]
#                                                      [--max-k 5]
#
# Walltime: ~30 seconds per chromosome at win10000.step2000 (16,500
#           windows). Iterate freely on threshold tuning.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# Argument parsing
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

# B01-style: all args named. CHROM is required; rest have defaults.
CHROM         <- NULL
OUT_JSON      <- NULL
INTERVAL_FILE <- NULL
Z_THRESHOLD   <- 2.5
MIN_L2_WIN    <- 5L
MERGE_GAP     <- 3L
PAD           <- 1L
MAX_K         <- 5L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom"            && i < length(args)) { CHROM         <- args[i + 1]; i <- i + 2L }
  else if (a == "--out-json"    && i < length(args)) { OUT_JSON      <- args[i + 1]; i <- i + 2L }
  else if (a == "--intervals"   && i < length(args)) { INTERVAL_FILE <- args[i + 1]; i <- i + 2L }
  else if (a == "--z-threshold" && i < length(args)) { Z_THRESHOLD   <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--min-l2-windows" && i < length(args)) { MIN_L2_WIN <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--merge-gap"   && i < length(args)) { MERGE_GAP     <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--pad"         && i < length(args)) { PAD           <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--max-k"       && i < length(args)) { MAX_K         <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}

if (is.null(CHROM)) {
  stop("Usage: Rscript STEP_TR_B_classify_theta.R --chrom <CHROM> [opts]")
}

# Configuration from 00_theta_config.sh
THETA_TSV_DIR   <- Sys.getenv("THETA_TSV_DIR",   unset = NA)
PESTPG_SCALE    <- Sys.getenv("PESTPG_SCALE",    unset = "win10000.step2000")
THETA_GRID_MODE <- Sys.getenv("THETA_GRID_MODE", unset = "native")
JSON_OUT_DIR    <- Sys.getenv("JSON_OUT_DIR",    unset = NA)

stopifnot(!is.na(THETA_TSV_DIR))
stopifnot(THETA_GRID_MODE %in% c("native", "dosage"))

# Resolve input TSV path
in_tsv_basename <- if (THETA_GRID_MODE == "native") {
  sprintf("theta_native.%s.%s.tsv.gz", CHROM, PESTPG_SCALE)
} else {
  sprintf("theta_dgrid.%s.tsv.gz", CHROM)
}
in_tsv <- file.path(THETA_TSV_DIR, in_tsv_basename)

if (!file.exists(in_tsv)) {
  stop("[STEP_TR_B] Input TSV missing: ", in_tsv,
       "\n  Run STEP_TR_A_compute_theta_matrices.R for this chromosome first.")
}

# Resolve output JSON path
if (is.null(OUT_JSON)) {
  if (is.na(JSON_OUT_DIR)) {
    OUT_JSON <- file.path(THETA_TSV_DIR, sprintf("%s_phase2_theta.json", CHROM))
  } else {
    chrom_dir <- file.path(JSON_OUT_DIR, CHROM)
    dir.create(chrom_dir, recursive = TRUE, showWarnings = FALSE)
    OUT_JSON <- file.path(chrom_dir, sprintf("%s_phase2_theta.json", CHROM))
  }
}

message("[STEP_TR_B] CHROM=", CHROM,
        "  scale=", PESTPG_SCALE,
        "  mode=", THETA_GRID_MODE)
message("[STEP_TR_B] input  = ", in_tsv)
message("[STEP_TR_B] output = ", OUT_JSON)
message("[STEP_TR_B] params: z=", Z_THRESHOLD, " min_l2=", MIN_L2_WIN,
        " gap=", MERGE_GAP, " pad=", PAD, " max_k=", MAX_K)

# =============================================================================
# Load and pivot to wide samples × windows matrix
# =============================================================================

message("[STEP_TR_B] Loading long-format TSV...")
t0 <- proc.time()
long_dt <- fread(in_tsv)
# `theta_pi` here is per-site (tP / nSites); see STEP_TR_A header comment.
# `tP_sum` is the raw pestPG window sum, preserved by STEP_TR_A so the
# atlas can offer a normalized/raw toggle for diagnostic purposes.
required <- c("sample", "chrom", "window_idx", "start_bp", "end_bp",
              "theta_pi", "n_sites")
missing <- setdiff(required, names(long_dt))
if (length(missing) > 0) {
  stop("[STEP_TR_B] Input TSV missing columns: ", paste(missing, collapse = ", "))
}
# tP_sum is optional — STEP_TR_A introduced it after the per-site fix
# (ANGSD issue #329). If the TSV was produced by an older STEP_TR_A and
# only carries `theta_pi`, the raw-sum track will be NA in the JSON.
has_tP_sum <- "tP_sum" %in% names(long_dt)
if (!has_tP_sum) {
  message("[STEP_TR_B]   note: TSV lacks tP_sum column (pre-fix STEP_TR_A);",
          " atlas raw-sum toggle will be empty.")
}
long_dt <- long_dt[chrom == CHROM]
message("[STEP_TR_B]   rows = ", nrow(long_dt))

# Window grid (canonical, from data — assumes STEP_TR_A wrote a consistent grid)
win_grid <- unique(long_dt[, .(window_idx, start_bp, end_bp)])
setkey(win_grid, window_idx)
n_win <- nrow(win_grid)
message("[STEP_TR_B]   windows = ", n_win)

# Sample order
sample_order <- sort(unique(long_dt$sample))
n_samp <- length(sample_order)
message("[STEP_TR_B]   samples = ", n_samp)

# Pivot to wide samples × windows matrix
samp_to_row <- setNames(seq_along(sample_order), sample_order)
win_idx_vals <- win_grid$window_idx
win_to_col <- setNames(seq_along(win_idx_vals), as.character(win_idx_vals))

# theta_mat = per-site θπ (tP / nSites) — used for ALL analysis
#   (PCA, envelopes, |Z|). This is the diversity-comparable estimate.
# tP_sum_mat = raw pestPG window-sum θπ — diagnostic-only, exposed in
#   the JSON so the atlas can offer a normalized/raw toggle. Useful for
#   spotting regions where the cohort lost callable sites (raw-sum dips
#   that the per-site track flattens out by construction).
theta_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                    dimnames = list(sample_order, NULL))
n_sites_mat <- matrix(NA_integer_, nrow = n_samp, ncol = n_win,
                      dimnames = list(sample_order, NULL))
tP_sum_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                     dimnames = list(sample_order, NULL))

rows <- samp_to_row[long_dt$sample]
cols <- win_to_col[as.character(long_dt$window_idx)]
good <- !is.na(rows) & !is.na(cols)
fill_idx <- cbind(rows[good], cols[good])
theta_mat[fill_idx]   <- long_dt$theta_pi[good]
n_sites_mat[fill_idx] <- long_dt$n_sites[good]
if (has_tP_sum) {
  tP_sum_mat[fill_idx] <- long_dt$tP_sum[good]
}

message("[STEP_TR_B]   matrix coverage: ",
        round(100 * sum(!is.na(theta_mat)) / length(theta_mat), 1), "%",
        "  in ", round((proc.time() - t0)[3], 1), "s")

# =============================================================================
# Per-window population metrics
# =============================================================================
# For each window w:
#   median_theta_pi[w] = chromosome-level robust centre across 226 samples
#   mad_theta_pi[w]    = robust dispersion (median absolute deviation × 1.4826)
#   iqr_theta_pi[w]    = inter-quartile range (cross-check)
#   max_abs_z[w]       = maximum |sample_dev| / mad_dev across samples,
#                        where sample_dev = theta_mat[s, w] - median_theta_pi[w]
#                        and mad_dev = MAD across all (sample, window) deviations
# =============================================================================

message("[STEP_TR_B] Computing per-window population metrics...")
t1 <- proc.time()

window_median <- apply(theta_mat, 2, median, na.rm = TRUE)
window_mad    <- apply(theta_mat, 2, function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 10) return(NA_real_)
  mad(v, constant = 1.4826)
})
window_iqr <- apply(theta_mat, 2, function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 10) return(NA_real_)
  IQR(v)
})

# Per-(sample, window) deviation from the per-window median, expressed in
# units of the chromosome-wide median deviation (a robust z-score). This
# is what the page-12 |Z| panel will render: one trace per sample, with
# the running max-|z| collapsed to a 1D profile for envelope detection.
dev_mat <- sweep(theta_mat, 2, window_median, FUN = "-")

# Pool deviations across all samples and windows to get a robust scale
all_devs <- as.vector(dev_mat)
all_devs <- all_devs[is.finite(all_devs)]
chrom_dev_mad <- if (length(all_devs) > 100) mad(all_devs, constant = 1.4826) else NA_real_
if (!is.finite(chrom_dev_mad) || chrom_dev_mad <= 0) chrom_dev_mad <- 1e-6

# Per-window |Z| summary: maximum absolute z-score across samples at this
# window. This is the L2-envelope detection input.
max_abs_z <- apply(dev_mat, 2, function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 10) return(NA_real_)
  max(abs(v / chrom_dev_mad))
})

# Top-10% mean |Z| as an alternative summary that is less sensitive to a
# single outlier sample (useful when an L2 contains a coordinated
# karyotype subgroup rather than one anomalous individual).
top10_abs_z <- apply(dev_mat, 2, function(v) {
  v <- abs(v[is.finite(v)] / chrom_dev_mad)
  if (length(v) < 10) return(NA_real_)
  k <- max(1L, ceiling(0.1 * length(v)))
  mean(sort(v, decreasing = TRUE)[seq_len(k)])
})

# Per-window n_sites_floor — useful as a QC flag for downstream consumers
n_sites_floor <- apply(n_sites_mat, 2, function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NA_integer_)
  as.integer(min(v))
})

message("[STEP_TR_B]   metrics computed in ", round((proc.time() - t1)[3], 1), "s")
message("[STEP_TR_B]   |Z|: median=", round(median(max_abs_z, na.rm = TRUE), 2),
        " p95=", round(quantile(max_abs_z, 0.95, na.rm = TRUE), 2),
        " p99=", round(quantile(max_abs_z, 0.99, na.rm = TRUE), 2))

# =============================================================================
# Per-window local PCA on (2·PAD + 1) θπ neighbourhood
# =============================================================================
# For each focal window w, stack columns [w-PAD .. w+PAD] of theta_mat
# (n_samp × (2·PAD + 1) values), run a local PCA, and keep:
#   pc1_loadings[w][1..n_samp], pc2_loadings[w][1..n_samp]
#   lambda_ratio[w] = lambda_1 / lambda_2  (one-dimensionality indicator)
#   lambda_1[w], lambda_2[w]
#
# Heteroscedastic-noise mitigation: each sample's contribution is weighted
# by sqrt(n_sites[s, w_focal]) / sqrt(median n_sites over the chrom), so
# samples with sparse callable coverage in a given window contribute less
# to the local covariance. This is standard heteroscedastic-PCA practice.
# =============================================================================

message("[STEP_TR_B] Running per-window local PCA (pad=", PAD, ")...")
t2 <- proc.time()

pc1_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                  dimnames = list(sample_order, NULL))
pc2_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                  dimnames = list(sample_order, NULL))
lambda_ratio_vec <- rep(NA_real_, n_win)
lambda_1_vec     <- rep(NA_real_, n_win)
lambda_2_vec     <- rep(NA_real_, n_win)

n_sites_chrom_median <- median(n_sites_mat, na.rm = TRUE)
if (!is.finite(n_sites_chrom_median) || n_sites_chrom_median <= 0) {
  n_sites_chrom_median <- 1
}

for (wi in seq_len(n_win)) {
  lo <- max(1L, wi - PAD)
  hi <- min(n_win, wi + PAD)
  block <- theta_mat[, lo:hi, drop = FALSE]
  # Drop rows (samples) that have any NA in the block
  ok <- complete.cases(block)
  if (sum(ok) < max(20L, ncol(block) + 2L)) next

  block_ok <- block[ok, , drop = FALSE]

  # Heteroscedastic weight from focal-window n_sites
  n_focal <- n_sites_mat[ok, wi]
  n_focal[!is.finite(n_focal) | n_focal <= 0] <- 1L
  weights <- sqrt(n_focal / n_sites_chrom_median)
  block_w <- block_ok * weights

  # Centre, then SVD (samples × neighbourhood is small: 226 × 3 typically)
  centred <- sweep(block_w, 2, colMeans(block_w), FUN = "-")
  sv <- tryCatch(svd(centred, nu = 2, nv = 0), error = function(e) NULL)
  if (is.null(sv)) next

  # First two left singular vectors give the per-sample PC loadings
  pc1_mat[ok, wi] <- as.numeric(sv$u[, 1])
  if (ncol(sv$u) >= 2) pc2_mat[ok, wi] <- as.numeric(sv$u[, 2])

  # Eigenvalues of the sample-side covariance: d^2 / (n - 1)
  d2 <- sv$d^2
  if (length(d2) >= 2 && d2[2] > 0) {
    lambda_ratio_vec[wi] <- d2[1] / d2[2]
    lambda_1_vec[wi] <- d2[1]
    lambda_2_vec[wi] <- d2[2]
  } else if (length(d2) >= 1) {
    lambda_1_vec[wi] <- d2[1]
  }
}

message("[STEP_TR_B]   local PCA computed in ", round((proc.time() - t2)[3], 1), "s")
message("[STEP_TR_B]   λ₁/λ₂ ratio: median=",
        round(median(lambda_ratio_vec, na.rm = TRUE), 2),
        " p95=", round(quantile(lambda_ratio_vec, 0.95, na.rm = TRUE), 2))

# =============================================================================
# L2 envelope detection from contiguous high-|Z| runs
# =============================================================================

message("[STEP_TR_B] Calling L2 envelopes (z>", Z_THRESHOLD,
        ", min=", MIN_L2_WIN, " windows, gap=", MERGE_GAP, ")...")

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

  # Merge adjacent envelopes within merge_gap
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

l2_envs <- call_envelopes(max_abs_z, Z_THRESHOLD, MIN_L2_WIN, MERGE_GAP)
if (nrow(l2_envs) > 0) {
  l2_envs[, `:=`(
    start_bp  = win_grid$start_bp[win_start],
    end_bp    = win_grid$end_bp[win_end],
    span_kb   = round((win_grid$end_bp[win_end] -
                       win_grid$start_bp[win_start]) / 1000, 1),
    peak_z    = vapply(seq_len(nrow(l2_envs)),
                       function(k) max(max_abs_z[l2_envs$win_start[k]:l2_envs$win_end[k]],
                                       na.rm = TRUE),
                       numeric(1)),
    mean_z    = vapply(seq_len(nrow(l2_envs)),
                       function(k) mean(max_abs_z[l2_envs$win_start[k]:l2_envs$win_end[k]],
                                        na.rm = TRUE),
                       numeric(1))
  )]
  l2_envs[, l2_id := paste0(CHROM, "_thpi_L2_", sprintf("%03d", seq_len(.N)))]
}

# L1 = merge L2s within a wider gap (defaults to 3× merge_gap)
l1_merge_gap <- 3L * MERGE_GAP
l1_envs <- if (nrow(l2_envs) >= 1) {
  if (nrow(l2_envs) == 1) {
    copy(l2_envs)[, `:=`(l1_id = paste0(CHROM, "_thpi_L1_001"),
                          n_l2  = 1L)][]
  } else {
    setorder(l2_envs, win_start)
    merged <- l2_envs[1, .(win_start, win_end, n_l2 = 1L)]
    for (k in seq.int(2L, nrow(l2_envs))) {
      if (l2_envs$win_start[k] - merged$win_end[nrow(merged)] <= l1_merge_gap) {
        merged[nrow(merged), win_end := l2_envs$win_end[k]]
        merged[nrow(merged), n_l2 := n_l2 + 1L]
      } else {
        merged <- rbind(merged, l2_envs[k, .(win_start, win_end, n_l2 = 1L)])
      }
    }
    merged[, `:=`(
      start_bp  = win_grid$start_bp[win_start],
      end_bp    = win_grid$end_bp[win_end],
      span_kb   = round((win_grid$end_bp[win_end] -
                         win_grid$start_bp[win_start]) / 1000, 1),
      l1_id     = paste0(CHROM, "_thpi_L1_", sprintf("%03d", seq_len(.N)))
    )]
    merged[]
  }
} else {
  data.table()
}

message("[STEP_TR_B]   L2 envelopes: ", nrow(l2_envs))
message("[STEP_TR_B]   L1 envelopes: ", nrow(l1_envs))
if (nrow(l2_envs) > 0) {
  message("[STEP_TR_B]   L2 size range: ", round(min(l2_envs$span_kb), 0), "–",
          round(max(l2_envs$span_kb), 0), " kb")
}

# =============================================================================
# Optional per-interval karyotype assignment
# =============================================================================

interval_assignments <- list()
if (!is.null(INTERVAL_FILE) && file.exists(INTERVAL_FILE)) {
  message("[STEP_TR_B] Loading intervals: ", INTERVAL_FILE)
  intv_dt <- fread(INTERVAL_FILE)
  if ("chrom" %in% names(intv_dt)) intv_dt <- intv_dt[chrom == CHROM]
  message("[STEP_TR_B]   intervals on ", CHROM, ": ", nrow(intv_dt))

  for (ii in seq_len(nrow(intv_dt))) {
    intv     <- intv_dt[ii]
    intv_id  <- intv$interval_id %||% sprintf("I%03d", ii)
    win_in   <- which(win_grid$start_bp >= intv$start_bp &
                      win_grid$end_bp   <= intv$end_bp)
    if (length(win_in) < 5) next

    sub_mat <- theta_mat[, win_in, drop = FALSE]
    sample_means <- rowMeans(sub_mat, na.rm = TRUE)
    valid <- is.finite(sample_means)
    if (sum(valid) < 10) next

    v <- sample_means[valid]
    sample_ids <- sample_order[valid]

    best_k <- 2L; best_sil <- -1; best_cl <- NULL
    for (k in 2:min(MAX_K, length(unique(round(v, 6))) - 1L)) {
      if (k >= length(v)) next
      cl <- tryCatch(kmeans(v, centers = k, nstart = 25, iter.max = 100),
                     error = function(e) NULL)
      if (is.null(cl)) next
      # Approximate silhouette
      sil_vals <- numeric(length(v))
      for (j in seq_along(v)) {
        own <- cl$cluster[j]
        own_members <- v[cl$cluster == own & seq_along(v) != j]
        a_j <- if (length(own_members) > 0) mean(abs(v[j] - own_members)) else 0
        other_dists <- vapply(setdiff(seq_len(k), own),
                              function(other) {
                                m <- v[cl$cluster == other]
                                if (length(m) > 0) mean(abs(v[j] - m)) else NA_real_
                              }, numeric(1))
        b_j <- if (length(other_dists) > 0) min(other_dists, na.rm = TRUE) else 0
        sil_vals[j] <- if (max(a_j, b_j) > 0) (b_j - a_j) / max(a_j, b_j) else 0
      }
      mean_sil <- mean(sil_vals)
      if (mean_sil > best_sil) { best_sil <- mean_sil; best_k <- k; best_cl <- cl }
    }

    if (!is.null(best_cl)) {
      # Label clusters by ascending centre
      centre_order <- order(best_cl$centers[, 1])
      class_labels <- character(best_k)
      if (best_k == 2) {
        class_labels[centre_order[1]] <- "LOW_DIV"
        class_labels[centre_order[2]] <- "HIGH_DIV"
      } else if (best_k == 3) {
        class_labels[centre_order[1]] <- "LOW_DIV"
        class_labels[centre_order[2]] <- "MID_DIV"
        class_labels[centre_order[3]] <- "HIGH_DIV"
      } else {
        for (j in seq_len(best_k)) {
          class_labels[centre_order[j]] <- sprintf("DIV_TIER_%d", j)
        }
      }
      assigned <- class_labels[best_cl$cluster]

      interval_assignments[[length(interval_assignments) + 1L]] <- list(
        interval_id      = intv_id,
        chrom            = CHROM,
        start_bp         = as.integer(intv$start_bp),
        end_bp           = as.integer(intv$end_bp),
        n_interval_windows = length(win_in),
        k                = as.integer(best_k),
        silhouette       = round(best_sil, 4),
        sample_ids       = sample_ids,
        sample_class     = assigned,
        sample_mean_div  = round(v, 6)
      )
    }
  }
  message("[STEP_TR_B]   intervals classified: ", length(interval_assignments))
}

# =============================================================================
# Build JSON layers for the page-12 atlas
# =============================================================================

message("[STEP_TR_B] Assembling JSON layers...")

# Helper: round and replace NA / non-finite with NULL
clean_numeric <- function(x, digits = 6) {
  out <- round(as.numeric(x), digits)
  out[!is.finite(out)] <- NA_real_
  out
}

# Layer: theta_pi_per_window (the dense data layer)
# Carries TWO per-sample tracks so the atlas can toggle modes:
#   theta_pi  — per-site density (tP / nSites). Default rendering.
#   tP_sum    — raw pestPG window sum. Diagnostic mode.
# All analysis layers below (envelopes, local PCA, |Z|) use theta_pi.
samples_block <- vector("list", n_samp)
for (s in seq_len(n_samp)) {
  entry <- list(
    sample_id = sample_order[s],
    theta_pi  = clean_numeric(theta_mat[s, ], 6),
    n_sites   = as.integer(n_sites_mat[s, ])
  )
  if (has_tP_sum) {
    # 4 digits is enough for diagnostic display; tP_sum has scale ~5
    entry$tP_sum <- clean_numeric(tP_sum_mat[s, ], 4)
  }
  samples_block[[s]] <- entry
}

# v4 turn 114: also build a flat row-major `values` vector (length n_w * n_s,
# row = window, col = sample). This matches the atlas's canonical Shape A
# contract for theta_pi_per_window.values (line ~32521 of Inversion_atlas.html
# detectSchemaAndLayers). It's redundant with samples[i].theta_pi (every value
# appears in both forms) but lets the atlas's vectorized readers grab the full
# matrix in one shot rather than walking the samples array. Cost: ~25% file
# size increase for theta_pi_per_window block — acceptable.
# (The atlas's case-block accepts EITHER form per the turn-114 dual-shape
# patch, so emitting both is purely a cleanness measure for new runs.)
values_flat <- numeric(n_win * n_samp)
for (wi in seq_len(n_win)) {
  values_flat[((wi - 1L) * n_samp + 1L):(wi * n_samp)] <-
    clean_numeric(theta_mat[, wi], 6)
}

theta_pi_per_window <- list(
  schema_version = 1L,
  layer          = "theta_pi_per_window",
  chrom          = CHROM,
  scale          = if (THETA_GRID_MODE == "native") PESTPG_SCALE else "dosage_grid",
  grid_mode      = THETA_GRID_MODE,
  # Available render modes for the page-12 toggle:
  #   "per_site"  → reads samples[].theta_pi (default; per-site θπ)
  #   "raw_sum"   → reads samples[].tP_sum   (diagnostic; pestPG window sum)
  # If has_tP_sum is FALSE, only "per_site" is available.
  available_modes = if (has_tP_sum) c("per_site", "raw_sum") else c("per_site"),
  default_mode    = "per_site",
  n_samples      = as.integer(n_samp),
  n_windows      = as.integer(n_win),
  # v4 turn 114: sample_ids = top-level string array of sample IDs in the
  # canonical row order used by `values`. Distinct from `samples` (which
  # is the array of per-sample-object entries, sample-major). Both are
  # emitted so atlas readers can pick whichever form they need.
  sample_ids     = sample_order,
  values         = clean_numeric(values_flat, 6),
  windows        = lapply(seq_len(n_win), function(wi) list(
    idx        = as.integer(win_grid$window_idx[wi]),
    start_bp   = as.integer(win_grid$start_bp[wi]),
    end_bp     = as.integer(win_grid$end_bp[wi]),
    n_sites_floor = as.integer(n_sites_floor[wi])
  )),
  samples        = samples_block
)

# Layer: theta_pi_local_pca
# v4 turn 114: emit BOTH `z_profile` (TR_B's original name, kept for
# backward compat with consumers / atlas-side caches built before turn 114)
# AND `z` (atlas-canonical contract name). They're identical numeric
# vectors — the atlas's case-block / recovery sweep accepts either, but
# emitting both means an atlas tab that reads `z` directly (without going
# through detectSchemaAndLayers's tolerance code path) just works.
theta_pi_local_pca <- list(
  schema_version = 1L,
  layer          = "theta_pi_local_pca",
  chrom          = CHROM,
  scale          = if (THETA_GRID_MODE == "native") PESTPG_SCALE else "dosage_grid",
  pad            = as.integer(PAD),
  n_samples      = as.integer(n_samp),
  n_windows      = as.integer(n_win),
  sample_order   = sample_order,
  pc1_loadings   = lapply(seq_len(n_win), function(wi) clean_numeric(pc1_mat[, wi], 6)),
  pc2_loadings   = lapply(seq_len(n_win), function(wi) clean_numeric(pc2_mat[, wi], 6)),
  lambda_1       = clean_numeric(lambda_1_vec, 6),
  lambda_2       = clean_numeric(lambda_2_vec, 6),
  lambda_ratio   = clean_numeric(lambda_ratio_vec, 4),
  z              = clean_numeric(max_abs_z, 4),     # atlas-canonical name
  z_profile      = clean_numeric(max_abs_z, 4),     # TR_B legacy name (alias)
  z_top10_mean   = clean_numeric(top10_abs_z, 4)
)

# Layer: theta_pi_envelopes
theta_pi_envelopes <- list(
  schema_version = 1L,
  layer          = "theta_pi_envelopes",
  chrom          = CHROM,
  z_threshold    = Z_THRESHOLD,
  min_l2_windows = as.integer(MIN_L2_WIN),
  merge_gap      = as.integer(MERGE_GAP),
  l2 = if (nrow(l2_envs) > 0) lapply(seq_len(nrow(l2_envs)), function(k) list(
    l2_id     = l2_envs$l2_id[k],
    win_start = as.integer(l2_envs$win_start[k]),
    win_end   = as.integer(l2_envs$win_end[k]),
    start_bp  = as.integer(l2_envs$start_bp[k]),
    end_bp    = as.integer(l2_envs$end_bp[k]),
    span_kb   = round(l2_envs$span_kb[k], 1),
    n_windows = as.integer(l2_envs$n_windows[k]),
    peak_z    = round(l2_envs$peak_z[k], 4),
    mean_z    = round(l2_envs$mean_z[k], 4)
  )) else list(),
  l1 = if (nrow(l1_envs) > 0) lapply(seq_len(nrow(l1_envs)), function(k) list(
    l1_id     = l1_envs$l1_id[k],
    win_start = as.integer(l1_envs$win_start[k]),
    win_end   = as.integer(l1_envs$win_end[k]),
    start_bp  = as.integer(l1_envs$start_bp[k]),
    end_bp    = as.integer(l1_envs$end_bp[k]),
    span_kb   = round(l1_envs$span_kb[k], 1),
    n_l2      = as.integer(l1_envs$n_l2[k])
  )) else list()
)

# Layer: tracks (per-window aggregates aligned with theta_pi_per_window's
# window grid, additive into the atlas's existing tracks dict)
tracks_layer <- list(
  theta_pi_median = list(
    min    = round(min(window_median, na.rm = TRUE), 6),
    max    = round(max(window_median, na.rm = TRUE), 6),
    values = clean_numeric(window_median, 6),
    pos_bp = as.integer(round((win_grid$start_bp + win_grid$end_bp) / 2))
  ),
  theta_pi_z = list(
    min    = round(min(max_abs_z, na.rm = TRUE), 4),
    max    = round(max(max_abs_z, na.rm = TRUE), 4),
    values = clean_numeric(max_abs_z, 4),
    pos_bp = as.integer(round((win_grid$start_bp + win_grid$end_bp) / 2))
  ),
  theta_pi_lambda_ratio = list(
    min    = round(min(lambda_ratio_vec, na.rm = TRUE), 4),
    max    = round(max(lambda_ratio_vec, na.rm = TRUE), 4),
    values = clean_numeric(lambda_ratio_vec, 4),
    pos_bp = as.integer(round((win_grid$start_bp + win_grid$end_bp) / 2))
  )
)

# Top-level JSON
# v4 turn 114: schema_version bumped 1 → 2. Now that the per-layer blocks
# emit atlas-canonical field names (theta_pi_per_window.values + .sample_ids,
# theta_pi_local_pca.z) alongside the legacy names, the JSON conforms to
# the v2 contract and the atlas's case-block detection path activates
# directly without falling through to inferLayersFromV1's recovery sweep.
# The atlas accepts either schema_version=1 OR =2 with these layer shapes
# (turn-114 dual-tolerance patch in detectSchemaAndLayers + inferLayersFromV1),
# so this bump is a cleanness measure, not a functional gate.
out_json_obj <- list(
  schema_version    = 2L,
  chrom             = CHROM,
  n_samples         = as.integer(n_samp),
  n_windows         = as.integer(n_win),
  scale             = if (THETA_GRID_MODE == "native") PESTPG_SCALE else "dosage_grid",
  grid_mode         = THETA_GRID_MODE,
  `_generated_at`   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  `_generator`      = "STEP_TR_B_classify_theta.R",
  `_layers_present` = c("theta_pi_per_window", "theta_pi_local_pca",
                        "theta_pi_envelopes", "tracks"),
  tracks                = tracks_layer,
  theta_pi_per_window   = theta_pi_per_window,
  theta_pi_local_pca    = theta_pi_local_pca,
  theta_pi_envelopes    = theta_pi_envelopes
)

if (length(interval_assignments) > 0) {
  out_json_obj$theta_pi_intervals <- interval_assignments
  out_json_obj$`_layers_present`  <- c(out_json_obj$`_layers_present`,
                                       "theta_pi_intervals")
}

# Write JSON (compact form; the atlas reader handles either pretty or compact)
write_json(out_json_obj, OUT_JSON,
           auto_unbox = TRUE, na = "null", pretty = FALSE, digits = NA)
fi <- file.info(OUT_JSON)
message("[STEP_TR_B] Wrote ", OUT_JSON,
        " (", round(fi$size / 1024 / 1024, 2), " MB)")

message("[STEP_TR_B] DONE — chrom=", CHROM)
