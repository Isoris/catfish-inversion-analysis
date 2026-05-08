#!/usr/bin/env Rscript

# =============================================================================
# STEP_C04b_snake3_ghsl_classify.R
#
# SNAKE 3 v6: LIGHT CLASSIFIER — Scoring + Karyotype + Decomposition
#
# PURPOSE:
#   Loads the pre-computed rolling divergence matrices from STEP_C04 (v6).
#   Runs in ~30 seconds per chromosome. Iterate 100 times while tuning.
#
# WHAT IT DOES:
#   A. METRICS (on rolling-smoothed data):
#      - Rank stability: Spearman rho between adjacent rolling windows
#      - Per-window distribution: mean, bimodality, contrast, tightness
#      - Per-window status: PASS/WEAK/FAIL (z-score against chr baseline)
#
#   B. KARYOTYPE CALLING (on rolling-smoothed ranks):
#      - Per-sample: runs of stable LOW (INV/INV) or HIGH (INV/nonINV)
#      - Much cleaner than raw-window calling because rolling kills noise
#
#   C. INTERVAL CLASSIFICATION (requires triangle/stair intervals):
#      - Per sample × per interval: extract rolling profile vector
#      - Whole-interval mean → rank samples → k-means (k=2..5) + silhouette
#      - Classify: INV/INV, HET, nonINV/nonINV per sample per interval
#
#   D. INTERVAL DECOMPOSITION (detects sub-systems):
#      - Per sample: changepoint detection on rolling profile within interval
#      - Group samples by changepoint location → sub-system membership
#      - Detects when an interval contains 2+ overlapping inversion systems
#
# OUTPUT:
#   snake3v6_window_track.tsv.gz       — per-window metrics (rolling-based)
#   snake3v6_karyotype_calls.tsv.gz    — per-sample stable runs
#   snake3v6_interval_genotypes.tsv.gz — per-sample × per-interval classification
#   snake3v6_interval_decomp.tsv.gz    — sub-system decomposition (if intervals given)
#   snake3v6_summary.tsv               — chromosome-level summary
#
# Usage:
#   Rscript STEP_C04b_snake3_ghsl_classify.R <matrices_dir> <outdir> \
#     [--chrom C_gar_LG01] \
#     [--scale 50] \
#     [--intervals <triangle_intervals.tsv.gz>] \
#     [--karyo_lo 0.15] [--karyo_hi 0.70] [--karyo_min_run 10] \
#     [--rank_window 5] \
#     [--max_k 5]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop(paste(
  "Usage: Rscript STEP_C04b_snake3_ghsl_classify.R <matrices_dir> <outdir> [opts]",
  "  <matrices_dir>  Dir with *.ghsl_v6_matrices.rds from STEP_C04",
  "  <outdir>        Output directory",
  sep = "\n"
))

matrices_dir <- args[1]
outdir       <- args[2]

CHROM_FILTER       <- NULL
SCALE              <- 50L        # which rolling scale to use (must match a saved scale)
INTERVAL_FILE      <- NULL       # triangle_intervals.tsv.gz
KARYO_LO           <- 0.15       # bottom quantile → INV/INV
KARYO_HI           <- 0.70       # top quantile → INV/nonINV
KARYO_MIN_RUN      <- 10L        # min consecutive windows for stable call
RANK_WINDOW        <- 5L         # smoothing for rank stability
MAX_K              <- 5L         # max clusters for interval classification

i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom" && i < length(args))        { CHROM_FILTER <- args[i + 1]; i <- i + 2L }
  else if (a == "--scale" && i < length(args))   { SCALE <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--intervals" && i < length(args)){ INTERVAL_FILE <- args[i + 1]; i <- i + 2L }
  else if (a == "--karyo_lo" && i < length(args)) { KARYO_LO <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--karyo_hi" && i < length(args)) { KARYO_HI <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--karyo_min_run" && i < length(args)) { KARYO_MIN_RUN <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--rank_window" && i < length(args))   { RANK_WINDOW <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--max_k" && i < length(args))         { MAX_K <- as.integer(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("================================================================")
message("[S3v6b] Snake 3 v6b: LIGHT CLASSIFIER")
message("================================================================")
message("[S3v6b] Matrices:   ", matrices_dir)
message("[S3v6b] Output:     ", outdir)
message("[S3v6b] Scale:      s", SCALE, " (~", SCALE * 5, " kb rolling)")
message("[S3v6b] Karyotype:  lo=", KARYO_LO, " hi=", KARYO_HI,
        " min_run=", KARYO_MIN_RUN)
if (!is.null(INTERVAL_FILE)) message("[S3v6b] Intervals:  ", INTERVAL_FILE)

# Load intervals if provided
intervals_dt <- NULL
if (!is.null(INTERVAL_FILE) && file.exists(INTERVAL_FILE)) {
  intervals_dt <- fread(INTERVAL_FILE)
  message("[S3v6b] Loaded ", nrow(intervals_dt), " intervals")
}

# =============================================================================
# A. ROLLING METRICS: rank stability, bimodality, contrast, scoring
# =============================================================================

compute_rolling_metrics <- function(roll_mat, rank_window = 5L,
                                    karyo_lo = 0.15, karyo_hi = 0.70) {
  n_samp <- nrow(roll_mat)
  n_win  <- ncol(roll_mat)

  # ── A1: Rank matrix (normalized 0-1 per window) ──
  rank_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                     dimnames = dimnames(roll_mat))
  for (wi in seq_len(n_win)) {
    vals <- roll_mat[, wi]
    valid <- !is.na(vals)
    if (sum(valid) < 10) next
    rank_mat[valid, wi] <- rank(vals[valid], ties.method = "average") / sum(valid)
  }

  # ── A2: Rank stability (Spearman rho between adjacent windows) ──
  pairwise_rho <- rep(NA_real_, n_win)
  for (wi in seq_len(n_win - 1)) {
    r1 <- rank_mat[, wi]
    r2 <- rank_mat[, wi + 1]
    valid <- !is.na(r1) & !is.na(r2)
    if (sum(valid) < 20) next
    pairwise_rho[wi] <- cor(r1[valid], r2[valid], method = "spearman")
  }

  # Smoothed rho
  smooth_rho <- rep(NA_real_, n_win)
  half <- floor(rank_window / 2)
  for (wi in seq_len(n_win)) {
    lo <- max(1, wi - half)
    hi <- min(n_win - 1, wi + half)
    rhos <- pairwise_rho[lo:hi]
    rhos <- rhos[!is.na(rhos)]
    if (length(rhos) >= 2) smooth_rho[wi] <- mean(rhos)
  }

  # ── A3: Per-window distribution metrics ──
  metrics <- data.table(
    window_idx     = seq_len(n_win),
    n_scored       = integer(n_win),
    div_mean       = numeric(n_win),
    div_median     = numeric(n_win),
    div_sd         = numeric(n_win),
    div_iqr        = numeric(n_win),
    div_skew       = numeric(n_win),
    div_bimodal    = numeric(n_win),
    n_low_div      = integer(n_win),
    n_high_div     = integer(n_win),
    n_mid_div      = integer(n_win),
    low_div_tightness  = numeric(n_win),
    high_div_tightness = numeric(n_win),
    mid_div_spread     = numeric(n_win)
  )

  for (wi in seq_len(n_win)) {
    vals <- roll_mat[, wi]
    valid <- !is.na(vals)
    n_valid <- sum(valid)
    metrics$n_scored[wi] <- n_valid
    if (n_valid < 20) next

    v <- vals[valid]
    metrics$div_mean[wi]   <- mean(v)
    metrics$div_median[wi] <- median(v)
    metrics$div_sd[wi]     <- sd(v)
    metrics$div_iqr[wi]    <- IQR(v)

    # Skewness
    m3 <- mean((v - mean(v))^3)
    s3 <- sd(v)^3
    metrics$div_skew[wi] <- if (s3 > 0) m3 / s3 else 0

    # Bimodality: count density peaks
    if (n_valid >= 30 && sd(v) > 0.001) {
      metrics$div_bimodal[wi] <- tryCatch({
        dens <- density(v, n = 64)
        dv <- diff(sign(diff(dens$y)))
        sum(dv == -2) + 1
      }, error = function(e) 1)
    } else {
      metrics$div_bimodal[wi] <- 1
    }

    # Quantile groups
    q_lo <- quantile(v, karyo_lo)
    q_hi <- quantile(v, karyo_hi)
    low_mask  <- v <= q_lo
    high_mask <- v >= q_hi
    mid_mask  <- v > q_lo & v < q_hi

    metrics$n_low_div[wi]  <- sum(low_mask)
    metrics$n_high_div[wi] <- sum(high_mask)
    metrics$n_mid_div[wi]  <- sum(mid_mask)

    metrics$low_div_tightness[wi]  <- if (sum(low_mask) >= 5) IQR(v[low_mask]) else NA_real_
    metrics$high_div_tightness[wi] <- if (sum(high_mask) >= 5) IQR(v[high_mask]) else NA_real_
    metrics$mid_div_spread[wi]     <- if (sum(mid_mask) >= 5) IQR(v[mid_mask]) else NA_real_
  }

  # ── A4: Z-score scoring ──
  metrics[, rank_stability := smooth_rho[seq_len(n_win)]]

  # Divergence contrast
  for (wi in seq_len(n_win)) {
    vals <- roll_mat[, wi]
    valid <- !is.na(vals)
    if (sum(valid) < 20) next
    v <- vals[valid]
    q_lo_v <- quantile(v, karyo_lo)
    q_hi_v <- quantile(v, karyo_hi)
    lo_vals <- v[v <= q_lo_v]
    hi_vals <- v[v >= q_hi_v]
    lo_m <- if (length(lo_vals) > 0) mean(lo_vals) else NA_real_
    hi_m <- if (length(hi_vals) > 0) mean(hi_vals) else NA_real_
    set(metrics, wi, "div_contrast", hi_m - lo_m)
  }
  if (!"div_contrast" %in% names(metrics)) metrics[, div_contrast := NA_real_]

  # Z-score against chromosome baseline
  rho_valid <- metrics$rank_stability[!is.na(metrics$rank_stability)]
  rho_med <- if (length(rho_valid) > 10) median(rho_valid) else 0.5
  rho_mad <- max(if (length(rho_valid) > 10) mad(rho_valid, constant = 1.4826) else 0.1, 0.02)

  dc_valid <- metrics$div_contrast[!is.na(metrics$div_contrast) & metrics$n_scored >= 20]
  dc_med <- if (length(dc_valid) > 10) median(dc_valid) else 0
  dc_mad <- max(if (length(dc_valid) > 10) mad(dc_valid, constant = 1.4826) else 0.1, 0.02)

  metrics[, rank_stability_z := fifelse(!is.na(rank_stability),
    (rank_stability - rho_med) / rho_mad, NA_real_)]
  metrics[, div_contrast_z := fifelse(!is.na(div_contrast),
    (div_contrast - dc_med) / dc_mad, NA_real_)]

  message("[S3v6b]   Baseline: rho_med=", round(rho_med, 4), " rho_mad=", round(rho_mad, 4),
          " contrast_med=", round(dc_med, 4), " contrast_mad=", round(dc_mad, 4))

  # Combined score
  metrics[, ghsl_v6_score := fifelse(
    !is.na(rank_stability_z) & !is.na(div_contrast_z) & n_scored >= 20,
    pnorm(rank_stability_z) * 0.40 +
    pnorm(div_contrast_z)   * 0.25 +
    fifelse(div_bimodal >= 2, 0.20, 0) +
    fifelse(!is.na(low_div_tightness) & low_div_tightness < 0.05, 0.15, 0),
    NA_real_
  )]

  # Status
  metrics[, ghsl_v6_status := fifelse(
    is.na(ghsl_v6_score), "FAIL",
    fifelse(ghsl_v6_score > 0.65 & rank_stability_z > 1.5, "PASS",
    fifelse(ghsl_v6_score > 0.50 & rank_stability_z > 0.5, "WEAK", "FAIL"))
  )]

  list(metrics = metrics, rank_mat = rank_mat, smooth_rho = smooth_rho)
}

# =============================================================================
# B. KARYOTYPE CALLING: per-sample runs on rolling ranks
# =============================================================================

call_karyotypes <- function(rank_mat, min_run = 10L,
                            karyo_lo = 0.15, karyo_hi = 0.70) {
  n_samples <- nrow(rank_mat)
  n_win     <- ncol(rank_mat)
  snames    <- rownames(rank_mat)
  karyo_calls <- list()

  for (si in seq_len(n_samples)) {
    ranks <- rank_mat[si, ]

    state <- rep("NA", n_win)
    state[!is.na(ranks) & ranks <= karyo_lo]                           <- "LOW"
    state[!is.na(ranks) & ranks >= karyo_hi]                           <- "HIGH"
    state[!is.na(ranks) & ranks > karyo_lo & ranks < karyo_hi]        <- "MID"

    rle_state   <- rle(state)
    run_ends    <- cumsum(rle_state$lengths)
    run_starts  <- c(1L, head(run_ends, -1) + 1L)

    for (ri in seq_along(rle_state$values)) {
      if (rle_state$values[ri] %in% c("LOW", "HIGH") &&
          rle_state$lengths[ri] >= min_run) {
        karyo_calls[[length(karyo_calls) + 1L]] <- data.table(
          sample_id    = snames[si],
          window_start = run_starts[ri],
          window_end   = run_ends[ri],
          n_windows    = rle_state$lengths[ri],
          call         = if (rle_state$values[ri] == "LOW") "INV_INV" else "INV_nonINV",
          mean_rank    = mean(ranks[run_starts[ri]:run_ends[ri]], na.rm = TRUE)
        )
      }
    }
  }

  if (length(karyo_calls) > 0) rbindlist(karyo_calls) else data.table()
}

# =============================================================================
# C. INTERVAL CLASSIFICATION: per-sample karyotype within each interval
# =============================================================================

classify_interval <- function(roll_mat, interval_windows, max_k = 5L) {
  # Extract sub-matrix for this interval's windows
  # roll_mat rows = samples, cols = windows
  # interval_windows = integer vector of window indices within the interval
  
  sub_mat <- roll_mat[, interval_windows, drop = FALSE]
  
  # Per-sample mean divergence across the interval
  sample_means <- rowMeans(sub_mat, na.rm = TRUE)
  valid <- !is.na(sample_means)
  
  if (sum(valid) < 10) {
    return(data.table(
      sample_id = rownames(roll_mat),
      interval_mean_div = sample_means,
      interval_class = "UNCLASSIFIED",
      interval_k = NA_integer_,
      silhouette = NA_real_
    ))
  }
  
  v <- sample_means[valid]
  sample_ids <- rownames(roll_mat)[valid]
  
  # Try k = 2..max_k, pick best silhouette
  best_k   <- 2L
  best_sil <- -1
  best_cl  <- NULL
  
  for (k in 2:min(max_k, length(unique(round(v, 6))) - 1)) {
    if (k >= length(v)) next
    cl <- tryCatch(kmeans(v, centers = k, nstart = 25, iter.max = 100),
                   error = function(e) NULL)
    if (is.null(cl)) next
    
    # Silhouette (simplified: mean within-cluster distance vs nearest-cluster)
    sil_vals <- numeric(length(v))
    for (i in seq_along(v)) {
      own_cluster <- cl$cluster[i]
      own_members <- v[cl$cluster == own_cluster & seq_along(v) != i]
      a_i <- if (length(own_members) > 0) mean(abs(v[i] - own_members)) else 0
      
      other_dists <- numeric(0)
      for (kk in seq_len(k)) {
        if (kk == own_cluster) next
        other_members <- v[cl$cluster == kk]
        if (length(other_members) > 0) {
          other_dists <- c(other_dists, mean(abs(v[i] - other_members)))
        }
      }
      b_i <- if (length(other_dists) > 0) min(other_dists) else 0
      sil_vals[i] <- if (max(a_i, b_i) > 0) (b_i - a_i) / max(a_i, b_i) else 0
    }
    
    mean_sil <- mean(sil_vals)
    if (mean_sil > best_sil) {
      best_sil <- mean_sil
      best_k   <- k
      best_cl  <- cl
    }
  }
  
  # Assign class labels based on cluster center ordering
  # Lowest center = INV_INV (both haplotypes locked → low divergence)
  # Highest center = INV_nonINV (haplotypes can't recombine → high divergence)
  # Middle = HET or SOUP
  #
  # Chat 14 (2026-04-18): k=4/k=5 labels changed from generic CLASS_N to
  # divergence-rank scheme. At k=4 the extra subgroup is an intermediate
  # divergence class; at k=5 HET sits in the middle as the balance
  # point. Part D (CUSUM changepoint clustering) is the sub-inversion-
  # boundary axis and reports its own n_clusters independently — two
  # axes, two outputs, not conflated.
  if (!is.null(best_cl)) {
    center_order <- order(best_cl$centers[, 1])
    class_labels <- character(best_k)
    if (best_k == 2) {
      class_labels[center_order[1]] <- "LOW_DIV"
      class_labels[center_order[2]] <- "HIGH_DIV"
    } else if (best_k == 3) {
      class_labels[center_order[1]] <- "INV_INV"
      class_labels[center_order[2]] <- "HET"
      class_labels[center_order[3]] <- "INV_nonINV"
    } else if (best_k == 4) {
      # Four divergence bins: INV_INV / INTER_LOW / INTER_HIGH / INV_nonINV
      class_labels[center_order[1]] <- "INV_INV"
      class_labels[center_order[2]] <- "INTER_LOW"
      class_labels[center_order[3]] <- "INTER_HIGH"
      class_labels[center_order[4]] <- "INV_nonINV"
    } else if (best_k == 5) {
      # Five bins with HET at the middle of the divergence gradient
      class_labels[center_order[1]] <- "INV_INV"
      class_labels[center_order[2]] <- "INTER_LOW"
      class_labels[center_order[3]] <- "HET"
      class_labels[center_order[4]] <- "INTER_HIGH"
      class_labels[center_order[5]] <- "INV_nonINV"
    } else {
      # k > 5 (shouldn't happen with MAX_K=5 default): rank-labelled
      for (j in seq_len(best_k)) {
        class_labels[center_order[j]] <- paste0("CLASS_", j)
      }
    }

    assigned_classes <- class_labels[best_cl$cluster]
  } else {
    assigned_classes <- rep("UNCLASSIFIED", sum(valid))
  }
  
  # Build output for ALL samples (including those with NA)
  out <- data.table(
    sample_id         = rownames(roll_mat),
    interval_mean_div = sample_means,
    interval_class    = "NA_DATA",
    interval_k        = best_k,
    silhouette        = best_sil
  )
  out[valid, interval_class := assigned_classes]
  
  out
}

# =============================================================================
# D. INTERVAL DECOMPOSITION: detect sub-systems via changepoint
# =============================================================================

decompose_interval <- function(roll_mat, interval_windows, sample_classes) {
  # For each sample, look at the divergence PROFILE across windows
  # within the interval. Find changepoints → where does divergence
  # shift significantly?
  #
  # Samples whose changepoints fall at the same window position
  # share the same sub-inversion boundary.
  #
  # Uses simple approach: sliding variance + threshold.
  # A changepoint is where the local mean shifts significantly.
  
  sub_mat <- roll_mat[, interval_windows, drop = FALSE]
  n_samp  <- nrow(sub_mat)
  n_win   <- ncol(sub_mat)
  snames  <- rownames(sub_mat)
  
  if (n_win < 10) return(data.table())
  
  # Per-sample changepoint detection (simple: biggest absolute diff in
  # cumulative sum deviation from mean)
  decomp_list <- list()
  
  for (si in seq_len(n_samp)) {
    profile <- sub_mat[si, ]
    valid   <- !is.na(profile)
    if (sum(valid) < 10) next
    
    # Fill NAs with local mean for changepoint detection
    profile_filled <- profile
    profile_filled[!valid] <- mean(profile, na.rm = TRUE)
    
    # CUSUM-based changepoint: cumulative sum of deviations from mean
    mu <- mean(profile_filled)
    cusum <- cumsum(profile_filled - mu)
    
    # The changepoint is where |cusum| is maximized
    # Multiple changepoints: find top N peaks in |cusum|
    abs_cusum <- abs(cusum)
    
    # Simple: find the single strongest changepoint
    cp_idx <- which.max(abs_cusum)
    cp_strength <- abs_cusum[cp_idx] / (n_win * sd(profile_filled))
    
    # Also compute: left-half mean vs right-half mean
    left_mean  <- mean(profile_filled[1:cp_idx], na.rm = TRUE)
    right_mean <- mean(profile_filled[(cp_idx + 1):n_win], na.rm = TRUE)
    asymmetry  <- right_mean - left_mean
    
    # Profile flatness: low variance = flat = consistent karyotype
    profile_var <- var(profile_filled)
    
    # Slope: linear trend across interval
    slope <- tryCatch(
      coef(lm(profile_filled ~ seq_len(n_win)))[2],
      error = function(e) NA_real_
    )
    
    decomp_list[[length(decomp_list) + 1L]] <- data.table(
      sample_id       = snames[si],
      changepoint_win = cp_idx,
      changepoint_pos = interval_windows[cp_idx],  # global window index
      cp_strength     = round(cp_strength, 4),
      left_mean       = round(left_mean, 4),
      right_mean      = round(right_mean, 4),
      asymmetry       = round(asymmetry, 4),
      profile_var     = round(profile_var, 6),
      profile_slope   = round(as.numeric(slope), 6)
    )
  }
  
  if (length(decomp_list) == 0) return(data.table())
  decomp_dt <- rbindlist(decomp_list)
  
  # Cluster samples by changepoint location
  # Samples with changepoints at similar positions share a sub-system boundary
  if (nrow(decomp_dt) >= 5) {
    cp_vals <- decomp_dt$changepoint_win
    
    # Only cluster if there's variation in changepoint location
    if (sd(cp_vals) > 2) {
      cp_cl <- tryCatch({
        # Try k=2 first (most common: one interval = two sub-systems)
        km <- kmeans(cp_vals, centers = 2, nstart = 25)
        # Check if the two centers are meaningfully separated
        center_gap <- abs(diff(sort(km$centers[, 1])))
        if (center_gap > n_win * 0.1) {  # at least 10% of interval apart
          km$cluster
        } else {
          rep(1L, length(cp_vals))  # not separated enough → one system
        }
      }, error = function(e) rep(1L, length(cp_vals)))
      
      decomp_dt[, cp_cluster := cp_cl]
    } else {
      decomp_dt[, cp_cluster := 1L]
    }
  } else {
    decomp_dt[, cp_cluster := 1L]
  }
  
  decomp_dt
}

# =============================================================================
# MAIN LOOP
# =============================================================================

rds_files <- sort(list.files(matrices_dir, pattern = "\\.ghsl_v6_matrices\\.rds$",
                             full.names = TRUE))
if (length(rds_files) == 0) stop("[S3v6b] FATAL: No .ghsl_v6_matrices.rds files in: ",
                                  matrices_dir)
message("[S3v6b] Found ", length(rds_files), " matrix files")

if (!is.null(CHROM_FILTER)) {
  rds_files <- rds_files[grepl(CHROM_FILTER, rds_files)]
}

all_metrics    <- list()
all_karyo      <- list()
all_interval   <- list()
all_decomp     <- list()
all_summary    <- list()

for (rds_f in rds_files) {
  t0 <- proc.time()
  message("\n================================================================")
  message("[S3v6b] Loading: ", basename(rds_f))

  mat_data <- readRDS(rds_f)
  chr         <- mat_data$chrom
  div_mat     <- mat_data$div_mat
  het_mat     <- mat_data$het_mat
  window_info <- mat_data$window_info
  snames      <- mat_data$sample_names
  n_win       <- ncol(div_mat)

  message("[S3v6b] ", chr, ": ", length(snames), " samples × ", n_win, " windows")

  # Select rolling scale
  scale_key <- paste0("s", SCALE)
  if (!scale_key %in% names(mat_data$rolling)) {
    avail <- paste(names(mat_data$rolling), collapse = ", ")
    message("[S3v6b] WARNING: scale '", scale_key, "' not found. Available: ", avail)
    message("[S3v6b]   Falling back to first available scale")
    scale_key <- names(mat_data$rolling)[1]
    message("[S3v6b]   Using: ", scale_key)
  }

  roll_mat <- mat_data$rolling[[scale_key]]
  message("[S3v6b] Using rolling scale: ", scale_key)

  # ── A: Rolling metrics ──
  message("[S3v6b] Computing rolling metrics...")
  result_A <- compute_rolling_metrics(roll_mat, rank_window = RANK_WINDOW,
                                       karyo_lo = KARYO_LO, karyo_hi = KARYO_HI)
  metrics  <- result_A$metrics
  rank_mat <- result_A$rank_mat

  n_pass <- sum(metrics$ghsl_v6_status == "PASS", na.rm = TRUE)
  n_weak <- sum(metrics$ghsl_v6_status == "WEAK", na.rm = TRUE)
  n_fail <- sum(metrics$ghsl_v6_status == "FAIL", na.rm = TRUE)
  message("[S3v6b]   PASS=", n_pass, " WEAK=", n_weak, " FAIL=", n_fail)

  # Add coords
  metrics[, `:=`(
    chrom            = chr,
    global_window_id = window_info$global_window_id[seq_len(n_win)],
    start_bp         = window_info$start_bp[seq_len(n_win)],
    end_bp           = window_info$end_bp[seq_len(n_win)],
    pos_mb           = window_info$pos_mb[seq_len(n_win)],
    rolling_scale    = scale_key
  )]

  # ── B: Karyotype calling ──
  message("[S3v6b] Calling karyotypes on rolling ranks...")
  karyo_dt <- call_karyotypes(rank_mat, min_run = KARYO_MIN_RUN,
                               karyo_lo = KARYO_LO, karyo_hi = KARYO_HI)
  if (nrow(karyo_dt) > 0) {
    karyo_dt[, chrom := chr]
    karyo_dt[, `:=`(
      start_bp = window_info$start_bp[window_start],
      end_bp   = window_info$end_bp[window_end],
      start_mb = round(window_info$start_bp[window_start] / 1e6, 2),
      end_mb   = round(window_info$end_bp[window_end] / 1e6, 2),
      rolling_scale = scale_key
    )]
    n_inv <- sum(karyo_dt$call == "INV_INV")
    n_non <- sum(karyo_dt$call == "INV_nonINV")
    n_uniq <- length(unique(karyo_dt$sample_id))
    message("[S3v6b]   Karyotypes: ", nrow(karyo_dt), " runs (",
            n_inv, " INV/INV, ", n_non, " INV/nonINV) across ",
            n_uniq, " samples")
  } else {
    message("[S3v6b]   No stable karyotype runs")
  }

  # ── C+D: Interval classification + decomposition ──
  chr_interval_dt <- data.table()
  chr_decomp_dt   <- data.table()

  if (!is.null(intervals_dt)) {
    chr_intervals <- intervals_dt[chrom == chr]
    if (nrow(chr_intervals) > 0) {
      message("[S3v6b] Classifying ", nrow(chr_intervals), " intervals...")

      for (ii in seq_len(nrow(chr_intervals))) {
        intv <- chr_intervals[ii]
        intv_id <- if ("interval_id" %in% names(intv)) intv$interval_id else paste0("I", ii)
        intv_start <- intv$start_bp
        intv_end   <- intv$end_bp

        # Find windows within this interval
        win_mask <- window_info$start_bp >= intv_start & window_info$end_bp <= intv_end
        win_idx  <- which(win_mask)

        if (length(win_idx) < 5) {
          message("[S3v6b]   Interval ", intv_id, ": only ", length(win_idx),
                  " windows, skipping")
          next
        }

        message("[S3v6b]   Interval ", intv_id, " (",
                round(intv_start / 1e6, 2), "-", round(intv_end / 1e6, 2),
                " Mb, ", length(win_idx), " windows)")

        # C: Classify
        cl_dt <- classify_interval(roll_mat, win_idx, max_k = MAX_K)
        cl_dt[, `:=`(
          chrom       = chr,
          interval_id = intv_id,
          interval_start_bp = intv_start,
          interval_end_bp   = intv_end,
          interval_start_mb = round(intv_start / 1e6, 2),
          interval_end_mb   = round(intv_end / 1e6, 2),
          n_interval_windows = length(win_idx),
          rolling_scale = scale_key
        )]

        # Report
        class_tab <- table(cl_dt$interval_class)
        class_str <- paste(paste0(names(class_tab), "=", class_tab), collapse = " ")
        message("[S3v6b]     k=", cl_dt$interval_k[1],
                " sil=", round(cl_dt$silhouette[1], 3),
                " | ", class_str)

        chr_interval_dt <- rbind(chr_interval_dt, cl_dt, fill = TRUE)

        # D: Decompose
        dec_dt <- decompose_interval(roll_mat, win_idx, cl_dt$interval_class)
        if (nrow(dec_dt) > 0) {
          dec_dt[, `:=`(
            chrom       = chr,
            interval_id = intv_id,
            interval_start_bp = intv_start,
            interval_end_bp   = intv_end
          )]
          # Map changepoint global window index to bp
          dec_dt[, changepoint_bp := window_info$start_bp[changepoint_pos]]
          dec_dt[, changepoint_mb := round(changepoint_bp / 1e6, 2)]

          n_clusters <- length(unique(dec_dt$cp_cluster))
          if (n_clusters > 1) {
            message("[S3v6b]     Decomposition: ", n_clusters,
                    " sub-systems detected (changepoint clusters)")
            for (cc in sort(unique(dec_dt$cp_cluster))) {
              cc_dt <- dec_dt[cp_cluster == cc]
              message("[S3v6b]       Cluster ", cc, ": n=", nrow(cc_dt),
                      " median_cp=", round(median(cc_dt$changepoint_mb), 2), " Mb",
                      " mean_asymm=", round(mean(cc_dt$asymmetry), 3))
            }
          } else {
            message("[S3v6b]     No sub-system decomposition (single changepoint cluster)")
          }

          chr_decomp_dt <- rbind(chr_decomp_dt, dec_dt, fill = TRUE)
        }
      }
    }
  }

  # ===========================================================================
  # Chat 14 (2026-04-18): Per-chromosome RDS emit
  # ===========================================================================
  # Three per-chrom RDS files are written to <outdir>/annot/ and
  # <outdir>/per_sample/. These are what downstream consumers (phase-4b
  # Tier-3, run_all.R 2d block scoring, utils/lib_ghsl_panel.R) read.
  #
  # Files:
  #   annot/<chr>.ghsl_v6.annot.rds      — thin per-window aggregates
  #                                         (one row per window, no sample
  #                                         dimension). Consumed by
  #                                         run_all.R 2d block scoring.
  #   annot/<chr>.ghsl_v6.karyotypes.rds — per-sample stable LOW/HIGH runs
  #                                         (one row per run). Consumed by
  #                                         lib_ghsl_confirmation.R for
  #                                         Tier-3 SPLIT detection.
  #   per_sample/<chr>.ghsl_v6.per_sample.rds — dense long-format panel
  #                                         (one row per sample×window).
  #                                         Carries divergence and rank
  #                                         at all rolling scales. The
  #                                         primary artifact for on-demand
  #                                         queries via utils/lib_ghsl_panel.R
  #                                         and reg$compute$ghsl_*.
  # ===========================================================================
  annot_outdir <- file.path(outdir, "annot")
  panel_outdir <- file.path(outdir, "per_sample")
  dir.create(annot_outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(panel_outdir, recursive = TRUE, showWarnings = FALSE)

  # --- annot RDS (thin, per-window) ---
  # `metrics` already has: window_idx, n_scored, div_mean, div_median,
  # div_sd, div_iqr, div_skew, div_bimodal, n_low_div, n_high_div,
  # n_mid_div, low_div_tightness, high_div_tightness, mid_div_spread,
  # rank_stability, rank_stability_z, div_contrast, div_contrast_z,
  # ghsl_v6_score, ghsl_v6_status, chrom, global_window_id, start_bp,
  # end_bp, pos_mb, rolling_scale.
  annot_rds <- file.path(annot_outdir, paste0(chr, ".ghsl_v6.annot.rds"))
  saveRDS(metrics, annot_rds)
  message("[S3v6b]   annot RDS: ", annot_rds,
          " (", nrow(metrics), " windows)")

  # --- karyotypes RDS (per-sample stable runs) ---
  # karyo_dt already has: sample_id, window_start, window_end, n_windows,
  # call, mean_rank, chrom, start_bp, end_bp, start_mb, end_mb,
  # rolling_scale.
  karyo_rds <- file.path(annot_outdir, paste0(chr, ".ghsl_v6.karyotypes.rds"))
  saveRDS(karyo_dt, karyo_rds)
  message("[S3v6b]   karyotypes RDS: ", karyo_rds,
          " (", nrow(karyo_dt), " runs)")

  # --- per-sample panel RDS (dense, long-format, multi-scale) ---
  # Build one row per (sample, window) with divergence + rank at every
  # available rolling scale, plus window-level metrics (same for all
  # samples at a given window; duplicated for join convenience).
  #
  # Available scales come from mat_data$rolling — which is every scale
  # the heavy engine (STEP_C04) was invoked with. Default is
  # s10/s20/s30/s40/s50/s100.
  scale_keys <- names(mat_data$rolling)
  n_samp_pan <- length(snames)

  # Quick helper: per-sample rank at a scale (normalized 0..1 per window,
  # like rank_mat but recomputed per scale). Cache one matrix per scale.
  compute_rank_mat_for_scale <- function(roll_mat_s) {
    rm_s <- matrix(NA_real_,
                   nrow = nrow(roll_mat_s),
                   ncol = ncol(roll_mat_s),
                   dimnames = dimnames(roll_mat_s))
    for (wi in seq_len(ncol(roll_mat_s))) {
      vals <- roll_mat_s[, wi]
      valid <- !is.na(vals)
      if (sum(valid) < 10) next
      rm_s[valid, wi] <- rank(vals[valid], ties.method = "average") /
                          sum(valid)
    }
    rm_s
  }

  # Build the long-format panel in two passes to avoid memory peaks:
  # first a compact list per scale, then cbind at the end.
  # Base columns — repeat per sample, common across scales.
  win_idx_vec <- rep(seq_len(n_win), each = n_samp_pan)
  samp_vec    <- rep(snames, times = n_win)
  # Expand window info to sample×window rows (window is the faster-
  # changing axis because each window's samples are contiguous).
  gwid_vec    <- window_info$global_window_id[win_idx_vec]
  start_vec   <- window_info$start_bp[win_idx_vec]
  end_vec     <- window_info$end_bp[win_idx_vec]
  posmb_vec   <- window_info$pos_mb[win_idx_vec]

  panel_dt <- data.table(
    sample_id         = samp_vec,
    chrom             = chr,
    window_idx        = win_idx_vec,
    global_window_id  = gwid_vec,
    start_bp          = start_vec,
    end_bp            = end_vec,
    pos_mb            = posmb_vec
  )

  # Per-scale divergence + rank columns
  # Layout note: matrix roll_mat_s is [sample × window], and we want
  # the long-format ordering where each window's samples are
  # contiguous — matches samp_vec / win_idx_vec above. As.vector on an
  # R matrix walks column-major (samples-within-windows) which is
  # exactly the order we want, so we can assign directly.
  for (sk in scale_keys) {
    rm_s <- mat_data$rolling[[sk]]
    if (is.null(rm_s)) next
    panel_dt[[paste0("div_roll_", sk)]] <- as.vector(rm_s)
    rank_s <- compute_rank_mat_for_scale(rm_s)
    panel_dt[[paste0("rank_in_cohort_", sk)]] <- as.vector(rank_s)
    # Rank-band at this scale using the same thresholds as Part B
    rb <- rep(NA_character_, length(panel_dt$sample_id))
    rv <- panel_dt[[paste0("rank_in_cohort_", sk)]]
    rb[!is.na(rv) & rv <= KARYO_LO]                      <- "LOW"
    rb[!is.na(rv) & rv >= KARYO_HI]                      <- "HIGH"
    rb[!is.na(rv) & rv > KARYO_LO & rv < KARYO_HI]       <- "MID"
    panel_dt[[paste0("rank_band_", sk)]] <- rb
    rm(rm_s, rank_s, rb, rv)
  }

  # in_stable_run / stable_run_call — per-(sample, window) lookup based
  # on Part B's karyo_dt. karyo_dt uses window indices; expand each run
  # into (sample, window) cells we can join against.
  panel_dt[, in_stable_run  := FALSE]
  panel_dt[, stable_run_call := NA_character_]
  if (nrow(karyo_dt) > 0) {
    # Build a long data.table of (sample_id, window_idx, call) for
    # every (sample, window) that sits inside a stable run.
    expand_run <- function(r) {
      data.table(
        sample_id  = r$sample_id,
        window_idx = seq.int(r$window_start, r$window_end),
        call       = r$call
      )
    }
    run_cells <- rbindlist(lapply(seq_len(nrow(karyo_dt)),
                                    function(i) expand_run(karyo_dt[i])))
    setkey(run_cells, sample_id, window_idx)
    # Join back onto panel
    panel_key <- panel_dt[, .(sample_id, window_idx)]
    setkey(panel_key, sample_id, window_idx)
    joined <- run_cells[panel_key]
    panel_dt[, in_stable_run   := !is.na(joined$call)]
    panel_dt[, stable_run_call := joined$call]
    rm(run_cells, panel_key, joined)
  }

  # Window-level columns (same value for all samples at a given window,
  # but carrying them in the panel saves a join for most consumers).
  panel_dt[, ghsl_v6_score   := metrics$ghsl_v6_score[win_idx_vec]]
  panel_dt[, ghsl_v6_status  := metrics$ghsl_v6_status[win_idx_vec]]
  panel_dt[, rank_stability  := metrics$rank_stability[win_idx_vec]]
  panel_dt[, div_contrast_z  := metrics$div_contrast_z[win_idx_vec]]
  panel_dt[, div_bimodal     := metrics$div_bimodal[win_idx_vec]]

  # Attach metadata (scale list, thresholds, timestamp, sample order)
  panel_meta <- list(
    chrom            = chr,
    n_samples        = n_samp_pan,
    n_windows        = n_win,
    scales_available = scale_keys,
    primary_scale    = scale_key,
    karyo_lo         = KARYO_LO,
    karyo_hi         = KARYO_HI,
    karyo_min_run    = KARYO_MIN_RUN,
    sample_order     = snames,
    window_info      = window_info,
    generated_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    classifier       = "STEP_C04b_snake3_ghsl_classify.R (v6, chat-14 patched)"
  )
  attr(panel_dt, "ghsl_panel_meta") <- panel_meta

  panel_rds <- file.path(panel_outdir, paste0(chr, ".ghsl_v6.per_sample.rds"))
  saveRDS(panel_dt, panel_rds)
  fsize <- round(file.info(panel_rds)$size / 1e6, 1)
  message("[S3v6b]   per-sample panel RDS: ", panel_rds,
          " (", nrow(panel_dt), " rows, ", fsize, " MB)")

  rm(panel_dt, win_idx_vec, samp_vec, gwid_vec, start_vec, end_vec,
     posmb_vec)
  invisible(gc(verbose = FALSE, full = FALSE))

  # Collect
  all_metrics[[length(all_metrics) + 1L]]    <- metrics
  all_karyo[[length(all_karyo) + 1L]]        <- karyo_dt
  if (nrow(chr_interval_dt) > 0)
    all_interval[[length(all_interval) + 1L]] <- chr_interval_dt
  if (nrow(chr_decomp_dt) > 0)
    all_decomp[[length(all_decomp) + 1L]]     <- chr_decomp_dt

  all_summary[[length(all_summary) + 1L]] <- data.table(
    chrom = chr, n_windows = n_win, rolling_scale = scale_key,
    n_pass = n_pass, n_weak = n_weak, n_fail = n_fail,
    n_karyo_calls = nrow(karyo_dt),
    n_intervals_classified = nrow(chr_interval_dt),
    n_decomp_records = nrow(chr_decomp_dt)
  )

  elapsed <- round((proc.time() - t0)[3], 1)
  message("[S3v6b] ", chr, " DONE in ", elapsed, "s")
}

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

message("\n[S3v6b] Writing outputs...")

metrics_dt  <- if (length(all_metrics) > 0) rbindlist(all_metrics, fill = TRUE) else data.table()
karyo_dt    <- if (length(all_karyo) > 0) rbindlist(all_karyo, fill = TRUE) else data.table()
intv_dt     <- if (length(all_interval) > 0) rbindlist(all_interval, fill = TRUE) else data.table()
decomp_dt   <- if (length(all_decomp) > 0) rbindlist(all_decomp, fill = TRUE) else data.table()
summ_dt     <- if (length(all_summary) > 0) rbindlist(all_summary) else data.table()

f_track <- file.path(outdir, "snake3v6_window_track.tsv.gz")
fwrite(metrics_dt, f_track, sep = "\t")
message("[S3v6b] Window track: ", f_track, " (", nrow(metrics_dt), " rows)")

f_karyo <- file.path(outdir, "snake3v6_karyotype_calls.tsv.gz")
fwrite(karyo_dt, f_karyo, sep = "\t")
message("[S3v6b] Karyotype calls: ", f_karyo, " (", nrow(karyo_dt), " rows)")

if (nrow(intv_dt) > 0) {
  f_intv <- file.path(outdir, "snake3v6_interval_genotypes.tsv.gz")
  fwrite(intv_dt, f_intv, sep = "\t")
  message("[S3v6b] Interval genotypes: ", f_intv, " (", nrow(intv_dt), " rows)")
}

if (nrow(decomp_dt) > 0) {
  f_decomp <- file.path(outdir, "snake3v6_interval_decomp.tsv.gz")
  fwrite(decomp_dt, f_decomp, sep = "\t")
  message("[S3v6b] Decomposition: ", f_decomp, " (", nrow(decomp_dt), " rows)")
}

f_summ <- file.path(outdir, "snake3v6_summary.tsv")
fwrite(summ_dt, f_summ, sep = "\t")

message("\n================================================================")
message("[DONE] Snake 3 v6b: Light Classifier")
message("================================================================")
message("  Track:      ", f_track)
message("  Karyotype:  ", f_karyo)
if (nrow(intv_dt) > 0) message("  Intervals:  ", file.path(outdir, "snake3v6_interval_genotypes.tsv.gz"))
if (nrow(decomp_dt) > 0) message("  Decomp:     ", file.path(outdir, "snake3v6_interval_decomp.tsv.gz"))
message("  Summary:    ", f_summ)

if (nrow(summ_dt) > 0) {
  message("\n  Totals: PASS=", sum(summ_dt$n_pass),
          " WEAK=", sum(summ_dt$n_weak),
          " FAIL=", sum(summ_dt$n_fail))
}

message("\n  Tune parameters and rerun:")
message("    --scale 20|50|100     Change rolling scale")
message("    --karyo_lo 0.15       Adjust INV/INV threshold")
message("    --karyo_hi 0.70       Adjust INV/nonINV threshold")
message("    --karyo_min_run 10    Min consecutive windows for call")
message("    --max_k 5             Max clusters for interval classification")
