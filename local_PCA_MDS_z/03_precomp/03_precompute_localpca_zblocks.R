#!/usr/bin/env Rscript

# =============================================================================
# 03_precompute_localpca_zblocks.R
#
# Per-chromosome precomputation for path 1 (local-PCA z-blocks). Takes the
# unified MDS RDS from step 02b and emits ONE precomp.rds per chromosome
# plus the suite of NN-smoothed similarity matrices that downstream L1/L2
# detectors and the atlas JSON exporter all consume.
#
# Parallelized across chromosomes via mclapply. Run ONCE after step 02b;
# downstream detect/plot/export steps run interactively against the outputs.
#
# Pipeline position:
#   02b mds_merge  -> 03 PRECOMPUTE  -> 04 detect_L1
#                                    -> 05 plot_L1
#                                    -> 06 detect_L2
#                                    -> 07 plot_L2
#                                    -> 08b atlas_json
#
# ── This is the path-1 precompute only ───────────────────────────────────
# Sibling discovery paths each have their own precompute (they share the
# script skeleton but read different feature matrices):
#     path 1  local-PCA z-blocks  (this script — dosage feature)
#     path 2  local-PCA theta_pi  (path_localpca_thetapi/, pestPG feature)
#     path 3  local-PCA GHSL      (path_localpca_GHSL/, Clair3-phased feature)
# Each path writes to its own scratch tree and produces its own atlas JSON.
#
# ── Inputs ────────────────────────────────────────────────────────────────
#   <step02b_outprefix>     positional 1: prefix passed to step 02b
#                           (the script appends .mds.rds and reads from there)
#   <outdir>                positional 2: output root
#   [--dosage_dir <dir>]    if given, the dosage CV (uses the SAME SNPs as
#                           the local PCA) is used as the het component
#                           of inv_likeness, otherwise het_contrast is used
#
# ── Outputs (in <outdir>/) ────────────────────────────────────────────────
#   precomp/<chr>.precomp.rds                 FULL precomp (no slim)
#       $dt       per-window data.table including:
#           position    : global_window_id, chrom, start_bp, end_bp,
#                         center_bp, n_snps
#           MDS         : MDS1..MDSk, MDS{1..k}_z, max_abs_z, max_z_axis
#           inv_likeness: inv_likeness, inv_dip_p, inv_het_contrast
#           band shape  : band_discreteness, diffuse_score,
#                         het_intermediacy, n_effective_clusters
#           het sources : pc1_bimodality, het_pc1_gap, het_mid_fraction,
#                         het_mid_variance, dosage_het_rate_*
#           local k=3   : local_delta12, local_entropy, local_ena,
#                         band1_frac, band2_frac, band3_frac
#           adaptive    : beta_pval, adaptive_seed, beta_alpha, beta_beta
#           seed        : seed_nn_dist
#           morphology  : flat_inv_score, spiky_inv_score, fragmentation_score,
#                         + jaggedness/run/peak/plateau/nbhood/block columns
#           **per-sample**: PC_1_Ind*, PC_2_Ind*  (NEEDED by 08b atlas JSON)
#       $sim_mat                       window x window similarity (raw)
#       $mds_mat, $chrom, $n_windows
#       $bg_continuity_quantiles
#
#   precomp/sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
#                              MDS-space k-NN smoothed similarity matrices
#                              Tune via env var:
#                                NN_SIM_SCALES="40,80,160,320" sbatch ...
#
#   window_dt.tsv.gz           genome-wide per-window scalar table
#   precomp_summary.tsv        per-chrom window counts + timing
#
# ── No more slim ──────────────────────────────────────────────────────────
# Earlier versions wrote both precomp.rds and a precomp.slim.rds (same
# without per-sample PC_*_* columns). Slim broke step 08b (the JSON
# exporter falls back to PC2 jitter when PC_2_* is absent). Decision:
# emit ONE full precomp.rds; steps 04-07 only touch window coords so the
# extra columns are harmless to them, and 08b gets real PC1+PC2.
#
# ── What's NOT in this precompute ─────────────────────────────────────────
# Per the path-1 separation:
#   * NO SV prior stamping        (was path 2 — handled at phase 4 catalog birth)
#   * NO Q-matrix / family_likeness / localQ_* / band_Q_*    (separate overlay)
#   * NO test_05 Fst              (phase-4 stamp)
#   * NO relatedness/kinship matrix loading (08a builds sample_metadata.tsv)
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: STEP_C01a_precompute.R v10.0 SLIM)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript STEP_C01a_precompute.R <step10_outprefix> <outdir> [--dosage_dir <dir>]")
}

step10_prefix <- args[1]
outdir        <- args[2]
precomp_dir   <- file.path(outdir, "precomp")
dir.create(precomp_dir, recursive = TRUE, showWarnings = FALSE)

# Optional args
dosage_dir <- NULL
i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--dosage_dir" && i < length(args)) { dosage_dir <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

# =============================================================================
# PARAMETERS
# =============================================================================

SEED_MDS_AXES    <- 5L    # number of MDS axes contributing to max_abs_z
SEED_NEIGHBOR_K  <- 3L    # k for seed_nn_dist

# =============================================================================
# LOAD MDS
# =============================================================================

mds_rds_file <- paste0(step10_prefix, ".mds.rds")
if (!file.exists(mds_rds_file)) stop("Missing: ", mds_rds_file)

message("[PRECOMP] v10.0 SLIM — local PCA z-outlier path only")
message("[PRECOMP] Loading ", mds_rds_file, " ...")
t_load <- proc.time()
mds_obj <- readRDS(mds_rds_file)
message("[PRECOMP] Loaded in ", round((proc.time() - t_load)[3], 1), "s")

per_chr <- mds_obj$per_chr
if (is.null(per_chr) || length(per_chr) == 0) stop("No per_chr data")
chroms <- names(per_chr)
message("[PRECOMP] ", length(chroms), " chromosomes")

# Extract sample names
sample_names <- NULL
for (chr_tmp in chroms) {
  dt_tmp <- as.data.table(per_chr[[chr_tmp]]$out_dt)
  pc1_cols <- grep("^PC_1_", names(dt_tmp), value = TRUE)
  if (length(pc1_cols) > 0) {
    sample_names <- sub("^PC_1_", "", pc1_cols)
    break
  }
}
if (!is.null(sample_names)) {
  message("[PRECOMP] Sample names: ", length(sample_names))
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

make_sim_mat <- function(dmat) {
  finite_vals <- dmat[is.finite(dmat)]
  dmax <- if (length(finite_vals) > 0) quantile(finite_vals, 0.95, na.rm = TRUE) else 1
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  sim <- 1 - pmin(dmat / dmax, 1)
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  sim
}

# =============================================================================
# DOSAGE-BASED HET RATE
# =============================================================================
# For each window, count sites per sample where dosage is het-like (0.3-1.7).
# This uses the SAME SNPs as the PCA (so it's biologically coupled to the
# local locus) but is computed from raw dosage values, not eigenvectors —
# so it is computationally independent from the PCA decomposition.
#
# A window with many samples showing high het rate = real heterozygosity
# = consistent with inversion het genotype.
#
# Returns: data.table with global_window_id, dosage_het_rate_median,
#          dosage_het_rate_sd, dosage_het_rate_cv (CV: high = bimodal het)

compute_dosage_het_rates <- function(dosage_dir, per_chr, chroms) {
  if (is.null(dosage_dir) || !dir.exists(dosage_dir)) return(data.table())
  all_rows <- list()

  for (chr in chroms) {
    dos_file <- file.path(dosage_dir, paste0(chr, ".dosage.tsv.gz"))
    if (!file.exists(dos_file)) {
      dos_file <- file.path(dosage_dir, paste0(chr, ".dosage.tsv"))
      if (!file.exists(dos_file)) next
    }

    dt <- as.data.table(per_chr[[chr]]$out_dt)
    if (nrow(dt) == 0 || !"start_bp" %in% names(dt)) next

    dos <- tryCatch(fread(dos_file, nThread = 4L), error = function(e) NULL)
    if (is.null(dos) || nrow(dos) < 100) {
      message("  [dosage] ", chr, ": cannot load or too few sites")
      next
    }

    pos_col <- intersect(c("POS", "pos", "position"), names(dos))
    if (length(pos_col) == 0) { message("  [dosage] ", chr, ": no POS column"); next }
    positions <- dos[[pos_col[1]]]
    dos_mat <- as.matrix(dos[, !..pos_col[1]])

    message("  [dosage] ", chr, ": ", nrow(dos_mat), " sites x ", ncol(dos_mat), " samples")

    for (wi in seq_len(nrow(dt))) {
      w_start <- dt$start_bp[wi]; w_end <- dt$end_bp[wi]
      wid <- dt$global_window_id[wi]

      site_idx <- which(positions >= w_start & positions <= w_end)
      if (length(site_idx) < 10) next

      w_dos <- dos_mat[site_idx, , drop = FALSE]

      # Per-sample het rate: fraction of sites with dosage in [0.3, 1.7]
      per_sample_het <- colMeans(w_dos >= 0.3 & w_dos <= 1.7, na.rm = TRUE)
      per_sample_het <- per_sample_het[is.finite(per_sample_het)]
      if (length(per_sample_het) < 20) next

      het_median <- median(per_sample_het)
      het_sd <- sd(per_sample_het)
      het_cv <- if (het_median > 0.01) het_sd / het_median else NA_real_

      # Bimodality of het rate distribution = inversion signal:
      # inversion → bimodal (het samples ~50%, hom samples ~10%) → high CV
      # background → unimodal → low CV
      all_rows[[length(all_rows) + 1]] <- data.table(
        chrom = chr, global_window_id = wid,
        dosage_het_rate_median = round(het_median, 4),
        dosage_het_rate_sd = round(het_sd, 4),
        dosage_het_rate_cv = round(het_cv %||% NA_real_, 4)
      )
    }
    message("  [dosage] ", chr, ": ", length(all_rows), " windows processed")
  }
  if (length(all_rows) > 0) rbindlist(all_rows) else data.table()
}

# =============================================================================
# INVERSION-LIKENESS (genome-wide pass)
# =============================================================================

compute_inv_likeness_all <- function(per_chr, chroms) {
  all_rows <- list()
  for (chr in chroms) {
    dt <- as.data.table(per_chr[[chr]]$out_dt)
    if (nrow(dt) == 0) next
    pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
    has_pc1 <- length(pc1_cols) > 0

    # ── CHROMOSOME-WIDE FIXED k=3 BANDS ──
    # Computed ONCE per chromosome from average PC1 across all windows.
    # Provides stable group labels (avoids per-window k-means label switching).
    chr_bands <- NULL
    if (has_pc1 && length(pc1_cols) >= 20) {
      pc1_mat <- as.matrix(dt[, ..pc1_cols])
      avg_pc1 <- colMeans(pc1_mat, na.rm = TRUE)
      valid <- is.finite(avg_pc1)
      if (sum(valid) >= 20) {
        avg_valid <- avg_pc1[valid]
        chr_km <- tryCatch(kmeans(avg_valid, centers = 3, nstart = 10),
                            error = function(e) NULL)
        if (!is.null(chr_km)) {
          co <- order(chr_km$centers[, 1])
          chr_bands <- integer(length(avg_valid))
          chr_bands[chr_km$cluster == co[1]] <- 1L  # low PC1
          chr_bands[chr_km$cluster == co[2]] <- 2L  # mid PC1 (het-like)
          chr_bands[chr_km$cluster == co[3]] <- 3L  # high PC1
          names(chr_bands) <- sub("^PC_1_", "", names(avg_valid))
          message("  [", chr, "] Fixed bands: ",
                  sum(chr_bands == 1), "/", sum(chr_bands == 2), "/", sum(chr_bands == 3),
                  " (band1/band2/band3)")
        }
      }
    }

    for (wi in seq_len(nrow(dt))) {
      wid <- dt$global_window_id[wi]
      dip_p <- NA_real_; het_contrast <- NA_real_

      # Source B (fixed bands)
      het_pc1_gap <- NA_real_
      het_mid_fraction <- NA_real_
      het_mid_variance <- NA_real_
      pc1_bimodality <- NA_real_
      local_delta12 <- NA_real_
      local_entropy <- NA_real_
      local_ena <- NA_real_
      band_sizes <- rep(NA_real_, 3)

      # Source C (filled by post-loop dosage merge)
      dosage_het_rate_median <- NA_real_
      dosage_het_rate_sd <- NA_real_
      dosage_het_rate_cv <- NA_real_

      km3 <- NULL
      co <- NULL

      if (has_pc1) {
        scores <- as.numeric(dt[wi, ..pc1_cols])
        scores <- scores[is.finite(scores)]
        if (length(scores) >= 20) {
          # Dip test for trimodality
          dip_p <- tryCatch({
            if (requireNamespace("diptest", quietly = TRUE)) {
              diptest::dip.test(scores)$p.value
            } else {
              km <- tryCatch(kmeans(scores, centers = 3, nstart = 3), error = function(e) NULL)
              if (!is.null(km)) {
                centers <- sort(km$centers[, 1])
                gaps <- diff(centers)
                within_sd <- sqrt(mean(km$withinss / km$size))
                if (within_sd > 0) 1 / (1 + (min(gaps) / within_sd)^2) else 1
              } else 1
            }
          }, error = function(e) NA_real_)

          # HET CONTRAST CHECK (Source A — per-window k=3)
          km3 <- tryCatch(kmeans(scores, centers = 3, nstart = 5), error = function(e) NULL)
          if (!is.null(km3)) {
            co <- order(km3$centers[, 1])
            var_per_group <- vapply(co, function(g) {
              vals <- scores[km3$cluster == g]
              if (length(vals) >= 3) var(vals) else NA_real_
            }, numeric(1))
            sizes <- km3$size[co]
            min_frac <- min(sizes) / length(scores)
            between_var <- var(km3$centers[, 1])
            within_var <- mean(var_per_group, na.rm = TRUE)
            if (is.finite(between_var) && is.finite(within_var) && within_var > 0) {
              het_contrast <- between_var / within_var
            }
            if (min_frac < 0.05) het_contrast <- het_contrast * 0.5

            # Source B: chromosome-wide fixed bands applied to this window
            if (!is.null(chr_bands)) {
              sample_ids <- sub("^PC_1_", "", pc1_cols)
              all_scores <- as.numeric(dt[wi, ..pc1_cols])
              names(all_scores) <- sample_ids

              matched <- intersect(names(all_scores), names(chr_bands))
              if (length(matched) >= 20) {
                s <- all_scores[matched]
                b <- chr_bands[matched]

                tab <- table(factor(b, levels = 1:3))
                band_sizes <- as.numeric(tab) / sum(tab)
                het_mid_fraction <- band_sizes[2]

                mid_s <- s[b == 2]
                het_mid_variance <- if (length(mid_s) >= 3) var(mid_s) else NA_real_

                band_means <- tapply(s, b, mean, na.rm = TRUE)
                if (length(band_means) >= 2) {
                  sorted_bm <- sort(band_means)
                  het_pc1_gap <- min(diff(sorted_bm))
                }

                outer_a <- s[b == 1]; outer_b <- s[b == 3]
                if (length(outer_a) >= 3 && length(outer_b) >= 3) {
                  mu_diff <- abs(mean(outer_a) - mean(outer_b))
                  pooled_sd <- sqrt((var(outer_a) + var(outer_b)) / 2)
                  pc1_bimodality <- if (pooled_sd > 0) mu_diff / pooled_sd else NA_real_
                }

                eps_e <- 1e-15
                bm_vec <- as.numeric(band_means[as.character(1:3)])
                if (all(is.finite(bm_vec)) && diff(range(bm_vec)) > 0) {
                  spread <- diff(range(bm_vec)) * 0.01
                  soft <- matrix(0, nrow = length(s), ncol = 3)
                  for (gi in 1:3) soft[, gi] <- 1 / (abs(s - bm_vec[gi]) + spread)
                  soft <- soft / rowSums(soft)
                  per_H <- -rowSums(soft * log(soft + eps_e))
                  local_entropy <- mean(per_H, na.rm = TRUE)
                  local_ena <- exp(local_entropy)
                  per_d12 <- apply(soft, 1, function(p) {
                    ps <- sort(p, decreasing = TRUE); ps[1] - ps[2]
                  })
                  local_delta12 <- mean(per_d12, na.rm = TRUE)
                }
              }
            } else {
              # Fallback: per-window km3
              band_sizes <- sizes / sum(sizes)
              het_mid_fraction <- sizes[2] / sum(sizes)
              mid_vals <- scores[km3$cluster == co[2]]
              het_mid_variance <- if (length(mid_vals) >= 3) var(mid_vals) else NA_real_
              sorted_centers <- sort(km3$centers[, 1])
              het_pc1_gap <- min(diff(sorted_centers))
            }
          }
        }
      }

      # ── BAND-SHAPE DESCRIPTORS ──
      band_discreteness <- NA_real_
      diffuse_score <- NA_real_
      het_intermediacy <- NA_real_
      n_effective_clusters <- NA_real_

      if (!is.null(km3) && has_pc1) {
        scores_all <- as.numeric(dt[wi, ..pc1_cols])
        scores_all <- scores_all[is.finite(scores_all)]
        if (length(scores_all) >= 20) {
          co_local <- order(km3$centers[, 1])
          sizes_local <- km3$size[co_local]
          centers_sorted <- sort(km3$centers[, 1])
          var_per_g <- vapply(co_local, function(g) {
            vals <- scores_all[km3$cluster == g]
            if (length(vals) >= 3) var(vals) else NA_real_
          }, numeric(1))

          # band_discreteness: between-group gap / within-group spread
          mean_within_sd <- sqrt(mean(var_per_g, na.rm = TRUE))
          min_gap <- min(diff(centers_sorted))
          if (is.finite(mean_within_sd) && mean_within_sd > 0) {
            band_discreteness <- min_gap / mean_within_sd
          }

          # diffuse_score: 1 - trimodality (high p_dip = unimodal/diffuse)
          min_frac_local <- min(sizes_local) / length(scores_all)
          dip_contrib <- if (is.finite(dip_p)) dip_p else 0.5
          size_penalty <- if (min_frac_local < 0.08) 0.3 else 0
          diffuse_score <- pmin(1, dip_contrib + size_penalty)

          # het_intermediacy: is mid band truly between outers?
          mid_mean <- centers_sorted[2]
          outer_midpoint <- (centers_sorted[1] + centers_sorted[3]) / 2
          full_gap <- centers_sorted[3] - centers_sorted[1]
          if (is.finite(full_gap) && full_gap > 0) {
            het_intermediacy <- 1 - pmin(1, abs(mid_mean - outer_midpoint) / (full_gap / 2))
          }

          # n_effective_clusters: gap-stat-proxy via BSS/TSS for k=2/3/4
          tss <- sum((scores_all - mean(scores_all))^2)
          if (tss > 0) {
            bss_k <- numeric(4)
            bss_k[1] <- 0
            bss_k[3] <- sum(km3$betweenss) / tss
            for (kk in c(2L, 4L)) {
              kmk <- tryCatch(kmeans(scores_all, centers = kk, nstart = 3),
                              error = function(e) NULL)
              bss_k[kk] <- if (!is.null(kmk)) sum(kmk$betweenss) / tss else 0
            }
            gains <- diff(bss_k)
            n_effective_clusters <- 1
            for (kk in 2:4) {
              if (gains[kk - 1] > 0.05) n_effective_clusters <- kk
            }
          }
        }
      }

      # ── INV_LIKENESS (composite) ──
      # 45% HET + 30% TRIMODALITY + 25% BAND DISCRETENESS
      # Pass 1: het component = het_contrast (PCA-derived).
      # Pass 2 (post-loop): if dosage CV available, replaces het component.
      s_dip <- if (is.finite(dip_p)) pmin(1, pmax(0, (0.5 - dip_p) / 0.45)) else 0
      s_het <- 0
      if (is.finite(het_contrast)) {
        s_het <- pmin(1, pmax(0, (het_contrast - 2) / 8))
      }
      s_discrete <- 0
      if (is.finite(band_discreteness)) {
        s_discrete <- pmin(1, pmax(0, (band_discreteness - 1) / 4))
      }
      inv_like <- 0.45 * s_het + 0.30 * s_dip + 0.25 * s_discrete

      all_rows[[length(all_rows) + 1]] <- data.table(
        chrom = chr, global_window_id = wid,
        start_bp = dt$start_bp[wi], end_bp = dt$end_bp[wi],
        inv_dip_p = round(dip_p, 4),
        inv_het_contrast = round(het_contrast, 2),
        inv_likeness = round(inv_like, 4),
        # Band-shape descriptors
        band_discreteness = round(band_discreteness %||% NA_real_, 4),
        diffuse_score = round(diffuse_score %||% NA_real_, 4),
        het_intermediacy = round(het_intermediacy %||% NA_real_, 4),
        n_effective_clusters = as.integer(n_effective_clusters %||% NA_integer_),
        # Multi-source het assessment
        pc1_bimodality = round(pc1_bimodality %||% NA_real_, 4),
        het_pc1_gap = round(het_pc1_gap %||% NA_real_, 4),
        het_mid_fraction = round(het_mid_fraction %||% NA_real_, 4),
        het_mid_variance = round(het_mid_variance %||% NA_real_, 6),
        # Local fixed-band metrics
        local_delta12 = round(local_delta12 %||% NA_real_, 4),
        local_entropy = round(local_entropy %||% NA_real_, 4),
        local_ena = round(local_ena %||% NA_real_, 4),
        # Band sizes
        band1_frac = round(band_sizes[1] %||% NA_real_, 4),
        band2_frac = round(band_sizes[2] %||% NA_real_, 4),
        band3_frac = round(band_sizes[3] %||% NA_real_, 4),
        # Dosage-based het (filled by post-loop merge if --dosage_dir provided)
        dosage_het_rate_median = round(dosage_het_rate_median %||% NA_real_, 4),
        dosage_het_rate_sd = round(dosage_het_rate_sd %||% NA_real_, 4),
        dosage_het_rate_cv = round(dosage_het_rate_cv %||% NA_real_, 4)
      )
    }
  }
  if (length(all_rows) > 0) rbindlist(all_rows) else data.table()
}

message("[PRECOMP] Computing inversion-likeness scores...")
t_il <- proc.time()
window_dt <- compute_inv_likeness_all(per_chr, chroms)
message("[PRECOMP] Inv-likeness: ", nrow(window_dt), " windows, ",
        sum(window_dt$inv_likeness >= 0.5, na.rm = TRUE), " with score >= 0.5 (",
        round((proc.time() - t_il)[3], 1), "s)")

# ── Dosage-based het rate (computationally independent from PCA) ──
dosage_het_dt <- data.table()
if (!is.null(dosage_dir)) {
  message("[PRECOMP] Computing dosage het rates from ", dosage_dir, " ...")
  t_dos <- proc.time()
  dosage_het_dt <- compute_dosage_het_rates(dosage_dir, per_chr, chroms)
  message("[PRECOMP] Dosage het: ", nrow(dosage_het_dt), " windows (",
          round((proc.time() - t_dos)[3], 1), "s)")

  if (nrow(dosage_het_dt) > 0) {
    window_dt[dosage_het_dt, on = c("chrom", "global_window_id"),
      `:=`(dosage_het_rate_median = i.dosage_het_rate_median,
           dosage_het_rate_sd = i.dosage_het_rate_sd,
           dosage_het_rate_cv = i.dosage_het_rate_cv)]

    n_filled <- sum(is.finite(window_dt$dosage_het_rate_cv))
    message("[PRECOMP] Dosage merge: ", n_filled, " windows got dosage values")

    # SECOND PASS: recompute inv_likeness with dosage CV as het component
    n_with_dosage <- sum(is.finite(window_dt$dosage_het_rate_cv) &
                          window_dt$dosage_het_rate_cv > 0)
    if (n_with_dosage > 0) {
      message("[PRECOMP] Upgrading inv_likeness for ", n_with_dosage,
              " windows with dosage het CV...")
      window_dt[is.finite(dosage_het_rate_cv) & dosage_het_rate_cv > 0,
        inv_likeness := {
          s_h <- pmin(1, pmax(0, (dosage_het_rate_cv - 0.1) / 0.6))
          s_d <- pmin(1, pmax(0, (0.5 - inv_dip_p) / 0.45))
          s_d <- ifelse(is.finite(s_d), s_d, 0)
          s_bd <- pmin(1, pmax(0, (band_discreteness - 1) / 4))
          s_bd <- ifelse(is.finite(s_bd), s_bd, 0)
          round(0.45 * s_h + 0.30 * s_d + 0.25 * s_bd, 4)
        }]
      message("[PRECOMP] Done: ", n_with_dosage, " windows upgraded to dosage-based inv_likeness")
    }
  }
} else {
  message("[PRECOMP] No --dosage_dir provided — dosage het rate columns will be NA")
}

# ── Track summary ──
n_inv <- sum(window_dt$inv_likeness >= 0.5, na.rm = TRUE)
n_dos <- sum(is.finite(window_dt$dosage_het_rate_median))
n_discrete <- sum(window_dt$band_discreteness > 2, na.rm = TRUE)
n_diffuse <- sum(window_dt$diffuse_score > 0.5, na.rm = TRUE)
message("[PRECOMP] ── TRACK SUMMARY ──")
message("  Inv-like >= 0.5:       ", n_inv, " / ", nrow(window_dt),
        "  (formula: 45% het + 30% trimodality + 25% band_discreteness)")
message("  Band discrete (>2):    ", n_discrete, " / ", nrow(window_dt))
message("  Diffuse (>0.5):        ", n_diffuse, " / ", nrow(window_dt))
message("  Dosage het computed:   ", n_dos, " / ", nrow(window_dt))

# ═══════════════════════════════════════════════════════════════════
# TEST 07: Beta adaptive threshold
# ═══════════════════════════════════════════════════════════════════
# Fit Beta(α,β) to per-chromosome inv_likeness distribution.
# Background windows follow the bulk; inversion windows sit in the
# extreme right tail. The p-value gives a chromosome-specific threshold.
#
# Adds: beta_pval, adaptive_seed, beta_alpha, beta_beta

beta_adaptive_pvalues <- function(scores) {
  valid <- scores[is.finite(scores) & scores > 0.001 & scores < 0.999]
  if (length(valid) < 50) {
    return(data.table(beta_pval = rep(NA_real_, length(scores)),
                       adaptive_seed = rep(FALSE, length(scores)),
                       beta_alpha = NA_real_, beta_beta = NA_real_))
  }
  fit <- tryCatch({
    if (requireNamespace("MASS", quietly = TRUE)) {
      MASS::fitdistr(valid, "beta", start = list(shape1 = 1, shape2 = 5))
    } else {
      m <- mean(valid); v <- var(valid)
      if (v >= m * (1 - m)) v <- m * (1 - m) * 0.99
      a <- m * (m * (1 - m) / v - 1)
      b <- (1 - m) * (m * (1 - m) / v - 1)
      list(estimate = c(shape1 = max(0.1, a), shape2 = max(0.1, b)))
    }
  }, error = function(e) list(estimate = c(shape1 = 1, shape2 = 5)))
  a <- fit$estimate["shape1"]; b <- fit$estimate["shape2"]
  pvals <- 1 - pbeta(scores, a, b)
  pvals[!is.finite(scores)] <- NA_real_
  data.table(beta_pval = pvals, adaptive_seed = pvals < 0.01 & is.finite(pvals),
             beta_alpha = round(a, 3), beta_beta = round(b, 3))
}

message("[PRECOMP] ── TEST 07: Beta adaptive threshold ──")
window_dt[, c("beta_pval", "adaptive_seed", "beta_alpha", "beta_beta") :=
              list(NA_real_, FALSE, NA_real_, NA_real_)]

for (chr_i in unique(window_dt$chrom)) {
  idx <- window_dt$chrom == chr_i
  scores_chr <- window_dt$inv_likeness[idx]
  bt <- beta_adaptive_pvalues(scores_chr)
  window_dt[idx, `:=`(
    beta_pval = bt$beta_pval,
    adaptive_seed = bt$adaptive_seed,
    beta_alpha = bt$beta_alpha,
    beta_beta = bt$beta_beta
  )]
  n_seeds <- sum(bt$adaptive_seed, na.rm = TRUE)
  message("  ", chr_i, ": Beta(", bt$beta_alpha[1], ",", bt$beta_beta[1],
          ") → ", n_seeds, " adaptive seeds (p<0.01)")
}

n_total_seeds <- sum(window_dt$adaptive_seed, na.rm = TRUE)
message("[PRECOMP] Total adaptive seeds: ", n_total_seeds, " / ", nrow(window_dt))

f_inv <- file.path(outdir, "window_dt.tsv.gz")
fwrite(window_dt, f_inv, sep = "\t")
message("[PRECOMP] → ", f_inv)

# =============================================================================
# PER-CHROMOSOME PRECOMPUTE (parallelized with mclapply)
# =============================================================================

N_CORES <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")), length(chroms))
if (N_CORES > 1) {
  message("[PRECOMP] Using ", N_CORES, " cores for parallel precompute")
}

precompute_one_chr <- function(chr) {
  t_chr <- proc.time()
  chr_obj <- per_chr[[chr]]
  if (is.null(chr_obj)) return(NULL)

  dt <- as.data.table(chr_obj$out_dt)
  setalloccol(dt)   # refresh selfref pointer after RDS deserialization
  dt <- dt[order(start_bp)]
  dmat <- chr_obj$dmat
  mds_cols <- grep("^MDS[0-9]+$", names(dt), value = TRUE)
  mds_mat <- as.matrix(dt[, ..mds_cols])

  # Merge ALL inv_likeness + het + adaptive metrics into per-chr dt
  if (nrow(window_dt) > 0 && "global_window_id" %in% names(dt)) {
    il_cols <- setdiff(names(window_dt), c(names(dt), "chrom"))
    il_cols <- c("global_window_id", il_cols)
    dt <- merge(dt, window_dt[, ..il_cols, with = FALSE],
                by = "global_window_id", all.x = TRUE)
    dt <- dt[order(start_bp)]
    message("[PRECOMP] ", chr, ": merged ", length(il_cols) - 1, " inv_like columns into dt")
  }

  # Robust z-scores (median/MAD) on each MDS axis
  # Use data.table::set() (not dt[[zc]] <- ...) — the [[<- assignment
  # invalidates the .internal.selfref pointer and triggers the harmless-
  # but-noisy "Invalid .internal.selfref detected" warning on the next :=.
  for (mc in mds_cols) {
    zc <- paste0(mc, "_z")
    vv <- dt[[mc]]
    med <- median(vv, na.rm = TRUE)
    mad_val <- mad(vv, na.rm = TRUE)
    if (is.finite(mad_val) && mad_val > 1e-10) {
      set(dt, j = zc, value = (vv - med) / mad_val)
    } else {
      sdev <- sd(vv, na.rm = TRUE)
      if (is.finite(sdev) && sdev > 0) {
        set(dt, j = zc, value = (vv - mean(vv, na.rm = TRUE)) / sdev)
      } else {
        set(dt, j = zc, value = 0)
      }
    }
  }

  z_cols <- grep("^MDS[0-9]+_z$", names(dt), value = TRUE)
  if (length(z_cols) > SEED_MDS_AXES) z_cols <- z_cols[seq_len(SEED_MDS_AXES)]
  if (length(z_cols) > 0) {
    dt[, max_abs_z := apply(.SD, 1, function(x) max(abs(x), na.rm = TRUE)), .SDcols = z_cols]
    # Which axis dominated? Useful for the JSON exporter and atlas pages.
    dt[, max_z_axis := apply(.SD, 1, function(x) {
      ax <- which.max(abs(x))
      if (length(ax) == 0) NA_integer_ else as.integer(ax)
    }), .SDcols = z_cols]
  } else {
    dt[, max_abs_z := 0]
    dt[, max_z_axis := NA_integer_]
  }

  # Align dimensions
  n_dt <- nrow(dt); n_dm <- nrow(dmat)
  if (n_dm != n_dt) {
    n <- min(n_dm, n_dt)
    dt <- dt[seq_len(n)]; dmat <- dmat[seq_len(n), seq_len(n), drop = FALSE]
    mds_mat <- mds_mat[seq_len(n), , drop = FALSE]
  }

  # Similarity matrix
  sim_mat <- make_sim_mat(dmat)

  # ── BACKGROUND CONTINUITY BASELINE ──
  bg_continuity <- numeric(0)
  chr_inv_med <- median(dt$inv_likeness[is.finite(dt$inv_likeness)])
  n_w <- nrow(dt)

  if (n_w >= 10 && is.finite(chr_inv_med)) {
    bg_pairs <- which(
      dt$inv_likeness[seq_len(n_w - 1)] < chr_inv_med &
      dt$inv_likeness[seq(2, n_w)] < chr_inv_med
    )
    if (length(bg_pairs) >= 20) {
      if (length(bg_pairs) > 2000) bg_pairs <- sort(sample(bg_pairs, 2000))
      bg_continuity <- vapply(bg_pairs, function(i) {
        j <- i + 1L
        if (is.finite(sim_mat[i, j])) sim_mat[i, j] else 0
      }, numeric(1))
    }
  }

  if (length(bg_continuity) >= 10) {
    bg_q <- quantile(bg_continuity, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
  } else {
    adj_sims <- vapply(seq_len(n_w - 1), function(i) sim_mat[i, i + 1L], numeric(1))
    bg_q <- quantile(adj_sims, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
  }

  message("[PRECOMP] ", chr, ": bg_continuity baseline: ",
          "q50=", round(bg_q["50%"], 3),
          " q85=", round(bg_q["85%"], 3),
          " q90=", round(bg_q["90%"], 3),
          " q95=", round(bg_q["95%"], 3),
          " (from ", length(bg_continuity), " bg pairs)")

  # Seed NN distances
  nn_dists <- vapply(seq_len(nrow(dmat)), function(i) {
    d <- dmat[i, ]; d[i] <- Inf; d <- d[is.finite(d)]
    if (length(d) == 0) return(Inf)
    mean(sort(d)[seq_len(min(SEED_NEIGHBOR_K, length(d)))], na.rm = TRUE)
  }, numeric(1))
  dmax_nn <- quantile(dmat[is.finite(dmat)], 0.95, na.rm = TRUE)
  if (!is.finite(dmax_nn) || dmax_nn == 0) dmax_nn <- 1
  dt[, seed_nn_dist := nn_dists / dmax_nn]

  # ═══════════════════════════════════════════════════════════════════════
  # MORPHOLOGY FEATURES
  # ═══════════════════════════════════════════════════════════════════════

  n_w <- nrow(dt)
  z_vec <- dt$max_abs_z

  if (n_w >= 5 && any(is.finite(z_vec) & z_vec > 0)) {
    z_vec[!is.finite(z_vec)] <- 0

    # 1a. Local jaggedness
    JAGGED_HALF <- 5L
    jagged_vec <- rep(NA_real_, n_w)
    for (i in seq_len(n_w)) {
      lo <- max(2L, i - JAGGED_HALF)
      hi <- min(n_w, i + JAGGED_HALF)
      local_z <- z_vec[lo:hi]
      if (length(local_z) >= 3) {
        jagged_vec[i] <- mean(abs(diff(local_z)), na.rm = TRUE)
      }
    }
    dt[, local_jaggedness := round(jagged_vec, 4)]

    # 1b. Elevated runs
    z_thresh <- max(quantile(z_vec, 0.75, na.rm = TRUE), 1.5)
    is_elevated <- z_vec >= z_thresh

    run_id <- integer(n_w)
    current_run <- 0L
    for (i in seq_len(n_w)) {
      if (is_elevated[i]) {
        if (i == 1L || !is_elevated[i - 1L]) current_run <- current_run + 1L
        run_id[i] <- current_run
      }
    }

    run_len_vec   <- rep(0L, n_w)
    run_mean_z    <- rep(NA_real_, n_w)
    run_sd_z      <- rep(NA_real_, n_w)
    run_start_vec <- rep(NA_integer_, n_w)
    run_end_vec   <- rep(NA_integer_, n_w)

    if (current_run > 0) {
      for (rid in seq_len(current_run)) {
        idx <- which(run_id == rid)
        rlen <- length(idx)
        rmean <- mean(z_vec[idx], na.rm = TRUE)
        rsd   <- if (rlen >= 2) sd(z_vec[idx], na.rm = TRUE) else 0
        rstart <- idx[1]; rend <- idx[rlen]
        run_len_vec[idx]   <- rlen
        run_mean_z[idx]    <- rmean
        run_sd_z[idx]      <- rsd
        run_start_vec[idx] <- rstart
        run_end_vec[idx]   <- rend
      }
    }

    dt[, local_run_len     := run_len_vec]
    dt[, local_run_mean_z  := round(run_mean_z, 4)]
    dt[, local_run_sd_z    := round(run_sd_z, 4)]
    dt[, run_start_idx     := run_start_vec]
    dt[, run_end_idx       := run_end_vec]
    dt[, candidate_length_windows := local_run_len]
    dt[, candidate_length_bp := fifelse(
      is.finite(run_start_idx) & is.finite(run_end_idx),
      as.numeric(dt$end_bp[pmin(run_end_vec, n_w)] - dt$start_bp[pmax(run_start_vec, 1L)]),
      NA_real_
    )]

    # 1c. Peak prominence
    PROM_HALF <- 20L
    prom_vec <- rep(NA_real_, n_w)
    for (i in seq_len(n_w)) {
      lo <- max(1L, i - PROM_HALF)
      hi <- min(n_w, i + PROM_HALF)
      bg_z <- z_vec[lo:hi]
      bg_level <- quantile(bg_z, 0.25, na.rm = TRUE)
      prom_vec[i] <- z_vec[i] - bg_level
    }
    dt[, local_peak_prominence := round(prom_vec, 4)]

    # 1d. Plateau flatness
    pf_vec <- rep(0, n_w)
    for (i in seq_len(n_w)) {
      rl <- run_len_vec[i]
      if (rl >= 3) {
        sd_term  <- 1 / (1 + (run_sd_z[i] %||% 0))
        jag_term <- 1 / (1 + (jagged_vec[i] %||% 0))
        pf_vec[i] <- rl * sd_term * jag_term
      }
    }
    pf_max <- max(pf_vec, na.rm = TRUE)
    if (is.finite(pf_max) && pf_max > 0) pf_vec <- pf_vec / pf_max
    dt[, plateau_flatness := round(pf_vec, 4)]

    # 2. Neighborhood support
    NBHOOD_RADII <- c(5L, 10L, 20L)
    for (rad in NBHOOD_RADII) {
      sup_vec <- rep(NA_real_, n_w)
      for (i in seq_len(n_w)) {
        lo <- max(1L, i - rad)
        hi <- min(n_w, i + rad)
        neighbors <- is_elevated[lo:hi]
        sup_vec[i] <- mean(neighbors, na.rm = TRUE)
      }
      col_name <- paste0("nbhood_support_", rad)
      set(dt, j = col_name, value = round(sup_vec, 4))
    }

    # 3. Sim_mat block morphology
    SIM_HALF <- 10L
    block_compact   <- rep(NA_real_, n_w)
    block_coherence <- rep(NA_real_, n_w)
    block_frag      <- rep(NA_real_, n_w)
    square_support  <- rep(NA_real_, n_w)
    chr_sim_med <- median(sim_mat[upper.tri(sim_mat)], na.rm = TRUE)

    for (i in seq_len(n_w)) {
      lo <- max(1L, i - SIM_HALF)
      hi <- min(n_w, i + SIM_HALF)
      bsz <- hi - lo + 1L
      if (bsz < 5) next
      block <- sim_mat[lo:hi, lo:hi, drop = FALSE]

      block_vals <- block[upper.tri(block)]
      if (length(block_vals) >= 3) {
        block_compact[i] <- mean(block_vals, na.rm = TRUE)
      }

      diag_dist <- abs(row(block) - col(block))
      near_diag <- block[diag_dist <= 2 & diag_dist > 0]
      far_diag  <- block[diag_dist > max(2, bsz %/% 3)]
      if (length(near_diag) >= 2 && length(far_diag) >= 2) {
        mn_near <- mean(near_diag, na.rm = TRUE)
        mn_far  <- mean(far_diag, na.rm = TRUE)
        block_coherence[i] <- if (mn_near > 0.01) mn_far / mn_near else 0
      }

      if (length(block_vals) >= 3) {
        bm <- mean(block_vals, na.rm = TRUE)
        bsd <- sd(block_vals, na.rm = TRUE)
        block_frag[i] <- if (bm > 0.01) bsd / bm else 0
      }

      if (is.finite(chr_sim_med) && length(block_vals) >= 3) {
        square_support[i] <- mean(block_vals > chr_sim_med, na.rm = TRUE)
      }
    }

    dt[, local_block_compactness   := round(block_compact, 4)]
    dt[, local_block_coherence     := round(block_coherence, 4)]
    dt[, local_block_fragmentation := round(block_frag, 4)]
    dt[, local_square_support      := round(square_support, 4)]

    # 4. Composite morphology scores
    clamp01 <- function(x) { x[!is.finite(x)] <- 0; pmin(1, pmax(0, x)) }
    na0     <- function(x) { x[!is.finite(x)] <- 0; x }

    # 4a. flat_inv_score
    fl_run    <- clamp01(log2(pmax(run_len_vec, 1)) / 5)
    fl_z      <- clamp01(na0(run_mean_z) / 4)
    fl_stab   <- clamp01(1 - na0(run_sd_z) / 2)
    fl_jag    <- clamp01(1 - na0(jagged_vec) / 2)
    fl_block  <- clamp01(na0(block_compact))
    fl_coher  <- clamp01(na0(block_coherence))
    fl_nofrag <- clamp01(1 - na0(block_frag))
    flat_raw <- (0.22 * fl_run + 0.16 * fl_z + 0.12 * fl_stab +
                 0.10 * fl_jag + 0.16 * fl_block + 0.12 * fl_coher +
                 0.12 * fl_nofrag)
    flat_raw[run_len_vec < 3] <- 0
    dt[, flat_inv_score := round(flat_raw, 4)]

    # 4b. spiky_inv_score
    sp_prom   <- clamp01(na0(prom_vec) / 3)
    sp_sup5   <- clamp01(na0(dt$nbhood_support_5))
    sp_block  <- clamp01(na0(block_compact))
    sp_nofrag <- clamp01(1 - na0(block_frag) * 0.5)
    sp_short  <- clamp01(fifelse(run_len_vec >= 3 & run_len_vec <= 12, 1,
                         fifelse(run_len_vec > 12, 0.5, 0.3)))
    sp_antiflat <- clamp01(1 - na0(dt$plateau_flatness) * 0.5)
    spiky_raw <- (0.28 * sp_prom + 0.17 * sp_sup5 + 0.17 * sp_block +
                  0.12 * sp_nofrag + 0.12 * sp_short + 0.14 * sp_antiflat)
    spiky_raw[z_vec < z_thresh * 0.5] <- 0
    dt[, spiky_inv_score := round(spiky_raw, 4)]

    # 4c. fragmentation_score
    fg_jag     <- clamp01(na0(jagged_vec) / 2)
    fg_frag    <- clamp01(na0(block_frag))
    fg_nocoher <- clamp01(1 - na0(block_coherence))
    sup5  <- na0(dt$nbhood_support_5)
    sup20 <- na0(dt$nbhood_support_20)
    fg_instab    <- clamp01(pmax(0, sup5 - sup20))
    fg_shortrun  <- clamp01(fifelse(run_len_vec >= 1 & run_len_vec <= 2, 1,
                            fifelse(run_len_vec == 3, 0.5, 0)))
    frag_raw <- (0.25 * fg_jag + 0.20 * fg_frag + 0.15 * fg_nocoher +
                 0.20 * fg_instab + 0.20 * fg_shortrun)
    frag_raw[z_vec < 1.0] <- 0
    dt[, fragmentation_score := round(frag_raw, 4)]

    message("[PRECOMP] ", chr, ": morphology features computed (",
            "flat>0.3: ", sum(dt$flat_inv_score > 0.3, na.rm = TRUE),
            ", spiky>0.3: ", sum(dt$spiky_inv_score > 0.3, na.rm = TRUE),
            ", frag>0.3: ", sum(dt$fragmentation_score > 0.3, na.rm = TRUE), ")")
  } else {
    morph_cols <- c("local_jaggedness", "local_run_len", "local_run_mean_z",
                    "local_run_sd_z", "run_start_idx", "run_end_idx",
                    "candidate_length_windows", "candidate_length_bp",
                    "local_peak_prominence", "plateau_flatness",
                    "nbhood_support_5", "nbhood_support_10", "nbhood_support_20",
                    "local_block_compactness", "local_block_coherence",
                    "local_block_fragmentation", "local_square_support",
                    "flat_inv_score", "spiky_inv_score", "fragmentation_score")
    for (mc in morph_cols) set(dt, j = mc, value = NA_real_)
    message("[PRECOMP] ", chr, ": too few windows for morphology (", n_w, ")")
  }

  # ── Save precomp RDS ──
  precomp <- list(
    dt = dt, sim_mat = sim_mat, mds_mat = mds_mat,
    chrom = chr, n_windows = nrow(dt),
    bg_continuity_quantiles = bg_q
  )
  rds_out <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
  saveRDS(precomp, rds_out)

  # ── NN-smoothed sim_mats at multiple scales ──
  # MDS-space k-nearest-neighbor smoothing. nn_birth (the coarsest scale at
  # which a block first appears) is the persistence indicator used by D02/D09.
  sim_dir <- file.path(precomp_dir, "sim_mats")
  dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sim_mat, file.path(sim_dir, paste0(chr, ".sim_mat_nn0.rds")))
  message("[PRECOMP] ", chr, ": saved sim_mat_nn0")

  nn_sim_scales <- as.integer(strsplit(
    Sys.getenv("NN_SIM_SCALES",
               "20,40,80,120,160,200,240,320"), ",")[[1]])
  nn_sim_scales <- nn_sim_scales[nn_sim_scales > 0 & is.finite(nn_sim_scales)]

  mds_cols_nn <- grep("^MDS[0-9]+$", names(dt), value = TRUE)
  if (length(mds_cols_nn) > 0 && length(nn_sim_scales) > 0) {
    mds_mat_nn <- as.matrix(dt[, ..mds_cols_nn])
    n_mds <- nrow(mds_mat_nn)
    for (k in nn_sim_scales) {
      k_use <- min(k, n_mds - 1L)
      if (k_use < 2L) {
        message("[PRECOMP] ", chr, ": skipping nn", k,
                " (n_windows=", n_mds, " too small)")
        next
      }
      t_nn <- proc.time()
      smoothed <- matrix(0, nrow = n_mds, ncol = ncol(mds_mat_nn))
      for (wi in seq_len(n_mds)) {
        d <- dmat[wi, ]; d[wi] <- Inf
        nn_idx <- order(d)[seq_len(k_use)]
        smoothed[wi, ] <- colMeans(
          mds_mat_nn[c(wi, nn_idx), , drop = FALSE], na.rm = TRUE)
      }
      nn_dmat <- as.matrix(dist(smoothed))
      nn_sim <- make_sim_mat(nn_dmat)
      saveRDS(nn_sim,
              file.path(sim_dir, paste0(chr, ".sim_mat_nn", k, ".rds")))
      elapsed_nn <- round((proc.time() - t_nn)[3], 1)
      message("[PRECOMP] ", chr, ": saved sim_mat_nn", k,
              " (k_use=", k_use, ", ", elapsed_nn, "s)")
    }
  }

  elapsed <- round((proc.time() - t_chr)[3], 1)
  message("[PRECOMP] ", chr, ": ", nrow(dt), " windows (", elapsed, "s)")

  # ── Per-chrom QC stats ──
  win_bp <- if (all(c("start_bp", "end_bp") %in% names(dt))) dt$end_bp - dt$start_bp else integer(0)
  win_kb_mean    <- if (length(win_bp) > 0) round(mean(win_bp) / 1000, 2) else NA_real_
  win_kb_median  <- if (length(win_bp) > 0) round(stats::median(win_bp) / 1000, 2) else NA_real_
  chrom_span_mb  <- if (length(win_bp) > 0) round((max(dt$end_bp) - min(dt$start_bp)) / 1e6, 3) else NA_real_
  windows_per_mb <- if (is.finite(chrom_span_mb) && chrom_span_mb > 0) round(nrow(dt) / chrom_span_mb, 2) else NA_real_

  inv_like_mean   <- if ("inv_likeness" %in% names(dt)) round(mean(dt$inv_likeness, na.rm = TRUE), 4) else NA_real_
  inv_like_median <- if ("inv_likeness" %in% names(dt)) round(stats::median(dt$inv_likeness, na.rm = TRUE), 4) else NA_real_
  inv_like_q95    <- if ("inv_likeness" %in% names(dt)) round(stats::quantile(dt$inv_likeness, 0.95, na.rm = TRUE, names = FALSE), 4) else NA_real_
  n_inv_like_030  <- if ("inv_likeness" %in% names(dt)) sum(dt$inv_likeness >= 0.30, na.rm = TRUE) else 0L
  n_inv_like_070  <- if ("inv_likeness" %in% names(dt)) sum(dt$inv_likeness >= 0.70, na.rm = TRUE) else 0L

  n_bd_high       <- if ("band_discreteness" %in% names(dt)) sum(dt$band_discreteness > 2, na.rm = TRUE) else 0L
  bd_mean         <- if ("band_discreteness" %in% names(dt)) round(mean(dt$band_discreteness, na.rm = TRUE), 4) else NA_real_
  diff_mean       <- if ("diffuse_score" %in% names(dt))     round(mean(dt$diffuse_score, na.rm = TRUE), 4) else NA_real_
  n_three_cluster <- if ("n_effective_clusters" %in% names(dt)) sum(dt$n_effective_clusters == 3L, na.rm = TRUE) else 0L

  dos_cov_n     <- if ("dosage_het_rate_median" %in% names(dt)) sum(is.finite(dt$dosage_het_rate_median)) else 0L
  dos_cov_frac  <- if (nrow(dt) > 0) round(dos_cov_n / nrow(dt), 3) else NA_real_

  n_seed_eligible <- if ("max_abs_z" %in% names(dt) && "inv_likeness" %in% names(dt)) {
    sum((dt$max_abs_z >= 2.5) | (dt$inv_likeness >= 0.7), na.rm = TRUE)
  } else NA_integer_

  message(sprintf("[PRECOMP]   %s  span=%s Mb  win=%d  mean_kb=%s  inv_like>=0.5=%d  seed_elig=%s",
    chr,
    if (is.finite(chrom_span_mb)) sprintf("%.1f", chrom_span_mb) else "NA",
    nrow(dt),
    if (is.finite(win_kb_mean)) sprintf("%.1f", win_kb_mean) else "NA",
    sum(dt$inv_likeness >= 0.5, na.rm = TRUE),
    if (is.finite(n_seed_eligible)) as.character(n_seed_eligible) else "NA"
  ))

  data.table(
    chrom = chr, n_windows = nrow(dt),
    window_kb_mean   = win_kb_mean,
    window_kb_median = win_kb_median,
    chrom_span_mb    = chrom_span_mb,
    windows_per_mb   = windows_per_mb,
    n_inv_like_050   = sum(dt$inv_likeness >= 0.5, na.rm = TRUE),
    n_z_above_2 = sum(dt$max_abs_z >= 2.0, na.rm = TRUE),
    n_z_above_3 = sum(dt$max_abs_z >= 3.0, na.rm = TRUE),
    median_max_z = round(median(dt$max_abs_z, na.rm = TRUE), 3),
    q95_max_z = round(quantile(dt$max_abs_z, 0.95, na.rm = TRUE), 3),
    bg_q50 = round(bg_q["50%"], 4),
    bg_q85 = round(bg_q["85%"], 4),
    bg_q90 = round(bg_q["90%"], 4),
    bg_q95 = round(bg_q["95%"], 4),
    n_bg_pairs = length(bg_continuity),
    n_flat_030 = sum(dt$flat_inv_score > 0.3, na.rm = TRUE),
    n_spiky_030 = sum(dt$spiky_inv_score > 0.3, na.rm = TRUE),
    n_frag_030 = sum(dt$fragmentation_score > 0.3, na.rm = TRUE),
    inv_like_mean = inv_like_mean,
    inv_like_median = inv_like_median,
    inv_like_q95 = inv_like_q95,
    n_inv_like_030 = n_inv_like_030,
    n_inv_like_070 = n_inv_like_070,
    band_discreteness_mean = bd_mean,
    n_band_discrete_gt2 = n_bd_high,
    diffuse_score_mean = diff_mean,
    n_three_cluster = n_three_cluster,
    n_dosage_covered = dos_cov_n,
    dosage_coverage_frac = dos_cov_frac,
    n_seed_eligible = n_seed_eligible,
    elapsed_sec = elapsed
  )
}

message("[PRECOMP] Processing ", length(chroms), " chromosomes...")
t_all <- proc.time()

if (N_CORES > 1) {
  summary_list <- parallel::mclapply(chroms, precompute_one_chr, mc.cores = N_CORES)
} else {
  summary_list <- lapply(chroms, precompute_one_chr)
}

summary_list <- summary_list[!vapply(summary_list, is.null, logical(1))]
summary_dt <- rbindlist(summary_list)
elapsed_total <- round((proc.time() - t_all)[3], 1)
fwrite(summary_dt, file.path(outdir, "precomp_summary.tsv"), sep = "\t")

message("\n[DONE] Precompute complete (v10.0 SLIM — local PCA z-outlier path only)")
message("  Chromosomes: ", length(chroms), " (", N_CORES, " cores, ", elapsed_total, "s total)")
message("  Total windows: ", sum(summary_dt$n_windows))
message("  Windows with inv_likeness >= 0.5: ", sum(summary_dt$n_inv_like_050))
message("  Windows with robust |z| >= 2.0: ", sum(summary_dt$n_z_above_2))
message("  Windows with robust |z| >= 3.0: ", sum(summary_dt$n_z_above_3))
message("  flat_inv_score > 0.3:      ", sum(summary_dt$n_flat_030))
message("  spiky_inv_score > 0.3:     ", sum(summary_dt$n_spiky_030))
message("  fragmentation_score > 0.3: ", sum(summary_dt$n_frag_030))

if ("n_inv_like_030" %in% names(summary_dt)) {
  message("  -- Inv-likeness distribution --")
  message(sprintf("    >=0.30:  %d  (total)", sum(summary_dt$n_inv_like_030)))
  message(sprintf("    >=0.50:  %d", sum(summary_dt$n_inv_like_050)))
  message(sprintf("    >=0.70:  %d", sum(summary_dt$n_inv_like_070)))
}

if ("n_dosage_covered" %in% names(summary_dt)) {
  tot_dos <- sum(summary_dt$n_dosage_covered, na.rm = TRUE)
  tot_win <- sum(summary_dt$n_windows, na.rm = TRUE)
  if (tot_dos > 0) {
    message(sprintf("  -- Dosage het coverage: %d / %d windows (%.1f%%)",
      tot_dos, tot_win, 100 * tot_dos / tot_win))
  }
}

message("\n  Precomp dir:  ", precomp_dir)
message("  Window TSV:   ", f_inv)
message("  Summary:      ", file.path(outdir, "precomp_summary.tsv"))
message("\n  Next: ")
message("    - JSON per chrom:   Rscript STEP_C01a_export_json.R <outdir> <json_outdir>")
message("    - Seeded regions:   sbatch LAUNCH_STEP_C01b_1_seeded_regions.slurm")
message("    - Block detection:  sbatch LAUNCH_PHASE_01C_block_detect.slurm")
