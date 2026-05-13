#!/usr/bin/env Rscript

# =============================================================================
# STEP_ZO_G_precompute.R
#
# Per-chromosome precomputation for the Z (local-PCA dosage) path. Reads
# the unified MDS RDS produced by ZO_F and emits ONE precomp.rds per
# chromosome plus the NN-smoothed similarity matrices that the L1/L2
# detectors and the atlas consume. Parallelized across chromosomes via
# mclapply.
#
# Pipeline position:
#   ZO_E (per-chrom MDS) -> ZO_F (merge to unified MDS RDS) ->
#   ZO_G (this script: precomp + sim_mats) ->
#   ZO_H / ZO_J (L1 / L2 stripe detection on sim_mat) ->
#   atlas (per-region sample PCA on demand)
#
# Inputs:
#   <step02b_outprefix>     positional 1: prefix from ZO_F
#                           (the script appends .mds.rds and reads from there)
#   <outdir>                positional 2: output root
#   [--dosage_dir <dir>]    NOT CURRENTLY CONSUMED (legacy flag retained for
#                           CLI compat). The dosage-het-rate pass that used
#                           it was removed 2026-05-13.
#
# Outputs (in <outdir>/):
#   precomp/<chr>.precomp.rds                 per-chrom precomp
#       $dt       per-window data.table:
#           position : global_window_id, chrom, start_bp, end_bp, mid_bp
#           MDS      : MDS1..MDSk, MDS{1..k}_z, max_abs_z, max_z_axis
#           seed     : seed_nn_dist
#           **per-sample**: PC_1_Ind*, PC_2_Ind* (atlas per-region PCA input)
#       $sim_mat                               window x window similarity (raw)
#       $mds_mat, $chrom, $n_windows
#
#   precomp/sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
#                              MDS-space k-NN smoothed similarity matrices
#                              Tune via env var:
#                                NN_SIM_SCALES="40,80,160,320" sbatch ...
#
#   precomp_summary.tsv        per-chrom window counts + timing
#
# Removed 2026-05-13:
#   - genome-wide inv_likeness pass (band_discreteness, diffuse_score,
#     pc1_bimodality, het_pc1_gap, n_effective_clusters, etc.)
#   - dosage_het_rate_* pass
#   - local-k3 features (local_delta12, local_entropy, local_ena)
#   - beta_adaptive p-values (adaptive_seed, beta_pval, beta_alpha, beta_beta)
#   - bg_continuity baseline
#   - morphology features (flat_inv_score, spiky_inv_score,
#     fragmentation_score, local_jaggedness, nbhood_support_*, etc.)
#   - window_dt.tsv.gz genome-wide table
# All were lostruct-style feature-engineering with zero downstream consumers
# in this codebase (L1/L2 read sim_mat; atlas does per-region sample PCA).
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
  n_w <- nrow(dt)

  # Seed NN distances
  nn_dists <- vapply(seq_len(nrow(dmat)), function(i) {
    d <- dmat[i, ]; d[i] <- Inf; d <- d[is.finite(d)]
    if (length(d) == 0) return(Inf)
    mean(sort(d)[seq_len(min(SEED_NEIGHBOR_K, length(d)))], na.rm = TRUE)
  }, numeric(1))
  dmax_nn <- quantile(dmat[is.finite(dmat)], 0.95, na.rm = TRUE)
  if (!is.finite(dmax_nn) || dmax_nn == 0) dmax_nn <- 1
  dt[, seed_nn_dist := nn_dists / dmax_nn]

  # ── Save precomp RDS ──
  precomp <- list(
    dt = dt, sim_mat = sim_mat, mds_mat = mds_mat,
    chrom = chr, n_windows = nrow(dt)
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

  message(sprintf("[PRECOMP]   %s  span=%s Mb  win=%d  mean_kb=%s",
    chr,
    if (is.finite(chrom_span_mb)) sprintf("%.1f", chrom_span_mb) else "NA",
    nrow(dt),
    if (is.finite(win_kb_mean)) sprintf("%.1f", win_kb_mean) else "NA"
  ))

  data.table(
    chrom = chr, n_windows = nrow(dt),
    window_kb_mean   = win_kb_mean,
    window_kb_median = win_kb_median,
    chrom_span_mb    = chrom_span_mb,
    windows_per_mb   = windows_per_mb,
    n_z_above_2 = sum(dt$max_abs_z >= 2.0, na.rm = TRUE),
    n_z_above_3 = sum(dt$max_abs_z >= 3.0, na.rm = TRUE),
    median_max_z = round(median(dt$max_abs_z, na.rm = TRUE), 3),
    q95_max_z = round(quantile(dt$max_abs_z, 0.95, na.rm = TRUE), 3),
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

message("\n[DONE] ZO_G precompute complete")
message("  Chromosomes: ", length(chroms), " (", N_CORES, " cores, ", elapsed_total, "s total)")
message("  Total windows: ", sum(summary_dt$n_windows))
message("  Windows with robust |z| >= 2.0: ", sum(summary_dt$n_z_above_2))
message("  Windows with robust |z| >= 3.0: ", sum(summary_dt$n_z_above_3))

message("\n  Precomp dir:  ", precomp_dir)
message("  Summary:      ", file.path(outdir, "precomp_summary.tsv"))
