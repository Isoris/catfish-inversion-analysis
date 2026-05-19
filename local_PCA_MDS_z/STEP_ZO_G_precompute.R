#!/usr/bin/env Rscript

# =============================================================================
# STEP_ZO_G_precompute.R  (v10, 2026-05-16, harmonized-layout cleanup)
#
# Per-chromosome precompute for the Z (local-PCA dosage) path. Reads the
# per-chrom MDS results from <mds_dir>/<chr>.mds_perchr.rds (written by
# ZO_E into the harmonized 03_mds/ slot, formerly <outdir>/tmp/) and emits
# ONE precomp.rds per chromosome plus the NN-smoothed similarity matrices
# that the L1/L2 detectors and the atlas consume. Parallelized across
# chromosomes via mclapply.
#
# v10 changes:
#   - Reads per-chrom MDS RDS from <mds_dir> (env PATH1_MDS / --mds_dir)
#     instead of the legacy <outdir>/tmp/ subdir.
#   - Writes per-chrom precomp directly into <outdir> (the harmonized
#     04_precomp/ slot) — no more <outdir>/precomp/ nesting.
#   - Does NOT emit a whole-genome combined RDS. Per-chrom only.
#
# Pipeline position (no merge step):
#   ZO_E (per-chrom MDS) -> ZO_G (this script: precomp + sim_mats) ->
#   ZO_H / ZO_J (L1 / L2 stripe detection on sim_mat) ->
#   atlas (per-region sample PCA on demand)
#
# Inputs:
#   <mds_prefix>            positional 1; legacy; only dirname() is used as
#                           a fallback to locate the per-chrom MDS RDS dir
#                           if --mds_dir is not given.
#   <outdir>                positional 2: precomp output root (04_precomp/)
#   [--mds_dir <dir>]       per-chrom MDS RDS dir from ZO_E (03_mds/).
#                           Defaults to env PATH1_MDS, then dirname(<mds_prefix>).
#   [--dosage_dir <dir>]    legacy flag; silently ignored (the dosage-het-rate
#                           pass was removed 2026-05-13).
#
# Outputs (in <outdir>/, the harmonized 04_precomp/ slot):
#   <chr>.precomp.rds                         per-chrom precomp
#       $dt       per-window data.table:
#           position : global_window_id, chrom, start_bp, end_bp, mid_bp
#           MDS      : MDS1..MDSk, MDS{1..k}_z, max_abs_z, max_z_axis
#           seed     : seed_nn_dist
#           **per-sample**: PC_1_Ind*, PC_2_Ind* (atlas per-region PCA input)
#       $sim_mat                               window x window similarity (raw)
#       $mds_mat, $chrom, $n_windows
#
#   sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
#                              MDS-space k-NN smoothed similarity matrices
#                              Tune via env var:
#                                NN_SIM_SCALES="40,80,160,320" sbatch ...
#
#   precomp_summary.tsv        per-chrom window counts + timing (TSV only;
#                              there is NO whole-genome combined RDS)
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
  stop("Usage: Rscript STEP_ZO_G_precompute.R <mds_prefix> <outdir>",
       " [--mds_dir <dir>]\n",
       "  <mds_prefix>  legacy positional; dirname() used as MDS dir fallback\n",
       "  <outdir>      precomp output root (harmonized 04_precomp slot)\n",
       "  --mds_dir     dir of <chr>.mds_perchr.rds from ZO_E (03_mds slot).\n",
       "                Defaults to env PATH1_MDS, then dirname(<mds_prefix>).")
}

step10_prefix <- args[1]
outdir        <- args[2]
# v10: precomp lands DIRECTLY in <outdir> (no nested precomp/). sim_mats
# remain in <outdir>/sim_mats/.
precomp_dir   <- outdir
dir.create(precomp_dir, recursive = TRUE, showWarnings = FALSE)

# Optional --mds_dir override; legacy --dosage_dir silently ignored.
mds_dir_arg <- NA_character_
i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--mds_dir"    && i < length(args)) { mds_dir_arg <- args[i + 1L]; i <- i + 2L }
  else if (a == "--dosage_dir" && i < length(args)) { i <- i + 2L }
  else { i <- i + 1L }
}

# =============================================================================
# PARAMETERS
# =============================================================================

SEED_MDS_AXES    <- 5L    # number of MDS axes contributing to max_abs_z
SEED_NEIGHBOR_K  <- 3L    # k for seed_nn_dist

# =============================================================================
# LOAD PER-CHROM MDS RESULTS
# =============================================================================
# v10: ZO_E now writes per-chrom <chr>.mds_perchr.rds directly into the
# harmonized 03_mds/ slot (was <outdir>/tmp/). Resolution order:
#   1. --mds_dir <dir>            explicit CLI override
#   2. $PATH1_MDS                 from harmonized config
#   3. dirname(<mds_prefix>)      fallback for legacy callers
#   4. <outdir>/tmp/              last-resort back-compat with v9 layouts

resolve_mds_dir <- function() {
  candidates <- c(
    if (!is.na(mds_dir_arg))             mds_dir_arg                      else NULL,
    if (nzchar(Sys.getenv("PATH1_MDS"))) Sys.getenv("PATH1_MDS")         else NULL,
    dirname(step10_prefix),
    file.path(outdir, "tmp")
  )
  for (cand in candidates) if (dir.exists(cand)) return(cand)
  stop("Could not locate per-chrom MDS dir. Tried: ",
       paste(candidates, collapse = " | "))
}
mds_dir <- resolve_mds_dir()

perchr_files <- sort(list.files(mds_dir, pattern = "\\.mds_perchr\\.rds$",
                                full.names = TRUE))
if (length(perchr_files) == 0L) stop("No .mds_perchr.rds files in ", mds_dir)

message("[PRECOMP] Loading ", length(perchr_files), " per-chrom MDS results from ", mds_dir)
t_load <- proc.time()
per_chr <- list()
for (f in perchr_files) {
  obj <- readRDS(f)
  if (!is.null(obj$skip) && obj$skip) {
    message("[PRECOMP] SKIP ", obj$chrom, ": ", obj$reason)
    next
  }
  chr <- obj$out_dt$chrom[1]
  per_chr[[chr]] <- list(out_dt = obj$out_dt, dmat = obj$dmat, mds = obj$mds)
}
message("[PRECOMP] Loaded in ", round((proc.time() - t_load)[3], 1), "s")

if (length(per_chr) == 0L) stop("No usable per-chrom MDS results")
chroms <- names(per_chr)
message("[PRECOMP] ", length(chroms), " chromosomes")

# Extract sample names from any chrom's per-sample PC columns
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

  # Vectorized max_abs_z + max_z_axis (much faster than apply(.SD, 1, ...))
  z_cols <- grep("^MDS[0-9]+_z$", names(dt), value = TRUE)
  if (length(z_cols) > SEED_MDS_AXES) z_cols <- z_cols[seq_len(SEED_MDS_AXES)]
  if (length(z_cols) > 0) {
    zmat <- abs(as.matrix(dt[, ..z_cols]))
    zmat[!is.finite(zmat)] <- 0
    dt[, max_abs_z := do.call(pmax, c(as.data.frame(zmat), list(na.rm = TRUE)))]
    dt[, max_z_axis := max.col(zmat, ties.method = "first")]
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

  # ── One-time per-chrom precomputations ─────────────────────────────────
  # Cache the 95th-percentile distance (used by the similarity transform) and
  # the per-window sort order of dmat (used by seed_nn_dist + every NN scale,
  # was previously re-sorted 8x).
  finite_d <- dmat[is.finite(dmat)]
  dmax <- if (length(finite_d) > 0L) quantile(finite_d, 0.95, na.rm = TRUE) else 1
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  rm(finite_d)

  n_w <- nrow(dmat)
  # sorted_idx[i, ] = indices of windows ranked by ascending dmat[i, ].
  # Computed once; reused for seed_nn_dist and every NN-smoothing scale.
  sorted_idx <- matrix(NA_integer_, nrow = n_w, ncol = n_w)
  for (i in seq_len(n_w)) {
    d <- dmat[i, ]; d[i] <- Inf
    sorted_idx[i, ] <- order(d)
  }

  # Base similarity matrix (uses cached dmax).
  sim_mat <- 1 - pmin(dmat / dmax, 1)
  sim_mat[!is.finite(sim_mat)] <- 0
  diag(sim_mat) <- 1

  # Seed NN distances (use cached sorted_idx; mean of top SEED_NEIGHBOR_K).
  K_SEED <- min(SEED_NEIGHBOR_K, n_w - 1L)
  nn_dists <- vapply(seq_len(n_w), function(i) {
    mean(dmat[i, sorted_idx[i, seq_len(K_SEED)]], na.rm = TRUE)
  }, numeric(1))
  dt[, seed_nn_dist := nn_dists / dmax]

  # ── Save precomp RDS ──
  precomp <- list(
    dt = dt, sim_mat = sim_mat, mds_mat = mds_mat,
    chrom = chr, n_windows = nrow(dt)
  )
  rds_out <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
  saveRDS(precomp, rds_out)

  # ── NN-smoothed sim_mats at multiple scales ────────────────────────────
  # nn0 is the unsmoothed base sim_mat (same as $sim_mat inside precomp.rds,
  # written here under the conventional name). nn{40,80,160,320} are the
  # scales L1 / L2 / atlas consume. Earlier defaults also wrote nn{20,120,
  # 200,240} which had zero readers. Override via NN_SIM_SCALES env var.
  sim_dir <- file.path(precomp_dir, "sim_mats")
  dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sim_mat, file.path(sim_dir, paste0(chr, ".sim_mat_nn0.rds")))
  message("[PRECOMP] ", chr, ": saved sim_mat_nn0")

  nn_sim_scales <- as.integer(strsplit(
    Sys.getenv("NN_SIM_SCALES", "40,80,160,320"), ",")[[1]])
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

      # Smooth MDS coords by averaging each window with its top-k neighbors
      # (neighbors picked by cached sorted_idx; no per-scale resort).
      smoothed <- matrix(0, nrow = n_mds, ncol = ncol(mds_mat_nn))
      for (wi in seq_len(n_mds)) {
        nn_idx <- sorted_idx[wi, seq_len(k_use)]
        smoothed[wi, ] <- colMeans(
          mds_mat_nn[c(wi, nn_idx), , drop = FALSE], na.rm = TRUE)
      }

      # Pairwise Euclidean distance via tcrossprod (BLAS) instead of
      # dist() + as.matrix(): D^2[i,j] = ||s_i||^2 + ||s_j||^2 - 2*<s_i,s_j>.
      # Much faster for moderate N (~5-10x for N=6000, K_MDS=5).
      sq_norms <- rowSums(smoothed * smoothed)
      nn_dsq <- outer(sq_norms, sq_norms, "+") - 2 * tcrossprod(smoothed)
      nn_dsq[nn_dsq < 0] <- 0
      nn_dmat <- sqrt(nn_dsq)

      # Similarity transform: linear rescaling against this scale's own 95th
      # percentile (matches the upstream sim_mat formula).
      nn_dmax <- quantile(nn_dmat[is.finite(nn_dmat) & nn_dmat > 0], 0.95,
                          na.rm = TRUE)
      if (!is.finite(nn_dmax) || nn_dmax == 0) nn_dmax <- 1
      nn_sim <- 1 - pmin(nn_dmat / nn_dmax, 1)
      nn_sim[!is.finite(nn_sim)] <- 0
      diag(nn_sim) <- 1

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
