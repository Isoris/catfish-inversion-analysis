#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_G_atlas_json.R   (drafts/v5)
# =============================================================================
# Atlas JSON exporter for the theta-pi local-PCA / MDS / CUSUM pipeline.
# Combines outputs from TR_B (precomp), TR_D (L2 envelopes), TR_E (carrier
# bands), TR_F (CUSUM) into one consolidated per-chrom JSON the page-12
# atlas reads.
#
# All emitted window indices are 0-INDEXED (atlas is JS). R-internal
# computation stays 1-indexed; only the JSON values are converted.
#
# Reads:  <OUTROOT>/03_per_chrom/<chr>/precomp.rds
#         <OUTROOT>/03_per_chrom/<chr>/L1_envelopes.tsv
#         <OUTROOT>/03_per_chrom/<chr>/L2_envelopes.tsv
#         <OUTROOT>/03_per_chrom/<chr>/carrier_assignments.tsv
#         <OUTROOT>/03_per_chrom/<chr>/cusum_per_sample.tsv.gz
#         <OUTROOT>/03_per_chrom/<chr>/cusum_boundary_dist.tsv
# Writes: <OUTROOT>/04_atlas_json/<chr>/<chr>_phase2_theta.json
# Usage:  Rscript STEP_TR_G_atlas_json.R --chrom <CHR>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
CHROM <- NULL
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom" && i < length(args)) { CHROM <- args[i + 1]; i <- i + 2L }
  else { i <- i + 1L }
}

OUTROOT      <- Sys.getenv("OUTROOT",      unset = NA)
PESTPG_SCALE <- Sys.getenv("PESTPG_SCALE", unset = "win10000.step2000")
JSON_OUT_DIR <- Sys.getenv("JSON_OUT_DIR", unset = file.path(OUTROOT, "04_atlas_json"))
stopifnot(!is.na(OUTROOT))

per_chrom_dir <- file.path(OUTROOT, "03_per_chrom")
chroms <- if (!is.null(CHROM)) CHROM else list.files(per_chrom_dir, pattern = "^C_gar_LG[0-9]+$")

clean <- function(x, digits = 6) {
  out <- round(as.numeric(x), digits); out[!is.finite(out)] <- NA_real_; out
}

for (chrom in chroms) {
  rds  <- file.path(per_chrom_dir, chrom, "precomp.rds")
  l1f  <- file.path(per_chrom_dir, chrom, "L1_envelopes.tsv")
  l2f  <- file.path(per_chrom_dir, chrom, "L2_envelopes.tsv")
  caf  <- file.path(per_chrom_dir, chrom, "carrier_assignments.tsv")
  csf  <- file.path(per_chrom_dir, chrom, "cusum_per_sample.tsv.gz")
  bdf  <- file.path(per_chrom_dir, chrom, "cusum_boundary_dist.tsv")
  if (!file.exists(rds)) { message("[TR_G] ", chrom, ": no precomp — skip"); next }

  precomp <- readRDS(rds)
  dt <- precomp$dt
  n_win  <- precomp$n_windows
  n_samp <- precomp$n_samples
  sample_order <- precomp$sample_order

  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  pc2_cols <- grep("^PC_2_", names(dt), value = TRUE)
  pc1_mat  <- as.matrix(dt[, ..pc1_cols])
  pc2_mat  <- if (length(pc2_cols) > 0) as.matrix(dt[, ..pc2_cols]) else matrix(NA_real_, n_win, 0)

  # --- theta_pi_per_window: rebuild from per_sample × per_window. We keep
  # per-window medians + per-window n_sites_floor on the dt; full per-sample
  # θπ values aren't in the precomp (kept off the .rds to save memory). If
  # this layer needs full per-sample θπ in the JSON, reload theta_native TSV.
  THETA_TSV_DIR <- Sys.getenv("THETA_TSV_DIR", unset = NA)
  long_dt <- if (!is.na(THETA_TSV_DIR)) {
    tsv <- file.path(THETA_TSV_DIR,
                     sprintf("theta_native.%s.%s.tsv.gz", chrom, PESTPG_SCALE))
    if (file.exists(tsv)) fread(tsv)[chrom == ..chrom] else NULL
  } else NULL

  per_window_block <- NULL
  if (!is.null(long_dt) && nrow(long_dt) > 0) {
    samp_to_row <- setNames(seq_along(sample_order), sample_order)
    win_to_col  <- setNames(seq_along(dt$window_idx), as.character(dt$window_idx))
    theta_mat   <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                          dimnames = list(sample_order, NULL))
    rows <- samp_to_row[long_dt$sample]
    cols <- win_to_col[as.character(long_dt$window_idx)]
    good <- !is.na(rows) & !is.na(cols)
    theta_mat[cbind(rows[good], cols[good])] <- long_dt$theta_pi[good]
    samples_block <- lapply(seq_len(n_samp), function(s) list(
      sample_id = sample_order[s],
      theta_pi  = clean(theta_mat[s, ], 6)
    ))
    values_flat <- as.numeric(t(theta_mat))   # row-major over windows
    per_window_block <- list(
      schema_version = 2L, layer = "theta_pi_per_window", chrom = chrom,
      scale = PESTPG_SCALE,
      n_samples = n_samp, n_windows = n_win,
      sample_ids = sample_order,
      values = clean(values_flat, 6),
      windows = lapply(seq_len(n_win), function(wi) list(
        idx = as.integer(dt$window_idx[wi]),
        start_bp = as.integer(dt$start_bp[wi]),
        end_bp   = as.integer(dt$end_bp[wi])
      )),
      samples = samples_block
    )
  }

  # --- theta_pi_local_pca (sim_mat + MDS + sign-aligned PC loadings) ---
  # Encode sim_mat
  if (precomp$sim_mat_format == "upper_triangle_float32") {
    M <- precomp$sim_mat
    ut_idx <- which(upper.tri(M, diag = TRUE), arr.ind = FALSE)
    sim_payload <- clean(M[ut_idx], 4)
    sim_band_half <- NULL
  } else {
    sim_payload <- clean(as.vector(t(precomp$sim_band)), 4)
    sim_band_half <- as.integer(precomp$sim_band_half)
  }
  local_pca_block <- list(
    schema_version       = 2L,
    layer                = "theta_pi_local_pca",
    chrom                = chrom,
    scale                = PESTPG_SCALE,
    n_samples            = n_samp,
    n_windows            = n_win,
    sample_order         = sample_order,
    pc1_loadings_aligned = lapply(seq_len(n_win), function(wi) clean(pc1_mat[wi, ], 6)),
    pc2_loadings_aligned = lapply(seq_len(n_win), function(wi) clean(pc2_mat[wi, ], 6)),
    lambda_1             = clean(dt$lambda_1, 6),
    lambda_2             = clean(dt$lambda_2, 6),
    lambda_ratio         = clean(dt$lambda_ratio, 4),
    z                    = clean(dt$max_abs_z, 4),
    z_top10_mean         = clean(dt$top10_abs_z, 4),
    mds_coords           = list(mds1 = clean(dt$MDS1, 6),
                                mds2 = clean(dt$MDS2, 6)),
    anchor_window_idx    = as.integer(dt$anchor_window_idx[1]),  # 0-indexed (set by TR_B)
    sim_mat_format       = precomp$sim_mat_format,
    sim_mat_band_half    = sim_band_half,
    sim_mat_n            = n_win,
    sim_mat              = sim_payload,
    `_qc`                = list(
      n_unflipped       = length(precomp$unflipped_windows),
      unflipped_windows = if (length(precomp$unflipped_windows) > 0)
                            as.integer(precomp$unflipped_windows) else integer(0)
    )
  )

  # --- theta_pi_envelopes (L1 + L2) — both 0-indexed in win_*_idx0 ---
  l1_dt <- if (file.exists(l1f)) fread(l1f) else data.table()
  l2_dt <- if (file.exists(l2f)) fread(l2f) else data.table()
  envs_block <- list(
    schema_version = 2L, layer = "theta_pi_envelopes", chrom = chrom,
    l1 = if (nrow(l1_dt) > 0) lapply(seq_len(nrow(l1_dt)), function(k) list(
      l1_id     = l1_dt$l1_id[k],
      win_start = as.integer(l1_dt$win_start_idx0[k]),
      win_end   = as.integer(l1_dt$win_end_idx0[k]),
      start_bp  = as.integer(l1_dt$start_bp[k]),
      end_bp    = as.integer(l1_dt$end_bp[k]),
      span_kb   = round(l1_dt$span_kb[k], 1),
      n_windows = as.integer(l1_dt$n_windows[k]),
      peak_z    = round(l1_dt$peak_z[k], 4),
      mean_z    = round(l1_dt$mean_z[k], 4)
    )) else list(),
    l2 = if (nrow(l2_dt) > 0) lapply(seq_len(nrow(l2_dt)), function(k) list(
      l2_id        = l2_dt$l2_id[k],
      candidate_id = l2_dt$candidate_id[k],
      l1_id        = l2_dt$l1_id[k],
      win_start    = as.integer(l2_dt$win_start_idx0[k]),
      win_end      = as.integer(l2_dt$win_end_idx0[k]),
      start_bp     = as.integer(l2_dt$start_bp[k]),
      end_bp       = as.integer(l2_dt$end_bp[k]),
      span_kb      = round(l2_dt$span_kb[k], 1),
      n_windows    = as.integer(l2_dt$n_windows[k]),
      peak_z       = round(l2_dt$peak_z[k], 4),
      mean_z       = round(l2_dt$mean_z[k], 4)
    )) else list()
  )

  # --- theta_pi_cusum (per-carrier breakpoints + boundary distribution) ---
  ca_dt <- if (file.exists(caf)) fread(caf) else data.table()
  cs_dt <- if (file.exists(csf)) fread(csf) else data.table()
  bd_dt <- if (file.exists(bdf)) fread(bdf) else data.table()

  cusum_candidates <- list()
  if (nrow(l2_dt) > 0L) {
    for (k in seq_len(nrow(l2_dt))) {
      cid <- l2_dt$candidate_id[k]
      ca_sub <- ca_dt[candidate_id == cid]
      cs_sub <- if (nrow(cs_dt) > 0) cs_dt[candidate_id == cid] else data.table()
      bd_sub <- if (nrow(bd_dt) > 0) bd_dt[candidate_id == cid] else data.table()

      bands <- list()
      for (b in unique(ca_sub$band)) {
        members <- ca_sub[band == b]
        cs_b    <- cs_sub[band == b]
        bd_b    <- bd_sub[band == b]
        bands[[b]] <- list(
          band            = b,
          n_members       = nrow(members),
          member_samples  = members$sample_id,
          per_sample      = lapply(seq_len(nrow(cs_b)), function(j) list(
            sample_id        = cs_b$sample_id[j],
            cp_bp            = as.integer(cs_b$cp_bp[j] %||% NA_integer_),
            strength         = round(cs_b$strength[j], 3),
            asymmetry        = as.integer(cs_b$asymmetry[j] %||% NA_integer_),
            informative      = as.logical(cs_b$informative[j]),
            cp_side_inferred = cs_b$cp_side_inferred[j]
          )),
          boundary_5_prime = if (nrow(bd_b[side == "5_prime"]) > 0) {
            r <- bd_b[side == "5_prime"]
            list(n_carriers = as.integer(r$n_carriers), median_bp = as.integer(r$median_bp),
                 iqr_kb = round(r$iqr_kb, 2), spread_class = r$spread_class,
                 peak_strength = round(r$peak_strength, 3))
          } else list(n_carriers = 0L),
          boundary_3_prime = if (nrow(bd_b[side == "3_prime"]) > 0) {
            r <- bd_b[side == "3_prime"]
            list(n_carriers = as.integer(r$n_carriers), median_bp = as.integer(r$median_bp),
                 iqr_kb = round(r$iqr_kb, 2), spread_class = r$spread_class,
                 peak_strength = round(r$peak_strength, 3))
          } else list(n_carriers = 0L)
        )
      }
      cusum_candidates[[length(cusum_candidates) + 1L]] <- list(
        candidate_id = cid,
        chrom        = chrom,
        l2_id        = l2_dt$l2_id[k],
        start_bp     = as.integer(l2_dt$start_bp[k]),
        end_bp       = as.integer(l2_dt$end_bp[k]),
        bands        = bands
      )
    }
  }

  cusum_block <- list(
    schema_version = 2L, layer = "theta_pi_cusum", chrom = chrom,
    n_candidates = length(cusum_candidates),
    candidates   = cusum_candidates
  )

  # --- tracks (per-window aggregates) ---
  tracks_block <- list(
    theta_pi_median = list(
      values = clean(dt$theta_pi_median, 6),
      pos_bp = as.integer(dt$mid_bp)
    ),
    theta_pi_z = list(
      values = clean(dt$max_abs_z, 4),
      pos_bp = as.integer(dt$mid_bp)
    ),
    theta_pi_lambda_ratio = list(
      values = clean(dt$lambda_ratio, 4),
      pos_bp = as.integer(dt$mid_bp)
    ),
    theta_pi_mds1 = list(values = clean(dt$MDS1, 6), pos_bp = as.integer(dt$mid_bp)),
    theta_pi_mds2 = list(values = clean(dt$MDS2, 6), pos_bp = as.integer(dt$mid_bp))
  )

  # --- top-level JSON ---
  layers <- c("theta_pi_local_pca", "theta_pi_envelopes", "theta_pi_cusum", "tracks")
  obj <- list(
    schema_version = 2L,
    chrom          = chrom,
    n_samples      = n_samp,
    n_windows      = n_win,
    scale          = PESTPG_SCALE,
    `_generated_at` = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    `_generator`    = "STEP_TR_G_atlas_json.R (v5)",
    `_layers_present` = layers,
    tracks                = tracks_block,
    theta_pi_local_pca    = local_pca_block,
    theta_pi_envelopes    = envs_block,
    theta_pi_cusum        = cusum_block
  )
  if (!is.null(per_window_block)) {
    obj$theta_pi_per_window <- per_window_block
    obj$`_layers_present` <- c("theta_pi_per_window", obj$`_layers_present`)
  }

  out_dir <- file.path(JSON_OUT_DIR, chrom)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_json <- file.path(out_dir, sprintf("%s_phase2_theta.json", chrom))
  write_json(obj, out_json, auto_unbox = TRUE, na = "null", pretty = FALSE, digits = NA)
  fi <- file.info(out_json)
  message(sprintf("[TR_G] %s: wrote %s (%.2f MB)",
                  chrom, out_json, fi$size / 1024 / 1024))
}

message("[TR_G] DONE")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a