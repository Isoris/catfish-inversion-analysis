#!/usr/bin/env Rscript
# =============================================================================
# STEP_GH_J_atlas_json.R
# =============================================================================
# Atlas JSON exporter for the GHSL pipeline. Mirrors the structure of TR_J
# (theta_pi) and ZO_M (z-blocks) but trimmed for GHSL:
#   - GHSL has no dense per-sample feature track (yet), so no per_window block
#     and no carrier / CUSUM block.
#   - Output is single-scale: just the coarse GHSL precomp produced by GH_C
#     at GHSL_PCA_SCALE_FULL (default s50), plus L1/L2 envelopes/boundaries.
#
# Layers emitted:
#   ghsl_local_pca   ← sim_mat thumbnail + MDS coords + per-window PC loadings
#                      + |Z| profile.
#   ghsl_envelopes   ← L1 + L2 envelopes + boundary records (from TR_D / TR_F
#                      run against the GHSL precomp).
#   tracks           ← per-window summary aggregates (div_median, max_abs_z,
#                      lambda_ratio, MDS1, MDS2, ...).
#
# All emitted window indices are 0-INDEXED (atlas is JS).
#
# CLI:
#   --chr <chrom>
#   --precomp_dir <GHSL_PRECOMP_DIR>
#   --l1_dir      <GHSL_L1_DIR>
#   --l2_dir      <GHSL_L2_DIR>
#   --outdir      <GHSL_JSON_DIR>    (per-chrom subdir created)
#   --scale       <s50>              (cosmetic label in JSON; default $GHSL_PCA_SCALE_FULL)
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

get_arg <- function(flag, default = NA_character_) {
  args <- commandArgs(trailingOnly = TRUE); i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}

chrom        <- get_arg("--chr")
precomp_dir  <- get_arg("--precomp_dir")
l1_dir       <- get_arg("--l1_dir")
l2_dir       <- get_arg("--l2_dir")
outdir       <- get_arg("--outdir")
scale_label  <- get_arg("--scale", Sys.getenv("GHSL_PCA_SCALE_FULL", unset = "s50"))

stopifnot(!is.na(chrom), !is.na(precomp_dir), !is.na(outdir))
if (is.na(l1_dir)) l1_dir <- file.path(dirname(precomp_dir), "04_L1_detect")
if (is.na(l2_dir)) l2_dir <- file.path(dirname(precomp_dir), "06_L2_detect")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a))) b else a

clean <- function(x, digits = 6) {
  out <- round(as.numeric(x), digits); out[!is.finite(out)] <- NA_real_; out
}
to_idx0 <- function(x) if (length(x) == 0L) integer(0) else as.integer(x) - 1L

precomp_rds <- file.path(precomp_dir, paste0(chrom, ".ghsl_precomp.rds"))
if (!file.exists(precomp_rds))
  stop(sprintf("[GH_J] %s: missing precomp at %s", chrom, precomp_rds))

precomp <- readRDS(precomp_rds); dt <- precomp$dt
n_win   <- as.integer(precomp$n_windows)
n_samp  <- as.integer(precomp$n_samples)
sample_order <- precomp$sample_order
npc     <- as.integer(precomp$npc %||% 2L)

# ── ghsl_local_pca: sim_mat + MDS + per-window PC loadings ────────────────────
pc_blocks <- vector("list", npc)
for (k in seq_len(npc)) {
  cols <- grep(paste0("^PC_", k, "_"), names(dt), value = TRUE)
  if (length(cols) == 0L) next
  pc_mat <- as.matrix(dt[, ..cols])
  pc_blocks[[k]] <- lapply(seq_len(n_win), function(wi) clean(pc_mat[wi, ], 6))
}

if (precomp$sim_mat_format == "upper_triangle_float32") {
  M <- precomp$sim_mat
  ut_idx <- which(upper.tri(M, diag = TRUE), arr.ind = FALSE)
  sim_payload   <- clean(M[ut_idx], 4)
  sim_band_half <- NULL
} else {
  sim_payload   <- clean(as.vector(t(precomp$sim_band)), 4)
  sim_band_half <- as.integer(precomp$sim_band_half)
}

mds_cols <- grep("^MDS[0-9]+$", names(dt), value = TRUE)
mds_coords <- setNames(
  lapply(mds_cols, function(c) clean(dt[[c]], 6)),
  tolower(mds_cols)
)

local_pca_block <- list(
  schema_version       = 1L,
  layer                = "ghsl_local_pca",
  chrom                = chrom, scale = scale_label,
  n_samples            = n_samp, n_windows = n_win,
  npc                  = npc,
  sample_order         = sample_order,
  pc_loadings_aligned  = pc_blocks,
  lambda               = lapply(seq_len(npc), function(k)
    clean(dt[[paste0("lambda_", k)]], 6)),
  lambda_ratio         = clean(dt$lambda_ratio, 4),
  z                    = clean(dt$max_abs_z, 4),
  z_direct             = clean(dt$div_z_direct, 4),
  max_z_axis           = as.integer(dt$max_z_axis),
  mds_coords           = mds_coords,
  anchor_window_idx    = as.integer(dt$anchor_window_idx[1]),
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

# ── ghsl_envelopes (L1 + L2) ──────────────────────────────────────────────────
read_or_empty <- function(path) if (file.exists(path)) fread(path) else data.table()

envs_to_list <- function(envs, level) {
  if (nrow(envs) == 0L) return(list())
  lapply(seq_len(nrow(envs)), function(k) {
    ws <- if ("win_start" %in% names(envs)) to_idx0(envs$win_start[k]) else NA_integer_
    we <- if ("win_end"   %in% names(envs)) to_idx0(envs$win_end[k])   else NA_integer_
    out <- list(
      level     = level,
      win_start = ws, win_end = we,
      start_bp  = as.integer(envs$start_bp[k] %||% NA_integer_),
      end_bp    = as.integer(envs$end_bp[k]   %||% NA_integer_),
      n_windows = as.integer(envs$n_windows[k] %||% NA_integer_)
    )
    for (extra in c("candidate_id", "l1_id", "l2_id", "mean_sim",
                    "blue_red_ratio", "boundary_score", "validation_status",
                    "status")) {
      if (extra %in% names(envs)) out[[extra]] <- envs[[extra]][k]
    }
    out
  })
}
bnds_to_list <- function(bnds, level) {
  if (nrow(bnds) == 0L) return(list())
  lapply(seq_len(nrow(bnds)), function(k) {
    out <- list(
      level     = level,
      boundary_w = if ("boundary_w" %in% names(bnds)) to_idx0(bnds$boundary_w[k]) else NA_integer_
    )
    for (extra in c("validation_status", "boundary_score", "grow_max_z",
                    "right_max_z", "left_max_z",
                    "right_frac_blue", "left_frac_blue",
                    "is_real", "status")) {
      if (extra %in% names(bnds)) out[[extra]] <- bnds[[extra]][k]
    }
    out
  })
}

l1_envs <- read_or_empty(file.path(l1_dir, paste0(chrom, ".L1_envelopes.tsv")))
l2_envs <- read_or_empty(file.path(l2_dir, paste0(chrom, ".L2_envelopes.tsv")))
l1_bnds <- read_or_empty(file.path(l1_dir, paste0(chrom, ".L1_boundaries.tsv")))
l2_bnds <- read_or_empty(file.path(l2_dir, paste0(chrom, ".L2_boundaries.tsv")))

envs_block <- list(
  schema_version = 1L, layer = "ghsl_envelopes", chrom = chrom,
  l1            = envs_to_list(l1_envs, "L1"),
  l1_boundaries = bnds_to_list(l1_bnds, "L1"),
  l2            = envs_to_list(l2_envs, "L2"),
  l2_boundaries = bnds_to_list(l2_bnds, "L2")
)

# ── tracks (per-window aggregates) ────────────────────────────────────────────
tracks_block <- list(
  ghsl_div_median   = list(values = clean(dt$div_median, 6),  pos_bp = as.integer(dt$mid_bp)),
  ghsl_z_mds        = list(values = clean(dt$max_abs_z, 4),   pos_bp = as.integer(dt$mid_bp)),
  ghsl_z_direct     = list(values = clean(dt$div_z_direct, 4),pos_bp = as.integer(dt$mid_bp)),
  ghsl_lambda_ratio = list(values = clean(dt$lambda_ratio, 4),pos_bp = as.integer(dt$mid_bp))
)
if ("MDS1" %in% names(dt))
  tracks_block$ghsl_mds1 <- list(values = clean(dt$MDS1, 6), pos_bp = as.integer(dt$mid_bp))
if ("MDS2" %in% names(dt))
  tracks_block$ghsl_mds2 <- list(values = clean(dt$MDS2, 6), pos_bp = as.integer(dt$mid_bp))
if ("sketch_novelty" %in% names(dt))
  tracks_block$ghsl_sketch_novelty <- list(values = clean(dt$sketch_novelty, 4),
                                            pos_bp = as.integer(dt$mid_bp))

layers <- c("ghsl_local_pca", "ghsl_envelopes", "tracks")

obj <- list(
  schema_version    = 1L,
  chrom             = chrom,
  n_samples         = n_samp,
  n_windows         = n_win,
  scale             = scale_label,
  `_generated_at`   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  `_generator`      = "STEP_GH_J_atlas_json.R",
  `_layers_present` = layers,
  tracks            = tracks_block,
  ghsl_local_pca    = local_pca_block,
  ghsl_envelopes    = envs_block
)

out_dir <- file.path(outdir, chrom)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_json <- file.path(out_dir, sprintf("%s_phase3_ghsl.json", chrom))
write_json(obj, out_json, auto_unbox = TRUE, na = "null", pretty = FALSE, digits = NA)
fi <- file.info(out_json)
message(sprintf("[GH_J] %s: %s (%.2f MB)", chrom, out_json, fi$size / 1024 / 1024))
message("[GH_J] DONE")
