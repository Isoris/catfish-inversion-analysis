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

chrom         <- get_arg("--chr")
precomp_dir   <- get_arg("--precomp_dir")
l1_dir        <- get_arg("--l1_dir")
l2_dir        <- get_arg("--l2_dir")
carriers_dir  <- get_arg("--carriers_dir")
cusum_dir     <- get_arg("--cusum_dir")
outdir        <- get_arg("--outdir")
scale_label   <- get_arg("--scale", Sys.getenv("GHSL_PCA_SCALE_FULL", unset = "s50"))

stopifnot(!is.na(chrom), !is.na(precomp_dir), !is.na(outdir))
if (is.na(l1_dir))       l1_dir       <- file.path(dirname(precomp_dir), "04_L1_detect")
if (is.na(l2_dir))       l2_dir       <- file.path(dirname(precomp_dir), "06_L2_detect")
# carriers/cusum dirs are optional — when missing, the carriers/cusum blocks
# are simply omitted (e.g. if you haven't run TR_H/TR_I yet).

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

# ── ghsl_cusum: per-L2-candidate band + per-sample CP records ─────────────────
# Mirrors theta_pi_cusum: candidates × bands, each band emits per_sample
# (one CP per carrier) and boundary_5_prime / boundary_3_prime consensus
# blocks (with the 80% sample-consensus fields from TR_I May 2026+).
cusum_block <- NULL
if (!is.na(carriers_dir) || !is.na(cusum_dir)) {
  ca_path <- if (!is.na(carriers_dir))
                file.path(carriers_dir, paste0(chrom, ".carrier_assignments.tsv"))
             else NA_character_
  cs_path <- if (!is.na(cusum_dir))
                file.path(cusum_dir, paste0(chrom, ".cusum_per_sample.tsv.gz"))
             else NA_character_
  bd_path <- if (!is.na(cusum_dir))
                file.path(cusum_dir, paste0(chrom, ".cusum_boundary_dist.tsv"))
             else NA_character_

  ca_dt <- if (!is.na(ca_path) && file.exists(ca_path)) fread(ca_path) else data.table()
  cs_dt <- if (!is.na(cs_path) && file.exists(cs_path)) fread(cs_path) else data.table()
  bd_dt <- if (!is.na(bd_path) && file.exists(bd_path)) fread(bd_path) else data.table()

  cusum_candidates <- list()
  if (nrow(l2_envs) > 0L && nrow(ca_dt) > 0L && "candidate_id" %in% names(l2_envs)) {
    for (k in seq_len(nrow(l2_envs))) {
      cid <- l2_envs$candidate_id[k]
      ca_sub <- ca_dt[candidate_id == cid]
      cs_sub <- if (nrow(cs_dt) > 0L) cs_dt[candidate_id == cid] else data.table()
      bd_sub <- if (nrow(bd_dt) > 0L) bd_dt[candidate_id == cid] else data.table()
      if (nrow(ca_sub) == 0L) next
      bands <- list()
      for (b in unique(ca_sub$band)) {
        members <- ca_sub[band == b]
        cs_b    <- cs_sub[band == b]
        bd_b    <- bd_sub[band == b]
        side_block <- function(which_side) {
          r <- bd_b[side == which_side]
          if (nrow(r) == 0L) return(list(n_carriers = 0L))
          out <- list(
            n_carriers     = as.integer(r$n_carriers),
            median_bp      = as.integer(r$median_bp),
            iqr_kb         = round(r$iqr_kb, 2),
            spread_class   = r$spread_class,
            peak_strength  = round(r$peak_strength, 3),
            consensus_cp_bp       = as.integer(r$consensus_cp_bp),
            consensus_strength    = round(r$consensus_strength, 3),
            consensus_informative = as.logical(r$consensus_informative)
          )
          for (col in c("n_carriers_supporting", "frac_carriers_supporting",
                        "consensus_tol_kb", "passes_frac_threshold")) {
            if (col %in% names(r)) out[[col]] <- r[[col]]
          }
          out
        }
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
          boundary_5_prime = side_block("5_prime"),
          boundary_3_prime = side_block("3_prime")
        )
      }
      cusum_candidates[[length(cusum_candidates) + 1L]] <- list(
        candidate_id = cid, chrom = chrom,
        start_bp     = as.integer(l2_envs$start_bp[k] %||% NA_integer_),
        end_bp       = as.integer(l2_envs$end_bp[k]   %||% NA_integer_),
        bands        = bands
      )
    }
  }
  cusum_block <- list(
    schema_version = 1L, layer = "ghsl_cusum", chrom = chrom,
    n_candidates   = length(cusum_candidates),
    candidates     = cusum_candidates
  )
}

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

# ── Harmonized schema v4 envelope (v10) ──────────────────────────────────────
# The original GHSL JSON layered everything under ghsl_local_pca /
# ghsl_envelopes / ghsl_cusum keys with schema_version=1. v4 merges those
# layered keys into the canonical z-style top-level (windows / sim_thumb /
# sim_scales / l{1,2}_envelopes / etc.) and keeps the legacy ghsl_*
# blocks as `extra` for back-compat with current consumers.
.find_lib_atlas <- function() {
  cand <- Sys.getenv("LIB_ATLAS_JSON", unset = "")
  if (nzchar(cand) && file.exists(cand)) return(cand)
  cand2 <- file.path(Sys.getenv("SCRIPT_DIR_SHARED", unset = ""),
                     "lib_atlas_json.R")
  if (nzchar(Sys.getenv("SCRIPT_DIR_SHARED")) && file.exists(cand2)) return(cand2)
  cur <- getwd()
  for (i in 1:8) {
    c3 <- file.path(cur, "_shared", "lib_atlas_json.R")
    if (file.exists(c3)) return(c3)
    cur <- dirname(cur)
  }
  stop("Could not find _shared/lib_atlas_json.R; set SCRIPT_DIR_SHARED or LIB_ATLAS_JSON")
}
source(.find_lib_atlas())

# Map layered GHSL blocks → harmonized top-level fields.
samples_json_gh <- lapply(sample_order, function(s) list(ind = s))
windows_json_gh <- lapply(seq_len(n_win), function(i) {
  rec <- list(
    window_idx    = as.integer(dt$window_idx[i]),
    start_bp      = as.integer(dt$start_bp[i]),
    end_bp        = as.integer(dt$end_bp[i]),
    mid_bp        = as.integer(dt$mid_bp[i]),
    max_abs_z     = clean(dt$max_abs_z[i], 4),
    z_direct      = clean(dt$div_z_direct[i], 4),
    lambda_ratio  = clean(dt$lambda_ratio[i], 4)
  )
  if ("MDS1" %in% names(dt)) rec$MDS1 <- clean(dt$MDS1[i], 6)
  if ("MDS2" %in% names(dt)) rec$MDS2 <- clean(dt$MDS2[i], 6)
  rec
})

# sim_thumb = the unsmoothed sim_mat as a flat vector (or NULL if banded).
sim_thumb_v   <- if (precomp$sim_mat_format == "upper_triangle_float32" &&
                     !is.null(precomp$sim_mat)) clean(as.vector(t(precomp$sim_mat)), 4) else NULL
sim_thumb_n_v <- if (!is.null(sim_thumb_v)) n_win else 0L

# Coerce envelopes to a v4-shaped list (already done by envs_to_list above).
out_dir <- file.path(outdir, chrom)
out_json <- file.path(out_dir, sprintf("%s_phase3_ghsl.json", chrom))

# carriers / cusum: GHSL has a single ghsl_cusum block today (no separate
# carriers). Pass it through as the cusum extension so the v4 envelope
# carries it without losing the existing field shape.
atlas_json_v4(
  pipeline          = "ghsl",
  chrom             = chrom,
  n_windows         = n_win,
  n_samples         = n_samp,
  samples           = samples_json_gh,
  windows           = windows_json_gh,
  sim_thumb         = sim_thumb_v,
  sim_thumb_n       = sim_thumb_n_v,
  sim_scales        = NULL,        # GHSL exporter emits single-scale today
  default_sim_scale = NULL,
  sim_q_lo          = NA_real_, sim_q_hi = NA_real_,
  z_clip            = NA_real_, z_max_min = NA_real_,
  z_column          = "max_abs_z",
  has_pc2           = FALSE,
  family_source     = "none",
  scale             = scale_label,
  l1_envelopes      = envs_block$l1,
  l1_boundaries     = envs_block$l1_boundaries,
  l2_envelopes      = envs_block$l2,
  l2_boundaries     = envs_block$l2_boundaries,
  tracks            = if (length(tracks_block) > 0) tracks_block else NULL,
  cusum             = cusum_block,
  generator         = "STEP_GH_J_atlas_json.R",
  out_path          = out_json,
  # Back-compat: keep the legacy layered GHSL blocks for current consumers.
  extra             = list(
    ghsl_local_pca = local_pca_block,
    ghsl_envelopes = envs_block,
    ghsl_cusum     = cusum_block
  )
)
fi <- file.info(out_json)
message(sprintf("[GH_J v4] %s: %s (%.2f MB)", chrom, out_json, fi$size / 1024 / 1024))
message("[GH_J] DONE")
