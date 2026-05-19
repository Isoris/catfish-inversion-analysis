#!/usr/bin/env Rscript
# =============================================================================
# _shared/lib_atlas_json.R — harmonized atlas JSON schema v4 (2026-05-16)
# =============================================================================
# Shared schema used by all three atlas exporters:
#   local_PCA_MDS_z/STEP_ZO_M_export_atlas_json.R        (pipeline = "z")
#   local_PCA_MDS_theta_pi/STEP_TR_J_atlas_json.R        (pipeline = "theta_pi")
#   local_PCA_MDS_GHSL/STEP_GH_J_atlas_json.R            (pipeline = "ghsl")
#
# Design: take z (schema v3) as the canonical envelope — flat top-level
# fields the atlas page-1 viewer already consumes — and add a small set of
# optional pipeline-specific extension blocks (per_window theta matrix for
# theta_pi; carriers/cusum blocks for theta_pi/GHSL; layered tracks for all).
#
# Top-level v4 keys (REQUIRED unless noted):
#
#   schema_version    : integer = 4
#   pipeline          : "z" | "theta_pi" | "ghsl"
#   chrom             : str
#   n_windows         : int
#   n_samples         : int
#   samples           : [{ind, cga, ancestry, family_id}, ...]
#                       (theta_pi/GHSL pipelines may omit identity-merge fields
#                       and emit just {ind}; family_source then = "none")
#   family_source     : "pairs" | "family" | "none"
#   theta_cutoff      : num | null
#   has_pc2           : bool
#   z_column          : str               name of the column carrying |Z|
#                                          ("max_abs_z" by default)
#   scale             : str | null        optional pipeline scale label
#                                          (e.g. theta_pi "win50000.step10000",
#                                          GHSL "s50"). null for z.
#
#   windows           : [{...per-window record (id, mid_bp, max_abs_z, ...)}]
#   sim_thumb         : flat numeric vector (length sim_thumb_n^2)
#   sim_thumb_n       : int
#   sim_scales        : {nn40: <vec>, nn80: <vec>, ...} keyed by NN_SIM_SCALES
#   default_sim_scale : "nn80" (typically)
#
#   sim_q_lo, sim_q_hi, z_clip, z_max_min : numeric color/clip controls
#
#   theta_range       : [lo, hi] | null     (z pipeline only; null otherwise)
#   has_theta         : bool                true iff theta_by_sample present
#
#   l1_envelopes      : [...] | null
#   l1_boundaries     : [...] | null
#   l2_envelopes      : [...] | null
#   l2_boundaries     : [...] | null
#   has_l1_envelopes, has_l2_envelopes, has_l1_boundaries, has_l2_boundaries : bool
#
#   tracks            : { <track_name>: {values: [...], pos_bp: [...]}, ... } | null
#                       Optional layered per-window aggregates.
#
# OPTIONAL extension blocks (NULL when absent):
#   per_window        : (theta_pi)  per-sample × dense-window theta_pi matrix
#                                    {matrix, sample_order, window_idx, scale}
#   grid_map          : (theta_pi)  coarse↔dense window-index lookup
#                                    {coarse_to_dense, dense_to_coarse}
#   carriers          : (theta_pi/GHSL) carrier_assignments block
#                                    {candidates: [{candidate_id, bands: [...]}]}
#   cusum             : (theta_pi/GHSL) cusum bands + per-sample CP records
#                                    {candidates: [{candidate_id, bands: [...]}]}
#
# META keys (always present):
#   _generated_at     : ISO-8601 timestamp
#   _generator        : script filename emitting the JSON
#   _layers_present   : sorted vector of which optional blocks are non-null
#                       e.g. ["per_window", "carriers", "cusum", "tracks"]
#
# Helper API
# ----------
# atlas_json_v4(..., pipeline = c("z","theta_pi","ghsl"), out_path,
#               per_window = NULL, grid_map = NULL,
#               carriers = NULL, cusum = NULL, tracks = NULL, ...)
#   - Validates the minimum required fields, fills NULLs to defaults,
#     stamps schema_version=4, _generated_at, _generator,
#     _layers_present, and writes via jsonlite::write_json(auto_unbox=TRUE,
#     na="null", digits=6).
# =============================================================================

suppressPackageStartupMessages({ library(jsonlite) })

ATLAS_JSON_SCHEMA_VERSION <- 4L

# Tiny helper: drop NULLs from a list (so JSON stays terse).
.atlas_drop_nulls <- function(x) x[!vapply(x, is.null, logical(1))]

# ---- public: write a harmonized v4 atlas JSON ------------------------------
atlas_json_v4 <- function(
  pipeline,                         # "z" | "theta_pi" | "ghsl"
  chrom,
  n_windows,
  n_samples,
  samples,                          # list of per-sample records
  windows,                          # list of per-window records
  sim_thumb,
  sim_thumb_n,
  sim_scales        = NULL,
  default_sim_scale = NULL,
  sim_q_lo          = NA_real_,
  sim_q_hi          = NA_real_,
  z_clip            = NA_real_,
  z_max_min         = NA_real_,
  z_column          = "max_abs_z",
  has_pc2           = FALSE,
  family_source     = "none",
  theta_cutoff      = NULL,
  scale             = NULL,
  theta_range       = NULL,
  has_theta         = FALSE,
  l1_envelopes      = NULL,
  l1_boundaries     = NULL,
  l2_envelopes      = NULL,
  l2_boundaries     = NULL,
  tracks            = NULL,
  per_window        = NULL,         # theta_pi-only extension
  grid_map          = NULL,         # theta_pi-only extension
  carriers          = NULL,         # theta_pi/GHSL extension
  cusum             = NULL,         # theta_pi/GHSL extension
  generator         = NA_character_,
  out_path          = NULL,
  extra             = NULL,         # named list of extra top-level keys (back-compat)
  pretty            = FALSE,
  digits            = 6
) {
  stopifnot(pipeline %in% c("z", "theta_pi", "ghsl"))
  if (is.null(out_path) || !nzchar(out_path))
    stop("atlas_json_v4: out_path is required")

  layers_present <- character(0)
  if (!is.null(per_window)) layers_present <- c(layers_present, "per_window")
  if (!is.null(grid_map))   layers_present <- c(layers_present, "grid_map")
  if (!is.null(carriers))   layers_present <- c(layers_present, "carriers")
  if (!is.null(cusum))      layers_present <- c(layers_present, "cusum")
  if (!is.null(tracks))     layers_present <- c(layers_present, "tracks")
  layers_present <- sort(layers_present)

  obj <- list(
    schema_version    = ATLAS_JSON_SCHEMA_VERSION,
    pipeline          = pipeline,
    chrom             = chrom,
    n_windows         = as.integer(n_windows),
    n_samples         = as.integer(n_samples),
    samples           = samples,
    family_source     = family_source,
    theta_cutoff      = theta_cutoff,
    has_pc2           = isTRUE(has_pc2),
    z_column          = z_column,
    scale             = scale,
    windows           = windows,
    sim_thumb         = sim_thumb,
    sim_thumb_n       = as.integer(sim_thumb_n),
    sim_scales        = sim_scales,
    default_sim_scale = default_sim_scale,
    sim_q_lo          = sim_q_lo,
    sim_q_hi          = sim_q_hi,
    z_clip            = z_clip,
    z_max_min         = z_max_min,
    theta_range       = theta_range,
    has_theta         = isTRUE(has_theta),
    l1_envelopes      = l1_envelopes,
    l1_boundaries     = l1_boundaries,
    l2_envelopes      = l2_envelopes,
    l2_boundaries     = l2_boundaries,
    has_l1_envelopes  = !is.null(l1_envelopes),
    has_l2_envelopes  = !is.null(l2_envelopes),
    has_l1_boundaries = !is.null(l1_boundaries),
    has_l2_boundaries = !is.null(l2_boundaries),
    tracks            = tracks,
    per_window        = per_window,
    grid_map          = grid_map,
    carriers          = carriers,
    cusum             = cusum,
    `_generated_at`   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    `_generator`      = if (is.na(generator)) "lib_atlas_json.R" else generator,
    `_layers_present` = layers_present
  )

  # Merge extra (caller-supplied legacy fields) without overwriting v4 keys.
  if (!is.null(extra)) {
    nx <- setdiff(names(extra), names(obj))
    for (k in nx) obj[[k]] <- extra[[k]]
  }

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(obj, out_path,
                       auto_unbox = TRUE, na = "null",
                       pretty = pretty, digits = digits)
  invisible(out_path)
}
