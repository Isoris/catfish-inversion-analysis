#!/usr/bin/env Rscript
# =============================================================================
# _shared/convert_atlas_json.R  (v10, 2026-05-16)
# =============================================================================
# Convert legacy atlas JSONs (any of the three pipelines' historical schemas)
# to the harmonized v4 schema defined in _shared/lib_atlas_json.R.
#
# Supported inputs:
#   z         schema_version=3 (flat top-level: chrom, samples, windows,
#             sim_thumb, sim_scales, l{1,2}_envelopes, l{1,2}_boundaries, ...)
#   theta_pi  schema_version=2 (layered: theta_pi_local_pca / _envelopes /
#             _cusum / _per_window / _grid_map + tracks)
#   ghsl      schema_version=1 (layered: ghsl_local_pca / _envelopes /
#             _cusum + tracks; scale field present)
#   v4        schema_version=4 (no-op pass-through; useful for batch
#             normalization runs)
#
# Output: harmonized v4 JSON (writes through atlas_json_v4()).
#
# Usage:
#   Rscript _shared/convert_atlas_json.R \
#       --in   <legacy.json>      \
#       --out  <harmonized_v4.json>
#
#   # Batch mode: convert every *.json under a directory tree
#   Rscript _shared/convert_atlas_json.R \
#       --in_dir  <legacy_root> \
#       --out_dir <v4_root> \
#       [--in_pattern "_phase[23]_(theta|ghsl)\\.json$|\\.atlas\\.json$"] \
#       [--inplace]                  # write next to source with .v4.json suffix
#
# The converter is best-effort: it preserves any unmapped legacy fields
# under `extra` so consumers that read non-canonical keys keep working.
# =============================================================================

suppressPackageStartupMessages({ library(jsonlite) })

# ---- locate lib_atlas_json.R ------------------------------------------------
.find_lib_atlas <- function() {
  cand <- Sys.getenv("LIB_ATLAS_JSON", unset = "")
  if (nzchar(cand) && file.exists(cand)) return(cand)
  here <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                   error = function(e) ".")
  same_dir <- file.path(dirname(here), "lib_atlas_json.R")
  if (file.exists(same_dir)) return(same_dir)
  cur <- if (nchar(here) > 1) dirname(here) else getwd()
  for (i in 1:8) {
    c2 <- file.path(cur, "_shared", "lib_atlas_json.R")
    if (file.exists(c2)) return(c2); cur <- dirname(cur)
  }
  stop("Could not locate _shared/lib_atlas_json.R")
}
source(.find_lib_atlas())

# ---- CLI --------------------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, .args)
  if (is.na(i) || i == length(.args)) return(default)
  .args[i + 1L]
}
.has_arg <- function(flag) !is.na(match(flag, .args))

IN       <- .get_arg("--in")
OUT      <- .get_arg("--out")
IN_DIR   <- .get_arg("--in_dir")
OUT_DIR  <- .get_arg("--out_dir")
PATTERN  <- .get_arg("--in_pattern",
                     default = "_phase[23]_(theta|ghsl)\\.json$|\\.atlas\\.json$")
INPLACE  <- .has_arg("--inplace")

if (is.na(IN) && is.na(IN_DIR))
  stop("Usage: --in <file> --out <file>   |   --in_dir <dir> --out_dir <dir>")

# ---- schema detection -------------------------------------------------------
detect_schema <- function(j) {
  sv <- if (is.null(j$schema_version)) NA_integer_ else as.integer(j$schema_version)
  if (!is.na(sv) && sv == 4L) return(list(pipeline = j$pipeline %||% "unknown",
                                          version = 4L))
  if (!is.null(j$pipeline) && j$pipeline %in% c("z", "theta_pi", "ghsl"))
    return(list(pipeline = j$pipeline, version = sv))
  # Layered schemas: detect by which `<pipeline>_*` blocks are present.
  if (any(startsWith(names(j), "ghsl_")))     return(list(pipeline = "ghsl",     version = sv))
  if (any(startsWith(names(j), "theta_pi_"))) return(list(pipeline = "theta_pi", version = sv))
  # Flat z schema: heuristic — has top-level windows + sim_thumb + samples.
  if (!is.null(j$windows) && !is.null(j$sim_thumb) && !is.null(j$samples))
    return(list(pipeline = "z", version = sv))
  list(pipeline = "unknown", version = sv)
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# ---- adapters per legacy schema ---------------------------------------------
adapt_z <- function(j) {
  # z v3 → v4: just stamp pipeline + bump schema_version (handled in writer).
  list(
    pipeline          = "z",
    chrom             = j$chrom,
    n_windows         = j$n_windows,
    n_samples         = j$n_samples,
    samples           = j$samples,
    windows           = j$windows,
    sim_thumb         = j$sim_thumb,
    sim_thumb_n       = j$sim_thumb_n,
    sim_scales        = j$sim_scales,
    default_sim_scale = j$default_sim_scale,
    sim_q_lo          = j$sim_q_lo %||% NA_real_,
    sim_q_hi          = j$sim_q_hi %||% NA_real_,
    z_clip            = j$z_clip   %||% NA_real_,
    z_max_min         = j$z_max_min %||% NA_real_,
    z_column          = j$z_column %||% "max_abs_z",
    has_pc2           = isTRUE(j$has_pc2),
    family_source     = j$family_source %||% "none",
    theta_cutoff      = j$theta_cutoff,
    scale             = NULL,
    theta_range       = j$theta_range,
    has_theta         = isTRUE(j$has_theta),
    l1_envelopes      = j$l1_envelopes,
    l1_boundaries     = j$l1_boundaries,
    l2_envelopes      = j$l2_envelopes,
    l2_boundaries     = j$l2_boundaries,
    tracks            = j$tracks
  )
}

adapt_theta_pi <- function(j) {
  lp <- j$theta_pi_local_pca %||% list()
  ev <- j$theta_pi_envelopes %||% list()
  cu <- j$theta_pi_cusum     %||% NULL
  pw <- j$theta_pi_per_window %||% NULL
  gm <- j$theta_pi_grid_map   %||% NULL

  sample_ids <- lp$sample_order %||% lapply(seq_len(j$n_samples %||% 0L), as.character)
  samples    <- lapply(sample_ids, function(s) list(ind = s))

  # Synthesise per-window records from lp (which holds z, mds_coords, lambda).
  n_win <- as.integer(j$n_windows %||% length(lp$z %||% integer(0)))
  z_vec  <- lp$z       %||% rep(NA_real_, n_win)
  zd_vec <- lp$z_direct %||% rep(NA_real_, n_win)
  lr_vec <- lp$lambda_ratio %||% rep(NA_real_, n_win)
  mds1   <- if (!is.null(lp$mds_coords$mds1)) lp$mds_coords$mds1 else rep(NA_real_, n_win)
  mds2   <- if (!is.null(lp$mds_coords$mds2)) lp$mds_coords$mds2 else rep(NA_real_, n_win)

  windows <- lapply(seq_len(n_win), function(i) list(
    window_idx   = i - 1L,
    max_abs_z    = z_vec[i],
    z_direct     = zd_vec[i],
    lambda_ratio = lr_vec[i],
    MDS1         = mds1[i],
    MDS2         = mds2[i]
  ))

  list(
    pipeline          = "theta_pi",
    chrom             = j$chrom,
    n_windows         = n_win,
    n_samples         = j$n_samples %||% length(sample_ids),
    samples           = samples,
    windows           = windows,
    sim_thumb         = lp$sim_mat,
    sim_thumb_n       = lp$sim_mat_n %||% n_win,
    sim_scales        = NULL,
    default_sim_scale = NULL,
    sim_q_lo          = NA_real_, sim_q_hi = NA_real_,
    z_clip            = NA_real_, z_max_min = NA_real_,
    z_column          = "max_abs_z",
    has_pc2           = FALSE,
    family_source     = "none",
    scale             = j$coarse_scale %||% j$scale,
    l1_envelopes      = ev$l1,
    l1_boundaries     = ev$l1_boundaries,
    l2_envelopes      = ev$l2,
    l2_boundaries     = ev$l2_boundaries,
    tracks            = j$tracks,
    per_window        = pw,
    grid_map          = gm,
    cusum             = cu,
    extra             = list(
      theta_pi_local_pca   = lp,
      theta_pi_envelopes   = ev,
      theta_pi_cusum       = cu,
      theta_pi_per_window  = pw,
      theta_pi_grid_map    = gm,
      coarse_scale         = j$coarse_scale,
      dense_scale          = j$dense_scale,
      n_windows_coarse     = j$n_windows_coarse,
      n_windows_dense      = j$n_windows_dense
    )
  )
}

adapt_ghsl <- function(j) {
  lp <- j$ghsl_local_pca %||% list()
  ev <- j$ghsl_envelopes %||% list()
  cu <- j$ghsl_cusum     %||% NULL

  sample_ids <- lp$sample_order %||% lapply(seq_len(j$n_samples %||% 0L), as.character)
  samples    <- lapply(sample_ids, function(s) list(ind = s))

  n_win <- as.integer(j$n_windows %||% length(lp$z %||% integer(0)))
  z_vec  <- lp$z       %||% rep(NA_real_, n_win)
  zd_vec <- lp$z_direct %||% rep(NA_real_, n_win)
  lr_vec <- lp$lambda_ratio %||% rep(NA_real_, n_win)
  mds1   <- if (!is.null(lp$mds_coords$mds1)) lp$mds_coords$mds1 else rep(NA_real_, n_win)
  mds2   <- if (!is.null(lp$mds_coords$mds2)) lp$mds_coords$mds2 else rep(NA_real_, n_win)

  windows <- lapply(seq_len(n_win), function(i) list(
    window_idx   = i - 1L,
    max_abs_z    = z_vec[i],
    z_direct     = zd_vec[i],
    lambda_ratio = lr_vec[i],
    MDS1         = mds1[i],
    MDS2         = mds2[i]
  ))

  list(
    pipeline          = "ghsl",
    chrom             = j$chrom,
    n_windows         = n_win,
    n_samples         = j$n_samples %||% length(sample_ids),
    samples           = samples,
    windows           = windows,
    sim_thumb         = lp$sim_mat,
    sim_thumb_n       = lp$sim_mat_n %||% n_win,
    sim_scales        = NULL,
    default_sim_scale = NULL,
    sim_q_lo          = NA_real_, sim_q_hi = NA_real_,
    z_clip            = NA_real_, z_max_min = NA_real_,
    z_column          = "max_abs_z",
    has_pc2           = FALSE,
    family_source     = "none",
    scale             = j$scale,
    l1_envelopes      = ev$l1,
    l1_boundaries     = ev$l1_boundaries,
    l2_envelopes      = ev$l2,
    l2_boundaries     = ev$l2_boundaries,
    tracks            = j$tracks,
    cusum             = cu,
    extra             = list(
      ghsl_local_pca = lp,
      ghsl_envelopes = ev,
      ghsl_cusum     = cu
    )
  )
}

# ---- driver -----------------------------------------------------------------
convert_one <- function(in_path, out_path) {
  j <- jsonlite::read_json(in_path, simplifyVector = FALSE)
  d <- detect_schema(j)
  if (d$pipeline == "unknown")
    stop("Cannot detect pipeline for ", in_path,
         " (schema_version=", d$version, ")")

  if (!is.na(d$version) && d$version == 4L) {
    message("[conv] ", in_path, " already v4 — copying through")
    args_list <- adapt_z(j)            # v4 z-style envelope = same shape
    args_list$pipeline <- d$pipeline
  } else if (d$pipeline == "z") {
    args_list <- adapt_z(j)
  } else if (d$pipeline == "theta_pi") {
    args_list <- adapt_theta_pi(j)
  } else if (d$pipeline == "ghsl") {
    args_list <- adapt_ghsl(j)
  } else {
    stop("Unhandled pipeline: ", d$pipeline)
  }

  args_list$generator <- paste0("convert_atlas_json.R<-", basename(in_path))
  args_list$out_path  <- out_path
  do.call(atlas_json_v4, args_list)
  message("[conv] ", in_path, " (", d$pipeline, "/v",
          d$version %||% "?", ") -> ", out_path)
}

if (!is.na(IN)) {
  if (is.na(OUT)) stop("--out is required when --in is given")
  convert_one(IN, OUT)
} else {
  if (!INPLACE && is.na(OUT_DIR))
    stop("--out_dir is required for batch mode (or pass --inplace)")
  files <- list.files(IN_DIR, pattern = PATTERN,
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0L)
    stop("No files matching '", PATTERN, "' under ", IN_DIR)
  message("[conv] batch: ", length(files), " files under ", IN_DIR)
  for (f in files) {
    rel <- if (INPLACE) {
      sub("\\.json$", ".v4.json", f)
    } else {
      file.path(OUT_DIR, sub(paste0("^", normalizePath(IN_DIR, mustWork = FALSE),
                                    "/?"), "",
                             normalizePath(f, mustWork = FALSE)))
    }
    dir.create(dirname(rel), recursive = TRUE, showWarnings = FALSE)
    tryCatch(convert_one(f, rel),
             error = function(e) message("[conv][ERR] ", f, ": ", conditionMessage(e)))
  }
}
