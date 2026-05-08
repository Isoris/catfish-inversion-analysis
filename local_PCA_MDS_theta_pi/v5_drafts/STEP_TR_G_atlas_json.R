#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_G_atlas_json.R   (drafts/v5)
# =============================================================================
# Atlas JSON exporter — combines outputs from TR_B (precomp), TR_C (L1),
# TR_D (L2), TR_E (carriers), TR_F (CUSUM) into one consolidated per-chrom
# JSON the page-12 atlas reads. Mirrors what local_PCA_MDS_z's 08b atlas
# JSON exporter does for the z-blocks pipeline.
#
# All emitted window indices are 0-INDEXED (atlas is JS).
#
# Reads (z-blocks-shaped layout produced by TR_B v5):
#   <PRECOMP_DIR>/<chr>.precomp.rds
#   <L1_DIR>/<chr>.L1_envelopes.tsv
#   <L1_DIR>/<chr>.L1_boundaries.tsv         (optional; included if present)
#   <L2_DIR>/<chr>.L2_envelopes.tsv
#   <L2_DIR>/<chr>.L2_boundaries.tsv         (optional)
#   <CARRIERS_DIR>/<chr>.carrier_assignments.tsv
#   <CUSUM_DIR>/<chr>.cusum_per_sample.tsv.gz
#   <CUSUM_DIR>/<chr>.cusum_boundary_dist.tsv
#   $THETA_TSV_DIR/theta_native.<chr>.<scale>.tsv.gz   (for theta_pi_per_window)
# Writes:
#   <JSON_OUT_DIR>/<chr>/<chr>_phase2_theta.json
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}

CHR          <- get_arg("--chr")
PRECOMP_DIR  <- get_arg("--precomp_dir")
L1_DIR       <- get_arg("--l1_dir")
L2_DIR       <- get_arg("--l2_dir")
CARRIERS_DIR <- get_arg("--carriers_dir")
CUSUM_DIR    <- get_arg("--cusum_dir")
JSON_OUT_DIR <- get_arg("--json_out_dir")

OUTROOT       <- Sys.getenv("OUTROOT",       unset = NA)
THETA_TSV_DIR <- Sys.getenv("THETA_TSV_DIR", unset = NA)
PESTPG_SCALE  <- Sys.getenv("PESTPG_SCALE",  unset = "win10000.step2000")
stopifnot(!is.na(OUTROOT))

if (is.na(PRECOMP_DIR))  PRECOMP_DIR  <- file.path(OUTROOT, "precomp")
if (is.na(L1_DIR))       L1_DIR       <- file.path(OUTROOT, "L1_detect")
if (is.na(L2_DIR))       L2_DIR       <- file.path(OUTROOT, "L2_detect")
if (is.na(CARRIERS_DIR)) CARRIERS_DIR <- file.path(OUTROOT, "carriers")
if (is.na(CUSUM_DIR))    CUSUM_DIR    <- file.path(OUTROOT, "cusum")
if (is.na(JSON_OUT_DIR)) JSON_OUT_DIR <- file.path(OUTROOT, "04_atlas_json")

chroms <- if (!is.na(CHR)) CHR else
  sub("\\.precomp\\.rds$", "",
      list.files(PRECOMP_DIR, pattern = "\\.precomp\\.rds$"))

clean <- function(x, digits = 6) {
  out <- round(as.numeric(x), digits); out[!is.finite(out)] <- NA_real_; out
}

# Convert 1-indexed win_start/win_end (z-blocks output) to 0-indexed for JSON.
to_idx0 <- function(x) if (length(x) == 0L) integer(0) else as.integer(x) - 1L

for (chrom in chroms) {
  rds <- file.path(PRECOMP_DIR, paste0(chrom, ".precomp.rds"))
  if (!file.exists(rds)) { message("[TR_G] ", chrom, ": no precomp — skip"); next }

  precomp <- readRDS(rds); dt <- precomp$dt
  n_win  <- precomp$n_windows; n_samp <- precomp$n_samples
  sample_order <- precomp$sample_order
  npc <- as.integer(precomp$npc %||% 2L)

  # Per-sample PC1..PCk arrays per window.
  pc_blocks <- vector("list", npc)
  for (k in seq_len(npc)) {
    cols <- grep(paste0("^PC_", k, "_"), names(dt), value = TRUE)
    if (length(cols) == 0) next
    pc_mat <- as.matrix(dt[, ..cols])
    pc_blocks[[k]] <- lapply(seq_len(n_win), function(wi) clean(pc_mat[wi, ], 6))
  }

  # ─ theta_pi_per_window (rebuild from theta_native TSV) ─
  per_window_block <- NULL
  if (!is.na(THETA_TSV_DIR)) {
    tsv <- file.path(THETA_TSV_DIR,
                     sprintf("theta_native.%s.%s.tsv.gz", chrom, PESTPG_SCALE))
    if (file.exists(tsv)) {
      long_dt <- fread(tsv); long_dt <- long_dt[chrom == ..chrom]
      if (nrow(long_dt) > 0) {
        samp_to_row <- setNames(seq_along(sample_order), sample_order)
        win_to_col  <- setNames(seq_along(dt$window_idx), as.character(dt$window_idx))
        theta_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                            dimnames = list(sample_order, NULL))
        rows <- samp_to_row[long_dt$sample]
        cols <- win_to_col[as.character(long_dt$window_idx)]
        good <- !is.na(rows) & !is.na(cols)
        theta_mat[cbind(rows[good], cols[good])] <- long_dt$theta_pi[good]
        per_window_block <- list(
          schema_version = 2L, layer = "theta_pi_per_window", chrom = chrom,
          scale = PESTPG_SCALE,
          n_samples = n_samp, n_windows = n_win,
          sample_ids = sample_order,
          values    = clean(as.numeric(t(theta_mat)), 6),         # row-major over windows
          windows   = lapply(seq_len(n_win), function(wi) list(
            idx = as.integer(dt$window_idx[wi]),
            start_bp = as.integer(dt$start_bp[wi]),
            end_bp   = as.integer(dt$end_bp[wi])
          )),
          samples = lapply(seq_len(n_samp), function(s) list(
            sample_id = sample_order[s],
            theta_pi  = clean(theta_mat[s, ], 6)
          ))
        )
      }
    }
  }

  # ─ theta_pi_local_pca ─
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
    chrom                = chrom, scale = PESTPG_SCALE,
    n_samples            = n_samp, n_windows = n_win,
    npc                  = npc,
    sample_order         = sample_order,
    pc_loadings_aligned  = pc_blocks,                # list of NPC matrices, sign-aligned
    lambda               = lapply(seq_len(npc), function(k)
      clean(dt[[paste0("lambda_", k)]], 6)),
    lambda_ratio         = clean(dt$lambda_ratio, 4),
    z                    = clean(dt$max_abs_z, 4),       # MDS-axis |Z| (z-blocks definition)
    z_theta_direct       = clean(dt$theta_z_direct, 4),  # θπ-direct |Z| (atlas alt track)
    max_z_axis           = as.integer(dt$max_z_axis),
    mds_coords           = list(
      mds1 = clean(dt$MDS1, 6), mds2 = clean(dt$MDS2, 6),
      mds3 = clean(dt$MDS3, 6), mds4 = clean(dt$MDS4, 6),
      mds5 = clean(dt$MDS5, 6)
    ),
    anchor_window_idx    = as.integer(dt$anchor_window_idx[1]),    # already 0-indexed
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

  # ─ theta_pi_envelopes (L1 + L2 from z-blocks-style detector outputs) ─
  read_or_empty <- function(path) if (file.exists(path)) fread(path) else data.table()
  l1_envs <- read_or_empty(file.path(L1_DIR, paste0(chrom, ".L1_envelopes.tsv")))
  l2_envs <- read_or_empty(file.path(L2_DIR, paste0(chrom, ".L2_envelopes.tsv")))
  l1_bnds <- read_or_empty(file.path(L1_DIR, paste0(chrom, ".L1_boundaries.tsv")))
  l2_bnds <- read_or_empty(file.path(L2_DIR, paste0(chrom, ".L2_boundaries.tsv")))

  # z-blocks 04/06 emit win_start/win_end as 1-indexed window positions.
  # Convert to 0-indexed for the atlas.
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

  envs_block <- list(
    schema_version = 2L, layer = "theta_pi_envelopes", chrom = chrom,
    l1            = envs_to_list(l1_envs, "L1"),
    l1_boundaries = bnds_to_list(l1_bnds, "L1"),
    l2            = envs_to_list(l2_envs, "L2"),
    l2_boundaries = bnds_to_list(l2_bnds, "L2")
  )

  # ─ theta_pi_cusum ─
  ca_dt <- read_or_empty(file.path(CARRIERS_DIR, paste0(chrom, ".carrier_assignments.tsv")))
  cs_dt <- if (file.exists(file.path(CUSUM_DIR, paste0(chrom, ".cusum_per_sample.tsv.gz"))))
              fread(file.path(CUSUM_DIR, paste0(chrom, ".cusum_per_sample.tsv.gz"))) else data.table()
  bd_dt <- read_or_empty(file.path(CUSUM_DIR, paste0(chrom, ".cusum_boundary_dist.tsv")))

  cusum_candidates <- list()
  if (nrow(l2_envs) > 0L && "candidate_id" %in% names(l2_envs)) {
    for (k in seq_len(nrow(l2_envs))) {
      cid <- l2_envs$candidate_id[k]
      ca_sub <- ca_dt[candidate_id == cid]
      cs_sub <- if (nrow(cs_dt) > 0) cs_dt[candidate_id == cid] else data.table()
      bd_sub <- if (nrow(bd_dt) > 0) bd_dt[candidate_id == cid] else data.table()
      bands <- list()
      for (b in unique(ca_sub$band)) {
        members <- ca_sub[band == b]
        cs_b    <- cs_sub[band == b]
        bd_b    <- bd_sub[band == b]
        side_block <- function(side) {
          r <- bd_b[side == ..side]
          if (nrow(r) == 0L) return(list(n_carriers = 0L))
          list(
            n_carriers     = as.integer(r$n_carriers),
            median_bp      = as.integer(r$median_bp),
            iqr_kb         = round(r$iqr_kb, 2),
            spread_class   = r$spread_class,
            peak_strength  = round(r$peak_strength, 3),
            consensus_cp_bp       = as.integer(r$consensus_cp_bp),
            consensus_strength    = round(r$consensus_strength, 3),
            consensus_informative = as.logical(r$consensus_informative)
          )
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
    schema_version = 2L, layer = "theta_pi_cusum", chrom = chrom,
    n_candidates   = length(cusum_candidates),
    candidates     = cusum_candidates
  )

  # ─ tracks ─
  tracks_block <- list(
    theta_pi_median = list(values = clean(dt$theta_pi_median, 6),
                            pos_bp = as.integer(dt$mid_bp)),
    theta_pi_z_mds  = list(values = clean(dt$max_abs_z, 4),
                            pos_bp = as.integer(dt$mid_bp)),
    theta_pi_z_direct = list(values = clean(dt$theta_z_direct, 4),
                              pos_bp = as.integer(dt$mid_bp)),
    theta_pi_lambda_ratio = list(values = clean(dt$lambda_ratio, 4),
                                  pos_bp = as.integer(dt$mid_bp)),
    theta_pi_mds1 = list(values = clean(dt$MDS1, 6), pos_bp = as.integer(dt$mid_bp)),
    theta_pi_mds2 = list(values = clean(dt$MDS2, 6), pos_bp = as.integer(dt$mid_bp))
  )

  layers <- c("theta_pi_local_pca", "theta_pi_envelopes", "theta_pi_cusum", "tracks")
  if (!is.null(per_window_block)) layers <- c("theta_pi_per_window", layers)

  obj <- list(
    schema_version = 2L,
    chrom = chrom, n_samples = n_samp, n_windows = n_win, scale = PESTPG_SCALE,
    `_generated_at`  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    `_generator`     = "STEP_TR_G_atlas_json.R (v5)",
    `_layers_present` = layers,
    tracks                = tracks_block,
    theta_pi_local_pca    = local_pca_block,
    theta_pi_envelopes    = envs_block,
    theta_pi_cusum        = cusum_block
  )
  if (!is.null(per_window_block)) obj$theta_pi_per_window <- per_window_block

  out_dir <- file.path(JSON_OUT_DIR, chrom)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_json <- file.path(out_dir, sprintf("%s_phase2_theta.json", chrom))
  write_json(obj, out_json, auto_unbox = TRUE, na = "null", pretty = FALSE, digits = NA)
  fi <- file.info(out_json)
  message(sprintf("[TR_G] %s: %s (%.2f MB)", chrom, out_json, fi$size / 1024 / 1024))
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
message("[TR_G] DONE")