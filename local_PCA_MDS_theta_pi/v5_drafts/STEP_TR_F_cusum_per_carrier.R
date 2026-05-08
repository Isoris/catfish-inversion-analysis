#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_F_cusum_per_carrier.R   (drafts/v5)
# =============================================================================
# Per-carrier CUSUM inside each L2 candidate. Wraps lib_persample_cusum.R
# (the kernel from spec_cusum/) and partitions the run by carrier band so
# we get a per-band breakpoint distribution per candidate.
#
# Pipeline position: after TR_B (precomp), TR_D (L2), TR_E (carriers).
# Each L2 candidate has a list of (sample_id → band) assignments from TR_E;
# this script slices the precomp's θπ matrix to the candidate windows ×
# the band's member samples, runs CUSUM, and writes:
#
#   <OUTROOT>/03_per_chrom/<chr>/cusum_per_sample.tsv.gz
#       Per-sample changepoints with band attribution. One row per
#       (sample × candidate). Columns mirror STEP_T05 + a `band` column
#       so downstream consumers (DC06, atlas) can group by carrier class.
#
#   <OUTROOT>/03_per_chrom/<chr>/cusum_boundary_dist.tsv
#       Per-candidate × per-band × per-side empirical distribution shape:
#         candidate_id chrom band side n_carriers median_bp iqr_kb
#         spread_class peak_strength
#       spread_class ∈ {tight, intermediate, ragged} per CUSUM_SPEC §5.2.
#
# Reads:   <OUTROOT>/03_per_chrom/<chr>/precomp.rds
#          <OUTROOT>/03_per_chrom/<chr>/L2_envelopes.tsv
#          <OUTROOT>/03_per_chrom/<chr>/carrier_assignments.tsv
#          v5_drafts/lib_persample_cusum.R   (the CUSUM kernel)
# Usage:   Rscript STEP_TR_F_cusum_per_carrier.R --chrom <CHR>
#                                                [--lib <path>]
#                                                [--bands "MID_DIV,HIGH_DIV"]
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
CHROM       <- NULL
LIB_PATH    <- NA_character_
BANDS_RUN   <- NA_character_  # comma-separated; default = all bands
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--chrom" && i < length(args)) { CHROM     <- args[i + 1]; i <- i + 2L }
  else if (a == "--lib"   && i < length(args)) { LIB_PATH  <- args[i + 1]; i <- i + 2L }
  else if (a == "--bands" && i < length(args)) { BANDS_RUN <- args[i + 1]; i <- i + 2L }
  else { i <- i + 1L }
}

OUTROOT <- Sys.getenv("OUTROOT", unset = NA)
stopifnot(!is.na(OUTROOT))

# Locate the CUSUM kernel.
if (is.na(LIB_PATH)) {
  guesses <- c(
    file.path(dirname(sys.frame(1)$ofile %||% "."), "lib_persample_cusum.R"),
    "lib_persample_cusum.R",
    "../v5_drafts/lib_persample_cusum.R",
    file.path(getwd(), "lib_persample_cusum.R")
  )
  hit <- guesses[file.exists(guesses)]
  if (length(hit) == 0L) stop("[TR_F] cannot locate lib_persample_cusum.R — pass --lib")
  LIB_PATH <- hit[1]
}
source(LIB_PATH)
message("[TR_F] cusum kernel: ", LIB_PATH)

per_chrom_dir <- file.path(OUTROOT, "03_per_chrom")
chroms <- if (!is.null(CHROM)) CHROM else list.files(per_chrom_dir, pattern = "^C_gar_LG[0-9]+$")
band_filter <- if (!is.na(BANDS_RUN)) strsplit(BANDS_RUN, ",")[[1]] else NULL

classify_spread <- function(iqr_bp) {
  if (!is.finite(iqr_bp)) return(NA_character_)
  if      (iqr_bp <  100000) "tight"
  else if (iqr_bp <  500000) "intermediate"
  else                       "ragged"
}

for (chrom in chroms) {
  rds <- file.path(per_chrom_dir, chrom, "precomp.rds")
  l2f <- file.path(per_chrom_dir, chrom, "L2_envelopes.tsv")
  caf <- file.path(per_chrom_dir, chrom, "carrier_assignments.tsv")
  if (!all(file.exists(c(rds, l2f, caf)))) {
    message("[TR_F] ", chrom, ": missing inputs — skip"); next
  }
  precomp <- readRDS(rds); dt <- precomp$dt
  l2_dt <- fread(l2f); ca_dt <- fread(caf)
  if (nrow(l2_dt) == 0L || nrow(ca_dt) == 0L) {
    fwrite(data.table(), file.path(per_chrom_dir, chrom, "cusum_per_sample.tsv.gz"),
           sep = "\t", compress = "gzip")
    fwrite(data.table(), file.path(per_chrom_dir, chrom, "cusum_boundary_dist.tsv"),
           sep = "\t")
    next
  }

  # Reconstruct the theta_pi matrix from PC_1_<sample> columns? No —
  # we need the actual θπ values. Read them from the long TSV.
  THETA_TSV_DIR <- Sys.getenv("THETA_TSV_DIR", unset = NA)
  PESTPG_SCALE  <- Sys.getenv("PESTPG_SCALE",  unset = "win10000.step2000")
  tsv <- file.path(THETA_TSV_DIR, sprintf("theta_native.%s.%s.tsv.gz", chrom, PESTPG_SCALE))
  if (!file.exists(tsv)) { message("[TR_F] ", chrom, ": missing TSV — skip"); next }
  long_dt <- fread(tsv)
  long_dt <- long_dt[chrom == ..chrom]
  sample_order <- precomp$sample_order
  win_idx <- dt$window_idx
  samp_to_row <- setNames(seq_along(sample_order), sample_order)
  win_to_col  <- setNames(seq_along(win_idx), as.character(win_idx))
  theta_mat <- matrix(NA_real_, nrow = length(sample_order), ncol = length(win_idx),
                      dimnames = list(sample_order, NULL))
  rows <- samp_to_row[long_dt$sample]
  cols <- win_to_col[as.character(long_dt$window_idx)]
  good <- !is.na(rows) & !is.na(cols)
  theta_mat[cbind(rows[good], cols[good])] <- long_dt$theta_pi[good]

  win_pos_bp <- dt$mid_bp

  per_sample_rows <- list()
  boundary_rows   <- list()

  for (k_row in seq_len(nrow(l2_dt))) {
    cid <- l2_dt$candidate_id[k_row]
    w_lo <- l2_dt$win_start[k_row]; w_hi <- l2_dt$win_end[k_row]
    if (w_hi - w_lo + 1L < 5L) next
    win_range_idx <- w_lo:w_hi
    win_pos_sub   <- win_pos_bp[win_range_idx]
    c_start <- l2_dt$start_bp[k_row]; c_end <- l2_dt$end_bp[k_row]

    ca_sub <- ca_dt[candidate_id == cid]
    if (nrow(ca_sub) == 0L) next
    bands_here <- if (is.null(band_filter)) unique(ca_sub$band) else
      intersect(band_filter, unique(ca_sub$band))

    for (b in bands_here) {
      members <- ca_sub[band == b, sample_id]
      if (length(members) < 5L) next
      M_sub <- theta_mat[members, win_range_idx, drop = FALSE]

      # ── (1) Per-sample CUSUM (existing) ─────────────────────────────────
      # Distribution: each carrier's individual breakpoint. Tells you whether
      # the boundary is tight, ragged, or bimodal across carriers.
      ps <- persample_cusum(M_sub, win_pos_sub, candidate_id = cid)
      ps[, `:=`(
        chrom              = chrom,
        band               = b,
        candidate_start_bp = c_start,
        candidate_end_bp   = c_end,
        dist_to_left_kb    = round(abs(cp_bp - c_start) / 1000, 2),
        dist_to_right_kb   = round(abs(cp_bp - c_end)   / 1000, 2),
        cp_side_inferred   = fifelse(is.na(cp_bp), NA_character_,
                                fifelse(abs(cp_bp - c_start) <= abs(cp_bp - c_end),
                                        "5_prime", "3_prime")),
        stream             = "theta"
      )]
      per_sample_rows[[length(per_sample_rows) + 1L]] <- ps

      # ── (2) Band-mean CUSUM (consensus boundary, sharper) ───────────────
      # Pool the band's samples by averaging θπ at each window, then run
      # CUSUM on the single 1D mean trace. Gives ONE breakpoint per band
      # per candidate, much sharper than per-sample because pooling reduces
      # noise by sqrt(n_band).
      mean_trace <- colMeans(M_sub, na.rm = TRUE)
      mean_mat   <- matrix(mean_trace, nrow = 1, dimnames = list(paste0("band:", b), NULL))
      consensus  <- persample_cusum(mean_mat, win_pos_sub, candidate_id = cid)
      consensus_cp_bp       <- as.integer(consensus$cp_bp[1])
      consensus_strength    <- round(consensus$strength[1], 3)
      consensus_asymmetry   <- as.integer(consensus$asymmetry[1])
      consensus_informative <- as.logical(consensus$informative[1])
      consensus_side <- if (is.na(consensus_cp_bp)) NA_character_
        else if (abs(consensus_cp_bp - c_start) <= abs(consensus_cp_bp - c_end)) "5_prime"
        else "3_prime"

      # Per-side boundary distribution: per-sample stats + consensus column.
      inf <- ps[informative == TRUE & is.finite(cp_bp)]
      for (side in c("5_prime", "3_prime")) {
        side_inf <- inf[cp_side_inferred == side]
        nc <- nrow(side_inf)
        # Did the band-mean CUSUM land on this side?
        cons_on_side <- !is.na(consensus_side) && consensus_side == side
        boundary_rows[[length(boundary_rows) + 1L]] <- data.table(
          candidate_id = cid, chrom = chrom, band = b, side = side,
          # Per-sample distribution (the manuscript-grade carrier spread):
          n_carriers   = nc,
          median_bp    = if (nc > 0L) as.integer(round(median(side_inf$cp_bp))) else NA_integer_,
          iqr_kb       = if (nc > 0L) round(IQR(side_inf$cp_bp) / 1000, 2) else NA_real_,
          spread_class = if (nc > 0L) classify_spread(IQR(side_inf$cp_bp)) else NA_character_,
          peak_strength = if (nc > 0L) round(max(side_inf$strength, na.rm = TRUE), 3) else NA_real_,
          # Band-mean consensus (the sharp boundary):
          consensus_cp_bp        = if (cons_on_side) consensus_cp_bp else NA_integer_,
          consensus_strength     = if (cons_on_side) consensus_strength else NA_real_,
          consensus_asymmetry    = if (cons_on_side) consensus_asymmetry else NA_integer_,
          consensus_informative  = if (cons_on_side) consensus_informative else NA,
          n_band_members         = length(members)
        )
      }
    }
  }

  ps_dt <- if (length(per_sample_rows) > 0) rbindlist(per_sample_rows, fill = TRUE) else data.table()
  bd_dt <- if (length(boundary_rows)   > 0) rbindlist(boundary_rows,   fill = TRUE) else data.table()
  fwrite(ps_dt, file.path(per_chrom_dir, chrom, "cusum_per_sample.tsv.gz"),
         sep = "\t", compress = "gzip")
  fwrite(bd_dt, file.path(per_chrom_dir, chrom, "cusum_boundary_dist.tsv"), sep = "\t")
  message(sprintf("[TR_F] %s: %d cusum rows across %d candidates × bands",
                  chrom, nrow(ps_dt), length(unique(ps_dt[, paste(candidate_id, band)]))))
}

message("[TR_F] DONE")

# Trivial null-coalesce so the lib-locate fallback compiles even when sys.frame
# is empty at script-execution time.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a