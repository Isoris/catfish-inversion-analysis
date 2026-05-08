#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_F_cusum_per_carrier.R   (drafts/v5)
# =============================================================================
# Per-carrier CUSUM inside each L2 candidate. Two outputs per (candidate ×
# band):
#
#   (1) PER-SAMPLE CUSUM — one breakpoint per carrier (the distribution).
#       Manuscript value: "carrier 1 breaks at 18.94 Mb, carrier 17 outlier
#       at 17.20 Mb, IQR = 70 kb tight 5' / 1.2 Mb ragged 3'."
#
#   (2) BAND-MEAN CUSUM — pool the band's samples into one mean θπ trace,
#       CUSUM the trace → ONE consensus breakpoint per band per side.
#       Sharper because pooling reduces noise by sqrt(n_band). Manuscript
#       value: "the HIGH_DIV band's 3' boundary is at 18.95 Mb."
#
# Both are stacked in cusum_boundary_dist.tsv.
#
# Reads:   <PRECOMP_DIR>/<chr>.precomp.rds
#          <L2_DIR>/<chr>.L2_envelopes.tsv
#          <CARRIERS_DIR>/<chr>.carrier_assignments.tsv
#          $THETA_TSV_DIR/theta_native.<chr>.<scale>.tsv.gz   (re-read for θπ values)
#          v5_drafts/lib_persample_cusum.R                     (CUSUM kernel)
# Writes:  <OUT_DIR>/<chr>.cusum_per_sample.tsv.gz
#          <OUT_DIR>/<chr>.cusum_boundary_dist.tsv
#
# Defaults (z-blocks-shaped layout from TR_B v5):
#   PRECOMP_DIR  = $OUTROOT/precomp
#   L2_DIR       = $OUTROOT/L2_detect
#   CARRIERS_DIR = $OUTROOT/carriers
#   OUT_DIR      = $OUTROOT/cusum
#
# Usage:
#   Rscript STEP_TR_F_cusum_per_carrier.R --chr <CHR>
#                                          [--lib <path/to/lib_persample_cusum.R>]
#                                          [--bands "MID_DIV,HIGH_DIV"]
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}

CHR          <- get_arg("--chr")
PRECOMP_DIR  <- get_arg("--precomp_dir")
L2_DIR       <- get_arg("--l2_dir")
CARRIERS_DIR <- get_arg("--carriers_dir")
OUT_DIR      <- get_arg("--out_dir")
LIB_PATH     <- get_arg("--lib")
BANDS_RUN    <- get_arg("--bands")

OUTROOT       <- Sys.getenv("OUTROOT",       unset = NA)
THETA_TSV_DIR <- Sys.getenv("THETA_TSV_DIR", unset = NA)
PESTPG_SCALE  <- Sys.getenv("PESTPG_SCALE",  unset = "win10000.step2000")
stopifnot(!is.na(OUTROOT), !is.na(THETA_TSV_DIR))

if (is.na(PRECOMP_DIR))  PRECOMP_DIR  <- file.path(OUTROOT, "precomp")
if (is.na(L2_DIR))       L2_DIR       <- file.path(OUTROOT, "L2_detect")
if (is.na(CARRIERS_DIR)) CARRIERS_DIR <- file.path(OUTROOT, "carriers")
if (is.na(OUT_DIR))      OUT_DIR      <- file.path(OUTROOT, "cusum")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Locate the CUSUM kernel.
if (is.na(LIB_PATH)) {
  guesses <- c(
    "lib_persample_cusum.R",
    "v5_drafts/lib_persample_cusum.R",
    file.path(getwd(), "lib_persample_cusum.R"),
    file.path(dirname(get_arg("--file=", "")), "lib_persample_cusum.R")
  )
  hit <- guesses[file.exists(guesses)]
  if (length(hit) == 0L) stop("[TR_F] cannot locate lib_persample_cusum.R — pass --lib")
  LIB_PATH <- hit[1]
}
source(LIB_PATH)
message("[TR_F] CUSUM kernel: ", LIB_PATH)

chroms <- if (!is.na(CHR)) CHR else
  sub("\\.precomp\\.rds$", "",
      list.files(PRECOMP_DIR, pattern = "\\.precomp\\.rds$"))

band_filter <- if (!is.na(BANDS_RUN)) strsplit(BANDS_RUN, ",")[[1]] else NULL

classify_spread <- function(iqr_bp) {
  if (!is.finite(iqr_bp)) return(NA_character_)
  if      (iqr_bp <  100000) "tight"
  else if (iqr_bp <  500000) "intermediate"
  else                       "ragged"
}

windows_in_range <- function(dt, l2_row) {
  if ("win_start" %in% names(l2_row) && is.finite(l2_row$win_start)) {
    seq.int(as.integer(l2_row$win_start), as.integer(l2_row$win_end))
  } else {
    which(dt$end_bp >= l2_row$start_bp & dt$start_bp <= l2_row$end_bp)
  }
}

for (chrom in chroms) {
  rds <- file.path(PRECOMP_DIR,  paste0(chrom, ".precomp.rds"))
  l2f <- file.path(L2_DIR,       paste0(chrom, ".L2_envelopes.tsv"))
  caf <- file.path(CARRIERS_DIR, paste0(chrom, ".carrier_assignments.tsv"))
  if (!all(file.exists(c(rds, l2f, caf)))) {
    message("[TR_F] ", chrom, ": missing inputs — skip"); next
  }
  precomp <- readRDS(rds); dt <- precomp$dt
  l2_dt <- fread(l2f); ca_dt <- fread(caf)
  if (nrow(l2_dt) == 0L || nrow(ca_dt) == 0L) {
    fwrite(data.table(), file.path(OUT_DIR, paste0(chrom, ".cusum_per_sample.tsv.gz")),
           sep = "\t", compress = "gzip")
    fwrite(data.table(), file.path(OUT_DIR, paste0(chrom, ".cusum_boundary_dist.tsv")),
           sep = "\t")
    next
  }
  if (!"candidate_id" %in% names(l2_dt)) {
    l2_dt[, candidate_id := paste0(chrom, "_C_", sprintf("%03d", seq_len(.N)))]
  }

  # Reload θπ values from theta_native TSV (matrix isn't on the precomp.rds).
  tsv <- file.path(THETA_TSV_DIR,
                   sprintf("theta_native.%s.%s.tsv.gz", chrom, PESTPG_SCALE))
  if (!file.exists(tsv)) { message("[TR_F] ", chrom, ": no TSV — skip"); next }
  long_dt <- fread(tsv); long_dt <- long_dt[chrom == ..chrom]
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
    in_win <- windows_in_range(dt, l2_dt[k_row])
    if (length(in_win) < 5L) next
    win_pos_sub <- win_pos_bp[in_win]
    c_start <- as.integer(min(dt$start_bp[in_win]))
    c_end   <- as.integer(max(dt$end_bp[in_win]))

    ca_sub <- ca_dt[candidate_id == cid]
    if (nrow(ca_sub) == 0L) next
    bands_here <- if (is.null(band_filter)) unique(ca_sub$band) else
      intersect(band_filter, unique(ca_sub$band))

    for (b in bands_here) {
      members <- ca_sub[band == b, sample_id]
      if (length(members) < 5L) next
      M_sub <- theta_mat[members, in_win, drop = FALSE]

      # (1) Per-sample CUSUM — distribution of carrier breakpoints.
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

      # (2) Band-mean CUSUM — pooled trace, ONE breakpoint, sharp consensus.
      mean_trace <- colMeans(M_sub, na.rm = TRUE)
      mean_mat   <- matrix(mean_trace, nrow = 1,
                           dimnames = list(paste0("band:", b), NULL))
      consensus <- persample_cusum(mean_mat, win_pos_sub, candidate_id = cid)
      cons_cp_bp     <- as.integer(consensus$cp_bp[1])
      cons_strength  <- round(consensus$strength[1], 3)
      cons_asym      <- as.integer(consensus$asymmetry[1])
      cons_inform    <- as.logical(consensus$informative[1])
      cons_side <- if (is.na(cons_cp_bp)) NA_character_
        else if (abs(cons_cp_bp - c_start) <= abs(cons_cp_bp - c_end)) "5_prime"
        else "3_prime"

      inf <- ps[informative == TRUE & is.finite(cp_bp)]
      for (side in c("5_prime", "3_prime")) {
        side_inf <- inf[cp_side_inferred == side]
        nc <- nrow(side_inf)
        cons_on_side <- !is.na(cons_side) && cons_side == side
        boundary_rows[[length(boundary_rows) + 1L]] <- data.table(
          candidate_id   = cid, chrom = chrom, band = b, side = side,
          n_carriers     = nc,
          median_bp      = if (nc > 0L) as.integer(round(median(side_inf$cp_bp))) else NA_integer_,
          iqr_kb         = if (nc > 0L) round(IQR(side_inf$cp_bp) / 1000, 2) else NA_real_,
          spread_class   = if (nc > 0L) classify_spread(IQR(side_inf$cp_bp)) else NA_character_,
          peak_strength  = if (nc > 0L) round(max(side_inf$strength, na.rm = TRUE), 3) else NA_real_,
          consensus_cp_bp        = if (cons_on_side) cons_cp_bp else NA_integer_,
          consensus_strength     = if (cons_on_side) cons_strength else NA_real_,
          consensus_asymmetry    = if (cons_on_side) cons_asym else NA_integer_,
          consensus_informative  = if (cons_on_side) cons_inform else NA,
          n_band_members         = length(members)
        )
      }
    }
  }

  ps_dt <- if (length(per_sample_rows) > 0) rbindlist(per_sample_rows, fill = TRUE) else data.table()
  bd_dt <- if (length(boundary_rows)   > 0) rbindlist(boundary_rows,   fill = TRUE) else data.table()
  fwrite(ps_dt, file.path(OUT_DIR, paste0(chrom, ".cusum_per_sample.tsv.gz")),
         sep = "\t", compress = "gzip")
  fwrite(bd_dt, file.path(OUT_DIR, paste0(chrom, ".cusum_boundary_dist.tsv")), sep = "\t")
  message(sprintf("[TR_F] %s: %d cusum rows across %d candidates",
                  chrom, nrow(ps_dt), length(unique(ps_dt$candidate_id))))
}

message("[TR_F] DONE")