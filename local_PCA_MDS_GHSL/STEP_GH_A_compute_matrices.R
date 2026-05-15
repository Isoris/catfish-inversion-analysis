#!/usr/bin/env Rscript
# =============================================================================
# STEP_GH_A_compute_matrices.R  (v2 — clean rewrite, 2026-05-12)
# =============================================================================
# GHSL stage A: per-sample × per-window haplotype-divergence matrices.
#
# Reads a Clair3 phased-SNP TSV for one chromosome, builds raw 5 kb (or
# whatever window grid the precomp dictates) divergence matrices, and
# produces multi-scale rolling aggregates for downstream local PCA / MDS
# (GH_B, GH_C) and the existing classifier (GH_B_classify).
#
# Definition (GHSL, Genome Heterozygosity by Sequence Length, after the
# Nature hybrid-potato 2024 paper):
#
#   For each sample s and window w:
#     div[s, w] = n_phased_het[s, w] / n_total[s, w]
#     het[s, w] = n_all_het[s, w]    / n_total[s, w]
#
#   Where:
#     - n_phased_het = HET sites with phase_gt matching "0|1" or "1|0"
#                      (positive evidence that hap1 ≠ hap2 within the
#                       sample at that site)
#     - n_all_het    = HET sites regardless of phase ("0|1", "1|0", "0/1")
#     - n_total      = all called variants for the sample in the window
#                      (phased hets + unphased hets + hom_var; hom_ref
#                       is filtered upstream by Clair3)
#
# Per-window denominator uses *all* variants (including hom_var) so the
# ratio is normalized by sequence length, as in the Nature paper. This
# matters for HOM samples: their dense hom_var contribution pushes the
# ratio toward 0, while HET carriers in unrecombining inversions push
# the ratio up. The "denominator confound" (carriers have systematically
# more variants inside inversions) is a known limitation, mitigated at
# downstream scales by rolling smoothing.
#
# v2 changes vs v6 / chat-14 paste:
#   1. SPEEDUP: compute_divergence_matrix is now a single data.table
#      group-by aggregation (no nested for(wi) for(si) loop). The "|"
#      match is grepl(..., fixed = TRUE), not regex. Expected: 100x+
#      end-to-end on a 5 kb-grid chromosome.
#   2. CORRECTNESS: rolling[["sK"]] is now the *ratio of sums*
#      (rolling_n_phased_het / rolling_n_total), not the mean of ratios.
#      Equivalent to weighting each raw window's ratio by its variant
#      count. Stable under sparse low-coverage windows.
#   3. NEW OUTPUTS: rolling_n_total[["sK"]] and rolling_n_phased_het[["sK"]]
#      stored alongside rolling[["sK"]]. Required by GH_B's local PCA
#      heteroscedastic weight (sqrt(n_tot / median(n_tot))) and useful
#      for QC overlays in the atlas.
#   4. CLEAN: sample list is read from --sample_list (one id per line),
#      not reverse-engineered from precomp.rds PC_1_* columns. Precomp
#      is read only for window grid (start_bp, end_bp).
#
# Output (one file per chromosome):
#   <outdir>/<chr>.ghsl_matrices.rds
#     $div_mat                  — [n_samp × n_win] raw n_phased_het/n_total
#     $het_mat                  — [n_samp × n_win] raw n_all_het/n_total
#     $n_total_mat              — [n_samp × n_win] integer denominator
#     $n_phased_het_mat         — [n_samp × n_win] integer numerator
#     $n_all_het_mat            — [n_samp × n_win] integer all-het count
#     $rolling                  — list of [n_samp × n_win] ratio-of-sums
#                                  matrices, one per scale ("s10","s50",…)
#     $rolling_het              — same shape, het-fraction
#     $rolling_n_total          — list of [n_samp × n_win] integer rolling
#                                  sums of n_total (per scale)
#     $rolling_n_phased_het     — same, rolling sums of n_phased_het
#     $window_info              — data.table: window_idx, start_bp, end_bp, mid_bp
#     $sample_names             — character [n_samp]
#     $chrom                    — character
#     $params                   — list of run parameters
#
# Usage:
#   Rscript STEP_GH_A_compute_matrices.R \
#     --precomp_dir <dir>            # window grid source (*.precomp.rds)
#     --ghsl_prep_dir <dir>          # <chr>.merged_phased_snps.tsv.gz inputs
#     --sample_list <file>           # one sample id per line
#     --outdir <dir>                 # output for <chr>.ghsl_matrices.rds
#     [--chrom <chr>]                # single chromosome (default: all in precomp_dir)
#     [--scales 10,50,100]           # rolling scales (default: 10,50,100)
#     [--min_total 3]                # min n_total per (sample, window) to score
#     [--qual_min 20]                # min QUAL filter on variants
#     [--gq_min 10]                  # min GQ filter on variants
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ── Argument parsing ─────────────────────────────────────────────────────────
PRECOMP_DIR   <- NA_character_
GHSL_PREP_DIR <- NA_character_
SAMPLE_LIST   <- NA_character_
OUTDIR        <- NA_character_
CHROM         <- NULL
SCALES        <- c(10L, 50L, 100L)
MIN_TOTAL     <- 3L
QUAL_MIN      <- 20
GQ_MIN        <- 10

args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--precomp_dir"   && i < length(args)) { PRECOMP_DIR   <- args[i + 1]; i <- i + 2L }
  else if (a == "--ghsl_prep_dir" && i < length(args)) { GHSL_PREP_DIR <- args[i + 1]; i <- i + 2L }
  else if (a == "--sample_list"   && i < length(args)) { SAMPLE_LIST   <- args[i + 1]; i <- i + 2L }
  else if (a == "--outdir"        && i < length(args)) { OUTDIR        <- args[i + 1]; i <- i + 2L }
  else if (a == "--chrom"         && i < length(args)) { CHROM         <- args[i + 1]; i <- i + 2L }
  else if (a == "--scales"        && i < length(args)) {
    SCALES <- as.integer(strsplit(args[i + 1], "[ ,]+")[[1]]); i <- i + 2L
  }
  else if (a == "--min_total"     && i < length(args)) { MIN_TOTAL <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--qual_min"      && i < length(args)) { QUAL_MIN  <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--gq_min"        && i < length(args)) { GQ_MIN    <- as.numeric(args[i + 1]); i <- i + 2L }
  else { i <- i + 1L }
}
stopifnot(!is.na(PRECOMP_DIR), !is.na(GHSL_PREP_DIR),
          !is.na(SAMPLE_LIST), !is.na(OUTDIR))
SCALES <- sort(unique(SCALES[SCALES > 0L]))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ── Sample list ──────────────────────────────────────────────────────────────
sample_names <- readLines(SAMPLE_LIST)
sample_names <- sample_names[nchar(sample_names) > 0]
n_samp <- length(sample_names)
message("[GH_A] cohort: ", n_samp, " samples")
message("[GH_A] scales: ", paste(SCALES, collapse = ","))
message("[GH_A] min_total=", MIN_TOTAL, ", qual_min=", QUAL_MIN, ", gq_min=", GQ_MIN)

# ── Chromosome list (from precomp filenames) ─────────────────────────────────
precomp_files <- sort(list.files(PRECOMP_DIR, pattern = "\\.precomp\\.rds$",
                                 full.names = TRUE))
if (length(precomp_files) == 0L) stop("[GH_A] no .precomp.rds files in ", PRECOMP_DIR)
chrom_from_file <- sub("\\.precomp\\.rds$", "", basename(precomp_files))
names(precomp_files) <- chrom_from_file
chroms <- if (is.null(CHROM)) chrom_from_file else intersect(chrom_from_file, CHROM)
if (length(chroms) == 0L) stop("[GH_A] --chrom ", CHROM, " not found in ", PRECOMP_DIR)
message("[GH_A] chroms: ", length(chroms))

# =============================================================================
# Core: vectorized per-window aggregation
# =============================================================================

# Inputs:
#   ghsl_dt    — data.table with columns: sample_id, pos, gt_class, phase_gt,
#                qual?, gq?  (filtered upstream)
#   win_grid   — data.table with columns: window_idx (0-based), start_bp, end_bp
#   n_samp     — length(sample_names)
#   sample_to_row — named integer: sample_id -> 1..n_samp
#
# Output: list of four matrices [n_samp × n_win]:
#   n_total_mat, n_phased_het_mat, n_all_het_mat, plus the derived
#   div_mat and het_mat. Missing (sample, window) cells stay NA (div, het)
#   or 0 (count matrices).
compute_divergence_matrices <- function(ghsl_dt, win_grid, n_samp, sample_to_row,
                                        min_total) {
  n_win <- nrow(win_grid)

  # Annotate per-variant booleans once (vectorized, no per-row R calls)
  ghsl_dt[, is_het    := tolower(gt_class) == "het"]
  ghsl_dt[, is_phased := grepl("|", phase_gt, fixed = TRUE)]
  # phased het = HET AND phase_gt has a literal "|"
  ghsl_dt[, is_phased_het := is_het & is_phased]

  # Bin each variant into a window by findInterval over start_bp.
  # win_grid is sorted by start_bp; off-grid variants (pos < first start or
  # pos > last end) become window_idx = NA and are dropped.
  starts <- win_grid$start_bp
  ends   <- win_grid$end_bp
  ghsl_dt[, wi := findInterval(pos, starts)]
  ghsl_dt <- ghsl_dt[wi >= 1L & wi <= n_win]
  # Discard variants past their assigned window's end_bp
  ghsl_dt <- ghsl_dt[pos <= ends[wi]]

  # Restrict to samples we know about
  ghsl_dt[, samp_row := sample_to_row[sample_id]]
  ghsl_dt <- ghsl_dt[!is.na(samp_row)]

  # Single group-by aggregation — replaces the nested for(wi) for(si) loop
  agg <- ghsl_dt[, .(
    n_total      = .N,
    n_all_het    = sum(is_het),
    n_phased_het = sum(is_phased_het)
  ), by = .(samp_row, wi)]

  # Pre-allocate output matrices
  n_total_mat      <- matrix(0L, nrow = n_samp, ncol = n_win,
                             dimnames = list(NULL, NULL))
  n_phased_het_mat <- matrix(0L, nrow = n_samp, ncol = n_win)
  n_all_het_mat    <- matrix(0L, nrow = n_samp, ncol = n_win)

  # Scatter the agg into matrices (one matrix-indexed assignment per matrix)
  rc <- cbind(agg$samp_row, agg$wi)
  n_total_mat[rc]      <- agg$n_total
  n_all_het_mat[rc]    <- agg$n_all_het
  n_phased_het_mat[rc] <- agg$n_phased_het

  # Derived ratios with min_total guard
  div_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win)
  het_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win)
  ok <- n_total_mat >= min_total
  div_mat[ok] <- n_phased_het_mat[ok] / n_total_mat[ok]
  het_mat[ok] <- n_all_het_mat[ok]    / n_total_mat[ok]

  list(div_mat          = div_mat,
       het_mat          = het_mat,
       n_total_mat      = n_total_mat,
       n_phased_het_mat = n_phased_het_mat,
       n_all_het_mat    = n_all_het_mat)
}

# =============================================================================
# Multi-scale rolling: ratio of rolling sums (NOT mean of ratios)
# =============================================================================

# For each scale K, returns four [n_samp × n_win] matrices:
#   rolling_n_total[[sK]]      — rolling SUM of n_total over K windows
#   rolling_n_phased_het[[sK]] — rolling SUM of n_phased_het
#   rolling_n_all_het[[sK]]    — rolling SUM of n_all_het
#   rolling_div[[sK]]          — rolling SUM(n_phased_het) / rolling SUM(n_total)
#   rolling_het[[sK]]          — rolling SUM(n_all_het)    / rolling SUM(n_total)
#
# data.table::frollsum operates column-by-column on a numeric vector. We
# transpose to (n_win × n_samp), apply frollsum to each column (which is
# fast C-level), then transpose back. This avoids the per-sample R loop
# that the v6 paste used and processes all 226 samples per scale in one
# data.table-style call.
compute_rolling <- function(n_total_mat, n_phased_het_mat, n_all_het_mat, scales) {
  n_samp <- nrow(n_total_mat); n_win <- ncol(n_total_mat)
  out_total   <- list()
  out_phased  <- list()
  out_allhet  <- list()
  out_div     <- list()
  out_het     <- list()
  # Transpose once — work column-major over windows
  tot_T <- t(n_total_mat)
  pha_T <- t(n_phased_het_mat)
  ahe_T <- t(n_all_het_mat)
  for (K in scales) {
    lab <- paste0("s", K)
    K_eff <- min(K, n_win)
    # frollsum: center-aligned, fill with NA at edges
    roll_tot <- vapply(seq_len(n_samp),
      function(s) frollsum(tot_T[, s], n = K_eff, align = "center", na.rm = TRUE),
      numeric(n_win))
    roll_pha <- vapply(seq_len(n_samp),
      function(s) frollsum(pha_T[, s], n = K_eff, align = "center", na.rm = TRUE),
      numeric(n_win))
    roll_ahe <- vapply(seq_len(n_samp),
      function(s) frollsum(ahe_T[, s], n = K_eff, align = "center", na.rm = TRUE),
      numeric(n_win))
    # Each frollsum returns length-n_win; vapply stacks → (n_win × n_samp).
    # Transpose back to (n_samp × n_win) for output.
    out_total[[lab]]  <- t(roll_tot)
    out_phased[[lab]] <- t(roll_pha)
    out_allhet[[lab]] <- t(roll_ahe)

    # Ratios — NA when denominator < K_eff (incomplete window) or = 0
    denom <- out_total[[lab]]
    div_K <- matrix(NA_real_, n_samp, n_win)
    het_K <- matrix(NA_real_, n_samp, n_win)
    ok <- is.finite(denom) & denom > 0
    div_K[ok] <- out_phased[[lab]][ok] / denom[ok]
    het_K[ok] <- out_allhet[[lab]][ok] / denom[ok]
    out_div[[lab]] <- div_K
    out_het[[lab]] <- het_K
    message(sprintf("[GH_A]   rolling s%d: median(n_tot per win)=%s, coverage=%.1f%%",
                    K,
                    formatC(median(denom[ok], na.rm = TRUE), big.mark = ",",
                            format = "d"),
                    100 * sum(ok) / length(denom)))
  }
  list(n_total = out_total, n_phased_het = out_phased, n_all_het = out_allhet,
       div = out_div, het = out_het)
}

# =============================================================================
# Main loop
# =============================================================================
sample_to_row <- setNames(seq_along(sample_names), sample_names)

for (chr in chroms) {
  t_chr <- proc.time()
  message("\n[GH_A] ===== ", chr, " =====")

  pc_file <- precomp_files[[chr]]
  pc <- readRDS(pc_file)
  if (is.null(pc) || is.null(pc$dt) || nrow(pc$dt) < 20L) {
    message("[GH_A] ", chr, ": precomp has <20 windows — skip"); next
  }
  # window_idx is a 0-based row-order index into the per-chrom window grid.
  # The slim Z-path precomp doesn't ship a column called window_idx
  # explicitly (it has window_index_chr, or just row order). Resolve to
  # whichever is available; fall back to seq_len()-1 if neither exists.
  pc_dt <- pc$dt
  if ("window_idx" %in% names(pc_dt)) {
    win_grid <- pc_dt[, .(window_idx, start_bp, end_bp)]
  } else if ("window_index_chr" %in% names(pc_dt)) {
    win_grid <- pc_dt[, .(window_idx = window_index_chr, start_bp, end_bp)]
  } else {
    win_grid <- pc_dt[, .(window_idx = seq_len(.N) - 1L, start_bp, end_bp)]
  }
  setorder(win_grid, start_bp)
  rm(pc, pc_dt); invisible(gc(verbose = FALSE))
  n_win <- nrow(win_grid)
  message("[GH_A] ", chr, ": ", n_win, " windows")

  ghsl_file <- file.path(GHSL_PREP_DIR, paste0(chr, ".merged_phased_snps.tsv.gz"))
  if (!file.exists(ghsl_file)) {
    message("[GH_A] ", chr, ": no ", basename(ghsl_file), " — skip"); next
  }

  t0 <- proc.time()
  ghsl_dt <- fread(ghsl_file)
  message(sprintf("[GH_A] %s: read %s rows in %.1fs",
                  chr,
                  formatC(nrow(ghsl_dt), big.mark = ","),
                  (proc.time() - t0)[3]))

  # QC filter (only on columns that exist in this file)
  if ("qual" %in% names(ghsl_dt)) ghsl_dt <- ghsl_dt[is.na(qual) | qual >= QUAL_MIN]
  if ("gq"   %in% names(ghsl_dt)) ghsl_dt <- ghsl_dt[is.na(gq)   | gq   >= GQ_MIN]

  # Aggregate
  t1 <- proc.time()
  agg <- compute_divergence_matrices(ghsl_dt, win_grid, n_samp, sample_to_row,
                                     MIN_TOTAL)
  message(sprintf("[GH_A] %s: divergence matrices (%d × %d) %.1fs",
                  chr, n_samp, n_win, (proc.time() - t1)[3]))
  rm(ghsl_dt); invisible(gc(verbose = FALSE))

  div_vals <- agg$div_mat[is.finite(agg$div_mat)]
  if (length(div_vals) > 0L) {
    message(sprintf("[GH_A] %s: div  min=%.4f med=%.4f max=%.4f n_scored=%s",
                    chr, min(div_vals), median(div_vals), max(div_vals),
                    formatC(length(div_vals), big.mark = ",")))
  }

  # Rolling (ratio of sums)
  t2 <- proc.time()
  roll <- compute_rolling(agg$n_total_mat, agg$n_phased_het_mat,
                          agg$n_all_het_mat, SCALES)
  message(sprintf("[GH_A] %s: rolling at %d scales %.1fs",
                  chr, length(SCALES), (proc.time() - t2)[3]))

  # Save
  window_info <- data.table(
    window_idx = win_grid$window_idx,
    start_bp   = win_grid$start_bp,
    end_bp     = win_grid$end_bp,
    mid_bp     = as.integer((win_grid$start_bp + win_grid$end_bp) / 2L)
  )
  out <- list(
    div_mat              = agg$div_mat,
    het_mat              = agg$het_mat,
    n_total_mat          = agg$n_total_mat,
    n_phased_het_mat     = agg$n_phased_het_mat,
    n_all_het_mat        = agg$n_all_het_mat,
    rolling              = roll$div,
    rolling_het          = roll$het,
    rolling_n_total      = roll$n_total,
    rolling_n_phased_het = roll$n_phased_het,
    rolling_n_all_het    = roll$n_all_het,
    window_info          = window_info,
    sample_names         = sample_names,
    chrom                = chr,
    params               = list(
      scales    = SCALES,
      min_total = MIN_TOTAL,
      qual_min  = QUAL_MIN,
      gq_min    = GQ_MIN,
      schema    = "GH_A v2 (ratio-of-sums rolling)"
    )
  )
  out_file <- file.path(OUTDIR, paste0(chr, ".ghsl_matrices.rds"))
  saveRDS(out, out_file)
  message(sprintf("[GH_A] %s: wrote %s (%.1f MB) in %.1fs total",
                  chr, out_file,
                  file.info(out_file)$size / 1e6,
                  (proc.time() - t_chr)[3]))
  rm(agg, roll, out); invisible(gc(verbose = FALSE))
}

message("[GH_A] DONE")
