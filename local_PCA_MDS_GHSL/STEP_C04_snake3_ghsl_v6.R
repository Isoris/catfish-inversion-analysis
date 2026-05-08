#!/usr/bin/env Rscript

# =============================================================================
# STEP_C04_snake3_ghsl_v6.R
#
# SNAKE 3 v6: HEAVY ENGINE — Compute + Save Divergence Matrices
#
# PURPOSE:
#   Run ONCE per chromosome. Loads 77M+ variants, computes per-sample
#   within-haplotype divergence at raw window resolution, then applies
#   rolling means at multiple configurable scales. Saves everything as
#   RDS for the light classifier (STEP_C04b) to iterate on in 30 seconds.
#
# WHAT IT COMPUTES:
#   1. Raw div_mat [N_samples × N_windows]: ghsl_div = n_phased_het / n_total
#   2. Raw het_mat [N_samples × N_windows]: het_div = n_all_het / n_total
#   3. Rolling means at configurable scales (default: 20, 50, 100 windows)
#   4. Metadata: n_sites_mat, n_phased_het_mat, window coords, sample names
#
# WHAT IT DOES NOT DO:
#   No scoring. No karyotype calling. No PASS/FAIL. That's STEP_C04b.
#
# OUTPUT:
#   <outdir>/<chr>.ghsl_v6_matrices.rds  — one RDS per chromosome containing:
#     $div_mat          — raw GHSL divergence [samples × windows]
#     $het_mat          — raw het rate [samples × windows]
#     $n_sites_mat      — total variants per sample × window
#     $n_phased_het_mat — phased het count per sample × window
#     $rolling           — named list of rolling div_mats by scale
#     $rolling_het       — named list of rolling het_mats by scale
#     $window_info      — data.table with window coords (start_bp, end_bp, etc.)
#     $sample_names     — character vector of sample IDs
#     $chrom            — chromosome name
#     $params           — list of parameters used
#
# Usage:
#   Rscript STEP_C04_snake3_ghsl_v6.R <precomp_dir> <ghsl_prep_dir> <outdir> \
#     [--chrom C_gar_LG01] [--scales 20,50,100] [--min_phased 3]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop(paste(
  "Usage: Rscript STEP_C04_snake3_ghsl_v6.R <precomp_dir> <ghsl_prep_dir> <outdir> [opts]",
  "  --chrom <chr>       Process single chromosome",
  "  --scales 20,50,100  Rolling window scales (comma-separated)",
  "  --min_phased 3      Min phased sites per sample per window",
  "  --qual_min 20       Min QUAL for variant filtering",
  "  --gq_min 10         Min GQ for variant filtering",
  sep = "\n"
))

precomp_dir   <- args[1]
ghsl_prep_dir <- args[2]
outdir        <- args[3]

# Defaults
chrom_filter     <- NULL
# Chat 14 (2026-04-18): expanded default scale ladder.
# The previous default (20, 50, 100) skipped the 100-200 kb range
# where sub-inversions and short recombinant tracts live. New ladder
# gives an evenly-spaced fine ladder for detail work (10..50) plus
# s100 for chromosome-overview plotting. All scales are computed from
# the same 5-kb base matrix so the added cost is small (a few
# frollmean passes).
#   s10  = ~50 kb rolling
#   s20  = ~100 kb
#   s30  = ~150 kb
#   s40  = ~200 kb
#   s50  = ~250 kb   (primary display/karyotype/score scale)
#   s100 = ~500 kb   (overview scale)
ROLLING_SCALES   <- c(10L, 20L, 30L, 40L, 50L, 100L)
MIN_PHASED_SITES <- 3L
QUAL_MIN         <- 20
GQ_MIN           <- 10

i <- 4L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom" && i < length(args)) {
    chrom_filter <- args[i + 1]; i <- i + 2L
  } else if (a == "--scales" && i < length(args)) {
    ROLLING_SCALES <- as.integer(strsplit(args[i + 1], ",")[[1]])
    i <- i + 2L
  } else if (a == "--min_phased" && i < length(args)) {
    MIN_PHASED_SITES <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--qual_min" && i < length(args)) {
    QUAL_MIN <- as.numeric(args[i + 1]); i <- i + 2L
  } else if (a == "--gq_min" && i < length(args)) {
    GQ_MIN <- as.numeric(args[i + 1]); i <- i + 2L
  } else { i <- i + 1L }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("================================================================")
message("[S3v6] Snake 3 v6: HEAVY ENGINE — Divergence Matrices")
message("================================================================")
message("[S3v6] Precomp:        ", precomp_dir)
message("[S3v6] GHSL prep:      ", ghsl_prep_dir)
message("[S3v6] Output:         ", outdir)
message("[S3v6] MIN_PHASED=", MIN_PHASED_SITES)
message("[S3v6] Rolling scales: ", paste(ROLLING_SCALES, collapse = ", "), " windows")

# =============================================================================
# LOAD PRECOMP (for window grid and sample names)
# =============================================================================

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("[S3v6] FATAL: No .precomp.rds files in: ", precomp_dir)
message("[S3v6] Found ", length(rds_files), " precomp RDS files")

# Chat 14 (2026-04-18) chunking fix:
# Previously the loop below eagerly readRDS'd every chromosome's
# precomp file into `precomp_list` BEFORE the --chrom filter was
# applied, wasting 1–3 GB RAM for single-chromosome runs. Now we
# build a filename→chrom index by peeking at each file's `chrom`
# field (cheap: readRDS + access + discard), apply --chrom filter
# to the index, and then readRDS the heavy matrices ONLY for the
# chromosomes we actually process, inside the main loop.

# Build index: filename -> chrom (without retaining matrix payloads)
precomp_index <- list()  # chrom -> file path
for (f in rds_files) {
  obj <- readRDS(f)
  chr_name <- obj$chrom
  precomp_index[[chr_name]] <- f
  rm(obj); invisible(gc(verbose = FALSE, full = FALSE))
}
chroms <- names(precomp_index)
if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)

# Get sample names by peeking at one file (first chrom we'll actually use)
precomp_sample_names <- NULL
if (length(chroms) > 0) {
  obj <- readRDS(precomp_index[[chroms[1]]])
  pc1_cols <- grep("^PC_1_", names(obj$dt), value = TRUE)
  if (length(pc1_cols) > 0) precomp_sample_names <- sub("^PC_1_", "", pc1_cols)
  rm(obj); invisible(gc(verbose = FALSE, full = FALSE))
}
if (is.null(precomp_sample_names)) stop("[S3v6] FATAL: No PC_1_ columns found")
n_samples <- length(precomp_sample_names)
sample_names <- precomp_sample_names

# Map Ind→CGA
if (grepl("^Ind[0-9]", sample_names[1])) {
  sf <- Sys.getenv("SAMPLES_IND", "")
  if (!nzchar(sf) || !file.exists(sf)) {
    base <- Sys.getenv("BASE", "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04")
    for (candidate in c(
      file.path(base, "het_roh/01_inputs_check/samples.ind"),
      file.path(base, "popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv")
    )) {
      if (file.exists(candidate)) { sf <- candidate; break }
    }
  }
  if (nzchar(sf) && file.exists(sf)) {
    real <- trimws(readLines(sf))
    real <- real[nzchar(real)]
    if (length(real) == n_samples) {
      sample_names <- real
      message("[S3v6] Sample mapping: ", precomp_sample_names[1], " -> ", sample_names[1],
              " (", n_samples, " samples)")
    }
  }
}

message("[S3v6] Samples: ", n_samples)
message("[S3v6] Chromosomes: ", length(chroms))

# =============================================================================
# STAGE 1: Compute per-sample within-haplotype divergence (from v5, verbatim)
# =============================================================================

compute_divergence_matrix <- function(ghsl_by_window, dt, n_win, available_samples) {
  n_samp <- length(available_samples)
  ghsl_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                     dimnames = list(available_samples, NULL))
  het_mat  <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                     dimnames = list(available_samples, NULL))
  n_sites_mat <- matrix(0L, nrow = n_samp, ncol = n_win,
                        dimnames = list(available_samples, NULL))
  n_phased_het_mat <- matrix(0L, nrow = n_samp, ncol = n_win,
                             dimnames = list(available_samples, NULL))

  for (wi in seq_len(n_win)) {
    win_dt <- ghsl_by_window[[as.character(wi)]]
    if (is.null(win_dt) || nrow(win_dt) == 0) next

    for (si in available_samples) {
      sv <- win_dt[sample_id == si]
      n_total <- nrow(sv)
      if (n_total < MIN_PHASED_SITES) next

      gt <- tolower(sv$gt_class)
      n_het_all <- sum(gt == "het")

      # Phased hets: phase_gt contains "|" (0|1 or 1|0)
      phased_het <- sv[gt == "het" & grepl("\\|", phase_gt)]
      n_phased_het <- nrow(phased_het)

      ghsl_mat[si, wi] <- n_phased_het / n_total
      het_mat[si, wi]  <- n_het_all / n_total
      n_sites_mat[si, wi] <- n_total
      n_phased_het_mat[si, wi] <- n_phased_het
    }

    if (wi %% 500 == 0) {
      n_scored <- sum(!is.na(ghsl_mat[, wi]))
      message("[S3v6]   Window ", wi, "/", n_win, " (", n_scored, " samples scored)")
    }
  }

  list(ghsl = ghsl_mat, het = het_mat, n_sites = n_sites_mat,
       n_phased_het = n_phased_het_mat)
}

# =============================================================================
# STAGE 2: Rolling aggregation at multiple scales
# =============================================================================

compute_rolling_matrices <- function(raw_mat, scales) {
  # For each scale K, compute per-sample rolling mean across K windows.
  # Uses data.table::frollmean for speed.
  # NA handling: na.rm = TRUE so partial windows at edges still produce values.
  #
  # Returns named list: rolling[["s20"]] = matrix same dims as raw_mat
  
  n_samp <- nrow(raw_mat)
  n_win  <- ncol(raw_mat)
  rolling <- list()

  for (K in scales) {
    label <- paste0("s", K)
    message("[S3v6]   Rolling mean: scale=", K, " windows (~", round(K * 5, 0), " kb)")

    roll_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                       dimnames = dimnames(raw_mat))

    for (si in seq_len(n_samp)) {
      row_vals <- raw_mat[si, ]
      # frollmean with align="center" for symmetric smoothing
      smoothed <- frollmean(row_vals, n = K, align = "center", na.rm = TRUE)
      roll_mat[si, ] <- smoothed
    }

    # Coverage check: how many cells went from NA to non-NA (rescue) or vice versa
    n_raw_valid  <- sum(!is.na(raw_mat))
    n_roll_valid <- sum(!is.na(roll_mat))
    message("[S3v6]     Coverage: raw=", round(100 * n_raw_valid / length(raw_mat), 1),
            "% → rolling=", round(100 * n_roll_valid / length(roll_mat), 1), "%")

    rolling[[label]] <- roll_mat
  }

  rolling
}

# =============================================================================
# MAIN LOOP
# =============================================================================

for (chr in chroms) {
  # Chat 14: load THIS chrom's precomp on demand (was: eager precomp_list)
  pc <- readRDS(precomp_index[[chr]])
  if (is.null(pc) || pc$n_windows < 20) {
    message("[S3v6] SKIP ", chr, ": only ", pc$n_windows %||% 0, " windows")
    rm(pc); invisible(gc(verbose = FALSE, full = FALSE))
    next
  }

  dt <- pc$dt
  n_win <- nrow(dt)
  message("\n================================================================")
  message("[S3v6] ======= ", chr, " (", n_win, " windows) =======")
  message("================================================================")

  # ── Load + pre-index ──
  ghsl_file <- file.path(ghsl_prep_dir, paste0(chr, ".merged_phased_snps.tsv.gz"))
  if (!file.exists(ghsl_file)) {
    message("[S3v6] SKIP ", chr, ": no merged phased SNPs file")
    next
  }

  t0 <- proc.time()
  message("[S3v6] Loading ", basename(ghsl_file), " ...")
  ghsl_dt <- fread(ghsl_file)
  message("[S3v6]   Raw: ", formatC(nrow(ghsl_dt), big.mark = ","), " variants in ",
          round((proc.time() - t0)[3], 1), "s")

  # QC filter
  if ("qual" %in% names(ghsl_dt)) ghsl_dt <- ghsl_dt[is.na(qual) | qual >= QUAL_MIN]
  if ("gq" %in% names(ghsl_dt)) ghsl_dt <- ghsl_dt[is.na(gq) | gq >= GQ_MIN]
  message("[S3v6]   After QC: ", formatC(nrow(ghsl_dt), big.mark = ","))

  available_samples <- intersect(unique(ghsl_dt$sample_id), sample_names)
  message("[S3v6]   Samples: ", length(available_samples), " / ", n_samples)

  if (length(available_samples) < 20) {
    message("[S3v6] SKIP ", chr, ": only ", length(available_samples), " samples with data")
    next
  }

  # Pre-index by window
  message("[S3v6]   Pre-indexing...")
  t_idx <- proc.time()
  ghsl_dt <- ghsl_dt[sample_id %in% available_samples]
  window_starts <- dt$start_bp
  window_ends   <- dt$end_bp
  ghsl_dt[, window_id := findInterval(pos, window_starts)]
  ghsl_dt[window_id > n_win, window_id := n_win]
  ghsl_dt[window_id < 1, window_id := 1L]
  ghsl_dt <- ghsl_dt[pos <= window_ends[window_id]]

  ghsl_by_window <- split(ghsl_dt, ghsl_dt$window_id)
  message("[S3v6]   Indexed: ", formatC(nrow(ghsl_dt), big.mark = ","),
          " variants in ", round((proc.time() - t_idx)[3], 1), "s")

  # ── STAGE 1: Raw divergence matrix ──
  message("[S3v6] Stage 1: Computing raw divergence matrix...")
  t1 <- proc.time()
  div_result <- compute_divergence_matrix(ghsl_by_window, dt, n_win, available_samples)
  elapsed1 <- round((proc.time() - t1)[3], 1)

  div_mat <- div_result$ghsl
  het_mat <- div_result$het

  n_scored <- sum(!is.na(div_mat))
  n_possible <- length(available_samples) * n_win
  message("[S3v6]   Raw divergence: ", formatC(n_scored, big.mark = ","),
          " / ", formatC(n_possible, big.mark = ","),
          " (", round(100 * n_scored / n_possible, 1), "%) in ", elapsed1, "s")

  # Quick stats
  ghsl_vals <- div_mat[!is.na(div_mat)]
  if (length(ghsl_vals) > 0) {
    message("[S3v6]   GHSL: min=", round(min(ghsl_vals), 4),
            " median=", round(median(ghsl_vals), 4),
            " max=", round(max(ghsl_vals), 4))
  }

  # ── STAGE 2: Rolling aggregation ──
  message("[S3v6] Stage 2: Rolling aggregation...")
  t2 <- proc.time()
  rolling_div <- compute_rolling_matrices(div_mat, ROLLING_SCALES)
  rolling_het <- compute_rolling_matrices(het_mat, ROLLING_SCALES)
  elapsed2 <- round((proc.time() - t2)[3], 1)
  message("[S3v6]   Rolling computed in ", elapsed2, "s")

  # ── STAGE 3: Save RDS ──
  message("[S3v6] Stage 3: Saving matrices...")

  window_info <- dt[, .(global_window_id, start_bp, end_bp)]
  window_info[, pos_mb := round((start_bp + end_bp) / 2e6, 4)]

  out_rds <- file.path(outdir, paste0(chr, ".ghsl_v6_matrices.rds"))
  saveRDS(list(
    div_mat          = div_mat,
    het_mat          = het_mat,
    n_sites_mat      = div_result$n_sites,
    n_phased_het_mat = div_result$n_phased_het,
    rolling          = rolling_div,
    rolling_het      = rolling_het,
    window_info      = window_info,
    sample_names     = available_samples,
    chrom            = chr,
    params           = list(
      min_phased_sites = MIN_PHASED_SITES,
      qual_min         = QUAL_MIN,
      gq_min           = GQ_MIN,
      rolling_scales   = ROLLING_SCALES
    )
  ), out_rds)

  fsize <- round(file.info(out_rds)$size / 1e6, 1)
  message("[S3v6]   Saved: ", out_rds, " (", fsize, " MB)")

  # Free memory (chat 14: include pc from on-demand readRDS)
  rm(ghsl_dt, ghsl_by_window, div_result, pc, dt)
  gc(verbose = FALSE)

  message("[S3v6] ", chr, " DONE in ", round((proc.time() - t0)[3], 1), "s total")
}

message("\n================================================================")
message("[DONE] Snake 3 v6: Heavy engine complete")
message("================================================================")
message("  Output: ", outdir, "/")
message("  Files:  <chr>.ghsl_v6_matrices.rds")
message("")
message("  Next: run STEP_C04b_snake3_ghsl_classify.R on these RDS files")
message("  That script loads in seconds — iterate on scoring/classification there.")
