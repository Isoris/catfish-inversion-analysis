#!/usr/bin/env Rscript

# =============================================================================
# inspect_ghsl_v6.R — sanity check on STEP_C04_snake3_ghsl_v6.R output
#
# Just reads <chr>.ghsl_v6_matrices.rds and tells you whether the GHSL run
# produced something sensible. No reclassification. No new architecture.
#
# Usage:
#   Rscript inspect_ghsl_v6.R <ghsl_v6_matrices.rds> [out_dir]
#
# What it prints:
#   - Matrix dimensions (samples × windows)
#   - Number / fraction of NA cells per scale
#   - Per-sample mean divergence range (low ≈ all-INV-likely; high ≈ recurrent het)
#   - Per-window mean divergence percentile breaks
#   - Coverage by chromosome position (count of non-NA windows)
#
# What it writes (optional, if [out_dir] given):
#   <out_dir>/<chr>_ghsl_summary.tsv     per-sample summary
#   <out_dir>/<chr>_ghsl_per_window.tsv  per-window mean / sd / n_samples
#   <out_dir>/<chr>_ghsl_heatmap_s30.png (if pheatmap available)
#
# This is read-only inspection — does not modify the RDS or write new
# matrices. Safe to run repeatedly.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript inspect_ghsl_v6.R <ghsl_v6_matrices.rds> [out_dir]")
}
rds_in  <- args[1]
out_dir <- if (length(args) >= 2) args[2] else NA_character_

if (!file.exists(rds_in)) stop("Input RDS not found: ", rds_in)

cat("[INSPECT] Loading ", rds_in, " ...\n", sep = "")
ghsl <- readRDS(rds_in)

cat("[INSPECT] Top-level fields: ", paste(names(ghsl), collapse = ", "), "\n", sep = "")

chrom        <- if (!is.null(ghsl$chrom))        ghsl$chrom        else "?"
sample_names <- if (!is.null(ghsl$sample_names)) ghsl$sample_names else character()
window_info  <- if (!is.null(ghsl$window_info))  as.data.table(ghsl$window_info) else data.table()
params       <- if (!is.null(ghsl$params))       ghsl$params       else list()

cat("[INSPECT] Chromosome: ", chrom, "\n", sep = "")
cat("[INSPECT] N samples:  ", length(sample_names), "\n", sep = "")
cat("[INSPECT] N windows:  ", nrow(window_info), "\n", sep = "")
if (length(params) > 0) {
  cat("[INSPECT] Params:\n")
  for (nm in names(params)) {
    val <- params[[nm]]
    val_str <- if (length(val) > 5) paste0(paste(head(val, 5), collapse = ","), ",...") else paste(val, collapse = ",")
    cat("           ", nm, " = ", val_str, "\n", sep = "")
  }
}

# Find all matrices: raw div_mat, het_mat, and rolling at each scale
mat_names <- grep("^(div_mat|het_mat|rolling_)", names(ghsl), value = TRUE)
if (length(mat_names) == 0) {
  cat("[INSPECT] WARNING: no div_mat / rolling_* fields found.\n")
  cat("[INSPECT] Available top-level fields are listed above.\n")
  quit(status = 1)
}

cat("\n[INSPECT] === Matrix-level summary ===\n")
summary_rows <- list()
for (mn in mat_names) {
  M <- ghsl[[mn]]
  if (!is.matrix(M) && !is.data.frame(M)) next
  if (is.data.frame(M)) M <- as.matrix(M)
  n_cells <- length(M)
  n_na    <- sum(is.na(M))
  frac_na <- round(n_na / n_cells, 4)
  vals <- M[is.finite(M)]
  q     <- if (length(vals) > 0) quantile(vals, c(0, 0.05, 0.5, 0.95, 1), na.rm = TRUE) else rep(NA_real_, 5)
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    matrix      = mn,
    nrow        = nrow(M),
    ncol        = ncol(M),
    frac_na     = frac_na,
    q05         = round(q[2], 4),
    median      = round(q[3], 4),
    q95         = round(q[4], 4)
  )
  cat(sprintf("  %-20s  %4d x %5d  NA=%.3f  q05=%.3f  med=%.3f  q95=%.3f\n",
              mn, nrow(M), ncol(M), frac_na,
              q[2], q[3], q[4]))
}
summary_dt <- rbindlist(summary_rows, fill = TRUE)

# Pick the most useful matrix for per-sample / per-window summary.
# Prefer rolling at s30 if present (the chat-14 default sub-inversion scale),
# else rolling_50, else rolling_20, else raw div_mat.
pick_matrix <- function(ghsl, prefs) {
  for (p in prefs) {
    nm <- paste0("rolling_", p)
    if (!is.null(ghsl[[nm]])) return(list(name = nm, M = ghsl[[nm]]))
  }
  if (!is.null(ghsl$div_mat)) return(list(name = "div_mat", M = ghsl$div_mat))
  return(NULL)
}
picked <- pick_matrix(ghsl, c(30, 50, 20, 40, 100, 10))
if (is.null(picked)) {
  cat("[INSPECT] WARNING: could not pick a representative matrix.\n")
  quit(status = 1)
}
cat("\n[INSPECT] Using ", picked$name, " for per-sample / per-window summary.\n", sep = "")

M <- picked$M
if (is.data.frame(M)) M <- as.matrix(M)

# Per-sample summary
sample_dt <- NULL
if (nrow(M) == length(sample_names)) {
  sample_dt <- data.table(
    sample      = sample_names,
    n_windows   = rowSums(is.finite(M)),
    mean_div    = round(rowMeans(M, na.rm = TRUE), 4),
    median_div  = round(apply(M, 1, function(x) median(x, na.rm = TRUE)), 4),
    sd_div      = round(apply(M, 1, function(x) sd(x, na.rm = TRUE)), 4)
  )
  cat("\n[INSPECT] Per-sample mean divergence (top 5 lowest, top 5 highest):\n")
  print(head(sample_dt[order(mean_div)], 5))
  cat("\n  ...\n\n")
  print(tail(sample_dt[order(mean_div)], 5))
}

# Per-window summary
window_dt <- NULL
if (ncol(M) == nrow(window_info) && nrow(window_info) > 0) {
  window_dt <- data.table(
    window_idx  = seq_len(ncol(M)),
    start_bp    = if ("start_bp"  %in% names(window_info)) window_info$start_bp  else NA_integer_,
    end_bp      = if ("end_bp"    %in% names(window_info)) window_info$end_bp    else NA_integer_,
    n_samples_scored = colSums(is.finite(M)),
    mean_div    = round(colMeans(M, na.rm = TRUE), 4),
    sd_div      = round(apply(M, 2, function(x) sd(x, na.rm = TRUE)), 4)
  )
  cat("\n[INSPECT] Per-window divergence quantiles (across windows):\n")
  q <- quantile(window_dt$mean_div, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)
  for (i in seq_along(q)) cat(sprintf("  %s : %.4f\n", names(q)[i], q[i]))

  if ("start_bp" %in% names(window_info)) {
    chrom_span <- range(c(window_info$start_bp, window_info$end_bp), na.rm = TRUE)
    cat(sprintf("\n[INSPECT] Chromosome span covered: %.1f – %.1f Mb\n",
                chrom_span[1] / 1e6, chrom_span[2] / 1e6))
  }
}

# Optional file outputs
if (!is.na(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(summary_dt, file.path(out_dir, paste0(chrom, "_ghsl_summary.tsv")), sep = "\t")
  if (!is.null(sample_dt)) {
    fwrite(sample_dt, file.path(out_dir, paste0(chrom, "_ghsl_per_sample.tsv")), sep = "\t")
  }
  if (!is.null(window_dt)) {
    fwrite(window_dt, file.path(out_dir, paste0(chrom, "_ghsl_per_window.tsv")), sep = "\t")
  }

  # Heatmap if pheatmap available
  if (requireNamespace("pheatmap", quietly = TRUE) && nrow(M) > 1 && ncol(M) > 1) {
    png_path <- file.path(out_dir, paste0(chrom, "_ghsl_heatmap_", picked$name, ".png"))
    Mv <- M
    Mv[is.na(Mv)] <- median(Mv, na.rm = TRUE)   # fill NAs with median for plotting only
    # If too many windows, downsample columns for the heatmap (keep file readable)
    if (ncol(Mv) > 1500) {
      step <- ceiling(ncol(Mv) / 1500)
      idx  <- seq(1L, ncol(Mv), by = step)
      Mv   <- Mv[, idx, drop = FALSE]
      cat("[INSPECT] Heatmap downsampled columns by factor ", step, "\n", sep = "")
    }
    rownames(Mv) <- sample_names[seq_len(nrow(Mv))]
    tryCatch({
      pheatmap::pheatmap(
        Mv,
        cluster_rows = TRUE, cluster_cols = FALSE,
        show_rownames = FALSE, show_colnames = FALSE,
        main = paste0(chrom, " — GHSL ", picked$name, " (rows = samples, cols = windows L→R)"),
        filename = png_path, width = 14, height = 8
      )
      cat("[INSPECT] Wrote heatmap: ", png_path, "\n", sep = "")
    }, error = function(e) {
      cat("[INSPECT] Heatmap failed: ", conditionMessage(e), "\n", sep = "")
    })
  } else {
    cat("[INSPECT] pheatmap not installed; skipping heatmap.\n")
    cat("           install.packages('pheatmap')   # if you want it next time\n")
  }

  cat("\n[INSPECT] Wrote summary files to ", out_dir, "\n", sep = "")
}

cat("\n[INSPECT] Done.\n")
