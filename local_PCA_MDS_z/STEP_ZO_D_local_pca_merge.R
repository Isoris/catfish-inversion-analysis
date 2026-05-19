#!/usr/bin/env Rscript

# =============================================================================
# STEP_ZO_D_local_pca_merge.R  (v10, 2026-05-16, harmonized layout)
#
# STAGE 2 (merge) — coordination step that runs ONCE after all 01b array
# tasks complete. Walks every per-chrom *_tmp.rds in <outdir>/ (was
# <outdir>/tmp/ in v9), assigns globally unique sequential window_id
# values across the whole genome, and rewrites each per-chrom file with
# the patched IDs. Fast (~seconds).
#
# Pipeline position:
#   ZO_C local_pca_compute  ->  ZO_D LOCAL_PCA_MERGE  ->  ZO_E MDS  -> ...
#
# ── Inputs ────────────────────────────────────────────────────────────────
#   --outdir   <dir>   the same outdir ZO_C wrote to (01_local_pca/).
#                      Reads <chr>.window_pca_tmp.rds and *.chr_meta.tsv.gz
#                      directly from <outdir>/ (no tmp/ subdir in v10).
#   --registry_dir <dir>   optional; where windows_master.tsv.gz lands
#                          (default = $PATH1_REGISTRY, fallback = outdir).
#
# ── Outputs ───────────────────────────────────────────────────────────────
# In <registry_dir>/ (harmonized 02_dense_registry/ slot):
#   windows_master.tsv.gz          THE master window registry (genome-wide)
#   windows_master_summary.tsv     per-chromosome summary
# In <outdir>/ (harmonized 01_local_pca/ slot):
#   <chr>.window_pca.rds           per-chrom PCA, finalized with window_ids
#   <chr>.window_pca.tsv.gz        per-chrom PCA loadings, finalized
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: STEP09b_stage2_merge_registry.R v8.3-parallel)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# PARSE ARGS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
outdir       <- NULL
registry_dir <- NA_character_

i <- 1L
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) {
    outdir <- args[i + 1]; i <- i + 2L
  } else if (args[i] == "--registry_dir" && i < length(args)) {
    registry_dir <- args[i + 1]; i <- i + 2L
  } else {
    i <- i + 1L
  }
}

if (is.null(outdir)) stop("Usage: Rscript STEP_ZO_D_local_pca_merge.R --outdir <dir> [--registry_dir <dir>]")

# v10: stage-1 files now live directly in <outdir> (no tmp/ subdir).
# Back-compat: if no _tmp.rds at top level but a tmp/ subdir is present,
# fall back to reading from <outdir>/tmp/.
src_dir <- outdir
if (length(list.files(src_dir, pattern = "\\.window_pca_tmp\\.rds$")) == 0L) {
  legacy <- file.path(outdir, "tmp")
  if (dir.exists(legacy)) src_dir <- legacy
  else stop("No stage-1 outputs found in ", outdir, " (or its tmp/ subdir).")
}
tmpdir <- src_dir   # keep variable name for back-compat below

# Registry dir: --registry_dir wins, then env PATH1_REGISTRY, else outdir.
if (is.na(registry_dir)) {
  env_reg <- Sys.getenv("PATH1_REGISTRY", unset = "")
  registry_dir <- if (nzchar(env_reg)) env_reg else outdir
}
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)

message("[ZO_D] ═══════ Master registry merge (v10 harmonized layout) ═══════")
message("[ZO_D] Reading from:   ", src_dir)
message("[ZO_D] Writing per-chr finals to: ", outdir)
message("[ZO_D] Writing registry to:      ", registry_dir)

# =============================================================================
# DISCOVER + SORT CHROMOSOMES
# =============================================================================

meta_files <- sort(list.files(tmpdir, pattern = "\\.chr_meta\\.tsv\\.gz$", full.names = TRUE))
if (length(meta_files) == 0) stop("No .chr_meta.tsv.gz files found in: ", tmpdir)

chroms <- sub("\\.chr_meta\\.tsv\\.gz$", "", basename(meta_files))
message("[STEP09b-S2] Found ", length(chroms), " chromosomes: ",
        paste(chroms[1:min(3, length(chroms))], collapse = ", "),
        if (length(chroms) > 3) paste0(", ... (", length(chroms), " total)") else "")

# =============================================================================
# PASS 1: READ ALL METAS, ASSIGN GLOBAL WINDOW_ID
# =============================================================================

t0 <- proc.time()[3]
global_window_id <- 0L
all_meta <- list()
all_summary <- list()

for (ci in seq_along(chroms)) {
  chr <- chroms[ci]
  mf <- meta_files[ci]

  # Read summary to check if chr had windows
  sumf <- file.path(tmpdir, paste0(chr, ".summary.tsv"))
  if (file.exists(sumf)) {
    sumdt <- fread(sumf)
    if (sumdt$n_windows[1] == 0L) {
      message("[STEP09b-S2] ", chr, ": 0 windows, skipping")
      next
    }
  }

  chr_meta <- fread(mf)
  n <- nrow(chr_meta)
  if (n == 0) next

  # Assign global window_id
  new_ids <- seq(global_window_id + 1L, global_window_id + n)
  chr_meta[, window_id := new_ids]
  global_window_id <- global_window_id + n

  all_meta[[chr]] <- chr_meta

  all_summary[[chr]] <- data.table(
    chrom = chr,
    n_snps = if (file.exists(sumf)) fread(sumf)$n_snps[1] else NA_integer_,
    n_windows = n,
    step_snps = chr_meta$step_snps[1],
    window_snps = chr_meta$window_snps[1],
    first_window_id = new_ids[1],
    last_window_id = new_ids[n]
  )

  message("[STEP09b-S2] ", chr, ": ", n, " windows → IDs ",
          new_ids[1], "–", new_ids[n])
}

# =============================================================================
# PASS 2: PATCH PER-CHR RDS + WRITE FINAL OUTPUTS
# =============================================================================

message("\n[STEP09b-S2] Patching per-chr RDS files...")

for (chr in names(all_meta)) {
  rds_tmp <- file.path(tmpdir, paste0(chr, ".window_pca_tmp.rds"))
  if (!file.exists(rds_tmp)) {
    message("[WARN] Missing tmp RDS for ", chr, " — skipping patch")
    next
  }

  obj <- readRDS(rds_tmp)
  chr_meta <- all_meta[[chr]]

  # Patch window_id into window_meta
  obj$window_meta <- chr_meta

  # Patch window_id into pca table
  pca <- as.data.table(obj$pca)
  pca[, window_id := chr_meta$window_id]
  pca <- pca[, c("window_id", setdiff(names(pca), "window_id")), with = FALSE]
  obj$pca <- as.data.frame(pca)

  # Write final STEP09-compatible RDS
  final_rds <- file.path(outdir, paste0("STEP09_", chr, ".window_pca.rds"))
  saveRDS(obj, final_rds)

  # Write final PCA table
  final_tsv <- file.path(outdir, paste0("STEP09_", chr, ".window_pca.tsv.gz"))
  fwrite(pca, final_tsv, sep = "\t")

  message("[STEP09b-S2] ", chr, " → ", basename(final_rds))
}

# =============================================================================
# WRITE MASTER REGISTRY
# =============================================================================

master <- if (length(all_meta) > 0) rbindlist(all_meta) else {
  data.table(window_id = integer(), chrom = character(), start_bp = integer())
}

summary_dt <- if (length(all_summary) > 0) rbindlist(all_summary) else {
  data.table(chrom = character())
}

# v10: registry artifacts land in registry_dir (02_dense_registry slot),
# per-chrom .window_pca.rds finalized files stay in outdir (01_local_pca).
f1 <- file.path(registry_dir, "windows_master.tsv.gz")
f2 <- file.path(registry_dir, "windows_master_summary.tsv")

fwrite(master, f1, sep = "\t")
fwrite(summary_dt, f2, sep = "\t")

elapsed <- round(proc.time()[3] - t0, 1)

message("\n[DONE] STEP09b Stage 2 — master registry complete (", elapsed, "s)")
message("  Master registry: ", f1, " (", nrow(master), " windows)")
message("  Summary: ", f2)
message("  Chromosomes: ", length(all_meta))
message("  Global window_id range: 1–", global_window_id)
