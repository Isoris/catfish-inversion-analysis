#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_D_merge_registry.R  (v1, 2026-05-12)
# =============================================================================
# Genome-wide finalization for the θπ path. Runs ONCE after all per-chrom
# STEP_TR_C runs have finished.
#
# This is NOT a structural equivalent of path-1's 01c + 02b. Path-1 merges
# because its MDS step (02a chunked_2x) pulls *background windows from
# other chromosomes* and therefore needs genome-wide window IDs to do its
# job. The θπ TR_C MDS is purely per-chromosome — no cross-chrom data
# dependency. So this script does ONLY the two things that genuinely
# benefit from a cross-chrom pass:
#
#   1. Assign a sequential genome-wide global_window_id across all chroms
#      and patch it into each <chr>.precomp.rds (and write a master TSV
#      registry that downstream code can join on by chrom + window_idx).
#   2. Cluster per-window MDS_outlier flags into candidate inversion
#      intervals with a gap-merge tolerance, producing a genome-wide
#      candidate_regions.tsv.gz that the downstream L1/L2 detectors and
#      atlas consume.
#
# It does NOT:
#   - Re-run anything heavy. The MDS / sim_mat / sketch outputs from TR_C
#     are untouched.
#   - Build a chunked-background MDS. That algorithm belongs in path-1
#     where it solves a real problem; in θπ it's noise.
#   - Touch local-mode (dense) precomp files (their MDS columns are NA
#     by design, so they have no MDS_outlier flags to cluster). Sketch
#     and sparse_edges already live in $SKETCH_DIR and don't need merge.
#
# Inputs:
#   --precomp_dir <dir>     directory holding <chr>.precomp.rds from TR_C
#                           (typically $PATH2_PRECOMP = 04_precomp/)
#   --outdir      <dir>     where to write merged registry artifacts.
#                           Default in v10: $PATH2_REGISTRY (02_dense_registry).
#                           v9 default was same as --precomp_dir.
#   --gap_bp      <int>     gap-merge tolerance for candidate clustering
#                           in bp (default 500_000)
#   --min_windows <int>     min windows per candidate region (default 3)
#   --z_thresh    <num>     |max_abs_z| threshold for "outlier" flag
#                           (default 3.0; matches TR_C's downstream
#                           expectation)
#   --patch_rds   <bool>    if "true" (default), rewrites each
#                           <chr>.precomp.rds with global_window_id
#                           added to dt; if "false", leaves RDS untouched
#                           and only writes the master TSV registry.
#
# Outputs (in --outdir):
#   windows_master.tsv.gz             genome-wide registry, columns:
#                                       global_window_id, chrom,
#                                       window_idx, start_bp, end_bp,
#                                       mid_bp, max_abs_z, is_outlier,
#                                       q_min, is_outlier_fdr
#   windows_master_summary.tsv        per-chrom row counts + ID ranges
#   candidate_regions.tsv.gz          candidate inversion intervals:
#                                       candidate_id, chrom, start_bp,
#                                       end_bp, center_bp, n_windows,
#                                       first_global_window_id,
#                                       last_global_window_id,
#                                       max_z_in_region, median_z_in_region
#   candidate_window_membership.tsv.gz region <-> window mapping
#   precomp_summary.tsv               concatenated TR_C per-chrom QC
#                                       (already exists from TR_C, but
#                                       gets a final consolidated copy)
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Args ────────────────────────────────────────────────────────────────────
PRECOMP_DIR <- NA_character_
OUTDIR      <- NA_character_
GAP_BP      <- 500000L
MIN_WINDOWS <- 3L
Z_THRESH    <- 3.0
PATCH_RDS   <- TRUE

args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--precomp_dir" && i < length(args)) { PRECOMP_DIR <- args[i + 1]; i <- i + 2L }
  else if (a == "--outdir"      && i < length(args)) { OUTDIR      <- args[i + 1]; i <- i + 2L }
  else if (a == "--gap_bp"      && i < length(args)) { GAP_BP      <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--min_windows" && i < length(args)) { MIN_WINDOWS <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--z_thresh"    && i < length(args)) { Z_THRESH    <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--fdr_q"       && i < length(args)) { i <- i + 2L }  # legacy no-op
  else if (a == "--patch_rds"   && i < length(args)) { PATCH_RDS   <- tolower(args[i + 1]) %in% c("true", "1", "yes"); i <- i + 2L }
  else { i <- i + 1L }
}
stopifnot(!is.na(PRECOMP_DIR))
# v10: default --outdir is now $PATH2_REGISTRY (02_dense_registry/), with
# fallback to precomp_dir for v9 back-compat.
if (is.na(OUTDIR)) {
  env_reg <- Sys.getenv("PATH2_REGISTRY", unset = "")
  OUTDIR <- if (nzchar(env_reg)) env_reg else PRECOMP_DIR
}
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

precomp_files <- sort(list.files(PRECOMP_DIR, pattern = "\\.precomp\\.rds$",
                                 full.names = TRUE))
if (length(precomp_files) == 0L) stop("[TR_D] no .precomp.rds in ", PRECOMP_DIR)

message(sprintf("[TR_D] precomp_dir=%s outdir=%s", PRECOMP_DIR, OUTDIR))
message(sprintf("[TR_D] z_thresh=%.2f patch_rds=%s",
                Z_THRESH, as.character(PATCH_RDS)))
message("[TR_D] found ", length(precomp_files), " precomp files")

# =============================================================================
# Pass 1: read each precomp's dt, assemble master registry with global IDs
# =============================================================================
t0 <- proc.time()[3]
global_id <- 0L
master_rows <- list()
summary_rows <- list()
chrom_order <- character()
local_mode_chroms <- character()

for (f in precomp_files) {
  pc <- readRDS(f)
  chr <- pc$chrom
  if (!is.null(pc$mode) && pc$mode == "local") {
    # Local mode has no MDS_outlier flags; skip from candidate clustering
    # but still include in registry for cross-chrom coordinate joins.
    local_mode_chroms <- c(local_mode_chroms, chr)
  }
  dt <- pc$dt
  n <- nrow(dt)
  if (n == 0L) { rm(pc); next }
  ids <- seq.int(global_id + 1L, global_id + n)
  global_id <- global_id + n

  reg <- data.table(
    global_window_id = ids,
    chrom            = chr,
    window_idx       = dt$window_idx,
    start_bp         = dt$start_bp,
    end_bp           = dt$end_bp,
    mid_bp           = dt$mid_bp,
    max_abs_z        = dt$max_abs_z,
    is_outlier       = is.finite(dt$max_abs_z) & dt$max_abs_z >= Z_THRESH,
    mode             = pc$mode %||% NA_character_
  )
  master_rows[[chr]] <- reg
  summary_rows[[chr]] <- data.table(
    chrom            = chr,
    n_windows        = n,
    first_global_id  = ids[1L],
    last_global_id   = ids[n],
    n_outlier        = sum(reg$is_outlier, na.rm = TRUE),
    median_max_z     = round(median(reg$max_abs_z, na.rm = TRUE), 3),
    q95_max_z        = round(quantile(reg$max_abs_z, 0.95, na.rm = TRUE), 3),
    mode             = pc$mode %||% NA_character_
  )
  chrom_order <- c(chrom_order, chr)
  message(sprintf("[TR_D] %s: %d windows -> global_id %d..%d (%d z-outliers; gap-bp merge disabled)",
                  chr, n, ids[1L], ids[n],
                  summary_rows[[chr]]$n_outlier))
  rm(pc, dt, reg); invisible(gc(verbose = FALSE))
}

master <- rbindlist(master_rows, fill = TRUE)
summary_dt <- rbindlist(summary_rows, fill = TRUE)

if (length(local_mode_chroms) > 0L) {
  message("[TR_D] note: ", length(local_mode_chroms),
          " local-mode chrom(s) included in registry but excluded from",
          " candidate clustering (no MDS axes): ",
          paste(local_mode_chroms, collapse = ","))
}

fwrite(master, file.path(OUTDIR, "windows_master.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(summary_dt, file.path(OUTDIR, "windows_master_summary.tsv"), sep = "\t")
message("[TR_D] wrote windows_master.tsv.gz (", nrow(master), " rows)")

# =============================================================================
# Pass 2 (optional): patch each .precomp.rds with global_window_id
# =============================================================================
if (PATCH_RDS) {
  message("[TR_D] patching per-chrom precomp.rds with global_window_id...")
  for (f in precomp_files) {
    pc <- readRDS(f)
    chr <- pc$chrom
    reg <- master_rows[[chr]]
    if (is.null(reg) || nrow(reg) != nrow(pc$dt)) {
      warning("[TR_D] ", chr, ": row count mismatch — skipping RDS patch")
      next
    }
    # Insert global_window_id as the first column of dt; cheap, no other change
    pc$dt[, global_window_id := reg$global_window_id]
    setcolorder(pc$dt, c("global_window_id",
                          setdiff(names(pc$dt), "global_window_id")))
    saveRDS(pc, f)
    rm(pc); invisible(gc(verbose = FALSE))
  }
  message("[TR_D] RDS patch complete")
} else {
  message("[TR_D] skipping RDS patch (--patch_rds false); registry-only mode")
}

# =============================================================================
# Pass 3: candidate-region clustering DISABLED 2026-05-13
# =============================================================================
# Legacy lostruct-style gap-bp merge of is_outlier flags. Not consumed by
# L1/L2 stripe detection (which reads sim_mat directly). In dense-outlier
# regimes the gap-bp walk over-merges across noise. Empty placeholder files
# are still written so any downstream tooling that auto-resolves the paths
# doesn't break.
# =============================================================================
cand_dt       <- data.table(candidate_id = integer(), chrom = character())
membership_dt <- data.table(candidate_id = integer(), global_window_id = integer())

fwrite(cand_dt, file.path(OUTDIR, "candidate_regions.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(membership_dt, file.path(OUTDIR, "candidate_window_membership.tsv.gz"),
       sep = "\t", compress = "gzip")

elapsed <- round(proc.time()[3] - t0, 1)
message("[TR_D] candidate regions: 0 (gap-bp merge disabled)")
message(sprintf("[TR_D] DONE in %.1fs", elapsed))
message("  registry:    ", file.path(OUTDIR, "windows_master.tsv.gz"))
