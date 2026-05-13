#!/usr/bin/env Rscript
# =============================================================================
# STEP_GH_D_merge_registry.R  (v1, 2026-05-12)
# =============================================================================
# Genome-wide finalization for the GHSL path. Runs ONCE after all per-chrom
# STEP_GH_C runs have finished.
#
# Direct analog of STEP_TR_D_merge_registry.R — same algorithm, same
# outputs, just operates on <chr>.ghsl_precomp.rds and writes to a
# GHSL-namespaced location so it coexists with the θπ merge artifacts.
#
# Same architectural caveat as TR_D: this is NOT a structural equivalent
# of path-1's 02b. Path-1 merges because its MDS step (02a chunked_2x)
# pulls *background windows from other chromosomes* and needs genome-wide
# window IDs to do that. GHSL's GH_C MDS is purely per-chromosome — no
# cross-chrom data dependency. So this script does ONLY:
#
#   1. Assign sequential global_window_id across all chroms and (optionally)
#      patch into each <chr>.ghsl_precomp.rds; always write the master
#      TSV registry that downstream code can join on.
#   2. Cluster MDS_outlier flags into candidate inversion regions.
#
# Local-mode (dense) GHSL precomps have NA MDS columns by design; they're
# included in the registry for cross-chrom coordinate joins but excluded
# from candidate clustering. Their sketch + sparse_edges live in
# $SKETCH_DIR and don't need merge.
#
# Inputs:
#   --precomp_dir <dir>     directory of <chr>.ghsl_precomp.rds from GH_C
#                           (typically $OUTROOT/ghsl_precomp_coarse)
#   --outdir      <dir>     where to write merged artifacts (default:
#                           same as precomp_dir)
#   --gap_bp      <int>     gap-merge tolerance (default 500_000)
#   --min_windows <int>     min windows per region (default 3)
#   --z_thresh    <num>     |max_abs_z| outlier threshold (default 3.0)
#   --patch_rds   <bool>    if "true" (default), patches each
#                           <chr>.ghsl_precomp.rds with global_window_id
#
# Outputs (in --outdir):
#   ghsl_windows_master.tsv.gz        adds columns q_min, is_outlier_fdr
#   ghsl_windows_master_summary.tsv   adds column n_outlier_fdr
#   ghsl_candidate_regions.tsv.gz
#   ghsl_candidate_window_membership.tsv.gz
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

`%||%` <- function(a, b) if (is.null(a)) b else a

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
if (is.na(OUTDIR)) OUTDIR <- PRECOMP_DIR
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

precomp_files <- sort(list.files(PRECOMP_DIR, pattern = "\\.ghsl_precomp\\.rds$",
                                 full.names = TRUE))
if (length(precomp_files) == 0L) stop("[GH_D] no .ghsl_precomp.rds in ", PRECOMP_DIR)

message(sprintf("[GH_D] precomp_dir=%s outdir=%s", PRECOMP_DIR, OUTDIR))
message(sprintf("[GH_D] z_thresh=%.2f patch_rds=%s",
                Z_THRESH, as.character(PATCH_RDS)))
message("[GH_D] found ", length(precomp_files), " precomp files")

# ── Pass 1: registry ────────────────────────────────────────────────────────
t0 <- proc.time()[3]
global_id <- 0L
master_rows <- list(); summary_rows <- list()
chrom_order <- character(); local_mode_chroms <- character()

for (f in precomp_files) {
  pc <- readRDS(f); chr <- pc$chrom
  if (!is.null(pc$mode) && pc$mode == "local") local_mode_chroms <- c(local_mode_chroms, chr)
  dt <- pc$dt; n <- nrow(dt)
  if (n == 0L) { rm(pc); next }
  ids <- seq.int(global_id + 1L, global_id + n); global_id <- global_id + n

  reg <- data.table(
    global_window_id = ids,
    chrom            = chr,
    window_idx       = dt$window_idx,
    start_bp         = dt$start_bp,
    end_bp           = dt$end_bp,
    mid_bp           = dt$mid_bp,
    max_abs_z        = dt$max_abs_z,
    is_outlier       = is.finite(dt$max_abs_z) & dt$max_abs_z >= Z_THRESH,
    mode             = pc$mode %||% NA_character_,
    scale            = pc$scale %||% NA_character_
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
    mode             = pc$mode %||% NA_character_,
    scale            = pc$scale %||% NA_character_
  )
  chrom_order <- c(chrom_order, chr)
  message(sprintf("[GH_D] %s: %d windows -> global_id %d..%d (%d z-outliers, scale=%s; gap-bp merge disabled)",
                  chr, n, ids[1L], ids[n],
                  summary_rows[[chr]]$n_outlier,
                  pc$scale %||% "?"))
  rm(pc, dt, reg); invisible(gc(verbose = FALSE))
}

master     <- rbindlist(master_rows, fill = TRUE)
summary_dt <- rbindlist(summary_rows, fill = TRUE)
if (length(local_mode_chroms) > 0L) {
  message("[GH_D] note: ", length(local_mode_chroms),
          " local-mode chrom(s) excluded from candidate clustering: ",
          paste(local_mode_chroms, collapse = ","))
}
fwrite(master, file.path(OUTDIR, "ghsl_windows_master.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(summary_dt, file.path(OUTDIR, "ghsl_windows_master_summary.tsv"), sep = "\t")
message("[GH_D] wrote ghsl_windows_master.tsv.gz (", nrow(master), " rows)")

# ── Pass 2: patch RDS ───────────────────────────────────────────────────────
if (PATCH_RDS) {
  message("[GH_D] patching per-chrom ghsl_precomp.rds with global_window_id...")
  for (f in precomp_files) {
    pc <- readRDS(f); chr <- pc$chrom
    reg <- master_rows[[chr]]
    if (is.null(reg) || nrow(reg) != nrow(pc$dt)) {
      warning("[GH_D] ", chr, ": row count mismatch — skipping RDS patch"); next
    }
    pc$dt[, global_window_id := reg$global_window_id]
    setcolorder(pc$dt, c("global_window_id",
                          setdiff(names(pc$dt), "global_window_id")))
    saveRDS(pc, f)
    rm(pc); invisible(gc(verbose = FALSE))
  }
  message("[GH_D] RDS patch complete")
}

# ── Pass 3: candidate regions DISABLED 2026-05-13 ───────────────────────────
# Legacy lostruct-style gap-bp merge of is_outlier flags. Not consumed by
# L1/L2 stripe detection (which reads sim_mat directly). Empty placeholder
# files still written so downstream path auto-resolution doesn't break.
cand_dt       <- data.table(candidate_id = integer(), chrom = character())
membership_dt <- data.table(candidate_id = integer(), global_window_id = integer())

fwrite(cand_dt, file.path(OUTDIR, "ghsl_candidate_regions.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(membership_dt, file.path(OUTDIR, "ghsl_candidate_window_membership.tsv.gz"),
       sep = "\t", compress = "gzip")

elapsed <- round(proc.time()[3] - t0, 1)
message("[GH_D] candidate regions: 0 (gap-bp merge disabled)")
message(sprintf("[GH_D] DONE in %.1fs", elapsed))
message("  registry:    ", file.path(OUTDIR, "ghsl_windows_master.tsv.gz"))
