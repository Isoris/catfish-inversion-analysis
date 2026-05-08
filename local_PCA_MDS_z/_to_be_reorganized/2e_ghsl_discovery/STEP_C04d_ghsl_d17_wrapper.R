#!/usr/bin/env Rscript
# =============================================================================
# STEP_C04d_ghsl_d17_wrapper.R
# =============================================================================
# Phase 2 / 2e_ghsl_discovery — adapter that runs the proven D17 boundary
# detector on the GHSL local-PCA outputs, producing L1/L2 boundaries and
# envelopes in the same TSV shape as the dosage pipeline.
#
# Why this exists
# ---------------
# STEP_C04c (turn 2) ships a |Z|-threshold envelope detector as SECONDARY.
# That's a 1D scan — sees peaks but doesn't see edges. The dosage pipeline
# uses a much sharper detector (D17): a diagonal cross-block scan that asks
# "are samples on the LEFT dissimilar from samples on the RIGHT at this
# position?", which is exactly what an inversion edge looks like. D17 has
# been calibrated on dosage and ships with a validated growing-W validator,
# adaptive thresholding, Ward-D2 segment-internal-homogeneity stratification,
# and STABLE_BLUE / EDGE / FAKE / DEMOTED status labels.
#
# D17's core statistic is signal-agnostic: per-diagonal Z-normalize sim_mat,
# then median over an upper-triangle WxW cross-block. That works on ANY
# sim_mat with self-similarity ~ 1 and pairwise similarity decreasing with
# pattern dissimilarity — which describes GHSL's |cor(pc1[i], pc1[j])|
# sim_mat exactly. The defaults were tuned for dosage but the adaptive
# thresholding (default mode) self-calibrates per chromosome from observed
# grow_max_z quantiles, so it auto-adapts to whatever the typical signal
# looks like for the input.
#
# This wrapper:
#  1. Reads <chr>.ghsl_v6_localpca.rds (STEP_C04c output)
#  2. Writes a temp precomp.rds matching D17's expected shape: pc$dt with
#     start_bp / end_bp columns (D17 also accepts start/end, window_start/end,
#     bp_start/bp_end — we use start_bp/end_bp by default).
#  3. Writes a temp sim_mat.rds (D17 accepts bare matrix or list with $sim_mat;
#     we use bare matrix for simplicity).
#  4. Invokes STEP_D17_multipass_L1_only_v7.R with the temp paths.
#  5. Invokes STEP_D17_multipass_L2_v8.R consuming the L1 envelopes.
#  6. Renames outputs from <chr_label>_d17L1_*.tsv to <chr>_ghsl_d17L1_*.tsv
#     so they don't collide with the dosage pipeline's outputs.
#
# Output TSVs (consumed natively by export_precomp_to_json_v3.R via the
# --l1_envelopes / --l2_envelopes / --l2_boundaries / --l1_boundaries
# flags):
#   <chr>_ghsl_d17L1_envelopes.tsv
#   <chr>_ghsl_d17L1_boundaries.tsv
#   <chr>_ghsl_d17L1_boundary_score_curve.tsv
#   <chr>_ghsl_d17L2_envelopes.tsv
#   <chr>_ghsl_d17L2_boundaries.tsv
#
# Asymmetry vs θπ: same wrapper exists for θπ (STEP_TR_C_theta_d17_wrapper.R).
# Both produce parallel L1/L2 TSVs. The page-3 atlas (GHSL) and page-12 atlas
# (θπ) overlay them on their respective |Z| traces.
#
# Note on detector hierarchy:
#  - For GHSL: D17 envelopes (this script) are PRIMARY architecturally
#    (geometry-based edge detection), STEP_C04b PASS-runs are PRIMARY
#    biologically (calibrated on real signal). They answer different
#    questions — keep both, don't reconcile, the atlas overlays both.
#    STEP_C04c's |Z|-threshold secondaries become tertiary when D17 is
#    available; in practice the atlas can render any subset based on
#    what's present.
#  - For θπ: D17 envelopes (the parallel wrapper) PLUS STEP_TR_B's
#    |Z|-threshold envelopes both ship as primary candidate sets. No
#    upstream production detector exists for θπ.
#
# Performance notes
# -----------------
# D17's main hot paths are O(N²), NOT O(N³) — earlier turn-2/3 thinking
# called this cubic in error, see correction in turn 6.
#   * Step 1 (diagonal mean/sd precompute): one pass over each diagonal,
#     Σ(N-d) ≈ N²/2 ops total.
#   * Grow validator: n_peaks · Σ W² across the grow series. n_peaks
#     scales with biology (~5-50 per chromosome), not N. Σ W² is the
#     cost driver — when W series is percent-of-N, this also scales O(N²).
#
# Empirical estimates (R wall time, single chromosome):
#   N=4,300  (GHSL 5kb on LG28):  ~5-10s   total
#   N=8,000:                        ~15-30s  total
#   N=16,500 (θπ fine grid):        ~60-120s total (estimated by N² scaling)
#
# For fine-grid θπ, if wall time is too long, pass a reduced grow series
# via --d17_l1_args. Default percent series is
#   0.001,0.005,0.01,0.02,0.04,0.07     (6 W values, largest 0.07·N)
# Faster alternatives:
#   --d17_l1_args "--boundary_grow_W_pct 0.005,0.02,0.05"   (~2.4× speedup)
#   --d17_l1_args "--boundary_grow_W_pct 0.005,0.02"        (much faster)
#
# Trade-off: fewer grow W values = less robust REAL/FAKE classification.
# The growing-W validator works because real boundaries stay blue at all
# W while fake boundaries drift positive at large W. Cutting the W
# series narrows that test. For first-pass exploration on unknown θπ
# data, fewer W values is fine.
#
# Usage
# -----
#   Rscript STEP_C04d_ghsl_d17_wrapper.R \
#     --localpca       <localpca_dir>/<CHR>.ghsl_v6_localpca.rds \
#     --d17_l1_script  /path/to/STEP_D17_multipass_L1_only_v7.R \
#     --d17_l2_script  /path/to/STEP_D17_multipass_L2_v8.R \
#     --outdir         <out_dir> \
#     [--chr           <chr_label>]                            \
#     [--keep_temp]                                            \
#     [--d17_l1_args   "<extra args to forward to L1 script>"] \
#     [--d17_l2_args   "<extra args to forward to L2 script>"] \
#     [--rscript       /path/to/Rscript]
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- CLI --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
has_flag <- function(flag) any(args == flag)

LOCALPCA      <- get_arg("--localpca")
D17_L1_SCRIPT <- get_arg("--d17_l1_script")
D17_L2_SCRIPT <- get_arg("--d17_l2_script")
OUTDIR        <- get_arg("--outdir", ".")
CHR_LABEL     <- get_arg("--chr", NULL)
KEEP_TEMP     <- has_flag("--keep_temp")
EXTRA_L1      <- get_arg("--d17_l1_args", "")
EXTRA_L2      <- get_arg("--d17_l2_args", "")
RSCRIPT       <- get_arg("--rscript", "Rscript")

if (is.null(LOCALPCA) || !file.exists(LOCALPCA)) {
  stop("[C04d-wrap] FATAL: --localpca <file>.ghsl_v6_localpca.rds is required")
}
if (is.null(D17_L1_SCRIPT) || !file.exists(D17_L1_SCRIPT)) {
  stop("[C04d-wrap] FATAL: --d17_l1_script must point at STEP_D17_multipass_L1_only_v7.R")
}
if (is.null(D17_L2_SCRIPT) || !file.exists(D17_L2_SCRIPT)) {
  stop("[C04d-wrap] FATAL: --d17_l2_script must point at STEP_D17_multipass_L2_v8.R")
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("================================================================")
message("[C04d-wrap] GHSL local-PCA → D17 boundary detector wrapper")
message("================================================================")
message("[C04d-wrap] localpca RDS : ", LOCALPCA)
message("[C04d-wrap] D17 L1 script: ", D17_L1_SCRIPT)
message("[C04d-wrap] D17 L2 script: ", D17_L2_SCRIPT)
message("[C04d-wrap] outdir       : ", OUTDIR)

# ---- Load STEP_C04c output --------------------------------------------------
message("[C04d-wrap] Reading STEP_C04c local-PCA RDS")
lp <- readRDS(LOCALPCA)
chrom <- if (!is.null(CHR_LABEL)) CHR_LABEL else lp$chrom
n_win <- lp$n_windows %||% nrow(lp$window_info)

if (is.null(chrom) || is.na(chrom)) {
  stop("[C04d-wrap] FATAL: cannot determine chrom (not in RDS, not in --chr)")
}
message("[C04d-wrap] chrom: ", chrom, "   n_windows: ", n_win)

# Sanity checks
if (is.null(lp$window_info) || nrow(lp$window_info) != n_win) {
  stop(sprintf("[C04d-wrap] FATAL: window_info shape mismatch (n=%d, want %d)",
               nrow(lp$window_info %||% data.frame()), n_win))
}
if (is.null(lp$sim_mat) || !is.matrix(lp$sim_mat)) {
  stop("[C04d-wrap] FATAL: lp$sim_mat is not a matrix (was STEP_C04c output corrupted?)")
}
if (nrow(lp$sim_mat) != n_win || ncol(lp$sim_mat) != n_win) {
  stop(sprintf("[C04d-wrap] FATAL: sim_mat is %dx%d but n_windows=%d",
               nrow(lp$sim_mat), ncol(lp$sim_mat), n_win))
}

# ---- Build D17-compatible precomp.rds ---------------------------------------
# D17 accepts pc$dt with start_bp/end_bp columns (line 308 of D17 L1 v7).
# We add the PC1/PC2 columns too in case any downstream consumer of the
# precomp expects them — D17 itself doesn't read them, but other tools
# (e.g. export_precomp_to_json_v3.R) do.
message("[C04d-wrap] Building D17-compatible precomp")
win_info <- lp$window_info
pc_dt <- data.table(
  start_bp = as.integer(win_info$start_bp),
  end_bp   = as.integer(win_info$end_bp),
  chrom    = chrom,
  lam_1    = as.numeric(lp$lambda_1 %||% rep(NA_real_, n_win)),
  lam_2    = as.numeric(lp$lambda_2 %||% rep(NA_real_, n_win)),
  z        = as.numeric(lp$z_profile %||% rep(NA_real_, n_win))
)

# Attach per-sample PC1 / PC2 columns (sign-aligned for downstream consumers)
sample_names <- lp$sample_names
if (!is.null(lp$pc1_loadings_aligned_mat)) {
  pc1_aligned <- lp$pc1_loadings_aligned_mat
  for (s_idx in seq_len(ncol(pc1_aligned))) {
    s_name <- sample_names[s_idx]
    pc_dt[[paste0("PC_1_", s_name)]] <- as.numeric(pc1_aligned[, s_idx])
  }
}
if (!is.null(lp$pc2_loadings_aligned_mat)) {
  pc2_aligned <- lp$pc2_loadings_aligned_mat
  for (s_idx in seq_len(ncol(pc2_aligned))) {
    s_name <- sample_names[s_idx]
    pc_dt[[paste0("PC_2_", s_name)]] <- as.numeric(pc2_aligned[, s_idx])
  }
}

precomp_obj <- list(
  dt          = pc_dt,
  chrom       = chrom,
  n_windows   = n_win,
  sample_names = sample_names,
  source      = "STEP_C04c local-PCA → C04d wrapper",
  generator   = "STEP_C04d_ghsl_d17_wrapper.R"
)

# Temp files in outdir (under a sub-folder so they're easy to delete)
temp_dir <- file.path(OUTDIR, paste0(".d17_temp_", chrom))
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

precomp_temp <- file.path(temp_dir, paste0(chrom, ".precomp.slim.rds"))
sim_mat_temp <- file.path(temp_dir, paste0(chrom, ".sim_mat.rds"))

message("[C04d-wrap] Writing temp precomp: ", precomp_temp)
saveRDS(precomp_obj, precomp_temp)

message("[C04d-wrap] Writing temp sim_mat: ", sim_mat_temp)
saveRDS(lp$sim_mat, sim_mat_temp)

# ---- Run D17 L1 -------------------------------------------------------------
# IMPORTANT: D17 by default does NOT run the boundary scan. The user must
# pass --boundary_scan to enable it. We always add this flag.
message("\n[C04d-wrap] === Invoking D17 L1 ===")
l1_cmd_parts <- c(
  shQuote(RSCRIPT),
  shQuote(D17_L1_SCRIPT),
  "--precomp",      shQuote(precomp_temp),
  "--sim_mat",      shQuote(sim_mat_temp),
  "--chr",          shQuote(chrom),
  "--outdir",       shQuote(OUTDIR),
  "--boundary_scan"
)
if (nchar(EXTRA_L1) > 0) {
  l1_cmd_parts <- c(l1_cmd_parts, EXTRA_L1)
}
l1_cmd <- paste(l1_cmd_parts, collapse = " ")
message("[C04d-wrap] L1 cmd: ", l1_cmd)
l1_status <- system(l1_cmd)
if (l1_status != 0) {
  stop("[C04d-wrap] FATAL: D17 L1 failed with exit code ", l1_status)
}

# Verify L1 outputs exist
l1_envelope_path  <- file.path(OUTDIR, paste0(chrom, "_d17L1_envelopes.tsv"))
l1_boundary_path  <- file.path(OUTDIR, paste0(chrom, "_d17L1_boundaries.tsv"))
l1_curve_path     <- file.path(OUTDIR, paste0(chrom, "_d17L1_boundary_score_curve.tsv"))
if (!file.exists(l1_envelope_path)) {
  stop("[C04d-wrap] FATAL: D17 L1 didn't produce expected envelope output: ",
       l1_envelope_path)
}

n_l1_env <- nrow(fread(l1_envelope_path))
n_l1_bnd <- if (file.exists(l1_boundary_path)) nrow(fread(l1_boundary_path)) else 0
message("[C04d-wrap] D17 L1 produced: ", n_l1_env, " L1 envelopes, ",
        n_l1_bnd, " L1 boundary peaks")

# ---- Run D17 L2 -------------------------------------------------------------
message("\n[C04d-wrap] === Invoking D17 L2 ===")
l2_cmd_parts <- c(
  shQuote(RSCRIPT),
  shQuote(D17_L2_SCRIPT),
  "--precomp",      shQuote(precomp_temp),
  "--sim_mat",      shQuote(sim_mat_temp),
  "--l1_catalogue", shQuote(l1_envelope_path),
  "--chr",          shQuote(chrom),
  "--outdir",       shQuote(OUTDIR)
)
if (nchar(EXTRA_L2) > 0) {
  l2_cmd_parts <- c(l2_cmd_parts, EXTRA_L2)
}
l2_cmd <- paste(l2_cmd_parts, collapse = " ")
message("[C04d-wrap] L2 cmd: ", l2_cmd)
l2_status <- system(l2_cmd)
if (l2_status != 0) {
  stop("[C04d-wrap] FATAL: D17 L2 failed with exit code ", l2_status)
}

l2_envelope_path <- file.path(OUTDIR, paste0(chrom, "_d17L2_envelopes.tsv"))
l2_boundary_path <- file.path(OUTDIR, paste0(chrom, "_d17L2_boundaries.tsv"))
n_l2_env <- if (file.exists(l2_envelope_path)) nrow(fread(l2_envelope_path)) else 0
n_l2_bnd <- if (file.exists(l2_boundary_path)) nrow(fread(l2_boundary_path)) else 0
message("[C04d-wrap] D17 L2 produced: ", n_l2_env, " L2 envelopes, ",
        n_l2_bnd, " L2 boundary peaks")

# ---- Rename outputs to ghsl-namespaced form ---------------------------------
# So they don't collide with dosage pipeline outputs in a shared outdir.
rename_outputs <- function(orig_path, ghsl_suffix) {
  if (!file.exists(orig_path)) return(invisible())
  new_path <- sub("_d17L([12])", paste0("_ghsl_d17L\\1"), orig_path)
  if (orig_path != new_path) {
    file.rename(orig_path, new_path)
    message("[C04d-wrap]   renamed: ", basename(orig_path),
            " → ", basename(new_path))
  }
}

message("\n[C04d-wrap] Renaming outputs to GHSL namespace")
rename_outputs(l1_envelope_path,  "ghsl")
rename_outputs(l1_boundary_path,  "ghsl")
rename_outputs(l1_curve_path,     "ghsl")
rename_outputs(l2_envelope_path,  "ghsl")
rename_outputs(l2_boundary_path,  "ghsl")

# Clean up temp files unless --keep_temp
if (!KEEP_TEMP) {
  message("[C04d-wrap] Cleaning temp dir: ", temp_dir)
  unlink(temp_dir, recursive = TRUE)
} else {
  message("[C04d-wrap] Keeping temp files (--keep_temp): ", temp_dir)
}

message("\n================================================================")
message("[DONE] GHSL D17 boundary detection complete for ", chrom)
message("================================================================")
message("  Outputs in: ", OUTDIR)
message("    <chr>_ghsl_d17L1_envelopes.tsv")
message("    <chr>_ghsl_d17L1_boundaries.tsv")
message("    <chr>_ghsl_d17L1_boundary_score_curve.tsv")
message("    <chr>_ghsl_d17L2_envelopes.tsv")
message("    <chr>_ghsl_d17L2_boundaries.tsv")
message("")
message("  Next: feed into export_precomp_to_json_v3.R via")
message("    --l1_envelopes  <chr>_ghsl_d17L1_envelopes.tsv")
message("    --l2_envelopes  <chr>_ghsl_d17L2_envelopes.tsv")
message("    --l1_boundaries <chr>_ghsl_d17L1_boundaries.tsv")
message("    --l2_boundaries <chr>_ghsl_d17L2_boundaries.tsv")
