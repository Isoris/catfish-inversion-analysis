#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_C_theta_d17_wrapper.R
# =============================================================================
# Phase 2 / 2f_theta_discovery — adapter that runs the proven D17 boundary
# detector on θπ local-PCA outputs, producing L1/L2 boundaries and envelopes
# in the same TSV shape as the dosage and GHSL pipelines.
#
# Why this exists
# ---------------
# STEP_TR_B v4 ships a |Z|-threshold envelope detector (PRIMARY for θπ since
# no upstream production θπ candidate detector exists yet). That's a 1D scan
# — sees peaks but doesn't see edges. For symmetry with GHSL's C04d wrapper,
# we ALSO run D17 on θπ to get geometry-based edge detection (cross-block
# scan: "are samples on the LEFT dissimilar from samples on the RIGHT here?").
#
# D17 needs a window×window sim_mat. STEP_TR_B v4 doesn't currently emit
# one — the v4 walk-back decision was "θπ doesn't need sim_mat for the
# 1D |Z| profile, which is sign-invariant by construction". That's still
# true for the |Z| profile, but D17 needs sim_mat for the cross-block scan.
# This wrapper computes sim_mat on the fly from STEP_TR_B's JSON output's
# pc1_loadings array.
#
# Output sim_mat size: at win10000.step2000 = ~16,500 windows on LG28, the
# sim_mat is ~1.1 GB float64 in R memory. Fine cluster-side (LANTA has
# plenty of RAM); never serialized to JSON for the browser. The D17 outputs
# are small TSVs that the browser handles trivially.
#
# This wrapper:
#  1. Reads STEP_TR_B's <chr>_phase2_theta.json
#  2. Extracts theta_pi_local_pca.pc1_loadings (n_windows × n_samples)
#  3. Computes sim_mat = |cor(pc1[i], pc1[j])| in R memory
#  4. Writes a temp precomp.rds matching D17's expected shape
#  5. Writes a temp sim_mat.rds (bare matrix)
#  6. Invokes STEP_D17_multipass_L1_only_v7.R with --boundary_scan
#  7. Invokes STEP_D17_multipass_L2_v8.R consuming L1 envelopes
#  8. Renames outputs to <chr>_theta_d17L*.tsv (avoid collision with dosage/GHSL)
#
# Output TSVs (consumed natively by export_precomp_to_json_v3.R via the
# --l1_envelopes / --l2_envelopes / --l2_boundaries / --l1_boundaries flags):
#   <chr>_theta_d17L1_envelopes.tsv
#   <chr>_theta_d17L1_boundaries.tsv
#   <chr>_theta_d17L1_boundary_score_curve.tsv
#   <chr>_theta_d17L2_envelopes.tsv
#   <chr>_theta_d17L2_boundaries.tsv
#
# These ship ALONGSIDE STEP_TR_B v4's |Z|-threshold envelopes (which stay
# as-is in the page-12 JSON's theta_pi_envelopes layer). The atlas overlays
# both — peak-based primaries from |Z|, edge-based secondaries from D17 —
# on the page-12 |Z| panel, distinguished by color.
#
# Asymmetry vs GHSL: same script (D17 boundary detector) runs on both
# streams. The conceptual difference (GHSL has STEP_C04b as upstream
# production candidate detector; θπ doesn't have an equivalent) is reflected
# only in the page-3 vs page-12 atlas's hierarchy of which layer is shown
# as primary vs secondary, not in the wrapper logic.
#
# Usage
# -----
#   Rscript STEP_TR_C_theta_d17_wrapper.R \
#     --theta_json     <chr>_phase2_theta.json \
#     --d17_l1_script  /path/to/STEP_D17_multipass_L1_only_v7.R \
#     --d17_l2_script  /path/to/STEP_D17_multipass_L2_v8.R \
#     --outdir         <out_dir> \
#     [--chr           <chr_label>]                            \
#     [--keep_temp]                                            \
#     [--d17_l1_args   "<extra args to forward to L1 script>"] \
#     [--d17_l2_args   "<extra args to forward to L2 script>"] \
#     [--rscript       /path/to/Rscript]
#
# Wall time: depends heavily on fine-grid N. Empirical estimate from
#   numpy benchmarks (R is ~5-10× slower at hot loops):
#     N=8,000   sim_mat ~15s (cor in C),  D17 L1 ~15-30s,   D17 L2 ~5-10s
#     N=16,500  sim_mat ~60s,             D17 L1 ~60-120s,  D17 L2 ~15-30s
#   Total at θπ fine grid: ~3-4 minutes per chromosome.
#
# If too slow on real data, the right fix is NOT optimizing D17 (validated
# code, don't touch) but tuning the grow W series via the wrapper's
# --d17_l1_args flag. See the GHSL wrapper's Performance notes section
# for concrete examples.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
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

THETA_JSON    <- get_arg("--theta_json")
D17_L1_SCRIPT <- get_arg("--d17_l1_script")
D17_L2_SCRIPT <- get_arg("--d17_l2_script")
OUTDIR        <- get_arg("--outdir", ".")
CHR_LABEL     <- get_arg("--chr", NULL)
KEEP_TEMP     <- has_flag("--keep_temp")
EXTRA_L1      <- get_arg("--d17_l1_args", "")
EXTRA_L2      <- get_arg("--d17_l2_args", "")
RSCRIPT       <- get_arg("--rscript", "Rscript")

if (is.null(THETA_JSON) || !file.exists(THETA_JSON)) {
  stop("[TR_C-wrap] FATAL: --theta_json <chr>_phase2_theta.json is required")
}
if (is.null(D17_L1_SCRIPT) || !file.exists(D17_L1_SCRIPT)) {
  stop("[TR_C-wrap] FATAL: --d17_l1_script must point at STEP_D17_multipass_L1_only_v7.R")
}
if (is.null(D17_L2_SCRIPT) || !file.exists(D17_L2_SCRIPT)) {
  stop("[TR_C-wrap] FATAL: --d17_l2_script must point at STEP_D17_multipass_L2_v8.R")
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("================================================================")
message("[TR_C-wrap] θπ → D17 boundary detector wrapper")
message("================================================================")
message("[TR_C-wrap] theta_json   : ", THETA_JSON)
message("[TR_C-wrap] D17 L1 script: ", D17_L1_SCRIPT)
message("[TR_C-wrap] D17 L2 script: ", D17_L2_SCRIPT)
message("[TR_C-wrap] outdir       : ", OUTDIR)

# ---- Read STEP_TR_B's JSON output -------------------------------------------
message("[TR_C-wrap] Reading STEP_TR_B JSON")
js <- fromJSON(THETA_JSON, simplifyVector = FALSE)
chrom <- if (!is.null(CHR_LABEL)) CHR_LABEL else js$chrom
n_samp <- as.integer(js$n_samples)
n_win  <- as.integer(js$n_windows)
sample_order <- unlist(js$theta_pi_local_pca$sample_order)

if (length(sample_order) != n_samp) {
  stop("[TR_C-wrap] FATAL: sample_order length ", length(sample_order),
       " != n_samples ", n_samp)
}
message("[TR_C-wrap] chrom=", chrom, "  n_samples=", n_samp,
        "  n_windows=", n_win)

# Window grid — the JSON's theta_pi_per_window layer carries pos_bp
ppw <- js$theta_pi_per_window
window_pos_bp <- as.integer(unlist(ppw$pos_bp %||% js$tracks$theta_pi_z$pos_bp))
if (length(window_pos_bp) != n_win) {
  stop("[TR_C-wrap] FATAL: window pos_bp length ", length(window_pos_bp),
       " != n_windows ", n_win)
}

# Reconstruct start_bp / end_bp from pos_bp (D17 wants both). We assume
# the windows are uniform; window size is the diff between consecutive
# pos_bp. Last window's end is pos + half-step.
window_step <- median(diff(window_pos_bp), na.rm = TRUE)
half_step <- floor(window_step / 2)
window_start_bp <- pmax(0L, window_pos_bp - half_step)
window_end_bp   <- window_pos_bp + (window_step - half_step) - 1L
message("[TR_C-wrap] inferred window step: ", window_step, " bp")

# ---- Reconstruct pc1 matrix from JSON ---------------------------------------
# theta_pi_local_pca$pc1_loadings is a list of n_windows arrays, each length
# n_samples. We need it as an n_samp × n_win matrix to compute sim_mat from
# columns.
message("[TR_C-wrap] Reconstructing pc1 matrix (n_samp × n_win) ...")
pc1_lol <- js$theta_pi_local_pca$pc1_loadings
if (length(pc1_lol) != n_win) {
  stop("[TR_C-wrap] FATAL: pc1_loadings has ", length(pc1_lol),
       " entries; expected ", n_win)
}
pc1_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                  dimnames = list(sample_order, NULL))
for (wi in seq_len(n_win)) {
  v <- unlist(pc1_lol[[wi]])
  if (length(v) == n_samp) pc1_mat[, wi] <- as.numeric(v)
}

# Same for pc2
pc2_lol <- js$theta_pi_local_pca$pc2_loadings
pc2_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                  dimnames = list(sample_order, NULL))
if (length(pc2_lol) == n_win) {
  for (wi in seq_len(n_win)) {
    v <- unlist(pc2_lol[[wi]])
    if (length(v) == n_samp) pc2_mat[, wi] <- as.numeric(v)
  }
}

n_finite_pc1 <- sum(is.finite(pc1_mat[, 1]))
message("[TR_C-wrap]   pc1 reconstructed: ", n_finite_pc1, " finite at window 1")

# ---- Compute sim_mat = |cor(pc1[, i], pc1[, j])| ----------------------------
# Uses R's cor() with pairwise.complete.obs. Sign-invariant by construction
# (absolute value). Diagonal forced to 1.
message("[TR_C-wrap] Computing sim_mat (n_win × n_win = ",
        n_win, "²) ...")
mem_estimate_mb <- (as.numeric(n_win)^2 * 8) / 1024 / 1024
message("[TR_C-wrap]   estimated peak memory for sim_mat: ",
        round(mem_estimate_mb, 0), " MB")
t_sim <- proc.time()
sim_mat <- suppressWarnings(abs(cor(pc1_mat, use = "pairwise.complete.obs")))
sim_mat[!is.finite(sim_mat)] <- 0
diag(sim_mat) <- 1
storage.mode(sim_mat) <- "double"
message("[TR_C-wrap]   sim_mat computed in ",
        round((proc.time() - t_sim)[3], 1), "s")

# ---- Build D17-compatible precomp.rds ---------------------------------------
message("[TR_C-wrap] Building D17-compatible precomp")
pc_dt <- data.table(
  start_bp = as.integer(window_start_bp),
  end_bp   = as.integer(window_end_bp),
  chrom    = chrom,
  lam_1    = as.numeric(unlist(js$theta_pi_local_pca$lambda_1)),
  lam_2    = as.numeric(unlist(js$theta_pi_local_pca$lambda_2)),
  z        = as.numeric(unlist(js$theta_pi_local_pca$z_profile))
)

# Attach per-sample PC1 / PC2 columns (raw, sign-ambiguous — matching what
# STEP_TR_B emits in v4; v5 retrofit will add aligned versions, but the
# wrapper works on whatever's in the JSON)
for (s_idx in seq_len(n_samp)) {
  s_name <- sample_order[s_idx]
  pc_dt[[paste0("PC_1_", s_name)]] <- as.numeric(pc1_mat[s_idx, ])
  pc_dt[[paste0("PC_2_", s_name)]] <- as.numeric(pc2_mat[s_idx, ])
}

precomp_obj <- list(
  dt           = pc_dt,
  chrom        = chrom,
  n_windows    = n_win,
  sample_names = sample_order,
  source       = "STEP_TR_B JSON → TR_C wrapper",
  generator    = "STEP_TR_C_theta_d17_wrapper.R"
)

temp_dir <- file.path(OUTDIR, paste0(".d17_temp_theta_", chrom))
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
precomp_temp <- file.path(temp_dir, paste0(chrom, ".theta_precomp.slim.rds"))
sim_mat_temp <- file.path(temp_dir, paste0(chrom, ".theta_sim_mat.rds"))

message("[TR_C-wrap] Writing temp precomp: ", precomp_temp)
saveRDS(precomp_obj, precomp_temp)

message("[TR_C-wrap] Writing temp sim_mat: ", sim_mat_temp)
saveRDS(sim_mat, sim_mat_temp)

# ---- Run D17 L1 -------------------------------------------------------------
message("\n[TR_C-wrap] === Invoking D17 L1 ===")
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
message("[TR_C-wrap] L1 cmd: ", l1_cmd)
l1_status <- system(l1_cmd)
if (l1_status != 0) {
  stop("[TR_C-wrap] FATAL: D17 L1 failed with exit code ", l1_status)
}

l1_envelope_path  <- file.path(OUTDIR, paste0(chrom, "_d17L1_envelopes.tsv"))
l1_boundary_path  <- file.path(OUTDIR, paste0(chrom, "_d17L1_boundaries.tsv"))
l1_curve_path     <- file.path(OUTDIR, paste0(chrom, "_d17L1_boundary_score_curve.tsv"))
if (!file.exists(l1_envelope_path)) {
  stop("[TR_C-wrap] FATAL: D17 L1 didn't produce expected envelope output: ",
       l1_envelope_path)
}

n_l1_env <- nrow(fread(l1_envelope_path))
n_l1_bnd <- if (file.exists(l1_boundary_path)) nrow(fread(l1_boundary_path)) else 0
message("[TR_C-wrap] D17 L1 produced: ", n_l1_env, " L1 envelopes, ",
        n_l1_bnd, " L1 boundary peaks")

# ---- Run D17 L2 -------------------------------------------------------------
message("\n[TR_C-wrap] === Invoking D17 L2 ===")
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
message("[TR_C-wrap] L2 cmd: ", l2_cmd)
l2_status <- system(l2_cmd)
if (l2_status != 0) {
  stop("[TR_C-wrap] FATAL: D17 L2 failed with exit code ", l2_status)
}

l2_envelope_path <- file.path(OUTDIR, paste0(chrom, "_d17L2_envelopes.tsv"))
l2_boundary_path <- file.path(OUTDIR, paste0(chrom, "_d17L2_boundaries.tsv"))
n_l2_env <- if (file.exists(l2_envelope_path)) nrow(fread(l2_envelope_path)) else 0
n_l2_bnd <- if (file.exists(l2_boundary_path)) nrow(fread(l2_boundary_path)) else 0
message("[TR_C-wrap] D17 L2 produced: ", n_l2_env, " L2 envelopes, ",
        n_l2_bnd, " L2 boundary peaks")

# ---- Rename outputs to theta-namespaced form --------------------------------
rename_outputs <- function(orig_path) {
  if (!file.exists(orig_path)) return(invisible())
  new_path <- sub("_d17L([12])", paste0("_theta_d17L\\1"), orig_path)
  if (orig_path != new_path) {
    file.rename(orig_path, new_path)
    message("[TR_C-wrap]   renamed: ", basename(orig_path),
            " → ", basename(new_path))
  }
}

message("\n[TR_C-wrap] Renaming outputs to θπ namespace")
rename_outputs(l1_envelope_path)
rename_outputs(l1_boundary_path)
rename_outputs(l1_curve_path)
rename_outputs(l2_envelope_path)
rename_outputs(l2_boundary_path)

if (!KEEP_TEMP) {
  message("[TR_C-wrap] Cleaning temp dir: ", temp_dir)
  unlink(temp_dir, recursive = TRUE)
} else {
  message("[TR_C-wrap] Keeping temp files (--keep_temp): ", temp_dir)
}

message("\n================================================================")
message("[DONE] θπ D17 boundary detection complete for ", chrom)
message("================================================================")
message("  Outputs in: ", OUTDIR)
message("    <chr>_theta_d17L1_envelopes.tsv")
message("    <chr>_theta_d17L1_boundaries.tsv")
message("    <chr>_theta_d17L1_boundary_score_curve.tsv")
message("    <chr>_theta_d17L2_envelopes.tsv")
message("    <chr>_theta_d17L2_boundaries.tsv")
message("")
message("  These ship ALONGSIDE STEP_TR_B v4's |Z|-threshold envelopes")
message("  (which stay in the page-12 JSON's theta_pi_envelopes layer).")
message("  The atlas overlays both on the |Z| panel.")
