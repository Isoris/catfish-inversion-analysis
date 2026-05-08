# ============================================================================
# v5_blocks.R — drop-in code chunks for porting STEP_TR_B v4 → v5
# ============================================================================
# Read SPEC_STEP_TR_B_v5.md first. This file has the actual R code blocks,
# clearly labeled, in the order they should appear in the patched v4 file.
#
# Each BLOCK is a self-contained section. Insert in numerical order at the
# locations described in SPEC §5.
# ============================================================================


# ============================================================================
# BLOCK 0 — Config knobs
# Insert near v4 lines 64–70 (alongside PAD, MAX_K, MIN_L2_WIN etc.)
# ============================================================================

SIM_BAND_HALF        <- 200L
SIM_N_FULL_THRESHOLD <- 6000L

# CLI overrides (insert in v4's existing CLI parsing while-loop ~lines 75–95):
# else if (a == "--sim-band-half"        && i < length(args)) { SIM_BAND_HALF        <- as.integer(args[i + 1]); i <- i + 2L }
# else if (a == "--sim-n-full-threshold" && i < length(args)) { SIM_N_FULL_THRESHOLD <- as.integer(args[i + 1]); i <- i + 2L }


# ============================================================================
# BLOCK 1 — Sim_mat construction
# INSERT AFTER v4 line 343 (after the local-PCA loop, before "L2 envelope").
# ============================================================================

# =============================================================================
# Window×window sim_mat from sign-invariant PC1 correlation
# =============================================================================
# v5 addition (was deferred in v4). sim_mat[i,j] = |cor(pc1_mat[, i], pc1_mat[, j])|.
# The absolute value handles eigenvector sign ambiguity by construction.
# Stored as banded ±SIM_BAND_HALF when n_windows > SIM_N_FULL_THRESHOLD, else
# full upper-triangle. Atlas reader symmetrizes / reflects on load.
# =============================================================================

message("[STEP_TR_B] Building window×window sim_mat (n=", n_win, ")...")
t_sim <- proc.time()

if (n_win <= SIM_N_FULL_THRESHOLD) {

  # ----- FULL dense (upper triangle stored) -----
  sim_mat <- matrix(NA_real_, n_win, n_win)
  for (i in seq_len(n_win)) {
    pc_i <- pc1_mat[, i]
    if (all(!is.finite(pc_i))) next
    for (j in seq.int(i, n_win)) {
      pc_j <- pc1_mat[, j]
      ok <- is.finite(pc_i) & is.finite(pc_j)
      if (sum(ok) < 10L) next
      r <- cor(pc_i[ok], pc_j[ok])
      if (is.finite(r)) sim_mat[i, j] <- abs(r)
    }
  }
  diag(sim_mat) <- 1
  sim_mat_format <- "upper_triangle_float32"
  sim_band <- NULL                                  # not used in this branch

} else {

  # ----- BANDED ±SIM_BAND_HALF (every window stores 2*BAND_HALF + 1 cells) -----
  n_band_cols <- 2L * SIM_BAND_HALF + 1L
  sim_band <- matrix(NA_real_, n_win, n_band_cols)
  # offset 'o' (1..n_band_cols) corresponds to column j = i - BAND_HALF - 1 + o
  for (i in seq_len(n_win)) {
    pc_i <- pc1_mat[, i]
    if (all(!is.finite(pc_i))) next
    j_lo <- max(1L,    i - SIM_BAND_HALF)
    j_hi <- min(n_win, i + SIM_BAND_HALF)
    for (j in j_lo:j_hi) {
      pc_j <- pc1_mat[, j]
      ok <- is.finite(pc_i) & is.finite(pc_j)
      if (sum(ok) < 10L) next
      r <- cor(pc_i[ok], pc_j[ok])
      if (is.finite(r)) {
        o <- j - i + SIM_BAND_HALF + 1L
        sim_band[i, o] <- abs(r)
      }
    }
  }
  # Diagonal is at offset SIM_BAND_HALF + 1; force it to 1.
  sim_band[, SIM_BAND_HALF + 1L] <- 1
  sim_mat_format <- paste0("banded_float32_pm", SIM_BAND_HALF)
  sim_mat <- NULL                                   # not used in this branch
}

message("[STEP_TR_B]   sim_mat built in ", round((proc.time() - t_sim)[3], 1), "s",
        " (format: ", sim_mat_format, ")")


# ============================================================================
# BLOCK 2 — MDS embedding
# INSERT AFTER BLOCK 1.
# ============================================================================

# =============================================================================
# 2D MDS embedding from 1 − sim_mat
# =============================================================================
# Reconstructs a full sim matrix (from upper triangle OR from band, filling
# off-band cells with chrom-median sim) and runs cmdscale(k=2). The full
# matrix is allocated transiently — for very large chroms (n_win > ~20,000)
# this can OOM on tight nodes; if so, downgrade to coarse-grid MDS (see SPEC §5.2).
# =============================================================================

message("[STEP_TR_B] Computing 2D MDS from 1 - sim_mat...")
t_mds <- proc.time()

# Reconstruct the full symmetric sim matrix transiently
if (sim_mat_format == "upper_triangle_float32") {
  sim_mat_full <- sim_mat
  # Mirror upper into lower
  sim_mat_full[lower.tri(sim_mat_full)] <- t(sim_mat_full)[lower.tri(sim_mat_full)]
} else {
  band_median <- median(sim_band, na.rm = TRUE)
  if (!is.finite(band_median)) band_median <- 0.0
  sim_mat_full <- matrix(band_median, n_win, n_win)
  for (i in seq_len(n_win)) {
    j_lo <- max(1L,    i - SIM_BAND_HALF)
    j_hi <- min(n_win, i + SIM_BAND_HALF)
    for (j in j_lo:j_hi) {
      o <- j - i + SIM_BAND_HALF + 1L
      v <- sim_band[i, o]
      if (is.finite(v)) sim_mat_full[i, j] <- v
    }
  }
}
diag(sim_mat_full) <- 1

mds_fit <- tryCatch(
  cmdscale(as.dist(1 - sim_mat_full), k = 2L),
  error = function(e) {
    message("[STEP_TR_B]   cmdscale FAILED: ", e$message)
    NULL
  }
)
rm(sim_mat_full); gc(verbose = FALSE)

if (is.null(mds_fit) || nrow(mds_fit) != n_win) {
  mds1 <- rep(NA_real_, n_win); mds2 <- rep(NA_real_, n_win)
} else {
  mds1 <- mds_fit[, 1]
  mds2 <- mds_fit[, 2]
}

message("[STEP_TR_B]   MDS computed in ", round((proc.time() - t_mds)[3], 1), "s")


# ============================================================================
# BLOCK 3 — Anchor-flip sign alignment
# INSERT AFTER BLOCK 2.
# ============================================================================

# =============================================================================
# Sign-aligned PC1/PC2 loadings via anchor-window flip
# =============================================================================
# Picks the max-|Z| window as the anchor and flips every other window's PC1/PC2
# to maximize correlation with the anchor. Renders the per-sample lines panel
# and PCA scatter coherently across windows.
# =============================================================================

anchor_idx <- which.max(max_abs_z)
if (length(anchor_idx) == 0L || !is.finite(max_abs_z[anchor_idx])) {
  med_z <- median(max_abs_z, na.rm = TRUE)
  if (!is.finite(med_z)) med_z <- 0
  anchor_idx <- which.min(abs(max_abs_z - med_z))[1]
}
anchor_idx <- as.integer(anchor_idx)
message("[STEP_TR_B] Anchor window for sign-alignment: idx=", anchor_idx,
        " (|Z|=", round(max_abs_z[anchor_idx], 2), ")")

pc1_aligned <- pc1_mat
pc2_aligned <- pc2_mat
unflipped_windows <- integer(0)

anchor_pc1 <- pc1_mat[, anchor_idx]
anchor_pc2 <- pc2_mat[, anchor_idx]

for (w in seq_len(n_win)) {
  if (w == anchor_idx) next

  ok1 <- is.finite(anchor_pc1) & is.finite(pc1_mat[, w])
  if (sum(ok1) < 10L) {
    unflipped_windows <- c(unflipped_windows, w)
    next
  }
  r1 <- cor(anchor_pc1[ok1], pc1_mat[ok1, w])
  if (is.finite(r1) && r1 < 0) pc1_aligned[, w] <- -pc1_mat[, w]

  ok2 <- is.finite(anchor_pc2) & is.finite(pc2_mat[, w])
  if (sum(ok2) >= 10L) {
    r2 <- cor(anchor_pc2[ok2], pc2_mat[ok2, w])
    if (is.finite(r2) && r2 < 0) pc2_aligned[, w] <- -pc2_mat[, w]
  }
}

n_unflipped <- length(unflipped_windows)
if (n_unflipped > 0L) {
  message("[STEP_TR_B]   ", n_unflipped, " windows un-flipped (insufficient overlap)")
}


# ============================================================================
# BLOCK 4 — Replacement for v4 lines 604–621 (the theta_pi_local_pca block)
# REPLACE that block with this. Schema bumps to 2.
# ============================================================================

# Layer: theta_pi_local_pca (v5: + sim_mat + MDS + sign-aligned loadings)
# v5 turn: schema_version 1 → 2. Adds sim_mat / mds_coords / pc[12]_loadings_aligned
# / anchor_window_idx, keeps all v4 fields including pc1_loadings (sign-ambiguous,
# kept for backward-compat consumers). Atlas's case-block recognises v2 by the
# new fields' presence.

# Encode sim_mat for JSON. clean_numeric drops to digits=4 (correlation precision
# ~0.0001 is plenty for visualization, halves wire size vs digits=6).
if (sim_mat_format == "upper_triangle_float32") {
  ut_idx     <- which(upper.tri(sim_mat, diag = TRUE), arr.ind = FALSE)
  sim_payload <- clean_numeric(sim_mat[ut_idx], 4)
} else {
  # banded: row-major over [n_win × (2*BAND_HALF + 1)]
  sim_payload <- clean_numeric(as.vector(t(sim_band)), 4)
}

theta_pi_local_pca <- list(
  schema_version       = 2L,
  layer                = "theta_pi_local_pca",
  chrom                = CHROM,
  scale                = if (THETA_GRID_MODE == "native") PESTPG_SCALE else "dosage_grid",
  pad                  = as.integer(PAD),
  n_samples            = as.integer(n_samp),
  n_windows            = as.integer(n_win),
  sample_order         = sample_order,

  # v4 (kept):
  pc1_loadings         = lapply(seq_len(n_win), function(wi) clean_numeric(pc1_mat[, wi], 6)),
  pc2_loadings         = lapply(seq_len(n_win), function(wi) clean_numeric(pc2_mat[, wi], 6)),
  lambda_1             = clean_numeric(lambda_1_vec, 6),
  lambda_2             = clean_numeric(lambda_2_vec, 6),
  lambda_ratio         = clean_numeric(lambda_ratio_vec, 4),
  z                    = clean_numeric(max_abs_z, 4),
  z_profile            = clean_numeric(max_abs_z, 4),
  z_top10_mean         = clean_numeric(top10_abs_z, 4),

  # v5 NEW:
  pc1_loadings_aligned = lapply(seq_len(n_win), function(wi) clean_numeric(pc1_aligned[, wi], 6)),
  pc2_loadings_aligned = lapply(seq_len(n_win), function(wi) clean_numeric(pc2_aligned[, wi], 6)),
  anchor_window_idx    = as.integer(anchor_idx - 1L),       # 0-indexed for atlas JS
  mds_coords           = list(
    mds1 = clean_numeric(mds1, 6),
    mds2 = clean_numeric(mds2, 6)
  ),
  sim_mat_format       = sim_mat_format,
  sim_mat_band_half    = if (sim_mat_format != "upper_triangle_float32") as.integer(SIM_BAND_HALF) else NULL,
  sim_mat_n            = as.integer(n_win),
  sim_mat              = sim_payload,
  `_qc`                = list(
    n_unflipped       = as.integer(n_unflipped),
    unflipped_windows = if (length(unflipped_windows) > 0)
                          as.integer(unflipped_windows - 1L)   # 0-indexed
                        else
                          integer(0)
  )
)


# ============================================================================
# Notes on indexing convention (read carefully)
# ----------------------------------------------------------------------------
# The atlas is JS — 0-indexed. R is 1-indexed.
# The JSON convention (v4 already follows this) is 0-indexed for everything
# the atlas reads positionally:
#   - anchor_window_idx → emit (anchor_idx - 1L)
#   - unflipped_windows → emit (unflipped_windows - 1L)
# The R-internal computations all stay 1-indexed. Only the emitted JSON values
# are 0-indexed. v4 already does this for L2/L1 envelope `win_start` / `win_end`
# (lines 633–648), so the convention is consistent.
# ============================================================================
