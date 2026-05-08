# STEP_TR_B v5 — design spec for porting from v4

**Status:** spec for self-implementation (you port; this doc tells you exactly what to add and where).
**Source:** Apr 30 2026 chat (`d0ac5b6e` — "Theta pi phase 2 preparation and scrubber handoff") § 2.3, plus May 2 chat (`487c7f04` — "Atlas handoff and theta pi analysis preparation") for the local-PCA loop generalization.

---

## 1. Goal

Make `STEP_TR_B_classify_theta.R` mirror `local_PCA_z` exactly: weighted local PCA → sim_mat → MDS → sign-aligned loadings → keep |Z|-based L1/L2 envelopes (still primary).

This walks back v4's "no sim_mat / no MDS" decision (recorded in v4's docstring lines 28–34). The reason for walking it back:

- Per-sample PC1/PC2 loadings are sign-ambiguous between windows. v4 emits raw, sign-ambiguous loadings → atlas's per-sample lines panel renders as random sign-flip noise, PCA scatter clusters jump randomly when scrolling.
- D17b-equivalent cluster-side validation needs per-window MDS coordinates. v4 doesn't emit them.

**What v4 got right that stays.** The 1D |Z|-profile (`max_abs_z`, `top10_abs_z`), the per-window local-PCA loop (heteroscedastic-weighted SVD, λ₁/λ₂), and the |Z|-threshold envelope detection are all correct. None of these change.

**What v5 adds:**

| Addition | Computed how | Stored where in JSON |
|---|---|---|
| `sim_mat[i, j]` | `\|cor(pc1_mat[, i], pc1_mat[, j])\|` over samples; absolute value handles sign ambiguity by construction | `theta_pi_local_pca.sim_mat` |
| `mds_coords` | `cmdscale(as.dist(1 − sim_mat), k = 2)` | `theta_pi_local_pca.mds_coords.{mds1, mds2}` |
| `pc1_loadings_aligned`, `pc2_loadings_aligned` | Anchor-flip: pick max-\|Z\| window as anchor; for every other window, flip its PC vector if `cor(pc1[w], pc1[anchor]) < 0` | `theta_pi_local_pca.pc1_loadings_aligned`, `pc2_loadings_aligned` |
| `anchor_window_idx` | Index of the anchor window | `theta_pi_local_pca.anchor_window_idx` |

`pc1_loadings` (sign-ambiguous) stays in the JSON — keeps backward compatibility for any consumer that already reads it. Atlas's lines panel will read the `_aligned` variant.

---

## 2. Sim_mat sizing decision (was open in prev chat)

Three options were on the table:

- (A) Coarser-grid sim_mat at win50000 (rebuild θπ matrix at coarser scale)
- (B) **Banded sim_mat ±k windows around diagonal** ← chosen
- (C) Int8 quantization of the full matrix

**Choice: banded float32 at ±200 windows.**

Reasoning:
- For LG28 at win10000.step2000 (16,500 windows): banded ±200 = 25 MB raw, fits the 120 MB JSON budget with margin. Full dense float32 is 1.1 GB raw / 260 MB gzip — over budget.
- ±200 windows = ±400 kb of flanking context per focal window at the win10000.step2000 scale. That's well beyond the inversion-flanking radius needed for sign alignment and local PCA-scatter rendering. Far off-diagonal pairs don't drive any consumer of this layer.
- Banded format keeps native float32 precision (no int8 quantization artefacts).
- Storage layout is simple: emit as `sim_mat_band[w][offset]` where `offset = 0..2·BAND_HALF`, and the corresponding column index is `j = w − BAND_HALF + offset`. Atlas reader symmetrizes / reflects on load (cells where `j` is out of range get `null`).

If `n_windows ≤ 6,000` (smaller chroms), the full dense matrix is only ~140 MB raw / ~50 MB gzip and fits comfortably. **For chroms below the threshold, emit full dense.** Use a `sim_mat_format` field to tell the atlas which variant is in the JSON:

- `"upper_triangle_float32"` — full dense, upper triangle only (for small n)
- `"banded_float32_±200"` — banded around diagonal
- `"banded_float32_±N"` — generic name, atlas reads `BAND_HALF` from `sim_mat_band_half`

Threshold: `n_windows ≤ N_FULL_THRESHOLD` (default 6000) → full upper-triangle. Otherwise → banded ±200.

This decision can be revisited later (int8 if 25 MB is still too much, coarser grid if the lines panel doesn't need ±200 of context). Banded is the path of least resistance.

---

## 3. Anchor-flip algorithm

Pseudocode:

```
1. Find anchor window:
     anchor_idx = which.max(max_abs_z)           # already computed in v4
     # Tie-break: if multiple equal-max, pick smallest idx.

2. Initialize aligned matrices:
     pc1_aligned <- pc1_mat   # copy
     pc2_aligned <- pc2_mat

3. For each window w ≠ anchor_idx:
     anchor_pc1 = pc1_mat[, anchor_idx]
     w_pc1      = pc1_mat[, w]
     # Use only the samples with finite values in BOTH (intersection):
     ok = is.finite(anchor_pc1) & is.finite(w_pc1)
     if (sum(ok) < 10) {
         # Not enough overlap to decide — leave un-flipped, mark in QC.
         next
     }
     r1 = cor(anchor_pc1[ok], w_pc1[ok])
     if (r1 < 0) pc1_aligned[, w] <- -pc1_mat[, w]

     # PC2 sign aligned independently to its own anchor (the same anchor window's PC2):
     r2 = cor(pc2_mat[ok, anchor_idx], pc2_mat[ok, w])
     if (is.finite(r2) && r2 < 0) pc2_aligned[, w] <- -pc2_mat[, w]
```

Notes:
- Independent sign-alignment for PC1 and PC2 is correct — they're orthogonal, their signs are independently arbitrary.
- The anchor's own row stays un-flipped by construction.
- Use the intersection of finite samples for the correlation. Don't impute, don't drop windows from the matrix — just leave un-flipped windows where the overlap is too thin (record their indices in a small `_qc.unflipped_windows` array if you want).
- Do NOT use `abs(cor) < threshold` to gate flips — small \|cor\| means weak alignment evidence, but flipping based on the sign of a small correlation is still correct (it makes the sign coherent with the anchor; a near-zero correlation just means the flip choice doesn't matter much).

---

## 4. Sim_mat construction algorithm

```
N <- n_win
BAND_HALF <- 200L                                  # configurable, default 200
N_FULL_THRESHOLD <- 6000L                          # configurable

if (N <= N_FULL_THRESHOLD) {
    # FULL dense, upper triangle only
    sim_mat <- matrix(NA_real_, N, N)
    for (i in seq_len(N)) {
        for (j in seq.int(i, N)) {
            ok <- is.finite(pc1_mat[, i]) & is.finite(pc1_mat[, j])
            if (sum(ok) < 10) next
            r <- cor(pc1_mat[ok, i], pc1_mat[ok, j])
            sim_mat[i, j] <- abs(r)
        }
    }
    diag(sim_mat) <- 1
    sim_mat_format <- "upper_triangle_float32"
} else {
    # BANDED, ±BAND_HALF
    sim_band <- matrix(NA_real_, N, 2L * BAND_HALF + 1L)
    # column 'offset' o = 0..2·BAND_HALF corresponds to j = i - BAND_HALF + o
    for (i in seq_len(N)) {
        for (o in seq_len(2L * BAND_HALF + 1L)) {
            j <- i - BAND_HALF - 1L + o
            if (j < 1L || j > N) next                 # leave NA → JSON null
            ok <- is.finite(pc1_mat[, i]) & is.finite(pc1_mat[, j])
            if (sum(ok) < 10) next
            r <- cor(pc1_mat[ok, i], pc1_mat[ok, j])
            sim_band[i, o] <- abs(r)
        }
    }
    sim_mat_format <- paste0("banded_float32_±", BAND_HALF)
}
```

**Performance.** Naïve loop is O(N² × n_samp). For LG28 with N=16,500, n_samp=226, the banded version does 16500 × 401 × 226 = 1.5 × 10⁹ ops. At ~10⁸/s for `cor()` calls in R (which is conservative — R's `cor` is C-vectorized, faster than this), that's ~15 s. **You can vectorize it for a 10× speedup** by computing `cor(pc1_mat)` (the full N×N correlation) chunk-wise — `cor(pc1_mat[, i:(i+99)])` returns a 100×N block — and reading off the diagonal band. Don't optimize until profiling shows it's a bottleneck; the naïve loop should be acceptable at ~30–60 s/chrom.

**Equivalent computation (vectorized chunk approach if needed):**

```r
# stdize pc1_mat columns once (treating columns as n_samp-length vectors)
# - centre each column by its mean over finite samples
# - scale each column so sum-of-squares = 1
# Then sim_mat[i,j] = abs( crossprod(pc1_std[, i], pc1_std[, j]) )
# which lets you compute the full correlation block by `crossprod(pc1_std)`.
# Handle NAs by zero-imputing AFTER centering & before normalizing
# (this is what `cor(use="pairwise.complete.obs")` does column-wise; for
# the full-matrix shortcut you do "complete.obs" which drops samples with
# any NA in any window — usually fine for theta_mat).
```

If you want exact `pairwise.complete.obs` semantics, the naïve loop is the simplest path. If you want speed, use `complete.obs` + `crossprod`. Up to you.

---

## 5. Diff plan against v4 — exact insertion points

The v4 file is 715 lines. The v5 retrofit adds **three blocks** and modifies **one block**:

### 5.1 Add: sim_mat construction (NEW BLOCK)

**Insert after line 343** (after the local-PCA loop, before "L2 envelope detection from contiguous high-|Z| runs"). About 60 lines.

This is a new section header:
```r
# =============================================================================
# Window×window sim_mat from sign-invariant PC1 correlation
# =============================================================================
# v5 addition (was deferred in v4). sim_mat[i,j] = |cor(pc1_mat[, i], pc1_mat[, j])|.
# The absolute value handles eigenvector sign ambiguity by construction.
# Stored as banded ±BAND_HALF when n_windows > N_FULL_THRESHOLD, else full
# upper-triangle. Atlas reader symmetrizes / reflects on load.
# =============================================================================
```

Followed by the construction code from §4 above. Add config knobs near the top of the script (where `PAD`, `MAX_K` etc. are declared at lines 64–70):

```r
SIM_BAND_HALF        <- 200L
SIM_N_FULL_THRESHOLD <- 6000L
```

Plus CLI overrides (the same `--sim-band-half` / `--sim-n-full-threshold` args, mirroring v4's CLI parsing block lines 75–95).

### 5.2 Add: MDS embedding (NEW BLOCK)

**Insert after sim_mat construction.** About 25 lines.

```r
# =============================================================================
# 2D MDS embedding from 1 − sim_mat
# =============================================================================
message("[STEP_TR_B] Computing 2D MDS from 1 - sim_mat ...")
t3 <- proc.time()

# For banded sim_mat, fill missing cells with their column-pairwise mean
# before cmdscale (cmdscale needs a full distance matrix). Alternative:
# convert NA → bg_distance (e.g. 0.5). The banded part captures the local
# geometry; far-off-diagonal cells just need to be "neutral" so they don't
# dominate.
if (sim_mat_format == "upper_triangle_float32") {
    # Symmetrize:
    sim_mat_full <- sim_mat
    sim_mat_full[lower.tri(sim_mat_full)] <- t(sim_mat_full)[lower.tri(sim_mat_full)]
} else {
    # Reconstruct full matrix from band, fill remainder with chrom-median sim
    sim_mat_full <- matrix(NA_real_, N, N)
    for (i in seq_len(N)) {
        for (o in seq_len(2L * SIM_BAND_HALF + 1L)) {
            j <- i - SIM_BAND_HALF - 1L + o
            if (j >= 1L && j <= N) sim_mat_full[i, j] <- sim_band[i, o]
        }
    }
    sim_mat_full[is.na(sim_mat_full)] <- median(sim_band, na.rm = TRUE)
}
diag(sim_mat_full) <- 1

mds_fit <- tryCatch(
    cmdscale(as.dist(1 - sim_mat_full), k = 2L),
    error = function(e) { message("[STEP_TR_B]   cmdscale failed: ", e$message); NULL }
)

if (is.null(mds_fit)) {
    mds1 <- rep(NA_real_, N); mds2 <- rep(NA_real_, N)
} else {
    mds1 <- mds_fit[, 1]; mds2 <- mds_fit[, 2]
}

message("[STEP_TR_B]   MDS computed in ", round((proc.time() - t3)[3], 1), "s")
```

**Memory note.** `sim_mat_full` at LG28 is 1 GB. If that's too much for the LANTA node in the SLURM array (32 GB / task currently), use a coarse-grid MDS instead: bin the windows into 100-wide bins, compute the bin-level sim_mat (165×165 for LG28), run MDS on that, and bilinearly interpolate back to per-window mds1/mds2. The atlas only uses MDS coords for the per-window scatter, which doesn't need full per-window resolution. **Defer this decision to the first dry-run on LG28 — try the full approach first, downgrade to coarse if you OOM.**

### 5.3 Add: anchor-flip sign alignment (NEW BLOCK)

**Insert after MDS.** About 30 lines.

```r
# =============================================================================
# Sign-aligned PC1/PC2 loadings via anchor-window flip
# =============================================================================
# v5 addition. Picks the max-|Z| window as the anchor and flips every other
# window's PC1/PC2 to maximize correlation with the anchor. This makes the
# atlas's per-sample lines panel and PCA scatter render coherently across
# windows instead of showing arbitrary sign flips.
# =============================================================================

anchor_idx <- which.max(max_abs_z)
if (length(anchor_idx) == 0L || !is.finite(max_abs_z[anchor_idx])) {
    # Fallback: median-|Z| window
    anchor_idx <- which.min(abs(max_abs_z - median(max_abs_z, na.rm = TRUE)))
}
anchor_idx <- as.integer(anchor_idx)
message("[STEP_TR_B] Anchor window for sign-alignment: idx=", anchor_idx,
        " (|Z|=", round(max_abs_z[anchor_idx], 2), ")")

pc1_aligned <- pc1_mat
pc2_aligned <- pc2_mat
n_unflipped <- 0L
unflipped_windows <- integer(0)

anchor_pc1 <- pc1_mat[, anchor_idx]
anchor_pc2 <- pc2_mat[, anchor_idx]

for (w in seq_len(n_win)) {
    if (w == anchor_idx) next
    ok1 <- is.finite(anchor_pc1) & is.finite(pc1_mat[, w])
    if (sum(ok1) < 10L) {
        n_unflipped <- n_unflipped + 1L
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

if (n_unflipped > 0L) {
    message("[STEP_TR_B]   ", n_unflipped, " windows un-flipped (insufficient overlap)")
}
```

### 5.4 Modify: JSON `theta_pi_local_pca` block

**Modify lines 604–621** (the `theta_pi_local_pca <- list(...)` block).

Add four fields, schema bumps to v5:

```r
theta_pi_local_pca <- list(
    schema_version       = 2L,                      # was 1L; bump for v5 fields
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
    anchor_window_idx    = anchor_idx,
    mds_coords           = list(
        mds1 = clean_numeric(mds1, 6),
        mds2 = clean_numeric(mds2, 6)
    ),
    sim_mat_format       = sim_mat_format,
    sim_mat_band_half    = if (sim_mat_format != "upper_triangle_float32") as.integer(SIM_BAND_HALF) else NULL,
    sim_mat_n            = as.integer(n_win),
    sim_mat              = if (sim_mat_format == "upper_triangle_float32") {
        # upper-triangle: emit as flat row-major over upper triangle including diag
        ut_idx <- which(upper.tri(sim_mat, diag = TRUE), arr.ind = FALSE)
        clean_numeric(sim_mat[ut_idx], 4)
    } else {
        # banded: emit row-major over the band matrix [N × (2·BAND_HALF + 1)]
        # row-major: row 0 = window 0's band, etc.
        clean_numeric(as.vector(t(sim_band)), 4)
    },
    `_qc`                = list(
        n_unflipped       = as.integer(n_unflipped),
        unflipped_windows = if (length(unflipped_windows) > 0) as.integer(unflipped_windows) else NULL
    )
)
```

Note: `digits = 4` for sim_mat (correlations don't need more precision than 0.0001 for visual rendering, and 4-digit cuts JSON size in half vs 6-digit).

### 5.5 Modify: top-level `_layers_present` (no change — unchanged)

The `_layers_present` array doesn't change; `theta_pi_local_pca` is still listed. The atlas's case-block reads schema_version=2 and detects the new fields by their presence (turn-114 dual-tolerance pattern).

### 5.6 Optional: write a small `_qc` summary at the top level

If you want post-hoc visibility into the sign-alignment quality across chroms, add at line ~692 (next to `_generated_at`):

```r
`_qc_signalign` = list(
    anchor_window_idx = anchor_idx,
    n_windows         = as.integer(n_win),
    n_unflipped       = as.integer(n_unflipped)
)
```

Not required.

---

## 6. Testing the port

Three checks before submitting the 28-chrom array:

### 6.1 Run on LG28 only (the validated prototype)

```bash
$RSCRIPT STEP_TR_B_classify_theta.R --chrom C_gar_LG28
```

Expected runtime increase vs v4: roughly +30–60 s (sim_mat construction + MDS).

### 6.2 Eyeball the JSON

```bash
jq '.schema_version' <out.json>          # should be 2
jq '._layers_present' <out.json>          # unchanged
jq '.theta_pi_local_pca | keys' <out.json>
# Expected new keys:
#   anchor_window_idx, mds_coords, pc1_loadings_aligned, pc2_loadings_aligned,
#   sim_mat, sim_mat_band_half, sim_mat_format, sim_mat_n, _qc
jq '.theta_pi_local_pca.anchor_window_idx' <out.json>
jq '.theta_pi_local_pca.mds_coords | length' <out.json>     # should be 2 (mds1, mds2)
jq '.theta_pi_local_pca.sim_mat | length' <out.json>
# For banded ±200 with N=16500: expect 16500 × 401 = 6,616,500
jq '.theta_pi_local_pca._qc.n_unflipped' <out.json>
# Should be 0–10. If much higher, investigate (probably means many windows
# have <10 finite-overlap samples with the anchor).
```

### 6.3 Sign-alignment sanity

```r
# In R, after running TR_B once on LG28:
o <- jsonlite::fromJSON("<chr>_phase2_theta.json", simplifyVector = FALSE)
pca <- o$theta_pi_local_pca
anchor <- pca$anchor_window_idx + 1L                     # R is 1-indexed
n     <- pca$n_windows
n_samp <- pca$n_samples

aligned <- do.call(cbind, lapply(pca$pc1_loadings_aligned,
                                 function(x) unlist(x, use.names = FALSE)))
unaligned <- do.call(cbind, lapply(pca$pc1_loadings,
                                   function(x) unlist(x, use.names = FALSE)))
# Correlations to anchor should be POSITIVE for aligned (after flip):
ank <- aligned[, anchor]
cors_after  <- sapply(seq_len(n), function(w) {
    ok <- is.finite(ank) & is.finite(aligned[, w])
    if (sum(ok) < 10) return(NA_real_)
    cor(ank[ok], aligned[ok, w])
})
cat("Median post-flip cor to anchor:", median(cors_after, na.rm = TRUE), "\n")
cat("Fraction negative post-flip:   ", mean(cors_after < 0, na.rm = TRUE), "\n")
# Expected: median > 0, fraction negative < 0.05
```

If the post-flip median correlation is negative or close to zero, the anchor flip didn't work — debug.

---

## 7. What does NOT change in the port

- The L1/L2 envelope detection logic (lines 348–429 of v4): unchanged. Stays primary.
- The `theta_pi_envelopes` JSON layer (lines 624–651): unchanged.
- The `tracks` layer (lines 655–674): unchanged.
- The `theta_pi_per_window` layer (lines 568–595): unchanged.
- STEP_TR_A, STEP_TR_C, STEP_TR_D, the launcher, the verifier, chrom.list, the config: all unchanged.

This means the rest of the bundle (incl. our `THETA_CONFIG_DIR` launcher fix) is independent of this porting work.

---

## 8. The chunk you have to start from (May 2 chat)

The May 2 chat has the **NPC-generalized local PCA loop** that goes a step beyond v4's hardcoded `nu = 2`. v4's loop already does what you need at NPC=2, so the May 2 generalization is **optional** — only worth pulling in if you want to expose PC3/PC4/PC5 to the atlas later. For the v5 retrofit goal (sim_mat + MDS + sign-alignment), v4's existing PCA loop is sufficient.

If you do want NPC>2 support, replace v4's lines 290–337 with the May 2 chunk:

```r
NPC <- 2L                                     # configurable; v5 keeps default 2

pc_mats <- vector("list", NPC)
for (k in seq_len(NPC)) {
    pc_mats[[k]] <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                           dimnames = list(sample_order, NULL))
}
lambda_mat       <- matrix(NA_real_, nrow = n_win, ncol = NPC)
lambda_ratio_vec <- rep(NA_real_, n_win)

n_sites_chrom_median <- median(n_sites_mat, na.rm = TRUE)
if (!is.finite(n_sites_chrom_median) || n_sites_chrom_median <= 0) {
    n_sites_chrom_median <- 1
}

for (wi in seq_len(n_win)) {
    lo <- max(1L, wi - PAD)
    hi <- min(n_win, wi + PAD)
    block <- theta_mat[, lo:hi, drop = FALSE]
    ok <- complete.cases(block)
    if (sum(ok) < max(20L, ncol(block) + 2L)) next
    block_ok <- block[ok, , drop = FALSE]

    n_focal <- n_sites_mat[ok, wi]
    n_focal[!is.finite(n_focal) | n_focal <= 0] <- 1L
    weights <- sqrt(n_focal / n_sites_chrom_median)
    block_w <- block_ok * weights

    centred <- sweep(block_w, 2, colMeans(block_w), FUN = "-")
    nu_req <- min(NPC, ncol(centred))
    sv <- tryCatch(svd(centred, nu = nu_req, nv = 0), error = function(e) NULL)
    if (is.null(sv)) next

    for (k in seq_len(min(nu_req, ncol(sv$u)))) {
        pc_mats[[k]][ok, wi] <- as.numeric(sv$u[, k])
    }

    d2 <- sv$d^2
    for (k in seq_len(min(NPC, length(d2)))) {
        lambda_mat[wi, k] <- d2[k]
    }
    if (length(d2) >= 2 && d2[2] > 0) {
        lambda_ratio_vec[wi] <- d2[1] / d2[2]
    }
}

# Convenience aliases for back-compat with the rest of the script:
lambda_1_vec <- lambda_mat[, 1]
lambda_2_vec <- if (NPC >= 2) lambda_mat[, 2] else rep(NA_real_, n_win)
pc1_mat      <- pc_mats[[1]]
pc2_mat      <- if (NPC >= 2) pc_mats[[2]] else matrix(NA_real_, nrow = n_samp,
                                                       ncol = n_win,
                                                       dimnames = list(sample_order, NULL))
```

Strictly optional for v5 retrofit; recommended for future NPC>2 support.

---

## 9. Effort estimate

- Sim_mat construction block: ~30 min to type out, ~30 min to debug on LG28.
- MDS block: ~15 min, may need 1 retry if memory issues on full matrix.
- Anchor-flip block: ~20 min.
- JSON layer modification: ~10 min.
- Atlas JSON sanity checks: ~30 min.

**Total: ~2-3 hours of focused work.** All blocks are self-contained additions; the only modification to existing code is the `theta_pi_local_pca` JSON list constructor.

---

## 10. References

- v4 file: `STEP_TR_B_classify_theta.R` lines 1–715 (already on LANTA, in the bundle).
- Apr 30 chat (`d0ac5b6e`) §2.3 — design rationale for v4 → v5.
- May 2 chat (`487c7f04`) — NPC-generalized local-PCA loop (optional pull).
- The dosage scrubber's analogous code in `local_PCA_z` lives under
  `phase_2_discovery/2a_local_pca/` on LANTA: `STEP_C01a_snake1_precompute.R` for the sim_mat construction reference,
  `STEP_C04c_ghsl_local_pca.R` for the GHSL parallel that uses the same `cor(pc1)`-based sim_mat
  and anchor-flip pattern (verified working in turn 2 of the Apr 30 plan).
