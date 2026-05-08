# HANDOFF — `2f_theta_discovery/` (path 2, awaiting refactor)

> **Status: working code, not yet refactored to the consolidated 01–09 skeleton.**
> This handoff documents what's here, how each script functions today, and
> how it would map into the canonical skeleton if/when refactored. Read this
> alongside `README_theta_pi_scaling.md` (operational notes) and
> `../../HANDOFF.md` (path 1 / canonical skeleton).
>
> Codebase: `inversion-popgen-toolkit` v8.5. The four active scripts in this
> folder ARE working on LANTA — don't break them. The refactor decision
> (do it or not) is a future-session call; see §5 below.

---

## 1. What this pipeline does

Path 2 = "θπ" — uses **per-sample per-window pairwise nucleotide diversity
θπ** (computed by ANGSD's `-doThetas` + `thetaStat do_stat` per sample) to
detect inversions and other structural rearrangements.

The core observation: inside an inverted region, both haplotypes accumulate
mutations independently (no recombination crosses the breakpoint), driving
**within-individual** sequence diversity up. θπ is a direct, sign-stable
quantitative measure of that — higher θπ means more diversity, lower means
less, with no eigenvector flip ambiguity.

This is **the simplest of the three paths algorithmically** because:
- θπ is sign-stable by construction → no need for sim_mat or MDS to
  recover sign-invariant similarity. The 1D |Z| profile on θπ values
  IS the canonical signal.
- The detection collapses to "find contiguous high-|Z| runs in 1D" —
  much simpler than path 1's 2D cross-block scan on a sim_mat.
- Local PCA on the θπ matrix is computed but plays a **secondary** role
  (descriptive only — pc1/pc2 loadings, λ ratio for atlas display).

The pipeline computes:
1. Per-sample × per-window θπ matrix from ANGSD pestPG files (heavy)
2. Per-window |Z| profile + window-local PCA (light, runs in seconds)
3. L1/L2 envelopes from contiguous high-|Z| runs (PRIMARY detection)
4. Optional D17 boundary detection on a sim_mat computed from local PCA
5. Page-12 atlas JSON

---

## 2. Current file inventory

```
2f_theta_discovery/
├── 00_theta_config.sh                     env vars (SAMPLE_LIST, paths, scale)
├── 00_sanity_check_pestPG_scaling.sh      preflight on pestPG existence
├── verify_theta_setup.sh                  7-check preflight (run before launching)
├── chrom.list                             28-line chromosome list for SLURM array
├── STEP_TR_A_compute_theta_matrices.R     heavy: pestPG → long TSV (~3-5 min/chr)
├── STEP_TR_B_classify_theta.R             light: TSV → JSON, ALL-IN-ONE (~30 s/chr)
├── STEP_TR_C_theta_d17_wrapper.R          optional: D17 detector on θπ pc1 sim_mat
├── STEP_TR_D_augment_theta_json.R         optional: post-process JSON to add D17 layer
├── LAUNCH_TR_theta_pi.slurm               array launcher (TR_A + TR_B per chrom)
└── README_theta_pi_scaling.md             operational notes (don't replace)
```

---

## 3. What each script does, in detail

### 3.1 `STEP_TR_A_compute_theta_matrices.R` — heavy precomputation

**Per chromosome.** Reads ANGSD pestPG files (one per sample × per chrom)
at `PESTPG_SCALE` (default `win10000.step2000`), emits a long-format TSV
with one row per (sample, window).

This is the **feature-build step** — analogous to path 1's
`01a_beagle_to_dosage.py`. It produces the θπ matrix that everything
downstream consumes.

**Window-grid policy** (`THETA_GRID_MODE` in `00_theta_config.sh`):
- `native` (default): use the pestPG grid as-is. On LG28 at
  `win10000.step2000` this yields ~16,500 windows per chrom. The atlas
  uses Int32Array lookup tables to translate between the dosage-grid
  cursor (state.cur) and the θπ-grid cursor (state.cur_thpi).
- `dosage`: nearest-midpoint join θπ values onto path 1's variable-bp
  window grid. Fallback if native-grid loadings prove too noisy.

**Inputs:**
- `${PESTPG_DIR}/${SAMPLE}.${PESTPG_SCALE}.pestPG` (226 samples × 28 chroms)
- `${SAMPLE_LIST}` (one CGA id per line)
- (dosage mode only) `${DOSAGE_WIN_BED_DIR}/${CHROM}/windows.bed`

**Outputs:**
- Native mode: `${THETA_TSV_DIR}/theta_native.<CHROM>.<SCALE>.tsv.gz`
- Dosage mode: `${THETA_TSV_DIR}/theta_dgrid.<CHROM>.tsv.gz`
- Columns: `sample chrom window_idx start_bp end_bp theta_pi n_sites`
- Walltime: ~3–5 min per chromosome on LANTA.

### 3.2 `STEP_TR_B_classify_theta.R` — the everything-script

**~30 s walltime per chromosome.** This is the script that does almost all
of the path-2 work in one pass. It:

1. Loads the long TSV from TR_A → builds samples × windows θπ matrix.
2. Computes per-window population metrics:
   median, MAD, IQR, **|Z| from per-sample deviations against chromosome
   baseline**.
3. Runs **window-by-window local PCA** over a (2·pad + 1) window
   neighborhood: PC1/PC2 loadings per sample, λ₁/λ₂ ratio, 1D-ness flag.
4. Calls **L2 envelopes** from contiguous |Z| > threshold runs (default
   threshold 2.5, min length 5 windows, merge gap 3).
5. Optionally calls **per-interval karyotype classes** via k-means on
   interval-mean θπ when `--intervals` is supplied (k=2..max_k with
   silhouette selection; mirrors GHSL's C04b output).
6. **Emits the page-12 atlas JSON directly** with these layers:
   - `theta_pi_per_window` — samples × windows θπ + n_sites
   - `theta_pi_local_pca` — PC1/PC2 loadings, λ ratios, |Z| profile
   - `theta_pi_envelopes` — L2 + L1 envelope coordinates (PRIMARY)
   - `tracks` — theta_pi_median, theta_pi_z, lambda_ratio

**This single script does the work of path 1's steps 01b + 01c + 02a + 02b
+ 03 + 04 + 06 + 08b combined.** Possible because:
- θπ is sign-stable → no MDS needed (path 1 needs MDS to denoise sign-
  ambiguous PCs)
- Detection works on the 1D |Z| profile → no sim_mat needed for primary
  detection (path 1 needs sim_mat for the 2D cross-block boundary scan)

Plot scripts (path 1's 05/07) don't exist for path 2 — the page-12 atlas
viewer renders directly from the JSON.

**Outputs:**
- `${JSON_OUT_DIR}/<CHROM>/<CHROM>_phase2_theta.json` (page-12 JSON)
- (intermediates inside the JSON, no separate TSVs)

### 3.3 `STEP_TR_C_theta_d17_wrapper.R` — optional D17 detector

**Optional.** Runs path 1's D17 boundary detector on a sim_mat computed
from TR_B's pc1_loadings — a "what if we used path 1's edge-based detector
on θπ?" exploration, for symmetry with GHSL's C04d.

Why optional: TR_B's |Z|-threshold detector IS the canonical PRIMARY for
θπ. There's no separate biologically-calibrated detector for θπ
(unlike GHSL, which has C04b's PASS-runs as biological PRIMARY). So D17
on θπ is a methodological cross-check, not a required ingredient.

**This wrapper currently invokes path 1's D17 scripts BY OLD NAME**:
- `STEP_D17_multipass_L1_only_v7.R` ← now `04_detect_L1_localpca_zblocks.R`
- `STEP_D17_multipass_L2_v8.R`      ← now `06_detect_L2_localpca_zblocks.R`

**This is broken in the consolidated layout** — same problem as GHSL's
C04d wrapper. Fix by updating the wrapper to call the new script names.

**Memory cost flagged in header:** at `win10000.step2000` ≈ 16,500 windows,
sim_mat ≈ 1.1 GB float64. Fine cluster-side; NEVER serialized to JSON
(only the small TSVs ship to browser).

**Outputs** (use old `_d17L*` suffix, will need renaming):
- `<chr>_theta_d17L1_envelopes.tsv` etc.

### 3.4 `STEP_TR_D_augment_theta_json.R` — JSON post-processor

**Optional.** Reads TR_B's `<chr>_phase2_theta.json` + TR_C's D17 TSVs,
adds a `theta_d17_envelopes` JSON layer, writes the augmented JSON in
place (or to a separate path).

This exists as a separate script (rather than being merged into TR_B)
because TR_B v4 already ships and works; TR_C is an additional path that
runs after TR_B. Keeping the dependency graph clean:

```
TR_A → TR_B → JSON(v1)
              ↓
TR_A → TR_B JSON → TR_C wrapper → D17 TSVs
                                       ↓
TR_B JSON + D17 TSVs → TR_D → JSON(v2 augmented)
```

---

## 4. Mapping to the canonical 01–09 skeleton

This pipeline does NOT decompose into the canonical skeleton — it's
**deliberately compressed** because θπ doesn't need most of the path-1
machinery.

| Canonical step                  | θπ equivalent                            |
|---------------------------------|------------------------------------------|
| 01a feature build               | `STEP_TR_A_compute_theta_matrices.R`     |
| 01b local PCA compute           | part of `STEP_TR_B_classify_theta.R`     |
| 01c local PCA merge             | (not needed — TR_B is per-chr)           |
| 02a MDS compute                 | (not needed — θπ is sign-stable)         |
| 02b MDS merge                   | (not needed)                             |
| 03 precompute (sim_mat)         | (not needed for primary detection)       |
| 04 detect L1                    | part of `STEP_TR_B` (|Z|-threshold)      |
| 05 plot L1                      | (not implemented — atlas renders JSON)   |
| 06 detect L2                    | part of `STEP_TR_B` (|Z|-threshold)      |
| 07 plot L2                      | (not implemented)                        |
| 08a sample metadata             | (consumes `_shared/sample_metadata.tsv`) |
| 08b atlas JSON                  | part of `STEP_TR_B` (emits directly)     |
| **Optional D17 cross-check**    | `STEP_TR_C` + `STEP_TR_D`                |

The mismatch: **TR_B does eight canonical steps in one script**. The
compression is justified — θπ doesn't need MDS or sim_mat for primary
detection, and the per-step splits in path 1 exist mainly because path 1's
intermediate artifacts (window_pca.rds, mds.rds, sim_mat.rds) are
**reused across many parameter sweeps**, while TR_B's output (the JSON)
is the only consumer of its intermediates.

---

## 5. Refactoring decision

**Option A — Refactor to canonical 01–09 skeleton.**

This would mean splitting `STEP_TR_B` into:
- `01b_local_pca_compute_thetapi.R` → computes per-window PC1/PC2 + saves precomp.rds
- `02_compute_z_profile_thetapi.R` → computes |Z| profile from samples × θπ matrix
  (no MDS — θπ is sign-stable)
- `04_detect_L1_thetapi.R` → finds contiguous high-|Z| L1 runs
- `06_detect_L2_thetapi.R` → finds L2 sub-runs inside each L1
- `08b_export_atlas_json_thetapi.R` → packs everything into JSON

Pros:
- Three discovery paths look symmetric on disk
- Each step's outputs (precomp, L1, L2 TSVs) become inspectable side-files
  rather than intermediate variables inside one script
- Power-user parameter sweeps on just step 04 or 06 don't require re-running
  the whole TR_B (saves ~20 s per iteration — modest gain)
- Optional D17 wrapper (TR_C) becomes redundant — you'd just have a `--mode d17`
  switch on `04_detect_L1_thetapi.R` to call the path-1 D17 logic instead of
  the |Z|-threshold logic

Cons:
- Significant surgery on a working ~30-second pipeline
- The compression in TR_B reflects a real architectural insight (θπ doesn't
  need sim_mat) — splitting it back out re-introduces overhead the design
  was avoiding
- Risk of regression on a pipeline that ships and works

**Option B — Lighter touch: keep current architecture, just rewire paths.**

What this looks like:
1. Move the folder to `local_PCA_thetapi/` (sibling of `local_PCA_z/`).
2. Rename the 4 R scripts to canonical names:
   - `STEP_TR_A_compute_theta_matrices.R` → `01_compute_theta_matrices.R`
   - `STEP_TR_B_classify_theta.R` → `02_classify_thetapi.R` (the everything-script)
   - `STEP_TR_C_theta_d17_wrapper.R` → `03_optional_d17_thetapi.R`
   - `STEP_TR_D_augment_theta_json.R` → `04_augment_atlas_json_thetapi.R`
3. Update TR_C wrapper to call new path-1 D17 script names
   (`04_detect_L1_localpca_zblocks.R` etc.).
4. Update output naming: `<chr>_theta_d17L1_envelopes.tsv` →
   `<chr>.thetapi_L1_envelopes.tsv`.
5. Have TR_B consume `_shared/sample_metadata.tsv` for the `samples`
   sub-layer instead of a separate `--samples` flag.
6. Update launchers + config to write to `path_localpca_thetapi/` scratch tree.

Pros:
- Achievable in ~30 minutes
- No risk of breaking working pipeline
- Names + scratch tree look like the others
- Honest about the architectural difference (θπ is naturally compressed)

Cons:
- 4 scripts instead of canonical 11 — visibly different from path 1
- The everything-script (`02_classify_thetapi.R`) is large and mixes
  concerns; future contributors might be confused why path 1 has 11 scripts
  and path 2 has 4

**Recommended: Option B.** The compression in `STEP_TR_B` is a feature,
not a bug — forcing it into 11 scripts to match path 1 would be
mechanical symmetry without algorithmic justification. Document the
asymmetry honestly in the per-path README.

---

## 6. Required input data on LANTA

When refactored to scratch-tree convention, this pipeline reads:
- `<scratch>/03_pestPG/<sample>.<scale>.pestPG` — ANGSD `-doThetas` output,
  one file per sample per chrom. **Currently produced by a pipeline OUTSIDE
  this toolkit** (`MODULE_3 / 02_run_heterozygosity.sh`). Not reproduced
  by `local_PCA_z`.
- `<scratch>/_shared/sample_list.tsv` — one CGA per line (or just use
  `_shared/sample_metadata.tsv` and project the `cga` column).

Outputs (when refactored) should land in:
- `<scratch>/path_localpca_thetapi/01_theta_matrices/theta_native.<chr>.tsv.gz`
- `<scratch>/path_localpca_thetapi/02_atlas_json/<chr>.atlas_thetapi.json`
- (if D17 used) `<scratch>/path_localpca_thetapi/03_d17_optional/<chr>.thetapi_L1_envelopes.tsv` etc.

---

## 7. What's NOT broken (don't fix)

- The "θπ doesn't need MDS or sim_mat for primary detection" architectural
  decision. It's correct.
- `STEP_TR_B`'s |Z|-threshold L1/L2 detector — defaults
  `Z_THRESHOLD=2.5`, `MIN_L2_WIN=5`, `MERGE_GAP=3`, `PAD=1` reflect tuning
  history.
- The native-grid window policy (16,500 windows per chrom on LG28) — the
  atlas Int32Array cursor lookup tables are designed for it.
- `verify_theta_setup.sh` 7-check preflight — saves debugging time on real
  runs.

---

## 8. Run order (current, working)

```bash
cd inversion_modules/phase_2_discovery/2f_theta_discovery
mamba activate assembly

# Preflight (~30 s)
bash verify_theta_setup.sh

# Single launcher does TR_A → TR_B per chrom
sbatch --array=1-28 LAUNCH_TR_theta_pi.slurm chrom.list

# Optional: add D17 cross-check layer to JSON
# (per chrom, after TR_B output exists)
Rscript STEP_TR_C_theta_d17_wrapper.R --chrom <CHR> ...
Rscript STEP_TR_D_augment_theta_json.R --theta_json ... --d17_l1_env ...
```

Per-chrom walltime: TR_A ≈ 3–5 min, TR_B ≈ 30 s. With 28 chroms in
parallel, total wall ≈ 5 min.

---

## 9. Open questions for the refactor session

1. **Native vs dosage window grid — is the dual mode still needed?** The
   header says dosage mode is a fallback if native proves too noisy. If
   native has been working in production, simplify by retiring dosage
   mode (one branch fewer to maintain).

2. **Is `STEP_TR_C_theta_d17_wrapper.R` actually used in production, or
   just an experiment?** If unused, retire. If used, the JSON layer
   `theta_d17_envelopes` should become canonical and TR_C/TR_D should
   merge into TR_B (one script emits everything, less plumbing).

3. **Shared sample metadata: TR_B currently reads samples internally; can
   it consume `_shared/sample_metadata.tsv` instead?** Same question as
   GHSL. Yes, this is the right move for cross-path consistency.

4. **The optional `--intervals` mode for k-means karyotype calling
   mirrors GHSL's C04b interval-class output. Should both paths use the
   same input format for `--intervals`?** Currently they each define
   their own. Harmonize.

Discuss before refactoring.
