# HANDOFF — `2e_ghsl_discovery/` (path 3, awaiting refactor)

> **Status: working code, not yet refactored to the consolidated 01–09 skeleton.**
> This handoff documents what's here, how each script functions today, and
> how it would map into the canonical skeleton if/when refactored. Read this
> alongside `README.md` (operational quickstart) and `../../HANDOFF.md`
> (path 1 / canonical skeleton).
>
> Codebase: `inversion-popgen-toolkit` v8.5. The five active scripts in this
> folder ARE working on LANTA — don't break them. The refactor decision
> (do it or not) is a future-session call; see §5 below.

---

## 1. What this pipeline does

GHSL = "Ghost Sequence-based Linkage" — uses **per-sample phased haplotype
divergence** within sliding windows to detect inversions and other long-range
structural rearrangements.

The core observation: inside a heterozygous inversion carrier, the two
phased haplotypes diverge dramatically (because they're effectively two
"different chromosomes" in the inverted region — no recombination crosses
the breakpoint). Outside inversions, phased haplotypes converge.

Per-sample within-haplotype divergence (measured by
`ghsl_div = n_phased_het / n_total`) is therefore a direct, sign-stable
indicator of inversion heterozygosity. No PCA sign-flip ambiguity, no MDS
projection; raw divergence values are interpretable on their own.

The pipeline computes:
1. Per-sample × per-window divergence matrix from Clair3-phased VCFs (heavy)
2. Per-window classification & per-sample karyotype calls (light, iterable)
3. Window-local PCA on the divergence matrix (mirrors path 1 architecturally)
4. D17 boundary detection on the local-PCA sim_mat (canonical detector)
5. Page-3 atlas JSON aggregating all four output streams

---

## 2. Current file inventory

```
2e_ghsl_discovery/
├── STEP_C04_snake3_ghsl_v6.R              heavy compute  (~1 hr/chr)
├── STEP_C04b_snake3_ghsl_classify.R       light classifier (~30 s/chr, iterable)
├── STEP_C04c_ghsl_local_pca.R             window-local PCA on div_mat (~5-10 min/chr)
├── STEP_C04d_ghsl_d17_wrapper.R           wraps D17 detector on local-PCA sim_mat
├── export_ghsl_to_json_v3.R               page-3 atlas JSON consolidator
├── LAUNCH_STEP_C04_ghsl_v6_compute.slurm
├── LAUNCH_STEP_C04b_ghsl_v6_classify.slurm
├── LAUNCH_STEP_C04cd_ghsl_enrichment.slurm
└── README.md                              operational quickstart (don't replace)
```

---

## 3. What each script does, in detail

### 3.1 `STEP_C04_snake3_ghsl_v6.R` — heavy compute (no local PCA)

**Per chromosome, runs once.** Reads Clair3-phased VCF (~77M variants
genome-wide), computes per-sample within-haplotype divergence at the
**5-kb base scale**, applies rolling means at multiple configurable scales
(default 10/20/30/40/50/100 windows = 50 kb / 100 kb / ... / 500 kb), saves
everything as RDS for the light classifier to iterate against in seconds.

**This is NOT local PCA.** This is the *feature-build step* — analogous
to path 1's `01a_beagle_to_dosage.py`. It produces the divergence matrix
that everything downstream consumes.

**Inputs:**
- `<precomp_dir>` — path 1 precomp directory (used for chrom window grid)
- `<ghsl_prep_dir>` — Clair3-phased per-sample VCFs (genotype-quality filtered)
- `<outdir>` — output directory

**Outputs** (one RDS per chromosome):
- `<chr>.ghsl_v6_matrices.rds`
  - `$div_mat`         — raw GHSL divergence [N_samples × N_windows]
  - `$het_mat`         — raw het rate [N_samples × N_windows]
  - `$n_sites_mat`     — total variants per sample × window
  - `$n_phased_het_mat`— phased het count per sample × window
  - `$rolling`          — list of rolling div_mats by scale (s10/s20/...s100)
  - `$rolling_het`      — list of rolling het_mats
  - `$window_info`     — data.table with start_bp/end_bp per window
  - `$sample_names, $chrom, $params`

### 3.2 `STEP_C04b_snake3_ghsl_classify.R` — GHSL-specific classifier

**Per chromosome, runs once after C04.** ~30 s walltime per chrom, designed
for fast iteration.

Computes:
- **Per-window metrics** (rolling-smoothed): rank stability between adjacent
  windows, mean, bimodality, contrast, tightness, PASS/WEAK/FAIL z-status
- **Per-sample karyotype calls**: runs of stable LOW (INV/INV) or HIGH
  (INV/nonINV) divergence, length-filtered
- **Per-interval classification** (when triangle/stair intervals are
  supplied): k-means on per-sample interval-mean divergence, k=2..5 with
  silhouette selection → INV/INV / HET / nonINV/nonINV
- **Sub-system decomposition**: changepoint detection on per-sample profile
  inside an interval, groups samples by changepoint location

**This step is GHSL-specific and has no analogue in path 1.** It produces
the **PRIMARY** GHSL candidate set (PASS-runs in the annot.rds) that is
biologically calibrated on real signal.

**Outputs:**
- `snake3v6_window_track.tsv.gz`        — per-window metrics
- `snake3v6_karyotype_calls.tsv.gz`     — per-sample stable runs
- `snake3v6_interval_genotypes.tsv.gz`  — per-sample × per-interval class
- `snake3v6_interval_decomp.tsv.gz`     — sub-system decomposition
- `snake3v6_summary.tsv`                — chromosome-level summary
- (separately written by C04 layer) `<chr>.ghsl_v6.annot.rds`,
  `.per_sample.rds`, `.karyotypes.rds` — consumed by JSON exporter

### 3.3 `STEP_C04c_ghsl_local_pca.R` — local PCA on div_mat

**Per chromosome.** This is the actual local-PCA step for path 3. Runs
window-by-window PCA on the **raw** div_mat (not the rolling-smoothed one;
the script header explains why), with heteroscedastic per-sample weighting
by `sqrt(n_phased_het / median)`.

Key design choices documented in the script header:
- **Sign alignment**: per-window PC1/PC2 are sign-ambiguous from SVD. The
  highest-|Z| window is picked as anchor; every other window's PC1 is
  flipped if `cor(pc1[i], pc1[anchor]) < 0`. Both raw and aligned loadings
  are emitted.
- **Sim_mat is dense, no banding/quantization** at GHSL's 5-kb base.
  LG28 ≈ 4,300 windows → sim_mat ≈ 74 MB raw / ~18 MB gzip — within the
  120 MB-per-chrom JSON budget.
- **MDS via cmdscale(1 - sim_mat, k = 2)** — same as path 1's MDS step.
- **Secondary L1/L2 envelopes** from |Z|-threshold runs in the z_profile.
  These are SECONDARY; the PRIMARY GHSL candidates come from C04b.

**This script is conceptually `01b_local_pca_compute + 02a_mds_compute + 03_precompute` rolled into ONE.** It does sliding-window PCA, sim_mat,
MDS, and emits a precomp-equivalent RDS in a single pass.

**Outputs:**
- `<chr>.ghsl_v6_localpca.rds`
  - `$pc1_loadings_mat, $pc2_loadings_mat` (raw, sign-ambiguous)
  - `$pc1_loadings_aligned_mat, $pc2_loadings_aligned_mat`
  - `$lambda_1, $lambda_2, $lambda_ratio`
  - `$z_profile, $z_top10_mean`
  - `$sim_mat` (full symmetric float64 matrix)
  - `$mds1, $mds2`
  - `$anchor_window_idx`
  - `$secondary_l2_envelopes, $secondary_l1_envelopes` (data.tables)
  - `$window_info, $sample_names, $chrom, $params`

### 3.4 `STEP_C04d_ghsl_d17_wrapper.R` — calls D17 on the GHSL sim_mat

Adapter that runs path 1's L1+L2 boundary detector (D17, now
`04_detect_L1_localpca_zblocks.R` and `06_detect_L2_localpca_zblocks.R` in
the consolidated layout) on the GHSL local-PCA sim_mat.

D17's core statistic is **signal-agnostic** — per-diagonal Z-normalize
sim_mat, then median over an upper-triangle WxW cross-block. That works
on ANY sim_mat with self-similarity ~ 1 and pairwise similarity decreasing
with pattern dissimilarity, which describes GHSL's `|cor(pc1[i], pc1[j])|`
sim_mat exactly. Adaptive thresholding self-calibrates per chromosome
from observed `grow_max_z` quantiles.

**This wrapper currently invokes path 1's D17 scripts BY OLD NAME**:
- `STEP_D17_multipass_L1_only_v7.R` ← now `04_detect_L1_localpca_zblocks.R`
- `STEP_D17_multipass_L2_v8.R`      ← now `06_detect_L2_localpca_zblocks.R`

**This is broken in the consolidated layout** — the wrapper hardcodes the
old paths. Fix by either:
1. Updating the wrapper to call the new script names, OR
2. Symlinking the new scripts to the old names, OR
3. Refactoring per §5.

**Outputs** (currently use the old `_d17L*` suffix; will need renaming
when this wrapper is refactored):
- `<chr>_ghsl_d17L1_envelopes.tsv`
- `<chr>_ghsl_d17L1_boundaries.tsv`
- `<chr>_ghsl_d17L1_boundary_score_curve.tsv`
- `<chr>_ghsl_d17L2_envelopes.tsv`
- `<chr>_ghsl_d17L2_boundaries.tsv`

### 3.5 `export_ghsl_to_json_v3.R` — page-3 atlas JSON

Consolidates GHSL outputs from FOUR sources into one page-3 JSON:

1. **C04 heavy outputs** (`.annot.rds`, `.per_sample.rds`, `.karyotypes.rds`):
   - `tracks` — per-window aggregated annotations
   - `ghsl_panel` — per-sample × per-window heatmap data
   - `ghsl_kstripes` — computed K=2..6 stripe assignments
   - `ghsl_karyotype_runs` — stable LOW/HIGH runs
   - `ghsl_envelopes` (PRIMARY) — PASS-runs from C04b
2. **C04c local-PCA** (`.ghsl_v6_localpca.rds`):
   - `ghsl_local_pca` — pc1/pc2 loadings + λ ratio + sim_mat
   - `ghsl_secondary_envelopes` — |Z|-threshold runs (SECONDARY)
3. **C04d D17 wrapper TSVs**:
   - `ghsl_d17_envelopes` — boundary-detector edges (TERTIARY)
4. **Optional sample metadata** (`samples.tsv` with `ind, cga, ancestry`)

**Maps to canonical step 08b** but with a much richer per-page layer set
than path 1's atlas JSON. The atlas page-3 viewer renders all three
envelope sources (PRIMARY / SECONDARY / TERTIARY) in different colors.

---

## 4. Mapping to the canonical 01–09 skeleton

This pipeline does NOT decompose cleanly into the canonical skeleton as
currently written. The mapping is:

| Canonical step                  | GHSL equivalent                          |
|---------------------------------|------------------------------------------|
| 01a feature build               | `STEP_C04_snake3_ghsl_v6.R` (heavy)      |
| 01b local PCA compute           | part of `STEP_C04c_ghsl_local_pca.R`     |
| 01c local PCA merge             | (not needed — C04c is already per-chr)   |
| 02a MDS compute                 | part of `STEP_C04c_ghsl_local_pca.R`     |
| 02b MDS merge                   | (not needed)                             |
| 03 precompute                   | part of `STEP_C04c_ghsl_local_pca.R`     |
| 04 detect L1                    | `STEP_C04d` (wraps path 1's D17 L1)      |
| 05 plot L1                      | (not implemented here — handled by atlas)|
| 06 detect L2                    | `STEP_C04d` (wraps path 1's D17 L2)      |
| 07 plot L2                      | (not implemented here)                   |
| 08a sample metadata             | (consumes `_shared/sample_metadata.tsv` from path 1) |
| 08b atlas JSON                  | `export_ghsl_to_json_v3.R`               |
| **GHSL-specific extras**        | `STEP_C04b` (classifier — no path-1 analogue) |

The mismatch: **C04c bundles three canonical steps into one script**, and
**C04b is GHSL-specific** with no path-1 analogue.

---

## 5. Refactoring decision

**Option A — Refactor to canonical 01–09 skeleton.**

Pros:
- Three discovery paths look symmetric on disk
- D17 wrapper hack goes away (each path runs its own canonical 04+06)
- Easier to maintain in the long term

Cons:
- Splitting `STEP_C04c_ghsl_local_pca.R` into 01b+02a+03 equivalents is
  non-trivial — the script has cohesion (sign-alignment uses MDS-derived
  anchor; sim_mat builds on aligned PCs)
- `STEP_C04b` doesn't fit the skeleton at all (it's the GHSL-specific
  biological calibrator); needs to live somewhere as a "GHSL extras" step
- Cost: ~1 full session of careful surgery. Risk of regression on a
  pipeline that currently works end-to-end on LANTA.

**Option B — Lighter touch: keep current architecture, just rewire paths.**

What this looks like:
1. Move the folder to `local_PCA_GHSL/` (sibling of `local_PCA_z/`).
2. Rename the 5 R scripts to canonical names that reflect what they do:
   - `STEP_C04_snake3_ghsl_v6.R` → `01_compute_ghsl_matrices.R`
   - `STEP_C04b_snake3_ghsl_classify.R` → `02_classify_ghsl.R` (GHSL-specific)
   - `STEP_C04c_ghsl_local_pca.R` → `03_local_pca_ghsl.R` (compute+MDS+precomp)
   - `STEP_C04d_ghsl_d17_wrapper.R` → `04_detect_L1L2_ghsl.R` (wraps D17)
   - `export_ghsl_to_json_v3.R` → `05_export_atlas_json_localpca_GHSL.R`
3. Update C04d to call the new path-1 D17 script names
   (`04_detect_L1_localpca_zblocks.R`, `06_detect_L2_localpca_zblocks.R`).
4. Update output naming: `<chr>_ghsl_d17L1_envelopes.tsv` → `<chr>.GHSL_L1_envelopes.tsv`.
5. Have the JSON exporter consume `_shared/sample_metadata.tsv` instead of `--samples`.
6. Update launchers to write to `path_localpca_GHSL/` scratch tree.

Pros:
- Achievable in ~30 minutes
- No risk of breaking working pipeline
- Names + scratch tree look like the others; users see symmetry
- Still admits the architectural difference (C04b is GHSL-specific)

Cons:
- 5 scripts instead of canonical 11 — different shape from path 1
- A future "compare paths 1/2/3 step-by-step" exercise has friction

**Recommended: Option B**, with Option A deferred until there's a clear
reason to take on the surgery cost (e.g. you want to add a step to all
three paths uniformly and the irregularity becomes painful).

---

## 6. Required input data on LANTA

When refactored to scratch-tree convention, this pipeline reads:
- `<scratch>/04_clair3_phased_GHSL/<chr>.<phased VCF>` — Clair3 phased
  per-sample VCF, genotype-quality filtered. **Currently produced by a
  pipeline OUTSIDE this toolkit** (Clair3 + Whatshap or similar). Not
  reproduced by `local_PCA_z`.
- `<scratch>/path_localpca_zblocks/04_precomp/precomp/` — used for chrom
  window grid (only). Not used as path-1 features per se, just to align
  windows.

Outputs (when refactored) should land in:
- `<scratch>/path_localpca_GHSL/01_ghsl_matrices/<chr>.ghsl_v6_matrices.rds`
- `<scratch>/path_localpca_GHSL/02_ghsl_classify/`
- `<scratch>/path_localpca_GHSL/03_ghsl_localpca/<chr>.ghsl_v6_localpca.rds`
- `<scratch>/path_localpca_GHSL/04_GHSL_L1L2/<chr>.GHSL_L1_envelopes.tsv` etc.
- `<scratch>/path_localpca_GHSL/05_atlas_json/<chr>.atlas_GHSL.json`

---

## 7. What's NOT broken (don't fix)

The following work today on LANTA and should not be touched without reason:

- D17 cross-block detector design (signal-agnostic, self-calibrating
  thresholds — proven on dosage, demonstrated to work on GHSL sim_mat)
- C04b's biological calibration (the classifier was hand-tuned on real
  data; the magic numbers `KARYO_LO=0.15`, `KARYO_HI=0.70`,
  `KARYO_MIN_RUN=10` reflect tuning history)
- The four-source JSON consolidation in v3 — atlas page-3 expects exactly
  these layers
- Sign-alignment via highest-|Z| anchor window (page-3 cursor coherence
  depends on this)

---

## 8. Run order (current, working)

```bash
cd inversion_modules/phase_2_discovery/2e_ghsl/

# Stage 1: heavy compute (~1 hr/chrom × 28 in parallel = ~1 hr wall)
sbatch LAUNCH_STEP_C04_ghsl_v6_compute.slurm

# Stage 2: light classifier (~30 s/chrom)
sbatch --dependency=afterok:<JOB1> LAUNCH_STEP_C04b_ghsl_v6_classify.slurm

# Stage 3: page-3 enrichment (C04c → C04d → export per chrom)
sbatch --dependency=afterok:<JOB2> LAUNCH_STEP_C04cd_ghsl_enrichment.slurm
```

Output per chromosome: `${GHSL_DIR}/json_out/<chr>/<chr>_phase2_ghsl.json`.

---

## 9. Open questions for the refactor session

These need YOUR input before any rename happens:

1. **Should `STEP_C04b` be in the canonical sibling `local_PCA_GHSL/` folder
   at all?** It does inversion classification, not local PCA. Could argue
   it belongs in a separate `phase_4_classify/` or `inversion_modules/karyotype/`
   — leaving `local_PCA_GHSL/` as a pure "compute + detect" pipeline.

2. **Is the secondary-envelope detector in C04c (|Z|-threshold) still
   needed once C04d (D17) is the canonical detector?** The header says
   they ship as different layers. If both are useful, fine. If D17
   strictly dominates, simplify.

3. **JSON layer names: harmonize with path 1?** Path 1 emits
   `l1_envelopes` / `l2_envelopes`; GHSL emits `ghsl_envelopes` /
   `ghsl_secondary_envelopes` / `ghsl_d17_envelopes`. Keeping per-path
   prefixes lets the atlas page-3 layer-toggle cleanly, but breaks
   script-level symmetry.

4. **Does C04 still need the `<precomp_dir>` argument**, or can it use
   only `<scratch>/04_clair3_phased_GHSL/` and the chromosome FASTA index?
   The dependency on path 1's precomp for window grid creates an awkward
   inter-path coupling.

Discuss before refactoring.
