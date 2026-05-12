# `local_PCA_MDS_GHSL/` — GHSL phased-haplotype discovery path

Path 3 of `catfish-inversion-analysis`. Same conceptual flow as
`local_PCA_MDS_z/` (dosage) and `local_PCA_MDS_theta_pi/` (θπ), but the
upstream feature is **per-window phased-haplotype divergence** computed
from merged phased Clair3 SNPs, not dosage or θπ.

For the high-level run sheet, see **HOW_TO_RUN_LG28.txt**. For path
symmetry across z / θπ / GHSL, see
[`../local_PCA_MDS_z/README.md`](../local_PCA_MDS_z/README.md) §14.

**Layout:** flat. Step scripts and SLURM launchers live at the top
level, prefixed `STEP_GH_` (A→E). The shared chrom list and config live
one level up.

---

## 1. Pipeline at a glance

```
postprocess_results/<chr>/<sample>/all_variants_with_phase.tsv
   │ (per-sample, per-chrom Clair3 postprocess output)
   ▼
merged phased SNPs               (STEP_GH_prep — bash, ~few min/chrom)
ghsl_prep/<chr>.merged_phased_snps.tsv.gz
   │
   ▼
divergence matrices              (STEP_GH_A — heavy, ~1 hr/chrom × 28 array)
   │
   ▼
classifier (PASS-runs +          (STEP_GH_B — ~30 s/chrom; PRIMARY biological
karyotypes + interval k-means     candidates + interval-CUSUM decomposition)
+ per-interval CUSUM)
   │
   ▼
local PCA precompute             (STEP_GH_C — ~5–10 min/chrom; sim_mat for D17 +
+ sim_mat + secondary             secondary |Z|-threshold envelopes)
|Z| envelopes
   │
   ▼
D17 detect_L1 + detect_L2        (STEP_GH_D — ~10–60 s/chrom; PRIMARY
                                  architectural candidates)
   │
   ▼
page-3 atlas JSON                (STEP_GH_E — packs all of the above into
                                  one <chr>_phase2_ghsl.json)
```

Five stages, lettered A–E, plus a pre-stage (`STEP_GH_prep_*`) that
merges the per-sample Clair3 postprocess output into per-chromosome
input for STEP_GH_A. Each lettered stage has either its own SLURM
launcher (`LAUNCH_STEP_GH_A_compute.slurm`,
`LAUNCH_STEP_GH_B_classify.slurm`) or is bundled together for a single
per-chrom run via `LAUNCH_STEP_GH_CDE_enrichment.slurm` (which does C
→ D → E in one job). prep → A → B → CDE is the canonical run order;
the prep step is one-time per cohort.

## 2. Why GHSL has an extra stage compared to z and θπ

Two parallel candidate streams ship from this folder, and the atlas
overlays both on page 3:

| Stream | Source | Authority |
|---|---|---|
| `ghsl_envelopes` (PRIMARY biological) | STEP_GH_B PASS-runs | calibrated on real signal — denominator-confound aware |
| `ghsl_d17_envelopes` (PRIMARY architectural) | STEP_GH_D D17 boundary scan | geometry-based edge detection on the sim_mat |
| `ghsl_secondary_envelopes` | STEP_GH_C |Z|-threshold scan | 1D fallback, kept for cross-check |

Path 1 (dosage) and path 2 (θπ) emit only one primary candidate stream
each — they don't have a calibrated upstream classifier. GHSL has both
because phased-haplotype divergence is biologically interpretable
(low div in inversion homozygotes, high div in heterozygotes) and the
classifier exploits that directly.

## 3. The stages

### `STEP_GH_prep_merged_phased_snps.sh` — pre-stage (one-time per cohort)

Merges per-sample Clair3 postprocess TSVs
(`postprocess_results/<chr>/<sample>/all_variants_with_phase.tsv`) into
one per-chromosome file
(`ghsl_prep/<chr>.merged_phased_snps.tsv.gz`) consumed by STEP_GH_A.
Filters to `IS_SIMPLE_BIALLELIC=TRUE && IS_SNP=TRUE`; keeps phase
information (`is_phased`, `phase_gt`, `ps`, `phase_block_id`,
`phase_tier`) plus QC fields (`qual`, `gq`, `dp`).

```bash
# Single chromosome
bash STEP_GH_prep_merged_phased_snps.sh <postprocess_dir> <ghsl_prep_outdir> C_gar_LG28

# All 28 in parallel as a SLURM array
sbatch --array=1-28 STEP_GH_prep_merged_phased_snps.sh <postprocess_dir> <ghsl_prep_outdir>
```

Run once when new Clair3 data arrives; STEP_GH_A onward reads from the
merged `ghsl_prep/` outputs.

### `STEP_GH_A_compute_matrices.R` — heavy engine (~1 hr/chrom, run once)

Input: `<chr>.merged_phased_snps.tsv.gz` (~77M variants total) + window
grid from the dosage precomp.

1. Stage 1: builds `div_mat[226 × N_windows]` (per-sample phased-het
   fraction) and `het_mat` (all-het fraction) at the 5-kb base scale.
2. Stage 2: applies `frollmean(align = "center")` at the configurable
   scale ladder (default `10, 20, 30, 40, 50, 100` windows ≈
   50/100/150/200/250/500 kb). All scales share the same base matrix so
   the cost is dominated by Stage 1.
3. Stage 3: writes one RDS per chromosome —
   `<chr>.ghsl_matrices.rds` containing raw matrices, all rolling
   matrices, window coords, sample names, and the param block used.

No scoring, no classification. Just data prep.

### `STEP_GH_B_classify.R` — light classifier (~30 s/chrom, iterate)

Reads the matrices RDS. Four independent stages:

- **A. Window metrics.** Rank stability (Spearman ρ between adjacent
  rolling windows), bimodality, contrast, tightness, z-scored against
  the chromosome baseline. Per-window status: `PASS / WEAK / FAIL`.
- **B. Karyotype calling.** Per-sample runs of stable LOW (INV/INV
  candidate) or HIGH (INV/nonINV candidate) rolling ranks.
- **C. Interval classification.** When `--intervals triangle_intervals.tsv.gz`
  is supplied, runs k-means with silhouette selection across k ∈ 2..5
  on per-sample interval-mean divergence. At k=3: INV/INV / HET /
  INV_nonINV per sample per interval.
- **D. Interval decomposition (per-sample CUSUM).** Per-sample
  changepoint detection inside each interval, clustered by changepoint
  position to detect whether one interval contains 2+ overlapping
  inversion systems. Emits per-sample asymmetry, slope, profile_var.

This is GHSL's analogue of θπ's `STEP_TR_F_cusum_per_carrier.R` — same
spirit (per-sample changepoint per candidate region), bundled into the
classifier rather than split into its own stage.

Outputs (TSVs, in `<outdir>/`):

- `ghsl_window_track.tsv.gz` — per-window metrics
- `ghsl_karyotype_calls.tsv.gz` — per-sample stable runs
- `ghsl_interval_genotypes.tsv.gz` — interval-level classification
  (when `--intervals` provided)
- `ghsl_interval_decomp.tsv.gz` — sub-system decomposition
  (when `--intervals` provided and changepoints separate cleanly)
- `ghsl_summary.tsv` — one row per chromosome

Per-chromosome RDS shards (consumed by STEP_GH_C / STEP_GH_E):

- `annot/<chr>.ghsl_annot.rds` — thin per-window aggregates
  (one row per window). Columns include `ghsl_score`, `ghsl_status`,
  `rank_stability`, `div_contrast_z`, `div_bimodal`, plus coordinates.
- `annot/<chr>.ghsl_karyotypes.rds` — per-sample stable LOW/HIGH runs
  (one row per run).
- `per_sample/<chr>.ghsl_per_sample.rds` — **dense long-format panel**,
  one row per (sample × window). Carries per-scale `div_roll_<s>`,
  `rank_in_cohort_<s>`, `rank_band_<s>` at every saved scale, plus
  `in_stable_run`, `stable_run_call`, and window-level `ghsl_score` /
  `ghsl_status`.

### `STEP_GH_C_precompute.R` — local PCA + sim_mat (~5–10 min/chrom)

Same role as θπ's `STEP_TR_B_v5_precompute.R` and dosage's
`STEP_ZO_G_precompute.R`: builds the precomp + sim_mat that the L1/L2
boundary detector consumes.

Reads `<chr>.ghsl_matrices.rds` and produces
`<chr>.ghsl_localpca.rds`:

- per-window pc1 / pc2 loadings (raw + sign-aligned to anchor window)
- λ₁, λ₂, λ_ratio (1D-ness indicator)
- robust |Z| profile of per-sample population deviations
- dense `sim_mat[N_windows × N_windows]` from `|cor(pc1[i], pc1[j])|`
- MDS coords (cmdscale of `1 - sim_mat`, k = 2)
- secondary `|Z|`-threshold L2 / L1 envelopes (cross-check layer)

Local PCA runs on the **raw** `div_mat` at 5-kb base — cross-sample
averaging across 226 samples reduces covariance noise to well below
biological signal; smoothed input would create artificial
autocorrelation. Pass `--smoothing-scale s50` to fall back to rolling
input if real data proves too noisy. Heteroscedastic weighting by
`sqrt(n_phased_het)` downweights samples with sparse phased calls per
window.

### `STEP_GH_D_detect_L1L2.R` — D17 boundary detector (~10–60 s/chrom)

Wraps the validated D17 multipass detector around STEP_GH_C's sim_mat
and emits L1 / L2 envelopes and boundaries in the same TSV shape as the
dosage and θπ pipelines:

```
<chr>_ghsl_d17L1_envelopes.tsv
<chr>_ghsl_d17L1_boundaries.tsv
<chr>_ghsl_d17L1_boundary_score_curve.tsv
<chr>_ghsl_d17L2_envelopes.tsv
<chr>_ghsl_d17L2_boundaries.tsv
```

D17's core statistic is signal-agnostic — per-diagonal Z-normalized
sim_mat, median over a WxW upper-triangle cross-block. It works on any
sim_mat with self-similarity ~ 1 and pairwise similarity decreasing
with pattern dissimilarity. Adaptive thresholding (default mode)
self-calibrates per chromosome from observed `grow_max_z` quantiles.

### `STEP_GH_E_atlas_json.R` — page-3 atlas JSON exporter

Consolidates 4 source RDSes + 4 D17 TSVs into one
`<chr>_phase2_ghsl.json`. Layers emitted:

| Layer | Source |
|---|---|
| `tracks` | STEP_GH_B annot RDS aggregates |
| `ghsl_panel` | STEP_GH_B per_sample RDS |
| `ghsl_kstripes` | computed K=2..6 stripe assigns |
| `ghsl_karyotype_runs` | STEP_GH_B karyotypes RDS |
| `ghsl_local_pca` | STEP_GH_C localpca RDS |
| `ghsl_envelopes` (PRIMARY biological) | STEP_GH_B annot PASS-runs |
| `ghsl_secondary_envelopes` | STEP_GH_C `z_profile` threshold |
| `ghsl_d17_envelopes` (PRIMARY architectural) | STEP_GH_D D17 TSVs |

## 4. Run order

Pre-stage (one-time per cohort), then three SLURM launchers, in this order:

```bash
# Pre-stage: merge per-sample Clair3 postprocess output → per-chrom TSV
sbatch --array=1-28 STEP_GH_prep_merged_phased_snps.sh \
       <postprocess_dir> ${GHSL_PREP_DIR}

# Stage A: heavy compute, ~1 hr/chrom × 28 in parallel (~1 hr wall)
sbatch LAUNCH_STEP_GH_A_compute.slurm

# Stage B: light classifier, ~30 s/chrom (single-job loop)
sbatch LAUNCH_STEP_GH_B_classify.slurm

# Stages C + D + E: page-3 enrichment, runs C → D → E per chrom
sbatch LAUNCH_STEP_GH_CDE_enrichment.slurm
```

After the CDE job, page-3 JSONs land at
`${GHSL_DIR}/json_out/<chr>/<chr>_phase2_ghsl.json` and can be drag-dropped
into the atlas alongside `<chr>_phase2_theta.json` to populate page 3 and
page 12.

For an interactive single-chromosome run sheet (LG28), see
**HOW_TO_RUN_LG28.txt**.

## 5. Defaults

- Heavy-engine scale ladder: `10,20,30,40,50,100` (windows). Override
  with `GHSL_SCALES=...`.
- Classifier primary scale: `s50` (250 kb rolling).
- Karyotype quantile cutoffs: `karyo_lo=0.15`, `karyo_hi=0.70`.
- Karyotype min run: `10` consecutive windows.
- Interval k-means: `max_k=5` (classifier) / `max_k=6` (atlas export).
- D17 threshold mode: adaptive (self-calibrating per chromosome).

## 6. Known caveats

- Per-sample divergence denominator is the per-sample variant count,
  which correlates with karyotype. Rolling smoothing at 50–100 windows
  absorbs much of the per-sample denominator variance, and Stage C's
  interval k-means uses interval means rather than per-sample ratios
  alone — which reduces the denominator's influence on final calls.
  Empirical question to revisit: correlation of `n_sites_mat` row-means
  against final interval classifications.
- Quantile cutoffs `0.15 / 0.70` in Stage B bias toward low-frequency
  inversions. For high-frequency inversions Stage B under-calls INV/INV
  anchors, but Stage C (k-means, no fixed quantiles) is unbiased — so
  the interval-level layer rescues high-frequency cases.
- Secondary `|Z|`-threshold envelopes (STEP_GH_C) are a 1D scan and can
  miss edges; D17 (STEP_GH_D) is the architecturally sharper detector.
  Both ship; the atlas overlays both for cross-checking.

## 7. Naming history

This folder previously used legacy identifiers `STEP_C04*`, `snake3`,
`v6`, and column prefixes `snake3v6_*` / `ghsl_v6_*`. All renamed in
this folder to the canonical `STEP_GH_*` (A→E) scheme matching
`local_PCA_MDS_z/` (`STEP_ZO_*`) and `local_PCA_MDS_theta_pi/`
(`STEP_TR_*`). Output column / RDS field names dropped the `_v6` /
`snake3v6_` prefixes accordingly.

`STEP_GH_E_atlas_json.R` retains a back-compat fallback for the
`ghsl_v6_status` / `ghsl_v6_score` column names so it still reads
already-generated RDS shards on the cluster — the column-name probe
prefers `ghsl_status` / `ghsl_score` first, then falls back. Any new
runs through STEP_GH_B emit the new names directly.
