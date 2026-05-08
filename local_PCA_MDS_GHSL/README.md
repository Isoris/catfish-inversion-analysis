# phase_2_discovery/2e_ghsl — GHSL haplotype contrast (Layer C), v6

*Terminology:* what older docs called "Snake 3" is this folder's GHSL
haplotype-contrast layer, which maps onto Layer C of the 4-layer
independence framework described in `inversion_modules/README.md`.
Filenames (e.g. `STEP_C04_snake3_ghsl_v6.R`) and output columns
(`snake3v6_*`) preserve the legacy identifier and are tracked for a
coordinated rename in `../2c_precomp/RENAMING.md` — not touched here.

**Status:** v6 installed 2026-04-18. v5 archived to
`_archive_superseded/2e_ghsl_v5/` (original README kept as
`README_v5.md` in that dir).

## Why v6 replaces v5

Two problems with v5:

1. **Per-window signal was too noisy.** The raw within-sample haplotype
   divergence at 5-kb window resolution bounced around enough that
   rank-stability scoring saturated on founder noise. PASS rates were
   dominated by chromosome-wide LD structure, not inversion signal.
2. **Iteration was impractical.** Every scoring-threshold tweak
   re-loaded 77M variants and recomputed the divergence matrix from
   scratch — about an hour per chromosome. Tuning scoring cutoffs
   required re-running the heavy compute, so iteration was effectively
   impossible.

v6 fixes both with a heavy/light split plus rolling-window smoothing.

## Two scripts, clean separation

### `STEP_C04_snake3_ghsl_v6.R` — heavy engine (~1 hr/chrom, run once)

- Stage 1 verbatim from v5: computes `div_mat[226 × N_windows]`
  (phased-het fraction) and `het_mat` (all-het fraction) from merged
  phased Clair3 SNPs.
- Stage 2 **new**: applies `frollmean(align = "center")` at
  configurable scales (**chat 14 default `10, 20, 30, 40, 50, 100`
  windows ≈ 50/100/150/200/250/500 kb** — a finer ladder than chat-13's
  `20, 50, 100`, chosen to resolve short recombinant tracts and
  sub-inversion blocks while keeping s100 as a chromosome-overview
  scale). All scales are computed from the same base matrix, so the
  extra scales add only a few minutes per chrom.
- Stage 3: saves everything as one RDS per chromosome —
  `<chr>.ghsl_v6_matrices.rds` containing raw matrices, all rolling
  matrices, window coords, sample names, and the param block used.
- No scoring, no classification, no PASS/FAIL. Just data prep.

### `STEP_C04b_snake3_ghsl_classify.R` — light classifier (~30 s, iterate)

Reads the v6 matrices RDS. Four independent stages:

- **A. Rolling metrics.** rank stability, bimodality, contrast,
  tightness, z-scored against chromosome baseline, PASS/WEAK/FAIL
  status per window.
- **B. Karyotype calling.** per-sample runs of stable LOW (INV/INV
  candidate) or HIGH (INV/nonINV candidate) rolling ranks. Same logic
  as v5 but on smoothed input, which makes the runs much cleaner.
- **C. Interval classification (new).** Given triangle intervals from
  `--intervals triangle_intervals.tsv.gz`, takes per-sample
  interval-mean divergence and runs k-means with silhouette selection
  across k in 2..5. At k=3 this produces INV/INV / HET / INV_nonINV
  per sample per interval. At k=2 it's LOW_DIV / HIGH_DIV (honest
  label — can't distinguish homozygote-INV from heterozygote with two
  clusters).
- **D. Interval decomposition (new).** Per-sample CUSUM changepoint
  inside each interval, clustered by changepoint position to detect
  whether one interval contains multiple overlapping inversion
  systems. Emits per-sample asymmetry, slope, profile_var.

## Usage

```bash
# Heavy, run once per chromosome:
Rscript STEP_C04_snake3_ghsl_v6.R \
  precomp_dir ghsl_prep_dir ghsl_v6_out \
  --chrom C_gar_LG01 \
  --scales 20,50,100

# Light, iterate freely:
Rscript STEP_C04b_snake3_ghsl_classify.R \
  ghsl_v6_out classify_out \
  --chrom C_gar_LG01 \
  --scale 50 \
  --intervals triangle_intervals.tsv.gz
```

Defaults: `scale=50` (250 kb rolling), `karyo_lo=0.15`, `karyo_hi=0.70`,
`karyo_min_run=10`, `rank_window=5`, `max_k=5`. Heavy-engine scale ladder
default (chat 14): `10,20,30,40,50,100`.

## Output files from the classifier

**Genome-wide TSVs (classifier sink, also used by figure scripts):**

- `snake3v6_window_track.tsv.gz` — per-window metrics
- `snake3v6_karyotype_calls.tsv.gz` — per-sample stable runs
- `snake3v6_interval_genotypes.tsv.gz` — interval-level classification
  (when `--intervals` provided)
- `snake3v6_interval_decomp.tsv.gz` — sub-system decomposition
  (when `--intervals` provided and changepoints separate cleanly)
- `snake3v6_summary.tsv` — one row per chromosome

**Per-chromosome RDS shards (chat 14 addition — primary downstream interface):**

- `annot/<chr>.ghsl_v6.annot.rds` — thin per-window aggregates
  (one row per window, no sample dimension). Columns include
  `ghsl_v6_score`, `ghsl_v6_status`, `rank_stability`, `div_contrast_z`,
  `div_bimodal`, plus coordinates. Consumed by `run_all.R` 2d block
  scoring.
- `annot/<chr>.ghsl_v6.karyotypes.rds` — per-sample stable LOW/HIGH
  runs (one row per run). Consumed by `lib_ghsl_confirmation.R` for
  Tier-3 SPLIT detection.
- `per_sample/<chr>.ghsl_v6.per_sample.rds` — **dense long-format panel**,
  one row per (sample × window). Carries per-scale `div_roll_<s>`,
  `rank_in_cohort_<s>`, `rank_band_<s>` at every scale present in the
  rolling matrix (default s10/s20/s30/s40/s50/s100), plus
  `in_stable_run`, `stable_run_call`, and window-level `ghsl_v6_score`/
  `ghsl_v6_status`. Has a `ghsl_panel_meta` attribute with sample
  order, window info, scales available, thresholds, and timestamp.
  This is the primary artifact for on-demand queries — see
  `utils/lib_ghsl_panel.R` and `reg$compute$ghsl_*`.

## Page-3 atlas enrichment (C04c / C04d / export_v3)

The C04 + C04b pipeline produces RDSes that 2d candidate detection
consumes. To populate the Inversion Atlas page 3 ("local PCA GHSL")
and its sibling diagnostic panels, three additional scripts run
*after* C04b:

| Script | Purpose | Output |
|---|---|---|
| `STEP_C04c_ghsl_local_pca.R` | Local PCA on the raw `div_mat` (5kb windows). Per-window pc1/pc2 loadings, lambda ratio, |Z| profile, sim_mat, MDS coords, sign-aligned loadings, |Z|-threshold L2/L1 envelopes. | `<chr>.ghsl_v6_localpca.rds` |
| `STEP_C04d_ghsl_d17_wrapper.R` | Wraps the validated D17 multipass detector around C04c's sim_mat. D17 is the same boundary-scan detector used on the dosage stream — diagonal cross-block test, adaptive thresholding, growing-W validator, Ward-D2 stratification. | `<chr>_ghsl_d17{L1,L2}_{envelopes,boundaries}.tsv` |
| `export_ghsl_to_json_v3.R` | Consolidates 4 source RDSes + 4 D17 TSVs into the page-3 atlas JSON. 8 layers: tracks, ghsl_panel, ghsl_kstripes, ghsl_karyotype_runs, ghsl_local_pca, ghsl_envelopes (PRIMARY = PASS-runs), ghsl_secondary_envelopes (|Z|), ghsl_d17_envelopes. | `<chr>_phase2_ghsl.json` |

Architecture: C04c uses RAW `div_mat` for local PCA (not rolling) —
cross-sample averaging across 226 samples reduces covariance noise to
~1/sqrt(226·pad) ≈ 0.04, well below biological signal; smoothed input
would create artificial autocorrelation between windows. A
`--smoothing-scale s50` fallback is available if real LG28 data
proves too noisy. Heteroscedastic weighting by `sqrt(n_phased_het)`
downweights samples with sparse phased calls per window.

D17 (C04d) is the same detector used on the dosage stream — its
core statistic is signal-agnostic (per-diagonal Z-normalized sim_mat,
median over WxW upper-triangle cross-block) and works on any sim_mat
with self-similarity ~ 1 and pairwise similarity decreasing with
pattern dissimilarity. Adaptive thresholding (default mode)
self-calibrates per chromosome.

### Run order

Three SLURM launchers, in this order:

```bash
# Stage 1: heavy compute, ~1 hr/chrom × 28 in parallel (~1 hr wall)
sbatch LAUNCH_STEP_C04_ghsl_v6_compute.slurm

# Stage 2: light classifier, ~30 s/chrom (single-job loop)
sbatch LAUNCH_STEP_C04b_ghsl_v6_classify.slurm

# Stage 3: page-3 enrichment, runs C04c → C04d → export per chrom
sbatch LAUNCH_STEP_C04cd_ghsl_enrichment.slurm
```

After Stage 3, page-3 JSONs land at `${GHSL_DIR}/json_out/<chr>/<chr>_phase2_ghsl.json` and can be drag-dropped into the atlas alongside `<chr>_phase2_theta.json` to populate page 3 and page 12.

## Wiring status

**Wired chat 14 (2026-04-18).** All four documented consumers read
v6 paths and columns:

- `phase_7_karyotype_groups/proposal/lib_ghsl_confirmation.R`
  reads v6 karyotype and annot RDS (keeps v5 fallback for transitional
  partial re-runs). Adds panel-based helper
  `ghsl_per_sample_panel_in_interval()` as an additive alternative to
  the run-overlap Tier-3 logic (run-overlap stays authoritative for
  SPLIT detection).
- `phase_2_discovery/2d_candidate_detection/run_all.R` reads v6 annot
  RDS with v6-native column names; v5 fallback retained for mixed-run
  cases.
- `phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
  reads `iv$ghsl_v6_score_max`.
- `phase_2_discovery/2d_candidate_detection/STEP_D05_ghsl_stability.R`
  (unused stub) — header updated to reference v6.

**On-demand query library wired at same time:**
`utils/lib_ghsl_panel.R` provides `load_ghsl_panel`, `ghsl_panel_range`,
`ghsl_panel_aggregate`, `ghsl_panel_subblock_scan` (for complex candidates
with soft boundaries giving N sub-blocks), `ghsl_panel_wide_matrix` (for
plotting), plus a handful of plot helpers. Registry integration via
`reg$compute$ghsl_at_candidate`, `ghsl_at_interval`, `ghsl_at_subblocks`,
`ghsl_at_block_subregions`, `ghsl_wide_matrix` — same pattern as the
existing `pairwise_stat` / `boundary_fst_profile` for FST/dxy live compute.

## Figures

The v5 figures script was archived because it reads v5's output layout.
v6 classifier outputs the same spirit of per-window metrics in a TSV,
so either a small update to the v5 figures script or a fresh v6
figures script is needed. Currently not included.

## What was kept from v5

- Stage 1 divergence compute (per-sample phased-het / total-variants
  per window). The denominator concern noted in the v5 audit
  (denominator = per-sample variant count, which correlates with
  karyotype) is still present structurally in v6. However, rolling
  smoothing at 50-100 windows absorbs a lot of the per-sample
  denominator variance, and Part C's interval classification uses
  k-means on interval means rather than per-sample ratios alone, which
  reduces the denominator's influence on final calls. Whether this is
  sufficient is an empirical question for the HPC calibration run —
  check correlation of `n_sites_mat` row-means against final interval
  classifications.
- Same karyotype-quantile cutoffs (0.15/0.70) in Part B stable-run
  calling. These bias toward low-frequency inversions; for high-
  frequency inversions Part B under-calls INV/INV anchors but Part C
  (k-means, no fixed quantiles) is unbiased.
