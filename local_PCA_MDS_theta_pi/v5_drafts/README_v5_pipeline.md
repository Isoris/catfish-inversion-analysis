# v5 pipeline drafts — theta-pi local-PCA + MDS + CUSUM

Six scripts. Each owns one stage. Drafts only — rearrange / split / merge
before commit.

## Stages and outputs

```
TR_A (existing, unchanged)
   └─ theta_native.<CHR>.<SCALE>.tsv.gz          per-chrom long TSV (sample × window θπ)

TR_B v5_precompute  ← REPLACES v4's all-in-one TR_B
   ├─ 03_per_chrom/<CHR>/precomp.rds              full precomp (dt + sim_mat + mds + bg_q + per-sample PC1/PC2)
   ├─ 03_per_chrom/<CHR>/sim_mat_nn{0,20,…,320}.rds   NN-smoothed sim_mats
   ├─ window_dt.tsv.gz                             genome-wide per-window scalar table
   └─ precomp_summary.tsv                          per-chrom QC

TR_C detect_L1   (cohort-level wide envelopes from |Z| of θπ deviation)
   └─ 03_per_chrom/<CHR>/L1_envelopes.tsv

TR_D detect_L2   (cohort-level tight envelopes inside each L1)
   ├─ 03_per_chrom/<CHR>/L2_envelopes.tsv
   └─ 02_mds/candidate_regions.tsv.gz              genome-wide rollup

TR_E classify_carriers   (k-means on per-sample PC1 inside each L2 → bands)
   └─ 03_per_chrom/<CHR>/carrier_assignments.tsv

TR_F cusum_per_carrier   (CUSUM per band per candidate; lib_persample_cusum.R kernel)
   ├─ 03_per_chrom/<CHR>/cusum_per_sample.tsv.gz   one row per (sample × candidate × band)
   └─ 03_per_chrom/<CHR>/cusum_boundary_dist.tsv   per-side spread classes (tight/intermediate/ragged)

TR_G atlas_json   (combine everything → page-12 JSON)
   └─ 04_atlas_json/<CHR>/<CHR>_phase2_theta.json
```

## Why both L1/L2 AND CUSUM (not "or")

L1/L2 and CUSUM detect different objects.

**L1 + L2 are cohort-level.** They scan max_abs_z of θπ deviation across all
226 samples and identify regions where the cohort signal is elevated. L1 is
the wide net (long runs, lenient z); L2 is the refined call inside each L1.
The output is a list of candidate regions — *where* to look. Direct port of
local_PCA_MDS_z/04_detect_L1 and 06_detect_L2 — same seed-and-grow logic,
same morphology gates (flat_inv_score, spiky_inv_score, fragmentation_score,
sim_mat block compactness), same beta-adaptive thresholding. Only difference
is the inv_likeness composite formula: dosage z-blocks uses
`45% het_contrast + 30% trimodality + 25% band_discreteness` (PC1 trimodality
indicators that don't apply to θπ); v5 substitutes
`50% normalized max_abs_z + 30% sim_mat block compactness + 20% λ-ratio`.
Everything else (the morphology features, beta p-values, NN-smoothed
sim_mats, seed-and-grow detector) is verbatim from z-blocks.

**CUSUM is per-carrier — TWO flavors.** Inside each L2 candidate, samples
are partitioned by their PC1 loading into bands (LOW_DIV / MID_DIV / HIGH_DIV
via k-means; TR_E). Then TR_F runs:

  1. *Per-sample CUSUM*: one breakpoint per (sample × candidate × band) →
     the carrier-spread distribution (tight / intermediate / ragged).
     Manuscript value: "carrier 1 breaks at 18.94 Mb, carrier 17 outlier at
     17.20 Mb → 3' boundary is bimodal".
  2. *Band-mean CUSUM*: pool the band's samples into one mean θπ trace
     (averaging reduces noise by sqrt(n_band)), CUSUM that single trace →
     ONE consensus breakpoint per (band × candidate × side). Sharp.
     Manuscript value: "the HIGH_DIV band's 3' boundary is at 18.95 Mb".

Both flavors are emitted side-by-side in `cusum_boundary_dist.tsv` and in the
atlas JSON's `theta_pi_cusum.candidates[].bands[].boundary_{5,3}_prime` block:
per-sample stats (`n_carriers`, `median_bp`, `iqr_kb`, `spread_class`) +
consensus stats (`consensus_cp_bp`, `consensus_strength`).

## Index convention

All emitted window indices in JSON are 0-indexed (atlas is JS).

R-internal stays 1-indexed. The conversion happens at JSON-emit time:
- TR_C / TR_D write both `win_start` (1-indexed, R-native) and
  `win_start_idx0` / `win_end_idx0` (0-indexed) to the L1/L2 TSVs.
  Downstream scripts pick whichever they need.
- TR_G reads `win_*_idx0` for the JSON.
- The v4 TR_B has the off-by-one issue noted in the audit (envelope coords
  emitted as 1-indexed); this v5 set fixes it at the source.

## Run order

```bash
source 00_theta_config.sh
# 1. theta matrices (already done if you ran v4):
$RSCRIPT STEP_TR_A_compute_theta_matrices.R --chrom C_gar_LG28
# 2. precompute (per-chrom or all):
$RSCRIPT v5_drafts/STEP_TR_B_v5_precompute.R --chrom C_gar_LG28
# 3-6. detect → classify → cusum → emit:
$RSCRIPT v5_drafts/STEP_TR_C_detect_L1.R         --chrom C_gar_LG28
$RSCRIPT v5_drafts/STEP_TR_D_detect_L2.R         --chrom C_gar_LG28
$RSCRIPT v5_drafts/STEP_TR_E_classify_carriers.R --chrom C_gar_LG28
$RSCRIPT v5_drafts/STEP_TR_F_cusum_per_carrier.R --chrom C_gar_LG28 --lib v5_drafts/lib_persample_cusum.R
$RSCRIPT v5_drafts/STEP_TR_G_atlas_json.R        --chrom C_gar_LG28
```

For the full 28-chrom run, drop `--chrom`; each script iterates `CHROM_LIST`
from the config (or scans `03_per_chrom/` for chroms that have a precomp).

## Knobs to tune

| Script | Knob | Default | Effect |
|---|---|---|---|
| TR_B | `--pad` | 1 | local-PCA neighbourhood half-width |
| TR_B | `--sim-band-half` | 200 | banded sim_mat half-width when n_win > threshold |
| TR_B | `--sim-n-full-threshold` | 6000 | n_win threshold to switch full ↔ banded |
| TR_C | `--z-l1` | 1.5 | lenient cohort-|Z| for L1 |
| TR_C | `--min-l1-windows` | 10 | min run length |
| TR_C | `--merge-gap` | 5 | merge L1 fragments |
| TR_D | `--z-l2` | 2.5 | strict cohort-|Z| for L2 inside each L1 |
| TR_D | `--min-l2-windows` | 5 | min run length |
| TR_E | `--max-k` | 3 | max k for k-means on PC1 |
| TR_F | `--bands` | (all) | comma-list of bands to run CUSUM on |

## Memory notes

- TR_B reconstructs the full sim matrix transiently for cmdscale + NN
  smoothing baseline. At LG28 (n_win=16,500) that's ~1 GB; if you OOM,
  the next iteration should add a coarse-grid MDS option (bin windows
  to 100-wide bins, MDS on bin-level sim, interpolate per-window).
- Per-chrom `precomp.rds` includes per-sample PC1/PC2 columns — at
  226 samples × 16,500 windows × 8 bytes × 2 PCs ≈ 60 MB raw. Each
  RDS gets gzipped on save (saveRDS default), typically <20 MB on disk.

## What's NOT in v5

- No `inv_likeness` composite (z-blocks-specific; rebuild from θπ
  signals if you want a cohort-level "inv-like" score).
- No GHSL / SV cross-stamping (Phase 4 catalog work).
- No interactive re-CUSUM on user-defined sample subsets (CUSUM_SPEC §8.2
  defers this to a future on-demand recomputation; cluster-side default
  is the source of truth for v1).