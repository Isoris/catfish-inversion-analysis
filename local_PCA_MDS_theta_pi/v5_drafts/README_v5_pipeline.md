# v5 pipeline drafts — theta-pi local-PCA + MDS + L1/L2 + CUSUM

## Architecture

The boundary-detection half (TR_C / TR_C_plot / TR_D / TR_D_plot) is a
**verbatim port** of the z-blocks scripts:

| v5 file                    | Source (z-blocks)                                  |
|----------------------------|----------------------------------------------------|
| `STEP_TR_C_detect_L1.R`    | `local_PCA_MDS_z/04_detect_L1/04_detect_L1_localpca_zblocks.R` |
| `STEP_TR_C_plot_L1.R`      | `local_PCA_MDS_z/05_plot_L1/05_plot_L1_localpca_zblocks.R`     |
| `STEP_TR_D_detect_L2.R`    | `local_PCA_MDS_z/06_detect_L2/06_detect_L2_localpca_zblocks.R` |
| `STEP_TR_D_plot_L2.R`      | `local_PCA_MDS_z/07_plot_L2/07_plot_L2_localpca_zblocks.R`     |

The detection logic (diagonal boundary scan, grow validator, Ward-adaptive
thresholding, quadrant validator) is signal-agnostic — it operates on
`precomp.rds` + `sim_mat_nn{N}.rds`, so any precomp shaped like the z-blocks
one will work. TR_B v5 produces exactly that shape.

## Stages

```
TR_A (existing, unchanged)
   └─ theta_native.<CHR>.<SCALE>.tsv.gz                per-chrom long TSV

TR_B v5_precompute   (REPLACES v4's all-in-one TR_B)
   ├─ precomp/<CHR>.precomp.rds                        z-blocks-shaped:
   │     dt: chrom, window_idx, start_bp, end_bp, mid_bp,
   │         theta_pi_median, theta_z_direct,
   │         MDS1..MDSk, MDS1_z..MDSk_z, max_abs_z, max_z_axis,
   │         lambda_1..lambda_NPC, lambda_ratio,
   │         anchor_window_idx, PC_1_<sample>..PC_NPC_<sample>
   │     sim_mat / sim_band, mds_mat, bg_continuity_quantiles,
   │     chrom, n_windows, n_samples, npc, k_mds,
   │     sample_order, unflipped_windows
   ├─ precomp/sim_mats/<CHR>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
   ├─ window_dt.tsv.gz                                 genome-wide rollup
   └─ precomp_summary.tsv

TR_C detect_L1            (verbatim 04_detect_L1; uses sim_mat_nn80 by default)
   ├─ L1_detect/<CHR>.L1_envelopes.tsv
   ├─ L1_detect/<CHR>.L1_boundaries.tsv
   └─ L1_detect/<CHR>.L1_score_curve.tsv

TR_C_plot                 (verbatim 05_plot_L1)
   └─ L1_plots/<CHR>.L1_overlay.pdf

TR_D detect_L2            (verbatim 06_detect_L2; sim_mat_nn40 inside each L1)
   ├─ L2_detect/<CHR>.L2_envelopes.tsv
   ├─ L2_detect/<CHR>.L2_boundaries.tsv
   ├─ L2_detect/<CHR>.L2_segment_stats.tsv
   ├─ L2_detect/<CHR>.L2_quadrant_validator.tsv
   └─ L2_detect/<CHR>.L2_quadrant_audit.tsv

TR_D_plot                 (verbatim 07_plot_L2)
   └─ L2_plots/<CHR>.L2_overlay.pdf

TR_E classify_carriers    (k-means on per-sample mean PC1 inside each L2)
   └─ carriers/<CHR>.carrier_assignments.tsv

TR_F cusum_per_carrier    (per-sample + band-mean CUSUM; lib_persample_cusum.R)
   ├─ cusum/<CHR>.cusum_per_sample.tsv.gz              one row per (sample × candidate × band)
   └─ cusum/<CHR>.cusum_boundary_dist.tsv              per-side, per-band: distribution + consensus

TR_G atlas_json           (combine everything → page-12 JSON)
   └─ 04_atlas_json/<CHR>/<CHR>_phase2_theta.json
```

## Why both L1/L2 (cohort) AND CUSUM (per-carrier)

Different objects.

**L1 + L2 (z-blocks ports)** are cohort-level. They scan the sim_mat for
boundary peaks (where adjacent regions are unusually separated) using
the diagonal cross-block scan, validate via grow + quadrant tests, and
partition the chromosome into segments. The output is the candidate
list — *where* to look. Identical algorithm to z-blocks; only the input
sim_mat differs (here built from `|cor(pc1[, i], pc1[, j])|`).

**CUSUM is per-carrier — TWO flavors.** Inside each L2 candidate, samples
are partitioned by their PC1 loading into bands (`LOW_DIV` / `MID_DIV` /
`HIGH_DIV` via k-means; TR_E). Then TR_F runs:

  1. *Per-sample CUSUM*: one breakpoint per (sample × candidate × band).
     Output is the carrier-spread distribution (`tight` / `intermediate` /
     `ragged` per IQR threshold).
  2. *Band-mean CUSUM*: pool the band's samples into one mean θπ trace,
     CUSUM that single trace → ONE consensus breakpoint per (band × side).
     Sharper because pooling reduces noise by sqrt(n_band).

Both stay in `cusum_boundary_dist.tsv` and the atlas JSON's
`theta_pi_cusum.candidates[].bands[].boundary_{5,3}_prime` block.

## Index convention

- z-blocks 04/06 emit `win_start` / `win_end` / `boundary_w` as 1-indexed
  (R-native). TR_G converts to 0-indexed at JSON-emit time so the atlas
  (JS) reads positionally without off-by-one.
- TR_B v5 stores `anchor_window_idx` as 0-indexed in the precomp `$dt`.
- Everything inside R stays 1-indexed.

## NPC (number of PCs to keep)

Default `NPC = 4` in TR_B v5 (was 2 in v4). Storage cost: `NPC × 226
samples × 16,500 windows × 8 bytes ≈ 30 MB per PC raw`, ~10 MB on disk
each after gzip. Bump via `--npc 5`. The atlas JSON layer's
`pc_loadings_aligned` is a list of length NPC.

`K_MDS = 5` (matches z-blocks `SEED_MDS_AXES`). MDS1..MDS5 + MDS1_z..MDS5_z
emitted on `$dt`; `max_abs_z` is `max(|MDSk_z|)` over the K_MDS axes
(the z-blocks definition that 04 expects).

A separate `theta_z_direct` column carries the v4-style θπ-direct |Z|
(per-sample dev from window cohort median). It's kept as an alt track
for the atlas; not used by 04/06.

## Run order

```bash
source 00_theta_config.sh

# Theta matrices (already done if you ran v4):
$RSCRIPT STEP_TR_A_compute_theta_matrices.R --chrom C_gar_LG28

# Precompute (per-chrom or all):
$RSCRIPT v5_drafts/STEP_TR_B_v5_precompute.R --chrom C_gar_LG28

# Boundary detection (z-blocks verbatim ports — outdir is up to you):
$RSCRIPT v5_drafts/STEP_TR_C_detect_L1.R \
    --precomp_dir $OUTROOT/precomp --chr C_gar_LG28 \
    --outdir      $OUTROOT/L1_detect

$RSCRIPT v5_drafts/STEP_TR_D_detect_L2.R \
    --precomp_dir $OUTROOT/precomp --chr C_gar_LG28 \
    --L1_dir      $OUTROOT/L1_detect \
    --outdir      $OUTROOT/L2_detect

# Plots (optional):
$RSCRIPT v5_drafts/STEP_TR_C_plot_L1.R \
    --precomp_dir $OUTROOT/precomp --L1_dir $OUTROOT/L1_detect \
    --chr C_gar_LG28 --outdir $OUTROOT/L1_plots

$RSCRIPT v5_drafts/STEP_TR_D_plot_L2.R \
    --precomp_dir $OUTROOT/precomp \
    --L1_dir $OUTROOT/L1_detect --L2_dir $OUTROOT/L2_detect \
    --chr C_gar_LG28 --outdir $OUTROOT/L2_plots

# Carriers + CUSUM + atlas JSON:
$RSCRIPT v5_drafts/STEP_TR_E_classify_carriers.R --chr C_gar_LG28
$RSCRIPT v5_drafts/STEP_TR_F_cusum_per_carrier.R --chr C_gar_LG28 \
    --lib v5_drafts/lib_persample_cusum.R
$RSCRIPT v5_drafts/STEP_TR_G_atlas_json.R --chr C_gar_LG28
```

For the full 28-chrom run, omit `--chrom` / `--chr` from each script and
the iterators pick up `CHROM_LIST` from the config (or scan
`<OUTROOT>/precomp/` for chroms that have a precomp).

## What changed since the previous draft

- **Dropped the morphology / `inv_likeness` / `beta_adaptive` machinery
  from TR_B.** Z-blocks 04/06 don't read those fields; they were dead
  weight. The precomp is now lean — only what the verbatim 04/06
  actually consume, plus the per-sample PC loadings the atlas needs.
- **TR_C / TR_D are now the actual z-blocks 04/06 scripts**, copied
  verbatim. They do diagonal boundary scan + grow validator (04) and
  segment-internal scan + Ward-adaptive thresholding + quadrant
  validator (06). Same algorithms as dosage; different input sim_mat.
- **Plots added**: TR_C_plot and TR_D_plot are 05 and 07 verbatim. Each
  emits a multi-page PDF (whole-chrom + per-segment zooms) with L1/L2
  envelope overlays.
- **TR_B `max_abs_z` is now the z-blocks definition** (max over K_MDS
  axis robust z-scores), so 04/06 find it where they expect. The
  v4-style θπ-direct |Z| is preserved as `theta_z_direct`.
- **NPC generalized** (default 4). Per-sample loadings `PC_1_<s>` ..
  `PC_NPC_<s>` are stored on `$dt`.
- **Output layout matches z-blocks**: `precomp/<chr>.precomp.rds` +
  `precomp/sim_mats/<chr>.sim_mat_nn{N}.rds`. 04/06's auto-resolution
  (`--precomp_dir <dir>`) finds them automatically.

## What's NOT in v5

- No GHSL / SV cross-stamping (Phase 4 catalog work).
- No interactive re-CUSUM on user-defined sample subsets (CUSUM_SPEC §8.2
  defers this to a future on-demand recomputation).
- The CUSUM kernel `lib_persample_cusum.R` is the same one in
  `spec_cusum/`; nothing new there.