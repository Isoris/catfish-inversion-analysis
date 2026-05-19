# `local_PCA_MDS_theta_pi/` — theta-pi local-PCA + MDS pipeline (path 2)

Path 2 of catfish-inversion-analysis. Sibling to `local_PCA_MDS_z/` (path 1,
dosage-based) and `local_PCA_MDS_GHSL/` (path 3, phased haplotypes).

## Architecture

The boundary-detection half (TR_D detect_L1, TR_E plot_L1, TR_F detect_L2,
TR_G plot_L2) is a **verbatim port** of the dosage path's z-blocks scripts:

| theta-pi step           | Source (dosage z-blocks)         |
|-------------------------|----------------------------------|
| `STEP_TR_D_detect_L1.R` | `STEP_ZO_H_detect_L1.R`          |
| `STEP_TR_E_plot_L1.R`   | `STEP_ZO_I_plot_L1.R`            |
| `STEP_TR_F_detect_L2.R` | `STEP_ZO_J_detect_L2.R`          |
| `STEP_TR_G_plot_L2.R`   | `STEP_ZO_K_plot_L2.R`            |

The detection logic (diagonal boundary scan, grow validator, Ward-adaptive
thresholding, quadrant validator) is signal-agnostic — it operates on
`precomp.rds` + `sim_mat_nn{N}.rds`, so any precomp shaped like the dosage
one will work. `STEP_TR_C_mds_compute.R` produces exactly that shape.

## Letter map (flat layout, B–J)

| Letter | Script | Purpose |
|--------|--------|---------|
| A | `STEP_TR_A_compute_theta_matrices.R` | reads pestPG → emits `theta_native.<chr>.<scale>.tsv.gz` (per-site `tP/nSites` + raw sum + n_sites) |
| **B** | `STEP_TR_B_local_pca_compute.R` | per-window heteroscedastic local PCA + anchor-flip → `<chr>.window_pca.rds` |
| **C** | `STEP_TR_C_mds_compute.R` | builds sim_mat → MDS → features → `<chr>.precomp.rds` (+ `sim_mats/` in full mode) |
| D | `STEP_TR_D_detect_L1.R` | chrom-wide boundary scan @ nn80 → L1 envelopes |
| E | `STEP_TR_E_plot_L1.R` | multi-page L1 overlay PDF |
| F | `STEP_TR_F_detect_L2.R` | per-L1-segment boundary scan @ nn40 → L2 envelopes |
| G | `STEP_TR_G_plot_L2.R` | multi-page L1+L2 overlay PDF |
| H | `STEP_TR_H_classify_carriers.R` | k-means on per-sample mean PC1 inside each L2 → carrier bands |
| I | `STEP_TR_I_cusum_per_carrier.R` | per-sample + band-mean CUSUM per (L2 envelope × band) |
| J | `STEP_TR_J_atlas_json.R` | combine into atlas JSON for the page-12 viewer |

Plus `lib_persample_cusum.R` — helper sourced by step I.

The retired predecessors (`STEP_TR_B_classify_theta.R` and friends) live in
`_legacy/` for reference. They are not part of the current chain.

## Dual-scale design (coarse + dense)

The MDS step's full-N×N reconstruction is the memory bottleneck. A
win10000.step2000 grid on LG28 gives ~16,500 windows → 2.18 GB sim_mat →
~10–15 GB peak during cmdscale. Splitting scales avoids this:

```
COARSE LAYER  (PESTPG_SCALE = win50000.step10000, ~5,000 windows on LG28)
   A  → B → C  --mode full   → precomp_coarse/<chr>.precomp.rds + sim_mats/
                                ↓
                                D detect_L1 → E plot_L1
                                F detect_L2 → G plot_L2
                                ↓
                                coarse L1/L2 envelopes (boundary resolution ~10–50 kb)

DENSE LAYER   (PESTPG_SCALE = win10000.step2000, ~16,500 windows on LG28)
   A  → B → C  --mode local  → precomp_dense/<chr>.precomp.rds (no sim_mats)
                                ↓
                                H classify_carriers (uses dense per-sample PC1)
                                I cusum_per_carrier (refines coarse boundaries)
                                ↓
                                refined boundaries + carrier bands

J atlas_json bundles coarse envelopes + dense refinements + carriers into the page-12 JSON.
```

`--mode full` does the full sim_mat reconstruction + cmdscale MDS + writes
the NN-smoothed sim_mats that L1/L2 detection consumes. Use at the coarse
scale where it's tractable.

`--mode local` skips full reconstruction and cmdscale entirely — just
banded sim_mat + per-window PC scores in the precomp. Use at the dense
scale for the per-sample / CUSUM / refinement workflow. **You cannot run
TR_D–G on a `--mode local` precomp** by design — there are no NN-smoothed
sim_mats to detect on.

**Dual-scale orchestration is wired** as of May 2026:

- `LAUNCH_TR_theta_pi.slurm` runs the full A→J chain at both scales per
  array task: coarse pass → dense pass → merge pass (H/I/J).
- `STEP_TR_H_classify_carriers.R` and `STEP_TR_I_cusum_per_carrier.R`
  default their `--precomp_dir` to `${OUTROOT}/precomp_dense` if it
  exists (falls back to `${OUTROOT}/precomp` for legacy single-scale).
- `STEP_TR_J_atlas_json.R` takes `--coarse_precomp_dir` and
  `--dense_precomp_dir` separately, builds a `theta_pi_grid_map` block
  with `coarse_to_dense[]` / `dense_to_coarse[]` window-index lookups
  (atlas uses these for cursor sync between heatmap clicks and dense
  per-sample line plots), and emits the per-window theta matrix at the
  DENSE scale while keeping the local_pca block (sim_mat thumbnails,
  MDS coords) at the COARSE scale.

Backward compat: TR_J still accepts the legacy `--precomp_dir` flag
(aliased to coarse) and falls back to single-scale mode if no dense
precomp exists.

## Knobs that matter

### `STEP_TR_C_mds_compute.R --mode full`
- `--kmds 5` (default) — number of MDS axes computed and stored. Detection
  only ever uses the first 5 (`SEED_MDS_AXES`); the old default of 20 wasted
  cmdscale work.
- `--sim-band-half 200` (default, in WINDOWS) — half-width of the banded
  sim_mat storage. Sets the maximum spatial scale at which window pairs
  retain real correlation; pairs farther apart get the band's median value
  during full reconstruction. At win50000.step10000 that's ±2 Mb of real
  correlation; at win10000.step2000 it'd be ±400 kb. Bump up if you suspect
  inversions wider than the band radius.
- `--sim-n-full-threshold 6000` — n_win above this stores banded, below
  stores full. Reconstruction to N×N is full either way in `--mode full`.

### `STEP_TR_C_mds_compute.R --mode local`
- `--sim-band-half` — same meaning, but storage is the only effect.
  No reconstruction is done.
- `--kmds`, `--sim-n-full-threshold` — ignored (no MDS performed).

### `STEP_TR_B_local_pca_compute.R`
- `--npc 4` (default; or env `LOCAL_PCA_NPC`) — number of PCs to keep.
  PC1 drives sim_mat; PC2 is a tiebreaker; PC3/PC4 carried through for
  the atlas's per-sample loading scatters.
- `--pad 1` (default) — local-PCA neighbourhood half-width.

## Outputs (under `$OUTROOT`)

| Path | Producer | Notes |
|------|----------|-------|
| `01_local_pca/theta_native.<chr>.<scale>.tsv.gz` | A | per-site theta_pi, raw tP_sum, n_sites |
| `01_local_pca/<chr>.window_pca.rds` | B | per-window PC scores + lambda + theta_z_direct (one per scale) |
| `precomp/<chr>.precomp.rds` | C `--mode full` (default subdir) | COARSE feature bundle + sim_mats |
| `precomp/sim_mats/<chr>.sim_mat_nn{0..320}.rds` | C `--mode full` | NN-smoothed sim_mats |
| `precomp_dense/<chr>.precomp.rds` | C `--mode local --out_subdir precomp_dense` | DENSE per-window features only |
| `L1_detect/<chr>.L1_*.tsv` | D | envelopes + boundaries + score curve (coarse scale) |
| `L1_plots/<chr>.L1_overlay.pdf` | E | |
| `L2_detect/<chr>.L2_*.tsv` | F | per-L1-segment partition (coarse scale) |
| `L2_plots/<chr>.L2_overlay.pdf` | G | |
| `carriers/<chr>.carrier_assignments.tsv` | H | k-means carrier band labels (uses dense PC1 inside coarse L2) |
| `cusum/<chr>.cusum_per_sample.tsv.gz` + `cusum/<chr>.cusum_boundary_dist.tsv` | I | per-sample CUSUM (dense theta inside coarse L2) |
| `04_atlas_json/<chr>/<chr>_phase2_theta.json` | J | merged dual-scale atlas bundle (coarse sim_mat + dense theta + grid_map) |

See [HOW_TO_RUN_LG28.txt](HOW_TO_RUN_LG28.txt) for canonical commands and
[README_theta_pi_scaling.md](README_theta_pi_scaling.md) for the per-site
vs. raw-sum θπ scaling caveat.

## Why path 2 is simpler than path 1

Path 1 (dosage) needs MDS because PCA eigenvectors are sign-ambiguous —
you can't directly compare PC1 of window i to PC1 of window j (they might
just be flipped). Path 1's MDS uses lostruct's sign-invariant angular
distance to recover meaningful similarity.

Path 2 (theta-pi) uses sign-stable per-window scalars (θπ is a magnitude,
no eigenvector flip ambiguity at the scalar level). The local PCA's PC1
still has a sign ambiguity per window, but the anchor-flip step in TR_B
resolves it: every window's PC1 is sign-aligned to the anchor window's
PC1 via correlation. So sim_mat is simply `|cor(PC1[,i], PC1[,j])|` —
no lostruct distance needed.
