# Per-step notes — what each script does, eats, and emits

Reference for the consolidated `local_PCA_MDS_z/` pipeline (path 1).

> Naming uses the May-2026 single-letter scheme (`STEP_ZO_B..M`). Older
> references to `01a..08b`, numbered subfolders (`01_dosage_pca/`,
> `02_mds/`, …), or `STEP_D17_*` / `STEP_A03_*_stage1` etc. are all
> retired — those names will not match anything in the current tree.

---

## Step B — `STEP_ZO_B_beagle_to_dosage.py`

Per-chrom: reads beagle.gz genotype likelihoods → emits dosage matrix +
sites file. Pre-existing converter, byte-identical to what was in the old
layout.

**Output (in `<SCRATCH>/02_dosage_sites/`):**
- `<chr>.dosage.tsv.gz` — site × sample dosage matrix
- `<chr>.sites.tsv.gz`  — per-site allele/position metadata

---

## Step C — `STEP_ZO_C_local_pca_compute.R`

STAGE 1 (compute): one SLURM array task per chromosome runs sliding-window
PCA on the dosage matrix.

**Defaults:** `winsize=100`, `step=20`, `npc=4`, includes full eigenvalue
spectrum (`scree_full`, `pve_full`) for downstream scree plots.

**Output (in `<SCRATCH>/path_localpca_zblocks/01_local_pca/tmp/`):**
- `<chr>.window_pca_tmp.rds` — per-chrom PCA + scree, `window_id = NA`
- `<chr>.chr_meta.tsv.gz`    — per-window metadata
- `<chr>.pca_table.tsv.gz`   — per-window PC loadings
- `<chr>.summary.tsv`        — per-chrom summary

---

## Step D — `STEP_ZO_D_local_pca_merge.R`

STAGE 2 (merge): single short job. Reads all per-chr
`tmp/*.chr_meta.tsv.gz`, assigns globally unique sequential `window_id`
values, patches each per-chr RDS with the real window_id.

**Output (in `<SCRATCH>/path_localpca_zblocks/02_dense_registry/`):**
- `windows_master.tsv.gz`       — THE master window registry
- `windows_master_summary.tsv`  — per-chrom counts
- `<chr>.window_pca.rds`        — per-chrom PCA, finalized
- `<chr>.window_pca.tsv.gz`     — per-chrom PCA table

---

## Step E — `STEP_ZO_E_mds_compute.R`

STAGE 1 (compute): one SLURM array task per **focal** chromosome. Loads
ALL per-chr `.window_pca.rds` (yes, all 28 — needed for chunked background
sampling), computes lostruct distance matrix on focal + background windows,
runs `cmdscale`, extracts focal-only coords.

**Defaults:** `mds_mode=chunked_2x`, `npc=4`, `mds_dims=20`, `z_thresh=3.0`.

The chunked mode samples `2× n_focal` background windows from non-focal
chromosomes, **excluding high-inv_likeness** ones (so other chromosomes'
inversions don't leak into "background"). `chunked_3x` and `chunked_4x` use
3× and 4× respectively. `chromosome` mode uses no background. `global` uses
everything.

**Output (in `<SCRATCH>/path_localpca_zblocks/03_mds/tmp/`):**
- `<focal_chr>.mds_perchr.rds` — per-focal-chr MDS coords + metadata

This is the heaviest step in the pipeline (~12h walltime per chromosome).

---

## Step F — `STEP_ZO_F_mds_merge.R`

STAGE 2 (merge): single short job. Stitches per-focal-chr MDS RDS files
into the unified `inversion_localpca.mds.rds` with `$per_chr` structure
(the format every downstream consumes).

Runs candidate-region detection (z>3 contiguous windows, gap-merging).

**Output (in `<SCRATCH>/path_localpca_zblocks/03_mds/`):**
- `inversion_localpca.mds.rds`                     — main RDS
- `inversion_localpca.window_mds.tsv.gz`           — per-window MDS coords
- `inversion_localpca.candidate_regions.tsv.gz`    — coarse candidates
- `inversion_localpca.candidate_window_membership.tsv.gz`
- `inversion_localpca.mds_mode_metadata.tsv`
- `inversion_localpca.mds_background_<chr>.txt`    — chunked-mode bg indices

---

## Step G — `STEP_ZO_G_precompute.R`

Path-1 only. Per-chrom precompute, parallelized via `mclapply`.

**Inputs:**
- `<step_F_outprefix>.mds.rds` (output of step F)
- `--dosage_dir <dir>` (optional; per-chrom dosage TSVs from step B;
  enables `dosage_het_rate_cv` as the het component of `inv_likeness`)

**Outputs (in `<SCRATCH>/path_localpca_zblocks/04_precomp/`):**
- `precomp/<chr>.precomp.rds` — full precomp, columns include:
  - position: `global_window_id`, `chrom`, `start_bp`, `end_bp`, `center_bp`, `n_snps`
  - MDS: `MDS1..MDSk`, `MDS1_z..MDSk_z`, `max_abs_z`, `max_z_axis`
  - inv_likeness: `inv_likeness`, `inv_dip_p`, `inv_het_contrast`
  - band-shape: `band_discreteness`, `diffuse_score`, `het_intermediacy`,
    `n_effective_clusters`
  - multi-source het: `pc1_bimodality`, `het_pc1_gap`, `het_mid_fraction`,
    `het_mid_variance`, `dosage_het_rate_{median,sd,cv}`
  - local fixed-band: `local_delta12`, `local_entropy`, `local_ena`,
    `band1_frac`, `band2_frac`, `band3_frac`
  - adaptive: `beta_pval`, `adaptive_seed`, `beta_alpha`, `beta_beta`
  - seed eligibility: `seed_nn_dist`
  - morphology (z-profile, neighbourhood, sim_mat-block)
  - **per-sample**: `PC_1_Ind*`, `PC_2_Ind*` (these are what step M needs!)
- `precomp/sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds`
- `window_dt.tsv.gz` — genome-wide per-window scalar table
- `precomp_summary.tsv` — per-chrom summary

**Precomp contract (must hold):**
- Present columns: `max_abs_z`, `max_z_axis`, `inv_likeness`,
  `band_discreteness`, `diffuse_score`, `het_intermediacy`,
  `n_effective_clusters`, `dosage_het_rate_cv`, `beta_pval`,
  `adaptive_seed`, `seed_nn_dist`, `flat_inv_score`, `spiky_inv_score`,
  `fragmentation_score`, `PC_1_Ind*`
- Absent columns (would mean v9.x contamination is back):
  `sv_inv_*`, `localQ_*_K08`, `family_likeness`, `test05_fst_*`

**Sanity check after running:**
```bash
Rscript -e 'x <- readRDS("<SCRATCH>/path_localpca_zblocks/04_precomp/precomp/C_gar_LG28.precomp.rds")
            want <- c("max_abs_z","max_z_axis","inv_likeness")
            missing <- setdiff(want, names(x$dt))
            pc_cols <- length(grep("^PC_1_", names(x$dt)))
            cat("missing dims:", length(missing), "  PC_1_* cols:", pc_cols, "\n")'
```
Expected: `missing dims: 0  PC_1_* cols: 226`.

**No more slim.** Earlier sessions used `precomp.slim.rds` (with `PC_*_*`
columns dropped). It broke step M (PC2 jitter fallback). Slim is dead.

---

## Step H — `STEP_ZO_H_detect_L1.R`

Boundary-driven L1 envelope discovery, one chromosome at a time. The
boundary scan slides a small WxW upper-triangle cross-block at offset G
along the diagonal of the sim_mat; each peak is validated by the
**growing-W cross-block z** validator.

**Default NN scale: 80.**

**Canonical run (defaults from history line 19553):**

```bash
Rscript STEP_ZO_H_detect_L1.R \
  --precomp_dir <precomp>/precomp \
  --chr         <chr> \
  --outdir      <outdir> \
  --boundary_scan TRUE \
  --boundary_validator_mode grow \
  --boundary_W 5 --boundary_offset 5 --boundary_min_dist 30
```

**Inputs:** auto-resolved from `--precomp_dir + --chr + --nn` (or pass
`--precomp <path>` and `--sim_mat <path>` to override).

**Outputs (in `<outdir>`):**
- `<chr>.L1_envelopes.tsv`     — L1 segments derived from STABLE_BLUE cuts
- `<chr>.L1_boundaries.tsv`    — every detected peak with validation status
  - `STABLE_BLUE` (real long-range separation, persists at all val_W)
  - `STABLE_RED` (rejected by grow validator)
  - `WEAK_BLUE`, `DECAYS`, `MARGINAL`, `EDGE`
- `<chr>.L1_score_curve.tsv`   — per-window boundary score curve

**Filter `--boundary_filter stable` in steps I/K shows only `STABLE_BLUE`.**

**Settled parameters (from earlier sessions):**
- `--boundary_min_dist 30` — sweep through 10/20/30 (history 19488→19510)
- `--boundary_validator_mode grow` — default; `both` adds perpendicular ray
- merge threshold raised to 3 in v7 (history 19541)

---

## Step I — `STEP_ZO_I_plot_L1.R`

Multi-page L1 PDF.
- Page 1: whole-chromosome heatmap (lower triangle = sim, upper = z) with
  L1 envelopes outlined and STABLE_BLUE boundaries marked with red squares
- Pages 2..(N+1): one zoomed page per L1 segment, with **local Z** computed
  from the segment's own sub-matrix (so colours are calibrated per segment)

**Default NN scale: 80.**

**Canonical run:**

```bash
Rscript STEP_ZO_I_plot_L1.R \
  --precomp_dir <precomp>/precomp \
  --L1_dir      <l1_outdir> \
  --chr         <chr> \
  --outdir      <plot_outdir> \
  --toggle_L1 yes --boundary_filter stable
```

**Output:** `<chr>.L1_overlay.pdf` (single multi-page PDF).

**Style locked-in over many sessions:**
- L1 envelope outlines: deep blue `#1F3A6E`
- L1 envelope label: deep blue with white halo via ggrepel
- Boundary squares filled with status color (STABLE_BLUE = red)

---

## Step J — `STEP_ZO_J_detect_L2.R`

Reads L1 envelopes, runs the same boundary scan INDEPENDENTLY INSIDE EACH
L1 SEGMENT, producing a finer L2 partition.

**Default NN scale: 40** (finer than L1's 80 because we're zooming inside
a segment — smoother sim_mats wash out sub-structure here).

**v8 features:**
- **Ward-adaptive thresholding**: each L1 segment gets its own
  (real_max, fake_min) thresholds derived from its `internal_homogeneity`
  and `architectural_complexity` priors. STRICT segments (looks like inside
  one inversion) get tighter bounds; LOOSE (real architectural transitions
  inside the L1 region) get wider bounds. `--ward_adaptive no` falls back
  to pooled thresholds.
- **Quadrant validator**: a 2D test that checks whether a candidate L2
  boundary inside an L1 segment is consistent with a true sub-block, by
  testing the four sim_mat quadrants the boundary divides. Enable with
  `--quadrant_validator yes`.

**Canonical run (history line 19585):**

```bash
Rscript STEP_ZO_J_detect_L2.R \
  --precomp_dir <precomp>/precomp \
  --L1_dir      <l1_outdir> \
  --chr         <chr> \
  --outdir      <l2_outdir> \
  --boundary_scan TRUE --boundary_validator_mode grow \
  --quadrant_validator yes \
  --weak_demote_score 0 \
  --quad_rescue_max_grow_z 1.5 \
  --quad_demote_on_fail yes \
  --quad_demote_drift_floor -1.0
```

**Outputs (in `<l2_outdir>`):**
- `<chr>.L2_envelopes.tsv`           — L2 sub-segments
- `<chr>.L2_boundaries.tsv`          — all detected peaks with status
- `<chr>.L2_segment_stats.tsv`       — per-L1-segment Ward stats + strictness
- `<chr>.L2_quadrant_validator.tsv`  — per-peak quadrant test (when enabled)
- `<chr>.L2_quadrant_audit.tsv`      — extra quadrant audit

---

## Step K — `STEP_ZO_K_plot_L2.R`

Multi-page PDF for the L1+L2 hierarchy.
- Page 1: same as step I page 1 (chromosome-wide nn80 with L1 boundaries)
- Pages 2..(N+1): per L1 segment, using **nn40** (`--sim_mat_l2`), with
  L2 boundaries inside that segment overlaid as red squares

**Default NN scales: 80 (chrom) + 40 (inside).**

**Canonical run:**

```bash
Rscript STEP_ZO_K_plot_L2.R \
  --precomp_dir <precomp>/precomp \
  --L1_dir      <l1_outdir> \
  --L2_dir      <l2_outdir> \
  --chr         <chr> \
  --outdir      <plot_outdir> \
  --boundary_filter stable
```

**Output:** `<chr>.L2_overlay.pdf`

---

## Step L — `STEP_ZO_L_build_sample_metadata.R`

Reconciles the three identity layers into ONE merged TSV, ONCE,
genome-wide.

**Three layers:**
1. **bamlist** — one CGA per line, in PC_1_Ind* order → maps `Ind0..IndN → CGAxxx`
2. **ngsRelate pairs** — 3-col `id1 id2 theta` → union-find with cutoff
   (default 0.177, Manichaikul 1st-degree) → `family_id`
3. **ancestry** — TSV with `cga<TAB>ancestry` from NGSadmix K=8 → `ancestry`

**Canonical run:**

```bash
Rscript STEP_ZO_L_build_sample_metadata.R \
  --bamlist     <path>/list_of_samples.tsv \
  --pairs       <path>/catfish_226_for_natora.txt \
  --theta_cutoff 0.177 \
  --ancestry    <path>/ngsadmix_K8_ancestry.tsv \
  --out         <SCRATCH>/_shared/sample_metadata.tsv
```

**Output (in `<SCRATCH>/_shared/sample_metadata.tsv`):**

| Column     | Description                                       |
|------------|---------------------------------------------------|
| `ind`      | `Ind0..IndN` (matches `PC_1_Ind*` column suffix)  |
| `cga`      | Real sample ID (after bamlist remap, no .bam ext) |
| `family_id`| Integer; isolates = -1                             |
| `ancestry` | Ancestry label (default "unknown")                 |

This same TSV will be reusable by paths 2 and 3 once they're refactored —
all three discovery paths use the same 226-sample cohort.

---

## Step M — `STEP_ZO_M_export_atlas_json.R`

The FINAL JSON. Embeds:
- multi-scale sim_mat thumbnails (`nn40`, `nn80`, `nn160`, `nn320`)
- per-distance local Z scores (matching the overlay PDF colorimetry)
- L1 envelopes + L1 boundaries
- L2 envelopes + L2 boundaries
- sample identity from `sample_metadata.tsv`

**Canonical run (with `--sample_metadata`):**

```bash
Rscript STEP_ZO_M_export_atlas_json.R \
  --precomp_dir     <SCRATCH>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir          <SCRATCH>/path_localpca_zblocks/05_L1 \
  --L2_dir          <SCRATCH>/path_localpca_zblocks/07_L2 \
  --chr             C_gar_LG28 \
  --sample_metadata <SCRATCH>/_shared/sample_metadata.tsv \
  --out             <SCRATCH>/path_localpca_zblocks/09_atlas_json/C_gar_LG28.atlas.json
```

**Legacy mode:** if `--sample_metadata` not given, falls back to the old
`--bamlist + --pairs + --samples` 3-flag interface for backward compat.

**Output:** `<chr>.atlas.json` — drag-and-drop into the atlas page-1
viewer / pca_scrubber.

---

## Glossary of file extensions

| Extension                                 | What it is                                       |
|-------------------------------------------|--------------------------------------------------|
| `<chr>.beagle.gz`                         | ANGSD GL output (input)                          |
| `<chr>.dosage.tsv.gz`                     | Step B output                                    |
| `<chr>.sites.tsv.gz`                      | Step B output                                    |
| `<chr>.window_pca.rds`                    | Step D output (top-npc eigvecs + full scree)     |
| `inversion_localpca.mds.rds`              | Step F output (`$per_chr` structure)             |
| `<chr>.precomp.rds`                       | Step G output (full, with PC_1_*/PC_2_* columns) |
| `<chr>.sim_mat_nn{0..320}.rds`            | Step G output                                    |
| `<chr>.L1_envelopes.tsv`                  | Step H output                                    |
| `<chr>.L1_boundaries.tsv`                 | Step H output                                    |
| `<chr>.L1_score_curve.tsv`                | Step H output                                    |
| `<chr>.L1_overlay.pdf`                    | Step I output                                    |
| `<chr>.L2_envelopes.tsv`                  | Step J output                                    |
| `<chr>.L2_boundaries.tsv`                 | Step J output                                    |
| `<chr>.L2_segment_stats.tsv`              | Step J output                                    |
| `<chr>.L2_quadrant_validator.tsv`         | Step J output (when --quadrant_validator yes)    |
| `<chr>.L2_overlay.pdf`                    | Step K output                                    |
| `sample_metadata.tsv`                     | Step L output (genome-wide one-off, in _shared/) |
| `<chr>.atlas.json`                        | Step M output (atlas/scrubber JSON)              |
