# `local_PCA_MDS_z/` — HANDOFF (flat layout, May 2026)

> Path 1 (local-PCA z-blocks) of catfish-inversion-analysis. Documents:
> 1. The **script tree** (this folder, flat layout)
> 2. The **scratch tree** on LANTA where outputs land
> 3. The **canonical command** for every step
> 4. The **conventions** (defaults, naming, prefix scheme)
> 5. **What's left** for future sessions
>
> All step scripts and SLURM launchers live flat at the top level of this
> folder, prefixed `STEP_ZO_<LETTER>_` (Z Outlier path, lettered B–M in
> pipeline order — A is reserved for any future upstream cohort-prep step).
> The shared config template lives at the repo root:
> `../00_inversion_config.sh.template`.

---

## 1. The two trees

### 1.1 Script tree (this repo)

Flat — no numbered subfolders. Each step's script and its SLURM launcher
live next to each other, both prefixed `STEP_ZO_<LETTER>_`.

```
catfish-inversion-analysis/
├── 00_inversion_config.sh.template          (shared across all 3 paths)
├── local_PCA_MDS_z/                         (this folder, path 1)
│   ├── STEP_ZO_B_beagle_to_dosage.py        (beagle.gz → dosage + sites)
│   ├── STEP_ZO_B_LAUNCH_beagle_to_dosage.slurm
│   ├── STEP_ZO_C_local_pca_compute.R        (sliding-window PCA, ARRAY per chrom)
│   ├── STEP_ZO_C_LAUNCH_local_pca_compute.slurm
│   ├── STEP_ZO_D_local_pca_merge.R          (assigns global window_ids)
│   ├── STEP_ZO_D_LAUNCH_local_pca_merge.slurm
│   ├── STEP_ZO_E_mds_compute.R              (lostruct + cmdscale, ARRAY per focal chrom)
│   ├── STEP_ZO_E_LAUNCH_mds_compute.slurm
│   ├── STEP_ZO_F_mds_merge.R                (assembles per_chr + candidate regions)
│   ├── STEP_ZO_F_LAUNCH_mds_merge.slurm
│   ├── STEP_ZO_G_precompute.R               (per-chrom features + NN sim_mats, mclapply)
│   ├── STEP_ZO_G_LAUNCH_precompute.slurm
│   ├── STEP_ZO_H_detect_L1.R                (default --nn 80)
│   ├── STEP_ZO_I_plot_L1.R                  (default --nn 80)
│   ├── STEP_ZO_J_detect_L2.R                (default --nn 40)
│   ├── STEP_ZO_K_plot_L2.R                  (default --nn 80 chrom + --nn_l2 40 inside)
│   ├── STEP_ZO_L_build_sample_metadata.R    (bamlist + ngsRelate + ancestry → ONE TSV)
│   ├── STEP_ZO_M_export_atlas_json.R        (per-chrom JSON for the scrubber)
│   ├── STEP_ZO_HtoM_run_one_chrom.sh        (driver chains H→M for one chrom)
│   ├── docs/PER_STEP_NOTES.md
│   ├── README.md
│   └── HANDOFF.md                           (this file)
├── local_PCA_MDS_thetapi/                  (path 2, sibling — STEP_TR_* prefix)
└── local_PCA_MDS_GHSL/                     (path 3, sibling — STEP_GH_* prefix)
```

### 1.2 Scratch tree on LANTA

The **scratch-tree folder names keep the original 01–09 numbering** —
those are output buckets, not script names. Only scripts changed.

```
${SCRATCH}/inversion_localpca_v8/
│
├── 01_beagle/                                ← ANGSD output (input)
│   └── <chr>.beagle.gz
│
├── 02_dosage_sites/                          ← step B output (SHARED upstream)
│   ├── <chr>.dosage.tsv.gz
│   └── <chr>.sites.tsv.gz
│
├── 03_pestPG/                                ← ANGSD -doThetas (path 2 input)
│   └── <chr>.pestPG  (or whatever ANGSD names it)
│
├── 04_clair3_phased_GHSL/                    ← Clair3 phased haplotypes (path 3 input)
│   └── <chr>.<phased haplotype output>
│
├── path_localpca_zblocks/                    ← path 1 outputs (THIS pipeline)
│   ├── 01_local_pca/                          (step C: tmp/<chr>.window_pca_tmp.rds)
│   ├── 02_dense_registry/                     (step D: <chr>.window_pca.rds, master.tsv.gz)
│   ├── 03_mds/                                (steps E/F: inversion_localpca.mds.rds)
│   ├── 04_precomp/                            (step G: precomp/<chr>.precomp.rds + sim_mats/)
│   ├── 05_L1/                                 (step H: <chr>.L1_envelopes/_boundaries/_score_curve.tsv)
│   ├── 06_L1_plots/                           (step I: <chr>.L1_overlay.pdf)
│   ├── 07_L2/                                 (step J: <chr>.L2_envelopes/_boundaries/_segment_stats/_quadrant_validator.tsv)
│   ├── 08_L2_plots/                           (step K: <chr>.L2_overlay.pdf)
│   └── 09_atlas_json/                         (step M: <chr>.atlas.json)
│
├── path_localpca_thetapi/                    ← path 2 outputs (SAME 01-09 layout)
├── path_localpca_GHSL/                       ← path 3 outputs (SAME 01-09 layout)
│
└── _shared/
    ├── sample_metadata.tsv                  ← step L output, consumed by all 3 path exporters
    ├── 00_inversion_config.sh
    └── reference/
        └── fClaHyb_Gar_LG.fa
```

**Why script-letter naming doesn't match scratch-folder numbering inside
each `path_*` folder:** the root has 4 upstream-shared folders eating the
first 4 numbers (`01_beagle`, `02_dosage_sites`, `03_pestPG`,
`04_clair3_phased_GHSL`). Inside `path_localpca_zblocks/` the numbering
restarts at `01_local_pca/`. The script-to-scratch mapping is:

| Script (letter) | Output folder                                  |
|-----------------|------------------------------------------------|
| `B`             | `02_dosage_sites/`  (root, shared)             |
| `C`             | `path_localpca_zblocks/01_local_pca/`          |
| `D`             | `path_localpca_zblocks/02_dense_registry/`     |
| `E`+`F`         | `path_localpca_zblocks/03_mds/`                |
| `G`             | `path_localpca_zblocks/04_precomp/`            |
| `H` detect_L1   | `path_localpca_zblocks/05_L1/`                 |
| `I` plot_L1     | `path_localpca_zblocks/06_L1_plots/`           |
| `J` detect_L2   | `path_localpca_zblocks/07_L2/`                 |
| `K` plot_L2     | `path_localpca_zblocks/08_L2_plots/`           |
| `L`             | `_shared/sample_metadata.tsv`                  |
| `M`             | `path_localpca_zblocks/09_atlas_json/`         |

Path 2 and path 3, when refactored, will have NO step-B equivalent (their
features come from upstream pipelines outside this toolkit), so they start
at `01_local_pca/` with no offset.

---

## 2. The data flow — copy/paste runnable on LANTA

All commands assume the canonical scratch tree above. The default config
file lives at `${SCRATCH}/inversion_localpca_v8/_shared/00_inversion_config.sh`.

### Step B — beagle.gz → dosage + sites

(Per-chrom; one task per chromosome line in `chrom.list`.)

Currently the launcher script `STEP_ZO_B_LAUNCH_beagle_to_dosage.slurm` is
**unchanged from before** — it predates the consolidated layout. Run it as
you have been; the output should land in `02_dosage_sites/`.

### Steps C/D — per-chrom dense local PCA (compute + merge)

```bash
# Stage 1: one chromosome per array task
sbatch --array=0-27 STEP_ZO_C_LAUNCH_local_pca_compute.slurm chrom.list

# Stage 2: merge (assigns global window_ids)
sbatch --dependency=afterok:<JOB_C> STEP_ZO_D_LAUNCH_local_pca_merge.slurm
```

Produces in `path_localpca_zblocks/02_dense_registry/`:
- `windows_master.tsv.gz` — the master window registry
- `<chr>.window_pca.rds` — per-chrom local PCA, top-`npc=4` eigvecs + full scree spectrum
- `<chr>.window_pca.tsv.gz`

### Steps E/F — lostruct distance + MDS (compute + merge)

```bash
# Stage 1: one focal chromosome per array task (~12h walltime)
sbatch --array=0-27 STEP_ZO_E_LAUNCH_mds_compute.slurm chrom.list

# Stage 2: merge into final mds.rds with $per_chr structure
sbatch --dependency=afterok:<JOB_E> STEP_ZO_F_LAUNCH_mds_merge.slurm
```

Produces in `path_localpca_zblocks/03_mds/`:
- `inversion_localpca.mds.rds` (with `$per_chr` field every downstream consumes)
- `inversion_localpca.window_mds.tsv.gz`
- `inversion_localpca.candidate_regions.tsv.gz`

Default mode: `chunked_2x` — each focal chrom is MDS'd against itself plus
2× background sampled from non-focal chromosomes (excluding high
inv_likeness windows so foreign inversions don't leak into "background").

### Step G — precompute (per-chrom features + NN sim_mats)

```bash
sbatch --dependency=afterok:<JOB_F> STEP_ZO_G_LAUNCH_precompute.slurm
```

Produces in `path_localpca_zblocks/04_precomp/`:
- `precomp/<chr>.precomp.rds` — full precomp (with `PC_1_*` and `PC_2_*` per-sample columns)
- `precomp/sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds`
- `window_dt.tsv.gz`, `precomp_summary.tsv`

**Tune NN scales saved:**
```bash
NN_SIM_SCALES="40,80,160,320" sbatch STEP_ZO_G_LAUNCH_precompute.slurm
```

### Steps H–M — interactive per chromosome

The bulk-compute steps are B–G. Steps H–M are tuned interactively
(parameter sweeps) and then rolled out per chromosome with canonical
defaults baked in. The driver:

```bash
bash STEP_ZO_HtoM_run_one_chrom.sh C_gar_LG28
```

Or call individual scripts with the canonical-mode flags:

```bash
# H — L1 detect (defaults: --nn 80, boundary_W 5, boundary_offset 5,
#                          boundary_min_dist 30, validator_mode grow)
Rscript STEP_ZO_H_detect_L1.R \
  --precomp_dir <SCRATCH>/path_localpca_zblocks/04_precomp/precomp \
  --chr         C_gar_LG28 \
  --outdir      <SCRATCH>/path_localpca_zblocks/05_L1 \
  --boundary_scan TRUE \
  --boundary_validator_mode grow \
  --boundary_W 5 --boundary_offset 5 --boundary_min_dist 30
```

```bash
# I — L1 plot (default --nn 80, --boundary_filter stable)
Rscript STEP_ZO_I_plot_L1.R \
  --precomp_dir <SCRATCH>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir      <SCRATCH>/path_localpca_zblocks/05_L1 \
  --chr         C_gar_LG28 \
  --outdir      <SCRATCH>/path_localpca_zblocks/06_L1_plots \
  --toggle_L1 yes --boundary_filter stable
```

```bash
# J — L2 detect (default --nn 40, with quadrant validator)
Rscript STEP_ZO_J_detect_L2.R \
  --precomp_dir <SCRATCH>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir      <SCRATCH>/path_localpca_zblocks/05_L1 \
  --chr         C_gar_LG28 \
  --outdir      <SCRATCH>/path_localpca_zblocks/07_L2 \
  --boundary_scan TRUE --boundary_validator_mode grow \
  --quadrant_validator yes \
  --weak_demote_score 0 \
  --quad_rescue_max_grow_z 1.5 \
  --quad_demote_on_fail yes \
  --quad_demote_drift_floor -1.0
```

```bash
# K — L2 plot (default --nn 80 chrom-wide + --nn_l2 40 inside-segment)
Rscript STEP_ZO_K_plot_L2.R \
  --precomp_dir <SCRATCH>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir      <SCRATCH>/path_localpca_zblocks/05_L1 \
  --L2_dir      <SCRATCH>/path_localpca_zblocks/07_L2 \
  --chr         C_gar_LG28 \
  --outdir      <SCRATCH>/path_localpca_zblocks/08_L2_plots \
  --boundary_filter stable
```

```bash
# L — build sample metadata (run ONCE, genome-wide)
Rscript STEP_ZO_L_build_sample_metadata.R \
  --bamlist    <path>/list_of_samples_one_per_line_same_bamfile_list.tsv \
  --pairs      <path>/catfish_226_for_natora.txt \
  --theta_cutoff 0.177 \
  --ancestry   <path>/ngsadmix_K8_ancestry.tsv \
  --out        <SCRATCH>/_shared/sample_metadata.tsv
```

```bash
# M — export atlas JSON per chromosome
Rscript STEP_ZO_M_export_atlas_json.R \
  --precomp_dir     <SCRATCH>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir          <SCRATCH>/path_localpca_zblocks/05_L1 \
  --L2_dir          <SCRATCH>/path_localpca_zblocks/07_L2 \
  --chr             C_gar_LG28 \
  --sample_metadata <SCRATCH>/_shared/sample_metadata.tsv \
  --out             <SCRATCH>/path_localpca_zblocks/09_atlas_json/C_gar_LG28.atlas.json
```

---

## 3. Conventions baked in

### 3.1 NN scale defaults

| Step    | Default NN | Override flag |
|---------|------------|---------------|
| H detect_L1 | 80         | `--nn N`        |
| I plot_L1   | 80         | `--nn N`        |
| J detect_L2 | 40         | `--nn N`        |
| K plot_L2   | 80 chrom-wide + 40 inside-segment | `--nn N` + `--nn_l2 N` |
| M atlas     | 40, 80, 160, 320 (multi-scale) | `--nn_list 40,80,160,320` |

L1 uses **nn80** because chromosome-wide it needs to suppress short-range
noise. L2 uses **nn40** because inside an L1 segment, finer resolution is
the goal. (History lines 19553 / 19556.)

### 3.2 Path resolution flags (canonical mode)

All four detect/plot scripts and step M accept:

- `--precomp_dir <dir>` — auto-resolves `<chr>.precomp.rds` and the right
  `<chr>.sim_mat_nn{N}.rds` from the directory + `--chr` + `--nn`
- `--L1_dir <dir>` — auto-resolves `<chr>.L1_envelopes.tsv` and `<chr>.L1_boundaries.tsv`
- `--L2_dir <dir>` — auto-resolves the L2 equivalents
- `--chr <label>` — chromosome to operate on

Power users can still pass `--precomp <path>`, `--sim_mat <path>`,
`--catalogue <path>` etc. directly to override auto-resolution. Useful for
parameter sweeps where you want a custom NN scale or a non-canonical L1
catalogue.

### 3.3 Artifact naming (no more `d17`)

| Old name (retired) | New name |
|---|---|
| `<chr>_d17L1_envelopes.tsv`         | `<chr>.L1_envelopes.tsv` |
| `<chr>_d17L1_boundaries.tsv`        | `<chr>.L1_boundaries.tsv` |
| `<chr>_d17L1_boundary_score_curve.tsv` | `<chr>.L1_score_curve.tsv` |
| `<chr>_d17L1_overlay.pdf`           | `<chr>.L1_overlay.pdf` |
| `<chr>_d17L2_envelopes.tsv`         | `<chr>.L2_envelopes.tsv` |
| `<chr>_d17L2_boundaries.tsv`        | `<chr>.L2_boundaries.tsv` |
| `<chr>_d17L2_segment_stats.tsv`     | `<chr>.L2_segment_stats.tsv` |
| `<chr>_d17L2_quadrant_validator.tsv` | `<chr>.L2_quadrant_validator.tsv` |
| `<chr>_d17L2_quadrant_audit.tsv`    | `<chr>.L2_quadrant_audit.tsv` |
| `<chr>_d17L2_overlay.pdf`           | `<chr>.L2_overlay.pdf` |

The `d17` prefix was carried over from `STEP_D17_*` session names. It
carried no information about what the file was. It's gone.

### 3.4 Slim precomp — KILLED

Earlier versions used `<chr>.precomp.slim.rds` (sample-level columns
dropped to keep file ~10 MB). Steps H–K only need window coordinates so
slim worked for them, but step M NEEDS `PC_1_*` and `PC_2_*` per-sample
columns (lost in slim → falls back to PC2 jitter).

**Decision: full precomp only.** All scripts read the same
`<chr>.precomp.rds` (~70-100 MB, includes per-sample PCs). Steps H–K
only touch window coords; the extra columns are harmless. Step M gets
real PC1+PC2.

The `prep_lg28_bundle.R` workaround that built slim copies is retired.

### 3.5 Identity reconciliation — separated and centralized

The three identity layers (bamlist remap, ngsRelate family graph, NGSadmix
ancestry) used to be mashed together inside the JSON exporter via three
separate flags (`--bamlist`, `--pairs`, `--samples`).

**New design: `STEP_ZO_L_build_sample_metadata.R`** consumes the three
independent inputs and produces ONE merged TSV (`sample_metadata.tsv`,
columns: `ind, cga, family_id, ancestry`). The exporter `M` then takes a
single `--sample_metadata` flag.

Benefits:
- Reconciliation logic lives in ONE place
- The merged TSV can be sanity-checked before paying JSON-build cost
- Same TSV reusable by paths 2 and 3 (all three discovery paths use the
  same 226 samples) — `_shared/sample_metadata.tsv` is the single source
  of truth

Legacy `--bamlist` + `--pairs` + `--samples` flags still work in M for
power users; if `--sample_metadata` is given it short-circuits the legacy
logic.

### 3.6 Naming history

Three earlier rename passes are baked into the current `STEP_ZO_<LETTER>_`
names:

- **stage1/stage2 → compute/merge** for the C/D and E/F pairs (the
  SLURM-array + post-array merge halves).
- **overlay → plot** for `I plot_L1` and `K plot_L2` (matches the
  `detect_*` naming).
- **numbered (01a..08b) → letters (B..M)** in the May-2026 pass, for
  consistency with the sibling theta-pi (`STEP_TR_B..G`) and band-pi
  (`STEP_PI_B..F`) modules. A is reserved for any future upstream
  cohort-prep step.

All folded into the current scheme. No script names in this repo still
carry the old `_stage1` / `_stage2` / `_overlay` / `_localpca_zblocks`
tags or the numeric `01a/01b/.../08b` prefixes.

---

## 4. The three discovery paths — symmetry

All three paths share this skeleton:

```
local_pca_compute  (per-chr ARRAY, sliding-window PCA on path-specific feature)
local_pca_merge    (global window_ids)
mds_compute + merge (lostruct + cmdscale, with chunked background)
precompute         (per-chrom features + NN sim_mats)
detect_L1          (chrom-wide z-block detection, nn80)
plot_L1            (multi-page PDF)
detect_L2          (per-L1-segment fine sub-block detection, nn40)
plot_L2            (multi-page PDF)
build_sample_metadata + export_atlas_json
                   (consume sample_metadata.tsv from _shared/, emit JSON)
```

Path 1 (this folder) has an extra `STEP_ZO_B_beagle_to_dosage.py` because
dosage is computed inside this toolkit. Paths 2 and 3 inherit their feature
matrices from upstream pipelines (ANGSD `-doThetas` for path 2; Clair3
phasing for path 3) and skip directly to `local_pca_compute` (their own
"step C equivalent").

When paths 2 and 3 are refactored later, they'll match the same structure,
just with different feature inputs feeding the local PCA. The MDS,
precompute, detect, plot, and atlas-JSON code is fully shared in design (and
maybe even partly in implementation — the L1/L2 detect logic operates on
sim_mats regardless of what the feature was).

---

## 5. What's left for future sessions

### High priority

1. **Path 2 (θπ) — finalize the v5 drafts at
   `local_PCA_MDS_thetapi/v5_drafts/`.** TR_C / TR_D / TR_C_plot /
   TR_D_plot are verbatim copies of `STEP_ZO_H..K` and should run on a
   theta-pi precomp out of the box. TR_E (carrier classification) and
   TR_F (per-sample + band-mean CUSUM) are exploratory — needs a
   real-data dry-run on LG28 before flattening into top-level
   `STEP_TR_*` scripts at the sibling folder root.

2. **`STEP_ZO_B_LAUNCH_beagle_to_dosage.slurm`** — predates the flatten;
   its config sourcing may still reference an old `_shared/` path, and its
   actual python invocation calls `STEP_A01_beagle_to_dosage_by_chr.py`
   from a different `${A2_STEPS}` directory rather than the local
   `STEP_ZO_B_beagle_to_dosage.py`. Update to source
   `../00_inversion_config.sh` (or your site-local copy of the template)
   AND switch the python call to the local script.

3. **Path 3 (GHSL) — flatten next.** Apply the same `STEP_GH_*` prefix
   scheme used here (B onwards, no numbers).

### Medium priority

4. **`STEP_ZO_L_build_sample_metadata.R` — test on real LANTA data.**
   Parses cleanly but has not been run end-to-end on the real bamlist +
   ngsRelate output. Verify against an existing manually-built
   sample_metadata before trusting it.

5. **Plug real values into `00_inversion_config.sh.template`** at the
   repo root. The template has the right shape but every path / env name
   / sample count is still a placeholder. `BASE`, `SCRATCH_ROOT`, and the
   `_shared/` location depend on your cluster setup.

### Low priority / deferred

6. **Genome-wide rollout SLURM array for steps H–M.** Currently steps
   H–K are interactive; the `STEP_ZO_HtoM_run_one_chrom.sh` driver
   runs them for one chrom at a time. A SLURM array launcher that runs
   the driver in parallel across 28 chromosomes would be useful for
   genome-wide atlas builds.

7. **`STEP_ZO_B_beagle_to_dosage.py` sanity pass.** Output naming
   (`<chr>.dosage.tsv.gz` + `<chr>.sites.tsv.gz`) matches the new
   convention, so it should be fine, but worth a quick review.

### Done in the May-2026 cleanup session

- Numbered subfolders (`01_dosage_pca/` … `08_atlas_json/`) removed; all
  scripts flattened to top level with the `STEP_ZO_` prefix.
- `99_launchers/` and `99_legacy/` removed; launchers also flattened
  alongside their scripts. `99_docs/` renamed to `docs/`.
- Long file-name suffix `_localpca_zblocks` dropped from script names
  (the prefix already says it).
- `00_inversion_config.sh.template` moved out to the repo root for
  sharing with sibling paths.
- Internal `${SCRIPT_DIR}/<old-folder>/<old-name>` references in every
  launcher patched to point at the new flat names.
- **Numeric step prefixes (01a..08b) standardized to letters (B..M)** for
  consistency with the sibling theta-pi (B–G) and band-pi (B–F) modules.

---

## 6. The full DAG

```
                         (per chromosome)
beagle.gz ──► B beagle_to_dosage ──► dosage + sites
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ C local_pca_compute (ARRAY)          │
                       │ D local_pca_merge   (single)         │
                       └──────────────────────────────────────┘
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ E mds_compute (ARRAY focal-chr)      │
                       │ F mds_merge   (single)               │
                       └──────────────────────────────────────┘
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ G precompute_localpca_zblocks        │
                       │   precomp/<chr>.precomp.rds          │
                       │   sim_mats/<chr>.sim_mat_nn{0..320}  │
                       └──────────────────────────────────────┘
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ H detect_L1 (nn80)                   │
                       │   <chr>.L1_envelopes.tsv             │
                       │   <chr>.L1_boundaries.tsv            │
                       │   <chr>.L1_score_curve.tsv           │
                       └──────────────────────────────────────┘
                                │                 │
                                │                 ▼
                                │        ┌──────────────────┐
                                │        │ I plot_L1        │
                                │        │ <chr>.L1_overlay │
                                │        └──────────────────┘
                                ▼
                       ┌──────────────────────────────────────┐
                       │ J detect_L2 (nn40, per L1 segment)   │
                       │   <chr>.L2_envelopes.tsv             │
                       │   <chr>.L2_boundaries.tsv            │
                       │   <chr>.L2_segment_stats.tsv         │
                       │   <chr>.L2_quadrant_validator.tsv    │
                       └──────────────────────────────────────┘
                                │                 │
                                │                 ▼
                                │        ┌──────────────────┐
                                │        │ K plot_L2        │
                                │        │ <chr>.L2_overlay │
                                │        └──────────────────┘
                                ▼
   bamlist ┐                    │
   pairs   ├─► L build_sample_metadata ──► sample_metadata.tsv
   ancestry┘  (genome-wide ONCE)           (in _shared/)
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ M export_atlas_json                  │
                       │   inputs: full precomp + L1+L2 +     │
                       │           sample_metadata.tsv +      │
                       │           sim_mats nn{40,80,160,320} │
                       │   output: <chr>.atlas.json           │
                       └──────────────────────────────────────┘
                                          │
                                          ▼
                                pca_scrubber_v3 / atlas page-1
```
