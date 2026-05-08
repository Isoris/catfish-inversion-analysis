# `local_PCA_z/` — HANDOFF (consolidated layout v1.0, May 2026)

> Single canonical folder for path 1 (local-PCA z-blocks) of the
> inversion-popgen-toolkit. This handoff documents:
> 1. The **script tree** (this folder)
> 2. The **scratch tree** on LANTA where outputs land
> 3. The **canonical command** for every step
> 4. The **conventions** (defaults, naming, directory structure)
> 5. **What's left** for future sessions
>
> Codebase: `inversion-popgen-toolkit` v8.5 / consolidated layout v1.0.

---

## 1. The two trees

### 1.1 Script tree (this repo)

```
local_PCA_z/
├── 01_dosage_pca/
│   ├── 01a_beagle_to_dosage.py            (beagle.gz → dosage + sites)
│   ├── 01b_local_pca_compute.R            (sliding-window PCA, ARRAY per chrom)
│   └── 01c_local_pca_merge.R              (assigns global window_ids)
├── 02_mds/
│   ├── 02a_mds_compute.R                  (lostruct + cmdscale, ARRAY per focal chrom)
│   └── 02b_mds_merge.R                    (assembles per_chr + candidate regions)
├── 03_precomp/
│   └── 03_precompute_localpca_zblocks.R   (per-chrom features + NN sim_mats, mclapply)
├── 04_detect_L1/
│   └── 04_detect_L1_localpca_zblocks.R    (default --nn 80)
├── 05_plot_L1/
│   └── 05_plot_L1_localpca_zblocks.R      (default --nn 80)
├── 06_detect_L2/
│   └── 06_detect_L2_localpca_zblocks.R    (default --nn 40)
├── 07_plot_L2/
│   └── 07_plot_L2_localpca_zblocks.R      (default --nn 80 chrom + --nn_l2 40 inside)
├── 08_atlas_json/
│   ├── 08a_build_sample_metadata.R        (bamlist + ngsRelate + ancestry → ONE TSV)
│   └── 08b_export_atlas_json_localpca_zblocks.R
├── 99_launchers/                           (SLURM wrappers + one-shot driver)
├── 99_legacy/                              (retired _legacy.R + monolithic versions)
├── 99_docs/                                (PER_STEP_NOTES.md, schema specs)
├── _to_be_reorganized/                     (path 2 + path 3 code, NOT YET refactored)
├── README.md
└── HANDOFF.md                              (this file)
```

### 1.2 Scratch tree on LANTA

```
${SCRATCH}/inversion_localpca_v8/
│
├── 01_beagle/                                ← ANGSD output (input)
│   └── <chr>.beagle.gz
│
├── 02_dosage_sites/                          ← step 01a output (SHARED upstream)
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
│   ├── 01_local_pca/                          (step 01b: tmp/<chr>.window_pca_tmp.rds)
│   ├── 02_dense_registry/                     (step 01c: <chr>.window_pca.rds, master.tsv.gz)
│   ├── 03_mds/                                (steps 02: inversion_localpca.mds.rds)
│   ├── 04_precomp/                            (step 03: precomp/<chr>.precomp.rds + sim_mats/)
│   ├── 05_L1/                                 (step 04: <chr>.L1_envelopes/_boundaries/_score_curve.tsv)
│   ├── 06_L1_plots/                           (step 05: <chr>.L1_overlay.pdf)
│   ├── 07_L2/                                 (step 06: <chr>.L2_envelopes/_boundaries/_segment_stats/_quadrant_validator.tsv)
│   ├── 08_L2_plots/                           (step 07: <chr>.L2_overlay.pdf)
│   └── 09_atlas_json/                         (step 08b: <chr>.atlas.json)
│
├── path_localpca_thetapi/                    ← path 2 outputs (SAME 01-09 layout)
├── path_localpca_GHSL/                       ← path 3 outputs (SAME 01-09 layout)
│
└── _shared/
    ├── sample_metadata.tsv                  ← step 08a output, consumed by all 3 path exporters
    ├── 00_inversion_config.sh
    └── reference/
        └── fClaHyb_Gar_LG.fa
```

**Why script-folder numbering doesn't match scratch-folder numbering inside
each `path_*` folder:** the root has 4 upstream-shared folders eating the
first 4 numbers (`01_beagle`, `02_dosage_sites`, `03_pestPG`,
`04_clair3_phased_GHSL`). Inside `path_localpca_zblocks/` the numbering
restarts at `01_local_pca/`. The script-to-scratch mapping is:

| Script         | Output folder                                  |
|----------------|------------------------------------------------|
| `01a`          | `02_dosage_sites/`  (root, shared)             |
| `01b`          | `path_localpca_zblocks/01_local_pca/`          |
| `01c`          | `path_localpca_zblocks/02_dense_registry/`     |
| `02a`+`02b`    | `path_localpca_zblocks/03_mds/`                |
| `03`           | `path_localpca_zblocks/04_precomp/`            |
| `04` detect_L1 | `path_localpca_zblocks/05_L1/`                 |
| `05` plot_L1   | `path_localpca_zblocks/06_L1_plots/`           |
| `06` detect_L2 | `path_localpca_zblocks/07_L2/`                 |
| `07` plot_L2   | `path_localpca_zblocks/08_L2_plots/`           |
| `08a`          | `_shared/sample_metadata.tsv`                  |
| `08b`          | `path_localpca_zblocks/09_atlas_json/`         |

Path 2 and path 3, when refactored, will have NO `01a` equivalent (their
features come from upstream pipelines outside this toolkit), so they start
at `01_local_pca/` with no offset.

---

## 2. The data flow — copy/paste runnable on LANTA

All commands assume the canonical scratch tree above. The default config
file lives at `${SCRATCH}/inversion_localpca_v8/_shared/00_inversion_config.sh`.

### Step 01 — beagle.gz → dosage + sites

(Per-chrom; one task per chromosome line in `chrom.list`.)

Currently the launcher script `01a_LAUNCH_beagle_to_dosage.slurm` is
**unchanged from before** — it predates the consolidated layout. Run it as
you have been; the output should land in `02_dosage_sites/`.

### Step 01b/c — per-chrom dense local PCA (compute + merge)

```bash
# Stage 1: one chromosome per array task
sbatch --array=0-27 99_launchers/01b_LAUNCH_local_pca_compute.slurm chrom.list

# Stage 2: merge (assigns global window_ids)
sbatch --dependency=afterok:<JOB01b> 99_launchers/01c_LAUNCH_local_pca_merge.slurm
```

Produces in `path_localpca_zblocks/02_dense_registry/`:
- `windows_master.tsv.gz` — the master window registry
- `<chr>.window_pca.rds` — per-chrom local PCA, top-`npc=4` eigvecs + full scree spectrum
- `<chr>.window_pca.tsv.gz`

### Step 02a/b — lostruct distance + MDS (compute + merge)

```bash
# Stage 1: one focal chromosome per array task (~12h walltime)
sbatch --array=0-27 99_launchers/02a_LAUNCH_mds_compute.slurm chrom.list

# Stage 2: merge into final mds.rds with $per_chr structure
sbatch --dependency=afterok:<JOB02a> 99_launchers/02b_LAUNCH_mds_merge.slurm
```

Produces in `path_localpca_zblocks/03_mds/`:
- `inversion_localpca.mds.rds` (with `$per_chr` field every downstream consumes)
- `inversion_localpca.window_mds.tsv.gz`
- `inversion_localpca.candidate_regions.tsv.gz`

Default mode: `chunked_2x` — each focal chrom is MDS'd against itself plus
2× background sampled from non-focal chromosomes (excluding high
inv_likeness windows so foreign inversions don't leak into "background").

### Step 03 — precompute (per-chrom features + NN sim_mats)

```bash
sbatch --dependency=afterok:<JOB02b> 99_launchers/03_LAUNCH_precompute_localpca_zblocks.slurm
```

Produces in `path_localpca_zblocks/04_precomp/`:
- `precomp/<chr>.precomp.rds` — full precomp (with `PC_1_*` and `PC_2_*` per-sample columns)
- `precomp/sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds`
- `window_dt.tsv.gz`, `precomp_summary.tsv`

**Tune NN scales saved:**
```bash
NN_SIM_SCALES="40,80,160,320" sbatch 99_launchers/03_LAUNCH_precompute_localpca_zblocks.slurm
```

### Steps 04–08b — interactive per chromosome

The bulk-compute steps are 01–03. Steps 04–08b are tuned interactively
(parameter sweeps) and then rolled out per chromosome with canonical
defaults baked in. The driver:

```bash
bash 99_launchers/04to08b_run_one_chrom.sh C_gar_LG28
```

Or call individual scripts with the canonical-mode flags:

```bash
# 04 — L1 detect (defaults: --nn 80, boundary_W 5, boundary_offset 5,
#                         boundary_min_dist 30, validator_mode grow)
Rscript 04_detect_L1/04_detect_L1_localpca_zblocks.R \
  --precomp_dir <SCRATCH>/path_localpca_zblocks/04_precomp/precomp \
  --chr         C_gar_LG28 \
  --outdir      <SCRATCH>/path_localpca_zblocks/05_L1 \
  --boundary_scan TRUE \
  --boundary_validator_mode grow \
  --boundary_W 5 --boundary_offset 5 --boundary_min_dist 30
```

```bash
# 05 — L1 plot (default --nn 80, --boundary_filter stable)
Rscript 05_plot_L1/05_plot_L1_localpca_zblocks.R \
  --precomp_dir <SCRATCH>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir      <SCRATCH>/path_localpca_zblocks/05_L1 \
  --chr         C_gar_LG28 \
  --outdir      <SCRATCH>/path_localpca_zblocks/06_L1_plots \
  --toggle_L1 yes --boundary_filter stable
```

```bash
# 06 — L2 detect (default --nn 40, with quadrant validator)
Rscript 06_detect_L2/06_detect_L2_localpca_zblocks.R \
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
# 07 — L2 plot (default --nn 80 chrom-wide + --nn_l2 40 inside-segment)
Rscript 07_plot_L2/07_plot_L2_localpca_zblocks.R \
  --precomp_dir <SCRATCH>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir      <SCRATCH>/path_localpca_zblocks/05_L1 \
  --L2_dir      <SCRATCH>/path_localpca_zblocks/07_L2 \
  --chr         C_gar_LG28 \
  --outdir      <SCRATCH>/path_localpca_zblocks/08_L2_plots \
  --boundary_filter stable
```

```bash
# 08a — build sample metadata (run ONCE, genome-wide)
Rscript 08_atlas_json/08a_build_sample_metadata.R \
  --bamlist    <path>/list_of_samples_one_per_line_same_bamfile_list.tsv \
  --pairs      <path>/catfish_226_for_natora.txt \
  --theta_cutoff 0.177 \
  --ancestry   <path>/ngsadmix_K8_ancestry.tsv \
  --out        <SCRATCH>/_shared/sample_metadata.tsv
```

```bash
# 08b — export atlas JSON per chromosome
Rscript 08_atlas_json/08b_export_atlas_json_localpca_zblocks.R \
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
| 04 detect_L1 | 80         | `--nn N`        |
| 05 plot_L1   | 80         | `--nn N`        |
| 06 detect_L2 | 40         | `--nn N`        |
| 07 plot_L2   | 80 chrom-wide + 40 inside-segment | `--nn N` + `--nn_l2 N` |
| 08b atlas    | 40, 80, 160, 320 (multi-scale) | `--nn_list 40,80,160,320` |

L1 uses **nn80** because chromosome-wide it needs to suppress short-range
noise. L2 uses **nn40** because inside an L1 segment, finer resolution is
the goal. (History lines 19553 / 19556.)

### 3.2 Path resolution flags (canonical mode)

All four detect/plot scripts and 08b accept:

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
dropped to keep file ~10 MB). Steps 04–07 only need window coordinates so
slim worked for them, but step 08b NEEDS `PC_1_*` and `PC_2_*` per-sample
columns (lost in slim → falls back to PC2 jitter).

**Decision: full precomp only.** All scripts read the same
`<chr>.precomp.rds` (~70-100 MB, includes per-sample PCs). Steps 04–07
only touch window coords; the extra columns are harmless. Step 08b gets
real PC1+PC2.

The `prep_lg28_bundle.R` workaround that built slim copies is retired.

### 3.5 Identity reconciliation — separated and centralized

The three identity layers (bamlist remap, ngsRelate family graph, NGSadmix
ancestry) used to be mashed together inside the JSON exporter via three
separate flags (`--bamlist`, `--pairs`, `--samples`).

**New design: `08a_build_sample_metadata.R`** consumes the three independent
inputs and produces ONE merged TSV (`sample_metadata.tsv`, columns:
`ind, cga, family_id, ancestry`). The exporter `08b` then takes a single
`--sample_metadata` flag.

Benefits:
- Reconciliation logic lives in ONE place
- The merged TSV can be sanity-checked before paying JSON-build cost
- Same TSV reusable by paths 2 and 3 (all three discovery paths use the
  same 226 samples) — `_shared/sample_metadata.tsv` is the single source
  of truth

Legacy `--bamlist` + `--pairs` + `--samples` flags still work in 08b for
power users; if `--sample_metadata` is given it short-circuits the legacy
logic.

### 3.6 Stage1/stage2 → compute/merge

Stage1+stage2 was infrastructure (parallelize across chromosomes via SLURM
array, then merge). Renamed for clarity:

- `01b_local_pca_compute.R` (was `01c_local_pca_stage1.R`)
- `01c_local_pca_merge.R`   (was `01d_local_pca_stage2.R`)
- `02a_mds_compute.R`       (was `02a_mds_stage1.R`)
- `02b_mds_merge.R`         (was `02b_mds_stage2.R`)

### 3.7 Overlay → plot

`05_overlay_L1.R` → `05_plot_L1_localpca_zblocks.R`
`07_overlay_L2.R` → `07_plot_L2_localpca_zblocks.R`

Same naming convention as the detect scripts, easier to grep.

---

## 4. The three discovery paths — symmetry

All three paths share this skeleton:

```
01_local_pca_compute  (per-chr ARRAY, sliding-window PCA on path-specific feature)
02_local_pca_merge    (global window_ids)
03_mds                (lostruct + cmdscale, with chunked background)
04_precompute         (per-chrom features + NN sim_mats)
05_detect_L1          (chrom-wide z-block detection, nn80)
06_plot_L1            (multi-page PDF)
07_detect_L2          (per-L1-segment fine sub-block detection, nn40)
08_plot_L2            (multi-page PDF)
09_atlas_json         (consume sample_metadata.tsv from _shared/, emit JSON)
```

Path 1 (this folder) has an extra `01a_beagle_to_dosage.py` because dosage
is computed inside this toolkit. Paths 2 and 3 inherit their feature
matrices from upstream pipelines (ANGSD `-doThetas` for path 2; Clair3
phasing for path 3) and skip directly to `01_local_pca_compute`.

When paths 2 and 3 are refactored later, they'll match the same structure,
just with different feature inputs feeding the local PCA. The MDS,
precompute, detect, plot, and atlas-JSON code is fully shared in design (and
maybe even partly in implementation — the L1/L2 detect logic operates on
sim_mats regardless of what the feature was).

---

## 5. What's left for future sessions

### High priority

1. **NEXT SESSION — reorganize path 2 (θπ).** The user will upload
   `2f_theta_discovery/` for a focused session. Read
   `_to_be_reorganized/2f_theta_discovery/HANDOFF.md` first — it
   recommends Option B (lighter touch: rename + rewire scratch paths,
   don't unpack the everything-script). The user also has an open
   question on whether to lift `04_detect_L1` / `06_detect_L2` /
   `05_plot_L1` / `07_plot_L2` into a `_shared_detect_plot/` folder
   callable from all three paths — see README §14.1. Discuss before
   touching anything.

2. **`01a_LAUNCH_beagle_to_dosage.slurm`** — still uses old paths, predates
   the consolidated layout. Update its config sourcing to match the other
   launchers.

3. **Reorganize path 3 (`_to_be_reorganized/2e_ghsl_discovery/`)** the same
   way. Likely more involved because `STEP_C04_snake3_ghsl_v6.R` looks like
   a monolithic compute script that may combine multiple canonical steps.

### Medium priority

4. **Sample metadata script (08a) — test on real LANTA data.** The script
   was newly written for this consolidation; it parses cleanly but has not
   been run end-to-end on the real bamlist + ngsRelate output. Verify
   against an existing manually-built sample_metadata before trusting it.

5. **Stub the `_shared/00_inversion_config.sh` template.** A skeleton
   exists in the bundle (`_shared/00_inversion_config.sh.template`), but
   the real values (paths, env names, sample counts) need to be plugged in
   from your live config.

### Low priority / deferred

6. **Genome-wide rollout SLURM array for steps 04-08b.** Currently steps
   04-07 are interactive; the `04to08b_run_one_chrom.sh` driver runs them
   for one chrom at a time. A SLURM array launcher that runs the driver
   in parallel across 28 chromosomes would be useful for genome-wide
   atlas builds.

7. **The legacy `01a_beagle_to_dosage.py` script should be reviewed.**
   Output naming (`<chr>.dosage.tsv.gz` + `<chr>.sites.tsv.gz`) matches
   the new convention, so it should be fine, but worth a sanity pass.

---

## 6. The full DAG

```
                         (per chromosome)
beagle.gz ──► 01a beagle_to_dosage ──► dosage + sites
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ 01b local_pca_compute (ARRAY)        │
                       │ 01c local_pca_merge   (single)       │
                       └──────────────────────────────────────┘
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ 02a mds_compute (ARRAY focal-chr)    │
                       │ 02b mds_merge   (single)             │
                       └──────────────────────────────────────┘
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ 03 precompute_localpca_zblocks       │
                       │   precomp/<chr>.precomp.rds          │
                       │   sim_mats/<chr>.sim_mat_nn{0..320}  │
                       └──────────────────────────────────────┘
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ 04 detect_L1 (nn80)                  │
                       │   <chr>.L1_envelopes.tsv             │
                       │   <chr>.L1_boundaries.tsv            │
                       │   <chr>.L1_score_curve.tsv           │
                       └──────────────────────────────────────┘
                                │                 │
                                │                 ▼
                                │        ┌──────────────────┐
                                │        │ 05 plot_L1       │
                                │        │ <chr>.L1_overlay │
                                │        └──────────────────┘
                                ▼
                       ┌──────────────────────────────────────┐
                       │ 06 detect_L2 (nn40, per L1 segment)  │
                       │   <chr>.L2_envelopes.tsv             │
                       │   <chr>.L2_boundaries.tsv            │
                       │   <chr>.L2_segment_stats.tsv         │
                       │   <chr>.L2_quadrant_validator.tsv    │
                       └──────────────────────────────────────┘
                                │                 │
                                │                 ▼
                                │        ┌──────────────────┐
                                │        │ 07 plot_L2       │
                                │        │ <chr>.L2_overlay │
                                │        └──────────────────┘
                                ▼
   bamlist ┐                    │
   pairs   ├─► 08a build_sample_metadata ─► sample_metadata.tsv
   ancestry┘  (genome-wide ONCE)           (in _shared/)
                                          │
                                          ▼
                       ┌──────────────────────────────────────┐
                       │ 08b export_atlas_json                │
                       │   inputs: full precomp + L1+L2 +     │
                       │           sample_metadata.tsv +      │
                       │           sim_mats nn{40,80,160,320} │
                       │   output: <chr>.atlas.json           │
                       └──────────────────────────────────────┘
                                          │
                                          ▼
                                pca_scrubber_v3 / atlas page-1
```
