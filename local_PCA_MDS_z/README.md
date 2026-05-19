# `local_PCA_MDS_z/` — detailed data flow

Path 1 of catfish-inversion-analysis (the dosage-based local-PCA z-blocks
pipeline). This README walks through the complete pipeline end-to-end: every
step's purpose, what it reads, what it writes, what knobs matter, and how the
pieces connect on disk.

For the high-level summary and quick-start commands, see **HANDOFF.md**.
For per-step parameter notes, see **`docs/PER_STEP_NOTES.md`**.

**Layout note:** all step scripts and SLURM launchers live flat at the top
level of this folder, prefixed `STEP_ZO_<LETTER>_` (Z Outlier path, lettered
B–M in pipeline order — A is reserved for any future upstream cohort-prep
step). The shared config template lives one level up, at the repo root
(`../00_inversion_config.sh.template`).

---

## Table of contents

1. [The big picture](#1-the-big-picture)
1.5. [Method — high level](#15-method)
2. [The two trees: scripts and scratch](#2-the-two-trees-scripts-and-scratch)
3. [Step B — beagle.gz to dosage](#3-step-b)
4. [Steps C/D — local PCA (compute + merge)](#4-steps-cd)
5. [Steps E/F — MDS (compute + merge)](#5-steps-ef)
6. [Step G — precompute](#6-step-g)
7. [Step H — detect_L1](#7-step-h)
8. [Step I — plot_L1](#8-step-i)
9. [Step J — detect_L2](#9-step-j)
10. [Step K — plot_L2](#10-step-k)
11. [Step L — build sample metadata](#11-step-l)
12. [Step M — export atlas JSON](#12-step-m)
13. [Glossary of artifacts](#13-glossary-of-artifacts)
14. [Three-path symmetry](#14-three-path-symmetry)

---

## 1. The big picture

The pipeline turns raw genotype likelihoods into **per-chromosome JSON files**
that the `pca_scrubber_v3` browser tool consumes. The end goal is interactive
inspection of every chromosome's chromatin-block / inversion structure across
the 226-sample cohort.

```
beagle.gz                    (input — ANGSD genotype likelihoods)
   │
   ▼
dosage + sites               (step B)
   │
   ▼
sliding-window local PCA     (steps C/D — array per chromosome)
   │
   ▼
lostruct distance + MDS      (steps E/F — array per focal chromosome)
   │
   ▼
precomp + NN sim_mats        (step G — features for downstream)
   │
   ▼
detect L1 (chromosome-wide)  (step H — boundary scan, nn80)
   │
   ├──► plot L1              (step I — multi-page PDF)
   ▼
detect L2 (per L1 segment)   (step J — sub-block scan, nn40)
   │
   ├──► plot L2              (step K — multi-page PDF)
   ▼
sample metadata              (step L — bamlist + ngsRelate + ancestry, ONCE)
   │
   ▼
atlas JSON                   (step M — per chromosome)
   │
   ▼
pca_scrubber_v3              (browser viewer)
```

**Three flavours of step:**

* **Heavy bulk compute** (steps C–G) — submitted as SLURM array jobs, run
  once per cohort, take many hours. Re-run when you change input data.
* **Interactive parameter-tuned** (steps H–K) — run from the command line
  per chromosome, fast (~1 min per chrom). Re-run while sweeping parameters.
* **One-shot identity merge + per-chrom export** (steps L/M) — L runs
  once genome-wide; M runs per chromosome.

**Two levels of structure:**

* **L1** — chromosome-wide. Where do major architectural blocks
  (e.g. inversions, large LD blocks) start and end? Detected on **nn80**
  smoothed similarity matrix.
* **L2** — sub-divisions of each L1 block. Inside each L1 segment, where are
  the finer sub-block transitions? Detected on **nn40** (finer scale; smoother
  matrices wash out the substructure here).

---

## 1.5 Method — high level <a name="15-method"></a>

This section explains the *algorithm*, not the file plumbing. The key
question: what's a "window" and how does it carry information from raw
genotype data through to the atlas?

### 1.5.1 Windows are SNP-count-defined, not bp-defined

Step C slides over the chromosome **counting SNPs**, not base pairs:

```
--winsize 100      # window = 100 SNPs
--step 20          # advance by 20 SNPs between windows (heavy overlap)
```

So in dense regions (many SNPs per kb), windows are physically narrow;
in sparse regions, windows span more bp to keep the SNP count constant.
Each window's **physical extent** (`start_bp`, `end_bp`) is recorded but
is a *consequence* of the SNP grid, not an input to it.

Why SNP-count windows: PCA quality depends on how many informative sites
you put into the covariance, not on the bp span. Holding SNPs constant
holds the per-window PCA noise floor constant.

The 80% overlap (step=20 with winsize=100) means adjacent windows share
80 of their 100 SNPs. This isn't redundancy — it's how you get smooth,
position-resolved tracks instead of jagged step functions. The smoothness
is necessary for the boundary-scan detectors in steps H and J to work.

For LG28 (~100,000 SNPs after filtering): you get ~5,000 windows at
winsize=100, step=20.

### 1.5.2 What each window IS, after step C

Each window is **NOT a single value**. After PCA on its 100×226 dosage
submatrix, each window carries:

| Per-window quantity | Shape | What it captures |
|---|---|---|
| Top eigenvalues `lam_1..lam_k` | 4 numbers | how concentrated the PCA is in low dimensions |
| Per-sample PC scores `PC_1_<sample>..PC_4_<sample>` | 226 × 4 = 904 numbers | each sample's position in the window's PC space |
| Full eigenvalue spectrum | 226 numbers | for downstream scree plots |
| Window coordinates | start_bp, end_bp, center_bp, n_snps | physical position |

So the "fingerprint" of a window is **226 × 4 = 904 PC scores** — one
4-vector per sample. Inversions show up as windows where heterozygous
inversion carriers cluster separately from the two homozygous classes
along PC1 (or PC1+PC2). Outside inversions, the cluster structure is
random sampling noise.

### 1.5.3 MDS reduces window-fingerprints to a 2D map (steps E/F)

Steps E/F compute pairwise **lostruct distances** between window
fingerprints. Lostruct's distance is angular: it measures how much the
top-K eigenvalue-weighted subspaces of two windows differ in orientation.
Two windows from inside the same inversion have similar PC subspaces →
small distance. A window inside an inversion vs. one outside → large
distance.

Then `cmdscale` projects the distance matrix to MDS coordinates. After
step F, each window IS a single point — but in MDS space, not in
genotype space:

```
window i  →  4-vector PCA fingerprint  →  lostruct dist to all others  →  MDS coord (x_i, y_i, ...)
```

The chunked-2x background sampling in step E is why focal-chromosome
inversions stand out: by including 2× background windows from non-focal
chromosomes (excluding their high-`inv_likeness` regions), the MDS
embedding has a meaningful "normal" cloud to be an outlier *against*.
Without background, you have nothing to be an outlier from.

### 1.5.4 The similarity matrix is the heart of everything downstream

Step G builds the per-chromosome similarity matrix:

```
sim_mat[i, j] = (some kernel of) MDS_distance(window_i, window_j)   for i, j in chromosome
```

`sim_mat` is N × N where N is the number of windows on this chromosome
(~5,000 for LG28). Self-similarity ≈ 1 on the diagonal; off-diagonal
similarity decreases with MDS distance.

**This single matrix encodes the architecture of the chromosome.** Big
inversions show up as large red blocks on the diagonal (windows inside
the inversion are all similar to each other). Inversion boundaries show
up as sharp transitions — the "blue cross-blocks" that step H detects.

NN smoothing (`sim_mat_nn{k}`) replaces each window's MDS coordinates
with the mean of itself + its k nearest neighbours in MDS space, then
recomputes similarity. Bigger k = more smoothing. nn80 is the canonical
L1 scale (chromosome-wide), nn40 is the canonical L2 scale (inside one
L1 segment).

**A note on the NN scales — honestly, we don't fully know.** ZO_G writes
`sim_mat_nn{0,40,80,160,320}` by default. Our working interpretation:

- **nn0**: unsmoothed raw `sim_mat`. Useful for cross-checking that the
  smoothing isn't washing out real structure.
- **nn40**: small-scale; preserves sub-block detail (good for L2 inside
  one called L1 segment, e.g. nested inversions, small structural
  variants, megabase-scale things).
- **nn80**: chromosome-wide L1 default; the level where one inversion
  shows up as one block and short-range noise is suppressed.
- **nn160**: bigger blocks. Probably the right scale for very large
  inversions or large LD blocks that span more than ~10 Mb.
- **nn320**: very big inversions — chromosome-arm-scale rearrangements
  if your organism has them.
- And what we *think* the low end (nn20, the old default we dropped)
  was picking up was: noise of family / pedigree LD on the scale of a
  handful of contiguous windows. We removed it because nothing
  downstream consumed it, but if you're working with a strongly
  family-structured sample you might want to add it back.

These interpretations are educated guesses, not anything we've
rigorously calibrated against ground truth. Different organisms, different
sample structures, different sequencing depths, different window sizes
will probably shift what each `k` actually captures. If you're adapting
this pipeline to another dataset and your L1 / L2 calls look weird, **try
overriding the NN scales** with `NN_SIM_SCALES="20,40,80,160,320"` or
similar before assuming anything else is broken. We genuinely don't have
enough comparative data to tell you the right values up front.

### 1.5.5 Detection: cross-block scan on the diagonal

Steps H and J use the **same algorithm**, just applied at different
scopes (chromosome vs inside one L1 segment):

1. Slide a small W×W upper-triangle cross-block at offset G along the
   diagonal of sim_mat. At position i:
   ```
   cross_block = sim_mat[(i-W+1):i, (i+G+1):(i+G+W)]
   boundary_score[i] = -median(z_normalize(cross_block))
   ```
   Bluer cross-block (more separation between regions) → higher score.
2. 1D peak detection on `boundary_score[i]`.
3. Validate each peak by re-evaluating at multiple W sizes (`grow_max_z`).
   REAL boundaries stay blue at all W; FAKE boundaries (a small blue
   cross inside one inversion) drift to positive z as W grows.

This algorithm is **signal-agnostic**: it doesn't care if sim_mat came
from dosage PCA, θπ PCA, or GHSL phased haplotypes — it only needs a
matrix where self-similarity ≈ 1 and pairwise similarity decreases with
pattern dissimilarity. That's why the same detect_L1 / detect_L2 scripts
are intended to be reused across all three discovery paths.

### 1.5.6 What flows where — the per-window column genealogy

| Quantity | Born at | Lives in | Used by |
|---|---|---|---|
| `start_bp`, `end_bp`, `center_bp` | C (per window from sites file) | every per-window TSV/RDS downstream | every step |
| `n_snps`, `step_snps`, `window_snps` | C | precomp `$dt`, plus chr_meta | provenance |
| `PC_1_<sample>..PC_4_<sample>` | C | window_pca.rds → carried into precomp `$dt` | step M atlas JSON (per-window per-sample scatters) |
| `lam_1..lam_4`, scree | C | window_pca.rds → carried into precomp | atlas scree plots |
| MDS coords `MDS1..MDSk` | E/F | mds.rds `$per_chr` → precomp `$dt` | sim_mat construction |
| `MDS{i}_z`, `max_abs_z`, `max_z_axis` | G | precomp `$dt` | candidate region detection, atlas Z-track |
| `sim_mat` (N×N) | G | precomp `$sim_mat` + sim_mats/<chr>.sim_mat_nn{k}.rds | steps H, I, J, K, M |
| `inv_likeness`, band shape, morphology | G | precomp `$dt` | downstream scoring (phase 4) |
| L1 envelopes / boundaries | H | TSVs in 05_L1/ | steps I, J, K, M |
| L2 envelopes / boundaries | J | TSVs in 07_L2/ | steps K, M |

Every per-window quantity flows through `start_bp` / `end_bp` as the
common thread — that's what allows steps H+ to align everything by
window position regardless of what feature the path is using.

### 1.5.7 Why path 2 (θπ) and path 3 (GHSL) are simpler architecturally

Path 1 needs **MDS + sim_mat** because PCA eigenvectors are sign-ambiguous
— you can't directly compare PC1 of window i to PC1 of window j (they
might just be flipped). The MDS step uses lostruct's sign-invariant
angular distance to recover meaningful similarity, then sim_mat gives
you the matrix the cross-block detector needs.

Paths 2 and 3 use sign-stable per-window scalars (θπ for path 2;
phased-haplotype divergence for path 3). They CAN compute a sim_mat
optionally (path 3 does, path 2 doesn't by default), but they don't
*need* MDS — the per-window scalar is itself a comparable quantity, and
the |Z| profile of those scalars across the chromosome is the canonical
signal. That's why path 2's pipeline is much shorter (one heavy + one
light script vs. path 1's six bulk-compute steps).

---

## 2. The two trees: scripts and scratch

These are kept conceptually separate.

### 2.1 Script tree (this repo)

Flat layout. Every step script and its SLURM launcher live at the top level,
prefixed `STEP_ZO_<LETTER>_`. The shared config template lives one level up
at the repo root.

```
catfish-inversion-analysis/
├── 00_inversion_config.sh.template          <- shared across all 3 paths
├── local_PCA_MDS_z/                         <- this folder (path 1)
│   ├── STEP_ZO_B_beagle_to_dosage.py
│   ├── STEP_ZO_B_LAUNCH_beagle_to_dosage.slurm
│   ├── STEP_ZO_C_local_pca_compute.R       <- array per chrom
│   ├── STEP_ZO_C_LAUNCH_local_pca_compute.slurm
│   ├── STEP_ZO_D_local_pca_merge.R         <- single short merge
│   ├── STEP_ZO_D_LAUNCH_local_pca_merge.slurm
│   ├── STEP_ZO_E_mds_compute.R             <- array per focal chrom (heaviest)
│   ├── STEP_ZO_E_LAUNCH_mds_compute.slurm
│   ├── STEP_ZO_F_mds_merge.R               <- single short merge
│   ├── STEP_ZO_F_LAUNCH_mds_merge.slurm
│   ├── STEP_ZO_G_precompute.R              <- mclapply across chroms
│   ├── STEP_ZO_G_LAUNCH_precompute.slurm
│   ├── STEP_ZO_H_detect_L1.R               <- default --nn 80
│   ├── STEP_ZO_I_plot_L1.R                 <- default --nn 80
│   ├── STEP_ZO_J_detect_L2.R               <- default --nn 40
│   ├── STEP_ZO_K_plot_L2.R                 <- default --nn 80 + --nn_l2 40
│   ├── STEP_ZO_L_build_sample_metadata.R   <- run ONCE (genome-wide)
│   ├── STEP_ZO_M_export_atlas_json.R       <- per chromosome
│   ├── STEP_ZO_HtoM_run_one_chrom.sh       <- driver chains H→M for one chrom
│   ├── docs/PER_STEP_NOTES.md
│   ├── HANDOFF.md
│   └── README.md                           <- this file
├── local_PCA_MDS_thetapi/                   <- path 2 (sibling, STEP_TR_*)
└── local_PCA_MDS_GHSL/                      <- path 3 (sibling, STEP_GH_*)
```

`STEP_ZO_` = "Z Outlier" — the dosage z-blocks discovery path. Sibling
folders use `STEP_TR_*` (theta-pi) and `STEP_GH_*` (GHSL) when they're
refactored to the same convention. Letters start at **B** (not A) because
A is reserved for upstream cohort-prep steps that may live one level up
in the inversion-analysis repo.

### 2.2 Scratch tree (LANTA filesystem)

Lives at `${SCRATCH}/inversion_localpca_v8/`. Holds outputs only. The
**scratch-tree folder names keep the original 01–09 numbering** — these
are output buckets, not script names, and downstream consumers (atlas,
sibling paths, registry) hardcode them. Only the script names changed.

```
inversion_localpca_v8/
│
├── 01_beagle/                            <- ANGSD output (input to step B)
│   └── <chr>.beagle.gz
│
├── 02_dosage_sites/                      <- step B output (SHARED upstream)
│   ├── <chr>.dosage.tsv.gz                  feeds path 1
│   └── <chr>.sites.tsv.gz                   feeds path 1
│
├── 03_pestPG/                            <- ANGSD -doThetas (path 2 input — future)
├── 04_clair3_phased_GHSL/                <- Clair3 phased haplotypes (path 3 input)
│
├── path_localpca_zblocks/                <- path 1 outputs (THIS pipeline)
│   ├── 01_local_pca/                          (step C: tmp/<chr>.window_pca_tmp.rds)
│   ├── 02_dense_registry/                     (step D: <chr>.window_pca.rds + master)
│   ├── 03_mds/                                (steps E/F: inversion_localpca.mds.rds)
│   ├── 04_precomp/                            (step G: precomp/ + sim_mats/)
│   ├── 05_L1/                                 (step H: <chr>.L1_*.tsv)
│   ├── 06_L1_plots/                           (step I: <chr>.L1_overlay.pdf)
│   ├── 07_L2/                                 (step J: <chr>.L2_*.tsv)
│   ├── 08_L2_plots/                           (step K: <chr>.L2_overlay.pdf)
│   └── 09_atlas_json/                         (step M: <chr>.atlas.json)
│
├── path_localpca_thetapi/                <- path 2 outputs (when refactored)
├── path_localpca_GHSL/                   <- path 3 outputs (when refactored)
│
└── _shared/
    ├── sample_metadata.tsv               <- step L output, used by ALL paths
    ├── 00_inversion_config.sh
    └── reference/fClaHyb_Gar_LG.fa
```

**Why script-letter naming doesn't match scratch-folder numbering:** the
root of the scratch tree has 4 upstream-shared folders that eat the first
4 numbers (`01_beagle`, `02_dosage_sites`, `03_pestPG`,
`04_clair3_phased_GHSL`). Inside `path_localpca_zblocks/` the numbering
restarts at `01_local_pca/` so the path folder is self-coherent. The
script-to-scratch mapping is in HANDOFF §1.2 if you need it.

---

## 3. Step B — beagle.gz → dosage <a name="3-step-b"></a>

**Script:** `STEP_ZO_B_beagle_to_dosage.py`
**Type:** per-chromosome (one task per chr line)
**Run on:** SLURM (`STEP_ZO_B_LAUNCH_beagle_to_dosage.slurm`)

### What it does

Converts ANGSD's beagle.gz genotype likelihoods to a flat dosage matrix
(site × sample) plus a sites file (per-site allele/position metadata). Most
downstream steps work on dosage, not on raw likelihoods.

### Reads

* `${SCRATCH}/01_beagle/<chr>.beagle.gz` (one chr per file)

### Writes

* `${SCRATCH}/02_dosage_sites/<chr>.dosage.tsv.gz` — sites × samples dosage
* `${SCRATCH}/02_dosage_sites/<chr>.sites.tsv.gz` — per-site metadata

### Why this output is in the SHARED root, not inside path_localpca_zblocks/

Because dosage feeds many things, not just inversion discovery. ngsRelate,
NGSadmix, Fst stamping at phase 4, etc. all consume the same dosage files.
Putting them at the scratch root makes them shareable.

---

## 4. Steps C / D — local PCA (compute + merge) <a name="4-steps-cd"></a>

### 4.1 The split

Local PCA on dosage is parallelized as **stage1 (compute) + stage2 (merge)**:

* **C compute** — ONE SLURM array task per chromosome. Each task is fully
  independent: reads its own dosage, runs sliding-window PCA, writes per-chr
  RDS with `window_id = NA`. ~6 hr walltime per chrom; 28 in parallel ≈ 6 hr
  total.
* **D merge** — ONE short job (~seconds) after the array completes. Reads
  every per-chr RDS in `tmp/`, assigns globally unique sequential
  `window_id` values across the genome, rewrites each per-chr RDS with the
  patched IDs.

The merge exists *only* because each per-chr task can't know how many
windows the chromosomes before it produced. That's the only piece of state
that can't be parallelized.

### 4.2 Step C — compute

**Script:** `STEP_ZO_C_local_pca_compute.R`
**Run on:** `sbatch --array=0-27 STEP_ZO_C_LAUNCH_local_pca_compute.slurm chrom.list`

#### What it does

For one chromosome:
1. Loads `<chr>.dosage.tsv.gz` (sites × samples).
2. Slides a window of `winsize=100` SNPs along the chromosome with `step=20`.
3. Runs PCA on the (samples × 100) submatrix at each window. Stores top
   `npc=4` eigenvectors per window AND the full eigenvalue spectrum
   (length n_samples) for downstream scree plots.

#### Reads

* `${SCRATCH}/02_dosage_sites/<chr>.dosage.tsv.gz`
* `${SCRATCH}/02_dosage_sites/<chr>.sites.tsv.gz`

#### Writes (in `${OUTDIR}/tmp/`)

* `<chr>.window_pca_tmp.rds` — per-chr PCA + scree, **`window_id = NA`**
* `<chr>.chr_meta.tsv.gz` — per-window metadata, local indices
* `<chr>.pca_table.tsv.gz` — per-window PCA loadings, `window_id = NA`
* `<chr>.summary.tsv` — one-row summary

#### Knobs

| Flag | Default | Notes |
|---|---|---|
| `--winsize` | 100 | SNPs per window |
| `--step` | 20 | SNP step between windows |
| `--npc` | 4 | top eigenvectors stored per window |

### 4.3 Step D — merge

**Script:** `STEP_ZO_D_local_pca_merge.R`
**Run on:** `sbatch --dependency=afterok:<JOB_C> STEP_ZO_D_LAUNCH_local_pca_merge.slurm`

#### What it does

1. Reads all `<outdir>/tmp/<chr>.chr_meta.tsv.gz` in chromosome sort order.
2. Counts windows per chromosome → assigns globally unique
   `window_id = 1, 2, ..., total_n_windows`.
3. Patches each per-chr `.window_pca_tmp.rds` with the real window_id.
4. Renames to `<chr>.window_pca.rds` (drops the `_tmp` suffix).
5. Writes `windows_master.tsv.gz` — the master genome-wide window registry.

#### Reads

* `${OUTDIR}/tmp/*.chr_meta.tsv.gz` (28 files)
* `${OUTDIR}/tmp/*.window_pca_tmp.rds` (28 files)

#### Writes (in `${OUTDIR}/`, which is `path_localpca_zblocks/02_dense_registry/`)

* `windows_master.tsv.gz` — master window registry
* `windows_master_summary.tsv` — per-chrom counts
* `<chr>.window_pca.rds` — finalized per-chr PCA RDS
* `<chr>.window_pca.tsv.gz` — finalized per-chr PCA loadings table

---

## 5. Steps E / F — MDS (compute + merge) <a name="5-steps-ef"></a>

### 5.1 The expensive step

Step E is the heaviest in the pipeline. For each focal chromosome, it
computes pairwise lostruct distances between **every pair of windows**
(focal + background). With ~1800 focal windows and 2× background = ~5400
windows, that's a 5400 × 5400 distance computation, O(n²). One focal
chromosome takes ~12 hours. 28 in parallel ≈ 12 hours wall.

### 5.2 Step E — compute

**Script:** `STEP_ZO_E_mds_compute.R`
**Run on:** `sbatch --array=0-27 STEP_ZO_E_LAUNCH_mds_compute.slurm chrom.list`

#### What it does (per focal chromosome)

1. Loads ALL 28 `<chr>.window_pca.rds` from step D (yes, every chrom — needed
   for chunked background sampling).
2. **Background sampling** depends on `--mds_mode`:
   * `chromosome` — focal chrom only, no background.
   * `global` — all 28 chroms together.
   * `chunked_2x/3x/4x` — focal + N×focal-window-count from non-focal chroms,
     **excluding high-`inv_likeness` windows**. The exclusion is critical:
     without it, inversions on other chromosomes leak into "background" and
     mask the focal-chrom inversion as "normal".
3. Computes lostruct distance on focal+background windows (O(n²) loop).
4. Runs `cmdscale` MDS in the higher-dim space.
5. Extracts focal-chrom-only MDS coordinates.
6. Writes per-focal-chr RDS to `tmp/`.

#### Reads

* `${RDS_DIR}/<chr>.window_pca.rds` for ALL 28 chroms

#### Writes

* `${OUTDIR}/tmp/<focal_chr>.mds_perchr.rds` — focal coords + per-axis z + bg metadata

#### Knobs

| Flag | Default | Notes |
|---|---|---|
| `--mds_mode` | `chunked_2x` | Canonical mode. 2× bg from other chroms. |
| `--npc` | 4 | Must match step C's `--npc` |
| `--mds_dims` | 20 | MDS dimensions retained |
| `--z_thresh` | 3.0 | Candidate region z threshold |

### 5.3 Step F — merge

**Script:** `STEP_ZO_F_mds_merge.R`
**Run on:** `sbatch --dependency=afterok:<JOB_E> STEP_ZO_F_LAUNCH_mds_merge.slurm`

#### What it does

1. Reads 28 per-focal-chr `mds_perchr.rds` from `tmp/`.
2. Assembles them into the unified `inversion_localpca.mds.rds` with the
   `$per_chr` structure every downstream consumer expects.
3. Detects candidate regions (contiguous z>3 windows, gap-merged at 500 kb).
4. Writes the various TSV side-files for inspection.

#### Reads

* `${TMPDIR_IN}/<chr>.mds_perchr.rds` (28 files)

#### Writes (at `${OUTPREFIX}.*`, where outprefix = `path_localpca_zblocks/03_mds/inversion_localpca`)

* `.mds.rds` — main object (`$dt`, `$per_chr`, `$candidate_regions`, `$mds_mode`)
* `.window_mds.tsv.gz` — per-window MDS coordinates
* `.candidate_regions.tsv.gz` — coarse z-outlier regions
* `.candidate_window_membership.tsv.gz` — region ↔ window mapping
* `.mds_mode_metadata.tsv` — run config
* `.mds_background_<chr>.txt` — per-focal-chr background indices (chunked modes)

---

## 6. Step G — precompute <a name="6-step-g"></a>

**Script:** `STEP_ZO_G_precompute.R`
**Run on:** `sbatch --dependency=afterok:<JOB_F> STEP_ZO_G_LAUNCH_precompute.slurm`

### What it does

This is the central "feature factory" of path 1. It takes the unified MDS
object from step F and emits per-chromosome:
1. A rich `precomp.rds` with dozens of per-window features.
2. A suite of NN-smoothed similarity matrices at multiple scales.

Parallelized across chromosomes via `mclapply` (28 cores, 1:1 with chromosomes).

### What's in the precomp `$dt` data.table

Every column is per-window. The interesting groups:

| Group | Columns |
|---|---|
| Position | `global_window_id`, `chrom`, `start_bp`, `end_bp`, `center_bp`, `n_snps` |
| MDS | `MDS1..MDSk`, `MDS1_z..MDSk_z`, `max_abs_z`, `max_z_axis` |
| inv_likeness | `inv_likeness`, `inv_dip_p`, `inv_het_contrast` |
| Band shape | `band_discreteness`, `diffuse_score`, `het_intermediacy`, `n_effective_clusters` |
| Het sources | `pc1_bimodality`, `het_pc1_gap`, `het_mid_fraction`, `het_mid_variance`, `dosage_het_rate_*` |
| Local k=3 | `local_delta12`, `local_entropy`, `local_ena`, `band1/2/3_frac` |
| Adaptive | `beta_pval`, `adaptive_seed`, `beta_alpha`, `beta_beta` |
| Seed | `seed_nn_dist` |
| Morphology | `flat_inv_score`, `spiky_inv_score`, `fragmentation_score`, plus jaggedness/run/peak/plateau/nbhood/block columns |
| **Per sample** | `PC_1_Ind*`, `PC_2_Ind*` ← **REQUIRED by step M** |

The per-sample `PC_*_Ind*` columns are why the slim hack was killed. Step
M plots per-window per-sample PC1/PC2 in the scrubber; without these
columns it silently falls back to PC2 jitter.

### NN-smoothed similarity matrices — the key abstraction

The raw similarity matrix `sim_mat` (window × window) is noisy. Smoothing
each window by averaging with its k nearest neighbours **in MDS space**
produces a denoised matrix at scale k. The script writes nine scales by
default:

| Scale | Use case |
|---|---|
| `nn0` | raw, no smoothing |
| `nn20` | very fine resolution |
| `nn40` | **L2 multipass** (per-segment finer resolution) |
| `nn80` | **L1 multipass** + chrom-wide overlay |
| `nn120` | inter-scale checks |
| `nn160` | atlas thumbnail |
| `nn200` | inter-scale checks |
| `nn240` | inter-scale checks |
| `nn320` | atlas thumbnail (coarsest) |

Tune the saved set with `NN_SIM_SCALES="40,80,160,320" sbatch ...`.

### Reads

* `${MDS_PREFIX}.mds.rds` — output of step F
* `${DOSAGE_DIR}/<chr>.dosage.tsv.gz` (optional, enables `dosage_het_rate_cv`)

### Writes (in `${OUTDIR}/`, which is `path_localpca_zblocks/04_precomp/`)

* `precomp/<chr>.precomp.rds` — full per-window features (with PC_1_*, PC_2_*)
* `precomp/sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds`
* `window_dt.tsv.gz` — genome-wide per-window scalar table
* `precomp_summary.tsv` — per-chrom window counts + timing

### Sanity check after step G

```bash
Rscript -e 'x <- readRDS("<scratch>/path_localpca_zblocks/04_precomp/precomp/C_gar_LG28.precomp.rds")
            want <- c("max_abs_z", "max_z_axis", "inv_likeness")
            missing <- setdiff(want, names(x$dt))
            pc_cols <- length(grep("^PC_1_", names(x$dt)))
            cat("missing dims:", length(missing),
                "  PC_1_* cols:", pc_cols, "\n")'
```

Expected: `missing dims: 0  PC_1_* cols: 226`.

---

## 7. Step H — detect_L1 <a name="7-step-h"></a>

**Script:** `STEP_ZO_H_detect_L1.R`
**Run on:** interactive, per chromosome (or via `STEP_ZO_HtoM_run_one_chrom.sh`).

### What it does — the boundary-scan algorithm

For each chromosome:
1. Loads `<chr>.precomp.rds` (window coords) and `<chr>.sim_mat_nn80.rds`.
2. **Boundary score curve.** Slides a small W × W upper-triangle cross-block
   at offset G along the diagonal of the sim_mat. At each window position `i`:
   ```
   cross = Z[(i-W+1):i, (i+G+1):(i+G+W)]
   boundary_score[i] = -median(cross)
   ```
   Bluer cross (more separation between regions) ⇒ higher score.
3. **1D peak detection** on `boundary_score[i]` with `--boundary_min_dist`
   neighbourhood and `--boundary_score_min` floor.
4. **Validation.** Each peak gets two independent tests:
   * **(a) Growing-W cross-block z** (default validator). Re-evaluates the
     cross-block at multiple W sizes (`5, 10, 15, 20, ...`) and tracks
     `grow_max_z = max over W of median(z)`. REAL boundaries stay blue at all
     W; FAKE boundaries (small blue cross inside one inversion) drift toward
     positive z as the cross-block reaches into surrounding red bulk.
   * **(b) Perpendicular rays** (diagnostic). 1D z values along
     `Z[i, i+d]` and `Z[i-d, i]`.
5. **Adaptive thresholding.** `real_max = min(q_real(grow_max_z), real_max_ceiling)`
   and `fake_min = max(q_fake(grow_max_z), fake_min_floor)`. Falls back to
   fixed thresholds when too few peaks.
6. **Partition.** STABLE_BLUE peaks become L1 segment cut points. Single-window
   segments (caused by adjacent-tie peaks) merge into the more-similar
   neighbour.

### Default NN scale: 80

Settled in earlier sessions (history line 19553). Chromosome-wide L1
detection needs to suppress short-range noise; nn80 is the right level.

### Inputs (canonical mode)

```bash
Rscript STEP_ZO_H_detect_L1.R \
  --precomp_dir <scratch>/path_localpca_zblocks/04_precomp/precomp \
  --chr         C_gar_LG28 \
  --outdir      <scratch>/path_localpca_zblocks/05_L1
```

The script auto-resolves:
* `--precomp` ← `<precomp_dir>/<chr>.precomp.rds`
* `--sim_mat` ← `<precomp_dir>/sim_mats/<chr>.sim_mat_nn80.rds` (default `--nn 80`)

To override (e.g. parameter sweep at nn120):
```bash
... --nn 120
```

Or the legacy explicit-paths mode:
```bash
... --precomp <path>.precomp.rds --sim_mat <path>.sim_mat_nn120.rds
```

### Canonical parameters (history line 19553)

```
--boundary_scan TRUE
--boundary_validator_mode grow
--boundary_W 5 --boundary_offset 5 --boundary_min_dist 30
```

`--boundary_min_dist 30` was settled after sweeping 10/20/30 (history
19488 → 19510). 30 gives the cleanest STABLE_BLUE distribution.

### Outputs (in `--outdir`)

| File | Contents |
|---|---|
| `<chr>.L1_envelopes.tsv` | L1 segments, one row per segment, status='ENVELOPE' |
| `<chr>.L1_boundaries.tsv` | every detected peak with validation status |
| `<chr>.L1_score_curve.tsv` | per-window boundary score curve |

Boundary statuses:
* `STABLE_BLUE` — real long-range separation, persists at all val_W
* `STABLE_RED` — rejected by grow validator
* `WEAK_BLUE`, `DECAYS`, `MARGINAL`, `EDGE` — various failure modes

The downstream filter `--boundary_filter stable` (used in plotting) shows
only `STABLE_BLUE`.

---

## 8. Step I — plot_L1 <a name="8-step-i"></a>

**Script:** `STEP_ZO_I_plot_L1.R`
**Run on:** interactive, per chromosome.

### What it does

Multi-page L1 overlay PDF, one chromosome:
* **Page 1.** Whole-chromosome heatmap with L1 envelopes outlined and
  STABLE_BLUE boundaries marked. Uses **GLOBAL** per-distance Z (computed
  once over the whole chromosome).
* **Pages 2..(N+1).** One zoomed page per L1 segment. Each shows that
  segment's window range with the same envelope outlines and boundary
  squares. Uses **LOCAL** per-distance Z (computed from the segment's own
  sub-matrix), so colours normalize to the segment's distribution and
  short-range contrast is visible without being washed out.

Heatmap encoding:
* **lower triangle** = similarity (5-color gradient, `#F8F8F8` → `#7E1F1F`)
* **upper triangle** = Z (diverging, `#2C5AA0` → `#FAFAFA` → `#B22222`)

### Default NN scale: 80 (matches step H)

### Inputs (canonical mode)

```bash
Rscript STEP_ZO_I_plot_L1.R \
  --precomp_dir <scratch>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir      <scratch>/path_localpca_zblocks/05_L1 \
  --chr         C_gar_LG28 \
  --outdir      <scratch>/path_localpca_zblocks/06_L1_plots \
  --toggle_L1 yes \
  --boundary_filter stable
```

### Outputs

* `<chr>.L1_overlay.pdf` — single multi-page PDF

---

## 9. Step J — detect_L2 <a name="9-step-j"></a>

**Script:** `STEP_ZO_J_detect_L2.R`
**Run on:** interactive, per chromosome.

### What it does

Reads the L1 envelope catalogue from step H, then runs the same
boundary-scan logic **independently inside each L1 segment**, producing a
finer L2 partition.

### Two key v8 features

#### (a) Ward-adaptive thresholding

For each L1 segment:
1. Computes Ward.D2 clustering on the segment's windows.
2. Computes intensity statistics (`mean_sim`, `blue_red_ratio`).
3. Derives two priors (chromosome-rank-normalized):
   * `internal_homogeneity` — high ⇒ looks like inside one inversion
   * `architectural_complexity` — high ⇒ real architectural transitions inside
4. `strictness = internal_homogeneity - architectural_complexity` ∈ [-1, +1].
5. Each segment gets its **own** `(real_max, fake_min)` thresholds by scaling
   the pooled baseline:
   * STRICT segments tighten (most peaks reject as FAKE)
   * LOOSE segments widen (more peaks pass as REAL)

Self-calibrating per chromosome — no hardcoded class table or magic
`mean_sim` cutoff. Pass `--ward_adaptive no` to fall back to pooled-only
thresholds.

#### (b) Quadrant validator

A 2D test. A candidate L2 boundary divides the sim_mat into four
quadrants. The validator checks whether those four quadrants behave
like a true sub-block:
* The two "inside" quadrants (above and below the diagonal at the
  candidate position) should be high-sim within themselves.
* The two "cross" quadrants should be low-sim.

Demote-on-fail rescues weak peaks via `grow_z`. Enable with
`--quadrant_validator yes`.

### Default NN scale: 40 (finer than L1's 80)

Inside an L1 segment we want finer-resolution sub-structure; smoother
sim_mats wash that out. The chromosome-wide L1 pass uses nn80 because it
needs to suppress short-range noise; L2 inside a segment wants the opposite.

### Inputs (canonical mode)

```bash
Rscript STEP_ZO_J_detect_L2.R \
  --precomp_dir <scratch>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir      <scratch>/path_localpca_zblocks/05_L1 \
  --chr         C_gar_LG28 \
  --outdir      <scratch>/path_localpca_zblocks/07_L2 \
  --boundary_scan TRUE \
  --boundary_validator_mode grow \
  --quadrant_validator yes \
  --weak_demote_score 0 \
  --quad_rescue_max_grow_z 1.5 \
  --quad_demote_on_fail yes \
  --quad_demote_drift_floor -1.0
```

These canonical parameters come from history line 19585.

### Outputs

| File | Contents |
|---|---|
| `<chr>.L2_envelopes.tsv` | L2 sub-segments |
| `<chr>.L2_boundaries.tsv` | all detected peaks with validation status |
| `<chr>.L2_segment_stats.tsv` | per-L1-segment Ward stats + strictness |
| `<chr>.L2_quadrant_validator.tsv` | per-peak quadrant test (when enabled) |
| `<chr>.L2_quadrant_audit.tsv` | extra audit info |

---

## 10. Step K — plot_L2 <a name="10-step-k"></a>

**Script:** `STEP_ZO_K_plot_L2.R`
**Run on:** interactive, per chromosome.

### What it does

Multi-page L1+L2 overlay PDF:
* **Page 1.** Whole chromosome (nn80 by default), L1 envelopes + L1
  boundaries. Same as step I page 1.
* **Pages 2..(N+1).** One zoomed page per L1 segment, drawn on **nn40** by
  default (controlled by `--sim_mat_l2`). Shows the segment's L1 outline
  plus all L2 boundaries inside, as red squares.

### Default NN scales: 80 chrom-wide + 40 inside

Override with `--nn` (chrom-wide) and `--nn_l2` (inside-segment) for sweeps.

### Inputs (canonical mode)

```bash
Rscript STEP_ZO_K_plot_L2.R \
  --precomp_dir <scratch>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir      <scratch>/path_localpca_zblocks/05_L1 \
  --L2_dir      <scratch>/path_localpca_zblocks/07_L2 \
  --chr         C_gar_LG28 \
  --outdir      <scratch>/path_localpca_zblocks/08_L2_plots \
  --boundary_filter stable
```

### Outputs

* `<chr>.L2_overlay.pdf` — multi-page PDF

---

## 11. Step L — build sample metadata <a name="11-step-l"></a>

**Script:** `STEP_ZO_L_build_sample_metadata.R`
**Run on:** ONCE, genome-wide (not per chromosome).

### Why this script exists

The old JSON exporter took three flags (`--bamlist`, `--pairs`, `--samples`)
and reconciled them inline. That meant the reconciliation logic ran 28 times
(once per chrom JSON) and lived buried inside the exporter. Bad.

The new design: ONE script merges the three identity layers into ONE TSV,
written to `_shared/sample_metadata.tsv`. The exporter takes a single
`--sample_metadata` flag.

Benefits:
* Reconciliation logic in one place.
* Sanity-check the merged TSV before paying the JSON-build cost.
* Same TSV reusable across paths 2 and 3 (all three discovery paths use
  the same 226-sample cohort).

### The three independent identity layers

1. **bamlist** — one CGA-ID per line, in PC_1_Ind* order. Maps
   `Ind0..IndN → CGA####`.
2. **ngsRelate pairs** — 3-col TSV `id1 id2 theta`. Union-find with cutoff
   (default 0.177, Manichaikul 1st-degree) builds the family graph;
   connected components get sequential `family_id`; isolates get -1.
3. **ancestry** — TSV with `cga<TAB>ancestry_label`. Default "unknown" if
   not provided (e.g. NGSadmix K=8 cluster assignment).

### Reads

* `--bamlist <path>` — one CGA per line
* `--pairs <path>` — ngsRelate output (3-col)
* `--ancestry <path>` — optional ancestry assignment (`cga, ancestry`)

### Writes

* `--out <path>` — one TSV with columns `ind, cga, family_id, ancestry`,
  rows aligned to PC_1_Ind* order. Canonical location:
  `${SCRATCH}/_shared/sample_metadata.tsv`.

### Run once

```bash
Rscript STEP_ZO_L_build_sample_metadata.R \
  --bamlist     <path>/list_of_samples.tsv \
  --pairs       <path>/catfish_226_for_natora.txt \
  --theta_cutoff 0.177 \
  --ancestry    <path>/ngsadmix_K8_ancestry.tsv \
  --out         <scratch>/_shared/sample_metadata.tsv
```

---

## 12. Step M — export atlas JSON <a name="12-step-m"></a>

**Script:** `STEP_ZO_M_export_atlas_json.R`
**Run on:** per chromosome, after step L has run once.

### What it does

For one chromosome, packs everything the `pca_scrubber_v3` browser tool
needs into a single JSON:

* multi-scale sim_mat thumbnails (nn40, nn80, nn160, nn320 by default)
* per-distance local Z scores (matches the PDF colorimetry)
* L1 envelopes + L1 boundaries
* L2 envelopes + L2 boundaries
* `samples[i]` with `cga`, `family_id`, `ancestry` (from sample_metadata.tsv)
* `default_sim_scale` (which scale to render initially)
* `family_source` ("sample_metadata" in canonical mode)

### Inputs (canonical mode)

```bash
Rscript STEP_ZO_M_export_atlas_json.R \
  --precomp_dir     <scratch>/path_localpca_zblocks/04_precomp/precomp \
  --L1_dir          <scratch>/path_localpca_zblocks/05_L1 \
  --L2_dir          <scratch>/path_localpca_zblocks/07_L2 \
  --chr             C_gar_LG28 \
  --sample_metadata <scratch>/_shared/sample_metadata.tsv \
  --out             <scratch>/path_localpca_zblocks/09_atlas_json/C_gar_LG28.atlas.json
```

The script auto-resolves precomp + sim_mats + L1/L2 TSVs from
`--precomp_dir` + `--L1_dir` + `--L2_dir` + `--chr`.

### Optional knobs

| Flag | Default | Notes |
|---|---|---|
| `--nn_list` | `40,80,160,320` | comma-separated NN scales to embed |
| `--thumb_n` | 200 | thumbnail downsample target |
| `--sim_q_lo` / `--sim_q_hi` | 0.05 / 0.95 | sim quantile clip |
| `--z_clip` | 5 | Z range clip for colormap |
| `--track_<label> <tsv>` | — | add per-window scalar track to scrubber |

### Outputs

* `<chr>.atlas.json` — drag-and-drop into `pca_scrubber_v3`

### Why the JSON is at the END

An earlier session created a JSON exporter that ran right after the
precompute, emitting a JSON without L1/L2 results. **Wrong placement.**
The atlas needs the multipass results to be useful, so the JSON has to be
built AFTER step K. That earlier exporter is RETIRED.

---

## 13. Glossary of artifacts <a name="13-glossary-of-artifacts"></a>

| Artifact | Producer | Consumer |
|---|---|---|
| `<chr>.beagle.gz` | ANGSD | step B |
| `<chr>.dosage.tsv.gz` + `<chr>.sites.tsv.gz` | step B | step C |
| `<chr>.window_pca.rds` | step D | step E |
| `windows_master.tsv.gz` | step D | (registry) |
| `inversion_localpca.mds.rds` | step F | step G |
| `<chr>.precomp.rds` | step G | steps H, I, J, K, M |
| `<chr>.sim_mat_nn{0..320}.rds` | step G | steps H, I, J, K, M |
| `<chr>.L1_envelopes.tsv` | step H | steps I, J, K, M |
| `<chr>.L1_boundaries.tsv` | step H | steps I, K, M |
| `<chr>.L1_score_curve.tsv` | step H | (debug) |
| `<chr>.L1_overlay.pdf` | step I | (visual review) |
| `<chr>.L2_envelopes.tsv` | step J | steps K, M |
| `<chr>.L2_boundaries.tsv` | step J | steps K, M |
| `<chr>.L2_segment_stats.tsv` | step J | (debug / Ward stats inspection) |
| `<chr>.L2_quadrant_validator.tsv` | step J | (debug) |
| `<chr>.L2_overlay.pdf` | step K | (visual review) |
| `sample_metadata.tsv` | step L | step M (× 28 chroms × 3 paths) |
| `<chr>.atlas.json` | step M | `pca_scrubber_v3` |

---

## 14. Three-path symmetry <a name="14-three-path-symmetry"></a>

This pipeline is path 1 of three sibling discovery paths. All three share
the same skeleton; they differ only in the upstream feature matrix.

| Path | Feature input | Sibling folder | Step prefix | Status |
|---|---|---|---|---|
| `localpca_zblocks` | dosage from beagle.gz | `local_PCA_MDS_z/` | `STEP_ZO_` | refactored ← THIS |
| `localpca_thetapi` | per-window θπ from ANGSD pestPG | `local_PCA_MDS_thetapi/` | `STEP_TR_` | partial — v5 drafts in `v5_drafts/` |
| `localpca_GHSL` | small phased haplotypes from Clair3 | `local_PCA_MDS_GHSL/` | `STEP_GH_` | working, awaiting flatten |

Each path:
* Has its own folder at the repo root, flat layout, prefixed with its own
  `STEP_<XX>_` series, lettered B onwards.
* Skips step B if its feature matrix comes from upstream (paths 2 and 3 do —
  ANGSD `-doThetas` for path 2, Clair3 phasing for path 3).
* Writes to its own `path_localpca_<feature>/` scratch tree on LANTA.
* Consumes the SAME `_shared/sample_metadata.tsv` from step L (single
  source of truth for the 226-sample identity layer).

The atlas viewer reconciles per-chromosome pages across the three paths so
each chromosome's view is a unified multi-path display.

### 14.1 Open design question — shared detect/plot scripts across paths

The boundary-scan algorithm in `STEP_ZO_H_detect_L1.R` and
`STEP_ZO_J_detect_L2.R` is **signal-agnostic** (cf. §1.5.5 above). It only
needs a square sim_mat where self-similarity ≈ 1 and pairwise similarity
decreases with pattern dissimilarity — which describes path 1's MDS-derived
sim_mat AND path 3's `|cor(pc1[i], pc1[j])|` sim_mat AND any sim_mat path 2
might compute.

Path 2's v5 drafts at `local_PCA_MDS_thetapi/v5_drafts/` already use this:
they copy `STEP_ZO_H` and `STEP_ZO_J` verbatim as `STEP_TR_C_detect_L1.R`
and `STEP_TR_D_detect_L2.R`, run them on a theta-pi-shaped precomp, and
the algorithm just works.

**The clean version of this** is to lift the 4 scripts (detect_L1,
detect_L2, plot_L1, plot_L2) into a `_shared_detect_plot/` folder at the
repo root, and have each per-path launcher pass the right canonical
parameters (`--nn`, `--boundary_W`, `--boundary_min_dist` etc.). The v5
drafts decided to copy verbatim instead because the dosage canonical
defaults aren't necessarily right for theta-pi yet.

The complication: **canonical parameters differ per path.**
- Path 1 dosage settled at `--nn 80` for L1, `--nn 40` for L2,
  `--boundary_W 5`, `--boundary_min_dist 30 windows` (history line 19553).
- Path 2 θπ has not been swept — windows are defined by ANGSD's bp grid
  (~16,500 windows on LG28 at win10000.step2000), not by SNP count, so
  the right `--nn` will be very different.
- Path 3 GHSL has not been swept either — 5-kb base scale, ~4,300 windows
  on LG28, no calibration of D17 parameters yet.

The proposal (deferred until path 2 is refactored): **one shared script
per concern, no path-bias in defaults, per-path launchers fix the
canonical values.** In particular `--nn` would lose its hard default and
become required, forcing each launcher to declare which scale it uses.
That makes the shared script genuinely path-neutral and the per-path
tuning explicit at the launcher level — the right place for it.

This is the structural change to discuss when you reorganize path 2.
