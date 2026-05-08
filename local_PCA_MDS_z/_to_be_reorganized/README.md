# `_to_be_reorganized/` — path 2 + path 3 awaiting refactor

This folder holds **working code** for the other two discovery paths. They
have not yet been refactored to follow the consolidated 01–09 skeleton
that path 1 (`local_PCA_z/`) now uses.

## What's here

```
_to_be_reorganized/
├── 2e_ghsl_discovery/           ← path 3: GHSL haplotype contrast (Clair3 phased)
│   ├── STEP_C04_snake3_ghsl_v6.R
│   ├── STEP_C04b_snake3_ghsl_classify.R
│   ├── STEP_C04c_ghsl_local_pca.R
│   ├── STEP_C04d_ghsl_d17_wrapper.R
│   ├── export_ghsl_to_json_v3.R
│   ├── LAUNCH_STEP_C04_ghsl_v6_compute.slurm
│   ├── LAUNCH_STEP_C04b_ghsl_v6_classify.slurm
│   ├── LAUNCH_STEP_C04cd_ghsl_enrichment.slurm
│   └── README.md
│
└── 2f_theta_discovery/          ← path 2: θπ diversity (ANGSD pestPG)
    ├── STEP_TR_A_compute_theta_matrices.R
    ├── STEP_TR_B_classify_theta.R
    ├── STEP_TR_C_theta_d17_wrapper.R
    ├── STEP_TR_D_augment_theta_json.R
    ├── 00_sanity_check_pestPG_scaling.sh
    ├── 00_theta_config.sh
    ├── LAUNCH_TR_theta_pi.slurm
    ├── chrom.list
    ├── verify_theta_setup.sh
    └── README_theta_pi_scaling.md
```

## Why they're parked here

These are working pipelines that have been validated end-to-end on LANTA.
Refactoring them to match the path-1 skeleton requires:

1. Reading every script to understand which canonical step (01–09) it
   corresponds to. Some scripts may need to be **split** (one current
   script doing both compute and merge) or **merged** (multiple scripts
   that together implement one canonical step).
2. Rewriting launchers to write into the new scratch tree
   (`path_localpca_thetapi/01_local_pca/` etc.).
3. Stripping cohort-specific naming conventions (`STEP_TR_*`, `STEP_C04*`,
   `D17` etc.) in favor of the canonical artifact names
   (`<chr>.L1_envelopes.tsv` etc.).

**Read each folder's `HANDOFF.md` for the full architectural mapping**:

- **`2e_ghsl_discovery/HANDOFF.md`** — path 3 (GHSL). Has 5 scripts; key
  finding is that `STEP_C04c_ghsl_local_pca.R` rolls 3 canonical steps
  (01b + 02a + 03) into one, and `STEP_C04b_snake3_ghsl_classify.R` is
  GHSL-specific with no path-1 analogue.

- **`2f_theta_discovery/HANDOFF.md`** — path 2 (θπ). Has 4 scripts; key
  finding is that `STEP_TR_B_classify_theta.R` deliberately compresses 8
  canonical steps into one, because θπ is sign-stable and doesn't need
  MDS or sim_mat for primary detection.

Both HANDOFFs recommend **Option B** (lighter touch — rename + rewire
paths, don't unpack monolithic scripts) over Option A (full refactor to
11-script canonical skeleton). That's an entire session per path, with
you walking me through each script's role. **Not something to do blind.**

## Target structure (when refactored)

Each path becomes a sibling folder of `local_PCA_z/`:

```
inversion-popgen-toolkit/
├── local_PCA_z/                ← path 1 (DONE, this consolidation)
├── local_PCA_thetapi/          ← path 2 (TO BUILD from 2f_theta_discovery)
└── local_PCA_GHSL/             ← path 3 (TO BUILD from 2e_ghsl_discovery)
```

With the same 01–09 skeleton:

```
local_PCA_<feature>/
├── 01_local_pca_compute.R     ← reads path-specific feature, runs sliding-window PCA
├── 02_local_pca_merge.R       ← global window IDs
├── 03_mds_compute.R + merge.R ← lostruct + cmdscale
├── 04_precompute_*.R          ← per-chrom features + NN sim_mats
├── 05_detect_L1_*.R           ← default --nn 80
├── 06_plot_L1_*.R
├── 07_detect_L2_*.R           ← default --nn 40
├── 08_plot_L2_*.R
├── 09_export_atlas_json_*.R   ← consumes _shared/sample_metadata.tsv
├── 99_launchers/
└── 99_docs/
```

Note: paths 2 and 3 don't need a `01a_beagle_to_dosage.py` equivalent —
their feature matrices come from upstream pipelines (ANGSD `-doThetas`
producing `<chr>.pestPG`, Clair3 phasing producing phased haplotypes)
that live OUTSIDE this toolkit. They start at `01_local_pca_compute.R`
reading directly from `<SCRATCH>/03_pestPG/` or
`<SCRATCH>/04_clair3_phased_GHSL/`.

## Order of attack (suggested)

1. **Path 2 (θπ) first.** Smaller, simpler. The 4 scripts (TR_A → TR_B →
   TR_C → TR_D) probably map ~1:1 to canonical steps with minimal splitting.
2. **Path 3 (GHSL) second.** More involved. `STEP_C04_snake3_ghsl_v6.R`
   looks like a monolithic compute that may combine local PCA + MDS +
   precompute. Will need careful unpacking.

## Inputs in the scratch tree

When refactored, these paths will read from:

| Path                 | Reads from                               |
|----------------------|------------------------------------------|
| `local_PCA_thetapi`  | `<SCRATCH>/03_pestPG/`                   |
| `local_PCA_GHSL`     | `<SCRATCH>/04_clair3_phased_GHSL/`       |

Both will write to their own sibling output trees
(`<SCRATCH>/path_localpca_thetapi/` and `<SCRATCH>/path_localpca_GHSL/`)
following the same 01–09 layout as path 1.

Both will consume the SAME `<SCRATCH>/_shared/sample_metadata.tsv` built
by `local_PCA_z/08_atlas_json/08a_build_sample_metadata.R`.

The atlas viewer drag-drops one JSON per path per chromosome (so `LG28`
gets three JSONs total — one per path — and the atlas reconciles them
into a unified per-chromosome page).
