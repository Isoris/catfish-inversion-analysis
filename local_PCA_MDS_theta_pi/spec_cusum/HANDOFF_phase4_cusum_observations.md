# Handoff — Per-sample CUSUM observations for boundary refinement
*Session ending 2026-05-01 · phase 4 / 4b θπ resolution + plans for 4a GHSL*

---

## What this work is for

In your hatchery cohort (226 *C. gariepinus*, ~9× WGS), an inversion candidate has many carriers, and **each carrier's curve on a per-sample × per-window matrix has its own changepoint position**. The aggregate cohort-mean signal gives one boundary estimate; the **distribution of per-carrier changepoint positions** gives:

- A boundary estimate that's empirical, not parametric (median/IQR of a 60-sample distribution rather than a model fit)
- The shape of the distribution itself — whether all carriers vote for one position (sharp), spread across a wide range (eroded), or split into two narrow clusters at different positions (your hatchery-LD multiplexing hypothesis)
- Per-carrier `(cp_bp, strength, asymmetry)` rows that can be **stacked across streams** (CUSUM-on-θπ + CUSUM-on-GHSL + ancestral-fragments) instead of pre-aggregated within-stream and then merged

The grain of this work is **observations per carrier per stream**, not per-stream summaries. The downstream consensus step (currently `phase_4c / 03_consensus_merge.R`, formerly `phase_6_breakpoint_refinement / 03_consensus_merge.R`) gets redesigned later to aggregate across the pooled per-carrier evidence.

---

## What's done

### `lib_persample_cusum.R` (the math)

**Path:** `phase_4_resolution/shared_lib/lib_persample_cusum.R`
**Size:** 149 lines
**Tests:** 9/9 unit tests pass via `test_lib_persample_cusum.R`

Pure utility, observation-only. Takes any sample × window numeric matrix + a `window_pos_bp` vector + a `candidate_id`. Returns a single `data.table` with one row per sample:

| Column | Meaning |
|---|---|
| `sample_id` | Row name from the matrix |
| `candidate_id` | Echoed from input |
| `cp_idx` | Window index of strongest changepoint in this sample's row |
| `cp_bp` | BP coordinate of that changepoint |
| `strength` | `max\|cumsum(centred)\| / sd(centred)` — scale-invariant CUSUM stat |
| `asymmetry` | Sign of the cumsum at cp (+1 = signal rises past, −1 = falls past) |
| `left_mean` / `right_mean` | Mean of the row before / after cp |
| `n_used` | Non-NA windows for this sample |
| `informative` | TRUE if `strength ≥ 2.0 × sqrt(n_used)` and `n_used ≥ 5` |

**Math:** one-sided CUSUM with mean-removal (Lancaster MATH337 §1.3). For row x of length W: subtract row mean, take cumsum of the centred series, find argmax of |cumsum| → that's the changepoint.

**The `informative` threshold is empirically calibrated**, not parametric. Under pure noise, `max|cumsum|/sd ~ sqrt(n)` (random walk variance scaling). Quantiles of `(max|cumsum|/sd)/sqrt(n)` from 2000 reps of pure rnorm: 0.79 / 1.29 / 1.59 / 1.93 at p=0.5/0.95/0.99/0.999. Threshold at 2.0 keeps <1/1000 noise rows; real mean-shift signal scores ~5× sqrt(n).

**What was deliberately removed in this session:**
- KDE peak finding on cohort cp_bp distribution (assumed Gaussian shape)
- Silverman bandwidth (also assumes Gaussian)
- Hartigan dip test for bimodality (uses a parametric null)
- Spread classification (`tight`/`moderate`/`ragged` thresholds at 100/500 kb)
- Per-mode aggregation

These all bake in distribution-shape priors that the actual biology may not match. The right place for shape interpretation is the consensus step, not the per-stream observation step. See "Why" section below.

### `STEP_T05_theta_cusum.R` (the θπ driver) ← runnable on LANTA

**Path:** `phase_4_resolution/4b_theta_resolution/STEP_T05_theta_cusum.R`
**Tests:** End-to-end test passes 7/7 via `test_step_t05_e2e.R` (synthetic JSON + real driver invocation + output assertions).

Reads a TR_B-output JSON (`<chrom>_phase2_theta.json` containing the `theta_pi_per_window` layer), reconstructs the sample × window matrix from the flat `values` array, slices per candidate, runs `persample_cusum`, writes per-sample TSV + per-candidate summary TSV.

**Run it on LANTA:**
```bash
Rscript phase_4_resolution/4b_theta_resolution/STEP_T05_theta_cusum.R \
    --json   /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/.../LG28_phase2_theta.json \
    --candidates /path/to/candidate_intervals.tsv \
    --out-dir  /path/to/cusum_output/LG28/ \
    --mode per_candidate
```

The candidate TSV needs columns `candidate_id`, `chrom`, `start_bp`, `end_bp` — same schema as `02_ancestral_fragments.R`. If you don't have one, use:
```bash
--mode whole_chrom
```
to run on the entire chromosome (one cp per sample, anywhere on the chrom) without a candidate list.

**Outputs:**
- `theta_cusum_per_sample.tsv.gz` — one row per (sample × candidate). 18 columns total. Stackable with `02_ancestral_fragments.R`'s per-sample TSV when you redesign consensus.
- `theta_cusum_summary.tsv` — one row per candidate. Empirical distribution shape (n_total / n_informative / cp_min / cp_q25 / cp_median / cp_q75 / cp_max / cp_iqr_kb / cp_mad_kb / n_left / n_right / n_asym_pos / n_asym_neg). **No KDE, no fits — pure observations.**

### Atlas plumbing turns 114-118

Already shipped in this session before the CUSUM work began. Atlas now correctly:
- Detects all 8 phase-2 enrichment layers (4 θπ: `theta_pi_per_window/local_pca/envelopes/d17_envelopes`; 4 GHSL: `ghsl_local_pca/envelopes/secondary_envelopes/d17_envelopes`)
- Recovers them via the recovery sweep
- Routes them through `inferLayersFromV1` for legacy schema_v1 inputs
- Merges them via the case block in `mergeEnrichmentLayers` (was the bug — `default: continue` silently dropped 8 layers)
- Registers them in `_SCHEMA_REGISTRY` so the schema badge popup shows all 10 phase-2e items

`STEP_TR_B_classify_theta.R` was modified to emit `theta_pi_per_window.values` as a flat row-major array + `sample_ids` as a top-level string array (atlas-canonical contract). Schema bumped to v2.

**Atlas test suite:** 30/30 passing as of session end.

---

## What's NOT done — explicit pickup list

### Highest priority: see real CUSUM output before any further design

The session's most important decision was to stop pre-aggregating in the CUSUM streams. **Before you build the GHSL driver or the atlas rendering, run T05 on a real LG28 θπ JSON.** Open `theta_cusum_summary.tsv` and look at the actual distribution stats. Open `theta_cusum_per_sample.tsv.gz` and look at the per-carrier cp_bp values for the LG28 candidate (60 / 106 / 60 karyotypes).

You need to see:
- Are 60 informative carriers showing up? (60 minor homozygotes + 106 hets = up to 166 potential carriers depending on which homozygote class is "minor")
- Does the cp_bp distribution look bimodal in real data, supporting your hatchery-LD insight?
- Does cp_side_inferred actually split sensibly between left/right boundaries?
- Is `strength` separation between carriers and non-carriers as clean as on synthetic data?

**Don't design the consensus step or the atlas visualization until you've looked at this output.** I can suggest things based on synthetic data; you have to look at the real data to know if the suggestions match biology.

### Next: GHSL driver `STEP_R02_ghsl_cusum.R`

Mirror of T05 with the GHSL matrix instead of θπ. Quentin said "for GHSL the hpc is back but I'm tired maybe do in few days bc its hard." Picking this up later.

The differences from T05:
- Source matrix: `ghsl_panel.div_roll` from a GHSL phase-2 JSON (or whatever the canonical export field is — check `STEP_C04b` or the GHSL export step)
- Window grid may differ from θπ's window grid (GHSL has its own sliding-window scheme)
- Otherwise identical: read JSON, build sample × window matrix, slice per candidate, call `persample_cusum`, decorate with the same column schema, write `ghsl_cusum_per_sample.tsv.gz` and `ghsl_cusum_summary.tsv`
- **Important:** keep the column schema identical to T05's so the two streams can be `rbind`ed for any future stacked-evidence consensus design. Just change `stream = "theta"` to `stream = "ghsl"`.

Estimated effort: 90 minutes including E2E test if the GHSL JSON schema is similar to TR_B's. Longer if the GHSL panel schema is different and needs translation.

### After both drivers run on real data: redesign consensus

Currently `phase_4c / 03_consensus_merge.R` (former `phase_6_breakpoint_refinement / 03_consensus_merge.R`) has these problems:

1. **Assumes one inversion = one start + one end per candidate.** With hatchery LD multiplexing two ancestral configurations, you might want two left and two right boundaries, or a flag indicating multimodality with both candidates reported.
2. **Pre-aggregates each stream before merging.** The 226 per-carrier observations from `02_ancestral_fragments.R` collapse to a single mode + CI before going into 03; CUSUM streams would do the same. Per-carrier richness is thrown away.
3. **`1.96/sqrt(n)` CI formula** assumes a normal sampling distribution of n=4-7 weighted point estimates that aren't even in commensurate units of evidence. The CI shape is dominated by the formula, not the data.
4. **SV (step37) didn't pan out empirically** but is still wired in at weight 0.5.

A redesigned consensus step would:
- Stack all per-carrier observations across streams: `02_ancestral_fragments_per_sample.tsv.gz` + `theta_cusum_per_sample.tsv.gz` + `ghsl_cusum_per_sample.tsv.gz`. Same grain (one row per sample × stream × candidate). Add a `stream` column to `02`'s output to match.
- For each candidate, look at the pooled cp_bp distribution. Empirical histogram. Modal position(s) by simple peak-counting on a fixed-bin histogram (no kernel choice).
- Drop the `1.96/sqrt(n)` CI. Use direct bootstrap on the observation pool (resample carriers with replacement, recompute mode) — same as `02`'s existing within-stream bootstrap, just on the stacked pool.
- Drop the SV stream until it's shown to add signal.
- Output: per-candidate TSV with `final_left_bp` / `final_right_bp` / CI per side **plus a flag** for multi-modal cases that emits multiple breakpoint candidates.

**Don't write this redesign until real-data CUSUM output is in.** The right design depends on what the actual distributions look like.

### Atlas rendering (page 3 + page 12 + boundaries page)

Quentin's design from the last image set (3 jpgs, image 1 = dosage heatmap with vertical stripe artifacts, image 2 = Plasma colormap legend "Regional het (dosage)", image 3 = per-sample heatmap with dashed grid lines):

**Three-row layout for page 3 (GHSL) and page 12 (θπ):**

1. **Top:** cohort CUSUM curve.
   - x = position
   - y = mean(|cumsum(centred row)|) across informative carriers, OR median, OR something simpler — needs real-data look first
   - Yellow vertical lines at "shape change" positions — TBD criterion (peaks of the curve? inflections? something simpler?)

2. **Middle:** representative sample lines.
   - Pick ~6-12 carriers via k-means or MDS on the per-sample CUSUM observations (`cp_bp`, `strength`, `asymmetry`)
   - Each line colored by group (Homo_1 / Het / Homo_2 = STEP21 karyotype)
   - Yellow tick on each line at that sample's `cp_bp`

3. **Bottom:** per-sample heatmap (image-3 style).
   - Rows = 226 samples, ordered by group then by `cp_bp`
   - x = window position
   - Color = signal value, **Plasma colormap matching image 2** (yellow = high diversity, dark purple = low)
   - Yellow tick per row at that sample's cp_bp
   - Horizontal black separators at group boundaries

**Boundaries page** = combined view: stack per-carrier observations from all streams (θπ CUSUM + GHSL CUSUM + ancestral fragments) as rows of cp_bp votes, shown as a stacked histogram or rug plot per stream, plus the consensus call as a vertical line.

**Don't build this yet** for the same reason — atlas visualization choices are easier to get right looking at real data.

---

## Why we walked back from the original CUSUM design

In the first version of `lib_persample_cusum.R`, I had it return three data.tables: `persample`, `modes`, `summary`. The `modes` table was KDE peak-finding on the cohort cp_bp distribution with Silverman bandwidth (50 kb floor) plus a Hartigan dip test for bimodality.

Quentin pushed back, correctly:

> "How can choose the statistical test of a distribution shape that we don't know what it will look like. To me your logic although sounds based on what we discussed its kind of not scientific like. Observe first empirically then correct the kernel."

The KDE + dip test approach assumes:
- The cohort cp_bp distribution is roughly Gaussian-mixture-like (so Silverman's rule gives a sensible bandwidth)
- Multimodality is the relevant binary distinction (one-mode vs many-mode)
- A 50 kb bandwidth floor is empirically motivated (it wasn't — I picked it from my prior, not from data)

In a hatchery cohort with founder LD, the actual distribution might be:
- Two narrow clusters (founder lineage A breakpoint vs founder lineage B breakpoint, both inherited as blocks)
- One narrow cluster + a long tail (one canonical breakpoint + recombinant-eroded carriers)
- Diffuse but unimodal (uniformly aged inversion, gene-conversion erosion)
- Discrete clusters at 2-3 specific positions (multiple ancestral haplotypes still segregating)

The shape determines what consensus means, and we don't know the shape until we look. So: gather observations, don't fit. The simplification removed:
- 4 helper functions (`.find_kde_modes`, `.assign_to_modes`, `.classify_spread`, `.dip_test`)
- 2 of the 3 returned data.tables (`modes`, `summary`)
- 4 magic-number constants (`KDE_BW_MIN_BP`, `SPREAD_TIGHT_BP`, `SPREAD_MODERATE_BP`, `MIN_SAMPLES_FOR_KDE`)

Lib went from 318 lines → 149 lines. Tests went from 17 checks (with 7 initial failures from un-tuned thresholds, then all passing after empirical calibration) → 9 checks all passing on the trimmed lib. End-to-end test added (7/7).

The empirical noise threshold for `informative` was kept because it's calibrated against pure noise data (not against a model), and because filtering noise rows out of the observation table is genuinely useful for downstream consumers (otherwise every carrier and every non-carrier emits a row with a `cp_bp`, and the non-carrier cp_bp is just where the random walk happened to peak). If you want to drop even the `informative` filter and let the consensus step apply its own, the lib still emits `strength` and `n_used` so any downstream filter is reproducible.

---

## File manifest

```
/home/claude/work/phase_4_resolution/
├── shared_lib/
│   ├── lib_persample_cusum.R              [149 lines, 9/9 unit tests pass]
│   └── test_lib_persample_cusum.R         [unit tests]
├── 4a_ghsl_resolution/                    [empty — STEP_R02 not yet built]
├── 4b_theta_resolution/
│   ├── STEP_T05_theta_cusum.R             [the θπ driver, ready to run on LANTA]
│   └── test_step_t05_e2e.R                [end-to-end test, 7/7 passes]
└── 4d_dual_clustering/                    [empty — DC06 not yet built]
```

Staged in `/mnt/user-data/outputs/`:
- `lib_persample_cusum.R`
- `test_lib_persample_cusum.R`
- `STEP_T05_theta_cusum.R`
- `test_step_t05_e2e.R`
- This handoff doc

Plus the earlier turn-114-to-118 outputs (already staged earlier in the session):
- `Inversion_atlas.html` (49,990 lines, 8 phase-2 enrichment layers wired through detectors / recovery / inferV1 / merge / registry)
- `STEP_TR_B_classify_theta.R` (v2 schema with flat `values` + `sample_ids`)
- `test_atlas_t114_*.js` through `test_atlas_t118_*.js` (atlas plumbing tests)
- `test_tr_b_harness.R`
- Earlier `HANDOFF_turns_114_to_118_phase2_enrichment_plumbing.md`

---

## Quick command reference

**Run unit tests on the lib:**
```bash
Rscript /home/claude/work/phase_4_resolution/shared_lib/test_lib_persample_cusum.R
```

**Run end-to-end test on the θπ driver:**
```bash
Rscript /home/claude/work/phase_4_resolution/4b_theta_resolution/test_step_t05_e2e.R
```

**Run T05 on real LANTA data:**
```bash
RSCRIPT=/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript
BASE=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04
JSON=$BASE/.../json_out/C_gar_LG28/C_gar_LG28_phase2_theta.json
CANDIDATES=$BASE/.../candidate_intervals.tsv
OUT=$BASE/.../cusum_output/C_gar_LG28/

$RSCRIPT $BASE/inversion_codebase_v8.5/phase_4_resolution/4b_theta_resolution/STEP_T05_theta_cusum.R \
    --json   $JSON \
    --candidates $CANDIDATES \
    --out-dir  $OUT \
    --mode per_candidate
```

If no candidate TSV exists yet:
```bash
$RSCRIPT ... --json $JSON --out-dir $OUT --mode whole_chrom
```

The driver locates the lib at `../shared_lib/lib_persample_cusum.R` relative to itself by default. Pass `--lib /absolute/path/to/lib_persample_cusum.R` to override.

---

## Constraints to remember

- 226 *C. gariepinus* hatchery cohort, **NOT F1 hybrid**, **NOT C. macrocephalus**. K clusters in this cohort reflect broodline structure, not species admixture. (See userMemories.)
- `MS_Inversions_North_african_catfish` is the manuscript these refinements feed.
- LANTA account/partition: `lt200308` / `compute`. Rscript at `/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript`.
- Quentin works across multiple chats simultaneously — version everything, write handoffs, don't assume the next session has full context.
- "Continue" from Quentin = next agreed item. Multiple "idk"s in a row = stop and document, don't push through.

---

End of session. Nothing on fire. Picking back up: run T05 on LG28, look at the real distribution, then decide what GHSL / atlas / consensus should look like.
