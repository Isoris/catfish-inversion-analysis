# STEP_TR_A + STEP_TR_B θπ scaling patch — Apr 30, 2026

## Why

ANGSD's `pestPG` `tP` column is the **sum of per-site θπ across sites with
data in the window**, not a per-site density. To compare windows or
samples you must divide by `nSites` (the per-window callable-site count,
last column of pestPG). Confirmed in:

- ANGSD GitHub issue #329 (Callithrix-omics, clairemerot, ANGSD owner)
- Korunes & Samuk 2021 (pixy) — the missing-data-aware denominator
- Direct trace of `thetaStat do_stat` source (`slice()` is a plain sum)

For our 226-sample pure *C. gariepinus* hatchery cohort at ~9× coverage,
`nSites` varies sample-to-sample due to coverage dropouts. Without the
division, samples with more callable sites in a window falsely look more
diverse — a coverage artifact masquerading as biology.

## Files in this bundle

| File | What |
|---|---|
| `00_sanity_check_pestPG_scaling.sh` | Run first on one real `.pestPG` to confirm the bug exists in your data. Takes 2 seconds. |
| `STEP_TR_A_compute_theta_matrices.R.patched` | Drop-in replacement. Adds per-site `theta_pi` (the divided value) + preserves `tP_sum` (raw window sum) as a new TSV column. |
| `STEP_TR_B_classify_theta.R.patched` | Drop-in replacement. Reads both columns, builds two matrices, emits both per-sample tracks to the page-12 JSON with a top-level `available_modes` toggle. **Analysis layers (envelopes, local PCA, \|Z\|) use the per-site track only**, since that's the diversity-comparable estimate. The raw-sum track is for diagnostic display only. |
| `STEP_TR_A_patch.diff` | Diff for code review / version control. |

## Workflow on LANTA

```bash
cd $CODEBASE/phase_2_discovery/2f_theta_discovery

# 1. Confirm bug exists in real data
bash 00_sanity_check_pestPG_scaling.sh \
  $PESTPG_DIR/CGA001.win10000.step2000.pestPG
# Expect: mean(tP) ≈ 5, mean(tP/nSites) ≈ 5e-4, VERDICT: SUM

# 2. Backup + drop in patched scripts
cp STEP_TR_A_compute_theta_matrices.R STEP_TR_A_compute_theta_matrices.R.bak.preTPfix
cp STEP_TR_B_classify_theta.R         STEP_TR_B_classify_theta.R.bak.preTPfix
cp /path/to/STEP_TR_A_compute_theta_matrices.R.patched STEP_TR_A_compute_theta_matrices.R
cp /path/to/STEP_TR_B_classify_theta.R.patched         STEP_TR_B_classify_theta.R

# 3. Re-run on LG28 only as part of the dry-run
source 00_theta_config.sh
$RSCRIPT STEP_TR_A_compute_theta_matrices.R C_gar_LG28
$RSCRIPT STEP_TR_B_classify_theta.R         C_gar_LG28
```

## What changes downstream

### TSV (STEP_TR_A output)

```
before:  sample chrom window_idx start_bp end_bp theta_pi  n_sites
after:   sample chrom window_idx start_bp end_bp theta_pi  tP_sum  n_sites
                                                  ↑
                                              now per-site (tP/nSites)
```

### JSON layer `theta_pi_per_window` (STEP_TR_B output)

New top-level fields:
```json
"available_modes": ["per_site", "raw_sum"],
"default_mode":    "per_site",
```

Each `samples[i]` entry now has both tracks:
```json
{
  "sample_id": "CGA001",
  "theta_pi":  [0.0023, 0.0019, ...],   // per-site, default
  "tP_sum":    [22.4,   18.7,   ...],   // raw sum, diagnostic
  "n_sites":   [9742,   9698,  ...]
}
```

### JSON layers `theta_pi_local_pca`, `theta_pi_envelopes`, `tracks`

Unchanged. All computed from per-site `theta_pi` only — that's the
diversity-comparable estimate. The raw-sum track exists for display, not
analysis.

## Atlas-side TODO (when page-12 renderers get wired)

The page-12 θπ panel should expose a small toggle widget:

```
[ per-site ] [ raw sum ]    ← reads available_modes from JSON
```

When `per_site`: render `samples[i].theta_pi` on the per-window scrubber
and the cohort heatmap.

When `raw_sum`: render `samples[i].tP_sum` on the same panels.

Envelopes, the |Z| ribbon, the local-PCA scatter, and the cross-method
concordance overlay **do not flip** — they always reference per-site.
Document this in the toggle's tooltip so reviewers don't get confused.

Diagnostic value of `raw_sum`: regions where the cohort lost callable
sites show up as pronounced dips in the raw-sum track that the per-site
track flattens (by construction). Useful for spotting low-mappability
regions, repeat-masked stretches, or regions where the cohort's filter
regime was unusually aggressive.

## Sanity check after patching

In R, after STEP_TR_A finishes on LG28:

```r
library(data.table)
dt <- fread("$THETA_TSV_DIR/theta_native.C_gar_LG28.win10000.step2000.tsv.gz")

# theta_pi should now be in the 1e-4 to 1e-2 range
summary(dt$theta_pi)

# tP_sum should be in the 1 to 50 range
summary(dt$tP_sum)

# Sanity: theta_pi == tP_sum / n_sites (where n_sites > 0)
all.equal(dt[n_sites > 0, theta_pi],
          dt[n_sites > 0, tP_sum / n_sites])
```
