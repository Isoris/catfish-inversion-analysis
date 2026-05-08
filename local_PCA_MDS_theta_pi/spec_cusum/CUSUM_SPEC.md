# CUSUM additions to GHSL and θπ pages — full specification

**Status:** spec recovered from the "phase 4 resolution" architecture chats
plus the Lancaster MATH337 changepoint detection notes
(`https://www.lancaster.ac.uk/...MATH337...`).
**Source mapping:**

- Architecture (R02 / T05 / DC06 / SC01): the other-Claude conversation,
  screenshots batch 1–2 of this session
- Math (CUSUM statistic, threshold, location estimate): Lancaster MATH337
  §1.3, screenshots batch 3 of this session
- Phase placement (4a_ghsl_resolution / 4b_theta_resolution / 4d_dual_clustering):
  screenshots batch 2 of this session

## 1. The biological insight

Boundary refinement at the **group / cohort level** (existing
phase 5 STEP_Q07d / Q07e) tells you "the inversion spans 14.0–19.5 Mb".
That's the manuscript number — but it's a single span, derived from
where Fst between groups returns to baseline.

CUSUM at the **per-sample / per-carrier level** tells you the
**distribution of boundary positions across carriers**. That's a
different biological object:

- A tight 5' boundary distribution (e.g. median 14.01 Mb, IQR 70 kb)
  with a ragged 3' (e.g. median 19.4 Mb, IQR 1.2 Mb) is consistent
  with **gene-conversion erosion at the 3' end and a sharp left
  breakpoint** — clean inversion mechanics
- A bimodal 3' distribution (45 carriers cluster at 18.95 Mb,
  9 carriers internal at 17.2 Mb) is consistent with **nested
  inversions or alternate-allele structure** — recurrent origin
  or heterogeneous breakpoints

These are **manuscript-grade observations** that group-mean methods
cannot produce. They land naturally on the Ancestry / shelf-QC page
when zoomed to a candidate region, populating the "boundaries panel"
with per-carrier evidence rather than a single span.

## 2. Three CUSUM scripts, one shared algorithm, three calibrations

The walked-back final count is **three CUSUM applications** (per the
other-Claude self-correction in the screenshots), each tuned for
what it consumes:

| Script | Phase folder | Input matrix | Signal flavor | Baseline | Threshold |
|---|---|---|---|---|---|
| `STEP_R02_ghsl_cusum.R` | `4a_ghsl_resolution/` | GHSL rolling-divergence (per-sample × per-window) | rate (fraction of phased het sites) | per-sample mean across candidate | calibrated for normalized fraction |
| `STEP_T05_theta_cusum.R` | `4b_theta_resolution/` | θπ rolling matrix (per-sample × per-window) | rate (per-bp diversity) | per-sample mean across candidate | calibrated for ANGSD-derived noise (includes invariant-site denominator) |
| `STEP_SC01_emit_subcandidates.R` cell cliff-walker | `4e_subcandidate_emission/` | per-cell mean tracks | difference (cell mean − rest-of-row mean) | zero (signal is already a difference) | median + 3×MAD of flank-Δ distribution |

R02 and T05 share the algorithm exactly. SC01's cell cliff-walker
**reuses the CUSUM kernel** but operates on a different signal
flavor (already-differenced rather than rate), so its threshold logic
differs even though the inner loop is identical.

**Architectural commitment:** write the CUSUM kernel **once** as a
shared utility (e.g. `phase_4_resolution/shared_lib/cusum_core.R`),
then R02, T05, and the SC01 cell-walker each call it with their own
calibration parameters. This avoids three independent reimplementations
that would drift apart.

## 3. The CUSUM math (canonical, from Lancaster MATH337 §1.3)

For a per-sample series y = (y₁, …, yₙ) within a candidate region:

### 3.1 The statistic at candidate position τ

```
C_τ = sqrt( τ(n−τ) / n ) · | ȳ₁:τ − ȳ_(τ+1):n |
```

where ȳ₁:τ is the mean of the pre-τ segment and ȳ_(τ+1):n is the
mean of the post-τ segment.

The sqrt term is a re-scaling so that under no-change, |C_τ| is
distributed as the absolute value of a standard normal scaled by σ.

### 3.2 Maximum and detection

```
C²_max = max_{τ ∈ {1, …, n−1}} C²_τ / σ²
```

Detect a changepoint if `C²_max > c` for chosen threshold c.

### 3.3 Estimate the changepoint location

```
τ̂ = argmax_{τ} C²_τ
```

This is the position that maximizes the CUSUM statistic (the
"peak" in the CUSUM trace). The strength of the change is then
`Δμ̂ = ȳ_(τ̂+1):n − ȳ_1:τ̂`.

### 3.4 O(n) implementation via running sums

The naive implementation is O(n²) because every τ recomputes both
segment means. The standard trick uses prefix sums:

```
S_n = Σᵢ yᵢ                  # total sum (compute once)
S = 0                         # running sum
C_max = 0
τ̂ = 0
FOR t = 1, ..., n−1:
    S += y_t                  # running pre-τ sum
    ȳ_1:t = S / t
    ȳ_(t+1):n = (S_n − S) / (n − t)
    C²_t = (t(n−t) / n) · (ȳ_1:t − ȳ_(t+1):n)²
    IF C²_t > C_max:
        C_max = C²_t
        τ̂ = t
RETURN τ̂, C_max
```

**This is O(n) per sample.** For 226 samples × ~16,500 windows on
LG28, the total compute is ~3.7M operations — trivial in JS or R.
Even per-candidate (typically ~1500 windows) it's ~340k operations.
Sub-second.

## 4. Calibration differences between R02 and T05

The kernel is identical; the calibration is per-script:

### 4.1 R02 (GHSL)

- **Input matrix:** rolling-divergence values per (sample, window).
  Values are bounded fractions of phased het sites, typically [0, 1].
- **Baseline:** per-sample mean across the candidate region.
  CUSUM walks deviations from that baseline.
- **σ² estimate:** sample-level robust variance (MAD²) computed
  on the candidate-region residuals (y − ȳ).
- **Threshold c:** quantile-based on the chromosomal null
  distribution of `C²_max` from non-candidate windows
  (i.e. compute `C²_max` on random non-candidate spans of
  matching length, take the 95th or 99th percentile as c).

### 4.2 T05 (θπ)

- **Input matrix:** rolling θπ values per (sample, window).
  Values are per-bp diversity rates from ANGSD — different noise
  structure (includes invariant-site denominator).
- **Baseline:** per-sample mean across the candidate region.
- **σ² estimate:** ANGSD-aware. Either MAD² of residuals
  (matches R02), or use ANGSD's per-site variance estimates
  if available. **Open question** — defer to first real run on
  LG28 to see whether MAD-based variance is sufficient.
- **Threshold c:** calibrated separately from R02 because the
  noise scale differs. Same procedure (chromosomal null
  distribution), but values won't match R02's.

### 4.3 SC01 cell cliff-walker

- **Input:** per-cell mean tracks where the signal is
  `cell_mean − rest-of-row_mean` (already a difference).
- **Baseline:** **zero**, not the cell's own mean. The signal is
  already centered.
- **Threshold:** `median + 3×MAD` of the flank-Δ distribution
  (not an absolute changepoint strength). This matches the
  "cliff" terminology — it's looking for a flank-to-flank
  drop, not a level shift.

The shared kernel takes parameters `(y, baseline_mode, sigma2_mode,
threshold_c)` and dispatches accordingly.

## 5. Output contracts (three TSVs per candidate, per signal)

This matches what the other-Claude already specified. Reproducing
verbatim from the screenshots so it's locked:

### 5.1 Per-sample changepoint — `<cid>_<signal>_cusum_persample.tsv`

```
sample_id    cp_position_bp    cp_strength    asymmetry
CGA_001      18_950_000        0.73           +0.42
CGA_007      18_945_000        0.69           +0.39
...
```

- `cp_position_bp` — argmax τ̂ converted to genomic bp
- `cp_strength` — normalized |CUSUM| at the peak (`C_max / sqrt(σ²)`
  or similar — bounded for cross-sample comparison)
- `asymmetry` — sign and magnitude of the mean shift Δμ̂.
  Positive = post-τ mean higher than pre-τ; negative = post-τ lower.

One row per sample. For non-carriers (or samples with no significant
changepoint), include a row with NA `cp_position_bp` and NA strength —
**don't drop them**, downstream code needs to know which samples had
no signal vs which weren't analyzed.

### 5.2 Per-candidate boundary distribution — `<cid>_<signal>_boundary_dist.tsv`

```
boundary_side   n_carriers   median_bp     iqr_bp    spread_class
5_prime         54           14_010_000    70_000    tight
3_prime         54           19_400_000    1_200_000 ragged
```

Computed from the per-sample changepoint distribution. The
`spread_class` field is a derived label:

- `tight` — IQR < 100 kb (single-origin clean inversion)
- `intermediate` — 100 kb ≤ IQR < 500 kb
- `ragged` — IQR ≥ 500 kb (gene-conversion-eroded, recurrent,
  or recombinant-rich)

Thresholds tunable; defaults from the other-Claude spec.

### 5.3 Per-candidate sub-system clusters — `<cid>_<signal>_subsystems.tsv`

```
cluster_id   n_samples   center_bp     span_bp   member_samples
1            45          18_950_000    150_000   CGA_001,CGA_007,...
2             9          17_200_000     80_000   CGA_032,CGA_054,...
```

From clustering changepoint positions across samples (k-means
or hierarchical with silhouette selection for k). Two clusters
in one boundary's distribution = nested inversion or
alternate-allele signature.

Computed for each boundary side independently (5' clusters,
3' clusters), reported together.

## 6. DC06 cross-signal concordance

`STEP_DC06_cusum_concordance.R` in `4d_dual_clustering/`. Loads R02
and T05 outputs, compares per-sample changepoint positions:

```
sample_id    ghsl_cp_bp    theta_cp_bp    delta_bp    concordance_class
CGA_001      18_950_000    18_960_000     10_000      concordant_50kb
CGA_007      18_945_000    18_990_000     45_000      concordant_50kb
CGA_032      17_200_000    18_900_000   1_700_000     ghsl_sharper
CGA_054      18_950_000    17_180_000   1_770_000     theta_sharper
...
```

Concordance classes:

- `concordant_50kb` — |delta_bp| ≤ 50 kb (manuscript-grade
  cross-validation)
- `concordant_500kb` — |delta_bp| ≤ 500 kb
- `ghsl_sharper` — GHSL detected a sharper boundary than θπ
  (suggestive of phasing-block-resolved boundary structure)
- `theta_sharper` — θπ detected a sharper boundary than GHSL
  (suggestive of θπ-resolved boundary structure not visible
  in phasing geometry)
- `discordant` — no agreement at any threshold

Reports the manuscript paragraph the other-Claude proposed verbatim:

> "Per-sample changepoint positions from GHSL and θπ analyses
> showed concordance within 50 kb for X of N carriers (Y%),
> supporting common biological boundaries. Discordant carriers
> split between samples where GHSL detected a sharper boundary
> than θπ (suggestive of phasing-block-resolved boundary structure)
> and the reverse (suggestive of θπ-resolved boundary structure
> not visible in phasing geometry)."

## 7. Where this lives in the atlas

Per ADR-13, **two new pages** receive CUSUM output:

### 7.1 GHSL page (existing, page 3)

Currently shows: chromosome-wide |Z| profile, sim_mat, candidate envelopes.

**Adds when zoomed to a candidate region:**

- **Per-carrier boundary panel** — for each candidate, render the
  5' and 3' boundary distributions as histograms or strip plots
  along the genomic axis, with median + IQR annotated
- **Sub-system cluster overlay** — color-code carriers by
  cluster_id from the subsystems TSV; show cluster centers
  as vertical reference lines
- **Spread-class badge** — `tight` / `intermediate` / `ragged`
  label per boundary side

This information is loaded only when the user has zoomed to a
candidate (lazy-load `LG##_ghsl_cusum.json` on candidate selection).

### 7.2 θπ page (existing, page 12)

Same structure as GHSL page additions, populated from T05 outputs.
Different file: `LG##_theta_cusum.json`.

### 7.3 Ancestry / shelf-QC page (new, per ADR-13)

The cross-signal view. Renders DC06's concordance grid plus
per-carrier overlays from both R02 and T05 simultaneously
(side-by-side or interleaved on the same genomic axis).

This is the page where the manuscript-grade observations live.

## 8. Computation: cluster-side vs JS-side

The user asked: **"could it be done in JS?"** Yes — but it shouldn't be.

### 8.1 What's cheap enough for JS

The CUSUM kernel itself is O(n) and runs in milliseconds on a
typical browser for one candidate. **In principle, JS could
recompute CUSUM on the fly when the user changes:**

- The candidate region (different zoom)
- The sample subset (lasso selection)
- The baseline mode (whole-cohort vs carrier-only)

That's actually useful — interactive re-CUSUM on user-defined
sample subsets is a real value-add over static cluster-side output.

### 8.2 What requires cluster-side computation

The **input matrices** themselves can't reasonably ship to the browser:

- GHSL rolling-divergence matrix: 226 samples × ~16,500 windows ×
  Float32 = ~15 MB per chromosome
- θπ rolling matrix: similar
- Per-chromosome × 28 chromosomes = ~420 MB

That's too much for client-side initial load, especially given the
atlas already loads multiple JSONs. **The cluster precomputes the
default CUSUM output (R02 / T05 TSVs) for the default sample-set
configuration.** If the user wants to re-CUSUM on a different
sample subset, the JS can reload the relevant matrix on-demand.

### 8.3 Recommended split

- **Cluster-side, always:** `STEP_R02_ghsl_cusum.R` and
  `STEP_T05_theta_cusum.R` produce the default output TSVs,
  packed into `LG##_ghsl_cusum.json` and `LG##_theta_cusum.json`
  as per-candidate layers
- **Cluster-side, on-demand:** if the user requests re-CUSUM for
  a custom sample subset, the cluster recomputes against the
  full matrix and streams a new TSV. Implementation detail —
  could be a lightweight web service or a "download these results
  with my custom subset" workflow
- **JS-side, always:** rendering of histograms, strip plots,
  cluster overlays, concordance grids from the loaded JSON.
  No CUSUM math in JS by default.
- **JS-side, optional later:** if it turns out cluster-on-demand
  is a hassle, ship a small per-candidate matrix slice
  (~1500 windows × 226 samples × Float32 = ~1.4 MB per candidate)
  with the JSON, and recompute CUSUM in JS when the user changes
  sample subset. Defer this until the workflow needs it.

**Default plan: pure cluster-side compute, JS just renders.**
Re-CUSUM-in-JS is a future add-on, not the v1 scope.

## 9. Build order (revised from previous spec)

The other-Claude's order from the screenshots:

1. **STEP_T01 + T02** — θπ rolling matrix (heavy engine, paired
   with classifier). Already in scope; not new here.
2. **STEP_R02 split from existing C04b** — ~1 day. The GHSL
   rolling-divergence is already computed in C04b; R02 just
   adds the CUSUM step.
3. **STEP_T05 mirrors R02** — ~1.5 days. Same kernel, θπ-specific
   calibration. Build immediately after R02 so they're paired.
4. **DC06 concordance** — ~0.5 day. Requires R02 and T05 outputs.
5. **`shared_lib/cusum_core.R`** — written as part of R02, then
   extracted and reused for T05 and SC01 cell-walker.
6. **Atlas-side rendering** — gated on real LG28 data and on
   the existing GHSL / θπ pages working. Don't start until
   then.

**~3 days of cluster-side R work** to produce all CUSUM outputs.
**~2–3 days of atlas-side JS work** to add the boundary panels
and concordance grid, after the JSON contracts are validated.

## 10. JSON layer additions

### 10.1 `LG##_ghsl_cusum.json` (new file)

```json
{
  "schema_version": "1.0",
  "chromosome": "C_gar_LG28",
  "signal": "ghsl",
  "candidates": [
    {
      "candidate_id": "LG28_C_001",
      "candidate_start_bp": 14000000,
      "candidate_end_bp": 19500000,
      "n_carriers": 54,
      "per_sample": [
        {"sample_id": "CGA_001", "cp_position_bp": 18950000,
         "cp_strength": 0.73, "asymmetry": 0.42, "side": "3_prime"},
        ...
      ],
      "boundary_dist": {
        "5_prime": {"n_carriers": 54, "median_bp": 14010000,
                    "iqr_bp": 70000, "spread_class": "tight"},
        "3_prime": {"n_carriers": 54, "median_bp": 19400000,
                    "iqr_bp": 1200000, "spread_class": "ragged"}
      },
      "subsystems": {
        "5_prime": [{"cluster_id": 1, "n_samples": 54,
                     "center_bp": 14010000, "span_bp": 70000,
                     "member_samples": ["CGA_001", ...]}],
        "3_prime": [{"cluster_id": 1, "n_samples": 45,
                     "center_bp": 18950000, "span_bp": 150000,
                     "member_samples": [...]},
                    {"cluster_id": 2, "n_samples": 9,
                     "center_bp": 17200000, "span_bp": 80000,
                     "member_samples": [...]}]
      }
    }
  ]
}
```

### 10.2 `LG##_theta_cusum.json` (new file)

Same schema, `signal: "theta"`. Different numerical values but
identical structure.

### 10.3 Future: concordance JSON

`LG##_concordance.json` packs DC06 output. Per ADR-13 the Ancestry
page is the consumer; defer schema until the page design is
locked.

## 11. Walk-back protection

Specific failure modes to guard against in future sessions:

- **"Let's just compute CUSUM in JS to skip the cluster step"** —
  reject. The input matrix is too big for default load. JS can
  recompute on user request, but cluster-side is the source of truth.
- **"Three CUSUM scripts is too many, let's unify"** — partial
  reject. The kernel IS unified (`shared_lib/cusum_core.R`); the
  scripts differ in calibration only. That's the right design;
  don't collapse the calibration logic into the kernel.
- **"Let's make CUSUM the discovery detector for ancestry-Q too"** —
  reject. CUSUM detects level shifts in **rate** signals. Q-regime
  detection (QR03) needs dispersion / bimodality / per-sample
  deviation, not level shifts. Different problem, different tool.
  See QR03_SPEC.md §1 for the rationale.
- **"The cluster-side TSVs are enough, skip the JSON exporter"** —
  reject. The atlas needs the JSON to render. The TSVs are the
  intermediate representation; JSON is the deliverable.

## 12. tl;dr

- **Three CUSUMs total**, sharing one kernel: R02 (GHSL),
  T05 (θπ), SC01 cell cliff-walker (sub-candidate emission).
- **The math is canonical** (Lancaster MATH337 §1.3) and runs
  in O(n) per sample.
- **Cluster-side compute by default** (input matrices too big for
  browser). JS just renders histograms, strip plots, cluster
  overlays, and concordance grids.
- **Outputs land on three pages:** GHSL page (page 3) and θπ
  page (page 12) gain per-candidate boundary panels;
  Ancestry / shelf-QC page (new) shows DC06 cross-signal
  concordance.
- **Build order:** R02 → T05 → DC06, ~3 days cluster-side work.
  Atlas rendering ~2–3 days after JSONs validate.
- **Don't build until LG28 dry-run produces real GHSL / θπ
  output.** Same gating as everything else this session.
