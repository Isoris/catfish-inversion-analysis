#!/usr/bin/env bash
# =============================================================================
# 00_theta_config.sh — central config for phase_2_discovery/2f_theta_discovery
# =============================================================================
# Single source of truth for the θπ-discovery stream of Phase 2. Every
# variable consumed by STEP_TR_A and STEP_TR_B lives here. Sibling configs
# follow the same convention:
#   phase_2_discovery/2a_local_pca/00_inversion_config.sh
#   phase_2_discovery/2e_ghsl_discovery/00_ghsl_config.sh
#   phase_2_discovery/2h_local_ancestry_regime/00_q_engine_config.sh
#
# Source this from any STEP_TR_* script:
#   source "$(dirname "${BASH_SOURCE[0]}")/00_theta_config.sh"
#
# Cohort: 226 pure C. gariepinus hatchery samples (MS_Inversions_North_african_catfish).
# DO NOT use this config for the F1 hybrid genome assembly cohort or the
# pure C. macrocephalus wild cohort — those are separate manuscripts.
# =============================================================================

set -euo pipefail

# ── Project paths ─────────────────────────────────────────────────────────
export BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
# Repo layout (post April-2026 reorganization):
#   ${BASE}/inversion-popgen-toolkit/inversion_modules/phase_2_discovery/2f_theta_discovery/
# CODEBASE points at the inversion_modules root (where phase_2_discovery/ lives).
# To use the legacy v8.5 layout, override on the command line:
#   CODEBASE=$BASE/inversion_codebase_v8.5 source 00_theta_config.sh
export CODEBASE="${CODEBASE:-${BASE}/inversion-popgen-toolkit/inversion_modules}"
export THETA_DISC_DIR="${CODEBASE}/phase_2_discovery/2f_theta_discovery"
export OUTROOT="${BASE}/phase_2_outputs/2f_theta_discovery"

# ── Reference / sample inputs (shared with sibling 2a/2e streams) ─────────
export REF="${BASE}/01_ref/Cgar_v2.fa"
export ANC="${BASE}/01_ref/Cgar_v2.fa"          # using ref as ancestor (folded SFS)
export SAMPLES_IND="${BASE}/01_inputs_check/samples.ind"
export SAMPLE_LIST="${BASE}/01_inputs_check/samples_226_pure_gariepinus.txt"
export BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
export CALLABLE_SITES="${BASE}/01_inputs_check/callable_sites.bed.gz"

# ── θπ source (already exists from MODULE_3 het/ROH stream) ───────────────
# These are produced upstream by 02_run_heterozygosity.sh per sample:
#   ANGSD -doSaf 1 -fold 1 -doThetas 1  →  *.thetas.idx
#   thetaStat do_stat -win <W> -step <S> -type 2  →  *.win<W>.step<S>.pestPG
#
# Per-sample SFS-folded θπ (NO n/(n-1) dependence — n=1 valid here).
# This is the right per-sample estimator. region_popstats.c is per-group
# only and undefined at n=1 — do NOT use it as the per-sample input.
#
# v4 scale choice — `win10000.step2000`:
#
# Earlier iterations (v0-v3) anchored θπ to the dosage scrubber's variable-bp
# window grid via nearest-midpoint join. v4 reverses this: θπ uses its own
# native bp-fixed grid at win10000.step2000 (16,500 windows per chrom on
# LG28). Rationale:
#
#   1. The dosage scrubber's sim_mat blow-up at fine scales (1.1 GB at
#      win10000) was the constraint driving v0-v3 toward the coarser grid.
#      θπ doesn't need a sim_mat — |Z|-based envelope detection on per-window
#      population deviations is a 1D scan (n_windows long, not n²). The
#      blow-up doesn't apply.
#   2. PCA structure recovery is between-sample-covariance-driven and
#      averaged across 226 samples, so per-sample θπ noise at finer scales
#      (~14% relative SE at 50-200 sites) doesn't significantly hurt
#      covariance recovery (noise reduces by 1/sqrt(226·pad) ≈ 0.04).
#   3. Per-sample loadings ARE noisier at finer scales (the lines panel and
#      K-means scatter can look uglier), but this is cosmetic; envelope
#      detection is still clean.
#   4. SNP-density variability (0.13-20+ SNPs/kb) doesn't affect θπ
#      stability — θπ depends on CALLABLE sites, not biSNPs. Regions where
#      the dosage scrubber wobbles (sparse SNPs, wide variable-bp windows)
#      are exactly where θπ at finer scale gives clean signal.
#
# Data-layer file size at this scale:
#   per-window matrix:  226 samples × 16,500 windows × 8 bytes = 30 MB raw
#   gzip:               ~6 MB
#
# To revert to the dosage grid (if the per-sample lines panel is too
# noisy on real data): set PESTPG_SCALE="win50000.step10000" and
# uncomment DOSAGE_WIN_BED_DIR usage in STEP_TR00b. One-line change.
export PESTPG_DIR="${BASE}/het_roh/02_heterozygosity/03_theta/multiscale"
export PESTPG_SCALE="win10000.step2000"
export PESTPG_GLOB="*.${PESTPG_SCALE}.pestPG"
export THETA_TSV_DIR="${BASE}/het_roh/02_heterozygosity/05_aggregated"
# After STEP_TR00b runs, the per-sample × per-θπ-window matrix lives at:
#   ${THETA_TSV_DIR}/theta_native.<CHROM>.<SCALE>.tsv.gz
# Columns: sample  chrom  window_idx  start_bp  end_bp  theta_pi  n_sites
# (window_idx is 0-based; this is the θπ-NATIVE grid, NOT the dosage grid)

# ── Dosage scrubber's window grid (kept for fallback / per-candidate Phase 4) ─
# Used by:
#   - The fallback path in STEP_TR00b if PESTPG_SCALE is set to
#     win50000.step10000 (dosage-grid-aligned θπ, the v3 design)
#   - Future Phase 4 per-candidate work where θπ is joined to dosage's
#     candidate intervals via bp-overlap
# Per-chromosome BED produced by 2a_local_pca/STEP09b_dense_window_registry.R.
# Each chromosome has its own subfolder containing windows.bed with columns:
#   chrom  start_bp  end_bp  window_idx  n_snps
export DOSAGE_WIN_BED_DIR="${CODEBASE}/phase_2_discovery/2a_local_pca"

# ── Cursor sync mechanism (v4 design) ─────────────────────────────────────
# Because θπ now uses its own native window grid (not the dosage grid), the
# atlas needs to translate between state.cur (dosage window index) and
# state.cur_thpi (θπ window index) on cursor changes. This is the v1
# lookup-table approach, resurrected for v4 — but cheap (~20 KB of two
# Int32Arrays) and the ONLY way to keep both panels' cursors visually
# synchronized when the grids genuinely differ in resolution.
#
# See cursor_sync_decision_v1.md for the mechanism. Atlas-side cost:
#   - state.thetaPi.dosageToThpi[dosage_idx] → θπ window idx
#   - state.thetaPi.thpiToDosage[thpi_idx]   → dosage window idx
#   - 6 lines added to setCur() to derive state.cur_thpi from state.cur
#   - Alt+→ keyboard handler for native θπ-step navigation
# v1's design notes apply directly; no new architectural decision needed.

# ── Window grid metadata (informational) ──────────────────────────────────
# Documents what scale we're at; STEP_TR01 emits these to the JSON layer.
export WIN_BP=10000
export STEP_BP=2000
export SCALE_LABEL="${PESTPG_SCALE}"   # e.g. "win10000.step2000"

# ── Local-PCA parameters (mirror dosage local PCA from sibling 2a) ────────
# Match the dosage scrubber's neighbourhood size so cross-stream candidate
# overlap is computed on comparable window scales.
export LOCAL_PCA_PAD=1                          # ±1 window neighbourhood for local PCA
export LOCAL_PCA_NPC=2                          # PC1 + PC2 sufficient for sim_mat
export SIM_MAT_METRIC="abs_cosine"              # |cos(loadings_i, loadings_j)|

# ── Envelope detection thresholds (mirror dosage envelopes) ───────────────
# Read the same values 2a/2d use; do NOT diverge without explicit reason.
export CONCORD_THR=0.85                         # L2 envelope concord threshold
export MERGE_THR=0.85                           # L2 → L1 merge threshold
export SILHOUETTE_MIN=0.45                      # L2 silhouette filter
export ENV_MIN_WINDOWS=5                        # min windows for L2 to be retained

# ── Compute / SLURM ───────────────────────────────────────────────────────
export RSCRIPT="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"
export MAMBA_ENV="assembly"
export SLURM_ACCOUNT="lt200308"
export SLURM_PARTITION="compute"
export SLURM_TIME="04:00:00"
export SLURM_MEM="32G"
export SLURM_CPUS=4
# 28 chromosomes — array goes 1..28, one chrom per task. Same pattern as
# MODULE_5A2's STEP09b/STEP10v2 chromosome arrays.
export SLURM_ARRAY="1-28"

# ── Chromosomes ───────────────────────────────────────────────────────────
# Pure C. gariepinus karyotype = 28 linkage groups (LG01..LG28).
export CHROM_LIST=("C_gar_LG01" "C_gar_LG02" "C_gar_LG03" "C_gar_LG04"
                   "C_gar_LG05" "C_gar_LG06" "C_gar_LG07" "C_gar_LG08"
                   "C_gar_LG09" "C_gar_LG10" "C_gar_LG11" "C_gar_LG12"
                   "C_gar_LG13" "C_gar_LG14" "C_gar_LG15" "C_gar_LG16"
                   "C_gar_LG17" "C_gar_LG18" "C_gar_LG19" "C_gar_LG20"
                   "C_gar_LG21" "C_gar_LG22" "C_gar_LG23" "C_gar_LG24"
                   "C_gar_LG25" "C_gar_LG26" "C_gar_LG27" "C_gar_LG28")

# ── Output JSON path (consumed by atlas page 12) ──────────────────────────
# STEP_TR_B writes one consolidated JSON per chromosome at:
#   ${JSON_OUT_DIR}/${CHROM}/${CHROM}_phase2_theta.json
# The file carries four atlas layers in one bundle: theta_pi_per_window,
# theta_pi_local_pca, theta_pi_envelopes, and a tracks contribution.
# Atlas precomp loader recognizes these layer keys at
# pca_scrubber_v4.html detectSchemaAndLayers() lines 24316–24327.
# No atlas-side change needed for layer detection.
export JSON_OUT_DIR="${OUTROOT}/json_out"

# Legacy per-layer dirs (used by the archived STEP_TR01 long-format
# emitter; kept for the popstats page and other non-page-12 consumers
# that still expect three separate files).
export OUT_PER_WINDOW_DIR="${OUTROOT}/01_per_window"
export OUT_LOCAL_PCA_DIR="${OUTROOT}/02_local_pca"
export OUT_ENVELOPES_DIR="${OUTROOT}/03_envelopes"

# ── Logging ───────────────────────────────────────────────────────────────
export LOG_DIR="${OUTROOT}/logs"

# ── Schema version (matches theta_pi_data_layer_spec_v0.md) ───────────────
export THETA_JSON_SCHEMA_VERSION=1

# ── Sanity: create output dirs if missing ─────────────────────────────────
mkdir -p "${JSON_OUT_DIR}" "${OUT_PER_WINDOW_DIR}" "${OUT_LOCAL_PCA_DIR}" \
         "${OUT_ENVELOPES_DIR}" "${LOG_DIR}"

# ── Print config on demand ────────────────────────────────────────────────
theta_config_print() {
  cat <<EOF
=== 2f_theta_discovery config ===
BASE             : ${BASE}
THETA_DISC_DIR   : ${THETA_DISC_DIR}
OUTROOT          : ${OUTROOT}
PESTPG_DIR       : ${PESTPG_DIR}
SCALE_LABEL      : ${SCALE_LABEL}
WIN_BP / STEP_BP : ${WIN_BP} / ${STEP_BP}
LOCAL_PCA_PAD    : ${LOCAL_PCA_PAD}
LOCAL_PCA_NPC    : ${LOCAL_PCA_NPC}
CONCORD_THR      : ${CONCORD_THR}
MERGE_THR        : ${MERGE_THR}
SILHOUETTE_MIN   : ${SILHOUETTE_MIN}
RSCRIPT          : ${RSCRIPT}
N_CHROMS         : ${#CHROM_LIST[@]}
SCHEMA_VERSION   : ${THETA_JSON_SCHEMA_VERSION}
=================================
EOF
}

# Optional: invoke ./00_theta_config.sh print  to dump config
if [[ "${1:-}" == "print" ]]; then
  theta_config_print
fi
