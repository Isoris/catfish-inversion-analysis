#!/usr/bin/env bash
# =============================================================================
# STEP_ZO_HtoM_run_one_chrom.sh  (v10, 2026-05-16, harmonized layout)
#
# Run steps L1_detect → L1_plot → L2_detect → L2_plot → atlas_json
# (formerly ZO_H/I/J/K + ZO_M) for ONE chromosome with canonical defaults.
# Pulls all paths from the harmonized config (00_inversion_config.sh v10):
#   - L1/L2 detect+plot scripts live in _shared/ (STEP_05/06/07/08_*.R)
#   - Plot stages honor MAKE_L1_PLOTS / MAKE_L2_PLOTS (default true)
#   - Precomp lives in 04_precomp/ (was 04_precomp_v2/precomp/ in v9)
#
# Usage:
#   bash STEP_ZO_HtoM_run_one_chrom.sh C_gar_LG28
#
# Override paths via env (rarely needed):
#   CONFIG, SAMPLE_META, MAKE_L1_PLOTS, MAKE_L2_PLOTS
# =============================================================================

set -euo pipefail

CHR="${1:?Provide chromosome label as argv1, e.g. C_gar_LG28}"

CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_localpca_v7/_shared/00_inversion_config.sh}"
if [[ -f "${CONFIG}" ]]; then
  set -a; source "${CONFIG}"; set +a
fi

# Defaults if config wasn't sourced (lets the script run standalone).
BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
SCRATCH_ROOT="${SCRATCH_ROOT:-${BASE}/inversion_localpca_v7}"
SCRIPT_DIR="${SCRIPT_DIR:-${BASE}/catfish-inversion-analysis/catfish-inversion-analysis/local_PCA_MDS_z}"
SCRIPT_DIR_SHARED="${SCRIPT_DIR_SHARED:-${BASE}/catfish-inversion-analysis/catfish-inversion-analysis/_shared}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"

PRECOMP_DIR="${PATH1_PRECOMP:-${SCRATCH_ROOT}/local_PCA_MDS_z/04_precomp}"
L1_DIR="${PATH1_L1:-${SCRATCH_ROOT}/local_PCA_MDS_z/05_L1}"
L1_PLOT_DIR="${PATH1_L1_PLOTS:-${SCRATCH_ROOT}/local_PCA_MDS_z/06_L1_plots}"
L2_DIR="${PATH1_L2:-${SCRATCH_ROOT}/local_PCA_MDS_z/07_L2}"
L2_PLOT_DIR="${PATH1_L2_PLOTS:-${SCRATCH_ROOT}/local_PCA_MDS_z/08_L2_plots}"
JSON_DIR="${PATH1_JSON:-${SCRATCH_ROOT}/local_PCA_MDS_z/09_atlas_json}"
SAMPLE_META="${SAMPLE_META:-${SHARED_DIR:-${SCRATCH_ROOT}/_shared}/sample_metadata.tsv}"
MAKE_L1_PLOTS="${MAKE_L1_PLOTS:-true}"
MAKE_L2_PLOTS="${MAKE_L2_PLOTS:-true}"
export MAKE_L1_PLOTS MAKE_L2_PLOTS

mkdir -p "${L1_DIR}" "${L1_PLOT_DIR}" "${L2_DIR}" "${L2_PLOT_DIR}" "${JSON_DIR}"

echo "============================================================"
echo "  Pipeline L1_detect → L1_plot → L2_detect → L2_plot → atlas"
echo "  chromosome      = ${CHR}"
echo "  precomp_dir     = ${PRECOMP_DIR}"
echo "  sample_meta     = ${SAMPLE_META}"
echo "  MAKE_L1_PLOTS   = ${MAKE_L1_PLOTS}"
echo "  MAKE_L2_PLOTS   = ${MAKE_L2_PLOTS}"
echo "============================================================"

# ── 05 L1 detect (default --nn 80) ──────────────────────────────────────────
echo ""
echo "=== _shared/STEP_05_L1_detect.R (nn80) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR_SHARED}/STEP_05_L1_detect.R" \
  --precomp_dir "${PRECOMP_DIR}" \
  --chr         "${CHR}" \
  --outdir      "${L1_DIR}" \
  --boundary_scan TRUE \
  --boundary_validator_mode grow \
  --boundary_W 5 --boundary_offset 5 --boundary_min_dist 30

# ── 06 L1 plot (gated by MAKE_L1_PLOTS) ─────────────────────────────────────
echo ""
echo "=== _shared/STEP_06_L1_plot.R (nn80) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR_SHARED}/STEP_06_L1_plot.R" \
  --precomp_dir     "${PRECOMP_DIR}" \
  --L1_dir          "${L1_DIR}" \
  --chr             "${CHR}" \
  --outdir          "${L1_PLOT_DIR}" \
  --toggle_L1       yes \
  --boundary_filter stable

# ── 07 L2 detect (default --nn 40) ──────────────────────────────────────────
echo ""
echo "=== _shared/STEP_07_L2_detect.R (nn40) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR_SHARED}/STEP_07_L2_detect.R" \
  --precomp_dir "${PRECOMP_DIR}" \
  --L1_dir      "${L1_DIR}" \
  --chr         "${CHR}" \
  --outdir      "${L2_DIR}" \
  --boundary_scan TRUE \
  --boundary_validator_mode grow \
  --quadrant_validator yes \
  --weak_demote_score 0 \
  --quad_rescue_max_grow_z 1.5 \
  --quad_demote_on_fail yes \
  --quad_demote_drift_floor -1.0

# ── 08 L2 plot (gated by MAKE_L2_PLOTS) ─────────────────────────────────────
echo ""
echo "=== _shared/STEP_08_L2_plot.R (nn80 chrom + nn40 inside) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR_SHARED}/STEP_08_L2_plot.R" \
  --precomp_dir     "${PRECOMP_DIR}" \
  --L1_dir          "${L1_DIR}" \
  --L2_dir          "${L2_DIR}" \
  --chr             "${CHR}" \
  --outdir          "${L2_PLOT_DIR}" \
  --boundary_filter stable

# ── 09 atlas JSON (needs sample_metadata.tsv from ZO_L; one-time genome-wide) ─
echo ""
echo "=== STEP_ZO_M_export_atlas_json.R ==="
if [[ ! -f "${SAMPLE_META}" ]]; then
  echo "[WARN] sample_metadata.tsv not found at ${SAMPLE_META}"
  echo "[WARN] Run STEP_ZO_L_build_sample_metadata.R once to create it (genome-wide)."
  echo "[WARN] Skipping JSON export for ${CHR}."
  exit 0
fi
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_ZO_M_export_atlas_json.R" \
  --precomp_dir     "${PRECOMP_DIR}" \
  --L1_dir          "${L1_DIR}" \
  --L2_dir          "${L2_DIR}" \
  --chr             "${CHR}" \
  --sample_metadata "${SAMPLE_META}" \
  --out             "${JSON_DIR}/${CHR}.atlas.json"

echo ""
echo "============================================================"
echo "  All steps complete for ${CHR}"
echo "  L1 envelopes : ${L1_DIR}/${CHR}.L1_envelopes.tsv"
echo "  L2 envelopes : ${L2_DIR}/${CHR}.L2_envelopes.tsv"
echo "  L1 plot      : ${L1_PLOT_DIR}/${CHR}.L1_overlay.pdf  (MAKE_L1_PLOTS=${MAKE_L1_PLOTS})"
echo "  L2 plot      : ${L2_PLOT_DIR}/${CHR}.L2_overlay.pdf  (MAKE_L2_PLOTS=${MAKE_L2_PLOTS})"
echo "  Atlas JSON   : ${JSON_DIR}/${CHR}.atlas.json"
echo "============================================================"
