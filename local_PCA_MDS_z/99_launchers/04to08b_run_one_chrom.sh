#!/usr/bin/env bash
# =============================================================================
# 04to08b_run_one_chrom.sh
#
# Run steps 04 → 05 → 06 → 07 → 08b for ONE chromosome with canonical defaults.
# All canonical parameters are baked in (history lines 19548 → 19650). The
# only required arg is the chromosome label.
#
# Usage:
#   bash 04to08b_run_one_chrom.sh C_gar_LG28
#
# Override paths via env vars (rarely needed):
#   SCRATCH_ROOT, SCRIPT_DIR, SAMPLE_META
#
# Each chromosome takes ~5-10 min interactive (no SLURM array needed; this
# script is meant to be called inside a single allocation or interactively).
# =============================================================================

set -euo pipefail

CHR="${1:?Provide chromosome label as argv1, e.g. C_gar_LG28}"

# Path resolution (defaults match the canonical scratch tree)
BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
SCRATCH_ROOT="${SCRATCH_ROOT:-${BASE}/inversion_localpca_v8}"
SCRIPT_DIR="${SCRIPT_DIR:-${BASE}/inversion-popgen-toolkit/local_PCA_z}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"

PRECOMP_DIR="${SCRATCH_ROOT}/path_localpca_zblocks/04_precomp/precomp"
L1_DIR="${SCRATCH_ROOT}/path_localpca_zblocks/05_L1"
L1_PLOT_DIR="${SCRATCH_ROOT}/path_localpca_zblocks/06_L1_plots"
L2_DIR="${SCRATCH_ROOT}/path_localpca_zblocks/07_L2"
L2_PLOT_DIR="${SCRATCH_ROOT}/path_localpca_zblocks/08_L2_plots"
JSON_DIR="${SCRATCH_ROOT}/path_localpca_zblocks/09_atlas_json"
SAMPLE_META="${SAMPLE_META:-${SCRATCH_ROOT}/_shared/sample_metadata.tsv}"

mkdir -p "${L1_DIR}" "${L1_PLOT_DIR}" "${L2_DIR}" "${L2_PLOT_DIR}" "${JSON_DIR}"

echo "============================================================"
echo "  Pipeline 04→05→06→07→08b   chromosome=${CHR}"
echo "  precomp_dir : ${PRECOMP_DIR}"
echo "  sample_meta : ${SAMPLE_META}"
echo "============================================================"

# ── Step 04: detect L1 (default --nn 80) ─────────────────────────────────────
echo ""
echo "=== 04 detect_L1 (nn80) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/04_detect_L1/04_detect_L1_localpca_zblocks.R" \
  --precomp_dir "${PRECOMP_DIR}" \
  --chr         "${CHR}" \
  --outdir      "${L1_DIR}" \
  --boundary_scan TRUE \
  --boundary_validator_mode grow \
  --boundary_W 5 --boundary_offset 5 --boundary_min_dist 30

# ── Step 05: plot L1 ─────────────────────────────────────────────────────────
echo ""
echo "=== 05 plot_L1 (nn80) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/05_plot_L1/05_plot_L1_localpca_zblocks.R" \
  --precomp_dir     "${PRECOMP_DIR}" \
  --L1_dir          "${L1_DIR}" \
  --chr             "${CHR}" \
  --outdir          "${L1_PLOT_DIR}" \
  --toggle_L1       yes \
  --boundary_filter stable

# ── Step 06: detect L2 (default --nn 40) ─────────────────────────────────────
echo ""
echo "=== 06 detect_L2 (nn40) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/06_detect_L2/06_detect_L2_localpca_zblocks.R" \
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

# ── Step 07: plot L2 ─────────────────────────────────────────────────────────
echo ""
echo "=== 07 plot_L2 (nn80 chrom + nn40 inside) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/07_plot_L2/07_plot_L2_localpca_zblocks.R" \
  --precomp_dir     "${PRECOMP_DIR}" \
  --L1_dir          "${L1_DIR}" \
  --L2_dir          "${L2_DIR}" \
  --chr             "${CHR}" \
  --outdir          "${L2_PLOT_DIR}" \
  --boundary_filter stable

# ── Step 08b: export atlas JSON (assumes 08a has been run once already) ──────
echo ""
echo "=== 08b export_atlas_json ==="
if [[ ! -f "${SAMPLE_META}" ]]; then
  echo "[WARN] sample_metadata.tsv not found at ${SAMPLE_META}"
  echo "[WARN] Run 08a_build_sample_metadata.R once to create it (genome-wide)."
  echo "[WARN] Skipping JSON export for ${CHR}."
  exit 0
fi
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/08_atlas_json/08b_export_atlas_json_localpca_zblocks.R" \
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
echo "  L1 plot      : ${L1_PLOT_DIR}/${CHR}.L1_overlay.pdf"
echo "  L2 plot      : ${L2_PLOT_DIR}/${CHR}.L2_overlay.pdf"
echo "  Atlas JSON   : ${JSON_DIR}/${CHR}.atlas.json"
echo "============================================================"
