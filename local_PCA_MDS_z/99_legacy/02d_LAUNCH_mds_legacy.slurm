#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_step10
#SBATCH -o logs/inv_step10.%j.out
#SBATCH -e logs/inv_step10.%j.err

set -euo pipefail

# =============================================================================
# LAUNCH_STEP10_window_mds_outliers.slurm  (v8.5)
#
# Legacy single-shot MDS outlier detection (superseded by STEP10v2 parallel).
# Kept for backward compatibility / comparison runs.
#
# Usage:
#   sbatch LAUNCH_STEP10_window_mds_outliers.slurm
# =============================================================================

source ~/.bashrc
mamba activate assembly

CONFIG="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/00_inversion_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

inv_init_dirs

A2_STEPS="${DISCOVERY2DIR}/steps"

"${RSCRIPT_BIN}" "${A2_STEPS}/STEP_B01_legacy_mds_perchr.R" \
  "${LPCA_DIR}" \
  "${MDS_PREFIX}" \
  "${NPC}" "${MDS_DIMS_DEFAULT}" "${Z_THRESH_DEFAULT}" \
  "${GAP_BP_CONTROL_DEFAULT}" "${MIN_WINDOWS_CONTROL_DEFAULT}"

echo "[DONE] STEP10 — $(date)"
