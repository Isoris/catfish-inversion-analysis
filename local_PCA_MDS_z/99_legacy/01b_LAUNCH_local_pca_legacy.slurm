#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_step09
#SBATCH -o logs/inv_step09.%A_%a.out
#SBATCH -e logs/inv_step09.%A_%a.err

set -euo pipefail

# =============================================================================
# LAUNCH_STEP09_local_pca_windows_by_chr.slurm  (v8.5)
#
# Per-chromosome local PCA (non-overlapping windows, legacy mode).
#
# Usage:
#   sbatch --array=0-27%8 LAUNCH_STEP09_local_pca_windows_by_chr.slurm chrom.list
# =============================================================================

source ~/.bashrc
mamba activate assembly

CONFIG="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/00_inversion_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

inv_init_dirs

CHROM_LIST="${1:?ERROR: provide chrom list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))

CHROM="$(sed -n "${LINE}p" "${CHROM_LIST}")"
if [[ -z "${CHROM}" ]]; then
  echo "[ERROR] Could not read chromosome from ${CHROM_LIST} line ${LINE}" >&2
  exit 1
fi

A2_STEPS="${DISCOVERY2DIR}/steps"

SITES="${DOSAGE_DIR}/${CHROM}.sites.tsv.gz"
DOSAGE="${DOSAGE_DIR}/${CHROM}.dosage.tsv.gz"
OUTPREFIX="${LPCA_DIR}/STEP09_${CHROM}"

[[ -x "${RSCRIPT_BIN}" ]] || { echo "[ERROR] Rscript not executable: ${RSCRIPT_BIN}" >&2; exit 1; }
[[ -s "${SITES}" ]] || { echo "[ERROR] Missing ${SITES}" >&2; exit 1; }
[[ -s "${DOSAGE}" ]] || { echo "[ERROR] Missing ${DOSAGE}" >&2; exit 1; }

if [[ -s "${OUTPREFIX}.window_meta.tsv.gz" && -s "${OUTPREFIX}.window_pca.rds" && -s "${OUTPREFIX}.window_pca.tsv.gz" ]]; then
  echo "[INFO] STEP09 outputs already exist for ${CHROM}, skipping"
  exit 0
fi

echo "[INFO] CHROM=${CHROM}"
echo "[INFO] Using Rscript: ${RSCRIPT_BIN}"
"${RSCRIPT_BIN}" -e 'cat("R=", R.version.string, "\n"); cat("Rscript=", Sys.which("Rscript"), "\n"); suppressPackageStartupMessages(library(RSpectra)); cat("RSpectra OK\n")'

"${RSCRIPT_BIN}" "${A2_STEPS}/STEP_A02_local_pca_windows_by_chr.R" \
  "${SITES}" \
  "${DOSAGE}" \
  "${OUTPREFIX}" \
  "${NPC}" \
  "${WINSIZE}"

echo "[DONE] ${CHROM}"
