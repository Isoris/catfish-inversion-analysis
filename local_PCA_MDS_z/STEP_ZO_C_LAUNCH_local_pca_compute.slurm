#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH -t 6:00:00
#SBATCH -J inv_C_compute
#SBATCH -o logs/inv_C_compute.%A_%a.out
#SBATCH -e logs/inv_C_compute.%A_%a.err

# =============================================================================
# STEP_ZO_C_LAUNCH_local_pca_compute.slurm   (v8.5 / consolidated layout v1.0)
#
# STAGE 1 (compute): one SLURM array task per chromosome. Sliding-window local
# PCA on dosage matrix. Writes per-chrom RDS with window_id = NA. The merge
# step (D) assigns global window_ids.
#
# Submit:
#   sbatch --array=0-27 STEP_ZO_C_LAUNCH_local_pca_compute.slurm chrom.list
#
# Output (path_localpca_zblocks/01_local_pca/tmp/):
#   <chr>.window_pca_tmp.rds, <chr>.chr_meta.tsv.gz, <chr>.pca_table.tsv.gz,
#   <chr>.summary.tsv
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="${CONFIG:-${SCRATCH:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}/inversion_localpca_v8/_shared/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a; source "${CONFIG}"; set +a
mkdir -p logs

SCRIPT_DIR="${SCRIPT_DIR:-${BASE}/inversion-popgen-toolkit/local_PCA_z}"
SCRATCH_ROOT="${SCRATCH_ROOT:-${BASE}/inversion_localpca_v8}"
DOSAGE_DIR="${DOSAGE_DIR:-${SCRATCH_ROOT}/02_dosage_sites}"
OUTDIR="${OUTDIR:-${SCRATCH_ROOT}/path_localpca_zblocks/01_local_pca}"

CHROM_LIST="${1:?ERROR: provide chrom.list as argv1}"
WINSIZE="${2:-100}"
STEP_SIZE="${3:-20}"
NPC_ARG="${4:-4}"

mkdir -p "${OUTDIR}/tmp"

IDX="${SLURM_ARRAY_TASK_ID:-0}"
LINE=$((IDX + 1))
CHROM="$(sed -n "${LINE}p" "${CHROM_LIST}")"
[[ -n "${CHROM}" ]] || { echo "[ERROR] cannot read chr from ${CHROM_LIST} line ${LINE}" >&2; exit 1; }

cat <<EOF
================================================================
  STEP_ZO_C local_pca_compute  —  ${CHROM}
  Array task : ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_MAX:-?}
  Dosage dir : ${DOSAGE_DIR}
  Out dir    : ${OUTDIR}
  winsize=${WINSIZE}  step=${STEP_SIZE}  npc=${NPC_ARG}
================================================================
EOF

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_ZO_C_local_pca_compute.R" \
  --dosage_dir "${DOSAGE_DIR}" \
  --outdir     "${OUTDIR}" \
  --chrom      "${CHROM}" \
  --winsize    "${WINSIZE}" \
  --step       "${STEP_SIZE}" \
  --npc        "${NPC_ARG}"

echo "[DONE] STEP_ZO_C — ${CHROM} — $(date)"
