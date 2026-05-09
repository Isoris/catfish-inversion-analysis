#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -J inv_E_compute
#SBATCH -o logs/inv_E_compute.%A_%a.out
#SBATCH -e logs/inv_E_compute.%A_%a.err

# =============================================================================
# STEP_ZO_E_LAUNCH_mds_compute.slurm  (v8.5 / consolidated layout v1.0)
#
# STAGE 1 (compute): one SLURM array task per FOCAL chromosome. Reads ALL
# per-chr window_pca.rds (needed for chunked background sampling), computes
# lostruct distance matrix, runs cmdscale, writes per-focal-chr MDS RDS to
# tmp/. Heaviest step in the pipeline (~12h walltime per chromosome).
#
# Submit:
#   sbatch --array=0-27 STEP_ZO_E_LAUNCH_mds_compute.slurm chrom.list
#
# Output (path_localpca_zblocks/03_mds/tmp/):
#   <focal_chr>.mds_perchr.rds  — per-focal-chr MDS coords
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
RDS_DIR="${RDS_DIR:-${SCRATCH_ROOT}/path_localpca_zblocks/02_dense_registry}"
OUTDIR="${OUTDIR:-${SCRATCH_ROOT}/path_localpca_zblocks/03_mds}"
OUTPREFIX="${OUTPREFIX:-inversion_localpca}"

# Defaults (canonical from history line 19553-area)
MDS_MODE="${MDS_MODE:-chunked_2x}"
NPC="${NPC:-4}"
MDS_DIMS="${MDS_DIMS:-20}"
Z_THRESH="${Z_THRESH:-3.0}"

CHROM_LIST="${1:?ERROR: provide chrom.list as argv1}"

mkdir -p "${OUTDIR}/tmp"

IDX="${SLURM_ARRAY_TASK_ID:-0}"
LINE=$((IDX + 1))
FOCAL_CHR="$(sed -n "${LINE}p" "${CHROM_LIST}")"
[[ -n "${FOCAL_CHR}" ]] || { echo "[ERROR] cannot read chr" >&2; exit 1; }

cat <<EOF
================================================================
  STEP_ZO_E mds_compute  —  focal=${FOCAL_CHR}
  Array task : ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_MAX:-?}
  RDS dir    : ${RDS_DIR}
  Out dir    : ${OUTDIR}
  mode=${MDS_MODE}  npc=${NPC}  mds_dims=${MDS_DIMS}  z_thresh=${Z_THRESH}
================================================================
EOF

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_ZO_E_mds_compute.R" \
  --rds_dir   "${RDS_DIR}" \
  --outdir    "${OUTDIR}" \
  --outprefix "${OUTPREFIX}" \
  --focal_chr "${FOCAL_CHR}" \
  --mds_mode  "${MDS_MODE}" \
  --npc       "${NPC}" \
  --mds_dims  "${MDS_DIMS}" \
  --z_thresh  "${Z_THRESH}"

echo "[DONE] STEP_ZO_E — focal=${FOCAL_CHR} — $(date)"
