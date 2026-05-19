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
# STEP_ZO_E_LAUNCH_mds_compute.slurm  (v10, 2026-05-16, harmonized layout)
#
# STAGE 1 (compute): one SLURM array task per FOCAL chromosome. Reads
# per-chr <chr>.window_pca.rds (from 01_local_pca/), computes the lostruct
# distance matrix, runs cmdscale, writes per-focal-chr MDS RDS to the
# 03_mds/ slot. Heaviest step in the pipeline (~12h walltime per chromosome
# in chunked_2x mode; ~4 min in chromosome mode).
#
# v10: per-chrom MDS RDS now lands directly in 03_mds/ — no tmp/ subdir.
# Downstream ZO_G discovers it via $PATH1_MDS or --mds_dir.
#
# Submit:
#   sbatch --array=0-27 STEP_ZO_E_LAUNCH_mds_compute.slurm chrom.list
#
# Output (in $PATH1_MDS = 03_mds/):
#   <focal_chr>.mds_perchr.rds  — per-focal-chr MDS coords
#   <focal_chr>.metadata.tsv    — per-focal-chr one-row metadata
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="${CONFIG:-${SCRATCH:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}/inversion_localpca_v7/_shared/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a; source "${CONFIG}"; set +a
mkdir -p logs

# Harmonized v10 paths from config.
RDS_DIR="${RDS_DIR:-${PATH1_LOCAL_PCA}}"   # was 02_dense_registry; in v10 the
                                            # finalized per-chr .window_pca.rds
                                            # files stay in 01_local_pca/.
OUTDIR="${OUTDIR:-${PATH1_MDS}}"            # 03_mds/ slot, no tmp/ subdir
OUTPREFIX="${OUTPREFIX:-${MDS_PREFIX_BASENAME:-inversion_localpca}}"

# Defaults (canonical from history line 19553-area)
MDS_MODE="${MDS_MODE:-chunked_2x}"
NPC="${NPC:-4}"
MDS_DIMS="${MDS_DIMS:-${MDS_DIMS_DEFAULT:-5}}"
Z_THRESH="${Z_THRESH:-${Z_THRESH_DEFAULT:-3.0}}"

CHROM_LIST="${1:?ERROR: provide chrom.list as argv1}"

mkdir -p "${OUTDIR}"

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
