#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -t 1:00:00
#SBATCH -J inv_02b_merge
#SBATCH -o logs/inv_02b_merge.%j.out
#SBATCH -e logs/inv_02b_merge.%j.err

# =============================================================================
# 02b_LAUNCH_mds_merge.slurm  (v8.5 / consolidated layout v1.0)
#
# STAGE 2 (merge): single short job. Stitches per-focal-chr MDS RDS files in
# tmp/ into the final inversion_localpca.mds.rds (with $per_chr structure
# every downstream script expects), plus candidate-region detection.
#
# Submit AFTER 02a stage1 array job completes:
#   sbatch --dependency=afterok:<JOBID_02a> 02b_LAUNCH_mds_merge.slurm
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
TMPDIR_IN="${TMPDIR_IN:-${SCRATCH_ROOT}/path_localpca_zblocks/03_mds/tmp}"
OUTPREFIX="${OUTPREFIX:-${SCRATCH_ROOT}/path_localpca_zblocks/03_mds/inversion_localpca}"

MDS_MODE="${MDS_MODE:-chunked_2x}"
GAP_BP="${GAP_BP:-500000}"
MIN_WINDOWS="${MIN_WINDOWS:-3}"

cat <<EOF
================================================================
  02b mds_merge
  tmp dir  : ${TMPDIR_IN}
  outprefix: ${OUTPREFIX}
  mode=${MDS_MODE}  gap_bp=${GAP_BP}  min_windows=${MIN_WINDOWS}
================================================================
EOF

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_ZO_02b_mds_merge.R" \
  --tmpdir      "${TMPDIR_IN}" \
  --outprefix   "${OUTPREFIX}" \
  --mds_mode    "${MDS_MODE}" \
  --gap_bp      "${GAP_BP}" \
  --min_windows "${MIN_WINDOWS}"

echo "[DONE] 02b — $(date)"
