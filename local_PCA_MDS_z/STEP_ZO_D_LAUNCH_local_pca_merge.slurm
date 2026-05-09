#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 1:00:00
#SBATCH -J inv_D_merge
#SBATCH -o logs/inv_D_merge.%j.out
#SBATCH -e logs/inv_D_merge.%j.err

# =============================================================================
# STEP_ZO_D_LAUNCH_local_pca_merge.slurm  (v8.5 / consolidated layout v1.0)
#
# STAGE 2 (merge): single short job. Reads the 28 per-chr RDS in tmp/, assigns
# globally unique window_ids, writes finalized per-chr <chr>.window_pca.rds +
# windows_master.tsv.gz (the master registry).
#
# Submit AFTER STEP_ZO_C stage1 array job completes:
#   sbatch --dependency=afterok:<JOBID_C> STEP_ZO_D_LAUNCH_local_pca_merge.slurm
#
# Output (path_localpca_zblocks/02_dense_registry/):
#   windows_master.tsv.gz, windows_master_summary.tsv,
#   <chr>.window_pca.rds, <chr>.window_pca.tsv.gz
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

# STEP_ZO_C wrote tmp/ inside 01_local_pca/. The merge READS from there and WRITES
# the final registry to 02_dense_registry/.
INDIR="${INDIR:-${SCRATCH_ROOT}/path_localpca_zblocks/01_local_pca}"
OUTDIR="${OUTDIR:-${SCRATCH_ROOT}/path_localpca_zblocks/02_dense_registry}"
mkdir -p "${OUTDIR}"

# The merge script expects --outdir to be the same dir that contains tmp/. So
# we symlink/copy the tmp from 01_local_pca/ into 02_dense_registry/ first.
if [[ ! -d "${OUTDIR}/tmp" ]]; then
  ln -sfn "${INDIR}/tmp" "${OUTDIR}/tmp"
fi

cat <<EOF
================================================================
  STEP_ZO_D local_pca_merge
  Reading tmp from : ${INDIR}/tmp
  Writing to       : ${OUTDIR}
================================================================
EOF

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_ZO_D_local_pca_merge.R" \
  --outdir "${OUTDIR}"

echo "[DONE] STEP_ZO_D — $(date)"
