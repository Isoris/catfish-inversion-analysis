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
# STEP_ZO_D_LAUNCH_local_pca_merge.slurm  (v10, 2026-05-16, harmonized layout)
#
# STAGE 2 (merge): single short job. Reads the per-chr <chr>.window_pca_tmp.rds
# from 01_local_pca/, assigns globally unique window_ids, writes finalized
# per-chr <chr>.window_pca.rds back into 01_local_pca/, and writes the
# windows_master.tsv.gz + windows_master_summary.tsv into 02_dense_registry/.
#
# v10: no <outdir>/tmp/ subdir anymore — the stage-1 *_tmp.rds files coexist
# with the finals in 01_local_pca/. The merge script takes --outdir for the
# per-chrom dir and --registry_dir for the registry artifacts.
#
# Submit AFTER STEP_ZO_C stage1 array job completes:
#   sbatch --dependency=afterok:<JOBID_C> STEP_ZO_D_LAUNCH_local_pca_merge.slurm
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="${CONFIG:-${SCRATCH:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}/inversion_localpca_v7/_shared/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a; source "${CONFIG}"; set +a
mkdir -p logs

SCRIPT_DIR="${SCRIPT_DIR:-${SCRIPT_DIR}}"
# Per-chrom dir = 01_local_pca slot; registry dir = 02_dense_registry slot.
OUTDIR="${OUTDIR:-${PATH1_LOCAL_PCA}}"
REGISTRY_DIR="${REGISTRY_DIR:-${PATH1_REGISTRY}}"
mkdir -p "${OUTDIR}" "${REGISTRY_DIR}"

cat <<EOF
================================================================
  STEP_ZO_D local_pca_merge (v10 harmonized layout)
  Stage-1 _tmp.rds source : ${OUTDIR}
  Per-chr finals          : ${OUTDIR}
  Registry artifacts      : ${REGISTRY_DIR}
================================================================
EOF

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_ZO_D_local_pca_merge.R" \
  --outdir       "${OUTDIR}" \
  --registry_dir "${REGISTRY_DIR}"

echo "[DONE] STEP_ZO_D — $(date)"
