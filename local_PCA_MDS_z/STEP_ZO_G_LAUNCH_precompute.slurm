#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 2:00:00
#SBATCH -J inv_G_precomp
#SBATCH -o logs/inv_G_precomp.%j.out
#SBATCH -e logs/inv_G_precomp.%j.err

# =============================================================================
# STEP_ZO_G_LAUNCH_precompute.slurm  (v10, 2026-05-16, harmonized layout)
#
# Per-chromosome precompute: builds <chr>.precomp.rds + the suite of
# NN-smoothed similarity matrices (<chr>.sim_mat_nn{0,...,320}.rds).
# Parallelised across chromosomes via mclapply. NO whole-genome
# combined RDS is written.
#
# Submit AFTER all ZO_E per-chrom MDS array tasks complete. ZO_G reads
# $PATH1_MDS/<chr>.mds_perchr.rds directly (was <OUTDIR>/tmp/ in v9).
#
# Output ($PATH1_PRECOMP = 04_precomp/):
#   <chr>.precomp.rds
#   sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
#   precomp_summary.tsv
#
# Sweet spot — 28 cores (1:1 with chromosomes), 237G memory cap. Real usage <80G.
# Tune the NN scales saved via:  NN_SIM_SCALES="40,80,160,320" sbatch ...
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="${CONFIG:-${SCRATCH:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}/inversion_localpca_v7/_shared/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a; source "${CONFIG}"; set +a
mkdir -p logs

# Harmonized v10 dirs from config.
MDS_DIR="${MDS_DIR:-${PATH1_MDS}}"
OUTDIR="${OUTDIR:-${PATH1_PRECOMP}}"
DOSAGE_DIR="${DOSAGE_DIR:-${SCRATCH_ROOT}/02_dosage_sites}"
# MDS_PREFIX is a legacy positional ZO_G still accepts; dirname() falls
# back to it only if --mds_dir resolution fails.
MDS_PREFIX="${MDS_PREFIX:-${PATH1_MDS_PREFIX}}"

# Sanity-check that ZO_E array tasks have written outputs to PATH1_MDS.
if [[ ! -d "${MDS_DIR}" ]] || [[ -z $(find "${MDS_DIR}" -maxdepth 1 -name "*.mds_perchr.rds" -print -quit) ]]; then
  echo "[G] missing ${MDS_DIR}/<chr>.mds_perchr.rds — run ZO_E first" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

cat <<EOF
================================================================
  STEP_ZO_G precompute (v10 harmonized layout)
  MDS dir         : ${MDS_DIR}
  Precomp out dir : ${OUTDIR}
  NN_SIM_SCALES   : ${NN_SIM_SCALES:-40,80,160,320}
  Cores           : ${SLURM_CPUS_PER_TASK:-1}
================================================================
EOF

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_ZO_G_precompute.R" \
  "${MDS_PREFIX}" \
  "${OUTDIR}" \
  --mds_dir "${MDS_DIR}"

cat <<EOF

  STEP_ZO_G complete — $(date)
  Precomp:  ${OUTDIR}/
  Sim mats: ${OUTDIR}/sim_mats/

  Sanity check:
    Rscript -e 'x <- readRDS("${OUTDIR}/C_gar_LG28.precomp.rds");
                pc_cols <- length(grep("^PC_1_", names(x\$dt)));
                cat("rows:", nrow(x\$dt), "  PC_1_* cols:", pc_cols,
                    "  sim_mat:", paste(dim(x\$sim_mat), collapse="x"), "\n")'

  Next: _shared/STEP_05_L1_detect.R, _shared/STEP_07_L2_detect.R.
================================================================
EOF
