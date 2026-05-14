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
# STEP_ZO_G_LAUNCH_precompute.slurm
#
# Per-chromosome precompute: builds <chr>.precomp.rds + the suite of NN-smoothed
# similarity matrices (<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds).
# Parallelised across chromosomes via mclapply.
#
# Submit AFTER all ZO_E per-chrom MDS array tasks complete. ZO_G reads
# <OUTDIR>/tmp/<chr>.mds_perchr.rds directly — the old ZO_F merge step is
# no longer needed.
#
# Output (path_localpca_zblocks/04_precomp/precomp/):
#   <chr>.precomp.rds
#   sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
#   ../precomp_summary.tsv
#
# Sweet spot — 28 cores (1:1 with chromosomes), 237G memory cap. Real usage <80G.
# Tune the NN scales saved via:  NN_SIM_SCALES="40,80,160,320" sbatch ...
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
MDS_PREFIX="${MDS_PREFIX:-${SCRATCH_ROOT}/path_localpca_zblocks/03_mds/inversion_localpca}"
OUTDIR="${OUTDIR:-${SCRATCH_ROOT}/path_localpca_zblocks/04_precomp}"
DOSAGE_DIR="${DOSAGE_DIR:-${SCRATCH_ROOT}/02_dosage_sites}"

# ZO_G reads per-chrom MDS results from <OUTDIR>/mds_perchr/<chr>.mds_perchr.rds
# (produced by ZO_E). Older runs wrote to <OUTDIR>/tmp/ — that's the back-compat
# fallback (handled inside ZO_G). Sanity-check at least one of those exists.
ZO_E_DIR=""
for cand in "${OUTDIR}/mds_perchr" "${OUTDIR}/tmp"; do
  if [[ -d "$cand" ]] && find "$cand" -maxdepth 1 -name "*.mds_perchr.rds" -print -quit | grep -q .; then
    ZO_E_DIR="$cand"; break
  fi
done
if [[ -z "$ZO_E_DIR" ]]; then
  echo "[G] missing <chr>.mds_perchr.rds in ${OUTDIR}/{mds_perchr,tmp}/ — run ZO_E first" >&2
  exit 1
fi
echo "[G] using per-chrom MDS dir: ${ZO_E_DIR}"

mkdir -p "${OUTDIR}"

cat <<EOF
================================================================
  STEP_ZO_G precompute
  ZO_E perchr   : ${ZO_E_DIR}
  Out dir       : ${OUTDIR}
  NN_SIM_SCALES : ${NN_SIM_SCALES:-20,40,80,120,160,200,240,320}
  Cores         : ${SLURM_CPUS_PER_TASK:-1}
================================================================
EOF

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_ZO_G_precompute.R" \
  "${MDS_PREFIX}" \
  "${OUTDIR}"

cat <<EOF

  STEP_ZO_G complete — $(date)
  Precomp:  ${OUTDIR}/precomp/
  Sim mats: ${OUTDIR}/precomp/sim_mats/

  Sanity check:
    Rscript -e 'x <- readRDS("${OUTDIR}/precomp/C_gar_LG28.precomp.rds");
                pc_cols <- length(grep("^PC_1_", names(x\$dt)));
                cat("rows:", nrow(x\$dt), "  PC_1_* cols:", pc_cols,
                    "  sim_mat:", paste(dim(x\$sim_mat), collapse="x"), "\n")'

  Next: STEP_ZO_H detect_L1, STEP_ZO_J detect_L2.
================================================================
EOF
