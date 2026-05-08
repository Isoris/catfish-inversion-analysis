#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=237G
#SBATCH -t 2:00:00
#SBATCH -J inv_03_precomp
#SBATCH -o logs/inv_03_precomp.%j.out
#SBATCH -e logs/inv_03_precomp.%j.err

# =============================================================================
# 03_LAUNCH_precompute_localpca_zblocks.slurm  (v8.5 / consolidated layout v1.0)
#
# Per-chromosome precompute: builds <chr>.precomp.rds + the suite of NN-smoothed
# similarity matrices (<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds).
# Parallelised across chromosomes via mclapply.
#
# Submit AFTER 02b mds_merge completes:
#   sbatch --dependency=afterok:<JOBID_02b> 03_LAUNCH_precompute_localpca_zblocks.slurm
#
# Output (path_localpca_zblocks/04_precomp/precomp/):
#   <chr>.precomp.rds                                           (full, no slim)
#   sim_mats/<chr>.sim_mat_nn{0,20,40,80,120,160,200,240,320}.rds
#   ../window_dt.tsv.gz, ../precomp_summary.tsv
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

[[ -f "${MDS_PREFIX}.mds.rds" ]] || { echo "[03] missing ${MDS_PREFIX}.mds.rds — run 02b first" >&2; exit 1; }

mkdir -p "${OUTDIR}"

cat <<EOF
================================================================
  03 precompute_localpca_zblocks  (v10.0 SLIM)
  MDS prefix    : ${MDS_PREFIX}
  Out dir       : ${OUTDIR}
  Dosage dir    : ${DOSAGE_DIR:-<not set>}
  NN_SIM_SCALES : ${NN_SIM_SCALES:-20,40,80,120,160,200,240,320}
  Cores         : ${SLURM_CPUS_PER_TASK:-1}
================================================================
EOF

DOSAGE_ARG=()
if [[ -n "${DOSAGE_DIR:-}" && -d "${DOSAGE_DIR}" ]]; then
  DOSAGE_ARG=(--dosage_dir "${DOSAGE_DIR}")
else
  echo "[03] WARNING: DOSAGE_DIR not set — dosage het CV columns will be NA"
fi

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/03_precomp/03_precompute_localpca_zblocks.R" \
  "${MDS_PREFIX}" \
  "${OUTDIR}" \
  "${DOSAGE_ARG[@]}"

cat <<EOF

  03 complete — $(date)
  Precomp:  ${OUTDIR}/precomp/
  Sim mats: ${OUTDIR}/precomp/sim_mats/

  Sanity check (full precomp, no slim — must include PC_1_*/PC_2_*):
    Rscript -e 'x <- readRDS("${OUTDIR}/precomp/C_gar_LG28.precomp.rds");
                want <- c("max_abs_z","max_z_axis","inv_likeness",
                          "band_discreteness","adaptive_seed");
                missing <- setdiff(want, names(x\$dt));
                pc_cols <- length(grep("^PC_1_", names(x\$dt)));
                cat("missing dims:", length(missing), "  PC_1_* cols:", pc_cols, "\n")'

  Next: 04_detect_L1, 05_plot_L1, 06_detect_L2, 07_plot_L2 (interactive),
        08a_build_sample_metadata + 08b_export_atlas_json
================================================================
EOF
