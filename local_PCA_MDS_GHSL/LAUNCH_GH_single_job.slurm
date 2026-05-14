#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH -t 4:00:00
#SBATCH -J ghsl_all
#SBATCH -o logs/ghsl_all.%j.out
#SBATCH -e logs/ghsl_all.%j.err

# =============================================================================
# LAUNCH_GH_single_job.slurm  (2026-05-14)
#
# GHSL phase-3 pipeline, single-job parallel layout (mirrors theta_pi).
# For each chromosome that has Clair3 postprocess data, runs:
#
#   PREP : STEP_GH_prep_merged_phased_snps.sh
#          (postprocess/<chr>/<sample>/all_variants_with_phase.tsv  ->
#           <chr>.merged_phased_snps.tsv.gz)
#   A    : STEP_GH_A_compute_matrices.R   (divergence matrices, multi-scale)
#   B    : STEP_GH_B_local_pca_compute.R  (per-window local PCA, scale s50)
#   C    : STEP_GH_C_mds_compute.R        (MDS + sim_mat + NN smoothing)
#   D    : STEP_GH_D_merge_registry.R     (per-chrom registry — runs once
#                                          after all chroms finish)
#
# Auto-detects which chroms have postprocess data in
# ${GHSL_POSTPROCESS_DIR}/<chr>/<sample>/ — useful while Clair3 is still
# running and only a subset of chroms is ready.
#
# Submit:
#   sbatch LAUNCH_GH_single_job.slurm
#
# Override at submit time:
#   N_PARALLEL=8 GHSL_SCALES=10,50 sbatch LAUNCH_GH_single_job.slurm
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ---- Central config (single explicit path, no search loop) ----------------
CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_localpca_v7/_shared/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
echo "[launcher] central config: ${CONFIG}"
set -a; source "${CONFIG}"; set +a
# Use SCRIPT_DIR_GHSL from the central config (the canonical GHSL script
# location). We deliberately do NOT derive from $(dirname BASH_SOURCE)
# because under sbatch, the launcher is copied to /var/spool/slurmd/job<ID>/
# and BASH_SOURCE points at the spool copy, not the real repo dir.
export SCRIPT_DIR="${SCRIPT_DIR_GHSL}"

# Required env from central config (must be exported for child R scripts).
export OUTROOT="${PATH3_ROOT}"
export GHSL_POSTPROCESS_DIR GHSL_PREP_DIR GHSL_MATRICES_DIR \
       GHSL_LOCALPCA_DIR GHSL_PRECOMP_DIR GHSL_SCALES \
       GHSL_PCA_SCALE_FULL GHSL_PCA_SCALE_DENSE \
       PRECOMP_DIR SAMPLES_IND RSCRIPT_BIN SCRIPT_DIR

mkdir -p logs "${GHSL_PREP_DIR}" "${GHSL_MATRICES_DIR}" "${GHSL_LOCALPCA_DIR}" \
         "${GHSL_PRECOMP_DIR}" "${PATH3_ROOT}/logs_perchrom"

# ---- BLAS thread pinning ---------------------------------------------------
# 64 CPU / 8 workers = 8 BLAS threads per worker.
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8

N_PARALLEL="${N_PARALLEL:-8}"

# ---- Auto-detect available chromosomes ------------------------------------
# A chrom is "available" if ${GHSL_POSTPROCESS_DIR}/<chr>/ exists and has
# at least one CGA*/ subdirectory.
if [[ ! -d "${GHSL_POSTPROCESS_DIR}" ]]; then
  echo "[ERR] postprocess dir missing: ${GHSL_POSTPROCESS_DIR}" >&2
  exit 1
fi

CHROMS=""
for chr_dir in "${GHSL_POSTPROCESS_DIR}"/C_gar_LG*; do
  [[ -d "$chr_dir" ]] || continue
  chr=$(basename "$chr_dir")
  n_samples=$(ls -d "$chr_dir"/CGA* 2>/dev/null | wc -l)
  if [[ $n_samples -gt 0 ]]; then
    CHROMS="${CHROMS}${chr}
"
  fi
done
CHROMS="$(echo "${CHROMS}" | awk 'NF>0')"
N_TOTAL="$(echo "${CHROMS}" | wc -l)"

echo "================================================================"
echo "  GHSL single-job — ${N_TOTAL} chroms with Clair3 data, ${N_PARALLEL} parallel"
echo "  Started:    $(date)"
echo "  POSTPROCESS: ${GHSL_POSTPROCESS_DIR}"
echo "  GHSL_PREP:   ${GHSL_PREP_DIR}"
echo "  OUTROOT:     ${OUTROOT}"
echo "  SCRIPT_DIR:  ${SCRIPT_DIR}"
echo "  Chroms:"
echo "${CHROMS}" | awk '{print "    "$0}'
echo "================================================================"

# ---- Per-chrom pipeline (PREP -> A -> B -> C) ------------------------------
PERCHROM_SCRIPT="$(mktemp /tmp/run_ghsl_chr_XXXXXX.sh)"
chmod +x "${PERCHROM_SCRIPT}"
trap "rm -f ${PERCHROM_SCRIPT}" EXIT

cat > "${PERCHROM_SCRIPT}" <<'PERCHROM_EOF'
#!/usr/bin/env bash
set -eo pipefail
CHROM="$1"
LOG="${PATH3_ROOT}/logs_perchrom/ghsl.${CHROM}.log"
exec >"${LOG}" 2>&1

echo "[$(date)] [${CHROM}] START"

# Skip-if-done: if precomp.rds already exists, this chrom is finished.
FINAL_PRECOMP="${GHSL_PRECOMP_DIR}/${CHROM}.ghsl_precomp.rds"
if [[ -s "${FINAL_PRECOMP}" ]]; then
  echo "[$(date)] [${CHROM}] SKIP — final precomp already present: ${FINAL_PRECOMP}"
  exit 0
fi

# Prep — merge per-sample Clair3 postprocess TSVs into per-chrom TSV.gz
PREP_OUT="${GHSL_PREP_DIR}/${CHROM}.merged_phased_snps.tsv.gz"
if [[ ! -s "${PREP_OUT}" ]]; then
  echo "[$(date)] [${CHROM}] === PREP ==="
  bash "${SCRIPT_DIR}/STEP_GH_prep_merged_phased_snps.sh" \
    "${GHSL_POSTPROCESS_DIR}" "${GHSL_PREP_DIR}" "${CHROM}"
else
  echo "[$(date)] [${CHROM}] prep already done, skipping"
fi

# A — divergence matrices
echo "[$(date)] [${CHROM}] === A compute_matrices ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_GH_A_compute_matrices.R" \
  --precomp_dir   "${PRECOMP_DIR}" \
  --ghsl_prep_dir "${GHSL_PREP_DIR}" \
  --sample_list   "${SAMPLES_IND}" \
  --outdir        "${GHSL_MATRICES_DIR}" \
  --chrom         "${CHROM}" \
  --scales        "${GHSL_SCALES}"

# B — local PCA per window
echo "[$(date)] [${CHROM}] === B local_pca ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_GH_B_local_pca_compute.R" \
  --matrices_dir "${GHSL_MATRICES_DIR}" \
  --outdir       "${GHSL_LOCALPCA_DIR}" \
  --scale        "${GHSL_PCA_SCALE_FULL}" \
  --chrom        "${CHROM}"

# C — MDS + sim_mat + NN smoothing (mode full)
echo "[$(date)] [${CHROM}] === C mds (mode=full) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_GH_C_mds_compute.R" \
  --chrom      "${CHROM}" \
  --mode       full \
  --out_subdir "03_precomp"

echo "[$(date)] [${CHROM}] DONE"
PERCHROM_EOF

# ---- Fan out via xargs -----------------------------------------------------
echo "[launcher] dispatching ${N_TOTAL} chroms × ${N_PARALLEL} parallel..."
echo "${CHROMS}" \
  | xargs -n1 -P"${N_PARALLEL}" -I{} bash -c "
      if ${PERCHROM_SCRIPT} {} ; then
        echo '[done] {}'
      else
        echo '[FAIL] {}  (see ${PATH3_ROOT}/logs_perchrom/ghsl.{}.log)'
      fi
    "

# ---- D — registry merge (once, after all per-chrom finished) ---------------
echo ""
echo "[$(date)] === D merge_registry (genome-wide) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_GH_D_merge_registry.R" \
  --precomp_dir "${GHSL_PRECOMP_DIR}" \
  --outdir      "${GHSL_PRECOMP_DIR}" \
  || echo "[WARN] D registry merge failed — partial state OK if some chroms FAILed"

echo ""
echo "================================================================"
echo "  GHSL pipeline complete. Finished: $(date)"
echo "================================================================"
echo "Outputs:"
echo "  prep:      ${GHSL_PREP_DIR}/<chr>.merged_phased_snps.tsv.gz"
echo "  matrices:  ${GHSL_MATRICES_DIR}/<chr>.ghsl_matrices.rds"
echo "  local_pca: ${GHSL_LOCALPCA_DIR}/<chr>.ghsl_window_pca.rds"
echo "  precomp:   ${GHSL_PRECOMP_DIR}/<chr>.ghsl_precomp.rds"
echo "  sim_mats:  ${GHSL_PRECOMP_DIR}/sim_mats/<chr>.sim_mat_nn{0,40,80,160,320}.rds"
echo ""
echo "Per-chrom logs: ${PATH3_ROOT}/logs_perchrom/ghsl.<chrom>.log"
echo "Re-run will SKIP chroms whose <chr>.ghsl_precomp.rds already exists."
