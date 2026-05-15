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
#   D    : STEP_TR_D_detect_L1.R          (L1 boundaries on GHSL sim_mat_nn80)
#   E    : STEP_TR_E_plot_L1.R            (L1 overlay PDF)
#   F    : STEP_TR_F_detect_L2.R          (L2 boundaries inside each L1)
#   G    : STEP_TR_G_plot_L2.R            (L2 overlay PDF)
#   J    : STEP_GH_J_atlas_json.R         (slim atlas JSON: sim_mat + L1/L2;
#                                          no carriers / no CUSUM — GHSL has
#                                          no dense per-sample track yet)
#   K    : STEP_GH_D_merge_registry.R     (genome-wide registry — runs once
#                                          after all per-chrom finish)
#
# TR_D/E/F/G are pipeline-agnostic and operate on any precomp + sim_mats with
# the standard window_dt schema, so we reuse them against GHSL outputs.
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
# GH_C reads its local-PCA input dir from $OUT_GH_LOCAL_PCA_DIR (default
# "01_ghsl_local_pca"). Our config writes GH_B output to GHSL_LOCALPCA_DIR
# (02_local_pca by default), so bridge the two names here.
export OUT_GH_LOCAL_PCA_DIR="${GHSL_LOCALPCA_DIR}"
# Launcher-local defaults for the L1/L2/atlas dirs. These are also defined in
# the central config, but we set them here too so the launcher is robust to
# running against an older config that doesn't know about them.
export GHSL_L1_DIR="${GHSL_L1_DIR:-${PATH3_ROOT}/04_L1_detect}"
export GHSL_L1_PLOTS_DIR="${GHSL_L1_PLOTS_DIR:-${PATH3_ROOT}/05_L1_plots}"
export GHSL_L2_DIR="${GHSL_L2_DIR:-${PATH3_ROOT}/06_L2_detect}"
export GHSL_L2_PLOTS_DIR="${GHSL_L2_PLOTS_DIR:-${PATH3_ROOT}/07_L2_plots}"
export GHSL_JSON_DIR="${GHSL_JSON_DIR:-${PATH3_ROOT}/08_atlas_json}"
# Same defensive default for SCRIPT_DIR_THETA — needed for TR_D/E/F/G calls.
export SCRIPT_DIR_THETA="${SCRIPT_DIR_THETA:-${REPO_ROOT:-$(dirname "${SCRIPT_DIR}")}/local_PCA_MDS_theta_pi}"
export GHSL_POSTPROCESS_DIR GHSL_PREP_DIR GHSL_MATRICES_DIR \
       GHSL_LOCALPCA_DIR GHSL_PRECOMP_DIR \
       GHSL_L1_DIR GHSL_L1_PLOTS_DIR GHSL_L2_DIR GHSL_L2_PLOTS_DIR GHSL_JSON_DIR \
       GHSL_SCALES \
       GHSL_PCA_SCALE_FULL GHSL_PCA_SCALE_DENSE \
       PRECOMP_DIR SAMPLES_IND RSCRIPT_BIN SCRIPT_DIR SCRIPT_DIR_THETA

mkdir -p logs "${GHSL_PREP_DIR}" "${GHSL_MATRICES_DIR}" "${GHSL_LOCALPCA_DIR}" \
         "${GHSL_PRECOMP_DIR}" \
         "${GHSL_L1_DIR}" "${GHSL_L1_PLOTS_DIR}" \
         "${GHSL_L2_DIR}" "${GHSL_L2_PLOTS_DIR}" \
         "${GHSL_JSON_DIR}" \
         "${PATH3_ROOT}/logs_perchrom"

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

# ---- Per-chrom pipeline (PREP -> A -> B -> C -> D -> E -> F -> G -> J) -----
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

# Skip-if-done: if the final atlas JSON exists, this chrom has run all the
# way through GH_C + L1/L2 detection + GH_J. Earlier outputs are reused
# automatically by the individual step skip-checks.
FINAL_JSON="${GHSL_JSON_DIR}/${CHROM}/${CHROM}_phase3_ghsl.json"
if [[ -s "${FINAL_JSON}" ]]; then
  echo "[$(date)] [${CHROM}] SKIP — final atlas JSON already present: ${FINAL_JSON}"
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
if [[ ! -s "${GHSL_PRECOMP_DIR}/${CHROM}.ghsl_precomp.rds" ]]; then
  echo "[$(date)] [${CHROM}] === C mds (mode=full) ==="
  "${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_GH_C_mds_compute.R" \
    --chrom      "${CHROM}" \
    --mode       full \
    --out_subdir "03_precomp"
else
  echo "[$(date)] [${CHROM}] precomp already done, skipping C"
fi

# TR scripts (D-G) operate directly on a precomp + sim_mats and are
# pipeline-agnostic — we point them at the GHSL precomp.
GHSL_PRECOMP_F="${GHSL_PRECOMP_DIR}/${CHROM}.ghsl_precomp.rds"
SIM_NN80="${GHSL_PRECOMP_DIR}/sim_mats/${CHROM}.sim_mat_nn80.rds"
SIM_NN40="${GHSL_PRECOMP_DIR}/sim_mats/${CHROM}.sim_mat_nn40.rds"

# D — detect L1 boundaries
echo "[$(date)] [${CHROM}] === D detect_L1 ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR_THETA}/STEP_TR_D_detect_L1.R" \
  --precomp "${GHSL_PRECOMP_F}" --sim_mat "${SIM_NN80}" \
  --chr "${CHROM}" --outdir "${GHSL_L1_DIR}" \
  --boundary_scan TRUE \
  --boundary_validator_mode grow \
  --boundary_W 5 --boundary_offset 5 --boundary_min_dist 30

# E — plot L1
echo "[$(date)] [${CHROM}] === E plot_L1 ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR_THETA}/STEP_TR_E_plot_L1.R" \
  --precomp "${GHSL_PRECOMP_F}" --sim_mat "${SIM_NN80}" \
  --catalogue  "${GHSL_L1_DIR}/${CHROM}.L1_envelopes.tsv" \
  --boundaries "${GHSL_L1_DIR}/${CHROM}.L1_boundaries.tsv" \
  --chr "${CHROM}" --outdir "${GHSL_L1_PLOTS_DIR}" \
  --toggle_L1 yes --boundary_filter stable

# F — detect L2 boundaries (inside each L1 envelope)
echo "[$(date)] [${CHROM}] === F detect_L2 ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR_THETA}/STEP_TR_F_detect_L2.R" \
  --precomp "${GHSL_PRECOMP_F}" --sim_mat "${SIM_NN40}" \
  --l1_catalogue "${GHSL_L1_DIR}/${CHROM}.L1_envelopes.tsv" \
  --chr "${CHROM}" --outdir "${GHSL_L2_DIR}" \
  --boundary_scan TRUE \
  --boundary_validator_mode grow \
  --quadrant_validator yes \
  --weak_demote_score 0 \
  --quad_rescue_max_grow_z 1.5 \
  --quad_demote_on_fail yes \
  --quad_demote_drift_floor -1.0

# G — plot L2
echo "[$(date)] [${CHROM}] === G plot_L2 ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR_THETA}/STEP_TR_G_plot_L2.R" \
  --precomp "${GHSL_PRECOMP_F}" --sim_mat "${SIM_NN80}" --sim_mat_l2 "${SIM_NN40}" \
  --l1_catalogue  "${GHSL_L1_DIR}/${CHROM}.L1_envelopes.tsv" \
  --l1_boundaries "${GHSL_L1_DIR}/${CHROM}.L1_boundaries.tsv" \
  --l2_envelopes  "${GHSL_L2_DIR}/${CHROM}.L2_envelopes.tsv" \
  --l2_boundaries "${GHSL_L2_DIR}/${CHROM}.L2_boundaries.tsv" \
  --chr "${CHROM}" --outdir "${GHSL_L2_PLOTS_DIR}" \
  --boundary_filter stable

# J — atlas JSON (slim, no carrier / no CUSUM — GHSL has no dense per-sample track yet)
echo "[$(date)] [${CHROM}] === J atlas_json ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_GH_J_atlas_json.R" \
  --chr         "${CHROM}" \
  --precomp_dir "${GHSL_PRECOMP_DIR}" \
  --l1_dir      "${GHSL_L1_DIR}" \
  --l2_dir      "${GHSL_L2_DIR}" \
  --outdir      "${GHSL_JSON_DIR}" \
  --scale       "${GHSL_PCA_SCALE_FULL}"

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

# ---- K — registry merge (once, after all per-chrom finished) ---------------
# (Script name is STEP_GH_D_merge_registry.R for historical reasons; renamed
# in the log to K to disambiguate from per-chrom D=detect_L1 above.)
echo ""
echo "[$(date)] === K merge_registry (genome-wide) ==="
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_GH_D_merge_registry.R" \
  --precomp_dir "${GHSL_PRECOMP_DIR}" \
  --outdir      "${GHSL_PRECOMP_DIR}" \
  || echo "[WARN] K registry merge failed — partial state OK if some chroms FAILed"

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
echo "  L1:        ${GHSL_L1_DIR}/<chr>.L1_{envelopes,boundaries,score_curve}.tsv"
echo "  L1 plot:   ${GHSL_L1_PLOTS_DIR}/<chr>.L1_overlay.pdf"
echo "  L2:        ${GHSL_L2_DIR}/<chr>.L2_{envelopes,boundaries,segment_stats,quadrant_audit}.tsv"
echo "  L2 plot:   ${GHSL_L2_PLOTS_DIR}/<chr>.L2_overlay.pdf"
echo "  atlas:     ${GHSL_JSON_DIR}/<chr>/<chr>_phase3_ghsl.json"
echo ""
ls -lh "${GHSL_JSON_DIR}"/*/*_phase3_ghsl.json 2>/dev/null | awk '{print "  "$NF}' || echo "  (no atlas JSONs yet)"
echo ""
echo "Per-chrom logs: ${PATH3_ROOT}/logs_perchrom/ghsl.<chrom>.log"
echo "Re-run will SKIP chroms whose <chr>_phase3_ghsl.json already exists."
