#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH -t 6:00:00
#SBATCH -J theta_pi_all
#SBATCH -o logs/theta_pi_all.%j.out
#SBATCH -e logs/theta_pi_all.%j.err

# =============================================================================
# LAUNCH_TR_theta_pi_single_job.slurm
#
# Drop-in alternative to LAUNCH_TR_theta_pi.slurm for sites that can't (or
# don't want to) use --array. Runs ALL chromosomes from chr.list inside a
# SINGLE job, fanned out N_PARALLEL=16 at a time via xargs (-P 16).
#
# Why 16?
#   64 CPU / 4 effective threads per worker ≈ 16 concurrent chroms.
#   128 GB / ~6-8 GB peak per worker (dense pass) ≈ 16 concurrent fits.
# With 16 workers all 28 chroms finish in roughly 2 batches → ~1-2 h wall.
# Tune via:
#   N_PARALLEL=12 sbatch LAUNCH_TR_theta_pi_single_job.slurm chr.list
#
# Each worker logs to ${OUTROOT}/logs/theta_pi.<chrom>.log so the main
# stdout stays readable. Use `tail -f` on a per-chrom log to see progress.
#
# Submit:
#   sbatch LAUNCH_TR_theta_pi_single_job.slurm chr.list
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ---- Configs (mirror LAUNCH_TR_theta_pi.slurm) -----------------------------
BASE_DIR="${BASE_DIR:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for cand in \
    "${BASE_DIR}/inversion_localpca_v7/_shared/00_inversion_config.sh" \
    "${BASE_DIR}/inversion-popgen-toolkit/00_inversion_config.sh" \
    "${BASE_DIR}/inversion-popgen-toolkit/inversion_modules/00_inversion_config.sh" \
    "${BASE_DIR}/inversion_codebase_v8.5/00_inversion_config.sh"; do
  if [[ -f "$cand" ]]; then CONFIG_CENTRAL="$cand"; break; fi
done
[[ -n "${CONFIG_CENTRAL:-}" ]] || { echo "no central config" >&2; exit 1; }
echo "[launcher] central config: ${CONFIG_CENTRAL}"
set -a; source "${CONFIG_CENTRAL}"; set +a

for cand in \
    "${SCRIPT_DIR}/00_theta_config.sh" \
    "${INVDIR:-${BASE_DIR}/inversion-popgen-toolkit/inversion_modules}/phase_2_discovery/2f_theta_discovery/00_theta_config.sh" \
    "${BASE_DIR}/inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery/00_theta_config.sh"; do
  if [[ -f "$cand" ]]; then CONFIG_THETA="$cand"; break; fi
done
[[ -n "${CONFIG_THETA:-}" ]] || { echo "no theta config" >&2; exit 1; }
echo "[launcher] theta config:   ${CONFIG_THETA}"
set -a; source "${CONFIG_THETA}"; set +a

if declare -F inv_init_dirs >/dev/null; then inv_init_dirs; fi
mkdir -p "${OUTROOT}" "${JSON_OUT_DIR}" "${LOG_DIR:-${OUTROOT}/logs}" logs

# ---- Scales ----------------------------------------------------------------
export COARSE_PESTPG_SCALE="${COARSE_PESTPG_SCALE:-win50000.step10000}"
export DENSE_PESTPG_SCALE="${DENSE_PESTPG_SCALE:-win10000.step2000}"

# ---- R scripts location ----------------------------------------------------
TR_STEPS="${SCRIPT_DIR}/steps"
[[ -d "${TR_STEPS}" ]] || TR_STEPS="${SCRIPT_DIR}"
export TR_STEPS
export OUTROOT JSON_OUT_DIR LOG_DIR
export COARSE_PESTPG_SCALE DENSE_PESTPG_SCALE
RSCRIPT="${RSCRIPT_BIN:-Rscript}"
export RSCRIPT

# ---- Limit BLAS threads inside each worker so we don't oversubscribe ------
# (64 CPU total / 16 workers = 4 BLAS threads per worker)
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4

# ---- Argv ------------------------------------------------------------------
CHROM_LIST="${1:?ERROR: provide chr.list as argv1 (one chrom per line, line 1 = header)}"
[[ -f "${CHROM_LIST}" ]] || { echo "[ERROR] chrom.list not found: ${CHROM_LIST}" >&2; exit 1; }

CHROMS="$(tail -n +2 "${CHROM_LIST}" | awk 'NF>0')"
N_TOTAL="$(echo "${CHROMS}" | wc -l)"
N_PARALLEL="${N_PARALLEL:-16}"

echo "================================================================"
echo "  θπ phase-2 single-job — ${N_TOTAL} chromosomes, ${N_PARALLEL} parallel"
echo "  Started:    $(date)"
echo "  OUTROOT:    ${OUTROOT}"
echo "  JSON out:   ${JSON_OUT_DIR}"
echo "  per-chrom logs: ${LOG_DIR}/theta_pi.<chrom>.log"
echo "================================================================"

# ---- Per-chrom pipeline (A → J) --------------------------------------------
# Written as an inline script so xargs can dispatch it cleanly.
PERCHROM_SCRIPT="$(mktemp /tmp/run_chr_XXXXXX.sh)"
chmod +x "${PERCHROM_SCRIPT}"
trap "rm -f ${PERCHROM_SCRIPT}" EXIT

cat > "${PERCHROM_SCRIPT}" <<'PERCHROM_EOF'
#!/usr/bin/env bash
set -eo pipefail
CHROM="$1"
LOG="${LOG_DIR:-${OUTROOT}/logs}/theta_pi.${CHROM}.log"
exec >"${LOG}" 2>&1

echo "[$(date)] [${CHROM}] START"

# Skip-if-done: if the final atlas JSON exists, this chrom is finished.
FINAL_JSON="${JSON_OUT_DIR}/${CHROM}/${CHROM}_phase2_theta.json"
if [[ -s "${FINAL_JSON}" ]]; then
  echo "[$(date)] [${CHROM}] SKIP — final JSON already present: ${FINAL_JSON}"
  exit 0
fi

# COARSE pass
export PESTPG_SCALE="${COARSE_PESTPG_SCALE}"
echo "[$(date)] [${CHROM}] === A coarse @ ${PESTPG_SCALE} ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_A_compute_theta_matrices.R" --chrom "${CHROM}"

echo "[$(date)] [${CHROM}] === B coarse ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_B_local_pca_compute.R"      --chrom "${CHROM}"

echo "[$(date)] [${CHROM}] === C coarse (mode=full) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_C_mds_compute.R"            --chrom "${CHROM}" --mode full

echo "[$(date)] [${CHROM}] === D detect_L1 ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_D_detect_L1.R" \
  --precomp_dir "${OUTROOT}/precomp" --chr "${CHROM}" --outdir "${OUTROOT}/L1_detect"

echo "[$(date)] [${CHROM}] === E plot_L1 ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_E_plot_L1.R" \
  --precomp_dir "${OUTROOT}/precomp" --L1_dir "${OUTROOT}/L1_detect" \
  --chr "${CHROM}" --outdir "${OUTROOT}/L1_plots"

echo "[$(date)] [${CHROM}] === F detect_L2 ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_F_detect_L2.R" \
  --precomp_dir "${OUTROOT}/precomp" --L1_dir "${OUTROOT}/L1_detect" \
  --chr "${CHROM}" --outdir "${OUTROOT}/L2_detect"

echo "[$(date)] [${CHROM}] === G plot_L2 ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_G_plot_L2.R" \
  --precomp_dir "${OUTROOT}/precomp" --L1_dir "${OUTROOT}/L1_detect" \
  --L2_dir "${OUTROOT}/L2_detect" --chr "${CHROM}" --outdir "${OUTROOT}/L2_plots"

# DENSE pass
export PESTPG_SCALE="${DENSE_PESTPG_SCALE}"
echo "[$(date)] [${CHROM}] === A dense @ ${PESTPG_SCALE} ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_A_compute_theta_matrices.R" --chrom "${CHROM}"

echo "[$(date)] [${CHROM}] === B dense ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_B_local_pca_compute.R"      --chrom "${CHROM}"

echo "[$(date)] [${CHROM}] === C dense (mode=local, banded) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_C_mds_compute.R" \
  --chrom "${CHROM}" --mode local --sim-band-half 100 --out_subdir precomp_dense

# MERGE pass
echo "[$(date)] [${CHROM}] === H classify_carriers ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_H_classify_carriers.R" --chr "${CHROM}"

echo "[$(date)] [${CHROM}] === I cusum_per_carrier ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_I_cusum_per_carrier.R" \
  --chr "${CHROM}" --lib "${TR_STEPS}/lib_persample_cusum.R"

echo "[$(date)] [${CHROM}] === J atlas_json ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_J_atlas_json.R" \
  --chr                 "${CHROM}" \
  --coarse_precomp_dir  "${OUTROOT}/precomp" \
  --dense_precomp_dir   "${OUTROOT}/precomp_dense" \
  --coarse_pestpg_scale "${COARSE_PESTPG_SCALE}" \
  --dense_pestpg_scale  "${DENSE_PESTPG_SCALE}"

echo "[$(date)] [${CHROM}] DONE"
PERCHROM_EOF

# ---- Fan out via xargs -----------------------------------------------------
echo "[launcher] dispatching ${N_TOTAL} chroms × ${N_PARALLEL} parallel..."
echo "${CHROMS}" \
  | xargs -n1 -P"${N_PARALLEL}" -I{} bash -c "
      if ${PERCHROM_SCRIPT} {} ; then
        echo '[done] {}'
      else
        echo '[FAIL] {}'
      fi
    "

echo ""
echo "================================================================"
echo "  All chromosomes processed.  Finished: $(date)"
echo "================================================================"
echo ""
echo "Final JSONs:"
ls -lh "${JSON_OUT_DIR}"/*/*_phase2_theta.json 2>/dev/null | awk '{print "  "$NF}' || echo "  (none found)"
echo ""
echo "Per-chrom logs: ${LOG_DIR}/theta_pi.<chrom>.log"
echo "Re-run will SKIP chroms with an existing _phase2_theta.json."
