#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH -t 2:00:00
#SBATCH -J inv_HtoM_all
#SBATCH -o logs/inv_HtoM_all.%j.out
#SBATCH -e logs/inv_HtoM_all.%j.err

# =============================================================================
# LAUNCH_ZO_HtoM_single_job.slurm
#
# Runs H/I/J/K/M (L1/L2 detection + plots + atlas JSON) for ALL chromosomes
# inside a SINGLE SLURM job, fanned out N_PARALLEL=16 at a time via xargs.
# Same pattern as LAUNCH_TR_theta_pi_single_job.slurm — use when SLURM
# arrays aren't available or when you got kicked off the login node mid-run.
#
# Each worker calls STEP_ZO_HtoM_run_one_chrom.sh for one chromosome.
# That wrapper internally checks for SAMPLE_META and skips ZO_M with a
# [WARN] if missing, completing H→K cleanly regardless.
#
# Why these sizes?
#   --cpus-per-task=64 / --mem=128G / -t 2:00:00
#   N_PARALLEL=16 (each worker pinned to 4 BLAS threads, 64/16=4)
#   Per-chrom cost ~ 2-3 min → 28 chroms / 16 parallel ≈ 5-10 min wall
#
# Tune:
#   N_PARALLEL=12 sbatch LAUNCH_ZO_HtoM_single_job.slurm
#   PRECOMP_DIR=/some/other/04_precomp_vN/precomp sbatch ...
#
# Submit:
#   sbatch LAUNCH_ZO_HtoM_single_job.slurm
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ---- Central config (single explicit path, no search loop) ----------------
CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_localpca_v7/_shared/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
echo "[launcher] central config: ${CONFIG}"
set -a; source "${CONFIG}"; set +a

# ---- Path resolution -------------------------------------------------------
# SCRIPT_DIR comes from the central config (set to ${REPO_ROOT}/local_PCA_MDS_z).
# We deliberately DO NOT use $(dirname "${BASH_SOURCE[0]}") here because under
# sbatch, BASH_SOURCE points at /var/spool/slurmd/job<JOBID>/ (where SLURM
# copied the launcher), not at the real script location.
export SCRATCH_ROOT="${SCRATCH_ROOT:-${BASE_DIR}/inversion_localpca_v7}"
export PATH1_ROOT="${PATH1_ROOT:-${SCRATCH_ROOT}/local_PCA_MDS_z}"
export PRECOMP_DIR="${PRECOMP_DIR:-${PATH1_ROOT}/04_precomp_v2/precomp}"
export SAMPLE_META="${SAMPLE_META:-${SCRATCH_ROOT}/_shared/sample_metadata.tsv}"

# ---- BLAS thread pinning ---------------------------------------------------
# 64 CPU / 16 workers = 4 BLAS threads per worker (no oversubscription)
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4

N_PARALLEL="${N_PARALLEL:-16}"

# ---- Logs ------------------------------------------------------------------
mkdir -p logs "${PATH1_ROOT}/logs_htom"

# ---- Sanity checks ---------------------------------------------------------
if ! ls "${PRECOMP_DIR}"/*.precomp.rds >/dev/null 2>&1; then
  echo "[ERR] no .precomp.rds files in ${PRECOMP_DIR}" >&2
  exit 1
fi
if [[ ! -x "${SCRIPT_DIR}/STEP_ZO_HtoM_run_one_chrom.sh" ]]; then
  echo "[ERR] wrapper not found or not executable: ${SCRIPT_DIR}/STEP_ZO_HtoM_run_one_chrom.sh" >&2
  exit 1
fi

CHROMS="$(ls "${PRECOMP_DIR}"/*.precomp.rds | sed 's|.*/||; s|\.precomp\.rds||')"
N_TOTAL="$(echo "${CHROMS}" | wc -l)"

echo "================================================================"
echo "  Z-path H→M single-job — ${N_TOTAL} chromosomes, ${N_PARALLEL} parallel"
echo "  Started:    $(date)"
echo "  SCRIPT_DIR:  ${SCRIPT_DIR}"
echo "  PRECOMP_DIR: ${PRECOMP_DIR}"
echo "  SAMPLE_META: ${SAMPLE_META}"
echo "  per-chrom logs: ${PATH1_ROOT}/logs_htom/htom.<chrom>.log"
echo "================================================================"

if [[ ! -f "${SAMPLE_META}" ]]; then
  echo "[WARN] SAMPLE_META not found at ${SAMPLE_META}"
  echo "[WARN] ZO_M will be skipped per chrom. H/I/J/K still complete."
fi

# ---- Fan out via xargs -----------------------------------------------------
echo "[launcher] dispatching ${N_TOTAL} chroms × ${N_PARALLEL} parallel..."
echo "${CHROMS}" \
  | xargs -n1 -P"${N_PARALLEL}" -I{} bash -c "
      LOG='${PATH1_ROOT}/logs_htom/htom.{}.log'
      if bash '${SCRIPT_DIR}/STEP_ZO_HtoM_run_one_chrom.sh' {} >\"\${LOG}\" 2>&1 ; then
        echo '[done] {}'
      else
        echo '[FAIL] {}  (see '\$LOG')'
      fi
    "

# ---- Scorecard -------------------------------------------------------------
echo ""
echo "================================================================"
echo "  All H→M chromosomes processed. Finished: $(date)"
echo "================================================================"
echo "Output counts by stage:"
for d in 05_L1 06_L1_plots 07_L2 08_L2_plots 09_atlas_json; do
  n=$(ls "${PATH1_ROOT}/${d}"/* 2>/dev/null | wc -l)
  echo "  ${PATH1_ROOT}/${d}: ${n} files"
done

# Highlight any FAIL chroms
fails=$(grep -E '^\[FAIL\]' /dev/null 2>/dev/null || true)
# (the actual FAIL lines were printed inline above; this is just a hint)
echo ""
echo "If any [FAIL] lines printed: tail ${PATH1_ROOT}/logs_htom/htom.<chrom>.log to diagnose."
