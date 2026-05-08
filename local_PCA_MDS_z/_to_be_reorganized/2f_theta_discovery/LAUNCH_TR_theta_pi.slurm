#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH -t 04:00:00
#SBATCH -J theta_pi
#SBATCH -o logs/theta_pi.%A_%a.out
#SBATCH -e logs/theta_pi.%A_%a.err

# =============================================================================
# LAUNCH_TR_theta_pi.slurm  (v8.5)
#
# One focal chromosome per array task. Mirrors LAUNCH_B01_mds_stage1.slurm
# conventions:
#   - reads chrom.list path from argv[1]
#   - sources the central inversion config (RSCRIPT_BIN, inv_init_dirs)
#   - sources 2f_theta_discovery's stream config (PESTPG_DIR, JSON_OUT_DIR…)
#   - SLURM_ARRAY_TASK_ID maps to chrom.list line N+1 (line 1 = header)
#   - calls STEP_TR_A then STEP_TR_B with --chrom <CHROM>
#
# θπ has no cross-chromosomal step (each chrom is fully independent), so
# there is no "stage 2" merge. The atlas loads one JSON per chromosome.
#
# Submit:
#   sbatch --array=1-28 LAUNCH_TR_theta_pi.slurm chrom.list
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ---- Configs ---------------------------------------------------------------
# Central inversion config (provides RSCRIPT_BIN, inv_init_dirs, INVDIR…)
# Try new repo layout first, fall back to legacy v8.5 location.
BASE_DIR="${BASE_DIR:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
for cand in \
    "${BASE_DIR}/inversion-popgen-toolkit/00_inversion_config.sh" \
    "${BASE_DIR}/inversion-popgen-toolkit/inversion_modules/00_inversion_config.sh" \
    "${BASE_DIR}/inversion_codebase_v8.5/00_inversion_config.sh"; do
  if [[ -f "$cand" ]]; then CONFIG_CENTRAL="$cand"; break; fi
done
if [[ -z "${CONFIG_CENTRAL:-}" ]]; then
  echo "Missing central config (tried new + legacy paths)" >&2
  exit 1
fi
echo "[launcher] using central config: ${CONFIG_CENTRAL}"

set -a
source "${CONFIG_CENTRAL}"
set +a

# Stream-specific config (provides PESTPG_DIR, JSON_OUT_DIR, OUTROOT, CHROM_LIST…)
# Resolve relative to this script's directory so it works regardless of repo
# layout. Falls back to INVDIR-relative for legacy compat.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for cand in \
    "${SCRIPT_DIR}/00_theta_config.sh" \
    "${INVDIR:-${BASE_DIR}/inversion-popgen-toolkit/inversion_modules}/phase_2_discovery/2f_theta_discovery/00_theta_config.sh" \
    "${BASE_DIR}/inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery/00_theta_config.sh"; do
  if [[ -f "$cand" ]]; then CONFIG_THETA="$cand"; break; fi
done
if [[ -z "${CONFIG_THETA:-}" ]]; then
  echo "Missing theta config (tried new + legacy paths)" >&2
  exit 1
fi
echo "[launcher] using theta config:   ${CONFIG_THETA}"

set -a
source "${CONFIG_THETA}"
set +a

# B01-style: ensure pipeline output directories exist
if declare -F inv_init_dirs >/dev/null; then
  inv_init_dirs
fi
mkdir -p "${OUTROOT}" "${JSON_OUT_DIR}" "${LOG_DIR:-${OUTROOT}/logs}" logs

# ---- Argv ------------------------------------------------------------------
CHROM_LIST="${1:?ERROR: provide chrom.list as argv1 (one chrom per line, line 1 = header)}"
[[ -f "${CHROM_LIST}" ]] || { echo "[ERROR] chrom.list not found: ${CHROM_LIST}" >&2; exit 1; }

# ---- Pick this task's chromosome -------------------------------------------
# B01 idiom: SLURM_ARRAY_TASK_ID=N reads line N+1 (skip header). For task 1,
# this reads line 2 of chrom.list. The companion chrom.list MUST therefore
# have a header line at line 1.
IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))

CHROM="$(sed -n "${LINE}p" "${CHROM_LIST}")"
if [[ -z "${CHROM}" ]]; then
  echo "[ERROR] Could not read chromosome from ${CHROM_LIST} line ${LINE}" >&2
  exit 1
fi

# ---- Where do the R scripts live? ------------------------------------------
# Mirror B01: scripts in <stream-dir>/steps/ if that exists, else in <stream-dir>/.
TR_STEPS="${THETA_CONFIG_DIR}/steps"
[[ -d "${TR_STEPS}" ]] || TR_STEPS="${THETA_CONFIG_DIR}"

# ---- Banner ----------------------------------------------------------------
echo "================================================================"
echo "  θπ phase-2 — ${CHROM}"
echo "  Array task: ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX:-?}"
echo "  PESTPG_DIR: ${PESTPG_DIR}"
echo "  PESTPG_SCALE: ${PESTPG_SCALE}"
echo "  JSON out:   ${JSON_OUT_DIR}/${CHROM}/"
echo "  JobID:      ${SLURM_JOB_ID}"
echo "  Started:    $(date)"
echo "================================================================"

# ---- Pipeline (no D17) ----------------------------------------------------
echo ""
echo "[$(date)] [${CHROM}] === STEP_TR_A precompute ==="
"${RSCRIPT_BIN:-Rscript}" "${TR_STEPS}/STEP_TR_A_compute_theta_matrices.R" \
  --chrom "${CHROM}"

echo ""
echo "[$(date)] [${CHROM}] === STEP_TR_B classifier ==="
"${RSCRIPT_BIN:-Rscript}" "${TR_STEPS}/STEP_TR_B_classify_theta.R" \
  --chrom "${CHROM}"

echo ""
echo "[$(date)] [${CHROM}] === DONE (no D17) ==="
ls -lh "${JSON_OUT_DIR}/${CHROM}/${CHROM}_phase2_theta.json"
