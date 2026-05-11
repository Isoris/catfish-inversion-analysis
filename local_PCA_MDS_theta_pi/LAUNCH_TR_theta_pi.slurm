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
# LAUNCH_TR_theta_pi.slurm  (dual-scale, A→J)
#
# One whole chromosome per array task. Runs BOTH scales of the theta-pi
# pipeline in sequence per chromosome:
#
#   COARSE pass  (win50000.step10000, ~5k windows on LG28)
#     A  theta_native at coarse scale
#     B  local PCA per window
#     C  MDS + features + sim_mats   (--mode full → precomp/)
#     D–G  L1/L2 detect + plots (boundary discovery)
#
#   DENSE pass   (win10000.step2000, ~16.5k windows on LG28)
#     A  theta_native at dense scale
#     B  local PCA per window
#     C  MDS + features banded only  (--mode local --out_subdir precomp_dense)
#
#   MERGE pass   (consumes coarse L1/L2 + dense per-window data)
#     H  classify_carriers   (k-means on dense per-sample PC1 inside coarse L2)
#     I  cusum_per_carrier   (dense per-sample CUSUM inside coarse L2)
#     J  atlas_json          (merged page-12 JSON: coarse sim_mat + dense theta tracks)
#
# All steps operate at whole-chromosome granularity. The "candidate" concept
# at the merge layer = the auto-detected L2 envelope from the coarse scan.
#
# θπ has no cross-chromosomal step (each chrom is fully independent), so
# there is no "stage 2" merge. The atlas loads one JSON per chromosome.
#
# Mirrors LAUNCH_B01_mds_stage1.slurm conventions: reads chrom.list path
# from argv[1], sources central inversion config + theta stream config,
# SLURM_ARRAY_TASK_ID maps to chrom.list line N+1 (line 1 = header).
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

# ---- Scale definitions -----------------------------------------------------
# Two scales used in this dual-scale launcher. Override via env if needed.
export COARSE_PESTPG_SCALE="${COARSE_PESTPG_SCALE:-win50000.step10000}"
export DENSE_PESTPG_SCALE="${DENSE_PESTPG_SCALE:-win10000.step2000}"

# ---- Where do the R scripts live? ------------------------------------------
# Mirror B01: scripts in <stream-dir>/steps/ if that exists, else in <stream-dir>/.
# Use SCRIPT_DIR (computed above from BASH_SOURCE) — robust to repo layout and
# matches the same logic used to locate CONFIG_THETA. Earlier revisions
# referenced an unset ${THETA_CONFIG_DIR} which crashed under `set -u`.
TR_STEPS="${SCRIPT_DIR}/steps"
[[ -d "${TR_STEPS}" ]] || TR_STEPS="${SCRIPT_DIR}"

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

# ---- Pipeline — dual-scale chain A→J --------------------------------------
RSCRIPT="${RSCRIPT_BIN:-Rscript}"

# ════════════════════════════════════════════════════════════════════════════
# COARSE pass (win50000.step10000) — full sim_mat + MDS + L1/L2 detection
# ════════════════════════════════════════════════════════════════════════════
export PESTPG_SCALE="${COARSE_PESTPG_SCALE}"

echo ""
echo "[$(date)] [${CHROM}] === A coarse — theta_native @ ${PESTPG_SCALE} ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_A_compute_theta_matrices.R" \
  --chrom "${CHROM}"

echo ""
echo "[$(date)] [${CHROM}] === B coarse — local PCA per window ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_B_local_pca_compute.R" \
  --chrom "${CHROM}"

echo ""
echo "[$(date)] [${CHROM}] === C coarse — MDS + features (mode=full) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_C_mds_compute.R" \
  --chrom "${CHROM}" \
  --mode  full

echo ""
echo "[$(date)] [${CHROM}] === D — detect_L1 (nn80) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_D_detect_L1.R" \
  --precomp_dir "${OUTROOT}/precomp" \
  --chr         "${CHROM}" \
  --outdir      "${OUTROOT}/L1_detect"

echo ""
echo "[$(date)] [${CHROM}] === E — plot_L1 ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_E_plot_L1.R" \
  --precomp_dir "${OUTROOT}/precomp" \
  --L1_dir      "${OUTROOT}/L1_detect" \
  --chr         "${CHROM}" \
  --outdir      "${OUTROOT}/L1_plots"

echo ""
echo "[$(date)] [${CHROM}] === F — detect_L2 (nn40 inside each L1) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_F_detect_L2.R" \
  --precomp_dir "${OUTROOT}/precomp" \
  --L1_dir      "${OUTROOT}/L1_detect" \
  --chr         "${CHROM}" \
  --outdir      "${OUTROOT}/L2_detect"

echo ""
echo "[$(date)] [${CHROM}] === G — plot_L2 ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_G_plot_L2.R" \
  --precomp_dir "${OUTROOT}/precomp" \
  --L1_dir      "${OUTROOT}/L1_detect" \
  --L2_dir      "${OUTROOT}/L2_detect" \
  --chr         "${CHROM}" \
  --outdir      "${OUTROOT}/L2_plots"

# ════════════════════════════════════════════════════════════════════════════
# DENSE pass (win10000.step2000) — banded sim_mat only, no MDS, no NN-smooths
# Switch PESTPG_SCALE so TR_A/B/C/I all read the dense theta_native files.
# ════════════════════════════════════════════════════════════════════════════
export PESTPG_SCALE="${DENSE_PESTPG_SCALE}"

echo ""
echo "[$(date)] [${CHROM}] === A dense — theta_native @ ${PESTPG_SCALE} ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_A_compute_theta_matrices.R" \
  --chrom "${CHROM}"

echo ""
echo "[$(date)] [${CHROM}] === B dense — local PCA per window ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_B_local_pca_compute.R" \
  --chrom "${CHROM}"

echo ""
echo "[$(date)] [${CHROM}] === C dense — banded sim_mat only (mode=local) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_C_mds_compute.R" \
  --chrom         "${CHROM}" \
  --mode          local \
  --sim-band-half 100 \
  --out_subdir    precomp_dense

# ════════════════════════════════════════════════════════════════════════════
# MERGE pass — H/I/J consume coarse L1/L2 envelopes + dense per-window data.
# PESTPG_SCALE stays at dense so TR_I reads the dense theta_native TSV.
# H and I auto-prefer ${OUTROOT}/precomp_dense if present.
# J reads BOTH precomps and emits the merged page-12 JSON with grid_map.
# ════════════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] [${CHROM}] === H — classify_carriers (k-means on dense PC1 per coarse L2) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_H_classify_carriers.R" \
  --chr "${CHROM}"

echo ""
echo "[$(date)] [${CHROM}] === I — per-sample + band-mean CUSUM (dense theta inside coarse L2) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_I_cusum_per_carrier.R" \
  --chr "${CHROM}" \
  --lib "${TR_STEPS}/lib_persample_cusum.R"

echo ""
echo "[$(date)] [${CHROM}] === J — atlas JSON (merged dual-scale page-12 bundle) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_J_atlas_json.R" \
  --chr                 "${CHROM}" \
  --coarse_precomp_dir  "${OUTROOT}/precomp" \
  --dense_precomp_dir   "${OUTROOT}/precomp_dense" \
  --coarse_pestpg_scale "${COARSE_PESTPG_SCALE}" \
  --dense_pestpg_scale  "${DENSE_PESTPG_SCALE}"

echo ""
echo "[$(date)] [${CHROM}] === DONE (dual-scale A→J) ==="
ls -lh "${JSON_OUT_DIR}/${CHROM}/${CHROM}_phase2_theta.json" 2>/dev/null \
  || echo "[WARN] expected JSON not found at ${JSON_OUT_DIR}/${CHROM}/${CHROM}_phase2_theta.json"
