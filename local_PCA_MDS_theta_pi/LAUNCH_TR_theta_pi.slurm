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
# LAUNCH_TR_theta_pi.slurm  (v10, 2026-05-16, harmonized layout — dual-scale, A→J)
#
# Harmonized v10 paths (from _shared/00_inversion_config.sh):
#   PATH2_LOCAL_PCA = 01_local_pca/         (TR_A / TR_B output)
#   PATH2_REGISTRY  = 02_dense_registry/    (TR_D registry artifacts)
#   PATH2_MDS       = 03_mds/               (reserved; TR_C embeds MDS in precomp)
#   PATH2_PRECOMP   = 04_precomp/           (TR_C --mode full)
#   PATH2_PRECOMP_DENSE = 04_precomp_dense/ (TR_C --mode local)
#   PATH2_L1 / L1_PLOTS = 05_L1 / 06_L1_plots
#   PATH2_L2 / L2_PLOTS = 07_L2 / 08_L2_plots
#   PATH2_CARRIERS  = 08_carriers/          (TR_H output)
#   PATH2_CUSUM     = 09_cusum/             (TR_I output)
#   PATH2_JSON      = 09_atlas_json/        (TR_J output)
#
# L1/L2 detect+plot calls now route to _shared/STEP_0{5,6,7,8}_*.R. Plot
# steps honor MAKE_L1_PLOTS / MAKE_L2_PLOTS (default true).
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
# Shared L1/L2 detect+plot scripts now live in _shared/ (one level up).
SCRIPT_DIR_SHARED="${SCRIPT_DIR_SHARED:-$(cd "${SCRIPT_DIR}/../_shared" 2>/dev/null && pwd || echo "${SCRIPT_DIR}/../_shared")}"

# Plot toggles — wired through to _shared/STEP_06 and _shared/STEP_08.
export MAKE_L1_PLOTS="${MAKE_L1_PLOTS:-true}"
export MAKE_L2_PLOTS="${MAKE_L2_PLOTS:-true}"

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
echo "[$(date)] [${CHROM}] === C coarse — MDS + features (mode=full → 04_precomp/) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_C_mds_compute.R" \
  --chrom      "${CHROM}" \
  --mode       full \
  --out_subdir "${SLOT_PRECOMP:-04_precomp}"

echo ""
echo "[$(date)] [${CHROM}] === D — _shared/STEP_05_L1_detect (nn80 → 05_L1/) ==="
"${RSCRIPT}" "${SCRIPT_DIR_SHARED}/STEP_05_L1_detect.R" \
  --precomp_dir "${PATH2_PRECOMP}" \
  --chr         "${CHROM}" \
  --outdir      "${PATH2_L1}"

echo ""
echo "[$(date)] [${CHROM}] === E — _shared/STEP_06_L1_plot (MAKE_L1_PLOTS=${MAKE_L1_PLOTS}) ==="
"${RSCRIPT}" "${SCRIPT_DIR_SHARED}/STEP_06_L1_plot.R" \
  --precomp_dir "${PATH2_PRECOMP}" \
  --L1_dir      "${PATH2_L1}" \
  --chr         "${CHROM}" \
  --outdir      "${PATH2_L1_PLOTS}"

echo ""
echo "[$(date)] [${CHROM}] === F — _shared/STEP_07_L2_detect (nn40 → 07_L2/) ==="
"${RSCRIPT}" "${SCRIPT_DIR_SHARED}/STEP_07_L2_detect.R" \
  --precomp_dir "${PATH2_PRECOMP}" \
  --L1_dir      "${PATH2_L1}" \
  --chr         "${CHROM}" \
  --outdir      "${PATH2_L2}"

echo ""
echo "[$(date)] [${CHROM}] === G — _shared/STEP_08_L2_plot (MAKE_L2_PLOTS=${MAKE_L2_PLOTS}) ==="
"${RSCRIPT}" "${SCRIPT_DIR_SHARED}/STEP_08_L2_plot.R" \
  --precomp_dir "${PATH2_PRECOMP}" \
  --L1_dir      "${PATH2_L1}" \
  --L2_dir      "${PATH2_L2}" \
  --chr         "${CHROM}" \
  --outdir      "${PATH2_L2_PLOTS}"

echo ""
echo "[$(date)] [${CHROM}] === TR_D — merge_registry (genome-wide candidates → 02_dense_registry/) ==="
# Optional in dual-scale per-chrom layout — safe to call per-chrom because
# TR_D scans precomp_dir each time. The launcher could also defer this to
# a single afterok dependency, but inline keeps the chain simple.
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_D_merge_registry.R" \
  --precomp_dir "${PATH2_PRECOMP}" \
  --outdir      "${PATH2_REGISTRY}"

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
echo "[$(date)] [${CHROM}] === C dense — banded sim_mat only (mode=local → 04_precomp_dense/) ==="
"${RSCRIPT}" "${TR_STEPS}/STEP_TR_C_mds_compute.R" \
  --chrom         "${CHROM}" \
  --mode          local \
  --sim-band-half 100 \
  --out_subdir    "${SLOT_PRECOMP_DENSE:-04_precomp_dense}"

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
  --coarse_precomp_dir  "${PATH2_PRECOMP}" \
  --dense_precomp_dir   "${PATH2_PRECOMP_DENSE}" \
  --coarse_pestpg_scale "${COARSE_PESTPG_SCALE}" \
  --dense_pestpg_scale  "${DENSE_PESTPG_SCALE}"

echo ""
echo "[$(date)] [${CHROM}] === DONE (dual-scale A→J) ==="
ls -lh "${JSON_OUT_DIR}/${CHROM}/${CHROM}_phase2_theta.json" 2>/dev/null \
  || echo "[WARN] expected JSON not found at ${JSON_OUT_DIR}/${CHROM}/${CHROM}_phase2_theta.json"
