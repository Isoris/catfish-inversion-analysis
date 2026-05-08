#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH -J inv_C04cd_ghsl
#SBATCH --array=1-28
#SBATCH -o logs/inv_C04cd_ghsl_%A_%a.out
#SBATCH -e logs/inv_C04cd_ghsl_%A_%a.err

# =============================================================================
# LAUNCH_STEP_C04cd_ghsl_enrichment.slurm
#
# Phase 2 / 2e_ghsl — page-3 atlas enrichment runner. Per chromosome, runs:
#   STEP_C04c_ghsl_local_pca   → <chr>.ghsl_v6_localpca.rds
#   STEP_C04d_ghsl_d17_wrapper → <chr>_ghsl_d17{L1,L2}_{envelopes,boundaries}.tsv
#   export_ghsl_to_json_v3     → <chr>_phase2_ghsl.json (page-3 atlas)
#
# Upstream prerequisites:
#   STEP_C04 v6 already ran  (matrices RDS + annot RDS + per_sample RDS + karyotypes RDS)
#   STEP_C04b classifier ran (annot RDS carries PASS/WEAK/FAIL)
#   D17 multipass scripts present at the configured paths (defaulted from inversion_config)
#
# Replaces the legacy run_ghsl_array.sh — same pipeline but config-driven and
# wired for the new inversion_modules/ layout.
#
# Usage:
#   sbatch LAUNCH_STEP_C04cd_ghsl_enrichment.slurm                 # array 1–28
#   CHR_OVERRIDE=C_gar_LG28 sbatch --array=1 LAUNCH_STEP_C04cd_ghsl_enrichment.slurm
#
# Tuning env overrides:
#   GHSL_PAD                 default 1   (local-PCA window padding)
#   GHSL_SMOOTHING_SCALE     default none (override with s50 etc. if raw too noisy)
#   GHSL_ANCHOR              default auto
#   GHSL_MIN_SAMPLES         default 50
#   GHSL_Z_THRESHOLD         default 2.5
#   GHSL_MIN_L2_WINDOWS      default 5
#   GHSL_MERGE_GAP           default 3
#   GHSL_PRIMARY_SCALE       default s50  (which rolling scale the export uses for sub-blocks)
#   GHSL_MAX_K               default 6
#   D17_L1_SCRIPT            default ${BASE}/utils/d17/STEP_D17_multipass_L1_only_v7.R
#   D17_L2_SCRIPT            default ${BASE}/utils/d17/STEP_D17_multipass_L2_v8.R
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# ── Central config (provides BASE, RSCRIPT_BIN, GHSL_DIR, JSON_OUT_DIR, ...) ──
CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit/inversion_modules/00_inversion_config.sh}"
# Fall back to legacy locations if the new path doesn't exist yet
if [[ ! -f "${CONFIG}" ]]; then
  for cand in \
      "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_modules/00_inversion_config.sh" \
      "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/00_inversion_config.sh"; do
    if [[ -f "$cand" ]]; then CONFIG="$cand"; break; fi
  done
fi
[[ -f "${CONFIG}" ]] || { echo "Missing config (tried new + legacy paths)" >&2; exit 1; }
echo "[C04cd] using config: ${CONFIG}"

set -a
source "${CONFIG}"
set +a

BRIDGE_SH="${BRIDGE_SH:-${BASE}/utils/pipeline_bridge.sh}"
[[ -f "${BRIDGE_SH}" ]] && source "${BRIDGE_SH}"

if declare -F inv_init_dirs >/dev/null; then
  inv_init_dirs
fi
mkdir -p logs

# ── Resolve chromosome from array index (or override) ─────────────────────
if [[ -n "${CHR_OVERRIDE:-}" ]]; then
  CHR="${CHR_OVERRIDE}"
else
  CHR="$(printf "C_gar_LG%02d" "${SLURM_ARRAY_TASK_ID}")"
fi

# ── Resolve dirs (mirror the canonical C04 / C04b launchers) ──────────────
GHSL_V6_DIR="${GHSL_V6_DIR:-${GHSL_DIR}/matrices}"
GHSL_LOCALPCA_DIR="${GHSL_LOCALPCA_DIR:-${GHSL_DIR}/localpca}"
GHSL_D17_DIR="${GHSL_D17_DIR:-${GHSL_DIR}/d17}"
JSON_OUT="${JSON_OUT:-${GHSL_DIR}/json_out}"
mkdir -p "${GHSL_LOCALPCA_DIR}" "${GHSL_D17_DIR}" "${JSON_OUT}"

# Required upstream RDSes
ANNOT_RDS="${GHSL_V6_DIR}/${CHR}.ghsl_v6.annot.rds"
PERSAMP_RDS="${GHSL_V6_DIR}/${CHR}.ghsl_v6.per_sample.rds"
KARYO_RDS="${GHSL_V6_DIR}/${CHR}.ghsl_v6.karyotypes.rds"
MATRICES_RDS="${GHSL_V6_DIR}/${CHR}.ghsl_v6_matrices.rds"
for f in "${MATRICES_RDS}" "${ANNOT_RDS}" "${PERSAMP_RDS}" "${KARYO_RDS}"; do
  [[ -f "$f" ]] || { echo "[C04cd] ERROR: missing upstream RDS — $f" >&2; exit 2; }
done

# D17 scripts (overridable via env)
D17_L1_SCRIPT="${D17_L1_SCRIPT:-${BASE}/utils/d17/STEP_D17_multipass_L1_only_v7.R}"
D17_L2_SCRIPT="${D17_L2_SCRIPT:-${BASE}/utils/d17/STEP_D17_multipass_L2_v8.R}"
[[ -f "${D17_L1_SCRIPT}" ]] || {
  echo "[C04cd] ERROR: D17 L1 script missing — ${D17_L1_SCRIPT}" >&2
  echo "[C04cd]        Set D17_L1_SCRIPT env var to point at it." >&2
  exit 2;
}
[[ -f "${D17_L2_SCRIPT}" ]] || {
  echo "[C04cd] ERROR: D17 L2 script missing — ${D17_L2_SCRIPT}" >&2
  echo "[C04cd]        Set D17_L2_SCRIPT env var to point at it." >&2
  exit 2;
}

# Tuning knobs (env-overridable)
PAD="${GHSL_PAD:-1}"
SMOOTHING_SCALE="${GHSL_SMOOTHING_SCALE:-none}"
ANCHOR="${GHSL_ANCHOR:-auto}"
MIN_SAMPLES="${GHSL_MIN_SAMPLES:-50}"
Z_THRESHOLD="${GHSL_Z_THRESHOLD:-2.5}"
MIN_L2_WINDOWS="${GHSL_MIN_L2_WINDOWS:-5}"
MERGE_GAP="${GHSL_MERGE_GAP:-3}"
PRIMARY_SCALE="${GHSL_PRIMARY_SCALE:-s50}"
MAX_K="${GHSL_MAX_K:-6}"

SELF_DIR="${BASE}/inversion_modules/phase_2_discovery/2e_ghsl"

echo "================================================================"
echo "  C04cd: GHSL page-3 enrichment — ${CHR}"
echo "  Array task:    ${SLURM_ARRAY_TASK_ID:-?} of ${SLURM_ARRAY_TASK_MAX:-?}"
echo "  Matrices RDS:  ${MATRICES_RDS}"
echo "  Annot RDS:     ${ANNOT_RDS}"
echo "  Per-sample RDS:${PERSAMP_RDS}"
echo "  Karyo RDS:     ${KARYO_RDS}"
echo "  Local-PCA out: ${GHSL_LOCALPCA_DIR}"
echo "  D17 out:       ${GHSL_D17_DIR}"
echo "  JSON out:      ${JSON_OUT}"
echo "  D17 L1 script: ${D17_L1_SCRIPT}"
echo "  D17 L2 script: ${D17_L2_SCRIPT}"
echo "  Knobs: pad=${PAD} smoothing=${SMOOTHING_SCALE} anchor=${ANCHOR}"
echo "         min_samples=${MIN_SAMPLES} z_thr=${Z_THRESHOLD}"
echo "         min_l2=${MIN_L2_WINDOWS} merge_gap=${MERGE_GAP}"
echo "         primary_scale=${PRIMARY_SCALE} max_k=${MAX_K}"
echo "  Started:       $(date)"
echo "================================================================"

# ── Stage 1: STEP_C04c local PCA ──────────────────────────────────────────
echo "[$(date)] [${CHR}] === STEP_C04c local PCA ==="
"${RSCRIPT_BIN}" "${SELF_DIR}/STEP_C04c_ghsl_local_pca.R" \
  --chrom         "${CHR}" \
  --ghsl_matrices "${MATRICES_RDS}" \
  --out_rds       "${GHSL_LOCALPCA_DIR}/${CHR}.ghsl_v6_localpca.rds" \
  --pad           "${PAD}" \
  --smoothing-scale "${SMOOTHING_SCALE}" \
  --anchor        "${ANCHOR}" \
  --min-samples-per-window "${MIN_SAMPLES}" \
  --z-threshold   "${Z_THRESHOLD}" \
  --min-l2-windows "${MIN_L2_WINDOWS}" \
  --merge-gap     "${MERGE_GAP}"

# ── Stage 2: STEP_C04d D17 wrapper ────────────────────────────────────────
echo "[$(date)] [${CHR}] === STEP_C04d D17 wrapper ==="
"${RSCRIPT_BIN}" "${SELF_DIR}/STEP_C04d_ghsl_d17_wrapper.R" \
  --chrom          "${CHR}" \
  --localpca       "${GHSL_LOCALPCA_DIR}/${CHR}.ghsl_v6_localpca.rds" \
  --out_dir        "${GHSL_D17_DIR}" \
  --d17_l1_script  "${D17_L1_SCRIPT}" \
  --d17_l2_script  "${D17_L2_SCRIPT}"

# ── Stage 3: export_ghsl_to_json_v3 ───────────────────────────────────────
echo "[$(date)] [${CHR}] === export_ghsl_to_json_v3 ==="
"${RSCRIPT_BIN}" "${SELF_DIR}/export_ghsl_to_json_v3.R" \
  --chrom          "${CHR}" \
  --annot_rds      "${ANNOT_RDS}" \
  --persamp_rds    "${PERSAMP_RDS}" \
  --karyo_rds      "${KARYO_RDS}" \
  --localpca_rds   "${GHSL_LOCALPCA_DIR}/${CHR}.ghsl_v6_localpca.rds" \
  --d17_l1_env     "${GHSL_D17_DIR}/${CHR}_ghsl_d17L1_envelopes.tsv" \
  --d17_l2_env     "${GHSL_D17_DIR}/${CHR}_ghsl_d17L2_envelopes.tsv" \
  --d17_l1_bnd     "${GHSL_D17_DIR}/${CHR}_ghsl_d17L1_boundaries.tsv" \
  --d17_l2_bnd     "${GHSL_D17_DIR}/${CHR}_ghsl_d17L2_boundaries.tsv" \
  --out_dir        "${JSON_OUT}" \
  --primary_scale  "${PRIMARY_SCALE}" \
  --max_k          "${MAX_K}"

echo "================================================================"
echo "  C04cd complete — ${CHR} — $(date)"
echo "  Local-PCA RDS: ${GHSL_LOCALPCA_DIR}/${CHR}.ghsl_v6_localpca.rds"
echo "  Phase-2 JSON:  ${JSON_OUT}/${CHR}/${CHR}_phase2_ghsl.json"
echo "================================================================"

ls -lh "${JSON_OUT}/${CHR}/${CHR}_phase2_ghsl.json"
