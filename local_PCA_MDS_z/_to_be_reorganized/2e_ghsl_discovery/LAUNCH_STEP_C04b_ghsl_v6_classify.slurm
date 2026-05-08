#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 0:30:00
#SBATCH -J inv_C04b_cls
#SBATCH -o logs/inv_C04b_cls_%j.out
#SBATCH -e logs/inv_C04b_cls_%j.err

# =============================================================================
# LAUNCH_STEP_C04b_ghsl_v6_classify.slurm
#
# Phase 2 / 2e — GHSL v6 LIGHT CLASSIFIER. Runs in ~30 seconds per chromosome.
# Iterate 100 times while tuning scoring cutoffs; no heavy recompute.
#
# Upstream:    STEP_C04 (matrices RDS)
# Downstream:  STEP_D05_ghsl_stability (phase 2d candidate detection)
#
# Usage:
#   sbatch LAUNCH_STEP_C04b_ghsl_v6_classify.slurm                  # all chroms
#   sbatch LAUNCH_STEP_C04b_ghsl_v6_classify.slurm C_gar_LG12       # single chrom
#
# Tuning env overrides:
#   GHSL_SCALE=50           rolling scale to use (must match a saved scale)
#   GHSL_KARYO_LO=0.15      bottom quantile → INV/INV
#   GHSL_KARYO_HI=0.70      top quantile → INV/nonINV
#   GHSL_KARYO_MIN_RUN=10   min consecutive windows for stable call
#   GHSL_MAX_K=5            max clusters for interval k-means
#   GHSL_INTERVALS=path     triangle_intervals.tsv.gz (enables per-interval class)
# =============================================================================

set -euo pipefail
source ~/.bashrc
mamba activate assembly

CONFIG="${CONFIG:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_modules/00_inversion_config.sh}"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

BRIDGE_SH="${BRIDGE_SH:-${BASE}/utils/pipeline_bridge.sh}"
[[ -f "${BRIDGE_SH}" ]] && source "${BRIDGE_SH}"

inv_init_dirs
mkdir -p logs

# ── Resolve dirs ─────────────────────────────────────────────────────────
MATRICES_DIR="${MATRICES_DIR:-${GHSL_DIR}/matrices}"
OUTDIR_RUN="${OUTDIR_RUN:-${GHSL_DIR}/classify}"
mkdir -p "${OUTDIR_RUN}"

[[ -d "${MATRICES_DIR}" ]] || {
  echo "[C04b] ERROR: matrices dir missing — ${MATRICES_DIR}" >&2
  echo "[C04b]        Run LAUNCH_STEP_C04_ghsl_v6_compute.slurm first." >&2
  exit 2;
}

# ── Optional single-chrom arg ────────────────────────────────────────────
CHR_ARG="${1:-}"
CHR_OPTS=""
if [[ -n "${CHR_ARG}" ]]; then
  CHR_OPTS="--chrom ${CHR_ARG}"
fi

# ── Tuning params (env with defaults) ────────────────────────────────────
SCALE="${GHSL_SCALE:-50}"
KARYO_LO="${GHSL_KARYO_LO:-0.15}"
KARYO_HI="${GHSL_KARYO_HI:-0.70}"
KARYO_MIN_RUN="${GHSL_KARYO_MIN_RUN:-10}"
RANK_WINDOW="${GHSL_RANK_WINDOW:-5}"
MAX_K="${GHSL_MAX_K:-5}"

INTERVAL_OPTS=""
if [[ -n "${GHSL_INTERVALS:-}" ]]; then
  [[ -f "${GHSL_INTERVALS}" ]] || {
    echo "[C04b] ERROR: GHSL_INTERVALS not found — ${GHSL_INTERVALS}" >&2; exit 2;
  }
  INTERVAL_OPTS="--intervals ${GHSL_INTERVALS}"
fi

echo "================================================================"
echo "  C04b: GHSL v6 LIGHT CLASSIFIER"
echo "  Chrom:         ${CHR_ARG:-all}"
echo "  Matrices dir:  ${MATRICES_DIR}"
echo "  Out dir:       ${OUTDIR_RUN}"
echo "  Scale:         ${SCALE}"
echo "  Karyo lo/hi:   ${KARYO_LO} / ${KARYO_HI}"
echo "  Min run:       ${KARYO_MIN_RUN}"
echo "  Rank window:   ${RANK_WINDOW}"
echo "  Max K:         ${MAX_K}"
echo "  Intervals:     ${GHSL_INTERVALS:-(none — Stages A+B only)}"
echo "  JobID:         ${SLURM_JOB_ID}"
echo "  Started:       $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_2_discovery/2e_ghsl/STEP_C04b_snake3_ghsl_classify.R" \
  "${MATRICES_DIR}" "${OUTDIR_RUN}" \
  ${CHR_OPTS} \
  --scale         "${SCALE}" \
  --karyo_lo      "${KARYO_LO}" \
  --karyo_hi      "${KARYO_HI}" \
  --karyo_min_run "${KARYO_MIN_RUN}" \
  --rank_window   "${RANK_WINDOW}" \
  --max_k         "${MAX_K}" \
  ${INTERVAL_OPTS}

echo "================================================================"
echo "  C04b complete — $(date)"
echo "  Outputs:"
echo "    ${OUTDIR_RUN}/snake3v6_window_track.tsv.gz"
echo "    ${OUTDIR_RUN}/snake3v6_karyotype_calls.tsv.gz"
echo "    ${OUTDIR_RUN}/snake3v6_summary.tsv"
[[ -n "${GHSL_INTERVALS:-}" ]] && echo "    ${OUTDIR_RUN}/snake3v6_interval_genotypes.tsv.gz"
echo "================================================================"
