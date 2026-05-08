#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -A lt200308
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -J inv_C04_ghsl
#SBATCH --array=1-28
#SBATCH -o logs/inv_C04_ghsl_%A_%a.out
#SBATCH -e logs/inv_C04_ghsl_%A_%a.err

# =============================================================================
# LAUNCH_STEP_C04_ghsl_v6_compute.slurm
#
# Phase 2 / 2e — GHSL v6 HEAVY ENGINE. Run ONCE per chromosome.
# Loads 77M+ phased Clair3 SNPs, computes within-sample haplotype divergence
# at raw window resolution, applies rolling means at configurable scales,
# saves matrices as RDS for the light classifier (STEP_C04b) to iterate on.
#
# Upstream:    prepare_ghsl_v3_merged.sh (ghsl_prep/<chr>.merged_phased_snps.tsv.gz)
#              STEP_C01a_precompute (precomp/<chr>.precomp.rds for window coords)
# Downstream:  STEP_C04b_snake3_ghsl_classify (light classifier, 30 s / chrom)
#
# Usage:
#   sbatch LAUNCH_STEP_C04_ghsl_v6_compute.slurm        # array 1–28, all chroms
#   CHR_OVERRIDE=C_gar_LG12 sbatch --array=1 LAUNCH_STEP_C04_ghsl_v6_compute.slurm
#     # single chrom — run as 1-task array, pass CHR_OVERRIDE for clarity
#
# Sizing: ~1h per chrom per the v6 README. 64G is comfortable for 226 samples
# × a single chromosome. 4 CPU — the rolling stage is not highly parallel.
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

# ── Resolve chromosome from array index (or override) ────────────────────
if [[ -n "${CHR_OVERRIDE:-}" ]]; then
  CHR="${CHR_OVERRIDE}"
else
  CHR="$(printf "C_gar_LG%02d" "${SLURM_ARRAY_TASK_ID}")"
fi

# ── Resolve dirs ─────────────────────────────────────────────────────────
PRECOMP_DIR_RUN="${PRECOMP_DIR_RUN:-${MDS_DIR}/snake_regions_multiscale/precomp}"
GHSL_PREP_DIR="${GHSL_PREP_DIR:-${INVDIR}/ghsl_prep}"
OUTDIR_RUN="${OUTDIR_RUN:-${GHSL_DIR}/matrices}"
mkdir -p "${OUTDIR_RUN}"

[[ -d "${PRECOMP_DIR_RUN}" ]] || {
  echo "[C04] ERROR: precomp dir missing — ${PRECOMP_DIR_RUN}" >&2; exit 2;
}
[[ -d "${GHSL_PREP_DIR}" ]] || {
  echo "[C04] ERROR: ghsl prep dir missing — ${GHSL_PREP_DIR}" >&2
  echo "[C04]        Run utils/prepare_ghsl_v3_merged.sh first." >&2
  exit 2;
}

# ── Optional tuning env overrides ────────────────────────────────────────
SCALES="${GHSL_SCALES:-10,20,30,40,50,100}"     # chat-14 expanded default ladder
MIN_PHASED="${GHSL_MIN_PHASED:-3}"
QUAL_MIN="${GHSL_QUAL_MIN:-20}"
GQ_MIN="${GHSL_GQ_MIN:-10}"

echo "================================================================"
echo "  C04: GHSL v6 HEAVY ENGINE — ${CHR}"
echo "  Array task: ${SLURM_ARRAY_TASK_ID:-?} of ${SLURM_ARRAY_TASK_MAX:-?}"
echo "  Precomp dir:   ${PRECOMP_DIR_RUN}"
echo "  GHSL prep dir: ${GHSL_PREP_DIR}"
echo "  Out dir:       ${OUTDIR_RUN}"
echo "  Scales:        ${SCALES}"
echo "  Min phased:    ${MIN_PHASED}"
echo "  QUAL / GQ min: ${QUAL_MIN} / ${GQ_MIN}"
echo "  JobID:         ${SLURM_JOB_ID}"
echo "  Started:       $(date)"
echo "================================================================"

"${RSCRIPT_BIN}" \
  "${BASE}/inversion_modules/phase_2_discovery/2e_ghsl/STEP_C04_snake3_ghsl_v6.R" \
  "${PRECOMP_DIR_RUN}" "${GHSL_PREP_DIR}" "${OUTDIR_RUN}" \
  --chrom       "${CHR}" \
  --scales      "${SCALES}" \
  --min_phased  "${MIN_PHASED}" \
  --qual_min    "${QUAL_MIN}" \
  --gq_min      "${GQ_MIN}"

echo "================================================================"
echo "  C04 complete — ${CHR} — $(date)"
echo "  Matrices RDS: ${OUTDIR_RUN}/${CHR}.ghsl_v6_matrices.rds"
echo "  Next:         sbatch LAUNCH_STEP_C04b_ghsl_v6_classify.slurm"
echo "================================================================"
