#!/usr/bin/env bash
# =============================================================================
# STEP_GH_prep_merged_phased_snps.sh
#
# GHSL pre-stage — merges per-sample Clair3 postprocess TSVs into one
# per-chromosome file consumed by STEP_GH_A_compute_matrices.R.
#
# Input:  <postprocess_dir>/<chr>/<sample>/all_variants_with_phase.tsv
# Output: <outdir>/<chr>.merged_phased_snps.tsv.gz
#         (the canonical ghsl_prep/ directory the GH_A launcher reads from)
#
# Filters: IS_SIMPLE_BIALLELIC=TRUE, IS_SNP=TRUE
# Keeps:   sample_id, pos, ref, alt, gt_class, is_phased, phase_gt,
#          ps, phase_block_id, phase_tier, qual, gq, dp
#
# Usage:
#   bash STEP_GH_prep_merged_phased_snps.sh <postprocess_dir> <outdir> [chrom]
#   bash STEP_GH_prep_merged_phased_snps.sh <postprocess_dir> <outdir> C_gar_LG28
#
# SLURM array (one task per chromosome):
#   sbatch --array=1-28 STEP_GH_prep_merged_phased_snps.sh <postprocess_dir> <outdir>
# =============================================================================

set -euo pipefail

PP_DIR="${1:?Usage: STEP_GH_prep_merged_phased_snps.sh <postprocess_dir> <outdir> [chrom]}"
OUTDIR="${2:?Usage: STEP_GH_prep_merged_phased_snps.sh <postprocess_dir> <outdir> [chrom]}"
CHROM="${3:-}"

mkdir -p "${OUTDIR}"

# Resolve chromosome from SLURM array or argument
if [[ -z "${CHROM}" ]] && [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  # Auto-detect chromosomes from postprocess directory
  CHROMS=($(ls -d "${PP_DIR}"/C_gar_LG* 2>/dev/null | xargs -n1 basename | sort))
  IDX=$((SLURM_ARRAY_TASK_ID - 1))
  if [[ ${IDX} -ge ${#CHROMS[@]} ]]; then
    echo "[SKIP] Task ${SLURM_ARRAY_TASK_ID} exceeds chromosome count ${#CHROMS[@]}"
    exit 0
  fi
  CHROM="${CHROMS[${IDX}]}"
fi

if [[ -z "${CHROM}" ]]; then
  # No chrom specified: run all sequentially
  echo "[GH_prep] Running all chromosomes sequentially..."
  for CHR_DIR in "${PP_DIR}"/C_gar_LG*; do
    CHR=$(basename "${CHR_DIR}")
    [[ -d "${CHR_DIR}" ]] || continue
    echo "[GH_prep] Processing ${CHR}..."
    bash "$0" "${PP_DIR}" "${OUTDIR}" "${CHR}"
  done
  echo "[GH_prep] All done."
  exit 0
fi

# Single chromosome processing
CHR_DIR="${PP_DIR}/${CHROM}"
if [[ ! -d "${CHR_DIR}" ]]; then
  echo "[ERROR] No directory: ${CHR_DIR}" >&2
  exit 1
fi

OUTFILE="${OUTDIR}/${CHROM}.merged_phased_snps.tsv.gz"
TMPFILE="${OUTDIR}/.tmp_${CHROM}_merged.tsv"

echo "[GH_prep] Chromosome: ${CHROM}"
echo "[GH_prep] Input: ${CHR_DIR}"
echo "[GH_prep] Output: ${OUTFILE}"

# Write header
echo -e "sample_id\tpos\tref\talt\tgt_class\tis_phased\tphase_gt\tps\tphase_block_id\tphase_tier\tqual\tgq\tdp" > "${TMPFILE}"

N_SAMPLES=0
N_VARIANTS=0

for SAMPLE_DIR in "${CHR_DIR}"/CGA*; do
  [[ -d "${SAMPLE_DIR}" ]] || continue
  SAMPLE=$(basename "${SAMPLE_DIR}")
  VFILE="${SAMPLE_DIR}/all_variants_with_phase.tsv"

  if [[ ! -f "${VFILE}" ]]; then
    echo "  [SKIP] ${SAMPLE}: no all_variants_with_phase.tsv"
    continue
  fi

  # Extract biallelic SNPs with key columns
  # Columns: CHROM=1, POS=2, REF=4, ALT1=5, QUAL=8, IS_SNP=13,
  #          IS_SIMPLE_BIALLELIC=15, GT_CLASS=20, GQ=21, DP=22,
  #          IS_PHASED=41, PHASE_GT=42, PS=43, PHASE_BLOCK_ID=44, PHASE_TIER=45
  N_NEW=$(awk -F'\t' -v sample="${SAMPLE}" '
    NR == 1 { next }
    $13 == "1" && $15 == "1" {
      print sample "\t" $2 "\t" $4 "\t" $5 "\t" $20 "\t" $41 "\t" $42 "\t" $43 "\t" $44 "\t" $45 "\t" $8 "\t" $21 "\t" $22
      n++
    }
    END { print n > "/dev/stderr" }
  ' "${VFILE}" 2>&1 >> "${TMPFILE}" | tail -1)

  N_SAMPLES=$((N_SAMPLES + 1))
  N_VARIANTS=$((N_VARIANTS + ${N_NEW:-0}))

  if (( N_SAMPLES % 50 == 0 )); then
    echo "  [GH_prep] ${N_SAMPLES} samples processed, ${N_VARIANTS} variants so far"
  fi
done

echo "[GH_prep] ${N_SAMPLES} samples, ${N_VARIANTS} variants total"

# Sort by pos then sample, compress
echo "[GH_prep] Sorting and compressing..."
(head -1 "${TMPFILE}"; tail -n +2 "${TMPFILE}" | sort -t$'\t' -k2,2n -k1,1) | gzip -c > "${OUTFILE}"
rm -f "${TMPFILE}"

echo "[GH_prep] Done: ${OUTFILE} ($(du -h "${OUTFILE}" | cut -f1))"
