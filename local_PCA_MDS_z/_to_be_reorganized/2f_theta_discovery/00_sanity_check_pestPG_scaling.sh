#!/usr/bin/env bash
# =============================================================================
# 00_sanity_check_pestPG_scaling.sh
# =============================================================================
# Run BEFORE patching STEP_TR_A. Confirms that the pestPG `tP` column is
# the per-window SUM (not per-site density) by checking magnitude against
# expected per-site θπ for a vertebrate.
#
# Expected scale on your data:
#   per-site θπ in a vertebrate hatchery cohort ≈ 1e-3 to 1e-2
#   window size = 10000 bp, nSites typically a few thousand
#   sum tP per window ≈ nSites × per-site θπ ≈ 2000–9000 × 1e-3 = 2–9
#
# What this script reports:
#   - first 5 windows from one sample's pestPG: WinCenter, tP, nSites
#   - mean(tP), mean(nSites), implied per-site = mean(tP)/mean(nSites)
#
# If mean(tP) is in the range 0.5 – 50: tP is a SUM (bug confirmed).
# If mean(tP) is in the range 1e-4 – 1e-2: tP is already a density.
#
# Usage: bash 00_sanity_check_pestPG_scaling.sh [path/to/one.pestPG]
# =============================================================================

set -euo pipefail

if [[ $# -ge 1 ]]; then
  PESTPG="$1"
else
  source 00_theta_config.sh 2>/dev/null || true
  PESTPG=$(ls "${PESTPG_DIR}"/*.${PESTPG_SCALE}.pestPG 2>/dev/null | head -1)
fi

if [[ -z "${PESTPG:-}" || ! -f "${PESTPG}" ]]; then
  echo "ERROR: no pestPG found. Pass one as argument."
  exit 1
fi

echo "Inspecting: ${PESTPG}"
echo
echo "Header line:"
head -1 "${PESTPG}"
echo
echo "First 5 data rows (just Chr / WinCenter / tP / nSites):"
grep -v '^#' "${PESTPG}" | head -5 | awk -F'\t' '{print $2, $3, $5, $14}'
echo
echo "Summary across LG28 (or whichever chr):"
grep -v '^#' "${PESTPG}" | awk -F'\t' -v chr="${CHROM:-LG28}" '
  $2 ~ chr {
    n++
    sum_tP += $5
    sum_ns += $14
    if ($5 > max_tP) max_tP = $5
    if ($14 > 0) {
      sum_per_site += $5 / $14
      n_per_site++
    }
  }
  END {
    if (n == 0) { print "no rows on chr "chr; exit 1 }
    printf "  n_windows           : %d\n", n
    printf "  mean(tP)            : %.6g\n", sum_tP / n
    printf "  max(tP)             : %.6g\n", max_tP
    printf "  mean(nSites)        : %.1f\n", sum_ns / n
    printf "  mean(tP / nSites)   : %.6e   <-- this is per-site θπ if tP is a sum\n", sum_per_site / n_per_site
    print  ""
    if (sum_tP / n > 0.1) {
      print "  VERDICT: tP looks like a SUM. Patch STEP_TR_A."
    } else {
      print "  VERDICT: tP looks like a per-site density already. Investigate before patching."
    }
  }'
