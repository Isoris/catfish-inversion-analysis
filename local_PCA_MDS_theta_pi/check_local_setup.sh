#!/usr/bin/env bash
set -uo pipefail
source /mnt/c/Users/quent/Desktop/catfish-inversion-analysis/local_PCA_MDS_theta_pi/00_theta_config.sh
echo
echo "=== printed ==="
theta_config_print
echo
echo "=== sanity ==="
if [ -d "$PESTPG_DIR" ]; then
  echo "OK   PESTPG_DIR    : $PESTPG_DIR"
else
  echo "FAIL PESTPG_DIR    : $PESTPG_DIR"
fi
n=$(find "$PESTPG_DIR" -maxdepth 1 -name "*.${PESTPG_SCALE}.pestPG" -type f 2>/dev/null | wc -l)
echo "     pestPG files at scale ${PESTPG_SCALE}: $n   [expect 226]"

if [ -f "$SAMPLE_LIST" ]; then
  nsamp=$(wc -l < "$SAMPLE_LIST")
  echo "OK   SAMPLE_LIST   : $SAMPLE_LIST   [$nsamp lines]"
else
  echo "FAIL SAMPLE_LIST   : $SAMPLE_LIST"
fi

[ -f "$REF" ]            && echo "OK   REF           : $REF" || echo "FAIL REF           : $REF"
[ -d "$THETA_TSV_DIR" ]  && echo "OK   THETA_TSV_DIR : $THETA_TSV_DIR" || echo "MAKE THETA_TSV_DIR : $THETA_TSV_DIR"
[ -d "$JSON_OUT_DIR" ]   && echo "OK   JSON_OUT_DIR  : $JSON_OUT_DIR"  || echo "MAKE JSON_OUT_DIR  : $JSON_OUT_DIR"
[ -x "$RSCRIPT" ]        && echo "OK   RSCRIPT       : $RSCRIPT"       || echo "WARN RSCRIPT       : $RSCRIPT [not executable]"

echo
echo "=== first pestPG header peek ==="
first=$(find "$PESTPG_DIR" -maxdepth 1 -name "*.${PESTPG_SCALE}.pestPG" -type f 2>/dev/null | head -1)
if [ -n "$first" ]; then
  echo "file: $first"
  head -2 "$first"
else
  echo "no pestPG file found at scale"
fi