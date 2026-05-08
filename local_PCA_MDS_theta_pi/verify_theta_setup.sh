#!/usr/bin/env bash
# =============================================================================
# verify_theta_setup.sh
# =============================================================================
# Pre-flight check before running TR_A / TR_B on LANTA.
#
# Verifies:
#   1. We're in the right directory (Atlas/_scripts/cluster_R/2f_theta_discovery/)
#   2. 00_theta_config.sh is sourceable and exports the expected vars
#   3. PESTPG_DIR contains 226 .pestPG files at the configured scale
#   4. STEP_TR_A has the per-site scaling fix (ANGSD issue #329)
#   5. STEP_TR_B has the turn-114 schema-v2 patch
#   6. The sanity check script confirms tP is a SUM (the bug is real in the data)
#
# Usage:
#   bash verify_theta_setup.sh
#
# Exits 0 if everything is ready to launch. Exits 1 with a clear message
# at the first failed check.
# =============================================================================

set -euo pipefail

PASS=0
FAIL=0
WARN=0

ok()   { echo "  [ ok ] $*"; PASS=$((PASS+1)); }
fail() { echo "  [FAIL] $*"; FAIL=$((FAIL+1)); }
warn() { echo "  [WARN] $*"; WARN=$((WARN+1)); }
hr()   { echo ""; echo "── $* ──"; }

echo "=========================================================================="
echo "  Theta-pi pre-flight check"
echo "  $(date)"
echo "=========================================================================="

# ── 1. Working directory ───────────────────────────────────────────────────
hr "1. Directory check"
PWD_NOW="$(pwd)"
if [[ "$PWD_NOW" =~ inversion_modules/phase_2_discovery/2f_theta_discovery$ ]]; then
  ok "in inversion_modules/phase_2_discovery/2f_theta_discovery"
elif [[ "$PWD_NOW" =~ 2f_theta_discovery$ ]]; then
  warn "in a 2f_theta_discovery dir but not under inversion_modules/phase_2_discovery/"
  warn "  current: $PWD_NOW"
  warn "  expected suffix: inversion_modules/phase_2_discovery/2f_theta_discovery"
  warn "  This is OK if the codebase has a different layout. Continuing..."
else
  fail "not in a 2f_theta_discovery directory."
  fail "  current: $PWD_NOW"
  fail "  cd to inversion_modules/phase_2_discovery/2f_theta_discovery first."
  exit 1
fi

# ── 2. Required files present ─────────────────────────────────────────────
hr "2. Required files present"
for f in 00_theta_config.sh \
         00_sanity_check_pestPG_scaling.sh \
         STEP_TR_A_compute_theta_matrices.R \
         STEP_TR_B_classify_theta.R \
         LAUNCH_TR_theta_pi.slurm \
         chrom.list; do
  if [[ -f "$f" ]]; then
    ok "$f"
  else
    fail "$f missing"
  fi
done

if [[ $FAIL -gt 0 ]]; then
  echo ""
  echo "FAIL: required files missing — check upload."
  exit 1
fi

# ── 3. Source the config ──────────────────────────────────────────────────
hr "3. Source 00_theta_config.sh"
# shellcheck disable=SC1091
source ./00_theta_config.sh
for v in PESTPG_DIR PESTPG_SCALE SAMPLE_LIST THETA_TSV_DIR JSON_OUT_DIR OUTROOT; do
  val="${!v:-}"
  if [[ -n "$val" ]]; then
    ok "$v = $val"
  else
    fail "$v is empty after sourcing 00_theta_config.sh"
  fi
done

if [[ $FAIL -gt 0 ]]; then
  echo ""
  echo "FAIL: config did not export required vars."
  exit 1
fi

# ── 4. pestPG inputs ──────────────────────────────────────────────────────
hr "4. pestPG inputs at PESTPG_DIR"
if [[ ! -d "$PESTPG_DIR" ]]; then
  fail "PESTPG_DIR does not exist: $PESTPG_DIR"
  fail "  → Upstream MODULE_3 / 02_run_heterozygosity.sh hasn't run yet."
  exit 1
fi
ok "PESTPG_DIR exists: $PESTPG_DIR"

n_pestpg=$(find "$PESTPG_DIR" -maxdepth 1 -name "*.${PESTPG_SCALE}.pestPG" -type f 2>/dev/null | wc -l)
if [[ "$n_pestpg" -ge 220 && "$n_pestpg" -le 230 ]]; then
  ok "found $n_pestpg pestPG files at scale $PESTPG_SCALE (expect ~226)"
elif [[ "$n_pestpg" -gt 0 ]]; then
  warn "found only $n_pestpg pestPG files (expected ~226). Some samples may be missing."
else
  fail "no pestPG files found at $PESTPG_DIR/*.${PESTPG_SCALE}.pestPG"
  fail "  → Check PESTPG_SCALE setting (currently: $PESTPG_SCALE)"
  fail "  → Or upstream MODULE_3 hasn't run for this scale."
  exit 1
fi

# Sample list
if [[ -f "$SAMPLE_LIST" ]]; then
  n_samp=$(grep -c -v '^$' "$SAMPLE_LIST" 2>/dev/null || echo 0)
  ok "SAMPLE_LIST has $n_samp samples"
else
  fail "SAMPLE_LIST not found: $SAMPLE_LIST"
fi

# ── 5. TR_A has the per-site fix ──────────────────────────────────────────
hr "5. STEP_TR_A has the per-site scaling fix (ANGSD GitHub issue #329)"
if grep -q "tP / nSites\|theta_pi_per_site\|GitHub issue #329" STEP_TR_A_compute_theta_matrices.R; then
  hits=$(grep -c "tP / nSites\|theta_pi_per_site\|GitHub issue #329" STEP_TR_A_compute_theta_matrices.R)
  ok "found $hits markers of the per-site fix in STEP_TR_A"
else
  fail "STEP_TR_A is the legacy version — per-site fix NOT applied"
  fail "  → tP would be treated as per-site density, but it's actually a window SUM."
  fail "  → Coverage artefacts will masquerade as biology."
  fail "  → Drop in the patched STEP_TR_A from this bundle."
  exit 1
fi

# ── 6. TR_B has the turn-114 schema-v2 patch ──────────────────────────────
hr "6. STEP_TR_B has the turn-114 atlas-canonical-names patch"
if grep -q "schema_version = 2L\|sample_ids" STEP_TR_B_classify_theta.R; then
  hits=$(grep -c "schema_version = 2L\|sample_ids\|values_flat\|^.*values  *=.*clean_numeric" STEP_TR_B_classify_theta.R)
  ok "found $hits markers of the turn-114 patch in STEP_TR_B"
else
  warn "STEP_TR_B is the pre-turn-114 version — atlas-canonical names NOT emitted"
  warn "  → JSON will have legacy field shape (.samples-as-objects, .z_profile only)."
  warn "  → Atlas detector handles this via dual-shape tolerance, but cleaner is better."
  warn "  → Drop in the patched STEP_TR_B from this bundle for new runs."
fi

# ── 7. pestPG sanity check (confirms the bug is real in real data) ────────
hr "7. pestPG sanity (confirms tP is a SUM, not per-site)"
sample_pestpg="$(find "$PESTPG_DIR" -maxdepth 1 -name "*.${PESTPG_SCALE}.pestPG" -type f 2>/dev/null | head -1)"
if [[ -z "${sample_pestpg:-}" ]]; then
  warn "no pestPG to sanity-check (skipped)"
else
  echo "  inspecting: $sample_pestpg"
  verdict=$(bash 00_sanity_check_pestPG_scaling.sh "$sample_pestpg" 2>&1 | grep "VERDICT:" | head -1)
  if echo "$verdict" | grep -q "SUM"; then
    ok "sanity verdict: tP is a SUM (matches expected; per-site fix in TR_A is necessary)"
  elif echo "$verdict" | grep -q "per-site density"; then
    warn "sanity verdict says tP is already per-site density"
    warn "  → unusual; investigate before running TR_A. Check ANGSD version."
  else
    warn "sanity check gave no clear verdict:"
    echo "  $verdict"
  fi
fi

# ── Summary ───────────────────────────────────────────────────────────────
echo ""
echo "=========================================================================="
echo "  SUMMARY: $PASS passed, $WARN warnings, $FAIL failed"
echo "=========================================================================="

if [[ $FAIL -gt 0 ]]; then
  echo ""
  echo "STATUS: NOT READY — fix the failures above before running."
  exit 1
fi

if [[ $WARN -gt 0 ]]; then
  echo ""
  echo "STATUS: ready with warnings."
  echo ""
  echo "Next steps:"
  echo "  # dry-run on LG28 (interactive, ~15 min):"
  echo "  RSCRIPT=\$(which Rscript)"
  echo "  \$RSCRIPT STEP_TR_A_compute_theta_matrices.R --chrom C_gar_LG28"
  echo "  \$RSCRIPT STEP_TR_B_classify_theta.R         --chrom C_gar_LG28"
  echo "  jq '.schema_version, ._layers_present, .n_windows' \\"
  echo "     \${JSON_OUT_DIR}/C_gar_LG28/C_gar_LG28_phase2_theta.json"
  echo ""
  echo "  # if LG28 looks right, full 28-chrom array:"
  echo "  sbatch --array=1-28 LAUNCH_TR_theta_pi.slurm chrom.list"
  exit 0
fi

echo ""
echo "STATUS: READY ✓"
echo ""
echo "Next steps:"
echo "  # dry-run on LG28 (interactive, ~15 min):"
echo "  RSCRIPT=\$(which Rscript)"
echo "  \$RSCRIPT STEP_TR_A_compute_theta_matrices.R --chrom C_gar_LG28"
echo "  \$RSCRIPT STEP_TR_B_classify_theta.R         --chrom C_gar_LG28"
echo "  jq '.schema_version, ._layers_present, .n_windows' \\"
echo "     \${JSON_OUT_DIR}/C_gar_LG28/C_gar_LG28_phase2_theta.json"
echo ""
echo "  # if LG28 looks right, full 28-chrom array:"
echo "  sbatch --array=1-28 LAUNCH_TR_theta_pi.slurm chrom.list"
