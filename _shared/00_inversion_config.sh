#!/usr/bin/env bash
# =============================================================================
# 00_inversion_config.sh — v9 (2026-05-14, post slim-pipeline cleanup)
#
# Source with:
#   set -a; source "${CONFIG}"; set +a   # auto-export
#
# Every variable below uses the ${VAR:-default} pattern so anything set in
# the environment BEFORE sourcing this file wins. This prevents the previous
# trap where SCRIPT_DIR got clobbered when a launcher had already captured
# its own location.
#
# Layout assumed:
#   ${BASE}/
#   ├── catfish-inversion-analysis/                  (outer repo — may be stale)
#   │   └── catfish-inversion-analysis/              (nested clone w/ latest fixes)
#   │       ├── local_PCA_MDS_z/
#   │       ├── local_PCA_MDS_theta_pi/
#   │       ├── local_PCA_MDS_GHSL/
#   │       └── _shared/
#   ├── inversion_localpca_v7/                       (data tree)
#   │   ├── 01_beagle/  02_dosage_sites/  03_theta_pi_pestPG/  04_clair3_phased_GHSL/
#   │   ├── _shared/00_inversion_config.sh           (this file)
#   │   ├── local_PCA_MDS_z/  local_PCA_MDS_theta_pi/  local_PCA_MDS_GHSL/
#   │   └── chr.list
#   └── het_roh/02_heterozygosity/03_theta/multiscale/   (pestPG files for theta_pi)
# =============================================================================

# ── Project root ─────────────────────────────────────────────────────────────
export BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── Results root ─────────────────────────────────────────────────────────────
export SCRATCH_ROOT="${SCRATCH_ROOT:-${BASE}/inversion_localpca_v7}"
export INVDIR="${INVDIR:-${SCRATCH_ROOT}}"     # legacy alias for back-compat

# ── Code root (nested repo with latest fixes) ────────────────────────────────
# NB: the OUTER repo (${BASE}/catfish-inversion-analysis/local_PCA_MDS_z) is
# the stale May-11 snapshot. The active scripts live in the NESTED clone.
# Override via env if you reorganize the repo layout in the future.
export REPO_ROOT="${REPO_ROOT:-${BASE}/catfish-inversion-analysis/catfish-inversion-analysis}"
export SCRIPT_DIR="${SCRIPT_DIR:-${REPO_ROOT}/local_PCA_MDS_z}"
export SCRIPT_DIR_THETA="${SCRIPT_DIR_THETA:-${REPO_ROOT}/local_PCA_MDS_theta_pi}"
export SCRIPT_DIR_GHSL="${SCRIPT_DIR_GHSL:-${REPO_ROOT}/local_PCA_MDS_GHSL}"

# ── Rscript binary (mambaforge assembly env) ─────────────────────────────────
export RSCRIPT_BIN="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"

# ── Reference / sample lists (stable, unchanged from v8.5) ───────────────────
export HETDIR="${HETDIR:-${BASE}/het_roh}"
export REF="${REF:-${BASE}/00-samples/fClaHyb_Gar_LG.fa}"
export REF_FAI="${REF_FAI:-${REF}.fai}"
export BAMLIST="${BAMLIST:-${HETDIR}/01_inputs_check/bamlist_qcpass.txt}"
export SAMPLES_IND="${SAMPLES_IND:-${HETDIR}/01_inputs_check/samples.ind}"
# SAMPLE_LIST is the alias the theta_pi R scripts (TR_A/B/...) read from
# Sys.getenv(). One sample id per line, in the same order as PC_*_Ind*.
export SAMPLE_LIST="${SAMPLE_LIST:-${SAMPLES_IND}}"

# ── Input data (renamed source dirs) ─────────────────────────────────────────
export BEAGLE_DIR="${BEAGLE_DIR:-${SCRATCH_ROOT}/01_beagle}"
export DOSAGE_DIR="${DOSAGE_DIR:-${SCRATCH_ROOT}/02_dosage_sites}"
export THETAPI_DIR="${THETAPI_DIR:-${SCRATCH_ROOT}/03_theta_pi_pestPG}"
export GHSL_DIR="${GHSL_DIR:-${SCRATCH_ROOT}/04_clair3_phased_GHSL}"

# pestPG files actually live in het_roh on Lanta (the THETAPI_DIR above is
# empty; this is the real source path used by STEP_TR_A).
export PESTPG_DIR="${PESTPG_DIR:-${BASE}/het_roh/02_heterozygosity/03_theta/multiscale}"

# ── Path 1 (z-blocks / dosage) output tree ───────────────────────────────────
export PATH1_ROOT="${PATH1_ROOT:-${SCRATCH_ROOT}/local_PCA_MDS_z}"
export PATH1_LOCAL_PCA="${PATH1_LOCAL_PCA:-${PATH1_ROOT}/01_local_pca}"
export PATH1_REGISTRY="${PATH1_REGISTRY:-${PATH1_ROOT}/02_dense_registry}"
export PATH1_MDS="${PATH1_MDS:-${PATH1_ROOT}/03_mds}"

# Default precomp dir is now 04_precomp_v2 (the post-slim-ZO_G output).
# 04_precomp (without _v2) is the legacy chunked_2x/inv_likeness output and
# is preserved for archeology.
export PATH1_PRECOMP="${PATH1_PRECOMP:-${PATH1_ROOT}/04_precomp_v2}"
export PRECOMP_DIR="${PRECOMP_DIR:-${PATH1_PRECOMP}/precomp}"
export PATH1_L1="${PATH1_L1:-${PATH1_ROOT}/05_L1}"
export PATH1_L1_PLOTS="${PATH1_L1_PLOTS:-${PATH1_ROOT}/06_L1_plots}"
export PATH1_L2="${PATH1_L2:-${PATH1_ROOT}/07_L2}"
export PATH1_L2_PLOTS="${PATH1_L2_PLOTS:-${PATH1_ROOT}/08_L2_plots}"
export PATH1_JSON="${PATH1_JSON:-${PATH1_ROOT}/09_atlas_json}"

# Path 1 MDS prefix (ZO_E writes <prefix>.mds_metadata.tsv; ZO_G uses dirname
# to locate the per-chrom mds_perchr/ or tmp/ directory.)
export MDS_PREFIX_BASENAME="${MDS_PREFIX_BASENAME:-inversion_localpca}"
export PATH1_MDS_PREFIX="${PATH1_MDS_PREFIX:-${PATH1_MDS}/${MDS_PREFIX_BASENAME}}"
export MDS_PREFIX="${MDS_PREFIX:-${PATH1_MDS_PREFIX}}"

# ── Path 2 (theta_pi) and Path 3 (GHSL) output trees ─────────────────────────
export PATH2_ROOT="${PATH2_ROOT:-${SCRATCH_ROOT}/local_PCA_MDS_theta_pi}"
export PATH3_ROOT="${PATH3_ROOT:-${SCRATCH_ROOT}/local_PCA_MDS_GHSL}"

# ── Genome-wide shared outputs ───────────────────────────────────────────────
export SHARED_DIR="${SHARED_DIR:-${SCRATCH_ROOT}/_shared}"
export SAMPLE_META="${SAMPLE_META:-${SHARED_DIR}/sample_metadata.tsv}"

# ── Logs ─────────────────────────────────────────────────────────────────────
export LOG_DIR="${LOG_DIR:-${SCRATCH_ROOT}/logs}"

# ── Pipeline parameters (post-2026-05-13 slim cleanup) ───────────────────────
export NPC="${NPC:-4}"            # per-window PCA depth
export WINSIZE="${WINSIZE:-100}"  # SNPs per window
export WINSTEP="${WINSTEP:-20}"   # window step (SNPs)

# MDS defaults — chunked_2x removed, only "chromosome" mode active.
# MDS_DIMS reduced from 20 to 5 to match TR_C / GH_C; the base sim_mat
# is built from dmat_focal not from MDS coords so K=5 captures top variance
# with margin without inflating compute.
export MDS_DIMS_DEFAULT="${MDS_DIMS_DEFAULT:-5}"
export Z_THRESH_DEFAULT="${Z_THRESH_DEFAULT:-3}"

# NN-smoothed similarity scales — only the 4 scales actually consumed by
# L1/L2/atlas (was 8; the others had zero readers and 1.2 GB/chrom of disk).
export NN_SIM_SCALES="${NN_SIM_SCALES:-40,80,160,320}"

# =============================================================================
# θπ pipeline parameters (merged from the old 00_theta_config.sh, 2026-05-14)
# =============================================================================
# OUTROOT is the legacy alias TR_A/B/C scripts read for the theta_pi output
# tree. Points at PATH2_ROOT.
export OUTROOT="${OUTROOT:-${PATH2_ROOT}}"

# pestPG scales — coarse for L1/L2 boundary discovery, dense for per-sample
# CUSUM inside called L2 envelopes.
export PESTPG_SCALE="${PESTPG_SCALE:-win10000.step2000}"
export PESTPG_GLOB="${PESTPG_GLOB:-*.${PESTPG_SCALE}.pestPG}"
export COARSE_PESTPG_SCALE="${COARSE_PESTPG_SCALE:-win50000.step10000}"
export DENSE_PESTPG_SCALE="${DENSE_PESTPG_SCALE:-win10000.step2000}"
export SCALE_LABEL="${SCALE_LABEL:-${PESTPG_SCALE}}"
export WIN_BP="${WIN_BP:-10000}"
export STEP_BP="${STEP_BP:-2000}"

# Intermediate per-sample × per-window matrices land in PATH2_ROOT/01_local_pca.
export THETA_TSV_DIR="${THETA_TSV_DIR:-${PATH2_ROOT}/01_local_pca}"

# Optional: dosage-grid window BED (links to Z path's 01_local_pca for
# joint-grid atlas rendering when needed).
export DOSAGE_WIN_BED_DIR="${DOSAGE_WIN_BED_DIR:-${PATH1_LOCAL_PCA}}"

# θπ local-PCA + sim_mat / L1/L2 thresholds.
export LOCAL_PCA_PAD="${LOCAL_PCA_PAD:-1}"
export LOCAL_PCA_NPC="${LOCAL_PCA_NPC:-4}"
export SIM_MAT_METRIC="${SIM_MAT_METRIC:-abs_cosine}"
export CONCORD_THR="${CONCORD_THR:-0.85}"
export MERGE_THR="${MERGE_THR:-0.85}"
export SILHOUETTE_MIN="${SILHOUETTE_MIN:-0.45}"
export ENV_MIN_WINDOWS="${ENV_MIN_WINDOWS:-5}"

# θπ output sub-trees + atlas JSON.
export OUT_PER_WINDOW_DIR="${OUT_PER_WINDOW_DIR:-${PATH2_ROOT}/01_local_pca}"
export OUT_LOCAL_PCA_DIR="${OUT_LOCAL_PCA_DIR:-${PATH2_ROOT}/01_local_pca}"
export OUT_ENVELOPES_DIR="${OUT_ENVELOPES_DIR:-${PATH2_ROOT}/03_per_chrom}"
export JSON_OUT_DIR="${JSON_OUT_DIR:-${PATH2_ROOT}/04_atlas_json}"
export THETA_JSON_SCHEMA_VERSION="${THETA_JSON_SCHEMA_VERSION:-1}"

# θπ R-script env vars (read via Sys.getenv() inside TR_A/B/C/H/I/J).
# Provide sensible defaults so the R scripts don't trip on NA when run
# without first sourcing 00_theta_config.sh.
export THETA_GRID_MODE="${THETA_GRID_MODE:-native}"      # native | dosage
export COHORT_ID="${COHORT_ID:-catfish_226}"             # used by sketch / sparse-edges paths
export SKETCH_DIR="${SKETCH_DIR:-${PATH2_ROOT}/03_dense_sketch}"
export CHROM_LIST="${CHROM_LIST:-${SCRATCH_ROOT}/chr.list}"
export N_CORES="${N_CORES:-${SLURM_CPUS_PER_TASK:-1}}"

# =============================================================================
# GHSL pipeline parameters (2026-05-14)
# =============================================================================
# Input: Clair3 postprocess per-sample TSVs at:
#   ${MODULE_4A_ROOT}/postprocess_results/<chr>/<sample>/all_variants_with_phase.tsv
# The prep step merges these into one <chr>.merged_phased_snps.tsv.gz per chrom.
export MODULE_4A_ROOT="${MODULE_4A_ROOT:-${BASE}/MODULE_4A_SNP_INDEL50_Clair3}"
export GHSL_POSTPROCESS_DIR="${GHSL_POSTPROCESS_DIR:-${MODULE_4A_ROOT}/postprocess_results}"

# GHSL output tree (subdirs under PATH3_ROOT).
export GHSL_PREP_DIR="${GHSL_PREP_DIR:-${PATH3_ROOT}/00_prep}"            # output of prep step
export GHSL_MATRICES_DIR="${GHSL_MATRICES_DIR:-${PATH3_ROOT}/01_matrices}" # output of GH_A
export GHSL_LOCALPCA_DIR="${GHSL_LOCALPCA_DIR:-${PATH3_ROOT}/02_local_pca}" # output of GH_B
export GHSL_PRECOMP_DIR="${GHSL_PRECOMP_DIR:-${PATH3_ROOT}/03_precomp}"    # output of GH_C
export GHSL_L1_DIR="${GHSL_L1_DIR:-${PATH3_ROOT}/04_L1_detect}"            # output of GH_D L1 detect (TR_D)
export GHSL_L1_PLOTS_DIR="${GHSL_L1_PLOTS_DIR:-${PATH3_ROOT}/05_L1_plots}" # output of GH_E plot_L1 (TR_E)
export GHSL_L2_DIR="${GHSL_L2_DIR:-${PATH3_ROOT}/06_L2_detect}"            # output of GH_F L2 detect (TR_F)
export GHSL_L2_PLOTS_DIR="${GHSL_L2_PLOTS_DIR:-${PATH3_ROOT}/07_L2_plots}" # output of GH_G plot_L2 (TR_G)
export GHSL_JSON_DIR="${GHSL_JSON_DIR:-${PATH3_ROOT}/08_atlas_json}"       # output of GH_J atlas JSON

# Multi-scale rolling aggregates (s10 = 10×raw-bp window, s50 = 50×, ...).
export GHSL_SCALES="${GHSL_SCALES:-10,20,30,40,50,100}"

# Default local-PCA scale (used by GH_B and GH_C when --scale not given).
export GHSL_PCA_SCALE_FULL="${GHSL_PCA_SCALE_FULL:-s50}"   # coarse, --mode full
export GHSL_PCA_SCALE_DENSE="${GHSL_PCA_SCALE_DENSE:-s10}" # dense, --mode local

# ── Helper functions ─────────────────────────────────────────────────────────
inv_timestamp() { date '+%F %T'; }
inv_log()  { echo "[$(inv_timestamp)] [INV] $*"; }
inv_err()  { echo "[$(inv_timestamp)] [INV] [ERROR] $*" >&2; }
inv_die()  { inv_err "$@"; exit 1; }

inv_check_file() {
  local f="$1" label="${2:-file}"
  [[ -s "$f" ]] || inv_die "Missing or empty ${label}: ${f}"
}

inv_init_dirs() {
  mkdir -p \
    "${SCRATCH_ROOT}" \
    "${BEAGLE_DIR}" "${DOSAGE_DIR}" "${THETAPI_DIR}" "${GHSL_DIR}" \
    "${PATH1_LOCAL_PCA}" "${PATH1_REGISTRY}" "${PATH1_MDS}" "${PATH1_PRECOMP}" \
    "${PATH1_L1}" "${PATH1_L1_PLOTS}" "${PATH1_L2}" "${PATH1_L2_PLOTS}" "${PATH1_JSON}" \
    "${PATH2_ROOT}" \
    "${THETA_TSV_DIR}" "${OUT_PER_WINDOW_DIR}" "${OUT_LOCAL_PCA_DIR}" \
    "${OUT_ENVELOPES_DIR}" "${JSON_OUT_DIR}" "${SKETCH_DIR}" \
    "${PATH2_ROOT}/precomp" "${PATH2_ROOT}/precomp_dense" \
    "${PATH2_ROOT}/L1_detect" "${PATH2_ROOT}/L1_plots" \
    "${PATH2_ROOT}/L2_detect" "${PATH2_ROOT}/L2_plots" \
    "${PATH3_ROOT}" \
    "${GHSL_PREP_DIR}" "${GHSL_MATRICES_DIR}" "${GHSL_LOCALPCA_DIR}" \
    "${GHSL_PRECOMP_DIR}" \
    "${GHSL_L1_DIR}" "${GHSL_L1_PLOTS_DIR}" \
    "${GHSL_L2_DIR}" "${GHSL_L2_PLOTS_DIR}" \
    "${GHSL_JSON_DIR}" \
    "${SHARED_DIR}" "${LOG_DIR}"
}
