#!/usr/bin/env bash
# =============================================================================
# 00_inversion_config.sh — v10 (2026-05-16, harmonized-layout cleanup)
#
# Source with:
#   set -a; source "${CONFIG}"; set +a   # auto-export
#
# Every variable below uses the ${VAR:-default} pattern so anything set in
# the environment BEFORE sourcing this file wins. This prevents the previous
# trap where SCRIPT_DIR got clobbered when a launcher had already captured
# its own location.
#
# Harmonized output layout (v10): all three pipelines (z, theta_pi, GHSL)
# share the SAME subdirectory numbering under their respective PATH<N>_ROOT.
# The shared L1/L2 detect+plot scripts in _shared/ are parameterized by
# --path1/--path2/--path3 and emit into the same numbered slots.
#
#   <PATH_ROOT>/
#     00_prep/            ← optional prep step (GHSL prep; theta_pi/z skip)
#     01_local_pca/       ← per-window local PCA (ZO_C / TR_B / GH_B)
#     02_dense_registry/  ← genome-wide window registry + merge
#                            (ZO_D / TR_D / GH_D)
#     03_mds/             ← per-chrom MDS results (ZO_E / TR_C / GH_C)
#                            One <chr>.mds_perchr.rds per chrom. NO tmp/.
#                            NO whole-genome combined RDS.
#     04_precomp/         ← per-chrom precomp.rds + sim_mats/
#                            (ZO_G / TR_C-full / GH_C-full). NO tmp/.
#                            NO whole-genome combined RDS.
#     04_precomp_dense/   ← optional dense-scale precomp for CUSUM/carriers
#                            (theta_pi / GHSL only)
#     05_L1/              ← L1 envelopes + boundaries (shared L1 detect)
#     06_L1_plots/        ← optional L1 plots (gated by MAKE_L1_PLOTS=true)
#     07_L2/              ← L2 envelopes + boundaries (shared L2 detect)
#     08_L2_plots/        ← optional L2 plots (gated by MAKE_L2_PLOTS=true)
#     08_carriers/        ← optional carrier assignments (theta_pi/GHSL)
#     09_atlas_json/      ← per-chrom atlas JSONs (harmonized schema v4)
#     09_cusum/           ← optional per-sample CUSUM (theta_pi/GHSL)
#     plot/               ← legacy alias; if MAKE_*_PLOTS=true and
#                            <PATH_ROOT>/plot exists, plots go there.
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
export SCRIPT_DIR_SHARED="${SCRIPT_DIR_SHARED:-${REPO_ROOT}/_shared}"

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

# =============================================================================
# Harmonized subdir slots — shared by ALL three pipelines (z / theta_pi / GHSL)
# Each path is a "slot" that resolves under PATH<N>_ROOT. The slot names are
# the SAME across pipelines so the shared L1/L2 scripts and the JSON exporter
# can locate inputs via <PATH_ROOT>/<SLOT>.
# =============================================================================
export SLOT_PREP="${SLOT_PREP:-00_prep}"
export SLOT_LOCAL_PCA="${SLOT_LOCAL_PCA:-01_local_pca}"
export SLOT_REGISTRY="${SLOT_REGISTRY:-02_dense_registry}"
export SLOT_MDS="${SLOT_MDS:-03_mds}"
export SLOT_PRECOMP="${SLOT_PRECOMP:-04_precomp}"
export SLOT_PRECOMP_DENSE="${SLOT_PRECOMP_DENSE:-04_precomp_dense}"
export SLOT_L1="${SLOT_L1:-05_L1}"
export SLOT_L1_PLOTS="${SLOT_L1_PLOTS:-06_L1_plots}"
export SLOT_L2="${SLOT_L2:-07_L2}"
export SLOT_L2_PLOTS="${SLOT_L2_PLOTS:-08_L2_plots}"
export SLOT_CARRIERS="${SLOT_CARRIERS:-08_carriers}"
export SLOT_JSON="${SLOT_JSON:-09_atlas_json}"
export SLOT_CUSUM="${SLOT_CUSUM:-09_cusum}"

# ── Plot toggles (facultative, default ON for back-compat) ───────────────────
# MAKE_L1_PLOTS / MAKE_L2_PLOTS control whether the shared L1/L2 detect
# scripts launch their downstream plotting passes. Set to "false" / "0" /
# "no" to skip. Independent toggles so you can skip L2 plots (slowest) but
# keep L1 plots, etc. Honored by:
#   _shared/STEP_05_L1_detect.R         (skips writing to SLOT_L1_PLOTS)
#   _shared/STEP_06_L1_plot.R           (refuses to run if false)
#   _shared/STEP_07_L2_detect.R         (skips writing to SLOT_L2_PLOTS)
#   _shared/STEP_08_L2_plot.R           (refuses to run if false)
export MAKE_L1_PLOTS="${MAKE_L1_PLOTS:-true}"
export MAKE_L2_PLOTS="${MAKE_L2_PLOTS:-true}"

# ── Path 1 (z-blocks / dosage) output tree ───────────────────────────────────
export PATH1_ROOT="${PATH1_ROOT:-${SCRATCH_ROOT}/local_PCA_MDS_z}"
export PATH1_LOCAL_PCA="${PATH1_LOCAL_PCA:-${PATH1_ROOT}/${SLOT_LOCAL_PCA}}"
export PATH1_REGISTRY="${PATH1_REGISTRY:-${PATH1_ROOT}/${SLOT_REGISTRY}}"
export PATH1_MDS="${PATH1_MDS:-${PATH1_ROOT}/${SLOT_MDS}}"

# v10: 04_precomp (was 04_precomp_v2/precomp/). The nested precomp/precomp/
# layer is gone — files land directly in <PATH1_ROOT>/04_precomp/.
# PRECOMP_DIR alias kept for back-compat with any old caller.
export PATH1_PRECOMP="${PATH1_PRECOMP:-${PATH1_ROOT}/${SLOT_PRECOMP}}"
export PRECOMP_DIR="${PRECOMP_DIR:-${PATH1_PRECOMP}}"
export PATH1_L1="${PATH1_L1:-${PATH1_ROOT}/${SLOT_L1}}"
export PATH1_L1_PLOTS="${PATH1_L1_PLOTS:-${PATH1_ROOT}/${SLOT_L1_PLOTS}}"
export PATH1_L2="${PATH1_L2:-${PATH1_ROOT}/${SLOT_L2}}"
export PATH1_L2_PLOTS="${PATH1_L2_PLOTS:-${PATH1_ROOT}/${SLOT_L2_PLOTS}}"
export PATH1_JSON="${PATH1_JSON:-${PATH1_ROOT}/${SLOT_JSON}}"

# Path 1 MDS prefix (ZO_E writes <prefix>.mds_metadata.tsv; ZO_G uses dirname
# to locate the per-chrom MDS dir.) v10: MDS files now live directly in
# 03_mds/<chr>.mds_perchr.rds — no tmp/ subdir.
export MDS_PREFIX_BASENAME="${MDS_PREFIX_BASENAME:-inversion_localpca}"
export PATH1_MDS_PREFIX="${PATH1_MDS_PREFIX:-${PATH1_MDS}/${MDS_PREFIX_BASENAME}}"
export MDS_PREFIX="${MDS_PREFIX:-${PATH1_MDS_PREFIX}}"

# ── Path 2 (theta_pi) output tree ────────────────────────────────────────────
export PATH2_ROOT="${PATH2_ROOT:-${SCRATCH_ROOT}/local_PCA_MDS_theta_pi}"
export PATH2_LOCAL_PCA="${PATH2_LOCAL_PCA:-${PATH2_ROOT}/${SLOT_LOCAL_PCA}}"
export PATH2_REGISTRY="${PATH2_REGISTRY:-${PATH2_ROOT}/${SLOT_REGISTRY}}"
export PATH2_MDS="${PATH2_MDS:-${PATH2_ROOT}/${SLOT_MDS}}"
export PATH2_PRECOMP="${PATH2_PRECOMP:-${PATH2_ROOT}/${SLOT_PRECOMP}}"
export PATH2_PRECOMP_DENSE="${PATH2_PRECOMP_DENSE:-${PATH2_ROOT}/${SLOT_PRECOMP_DENSE}}"
export PATH2_L1="${PATH2_L1:-${PATH2_ROOT}/${SLOT_L1}}"
export PATH2_L1_PLOTS="${PATH2_L1_PLOTS:-${PATH2_ROOT}/${SLOT_L1_PLOTS}}"
export PATH2_L2="${PATH2_L2:-${PATH2_ROOT}/${SLOT_L2}}"
export PATH2_L2_PLOTS="${PATH2_L2_PLOTS:-${PATH2_ROOT}/${SLOT_L2_PLOTS}}"
export PATH2_CARRIERS="${PATH2_CARRIERS:-${PATH2_ROOT}/${SLOT_CARRIERS}}"
export PATH2_CUSUM="${PATH2_CUSUM:-${PATH2_ROOT}/${SLOT_CUSUM}}"
export PATH2_JSON="${PATH2_JSON:-${PATH2_ROOT}/${SLOT_JSON}}"

# ── Path 3 (GHSL) output tree ────────────────────────────────────────────────
export PATH3_ROOT="${PATH3_ROOT:-${SCRATCH_ROOT}/local_PCA_MDS_GHSL}"
export PATH3_PREP="${PATH3_PREP:-${PATH3_ROOT}/${SLOT_PREP}}"
export PATH3_LOCAL_PCA="${PATH3_LOCAL_PCA:-${PATH3_ROOT}/${SLOT_LOCAL_PCA}}"
export PATH3_REGISTRY="${PATH3_REGISTRY:-${PATH3_ROOT}/${SLOT_REGISTRY}}"
export PATH3_MDS="${PATH3_MDS:-${PATH3_ROOT}/${SLOT_MDS}}"
export PATH3_PRECOMP="${PATH3_PRECOMP:-${PATH3_ROOT}/${SLOT_PRECOMP}}"
export PATH3_PRECOMP_DENSE="${PATH3_PRECOMP_DENSE:-${PATH3_ROOT}/${SLOT_PRECOMP_DENSE}}"
export PATH3_L1="${PATH3_L1:-${PATH3_ROOT}/${SLOT_L1}}"
export PATH3_L1_PLOTS="${PATH3_L1_PLOTS:-${PATH3_ROOT}/${SLOT_L1_PLOTS}}"
export PATH3_L2="${PATH3_L2:-${PATH3_ROOT}/${SLOT_L2}}"
export PATH3_L2_PLOTS="${PATH3_L2_PLOTS:-${PATH3_ROOT}/${SLOT_L2_PLOTS}}"
export PATH3_CARRIERS="${PATH3_CARRIERS:-${PATH3_ROOT}/${SLOT_CARRIERS}}"
export PATH3_CUSUM="${PATH3_CUSUM:-${PATH3_ROOT}/${SLOT_CUSUM}}"
export PATH3_JSON="${PATH3_JSON:-${PATH3_ROOT}/${SLOT_JSON}}"

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

# Intermediate per-sample × per-window matrices land in PATH2_LOCAL_PCA.
export THETA_TSV_DIR="${THETA_TSV_DIR:-${PATH2_LOCAL_PCA}}"

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

# θπ output sub-trees — now aliased to the harmonized slots.
export OUT_PER_WINDOW_DIR="${OUT_PER_WINDOW_DIR:-${PATH2_LOCAL_PCA}}"
export OUT_LOCAL_PCA_DIR="${OUT_LOCAL_PCA_DIR:-${PATH2_LOCAL_PCA}}"
export OUT_ENVELOPES_DIR="${OUT_ENVELOPES_DIR:-${PATH2_L1}}"   # L1 envelopes
export JSON_OUT_DIR="${JSON_OUT_DIR:-${PATH2_JSON}}"
# Harmonized atlas JSON schema v4 (unified across z / theta_pi / GHSL).
# Use the converter (_shared/convert_atlas_json.R) to migrate older outputs.
export THETA_JSON_SCHEMA_VERSION="${THETA_JSON_SCHEMA_VERSION:-4}"
export ATLAS_JSON_SCHEMA_VERSION="${ATLAS_JSON_SCHEMA_VERSION:-4}"

# θπ R-script env vars (read via Sys.getenv() inside TR_A/B/C/H/I/J).
# Provide sensible defaults so the R scripts don't trip on NA when run
# without first sourcing 00_theta_config.sh.
export THETA_GRID_MODE="${THETA_GRID_MODE:-native}"      # native | dosage
export COHORT_ID="${COHORT_ID:-catfish_226}"             # used by sketch / sparse-edges paths
# Sketch dir co-located with PATH2_PRECOMP_DENSE for the dense-mode TR_C run.
export SKETCH_DIR="${SKETCH_DIR:-${PATH2_PRECOMP_DENSE}/sketch}"
export CHROM_LIST="${CHROM_LIST:-${SCRATCH_ROOT}/chr.list}"
export N_CORES="${N_CORES:-${SLURM_CPUS_PER_TASK:-1}}"

# =============================================================================
# GHSL pipeline parameters (2026-05-14, harmonized 2026-05-16)
# =============================================================================
# Input: Clair3 postprocess per-sample TSVs at:
#   ${MODULE_4A_ROOT}/postprocess_results/<chr>/<sample>/all_variants_with_phase.tsv
# The prep step merges these into one <chr>.merged_phased_snps.tsv.gz per chrom.
export MODULE_4A_ROOT="${MODULE_4A_ROOT:-${BASE}/MODULE_4A_SNP_INDEL50_Clair3}"
export GHSL_POSTPROCESS_DIR="${GHSL_POSTPROCESS_DIR:-${MODULE_4A_ROOT}/postprocess_results}"

# Legacy GHSL_*_DIR aliases now point at the harmonized PATH3_* slots. The
# numbering matches z / theta_pi exactly — there is no 00_prep/01_matrices/
# 02_local_pca split anymore; GHSL prep lives in 00_prep, GH_A writes the
# per-window PCA inputs into 01_local_pca alongside GH_B's output, GH_C
# writes per-chrom MDS+precomp into 03_mds and 04_precomp, etc.
export GHSL_PREP_DIR="${GHSL_PREP_DIR:-${PATH3_PREP}}"            # output of prep step
export GHSL_MATRICES_DIR="${GHSL_MATRICES_DIR:-${PATH3_LOCAL_PCA}}" # output of GH_A
export GHSL_LOCALPCA_DIR="${GHSL_LOCALPCA_DIR:-${PATH3_LOCAL_PCA}}" # output of GH_B
export GHSL_MDS_DIR="${GHSL_MDS_DIR:-${PATH3_MDS}}"                # output of GH_C (mds_perchr)
export GHSL_PRECOMP_DIR="${GHSL_PRECOMP_DIR:-${PATH3_PRECOMP}}"           # output of GH_C (coarse, mode=full)
export GHSL_PRECOMP_DENSE_DIR="${GHSL_PRECOMP_DENSE_DIR:-${PATH3_PRECOMP_DENSE}}" # output of GH_C (dense, mode=local) — for carriers/CUSUM
export GHSL_L1_DIR="${GHSL_L1_DIR:-${PATH3_L1}}"                          # output of shared L1 detect
export GHSL_L1_PLOTS_DIR="${GHSL_L1_PLOTS_DIR:-${PATH3_L1_PLOTS}}"        # output of shared L1 plot
export GHSL_L2_DIR="${GHSL_L2_DIR:-${PATH3_L2}}"                          # output of shared L2 detect
export GHSL_L2_PLOTS_DIR="${GHSL_L2_PLOTS_DIR:-${PATH3_L2_PLOTS}}"        # output of shared L2 plot
export GHSL_CARRIERS_DIR="${GHSL_CARRIERS_DIR:-${PATH3_CARRIERS}}"        # output of carriers
export GHSL_CUSUM_DIR="${GHSL_CUSUM_DIR:-${PATH3_CUSUM}}"                 # output of cusum
export GHSL_JSON_DIR="${GHSL_JSON_DIR:-${PATH3_JSON}}"                    # output of GH_J atlas JSON

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

# Truthiness helper for MAKE_L1_PLOTS / MAKE_L2_PLOTS env vars. Accepts
# true|1|yes|y|on (case-insensitive). Everything else is false.
inv_is_true() {
  case "${1,,}" in
    1|true|yes|y|on) return 0 ;;
    *)               return 1 ;;
  esac
}

inv_init_dirs() {
  mkdir -p \
    "${SCRATCH_ROOT}" \
    "${BEAGLE_DIR}" "${DOSAGE_DIR}" "${THETAPI_DIR}" "${GHSL_DIR}" \
    "${PATH1_LOCAL_PCA}" "${PATH1_REGISTRY}" "${PATH1_MDS}" "${PATH1_PRECOMP}" \
    "${PATH1_L1}" "${PATH1_L1_PLOTS}" "${PATH1_L2}" "${PATH1_L2_PLOTS}" "${PATH1_JSON}" \
    "${PATH2_ROOT}" \
    "${PATH2_LOCAL_PCA}" "${PATH2_REGISTRY}" "${PATH2_MDS}" \
    "${PATH2_PRECOMP}" "${PATH2_PRECOMP_DENSE}" \
    "${PATH2_L1}" "${PATH2_L1_PLOTS}" "${PATH2_L2}" "${PATH2_L2_PLOTS}" \
    "${PATH2_CARRIERS}" "${PATH2_CUSUM}" "${PATH2_JSON}" \
    "${SKETCH_DIR}" \
    "${PATH3_ROOT}" \
    "${PATH3_PREP}" "${PATH3_LOCAL_PCA}" "${PATH3_REGISTRY}" "${PATH3_MDS}" \
    "${PATH3_PRECOMP}" "${PATH3_PRECOMP_DENSE}" \
    "${PATH3_L1}" "${PATH3_L1_PLOTS}" "${PATH3_L2}" "${PATH3_L2_PLOTS}" \
    "${PATH3_CARRIERS}" "${PATH3_CUSUM}" "${PATH3_JSON}" \
    "${SHARED_DIR}" "${LOG_DIR}"
}
