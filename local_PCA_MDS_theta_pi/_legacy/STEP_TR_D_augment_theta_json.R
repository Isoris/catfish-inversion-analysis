#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_D_augment_theta_json.R
# =============================================================================
# Phase 2 / 2f_theta_discovery — page-12 JSON post-processor.
#
# Reads STEP_TR_B's <chr>_phase2_theta.json output and the θπ D17 wrapper TSVs
# (STEP_TR_C output). Adds the theta_d17_envelopes layer, writes the augmented
# JSON in place (or to a separate output if --out_json is given).
#
# Why a separate post-processor instead of modifying STEP_TR_B?
#   STEP_TR_B v4 ships and works. Its envelope detection (|Z|-threshold) is
#   PRIMARY for θπ. The D17 wrapper is an additional detector path that
#   ships its outputs as TSVs. Loading those TSVs into a JSON layer is a
#   pure shape-transformation that runs after both upstreams have produced
#   their outputs. Keeping it separate makes the dependency graph clean:
#
#     STEP_TR_A → STEP_TR_B → STEP_TR_B's JSON output
#                     ↓
#     STEP_TR_A → STEP_TR_B JSON → STEP_TR_C wrapper → D17 TSVs
#                                                          ↓
#     STEP_TR_B's JSON + D17 TSVs → STEP_TR_D → augmented JSON
#
# Symmetric to how export_ghsl_to_json_v3.R consumes both STEP_C04c outputs
# AND STEP_C04d wrapper TSVs in one consolidated emission.
#
# Usage
# -----
#   Rscript STEP_TR_D_augment_theta_json.R \
#     --theta_json    <chr>_phase2_theta.json \
#     --d17_l1_env    <chr>_theta_d17L1_envelopes.tsv \
#     --d17_l2_env    <chr>_theta_d17L2_envelopes.tsv \
#     --d17_l1_bnd    <chr>_theta_d17L1_boundaries.tsv \
#     --d17_l2_bnd    <chr>_theta_d17L2_boundaries.tsv \
#     [--out_json     <path>]    # default: overwrite --theta_json in place
#
# All --d17_* flags are individually optional; missing files just produce
# an empty list for that sub-layer.
#
# If neither d17 TSV is provided, the layer is omitted entirely (the JSON
# is written back unchanged).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- CLI --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}

THETA_JSON  <- get_arg("--theta_json")
D17_L1_ENV  <- get_arg("--d17_l1_env",  NULL)
D17_L2_ENV  <- get_arg("--d17_l2_env",  NULL)
D17_L1_BND  <- get_arg("--d17_l1_bnd",  NULL)
D17_L2_BND  <- get_arg("--d17_l2_bnd",  NULL)
OUT_JSON    <- get_arg("--out_json",    THETA_JSON)

if (is.null(THETA_JSON) || !file.exists(THETA_JSON)) {
  stop("[TR_D] FATAL: --theta_json <chr>_phase2_theta.json is required")
}

cat("================================================================\n")
cat("[TR_D] Augmenting θπ JSON with D17 envelopes\n")
cat("================================================================\n")
cat("[TR_D] input  : ", THETA_JSON, "\n", sep = "")
cat("[TR_D] output : ", OUT_JSON,   "\n", sep = "")
cat("[TR_D] d17 L1 env: ", D17_L1_ENV %||% "(skipped)", "\n", sep = "")
cat("[TR_D] d17 L2 env: ", D17_L2_ENV %||% "(skipped)", "\n", sep = "")
cat("[TR_D] d17 L1 bnd: ", D17_L1_BND %||% "(skipped)", "\n", sep = "")
cat("[TR_D] d17 L2 bnd: ", D17_L2_BND %||% "(skipped)", "\n", sep = "")

# ---- Helpers ----------------------------------------------------------------
round4 <- function(x) round(x, 4)

load_d17_tsv <- function(path, kind) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  cat("[TR_D]   loading ", kind, ": ", basename(path), "\n", sep = "")
  fread(path)
}

d17_env_to_records <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(list())
  lapply(seq_len(nrow(dt)), function(i) {
    r <- as.list(dt[i])
    list(
      candidate_id = as.character(r$candidate_id),
      start_w      = as.integer(r$start_w),
      end_w        = as.integer(r$end_w),
      start_bp     = as.integer(r$start_bp),
      end_bp       = as.integer(r$end_bp),
      n_windows    = as.integer(r$n_windows),
      mean_sim     = round4(as.numeric(r$mean_sim    %||% NA)),
      density_p70  = round4(as.numeric(r$density_p70 %||% NA)),
      status       = as.character(r$status %||% NA),
      parent_l1_id = if ("parent_l1_id" %in% names(dt))
                       as.character(r$parent_l1_id) else NULL
    )
  })
}

d17_bnd_to_records <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(list())
  lapply(seq_len(nrow(dt)), function(i) {
    r <- as.list(dt[i])
    out <- list(
      boundary_idx      = as.character(r$boundary_idx),
      boundary_w        = as.integer(r$boundary_w),
      boundary_bp       = as.integer(r$boundary_bp),
      validation_status = as.character(r$validation_status),
      boundary_score    = round4(as.numeric(r$boundary_score))
    )
    if ("grow_max_z"   %in% names(dt))
      out$grow_max_z   <- round4(as.numeric(r$grow_max_z))
    if ("parent_l1_id" %in% names(dt))
      out$parent_l1_id <- as.character(r$parent_l1_id)
    out
  })
}

# ---- Load existing JSON -----------------------------------------------------
cat("\n[TR_D] Loading existing θπ JSON\n")
js <- fromJSON(THETA_JSON, simplifyVector = FALSE)
cat("[TR_D]   chrom: ", js$chrom %||% "(unknown)",
    "  layers: ", paste(unlist(js$`_layers_present`), collapse = ", "),
    "\n", sep = "")

# ---- Load D17 TSVs ----------------------------------------------------------
d17_l1_env <- load_d17_tsv(D17_L1_ENV, "L1 envelopes")
d17_l2_env <- load_d17_tsv(D17_L2_ENV, "L2 envelopes")
d17_l1_bnd <- load_d17_tsv(D17_L1_BND, "L1 boundaries")
d17_l2_bnd <- load_d17_tsv(D17_L2_BND, "L2 boundaries")

n_inputs <- sum(c(!is.null(d17_l1_env), !is.null(d17_l2_env),
                  !is.null(d17_l1_bnd), !is.null(d17_l2_bnd)))
if (n_inputs == 0) {
  cat("[TR_D] No D17 inputs provided — passing JSON through unchanged\n")
} else {
  cat("[TR_D] Building theta_d17_envelopes layer\n")
  d17_layer <- list(
    schema_version = 2L,
    layer          = "theta_d17_envelopes",
    chrom          = js$chrom,
    source         = "STEP_TR_C wrapper → D17 cross-block boundary detector",
    detector       = "STEP_D17_multipass_L1_only_v7.R + STEP_D17_multipass_L2_v8.R",
    l1_envelopes   = d17_env_to_records(d17_l1_env),
    l2_envelopes   = d17_env_to_records(d17_l2_env),
    l1_boundaries  = d17_bnd_to_records(d17_l1_bnd),
    l2_boundaries  = d17_bnd_to_records(d17_l2_bnd)
  )
  cat("[TR_D]   theta_d17_envelopes: L1=", length(d17_layer$l1_envelopes),
      " L2=",  length(d17_layer$l2_envelopes),
      " bL1=", length(d17_layer$l1_boundaries),
      " bL2=", length(d17_layer$l2_boundaries), "\n", sep = "")

  js$theta_d17_envelopes <- d17_layer

  # Update layers_present
  layers <- unlist(js$`_layers_present`)
  if (!"theta_d17_envelopes" %in% layers) {
    layers <- c(layers, "theta_d17_envelopes")
  }
  js$`_layers_present` <- layers
  js$`_augmented_at` <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  js$`_augmenter`    <- "STEP_TR_D_augment_theta_json.R"
}

# ---- Write ------------------------------------------------------------------
cat("\n[TR_D] Writing augmented JSON\n")
write_json(js, OUT_JSON,
           auto_unbox = TRUE, na = "null", pretty = FALSE, digits = NA)
fi <- file.info(OUT_JSON)
cat("[TR_D] DONE — ", OUT_JSON, " (",
    round(fi$size / 1024 / 1024, 2), " MB, ",
    if (n_inputs > 0) "augmented" else "unchanged", ")\n", sep = "")
