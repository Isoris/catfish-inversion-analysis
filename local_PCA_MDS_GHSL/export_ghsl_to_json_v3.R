#!/usr/bin/env Rscript
# =============================================================================
# export_ghsl_to_json_v3.R
# =============================================================================
# Phase 2 / 2e_ghsl_discovery — page-3 atlas JSON exporter.
#
# Consolidates GHSL outputs from FOUR sources into one page-3 JSON:
#
#   1. STEP_C04   v6 heavy:    <chr>.ghsl_v6.annot.rds         → tracks +
#                              <chr>.ghsl_v6.per_sample.rds      ghsl_panel +
#                              <chr>.ghsl_v6.karyotypes.rds      ghsl_kstripes +
#                                                                ghsl_karyotype_runs +
#                                                                ghsl_envelopes
#                                                                (PASS-runs)
#   2. STEP_C04c  local-PCA:   <chr>.ghsl_v6_localpca.rds      → ghsl_local_pca +
#                                                                ghsl_secondary_envelopes
#                                                                (|Z|-threshold)
#   3. STEP_C04d  D17 wrapper: <chr>_ghsl_d17L1_envelopes.tsv  → ghsl_d17_envelopes
#                              <chr>_ghsl_d17L2_envelopes.tsv    (boundary-detector)
#                              <chr>_ghsl_d17L1_boundaries.tsv
#                              <chr>_ghsl_d17L2_boundaries.tsv
#   4. (optional) sample meta: <samples>.tsv (ind, cga, ancestry)
#
# Output: <out_dir>/<chr>/<chr>_phase2_ghsl.json (renamed from legacy
# _phase4a_ghsl.json per the Phase 2 reshape; back-compat handled in the
# atlas's loader).
#
# Layer inventory:
#   Layer name                  | Source                          | Status
#   ---------------------------- | -------------------------------- | -------
#   tracks                       | STEP_C04 annot RDS aggregates    | existing
#   ghsl_panel                   | STEP_C04 per_sample RDS          | existing
#   ghsl_kstripes                | computed K=2..6 stripe assigns   | existing
#   ghsl_karyotype_runs          | STEP_C04 karyotypes RDS          | existing
#   ghsl_local_pca               | STEP_C04c localpca RDS           | NEW turn 2/3
#   ghsl_envelopes (PRIMARY)     | STEP_C04 annot PASS-runs         | NEW turn 3
#   ghsl_secondary_envelopes     | STEP_C04c localpca z_profile     | NEW turn 3
#   ghsl_d17_envelopes           | STEP_C04d D17 TSVs               | NEW turn 5
#
# Asymmetry vs θπ: page 12 (STEP_TR_B emits its own JSON directly) reads only
# theta_pi_per_window + theta_pi_local_pca + theta_pi_envelopes (|Z|) and
# now also theta_d17_envelopes (parallel D17 wrapper output). Same multi-
# layer architecture, different layer set.
#
# Usage
# -----
#   Rscript export_ghsl_to_json_v3.R \
#     --chrom         C_gar_LG28 \
#     --annot_rds     <ghsl_v6_dir>/<chrom>.ghsl_v6.annot.rds \
#     --persamp_rds   <ghsl_v6_dir>/<chrom>.ghsl_v6.per_sample.rds \
#     --karyo_rds     <ghsl_v6_dir>/<chrom>.ghsl_v6.karyotypes.rds \
#     --localpca_rds  <localpca_dir>/<chrom>.ghsl_v6_localpca.rds \
#     --d17_l1_env    <d17_dir>/<chrom>_ghsl_d17L1_envelopes.tsv \
#     --d17_l2_env    <d17_dir>/<chrom>_ghsl_d17L2_envelopes.tsv \
#     --d17_l1_bnd    <d17_dir>/<chrom>_ghsl_d17L1_boundaries.tsv \
#     --d17_l2_bnd    <d17_dir>/<chrom>_ghsl_d17L2_boundaries.tsv \
#     --samples       samples.tsv \
#     --out_dir       <json_out_dir> \
#     [--primary_scale s50]    \
#     [--max_k 6]              \
#     [--sim_mat_thumb_n 200]
#
# All --d17_* flags are OPTIONAL: if missing, ghsl_d17_envelopes layer is
# omitted from the JSON. Same for --localpca_rds (turn 2 may not have run
# yet on all chromosomes when running the exporter genome-wide).
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

CHROM         <- get_arg("--chrom")
ANNOT_RDS     <- get_arg("--annot_rds")
PERSAMP_RDS   <- get_arg("--persamp_rds")
KARYO_RDS     <- get_arg("--karyo_rds")
LOCALPCA_RDS  <- get_arg("--localpca_rds", NULL)         # optional — turn 2 may not have run
D17_L1_ENV    <- get_arg("--d17_l1_env",   NULL)         # optional — turn 4 may not have run
D17_L2_ENV    <- get_arg("--d17_l2_env",   NULL)
D17_L1_BND    <- get_arg("--d17_l1_bnd",   NULL)
D17_L2_BND    <- get_arg("--d17_l2_bnd",   NULL)
SAMPLES_TSV   <- get_arg("--samples",      NULL)
OUT_DIR       <- get_arg("--out_dir",      ".")
PRIMARY_SCALE <- get_arg("--primary_scale", "s50")
MAX_K         <- as.integer(get_arg("--max_k", "6"))
SIM_THUMB_N   <- as.integer(get_arg("--sim_mat_thumb_n", "200"))

stopifnot(!is.null(CHROM), !is.null(ANNOT_RDS), file.exists(ANNOT_RDS))
stopifnot(!is.null(PERSAMP_RDS), file.exists(PERSAMP_RDS))
stopifnot(!is.null(KARYO_RDS),   file.exists(KARYO_RDS))

dir.create(file.path(OUT_DIR, CHROM), recursive = TRUE, showWarnings = FALSE)
OUT_JSON <- file.path(OUT_DIR, CHROM, paste0(CHROM, "_phase2_ghsl.json"))

cat("================================================================\n")
cat("[GHSL_EXPORT v3] Building page-3 JSON for ", CHROM, "\n", sep = "")
cat("================================================================\n")
cat("[GHSL_EXPORT v3] annot:    ", ANNOT_RDS,    "\n", sep = "")
cat("[GHSL_EXPORT v3] persamp:  ", PERSAMP_RDS,  "\n", sep = "")
cat("[GHSL_EXPORT v3] karyo:    ", KARYO_RDS,    "\n", sep = "")
cat("[GHSL_EXPORT v3] localpca: ", LOCALPCA_RDS %||% "(skipped)", "\n", sep = "")
cat("[GHSL_EXPORT v3] d17 L1:   ", D17_L1_ENV   %||% "(skipped)", "\n", sep = "")
cat("[GHSL_EXPORT v3] d17 L2:   ", D17_L2_ENV   %||% "(skipped)", "\n", sep = "")
cat("[GHSL_EXPORT v3] samples:  ", SAMPLES_TSV  %||% "(none)",    "\n", sep = "")
cat("[GHSL_EXPORT v3] output:   ", OUT_JSON,     "\n", sep = "")

# ---- Helpers ----------------------------------------------------------------
clean_numeric <- function(x, digits = 6) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  round(x, digits)
}
round4 <- function(x) round(x, 4)

# ---- Load STEP_C04 v6 annot RDS ---------------------------------------------
# Carries per-window stats: ghsl_v6_score, ghsl_v6_status (PASS/WEAK/FAIL),
# div_median, rank_stability, plus window_info coordinates. STEP_C04 produces
# this; don't refactor.
cat("\n[GHSL_EXPORT v3] Loading annot RDS\n")
annot <- readRDS(ANNOT_RDS)
if (is.data.frame(annot)) {
  annot_dt <- as.data.table(annot)
} else if (is.list(annot) && !is.null(annot$dt)) {
  annot_dt <- as.data.table(annot$dt)
} else if (is.list(annot) && !is.null(annot$annot)) {
  annot_dt <- as.data.table(annot$annot)
} else {
  stop("[GHSL_EXPORT v3] annot RDS structure not recognized")
}
n_windows <- nrow(annot_dt)
cat("[GHSL_EXPORT v3]   annot: ", n_windows, " windows\n", sep = "")

# Probe column names defensively (production v6 may rename across versions)
status_col <- intersect(c("ghsl_v6_status", "ghsl_status", "status"),
                        names(annot_dt))[1]
score_col  <- intersect(c("ghsl_v6_score", "ghsl_score", "score"),
                        names(annot_dt))[1]
if (is.na(status_col)) {
  cat("[GHSL_EXPORT v3]   warning: no status column in annot RDS — ",
      "envelopes layer will be empty\n", sep = "")
}

# ---- Load STEP_C04 v6 per_sample RDS ----------------------------------------
# Carries div_roll[scale][n_samples × n_windows] — the rolling-smoothed
# divergence panel rendered by the page-3-bis K-stripe heatmap and used by
# the atlas's lines-color-mode picker for per-sample GHSL coloring.
cat("[GHSL_EXPORT v3] Loading per_sample RDS\n")
ps <- readRDS(PERSAMP_RDS)
n_samples    <- ps$n_samples %||% length(ps$sample_names)
sample_names <- ps$sample_names

# ---- Load STEP_C04 v6 karyotypes RDS ----------------------------------------
cat("[GHSL_EXPORT v3] Loading karyotypes RDS\n")
karyo <- readRDS(KARYO_RDS)

# ---- Sample identity layer --------------------------------------------------
sample_meta <- data.table(ind = sample_names, cga = sample_names,
                          ancestry = "unknown")
if (!is.null(SAMPLES_TSV) && file.exists(SAMPLES_TSV)) {
  cat("[GHSL_EXPORT v3] Loading sample metadata: ", SAMPLES_TSV, "\n", sep = "")
  smeta <- fread(SAMPLES_TSV, header = TRUE)
  setnames(smeta, tolower(names(smeta)))
  if ("ind" %in% names(smeta)) {
    sample_meta <- merge(sample_meta[, .(ind)], smeta, by = "ind",
                         all.x = TRUE, sort = FALSE)
  } else if ("cga" %in% names(smeta)) {
    sample_meta[, cga := ind]
    sample_meta <- merge(sample_meta, smeta, by = "cga",
                         all.x = TRUE, sort = FALSE)
  }
  if (!"ancestry" %in% names(sample_meta)) sample_meta[, ancestry := "unknown"]
  if (!"cga"      %in% names(sample_meta)) sample_meta[, cga      := ind]
  sample_meta[is.na(cga),      cga      := ind]
  sample_meta[is.na(ancestry), ancestry := "unknown"]
}
sample_meta <- sample_meta[match(sample_names, ind)]

# =============================================================================
# Layer: tracks
# =============================================================================
# Per-window aggregates additive into the cross-page chromosome track strip.
# Same shape as the dosage scrubber's tracks.
build_tracks_layer <- function(annot_dt) {
  tracks <- list()
  pos_bp_centers <- as.integer(round(
    (annot_dt$start_bp + annot_dt$end_bp) / 2))
  add_track <- function(label, values, digits = 4) {
    finite <- values[is.finite(values)]
    if (length(finite) == 0L) return(invisible())
    tracks[[label]] <<- list(
      values = clean_numeric(values, digits),
      pos_bp = pos_bp_centers,
      min    = round(min(finite), digits),
      max    = round(max(finite), digits),
      mean   = round(mean(finite), digits)
    )
  }
  if (!is.null(score_col)  && score_col  %in% names(annot_dt))
    add_track("ghsl_score",        annot_dt[[score_col]])
  if ("div_median" %in% names(annot_dt))
    add_track("ghsl_div_median",   annot_dt$div_median, digits = 6)
  if ("rank_stability" %in% names(annot_dt))
    add_track("ghsl_rank_stability", annot_dt$rank_stability)
  tracks
}
tracks_layer <- build_tracks_layer(annot_dt)
cat("[GHSL_EXPORT v3]   tracks: ", length(tracks_layer),
    " per-window aggregates\n", sep = "")

# =============================================================================
# Layer: ghsl_panel  (existing v2 contract — div_roll matrices per scale)
# =============================================================================
build_panel_layer <- function(ps, primary_scale) {
  div_roll <- ps$div_roll %||% ps$rolling
  if (is.null(div_roll)) {
    warning("[GHSL_EXPORT v3] per_sample RDS has no div_roll/rolling matrices")
    return(NULL)
  }
  # Pre-aggregate cohort-mean divergence at the primary scale (used by
  # the atlas's track derivation; cheap and avoids the browser having to
  # average across samples each render).
  primary_mat <- div_roll[[primary_scale]] %||%
                 div_roll[[paste0("s", primary_scale)]]
  div_median <- if (!is.null(primary_mat))
    apply(primary_mat, 2, function(v) median(v, na.rm = TRUE))
  else NULL

  list(
    schema_version = 2L,
    primary_scale  = primary_scale,
    n_samples      = as.integer(n_samples),
    n_windows      = as.integer(n_windows),
    sample_names   = sample_names,
    start_bp       = as.integer(annot_dt$start_bp),
    end_bp         = as.integer(annot_dt$end_bp),
    div_roll       = lapply(div_roll, function(m) {
      m_cleaned <- m
      m_cleaned[!is.finite(m_cleaned)] <- NA_real_
      round(m_cleaned, 6)
    }),
    div_median     = if (!is.null(div_median)) clean_numeric(div_median, 6) else NULL
  )
}
panel_layer <- build_panel_layer(ps, PRIMARY_SCALE)
cat("[GHSL_EXPORT v3]   ghsl_panel: ", n_samples, " samples × ",
    n_windows, " windows × ",
    if (!is.null(panel_layer)) length(panel_layer$div_roll) else 0,
    " scales\n", sep = "")

# =============================================================================
# Layer: ghsl_kstripes  (existing v2 contract — K=2..MAX_K stripe assignments)
# =============================================================================
build_kstripes_layer <- function(ps, primary_scale, max_k) {
  div_roll <- ps$div_roll %||% ps$rolling
  primary_mat <- div_roll[[primary_scale]] %||%
                 div_roll[[paste0("s", primary_scale)]]
  if (is.null(primary_mat)) return(NULL)
  # Per-sample mean rank across windows (chromosome-wide). Samples are
  # then K-binned by quantile cut on this rank — this is the atlas's
  # cross-page sample coloring source.
  rank_mat <- apply(primary_mat, 2,
                    function(v) rank(v, na.last = "keep", ties.method = "average"))
  sample_mean_rank <- rowMeans(rank_mat, na.rm = TRUE)
  # rank-mean-normalize to 0..1 so K=2 / K=3 / K=4 quantile cuts are stable
  ranks_norm <- (rank(sample_mean_rank,
                      na.last = "keep", ties.method = "average") - 1) /
                (n_samples - 1)
  by_k <- list()
  for (K in 2:max_k) {
    cuts <- seq.int(0, 1, length.out = K + 1L)
    stripe <- as.integer(cut(ranks_norm, breaks = cuts,
                             include.lowest = TRUE, labels = FALSE))
    n_per <- as.integer(table(factor(stripe, levels = seq_len(K))))
    by_k[[as.character(K)]] <- list(
      stripe_per_sample = stripe,
      n_per_stripe      = n_per
    )
  }
  list(
    schema_version    = 2L,
    primary_scale     = primary_scale,
    n_samples         = as.integer(n_samples),
    sample_mean_rank  = clean_numeric(sample_mean_rank, 4),
    by_k              = by_k
  )
}
kstripes_layer <- build_kstripes_layer(ps, PRIMARY_SCALE, MAX_K)
if (!is.null(kstripes_layer))
  cat("[GHSL_EXPORT v3]   ghsl_kstripes: K=2..", MAX_K,
      " stripe assignments\n", sep = "")

# =============================================================================
# Layer: ghsl_karyotype_runs  (existing v2 contract)
# =============================================================================
build_karyo_runs_layer <- function(karyo, n_samples) {
  if (is.null(karyo)) return(NULL)
  # Karyotypes RDS may carry a data.table of stable LOW/HIGH runs per sample;
  # exact field name varies. Probe for a few common ones.
  runs <- karyo$karyotype_runs %||% karyo$runs %||% karyo$dt
  if (is.null(runs) || (is.data.frame(runs) && nrow(runs) == 0L)) return(NULL)
  runs_dt <- as.data.table(runs)
  records <- lapply(seq_len(nrow(runs_dt)), function(i) {
    r <- as.list(runs_dt[i])
    out <- list(
      sample_name = as.character(r$sample_name %||% r$sample %||% NA),
      start_bp    = as.integer(r$start_bp),
      end_bp      = as.integer(r$end_bp),
      class       = as.character(r$class %||% r$status %||% "UNKNOWN")
    )
    if (!is.null(r$n_windows)) out$n_windows <- as.integer(r$n_windows)
    if (!is.null(r$mean_score)) out$mean_score <- round4(as.numeric(r$mean_score))
    out
  })
  list(
    schema_version = 2L,
    n_runs         = nrow(runs_dt),
    runs           = records
  )
}
karyo_runs_layer <- build_karyo_runs_layer(karyo, n_samples)
if (!is.null(karyo_runs_layer))
  cat("[GHSL_EXPORT v3]   ghsl_karyotype_runs: ", karyo_runs_layer$n_runs,
      " stable LOW/HIGH runs\n", sep = "")

# =============================================================================
# Layer: ghsl_envelopes  (PRIMARY — STEP_C04b PASS-runs, NEW turn 3)
# =============================================================================
# Lift contiguous PASS-status windows from annot_dt into L2 candidate runs.
# Compute n_pass / n_weak / n_fail / mean_score / peak_score per run, plus
# bp coordinates from annot_dt.
build_primary_envelopes_layer <- function(annot_dt, status_col, score_col) {
  if (is.na(status_col)) return(NULL)
  status_vec <- annot_dt[[status_col]]
  is_pass <- status_vec == "PASS" & !is.na(status_vec)
  if (!any(is_pass)) {
    return(list(
      schema_version = 2L,
      layer          = "ghsl_envelopes",
      chrom          = CHROM,
      source         = paste0("STEP_C04b annot RDS (", status_col, " == 'PASS')"),
      l2_envelopes   = list(),
      l1_envelopes   = list()
    ))
  }
  # Run-length-encode contiguous PASS regions, allowing single-window WEAK/FAIL
  # gaps without breaking the run (matching the dosage pipeline's L2 merge logic)
  MERGE_GAP <- 3L
  runs <- rle(is_pass)
  # Identify run starts/ends in window indices
  ends   <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1) + 1L)
  pass_idx <- which(runs$values)
  if (length(pass_idx) == 0L) {
    return(list(
      schema_version = 2L, layer = "ghsl_envelopes", chrom = CHROM,
      source = "STEP_C04b annot RDS (no PASS windows found)",
      l2_envelopes = list(), l1_envelopes = list()
    ))
  }
  l2_dt <- data.table(
    win_start = starts[pass_idx],
    win_end   = ends[pass_idx],
    n_windows = runs$lengths[pass_idx]
  )
  # Merge L2 runs separated by ≤ MERGE_GAP windows of non-PASS
  setorder(l2_dt, win_start)
  if (nrow(l2_dt) > 1) {
    merged <- l2_dt[1]
    for (k in seq.int(2L, nrow(l2_dt))) {
      if (l2_dt$win_start[k] - merged$win_end[nrow(merged)] <= MERGE_GAP) {
        merged[nrow(merged), win_end := l2_dt$win_end[k]]
        merged[nrow(merged), n_windows := win_end - win_start + 1L]
      } else {
        merged <- rbind(merged, l2_dt[k])
      }
    }
    l2_dt <- merged
  }
  # Annotate each L2 with PASS/WEAK/FAIL counts + score stats
  score_vec <- if (!is.na(score_col)) annot_dt[[score_col]] else rep(NA_real_, n_windows)
  l2_dt[, `:=`(
    start_bp   = annot_dt$start_bp[win_start],
    end_bp     = annot_dt$end_bp[win_end],
    span_kb    = round((annot_dt$end_bp[win_end] -
                        annot_dt$start_bp[win_start]) / 1000, 1),
    n_pass     = vapply(seq_len(.N), function(k) {
      sub <- status_vec[win_start[k]:win_end[k]]
      sum(sub == "PASS", na.rm = TRUE)
    }, integer(1)),
    n_weak     = vapply(seq_len(.N), function(k) {
      sub <- status_vec[win_start[k]:win_end[k]]
      sum(sub == "WEAK", na.rm = TRUE)
    }, integer(1)),
    n_fail     = vapply(seq_len(.N), function(k) {
      sub <- status_vec[win_start[k]:win_end[k]]
      sum(sub == "FAIL", na.rm = TRUE)
    }, integer(1)),
    mean_score = vapply(seq_len(.N), function(k) {
      sub <- score_vec[win_start[k]:win_end[k]]
      mean(sub, na.rm = TRUE)
    }, numeric(1)),
    peak_score = vapply(seq_len(.N), function(k) {
      sub <- score_vec[win_start[k]:win_end[k]]
      max(sub, na.rm = TRUE)
    }, numeric(1)),
    l2_id      = paste0(CHROM, "_ghsl_L2_", sprintf("%03d", seq_len(.N)))
  )]
  # L1 = merge L2s within 9 windows
  L1_GAP <- 9L
  l1_dt <- if (nrow(l2_dt) == 0) data.table()
           else if (nrow(l2_dt) == 1) l2_dt[, .(win_start, win_end, n_l2 = 1L)]
           else {
             merged <- l2_dt[1, .(win_start, win_end, n_l2 = 1L)]
             for (k in seq.int(2L, nrow(l2_dt))) {
               if (l2_dt$win_start[k] - merged$win_end[nrow(merged)] <= L1_GAP) {
                 merged[nrow(merged), win_end := l2_dt$win_end[k]]
                 merged[nrow(merged), n_l2 := n_l2 + 1L]
               } else {
                 merged <- rbind(merged, l2_dt[k, .(win_start, win_end, n_l2 = 1L)])
               }
             }
             merged
           }
  if (nrow(l1_dt) > 0) {
    l1_dt[, `:=`(
      start_bp = annot_dt$start_bp[win_start],
      end_bp   = annot_dt$end_bp[win_end],
      span_kb  = round((annot_dt$end_bp[win_end] -
                        annot_dt$start_bp[win_start]) / 1000, 1),
      l1_id    = paste0(CHROM, "_ghsl_L1_", sprintf("%03d", seq_len(.N)))
    )]
  }
  list(
    schema_version = 2L,
    layer          = "ghsl_envelopes",
    chrom          = CHROM,
    source         = paste0("STEP_C04b annot RDS (", status_col, " == 'PASS')"),
    z_threshold_equivalent = NA,
    l2_envelopes   = lapply(seq_len(nrow(l2_dt)), function(k) {
      list(l2_id = l2_dt$l2_id[k],
           win_start = as.integer(l2_dt$win_start[k]),
           win_end   = as.integer(l2_dt$win_end[k]),
           start_bp  = as.integer(l2_dt$start_bp[k]),
           end_bp    = as.integer(l2_dt$end_bp[k]),
           span_kb   = l2_dt$span_kb[k],
           n_windows = as.integer(l2_dt$n_windows[k]),
           n_pass    = as.integer(l2_dt$n_pass[k]),
           n_weak    = as.integer(l2_dt$n_weak[k]),
           n_fail    = as.integer(l2_dt$n_fail[k]),
           mean_score = round4(l2_dt$mean_score[k]),
           peak_score = round4(l2_dt$peak_score[k]))
    }),
    l1_envelopes   = lapply(seq_len(nrow(l1_dt)), function(k) {
      list(l1_id = l1_dt$l1_id[k],
           win_start = as.integer(l1_dt$win_start[k]),
           win_end   = as.integer(l1_dt$win_end[k]),
           start_bp  = as.integer(l1_dt$start_bp[k]),
           end_bp    = as.integer(l1_dt$end_bp[k]),
           span_kb   = l1_dt$span_kb[k],
           n_l2      = as.integer(l1_dt$n_l2[k]))
    })
  )
}
envelopes_layer <- build_primary_envelopes_layer(annot_dt, status_col, score_col)
if (!is.null(envelopes_layer))
  cat("[GHSL_EXPORT v3]   ghsl_envelopes (PRIMARY): ",
      length(envelopes_layer$l2_envelopes), " L2 + ",
      length(envelopes_layer$l1_envelopes), " L1\n", sep = "")

# =============================================================================
# Layer: ghsl_local_pca + ghsl_secondary_envelopes  (turn 2/3, optional)
# =============================================================================
local_pca_layer <- NULL
sec_env_layer   <- NULL
if (!is.null(LOCALPCA_RDS) && file.exists(LOCALPCA_RDS)) {
  cat("[GHSL_EXPORT v3] Loading STEP_C04c localpca RDS\n")
  lp <- readRDS(LOCALPCA_RDS)

  # Pack sim_mat as upper triangle to halve serialization cost
  pack_upper_triangle <- function(sim) {
    n <- nrow(sim)
    out <- numeric(n * (n + 1) / 2)
    pos <- 1L
    for (i in seq_len(n)) {
      ll <- n - i + 1L
      out[pos:(pos + ll - 1L)] <- sim[i, i:n]
      pos <- pos + ll
    }
    out
  }
  sim_packed <- pack_upper_triangle(lp$sim_mat)

  # Per-window loadings as list-of-arrays (n_windows entries, each length n_samples)
  pc1_aligned <- lp$pc1_loadings_aligned_mat
  pc2_aligned <- lp$pc2_loadings_aligned_mat
  pc1_raw     <- lp$pc1_loadings_mat
  pc2_raw     <- lp$pc2_loadings_mat

  loadings_to_lol <- function(mat) {
    if (is.null(mat)) return(NULL)
    lapply(seq_len(ncol(mat)), function(wi) clean_numeric(mat[, wi], 6))
  }

  local_pca_layer <- list(
    schema_version = 2L,
    layer          = "ghsl_local_pca",
    chrom          = CHROM,
    scale          = "5kb",
    pad            = as.integer(lp$params$pad %||% 1L),
    n_samples      = as.integer(lp$n_samples),
    n_windows      = as.integer(lp$n_windows),
    sample_order   = lp$sample_names,
    pc1_loadings           = loadings_to_lol(pc1_raw),
    pc2_loadings           = loadings_to_lol(pc2_raw),
    pc1_loadings_aligned   = loadings_to_lol(pc1_aligned),
    pc2_loadings_aligned   = loadings_to_lol(pc2_aligned),
    lambda_1       = clean_numeric(lp$lambda_1, 6),
    lambda_2       = clean_numeric(lp$lambda_2, 6),
    lambda_ratio   = clean_numeric(lp$lambda_ratio, 4),
    z_profile      = clean_numeric(lp$z_profile, 4),
    z_top10_mean   = clean_numeric(lp$z_top10_mean, 4),
    mds_coords     = list(
      mds1 = clean_numeric(lp$mds1, 6),
      mds2 = clean_numeric(lp$mds2, 6)
    ),
    anchor_window_idx = as.integer(lp$anchor_window_idx),
    sim_mat_format    = "upper_triangle_float32",
    sim_mat_n         = as.integer(nrow(lp$sim_mat)),
    sim_mat           = clean_numeric(sim_packed, 4),
    `_compute_meta`   = list(
      weighting       = "n_phased_het",
      smoothing_input = lp$params$smoothing_scale %||% "none"
    )
  )
  cat("[GHSL_EXPORT v3]   ghsl_local_pca: ", lp$n_samples, " samples × ",
      lp$n_windows, " windows; sim_mat upper-triangle packed (",
      length(sim_packed), " floats)\n", sep = "")

  # Secondary envelopes from local-PCA z_profile threshold
  if (!is.null(lp$secondary_l2_envelopes) && nrow(lp$secondary_l2_envelopes) > 0) {
    sec_l2_records <- lapply(seq_len(nrow(lp$secondary_l2_envelopes)), function(k) {
      r <- lp$secondary_l2_envelopes[k]
      list(sec_id = r$sec_id, win_start = as.integer(r$win_start),
           win_end = as.integer(r$win_end),
           start_bp = as.integer(r$start_bp), end_bp = as.integer(r$end_bp),
           span_kb = r$span_kb, n_windows = as.integer(r$n_windows),
           peak_z = round4(r$peak_z), mean_z = round4(r$mean_z))
    })
  } else sec_l2_records <- list()
  if (!is.null(lp$secondary_l1_envelopes) && nrow(lp$secondary_l1_envelopes) > 0) {
    sec_l1_records <- lapply(seq_len(nrow(lp$secondary_l1_envelopes)), function(k) {
      r <- lp$secondary_l1_envelopes[k]
      list(sec_id = r$sec_id, win_start = as.integer(r$win_start),
           win_end = as.integer(r$win_end),
           start_bp = as.integer(r$start_bp), end_bp = as.integer(r$end_bp),
           span_kb = r$span_kb, n_l2 = as.integer(r$n_l2))
    })
  } else sec_l1_records <- list()
  sec_env_layer <- list(
    schema_version = 2L,
    layer          = "ghsl_secondary_envelopes",
    chrom          = CHROM,
    source         = "STEP_C04c local-PCA z_profile threshold",
    z_threshold    = as.numeric(lp$params$z_threshold %||% 2.5),
    min_l2_windows = as.integer(lp$params$min_l2_windows %||% 5L),
    merge_gap      = as.integer(lp$params$merge_gap %||% 3L),
    secondary_l2_envelopes = sec_l2_records,
    secondary_l1_envelopes = sec_l1_records
  )
  cat("[GHSL_EXPORT v3]   ghsl_secondary_envelopes: ",
      length(sec_l2_records), " L2 + ",
      length(sec_l1_records), " L1\n", sep = "")
}

# =============================================================================
# Layer: ghsl_d17_envelopes  (turn 5, optional — from D17 wrapper TSVs)
# =============================================================================
d17_layer <- NULL
load_d17_tsv <- function(path, kind) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  cat("[GHSL_EXPORT v3] Loading D17 ", kind, " TSV: ",
      basename(path), "\n", sep = "")
  fread(path)
}
d17_l1_env <- load_d17_tsv(D17_L1_ENV, "L1 envelopes")
d17_l2_env <- load_d17_tsv(D17_L2_ENV, "L2 envelopes")
d17_l1_bnd <- load_d17_tsv(D17_L1_BND, "L1 boundaries")
d17_l2_bnd <- load_d17_tsv(D17_L2_BND, "L2 boundaries")
if (!is.null(d17_l1_env) || !is.null(d17_l2_env) ||
    !is.null(d17_l1_bnd) || !is.null(d17_l2_bnd)) {
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
        mean_sim     = round4(as.numeric(r$mean_sim %||% NA)),
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
      if ("grow_max_z" %in% names(dt))
        out$grow_max_z <- round4(as.numeric(r$grow_max_z))
      if ("parent_l1_id" %in% names(dt))
        out$parent_l1_id <- as.character(r$parent_l1_id)
      out
    })
  }
  d17_layer <- list(
    schema_version = 2L,
    layer          = "ghsl_d17_envelopes",
    chrom          = CHROM,
    source         = "STEP_C04d D17 cross-block boundary detector",
    detector       = "STEP_D17_multipass_L1_only_v7.R + STEP_D17_multipass_L2_v8.R",
    l1_envelopes   = d17_env_to_records(d17_l1_env),
    l2_envelopes   = d17_env_to_records(d17_l2_env),
    l1_boundaries  = d17_bnd_to_records(d17_l1_bnd),
    l2_boundaries  = d17_bnd_to_records(d17_l2_bnd)
  )
  cat("[GHSL_EXPORT v3]   ghsl_d17_envelopes: L1=",
      length(d17_layer$l1_envelopes), " L2=",
      length(d17_layer$l2_envelopes), " bL1=",
      length(d17_layer$l1_boundaries), " bL2=",
      length(d17_layer$l2_boundaries), "\n", sep = "")
}

# =============================================================================
# Top-level JSON
# =============================================================================
out <- list(
  schema_version  = 2L,
  chrom           = CHROM,
  n_windows       = as.integer(n_windows),
  n_samples       = as.integer(n_samples),
  `_generated_at` = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  `_generator`    = "export_ghsl_to_json_v3.R",
  samples = lapply(seq_len(nrow(sample_meta)), function(i) {
    list(ind = sample_meta$ind[i], cga = sample_meta$cga[i],
         ancestry = sample_meta$ancestry[i])
  })
)

layers_present <- character(0)
if (length(tracks_layer) > 0) {
  out$tracks <- tracks_layer
  layers_present <- c(layers_present, "tracks")
}
if (!is.null(panel_layer)) {
  out$ghsl_panel <- panel_layer
  layers_present <- c(layers_present, "ghsl_panel")
}
if (!is.null(kstripes_layer)) {
  out$ghsl_kstripes <- kstripes_layer
  layers_present <- c(layers_present, "ghsl_kstripes")
}
if (!is.null(karyo_runs_layer)) {
  out$ghsl_karyotype_runs <- karyo_runs_layer
  layers_present <- c(layers_present, "ghsl_karyotype_runs")
}
if (!is.null(envelopes_layer)) {
  out$ghsl_envelopes <- envelopes_layer
  layers_present <- c(layers_present, "ghsl_envelopes")
}
if (!is.null(local_pca_layer)) {
  out$ghsl_local_pca <- local_pca_layer
  layers_present <- c(layers_present, "ghsl_local_pca")
}
if (!is.null(sec_env_layer)) {
  out$ghsl_secondary_envelopes <- sec_env_layer
  layers_present <- c(layers_present, "ghsl_secondary_envelopes")
}
if (!is.null(d17_layer)) {
  out$ghsl_d17_envelopes <- d17_layer
  layers_present <- c(layers_present, "ghsl_d17_envelopes")
}
out$`_layers_present` <- layers_present

# ---- Write ------------------------------------------------------------------
cat("\n[GHSL_EXPORT v3] Layers present: ", paste(layers_present, collapse = ", "),
    "\n", sep = "")
cat("[GHSL_EXPORT v3] Serializing JSON...\n", sep = "")
t_ser <- proc.time()
write_json(out, OUT_JSON,
           auto_unbox = TRUE, na = "null", pretty = FALSE, digits = NA)
cat("[GHSL_EXPORT v3]   ser time: ",
    round((proc.time() - t_ser)[3], 1), "s\n", sep = "")

fi <- file.info(OUT_JSON)
cat("[GHSL_EXPORT v3] DONE — ", OUT_JSON, " (",
    round(fi$size / 1024 / 1024, 2), " MB)\n", sep = "")
