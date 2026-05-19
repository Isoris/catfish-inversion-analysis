#!/usr/bin/env Rscript

# =============================================================================
# 08b_export_atlas_json_localpca_zblocks.R
#
# Per-chromosome atlas JSON exporter. Consumes the full precomp.rds, the
# L1+L2 envelope/boundary TSVs from steps 04+06, the multi-scale NN
# similarity matrices from step 03, and the merged sample_metadata.tsv
# from step 08a. Emits ONE <chr>.atlas.json drag-droppable into the
# pca_scrubber / atlas page-1 viewer.
#
# Pipeline position:
#   03 precompute      -| 
#   04 detect_L1       -|
#   05 plot_L1         -|->  08b ATLAS_JSON  ->  pca_scrubber_v3
#   06 detect_L2       -|
#   07 plot_L2         -|
#   08a sample_metadata-|
#
# ── JSON contents ─────────────────────────────────────────────────────────
#   * multi-scale sim_mat thumbnails        (nn40, nn80, nn160, nn320 by default)
#   * per-distance local Z scores           (matches the PDF colorimetry)
#   * L1 envelopes + L1 boundaries          (full validation status)
#   * L2 envelopes + L2 boundaries          (segment_stats + quadrant info)
#   * samples[i].cga / family_id / ancestry (from sample_metadata.tsv)
#   * default_sim_scale                     (which scale to render initially)
#   * family_source                         "sample_metadata" in canonical mode
#
# ── Inputs (canonical mode) ───────────────────────────────────────────────
#   --precomp_dir     <dir>    auto-resolves precomp + sim_mats nn{40,80,160,320}
#   --L1_dir          <dir>    auto-resolves <chr>.L1_envelopes/_boundaries.tsv
#   --L2_dir          <dir>    auto-resolves <chr>.L2_envelopes/_boundaries.tsv
#   --chr             <label>  chromosome to export
#   --sample_metadata <tsv>    pre-merged identity TSV from step 08a, columns:
#                                ind, cga, family_id, ancestry
#   --out             <path>   output JSON (e.g. <chr>.atlas.json)
#   [--nn_list 40,80,160,320]  comma-separated NN scales to embed
#
# ── Inputs (legacy / power-user mode) ─────────────────────────────────────
#   --precomp <path>           --sim_mat_nn40 <path> ... --sim_mat_nn320 <path>
#   --l1_envelopes <path>      --l1_boundaries <path>
#   --l2_envelopes <path>      --l2_boundaries <path>
#   --bamlist <path>           --pairs <path> [--theta_cutoff 0.177]
#   --samples <path>           --family <path>
#   --track_<label> <tsv>      add a per-window scalar track to the scrubber
#                              (TSV cols: chrom + start_bp/end_bp/value
#                              OR chrom + window_idx + value
#                              OR chrom + center_bp + value)
#
# ── Sample metadata ───────────────────────────────────────────────────────
# Canonical mode: --sample_metadata short-circuits the legacy reconciliation
# (--bamlist + --pairs/--family + --samples). The merged TSV from 08a is the
# single source of truth. Same TSV is reusable across paths 1/2/3.
# Legacy 3-flag mode still works if --sample_metadata is not given.
#
# ── Knob defaults ─────────────────────────────────────────────────────────
#   --thumb_n     200     thumbnail downsample target
#   --sim_q_lo    0.05    sim quantile clip low
#   --sim_q_hi    0.95    sim quantile clip high
#   --z_clip      5       Z range clip for colormap
#   --z_max_min   2.5     min |Z| to count a window as 'lit'
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: export_precomp_to_json_v3.R)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- CLI ---------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) == 0) return(default); args[i + 1]
}
PRECOMP    <- get_arg("--precomp")
SIMMAT     <- get_arg("--sim_mat", NULL)        # legacy: unlabeled default
SAMPLES    <- get_arg("--samples", NULL)
BAMLIST    <- get_arg("--bamlist", NULL)        # one CGA per line, precomp order
FAMILY     <- get_arg("--family", NULL)         # 2-col TSV: cga<TAB>family_id
PAIRS      <- get_arg("--pairs", NULL)          # 3-col: id1 id2 theta
THETA_CUT  <- as.numeric(get_arg("--theta_cutoff", "0.177"))
THETA      <- get_arg("--theta", NULL)
L1_ENV     <- get_arg("--l1_envelopes",  NULL)
L2_ENV     <- get_arg("--l2_envelopes",  NULL)
L2_BND     <- get_arg("--l2_boundaries", NULL)
L1_BND     <- get_arg("--l1_boundaries", NULL)
OUT        <- get_arg("--out", "precomp_scrubber.json")
THUMB_N    <- as.integer(get_arg("--thumb_n", "200"))
SIM_Q_LO   <- as.numeric(get_arg("--sim_q_lo", "0.05"))
SIM_Q_HI   <- as.numeric(get_arg("--sim_q_hi", "0.95"))
Z_CLIP     <- as.numeric(get_arg("--z_clip", "5"))
Z_MAX_MIN  <- as.numeric(get_arg("--z_max_min", "2.5"))

# ---- Canonical-mode flags ---------------------------------------------------
# New v3.1 — cleaner alternative to passing every file path explicitly.
#
# --sample_metadata <tsv>    : pre-merged identity TSV from 08a_build_sample_metadata.R
#                              (columns: ind, cga, family_id, ancestry). Replaces
#                              the trio of --bamlist + --pairs + --samples.
# --precomp_dir <dir>        : auto-resolves --precomp + --sim_mat_nn{40,80,160,320}
# --chr <chr_label>          : the chromosome to export (e.g. C_gar_LG28)
# --L1_dir / --L2_dir <dir>  : auto-resolve L1/L2 envelopes + boundaries TSVs
#
# Power users can still pass explicit per-file flags to override.
SAMPLE_META <- get_arg("--sample_metadata", NULL)
PRECOMP_DIR <- get_arg("--precomp_dir",     NULL)
CHR         <- get_arg("--chr",             NULL)
L1_DIR      <- get_arg("--L1_dir",          NULL)
L2_DIR      <- get_arg("--L2_dir",          NULL)
NN_LIST     <- strsplit(get_arg("--nn_list", "40,80,160,320"), ",")[[1]]

if (is.null(PRECOMP) && !is.null(PRECOMP_DIR) && !is.null(CHR)) {
  cand <- file.path(PRECOMP_DIR, paste0(CHR, ".precomp.rds"))
  if (file.exists(cand)) { PRECOMP <- cand; message("[export v3] resolved --precomp: ", cand) }
}
if (!is.null(PRECOMP_DIR) && !is.null(CHR)) {
  for (nn in NN_LIST) {
    label <- paste0("nn", trimws(nn))
    if (is.null(sim_scales_user[[label]])) {
      cand <- file.path(PRECOMP_DIR, "sim_mats",
                        paste0(CHR, ".sim_mat_", label, ".rds"))
      if (file.exists(cand)) {
        # injected lazily; the discovery loop below doesn't see CLI flags so
        # we stash them into sim_scales_user directly.
        # (sim_scales_user is initialized below — pre-create here.)
      }
    }
  }
}
if (is.null(L1_ENV) && !is.null(L1_DIR) && !is.null(CHR)) {
  cand <- file.path(L1_DIR, paste0(CHR, ".L1_envelopes.tsv"))
  if (file.exists(cand)) { L1_ENV <- cand; message("[export v3] resolved --l1_envelopes: ", cand) }
}
if (is.null(L1_BND) && !is.null(L1_DIR) && !is.null(CHR)) {
  cand <- file.path(L1_DIR, paste0(CHR, ".L1_boundaries.tsv"))
  if (file.exists(cand)) { L1_BND <- cand; message("[export v3] resolved --l1_boundaries: ", cand) }
}
if (is.null(L2_ENV) && !is.null(L2_DIR) && !is.null(CHR)) {
  cand <- file.path(L2_DIR, paste0(CHR, ".L2_envelopes.tsv"))
  if (file.exists(cand)) { L2_ENV <- cand; message("[export v3] resolved --l2_envelopes: ", cand) }
}
if (is.null(L2_BND) && !is.null(L2_DIR) && !is.null(CHR)) {
  cand <- file.path(L2_DIR, paste0(CHR, ".L2_boundaries.tsv"))
  if (file.exists(cand)) { L2_BND <- cand; message("[export v3] resolved --l2_boundaries: ", cand) }
}

stopifnot(!is.null(PRECOMP), file.exists(PRECOMP))

# Discover --sim_mat_<label> flags. Each provides one named scale.
# Example: --sim_mat_nn40 file.rds  -> scales[["nn40"]] <- file.rds
sim_scales_user <- list()
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (grepl("^--sim_mat_", a) && i < length(args)) {
    label <- sub("^--sim_mat_", "", a)
    sim_scales_user[[label]] <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

# Auto-resolve sim_mats from --precomp_dir + --chr + --nn_list (default: nn40,80,160,320)
if (!is.null(PRECOMP_DIR) && !is.null(CHR)) {
  for (nn in NN_LIST) {
    label <- paste0("nn", trimws(nn))
    if (is.null(sim_scales_user[[label]])) {
      cand <- file.path(PRECOMP_DIR, "sim_mats",
                        paste0(CHR, ".sim_mat_", label, ".rds"))
      if (file.exists(cand)) {
        sim_scales_user[[label]] <- cand
        message("[export v3] resolved --sim_mat_", label, ": ", cand)
      }
    }
  }
}

# ---- Sample metadata: --sample_metadata replaces --bamlist + --pairs + --samples ----
# If --sample_metadata is given, populate the legacy variables from the merged TSV.
if (!is.null(SAMPLE_META) && file.exists(SAMPLE_META)) {
  message("[export v3] Loading sample_metadata: ", SAMPLE_META)
  smeta <- fread(SAMPLE_META, header = TRUE)
  setnames(smeta, tolower(names(smeta)))
  required <- c("ind", "cga", "family_id", "ancestry")
  miss <- setdiff(required, names(smeta))
  if (length(miss) > 0)
    stop("[export v3] sample_metadata TSV missing required columns: ",
         paste(miss, collapse = ", "))
  # Stash for use later in the script. The legacy SAMPLES/BAMLIST/FAMILY paths
  # below will be skipped because we'll inject sample_meta directly.
  message("[export v3]   loaded ", nrow(smeta), " samples from sample_metadata TSV")
  attr(smeta, "from_sample_metadata") <- TRUE
  SAMPLE_META_DT <- smeta
} else {
  SAMPLE_META_DT <- NULL
}

# Discover --track_<label> flags. Each provides one per-window scalar track.
# Example: --track_theta theta.tsv -> tracks_user[["theta"]] <- "theta.tsv"
# Each TSV is expected to have either:
#   (a) columns chrom, start_bp, end_bp, value          (range form, more flexible)
#   (b) columns chrom, window_idx (1-based), value      (index form, exact match)
# Aligned to the precomp's window grid by nearest-center-bp matching.
tracks_user <- list()
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (grepl("^--track_", a) && i < length(args)) {
    label <- sub("^--track_", "", a)
    tracks_user[[label]] <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

message("[export v3] Loading precomp: ", PRECOMP)
pc    <- readRDS(PRECOMP)
dt    <- as.data.table(pc$dt)
chrom <- pc$chrom %||% dt$chrom[1]
n_win <- pc$n_windows %||% nrow(dt)

# ---- PC columns --------------------------------------------------------------
pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
pc2_cols <- grep("^PC_2_", names(dt), value = TRUE)
if (length(pc1_cols) == 0) stop("No PC_1_* columns in precomp$dt")
sample_ids <- sub("^PC_1_", "", pc1_cols)
n_samp     <- length(sample_ids)
have_pc2   <- length(pc2_cols) == length(pc1_cols)
message("[export v3] Windows: ", n_win, "   Samples: ", n_samp,
        "   PC2 per-sample: ", if (have_pc2) "yes" else "no (slim precomp; will jitter)")
# Print available dt columns once (helps Z-column debugging)
non_pc <- setdiff(names(dt), c(pc1_cols, pc2_cols))
message("[export v3] dt columns (non-PC): ", paste(non_pc, collapse = ", "))

# ---- Sample metadata ---------------------------------------------------------
sample_meta <- data.table(ind = sample_ids, cga = sample_ids, ancestry = "unknown")

# NEW (canonical path): if --sample_metadata was given, use it directly and
# skip the legacy --bamlist / --pairs / --samples reconciliation logic below.
# The merged TSV from 08a already has ind, cga, family_id, ancestry.
USE_MERGED_META <- !is.null(SAMPLE_META_DT)
if (USE_MERGED_META) {
  # Align rows to precomp's PC_1_Ind* order via the ind column
  sample_meta <- merge(sample_meta[, .(ind)], SAMPLE_META_DT,
                       by = "ind", all.x = TRUE, sort = FALSE)
  if (!"family_id" %in% names(sample_meta)) sample_meta[, family_id := -1L]
  if (!"ancestry"  %in% names(sample_meta)) sample_meta[, ancestry  := "unknown"]
  if (!"cga"       %in% names(sample_meta)) sample_meta[, cga       := ind]
  sample_meta[is.na(cga),       cga       := ind]
  sample_meta[is.na(ancestry),  ancestry  := "unknown"]
  sample_meta[is.na(family_id), family_id := -1L]
  message("[export v3]   sample_metadata applied: ",
          sum(sample_meta$family_id >= 0L), " in families, ",
          sum(sample_meta$family_id <  0L), " isolates")
}

if (!USE_MERGED_META && !is.null(SAMPLES) && file.exists(SAMPLES)) {
  message("[export v3] Loading sample metadata: ", SAMPLES)
  smeta <- fread(SAMPLES, header = TRUE)
  setnames(smeta, tolower(names(smeta)))
  if ("ind" %in% names(smeta)) {
    sample_meta <- merge(sample_meta[, .(ind)], smeta, by = "ind",
                         all.x = TRUE, sort = FALSE)
  } else if ("cga" %in% names(smeta)) {
    sample_meta[, cga := ind]
    sample_meta <- merge(sample_meta, smeta, by = "cga", all.x = TRUE, sort = FALSE)
  }
  if (!"ancestry" %in% names(sample_meta)) sample_meta[, ancestry := "unknown"]
  if (!"cga" %in% names(sample_meta))      sample_meta[, cga := ind]
  sample_meta[is.na(cga),      cga      := ind]
  sample_meta[is.na(ancestry), ancestry := "unknown"]
}

# ---- Bamlist remap (Ind0..IndN -> CGAxxx) -----------------------------------
# Ported from STEP_D16 v3.2. The precomp's PC_1_Ind* columns use sequential
# indices in bamfile order. The bamlist file is one CGA ID per line, in the
# SAME order. Row 1 of bamlist = Ind0, row 2 = Ind1, etc.
if (!USE_MERGED_META && !is.null(BAMLIST) && file.exists(BAMLIST)) {
  message("[export v3] Loading bamlist for Ind -> CGA remap: ", BAMLIST)
  bam_ids <- readLines(BAMLIST)
  bam_ids <- bam_ids[nzchar(trimws(bam_ids))]
  if (length(bam_ids) != n_samp) {
    warning(sprintf(
      "[export v3] bamlist has %d entries but precomp has %d samples — using min",
      length(bam_ids), n_samp))
  }
  ncommon <- min(length(bam_ids), n_samp)
  ind_to_cga <- setNames(bam_ids[seq_len(ncommon)],
                          paste0("Ind", seq_len(ncommon) - 1L))
  # Trim ".bam" or full paths from bamlist entries — keep just the leaf basename
  # without extension. STEP_D16 leaves them as-is (so e.g. CGA0001 stays CGA0001),
  # but if your bamlist is full bam paths, this lets the scrubber show clean labels.
  ind_to_cga <- vapply(ind_to_cga, function(x) {
    base <- basename(x)
    sub("\\.(bam|cram|sam)$", "", base, ignore.case = TRUE)
  }, character(1))
  sample_meta[, cga := ind_to_cga[ind]]
  sample_meta[is.na(cga), cga := ind]
  message("[export v3]   remapped ", sum(!is.na(ind_to_cga[sample_meta$ind])),
          " / ", n_samp, " sample IDs to CGA")
  # Preview first 5 mappings so user can sanity-check the order
  preview <- paste0(head(sample_meta$ind, 5), " -> ", head(sample_meta$cga, 5))
  message("[export v3]   preview: ", paste(preview, collapse = " | "))
}

# ---- Family info (--pairs preferred, fallback to --family) ------------------
# Ported from STEP_D16 v3.2.
# Two routes, priority order:
#   1. --pairs + --theta_cutoff: union-find on the relatedness graph
#   2. --family: natora-style 2-col TSV
# Both produce family_id per sample. Unmatched -> -1.
build_components_from_pairs <- function(pairs_path, theta_thr, all_ids) {
  message("[export v3]   reading pairs file: ", pairs_path)
  first_line <- readLines(pairs_path, n = 1L)
  parts <- strsplit(first_line, "[\t ]+")[[1]]
  has_header <- !suppressWarnings(is.finite(as.numeric(parts[length(parts)])))
  pairs_dt <- fread(pairs_path, header = has_header,
                    col.names = c("id1", "id2", "theta"))
  pairs_dt[, theta := as.numeric(theta)]
  message("[export v3]   total pairs: ", nrow(pairs_dt),
          "  (header_detected=", has_header, ")")
  qs <- quantile(pairs_dt$theta, c(0.5, 0.9, 0.95, 0.99, 1.0), na.rm = TRUE)
  message(sprintf("[export v3]   theta dist: median=%.4f p90=%.4f p95=%.4f p99=%.4f max=%.4f",
                  qs[1], qs[2], qs[3], qs[4], qs[5]))
  thr <- theta_thr
  pairs_dt <- pairs_dt[theta >= thr]
  message("[export v3]   pairs above theta_cutoff=", thr, ": ", nrow(pairs_dt))

  all_ids_in_pairs <- unique(c(pairs_dt$id1, pairs_dt$id2))
  universe <- unique(c(all_ids, all_ids_in_pairs))
  N <- length(universe)
  id2idx <- setNames(seq_len(N), universe)
  parent <- seq_len(N)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[parent[x]] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  union_xy <- function(a, b) {
    ra <- find_root(a); rb <- find_root(b)
    if (ra != rb) parent[ra] <<- rb
  }
  if (nrow(pairs_dt) > 0) {
    a_idx <- id2idx[pairs_dt$id1]
    b_idx <- id2idx[pairs_dt$id2]
    valid <- !is.na(a_idx) & !is.na(b_idx)
    a_idx <- a_idx[valid]; b_idx <- b_idx[valid]
    for (i in seq_along(a_idx)) union_xy(a_idx[i], b_idx[i])
  }
  roots <- vapply(seq_len(N), find_root, integer(1))
  comp_size <- table(roots)
  ordered_roots <- as.integer(names(sort(comp_size, decreasing = TRUE)))
  fam_id_map <- setNames(seq_along(ordered_roots), ordered_roots)
  family_id <- fam_id_map[as.character(roots)]
  out <- data.table(cga = universe, family_id = as.integer(family_id))
  out <- out[cga %in% all_ids]
  out
}

family_source <- "none"
if (USE_MERGED_META) {
  family_source <- "sample_metadata"
  message("[export v3]   family_id taken from sample_metadata TSV")
} else if (!is.null(PAIRS) && file.exists(PAIRS)) {
  message("[export v3] Loading family info from pairs file (theta>=", THETA_CUT, ")")
  fam_dt <- build_components_from_pairs(PAIRS, THETA_CUT, all_ids = sample_meta$cga)
  sample_meta <- merge(sample_meta, fam_dt, by = "cga", all.x = TRUE, sort = FALSE)
  sample_meta[is.na(family_id), family_id := -1L]
  family_source <- "pairs"
  fam_size <- sample_meta[family_id != -1L, .N, by = family_id][order(-N)]
  n_hub <- sum(fam_size$N >= 4L)
  message("[export v3]   total components: ", nrow(fam_size),
          "  hubs (>=4 samples): ", n_hub)
  if (n_hub > 0) {
    msg <- paste0("F", fam_size$family_id[1:min(5, n_hub)],
                  "(n=", fam_size$N[1:min(5, n_hub)], ")",
                  collapse = ", ")
    message("[export v3]   biggest: ", msg)
  }
} else if (!is.null(FAMILY) && file.exists(FAMILY)) {
  message("[export v3] Loading family info from --family file (natora-style)")
  fam_dt <- fread(FAMILY, header = FALSE, col.names = c("cga", "family_id"))
  sample_meta <- merge(sample_meta, fam_dt, by = "cga", all.x = TRUE, sort = FALSE)
  sample_meta[is.na(family_id), family_id := -1L]
  family_source <- "family"
  message("[export v3]   ", length(unique(sample_meta$family_id)),
          " distinct family IDs (incl. unmatched=-1)")
} else {
  sample_meta[, family_id := -1L]
}

# Re-align sample_meta to precomp Ind order (merges may have shuffled it)
sample_meta <- sample_meta[match(sample_ids, ind)]


# ---- Z column ----------------------------------------------------------------
z_col <- NULL
for (cand in c("robust_z", "z_robust", "z", "mds_z_robust", "mds_z1_robust",
               "z_pc1", "robust_z_pc1", "lam_1_z", "z_lam1")) {
  if (cand %in% names(dt)) { z_col <- cand; break }
}
if (is.null(z_col)) {
  # Fall back to lam_1 (first-eigenvalue magnitude per window). This is a
  # sensible per-window "structure strength" signal — same as what most
  # local-PCA pipelines plot. The scrubber's Z-track ends up showing
  # |lam_1| in this case, with the threshold lines (1.2 / 1.8 / 2.5) being
  # less meaningful but the relative shape preserved.
  if ("lam_1" %in% names(dt)) {
    z_col <- "lam_1"
    message("[export v3] No labeled Z column found; falling back to lam_1 (first-eigenvalue magnitude)")
  } else {
    message("[export v3] No Z column or lam_1 found, setting Z to 0 (Z-track will be flat)")
    dt[, .export_z := 0]; z_col <- ".export_z"
  }
} else {
  message("[export v3] Z column: ", z_col)
}

# ---- Theta alignment ---------------------------------------------------------
theta_by_sample <- NULL
theta_range     <- NULL
if (!is.null(THETA) && file.exists(THETA)) {
  message("[export v3] Loading theta track: ", THETA)
  theta <- fread(THETA)
  theta <- theta[chrom == chrom]
  samples_in_theta <- intersect(sample_meta$cga, unique(theta$sample))
  message("[export v3] Theta samples matched: ", length(samples_in_theta), " / ", n_samp)

  win_centers <- as.integer((dt$start_bp + dt$end_bp) / 2)
  theta_by_sample <- matrix(NA_real_, nrow = n_win, ncol = n_samp)
  colnames(theta_by_sample) <- sample_meta$cga

  for (s in samples_in_theta) {
    ts <- theta[sample == s]
    if (nrow(ts) < 2) next
    setorder(ts, win_center_bp)
    spacing <- median(diff(ts$win_center_bp), na.rm = TRUE)
    max_gap <- 2 * spacing
    idx     <- findInterval(win_centers, ts$win_center_bp, all.inside = TRUE)
    left    <- pmax(1L, idx)
    right   <- pmin(nrow(ts), idx + 1L)
    dL      <- abs(win_centers - ts$win_center_bp[left])
    dR      <- abs(win_centers - ts$win_center_bp[right])
    nearest <- ifelse(dL <= dR, left, right)
    dmin    <- pmin(dL, dR)
    vals    <- ts$theta_pi_persite[nearest]
    vals[dmin > max_gap] <- NA_real_
    col_idx <- match(s, sample_meta$cga)
    theta_by_sample[, col_idx] <- vals
  }
  theta_range <- range(theta_by_sample, na.rm = TRUE, finite = TRUE)
  message(sprintf("[export v3] Theta range: [%.3e, %.3e]", theta_range[1], theta_range[2]))
}

# ---- Per-window scalar tracks (--track_<label> ...) -------------------------
# Each track is a per-window numeric value. The TSV is expected to have ONE of:
#   (a) chrom, start_bp, end_bp, value               # range form (recommended)
#   (b) chrom, window_idx, value                     # index form (1-based)
#   (c) chrom, center_bp, value                      # center form
# We align by nearest center_bp to the precomp's window grid.
# Output: tracks_json = { label1 = {values: [...], min, max, mean, n_finite}, ... }
tracks_json <- list()
if (length(tracks_user) > 0) {
  win_centers <- as.integer((dt$start_bp + dt$end_bp) / 2)
  win_spacing <- median(diff(win_centers), na.rm = TRUE)
  for (label in names(tracks_user)) {
    track_path <- tracks_user[[label]]
    if (!file.exists(track_path)) {
      warning(sprintf("[export v3] track '%s': file not found: %s — skipping",
                      label, track_path))
      next
    }
    message(sprintf("[export v3] Loading track '%s' from %s", label, track_path))
    tdt <- tryCatch(fread(track_path), error = function(e) NULL)
    if (is.null(tdt) || nrow(tdt) == 0) {
      warning(sprintf("[export v3] track '%s': empty or unreadable", label))
      next
    }
    setnames(tdt, tolower(names(tdt)))
    # Filter to current chromosome if a chrom column exists
    if ("chrom" %in% names(tdt)) tdt <- tdt[chrom == chrom]
    if (nrow(tdt) == 0) {
      warning(sprintf("[export v3] track '%s': no rows match chrom=%s", label, chrom))
      next
    }
    # Determine form
    has_range  <- all(c("start_bp", "end_bp") %in% names(tdt)) && "value" %in% names(tdt)
    has_idx    <- "window_idx" %in% names(tdt) && "value" %in% names(tdt)
    has_center <- "center_bp" %in% names(tdt) && "value" %in% names(tdt)
    values_vec <- rep(NA_real_, n_win)
    if (has_idx) {
      # 1-based window indices: directly assign
      idx <- as.integer(tdt$window_idx)
      ok  <- idx >= 1 & idx <= n_win
      values_vec[idx[ok]] <- as.numeric(tdt$value[ok])
      message(sprintf("[export v3]   '%s': index form, %d / %d windows filled",
                      label, sum(ok & is.finite(values_vec[idx[ok]])), n_win))
    } else if (has_range) {
      # Range form: align by nearest center
      tcenter <- as.integer((tdt$start_bp + tdt$end_bp) / 2)
      setorder(tdt, start_bp)
      tcenter <- as.integer((tdt$start_bp + tdt$end_bp) / 2)
      idx <- findInterval(win_centers, tcenter, all.inside = TRUE)
      left  <- pmax(1L, idx)
      right <- pmin(nrow(tdt), idx + 1L)
      dL <- abs(win_centers - tcenter[left])
      dR <- abs(win_centers - tcenter[right])
      nearest <- ifelse(dL <= dR, left, right)
      dmin <- pmin(dL, dR)
      vals <- as.numeric(tdt$value[nearest])
      vals[dmin > 2 * win_spacing] <- NA_real_
      values_vec <- vals
      message(sprintf("[export v3]   '%s': range form, %d / %d windows filled",
                      label, sum(is.finite(values_vec)), n_win))
    } else if (has_center) {
      tcenter <- as.integer(tdt$center_bp)
      setorder(tdt, center_bp)
      tcenter <- as.integer(tdt$center_bp)
      idx <- findInterval(win_centers, tcenter, all.inside = TRUE)
      left  <- pmax(1L, idx)
      right <- pmin(nrow(tdt), idx + 1L)
      dL <- abs(win_centers - tcenter[left])
      dR <- abs(win_centers - tcenter[right])
      nearest <- ifelse(dL <= dR, left, right)
      dmin <- pmin(dL, dR)
      vals <- as.numeric(tdt$value[nearest])
      vals[dmin > 2 * win_spacing] <- NA_real_
      values_vec <- vals
      message(sprintf("[export v3]   '%s': center form, %d / %d windows filled",
                      label, sum(is.finite(values_vec)), n_win))
    } else {
      warning(sprintf(
        "[export v3] track '%s': TSV must have one of (start_bp+end_bp+value), (window_idx+value), or (center_bp+value). Got cols: %s",
        label, paste(names(tdt), collapse = ", ")))
      next
    }
    finite_vec <- values_vec[is.finite(values_vec)]
    tracks_json[[label]] <- list(
      values   = signif(values_vec, 5),
      min      = if (length(finite_vec) > 0) as.numeric(min(finite_vec)) else NA_real_,
      max      = if (length(finite_vec) > 0) as.numeric(max(finite_vec)) else NA_real_,
      mean     = if (length(finite_vec) > 0) as.numeric(mean(finite_vec)) else NA_real_,
      n_finite = length(finite_vec)
    )
  }
  if (length(tracks_json) > 0) {
    message(sprintf("[export v3] Tracks embedded: %s",
                    paste(names(tracks_json), collapse = ", ")))
  }
}

# ---- Build per-window arrays -------------------------------------------------
message("[export v3] Building per-window PC arrays ...")
pc1_mat <- as.matrix(dt[, ..pc1_cols])
storage.mode(pc1_mat) <- "double"
if (have_pc2) {
  pc2_mat <- as.matrix(dt[, ..pc2_cols])
  storage.mode(pc2_mat) <- "double"
} else {
  pc2_mat <- NULL
  # Per-sample stable y-jitter for the PCA scatter when no real PC2 exists.
  # Same idea as D16's strip plot: deterministic per-sample, identical
  # across windows, so each sample sits at the SAME y everywhere — easy to
  # follow visually. Scaled to roughly match a typical PC2 range.
  set.seed(7L)
  pc2_jitter <- runif(n_samp, min = -0.10, max = 0.10)
}
round4  <- function(x) round(x, 4)

windows <- vector("list", n_win)
for (i in seq_len(n_win)) {
  pc2_vec <- if (have_pc2) round4(pc2_mat[i, ]) else round4(pc2_jitter)
  w <- list(
    idx       = i - 1L,
    start_bp  = as.integer(dt$start_bp[i]),
    end_bp    = as.integer(dt$end_bp[i]),
    center_mb = round(((dt$start_bp[i] + dt$end_bp[i]) / 2) / 1e6, 4),
    z         = round4(dt[[z_col]][i]),
    lam1      = round4(dt$lam_1[i] %||% NA_real_),
    lam2      = round4(dt$lam_2[i] %||% NA_real_),
    pc1       = round4(pc1_mat[i, ]),
    pc2       = pc2_vec
  )
  if (!is.null(theta_by_sample)) w$theta <- signif(theta_by_sample[i, ], 4)
  windows[[i]] <- w
  if (i %% 500 == 0) message("  window ", i, " / ", n_win)
}

# ---- Sim_mat thumbnails (multi-scale, sim + local-Z, matches L2 overlay style) ---
# We support multiple smoothing scales. Each gets its own sim thumbnail and
# its own local-Z thumbnail, computed exactly like STEP_D17c_overlay_plot
# (compute_local_z + q05/q95 quantile anchors).

# Compute per-distance local Z (ported from STEP_D17c_overlay_plot_L2_v8.R)
compute_local_z <- function(sm, z_clip = 5) {
  Nl  <- nrow(sm)
  out <- matrix(NA_real_, Nl, Nl)
  for (d in 0:(Nl - 1L)) {
    if (d == 0L) {
      v <- diag(sm); ii <- seq_len(Nl); jj <- ii
    } else {
      ii <- seq.int(1L, Nl - d); jj <- ii + d
      v  <- sm[cbind(ii, jj)]
    }
    okv <- v[is.finite(v)]
    if (length(okv) < 5L) next
    mu <- mean(okv); sg <- sd(okv)
    if (!is.finite(sg) || sg < 1e-9) next
    z <- (v - mu) / sg
    out[cbind(ii, jj)] <- z
    out[cbind(jj, ii)] <- z
  }
  out[out >  z_clip] <-  z_clip
  out[out < -z_clip] <- -z_clip
  out
}

build_thumb <- function(sm, label) {
  if (is.null(sm) || !is.matrix(sm)) return(NULL)
  N    <- nrow(sm)
  step <- max(1, floor(N / THUMB_N))
  idx  <- seq(1, N, by = step)[seq_len(min(THUMB_N, ceiling(N / step)))]
  sub  <- sm[idx, idx, drop = FALSE]
  message(sprintf("[export v3]   %s: %dx%d full -> %dx%d thumbnail",
                  label, N, N, length(idx), length(idx)))
  z_sub <- compute_local_z(sub, z_clip = Z_CLIP)
  # Quantile anchors for the similarity colormap (q05/q95 by default)
  sv_finite <- sub[is.finite(sub)]
  q_lo <- if (length(sv_finite) > 0L)
            as.numeric(quantile(sv_finite, SIM_Q_LO, na.rm = TRUE)) else 0.0
  q_hi <- if (length(sv_finite) > 0L)
            as.numeric(quantile(sv_finite, SIM_Q_HI, na.rm = TRUE)) else 1.0
  if (!is.finite(q_lo)) q_lo <- 0.0
  if (!is.finite(q_hi) || q_hi <= q_lo) q_hi <- min(1.0, q_lo + 0.1)
  z_max_local <- max(abs(z_sub), na.rm = TRUE)
  if (!is.finite(z_max_local) || z_max_local < 1e-3) z_max_local <- 1
  z_max_local <- min(z_max_local, Z_CLIP)
  z_max_local <- max(z_max_local, Z_MAX_MIN)
  list(
    label    = label,
    n        = length(idx),
    sim      = round(as.vector(t(sub)), 3),
    z        = round(as.vector(t(z_sub)), 3),
    q_lo     = round(q_lo, 4),
    q_hi     = round(q_hi, 4),
    z_max    = round(z_max_local, 4)
  )
}

# Collect scale -> file path mapping. Priority order: explicit --sim_mat_*
# flags, then the legacy --sim_mat (label "default"), then pc$sim_mat
# (label "precomp").
scale_map <- list()
for (lab in names(sim_scales_user)) scale_map[[lab]] <- sim_scales_user[[lab]]
if (length(scale_map) == 0L && !is.null(SIMMAT) && file.exists(SIMMAT)) {
  scale_map[["default"]] <- SIMMAT
}

sim_scales_json <- list()
default_scale <- NULL

# pc$sim_mat goes in as "precomp" if present
if (!is.null(pc$sim_mat) && is.matrix(pc$sim_mat)) {
  message("[export v3] Building thumbnail for pc$sim_mat (label: precomp)")
  th <- build_thumb(pc$sim_mat, "precomp")
  if (!is.null(th)) sim_scales_json[["precomp"]] <- th
  default_scale <- "precomp"
}

# Each user-supplied scale
for (lab in names(scale_map)) {
  path <- scale_map[[lab]]
  if (!file.exists(path)) {
    message("[export v3]   skip ", lab, ": file not found: ", path)
    next
  }
  message("[export v3] Loading sim_mat scale '", lab, "' from ", path)
  sm_obj <- readRDS(path)
  # Support both bare matrix and a list with $sim_mat
  sm <- if (is.matrix(sm_obj)) sm_obj
        else if (!is.null(sm_obj$sim_mat) && is.matrix(sm_obj$sim_mat)) sm_obj$sim_mat
        else NULL
  if (is.null(sm)) {
    message("[export v3]   skip ", lab, ": cannot find a similarity matrix in RDS")
    next
  }
  if (nrow(sm) != n_win) {
    message("[export v3]   warn ", lab, ": nrow(sm)=", nrow(sm),
            " differs from precomp n_win=", n_win, "  -- including anyway")
  }
  th <- build_thumb(sm, lab)
  if (!is.null(th)) sim_scales_json[[lab]] <- th
  if (is.null(default_scale)) default_scale <- lab
}

# Pick a sensible default if the user has multiple scales: prefer nn40 ->
# nn80 -> nn160 -> nn20 -> nn320 -> first available
if (!is.null(default_scale)) {
  for (preferred in c("nn40", "nn80", "nn160", "nn20", "nn320")) {
    if (!is.null(sim_scales_json[[preferred]])) {
      default_scale <- preferred; break
    }
  }
}

# Backward-compat top-level keys: keep the old single-thumbnail keys
# pointing at the default scale, so older scrubber versions still render.
sim_thumb   <- if (!is.null(default_scale)) sim_scales_json[[default_scale]]$sim else NULL
sim_thumb_n <- if (!is.null(default_scale)) sim_scales_json[[default_scale]]$n   else 0L
message("[export v3] sim_mat scales embedded: ",
        if (length(sim_scales_json) == 0L) "(none)"
        else paste(names(sim_scales_json), collapse = ", "),
        if (!is.null(default_scale)) paste0("  (default: ", default_scale, ")") else "")


# ---- L1 envelopes ------------------------------------------------------------
# Schema (from STEP_D17_multipass_L1_only_v7.R):
#   chr, candidate_id, start_w, end_w, start_bp, end_bp, n_windows, scale_W,
#   mean_sim, density_p70, density_adaptive, threshold, trim_*, status
load_envelopes <- function(path, kind) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  message("[export v3] Loading ", kind, " envelopes: ", path)
  dt <- fread(path)
  if (nrow(dt) == 0) return(NULL)
  message("[export v3]   ", kind, " envelopes: ", nrow(dt))
  dt
}
l1_env_dt <- load_envelopes(L1_ENV, "L1")
l2_env_dt <- load_envelopes(L2_ENV, "L2")

env_to_list <- function(dt, want_parent = FALSE) {
  if (is.null(dt)) return(NULL)
  out <- vector("list", nrow(dt))
  for (i in seq_len(nrow(dt))) {
    rec <- list(
      candidate_id = as.character(dt$candidate_id[i]),
      start_w      = as.integer(dt$start_w[i]),
      end_w        = as.integer(dt$end_w[i]),
      start_bp     = as.integer(dt$start_bp[i]),
      end_bp       = as.integer(dt$end_bp[i]),
      n_windows    = as.integer(dt$n_windows[i]),
      mean_sim     = round4(as.numeric(dt$mean_sim[i])),
      density_p70  = round4(as.numeric(dt$density_p70[i])),
      status       = as.character(dt$status[i])
    )
    if (want_parent && "parent_l1_id" %in% names(dt)) {
      rec$parent_l1_id <- as.character(dt$parent_l1_id[i])
    }
    out[[i]] <- rec
  }
  out
}
l1_envelopes_json <- env_to_list(l1_env_dt, want_parent = FALSE)
l2_envelopes_json <- env_to_list(l2_env_dt, want_parent = TRUE)

# ---- L2 boundaries -----------------------------------------------------------
# Schema (from STEP_D17_multipass_L2_v8.R): we keep a tight set of columns —
# the rest are diagnostic and bloat the JSON.
load_boundaries <- function(path, kind) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  message("[export v3] Loading ", kind, " boundaries: ", path)
  dt <- fread(path)
  if (nrow(dt) == 0) return(NULL)
  message("[export v3]   ", kind, " boundaries: ", nrow(dt),
          "   STABLE_BLUE: ", sum(dt$validation_status == "STABLE_BLUE", na.rm = TRUE))
  dt
}
l2_bnd_dt <- load_boundaries(L2_BND, "L2")
l1_bnd_dt <- load_boundaries(L1_BND, "L1")

l2_boundaries_json <- NULL
if (!is.null(l2_bnd_dt)) {
  l2_boundaries_json <- vector("list", nrow(l2_bnd_dt))
  has_quad   <- "quad_verdict"      %in% names(l2_bnd_dt)
  has_qdemo  <- "quadrant_demoted"  %in% names(l2_bnd_dt)
  has_grow   <- "grow_max_z"        %in% names(l2_bnd_dt)
  has_dedup  <- "dedup_dropped"     %in% names(l2_bnd_dt)
  for (i in seq_len(nrow(l2_bnd_dt))) {
    rec <- list(
      boundary_idx      = as.character(l2_bnd_dt$boundary_idx[i]),
      parent_l1_id      = as.character(l2_bnd_dt$parent_l1_id[i]),
      boundary_w        = as.integer(l2_bnd_dt$boundary_w[i]),
      boundary_bp       = as.integer(l2_bnd_dt$boundary_bp[i]),
      validation_status = as.character(l2_bnd_dt$validation_status[i]),
      boundary_score    = round4(as.numeric(l2_bnd_dt$boundary_score[i]))
    )
    if (has_grow)  rec$grow_max_z       <- round4(as.numeric(l2_bnd_dt$grow_max_z[i]))
    if (has_quad)  rec$quad_verdict     <- as.character(l2_bnd_dt$quad_verdict[i])
    if (has_qdemo) rec$quadrant_demoted <- isTRUE(as.logical(l2_bnd_dt$quadrant_demoted[i]))
    if (has_dedup) rec$dedup_dropped    <- isTRUE(as.logical(l2_bnd_dt$dedup_dropped[i]))
    l2_boundaries_json[[i]] <- rec
  }
}

l1_boundaries_json <- NULL
if (!is.null(l1_bnd_dt)) {
  l1_boundaries_json <- vector("list", nrow(l1_bnd_dt))
  has_grow <- "grow_max_z" %in% names(l1_bnd_dt)
  for (i in seq_len(nrow(l1_bnd_dt))) {
    rec <- list(
      boundary_idx      = as.character(l1_bnd_dt$boundary_idx[i]),
      boundary_w        = as.integer(l1_bnd_dt$boundary_w[i]),
      boundary_bp       = as.integer(l1_bnd_dt$boundary_bp[i]),
      validation_status = as.character(l1_bnd_dt$validation_status[i]),
      boundary_score    = round4(as.numeric(l1_bnd_dt$boundary_score[i]))
    )
    if (has_grow) rec$grow_max_z <- round4(as.numeric(l1_bnd_dt$grow_max_z[i]))
    l1_boundaries_json[[i]] <- rec
  }
}

# ---- Assemble & write (v10 / schema v4) -------------------------------------
# Emit harmonized schema v4 via _shared/lib_atlas_json.R. Pipeline tag = "z".
# Optional per_window / carriers / cusum stay NULL (z path has no carriers).
.find_lib_atlas <- function() {
  cand <- Sys.getenv("LIB_ATLAS_JSON", unset = "")
  if (nzchar(cand) && file.exists(cand)) return(cand)
  cand <- file.path(Sys.getenv("SCRIPT_DIR_SHARED", unset = ""),
                   "lib_atlas_json.R")
  if (nzchar(Sys.getenv("SCRIPT_DIR_SHARED")) && file.exists(cand)) return(cand)
  cur <- tryCatch(normalizePath(getwd(), mustWork = FALSE),
                  error = function(e) ".")
  for (i in 1:8) {
    c2 <- file.path(cur, "_shared", "lib_atlas_json.R")
    if (file.exists(c2)) return(c2)
    cur <- dirname(cur)
  }
  stop("Could not find _shared/lib_atlas_json.R; set SCRIPT_DIR_SHARED or LIB_ATLAS_JSON")
}
source(.find_lib_atlas())

samples_json <- lapply(seq_len(nrow(sample_meta)), function(i) {
  list(ind = sample_meta$ind[i], cga = sample_meta$cga[i],
       ancestry = sample_meta$ancestry[i],
       family_id = as.integer(sample_meta$family_id[i] %||% -1L))
})

message("[export v4] Writing harmonized JSON: ", OUT)
atlas_json_v4(
  pipeline          = "z",
  chrom             = chrom,
  n_windows         = n_win,
  n_samples         = n_samp,
  samples           = samples_json,
  windows           = windows,
  sim_thumb         = sim_thumb,
  sim_thumb_n       = sim_thumb_n,
  sim_scales        = sim_scales_json,
  default_sim_scale = default_scale,
  sim_q_lo          = SIM_Q_LO,
  sim_q_hi          = SIM_Q_HI,
  z_clip            = Z_CLIP,
  z_max_min         = Z_MAX_MIN,
  z_column          = z_col,
  has_pc2           = have_pc2,
  family_source     = family_source,
  theta_cutoff      = if (family_source == "pairs") THETA_CUT else NULL,
  scale             = NULL,
  theta_range       = if (!is.null(theta_range)) as.numeric(theta_range) else NULL,
  has_theta         = !is.null(theta_by_sample),
  l1_envelopes      = l1_envelopes_json,
  l1_boundaries     = l1_boundaries_json,
  l2_envelopes      = l2_envelopes_json,
  l2_boundaries     = l2_boundaries_json,
  tracks            = if (length(tracks_json) > 0) tracks_json else NULL,
  generator         = "STEP_ZO_M_export_atlas_json.R",
  out_path          = OUT
)
fs <- file.info(OUT)$size
message(sprintf("[export v4] Done. Output size: %.1f MB", fs / 1e6))
