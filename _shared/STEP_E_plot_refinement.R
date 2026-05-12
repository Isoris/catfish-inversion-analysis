#!/usr/bin/env Rscript
# =============================================================================
# STEP_E_plot_refinement.R  (v1, 2026-05-12)
# =============================================================================
# Multi-page L1 + refinement overlay PDF for one chromosome.
#
#   Page 1            whole-chromosome heatmap with L1 envelopes outlined
#                     and L1 boundaries marked. GLOBAL per-distance Z.
#                     (Identical to 05_plot_L1_localpca_zblocks.R page 1.)
#   Pages 2..(N+1)    one per D17 boundary that was refined by STEP_E_refine.
#                     Three stacked panels:
#                       A: zoomed sim_mat (lower=sim, upper=Z), with glyphs
#                          on the upper-triangle diagonal at each refined
#                          node bp, and arcs between glyphs for edges.
#                       B: nodes summary table
#                       C: edges summary table
#
# Glyph encoding:
#   - position : (anchor_w, anchor_w) on upper-triangle diagonal,
#                anchor_w = window index closest to node bp_median
#   - size     : sqrt(n_samples), capped to a readable range
#   - fill     : mean feature value (theta_pi or GHSL div) at node bp,
#                averaged over contributing samples (continuous gradient).
#                If --band_tsv provided, fill instead = band-majority color.
#   - text     : n_samples integer, white on a dark halo
#
# Arc encoding (geom_curve in upper triangle):
#   - paired_flanks : thick deep-red solid line
#   - co_event      : medium gray dashed line
#   - weak_overlap  : thin pale dotted line
#
# Path-agnostic via --label {theta,ghsl}. Loaders match STEP_E_refine.
#
# Inputs:
#   --precomp_dir <dir>          for window grid + sim_mat (auto-resolves)
#   --L1_dir <dir>               for catalogue + boundaries from D17
#                                  (envelopes and boundaries TSVs, with the
#                                  same <label>_<chr>.L1_* prefix that
#                                  STEP_D17_boundary_detect_L1.R writes)
#   --refine_dir <dir>           directory holding STEP_E_refine outputs:
#                                  <label>_refinement_summary.tsv
#                                  <label>_<chr>_<boundary_idx>.nodes.tsv
#                                  <label>_<chr>_<boundary_idx>.edges.tsv
#                                  <label>_<chr>_<boundary_idx>.sample_peaks.tsv
#   --feature_source <dir>       same as in STEP_E_refine: theta TSV dir or
#                                  ghsl_matrices dir
#   --label theta|ghsl           selects loader + per-path defaults
#   --scale <str>                feature scale that refinement used
#                                  theta: "win10000.step2000"
#                                  ghsl:  "raw" / "s10" / "s50"
#   --chr <CHR>                  chromosome
#   --outdir <dir>               output dir
#
# Optional:
#   --precomp_suffix <str>       "precomp.rds" (theta, default) or "ghsl_precomp.rds"
#   --nn <int>                   sim_mat nn scale (default 80)
#   --band_tsv <path>            optional: sample_id -> band_id TSV. If
#                                provided, glyph fill = band majority color.
#   --window_kb <num>            refinement halfwidth (must match STEP_E_refine
#                                value; default 200 theta / 100 ghsl)
#   --page_size <num>            default 14
#
# Output:
#   <label>_<chr>.refinement_overlay.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggnewscale)
  library(scales)
  library(patchwork)
  library(grid)
  library(gridExtra)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# ── Args ────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

PRECOMP_DIR    <- get_arg("--precomp_dir")
L1_DIR         <- get_arg("--L1_dir")
REFINE_DIR     <- get_arg("--refine_dir")
FEATURE_SOURCE <- get_arg("--feature_source")
LABEL          <- get_arg("--label", "theta")
SCALE          <- get_arg("--scale",
                          if (get_arg("--label", "theta") == "ghsl") "raw" else "win10000.step2000")
CHR            <- get_arg("--chr")
OUTDIR         <- get_arg("--outdir", ".")
PRECOMP_SUFFIX <- get_arg("--precomp_suffix",
                          if (LABEL == "ghsl") "ghsl_precomp.rds" else "precomp.rds")
NN             <- as.integer(get_arg("--nn", "80"))
BAND_TSV       <- get_arg("--band_tsv")
WINDOW_KB      <- as.numeric(get_arg("--window_kb",
                                     if (LABEL == "ghsl") "100" else "200"))
PAGE_SIZE      <- as.numeric(get_arg("--page_size", "14"))
Z_CLIP         <- as.numeric(get_arg("--z_clip", "5"))

stopifnot(!is.na(PRECOMP_DIR), !is.na(L1_DIR), !is.na(REFINE_DIR),
          !is.na(FEATURE_SOURCE), !is.na(CHR))
stopifnot(LABEL %in% c("theta", "ghsl"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

LP <- paste0("[", LABEL, ".plot_refine] ")
cat(LP, "chr=", CHR, " label=", LABEL, " scale=", SCALE, "\n", sep = "")

# ── Resolve files ───────────────────────────────────────────────────────────
precomp_f <- file.path(PRECOMP_DIR, paste0(CHR, ".", PRECOMP_SUFFIX))
sim_mat_f <- file.path(PRECOMP_DIR, "sim_mats", sprintf("%s.sim_mat_nn%d.rds", CHR, NN))
cat_f     <- file.path(L1_DIR, sprintf("%s_%s.L1_envelopes.tsv",  LABEL, CHR))
bnd_f     <- file.path(L1_DIR, sprintf("%s_%s.L1_boundaries.tsv", LABEL, CHR))
summary_f <- file.path(REFINE_DIR, sprintf("%s_refinement_summary.tsv", LABEL))

for (f in c(precomp_f, sim_mat_f, cat_f, bnd_f, summary_f)) {
  if (!file.exists(f)) stop(LP, "missing input: ", f)
}
cat(LP, "precomp=", precomp_f, "\n", sep = "")
cat(LP, "sim_mat=", sim_mat_f, "\n", sep = "")
cat(LP, "L1 envelopes=", cat_f, "\n", sep = "")
cat(LP, "L1 boundaries=", bnd_f, "\n", sep = "")
cat(LP, "refinement summary=", summary_f, "\n", sep = "")

# ── Load precomp + sim_mat ──────────────────────────────────────────────────
pc <- readRDS(precomp_f)
dt_pc <- as.data.table(pc$dt)
n_win <- nrow(dt_pc)
win_starts <- if ("start_bp" %in% names(dt_pc)) dt_pc$start_bp else dt_pc$start
win_ends   <- if ("end_bp"   %in% names(dt_pc)) dt_pc$end_bp   else dt_pc$end
win_mids   <- as.integer((win_starts + win_ends) / 2L)
cat(LP, "N windows: ", n_win, "\n", sep = "")

sm_obj <- readRDS(sim_mat_f)
sim_mat <- if (is.matrix(sm_obj)) sm_obj
           else if (is.list(sm_obj) && !is.null(sm_obj$sim_mat)) sm_obj$sim_mat
           else if (is.list(sm_obj) && length(sm_obj) == 1L) sm_obj[[1]]
           else stop(LP, "sim_mat structure not recognized")
storage.mode(sim_mat) <- "double"
stopifnot(nrow(sim_mat) == n_win, ncol(sim_mat) == n_win)

# ── Load catalog + boundaries ───────────────────────────────────────────────
cat_dt <- fread(cat_f)
cat_dt <- cat_dt[status == "ENVELOPE" & !is.na(start_w) & !is.na(end_w)]
bdt    <- fread(bnd_f)
if ("validation_status" %in% names(bdt)) bdt <- bdt[validation_status == "STABLE_BLUE"]
cat(LP, "envelopes=", nrow(cat_dt), " boundaries=", nrow(bdt), "\n", sep = "")

# ── Load refinement summary, filter to this chrom ───────────────────────────
summary_dt <- fread(summary_f)
summary_dt <- summary_dt[chr == CHR]
cat(LP, "refined boundaries on this chrom: ", nrow(summary_dt), "\n", sep = "")

# ── Load feature data (same loaders as STEP_E_refine) ───────────────────────
load_feat_theta <- function(chr, feature_dir, scale_label) {
  tsv <- file.path(feature_dir, sprintf("theta_native.%s.%s.tsv.gz", chr, scale_label))
  if (!file.exists(tsv)) stop("[loader theta] missing ", tsv)
  long_dt <- fread(tsv)[chrom == chr]
  win_grid <- unique(long_dt[, .(window_idx, start_bp, end_bp)])
  setorder(win_grid, window_idx)
  win_grid[, mid_bp := as.integer((start_bp + end_bp) / 2L)]
  nw <- nrow(win_grid); sn <- sort(unique(long_dt$sample))
  samp_to_row <- setNames(seq_along(sn), sn)
  win_to_col  <- setNames(seq_len(nw), as.character(win_grid$window_idx))
  feat_mat <- matrix(NA_real_, nrow = length(sn), ncol = nw,
                     dimnames = list(sn, NULL))
  rows <- samp_to_row[long_dt$sample]; cols <- win_to_col[as.character(long_dt$window_idx)]
  good <- !is.na(rows) & !is.na(cols)
  feat_mat[cbind(rows[good], cols[good])] <- long_dt$theta_pi[good]
  list(feat_mat = feat_mat, win_mids = win_grid$mid_bp, sample_names = sn)
}
load_feat_ghsl <- function(chr, feature_dir, scale_label) {
  rds <- file.path(feature_dir, sprintf("%s.ghsl_matrices.rds", chr))
  gm <- readRDS(rds)
  feat_mat <- if (scale_label == "raw") gm$div_mat else gm$rolling[[scale_label]]
  if (is.null(feat_mat)) stop("[loader ghsl] scale '", scale_label, "' not found")
  list(feat_mat = feat_mat, win_mids = gm$window_info$mid_bp,
       sample_names = gm$sample_names)
}
fd <- if (LABEL == "theta") load_feat_theta(CHR, FEATURE_SOURCE, SCALE)
      else                  load_feat_ghsl( CHR, FEATURE_SOURCE, SCALE)
feat_mat <- fd$feat_mat
feat_mids <- fd$win_mids
feat_samples <- fd$sample_names

# Optional band assignment
band_dt <- NULL
band_palette <- c("1" = "#3B7BBF", "2" = "#E07B3F", "3" = "#7AB36F",
                  "4" = "#9F70B2", "5" = "#C46E8B", "6" = "#666666")
if (!is.na(BAND_TSV) && file.exists(BAND_TSV)) {
  band_dt <- fread(BAND_TSV)
  cat(LP, "band assignments loaded: ", nrow(band_dt), " samples\n", sep = "")
}

# ── Per-diagonal Z (global) ────────────────────────────────────────────────
compute_global_Z <- function(sm, n, z_clip) {
  z <- matrix(0, n, n)
  for (d in 0:(n - 1L)) {
    idx_i <- seq_len(n - d); idx_j <- idx_i + d
    idx_u <- cbind(idx_i, idx_j); idx_l <- cbind(idx_j, idx_i)
    v_u <- sm[idx_u]; v_l <- sm[idx_l]
    vals <- if (d == 0L) v_u else c(v_u, v_l)
    mu <- mean(vals, na.rm = TRUE); sg <- sd(vals, na.rm = TRUE)
    if (!is.finite(sg) || sg < 1e-9) sg <- 1
    z[idx_u] <- (v_u - mu) / sg
    if (d > 0L) z[idx_l] <- (v_l - mu) / sg
  }
  z[z >  z_clip] <-  z_clip; z[z < -z_clip] <- -z_clip
  z
}
compute_local_Z <- function(sm, z_clip) {
  nl <- nrow(sm); out <- matrix(NA_real_, nl, nl)
  for (d in 0:(nl - 1L)) {
    if (d == 0L) { v <- diag(sm); ii <- seq_len(nl); jj <- ii }
    else { ii <- seq.int(1L, nl - d); jj <- ii + d; v <- sm[cbind(ii, jj)] }
    okv <- v[is.finite(v)]
    if (length(okv) < 5L) next
    mu <- mean(okv); sg <- sd(okv); if (!is.finite(sg) || sg < 1e-9) next
    z <- (v - mu) / sg
    out[cbind(ii, jj)] <- z; out[cbind(jj, ii)] <- z
  }
  out[out >  z_clip] <-  z_clip; out[out < -z_clip] <- -z_clip
  out
}

cat(LP, "computing global Z\n", sep = "")
z_mat <- compute_global_Z(sim_mat, n_win, Z_CLIP)

# ── Sim layer & Z layer builder for a window range ──────────────────────────
build_layers <- function(w_lo, w_hi, use_local_z) {
  abs_idx <- w_lo:w_hi
  sm_sub <- sim_mat[abs_idx, abs_idx]
  z_sub  <- if (use_local_z) compute_local_Z(sm_sub, Z_CLIP) else z_mat[abs_idx, abs_idx]
  nl <- length(abs_idx)
  ii <- rep(abs_idx, times = nl); jj <- rep(abs_idx, each = nl)
  lo_mask <- jj <= ii; up_mask <- jj > ii
  list(
    sim = data.table(i = ii[lo_mask], j = jj[lo_mask],
                     v = as.numeric(sm_sub)[lo_mask]),
    z   = data.table(i = ii[up_mask], j = jj[up_mask],
                     v = as.numeric(z_sub)[up_mask])
  )
}

build_axis <- function(w_lo, w_hi) {
  br <- pretty(c(w_lo, w_hi), n = 8L)
  br <- unique(as.integer(round(br))); br <- br[br >= w_lo & br <= w_hi]
  lab <- sprintf("%.2f", win_starts[pmin(pmax(br, 1L), n_win)] / 1e6)
  list(breaks = br, labels = lab)
}

# Color limits for sim layer
all_sim <- as.numeric(sim_mat); all_sim <- all_sim[is.finite(all_sim)]
sim_q_lo <- as.numeric(quantile(all_sim, 0.05, na.rm = TRUE))
sim_q_hi <- as.numeric(quantile(all_sim, 0.98, na.rm = TRUE))

# ─────────────────────────────────────────────────────────────────────────────
# Page 1: whole-chromosome overview (same recipe as L1 plot)
# ─────────────────────────────────────────────────────────────────────────────
build_page1 <- function() {
  cat(LP, "building page 1 (whole chrom)\n", sep = "")
  layers <- build_layers(1L, n_win, use_local_z = FALSE)
  axes <- build_axis(1L, n_win)

  p <- ggplot() +
    geom_raster(data = layers$sim, aes(x = i, y = j, fill = v)) +
    scale_fill_gradientn(
      colours = c("#F8F8F8","#A8DBC2","#F2DC78","#E08838","#7E1F1F"),
      values  = scales::rescale(c(0, sim_q_lo, (sim_q_lo + sim_q_hi)/2, sim_q_hi, 1)),
      limits = c(0, 1), oob = scales::squish, name = "Similarity"
    ) +
    ggnewscale::new_scale_fill() +
    geom_raster(data = layers$z, aes(x = i, y = j, fill = v)) +
    scale_fill_gradient2(
      low = "#2C5AA0", mid = "#FAFAFA", high = "#B22222", midpoint = 0,
      limits = c(-Z_CLIP, Z_CLIP), oob = scales::squish,
      name = sprintf("Z (nn%d)", NN)
    )

  if (nrow(cat_dt) > 0L) {
    cdv <- copy(cat_dt)
    cdv[, draw_start := pmax(start_w, 1L)]
    cdv[, draw_end   := pmin(end_w,   n_win)]
    p <- p + geom_rect(
      data = cdv,
      aes(xmin = draw_start - 0.5, xmax = draw_end + 0.5,
          ymin = draw_start - 0.5, ymax = draw_end + 0.5),
      colour = "#1F3A6E", fill = NA, linewidth = 0.6, inherit.aes = FALSE
    )
  }

  if (nrow(bdt) > 0L) {
    bv <- copy(bdt)
    bW <- if ("boundary_W" %in% names(bv)) as.integer(bv$boundary_W[1]) else 5L
    bG <- if ("boundary_offset" %in% names(bv)) as.integer(bv$boundary_offset[1]) else 5L
    shift <- as.integer(ceiling((bG + 1L) / 2L))
    bv[, anchor_w := boundary_w + shift]
    bv[, sq_xmin := anchor_w - bW + 0.5]
    bv[, sq_xmax := anchor_w + 0.5]
    bv[, sq_ymin := anchor_w - 0.5]
    bv[, sq_ymax := anchor_w + bW - 0.5]
    bv <- bv[sq_xmax <= n_win + 0.5 & sq_ymax <= n_win + 0.5 &
             sq_xmin >= 0.5 & sq_ymin >= 0.5]
    if (nrow(bv) > 0L) {
      p <- p + ggnewscale::new_scale_fill() +
        geom_rect(data = bv,
                  aes(xmin = sq_xmin, xmax = sq_xmax,
                      ymin = sq_ymin, ymax = sq_ymax),
                  fill = "#C8102E", colour = NA, inherit.aes = FALSE) +
        geom_text(data = bv,
                  aes(x = anchor_w + (n_win * 0.005), y = anchor_w,
                      label = boundary_idx),
                  hjust = 0, size = 1.2, colour = "#7A0000", fontface = "bold",
                  inherit.aes = FALSE)
    }
  }

  p + coord_equal(xlim = c(0.5, n_win + 0.5), ylim = c(0.5, n_win + 0.5), expand = FALSE) +
    scale_x_continuous(name = "window", breaks = axes$breaks,
                       sec.axis = sec_axis(~ ., breaks = axes$breaks,
                                           labels = axes$labels, name = "Mb"),
                       expand = c(0, 0)) +
    scale_y_continuous(name = "window", breaks = axes$breaks,
                       sec.axis = sec_axis(~ ., breaks = axes$breaks,
                                           labels = axes$labels, name = "Mb"),
                       expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          plot.title = element_text(size = 12, face = "bold")) +
    labs(title = sprintf("%s | %s | whole chromosome | L1 envelopes (%d), boundaries (%d)",
                         CHR, LABEL, nrow(cat_dt), nrow(bdt)))
}

# ─────────────────────────────────────────────────────────────────────────────
# Refinement page (one per refined D17 boundary)
# ─────────────────────────────────────────────────────────────────────────────

bp_to_window <- function(bp) {
  # Returns the window index in the PRECOMP grid closest to bp
  idx <- which.min(abs(win_mids - bp))
  as.integer(idx)
}

# Find feature value at a given bp for a given sample list. Returns NA if
# the bp is outside the feature grid or no samples have values.
feature_mean_at_bp <- function(bp, sample_ids) {
  fidx <- which.min(abs(feat_mids - bp))
  if (length(fidx) == 0L) return(NA_real_)
  rows <- match(sample_ids, feat_samples)
  rows <- rows[!is.na(rows)]
  if (length(rows) == 0L) return(NA_real_)
  vals <- feat_mat[rows, fidx]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) return(NA_real_)
  mean(vals)
}

# Find band majority for a sample list (returns the most common band_id and
# its fraction). Returns list(band = NA, frac = NA) if no band_tsv.
band_majority <- function(sample_ids) {
  if (is.null(band_dt)) return(list(band = NA_character_, frac = NA_real_))
  bands <- band_dt[sample_id %in% sample_ids]$band_id
  bands <- bands[!is.na(bands)]
  if (length(bands) == 0L) return(list(band = NA_character_, frac = NA_real_))
  tab <- sort(table(bands), decreasing = TRUE)
  list(band = as.character(names(tab)[1]),
       frac = unname(as.numeric(tab[1])) / length(bands))
}

build_refinement_page <- function(row_summary) {
  bdy_idx <- row_summary$boundary_idx
  bdy_bp  <- as.integer(row_summary$bp_boundary)
  win_lo_bp <- bdy_bp - WINDOW_KB * 1000L
  win_hi_bp <- bdy_bp + WINDOW_KB * 1000L
  w_lo <- max(1L, bp_to_window(win_lo_bp))
  w_hi <- min(n_win, bp_to_window(win_hi_bp))
  if (w_hi - w_lo < 4L) {
    cat(LP, "  ", bdy_idx, ": refinement window too small, skip page\n", sep = "")
    return(NULL)
  }

  # Load node and edge tables (file names match STEP_E_refine output naming)
  nodes_f <- file.path(REFINE_DIR, sprintf("%s_%s_%s.nodes.tsv", LABEL, CHR, bdy_idx))
  edges_f <- file.path(REFINE_DIR, sprintf("%s_%s_%s.edges.tsv", LABEL, CHR, bdy_idx))
  if (!file.exists(nodes_f)) {
    cat(LP, "  ", bdy_idx, ": missing nodes file, skip\n", sep = ""); return(NULL)
  }
  nodes <- fread(nodes_f)
  edges <- if (file.exists(edges_f)) fread(edges_f) else data.table()
  if (nrow(nodes) == 0L) {
    cat(LP, "  ", bdy_idx, ": empty nodes table, skip\n", sep = ""); return(NULL)
  }

  # Map each node to a window index for plotting
  nodes[, anchor_w := vapply(bp_median, bp_to_window, integer(1))]

  # Compute glyph fill: mean feature, or band majority color
  nodes[, mean_feat := NA_real_]
  nodes[, band_id   := NA_character_]
  nodes[, band_frac := NA_real_]
  for (k in seq_len(nrow(nodes))) {
    samp_ids <- strsplit(nodes$samples_csv[k], ",")[[1]]
    nodes$mean_feat[k] <- feature_mean_at_bp(nodes$bp_median[k], samp_ids)
    bm <- band_majority(samp_ids)
    nodes$band_id[k]   <- bm$band
    nodes$band_frac[k] <- bm$frac
  }

  # Glyph fill color: prefer band majority if available, else mean_feat gradient
  use_bands <- !is.null(band_dt) && any(!is.na(nodes$band_id))

  # Build matrix layers (local Z for sharper per-window contrast)
  layers <- build_layers(w_lo, w_hi, use_local_z = TRUE)
  axes <- build_axis(w_lo, w_hi)
  span_w <- w_hi - w_lo + 1L

  # PANEL A: matrix + glyphs + arcs ──────────────────────────────────────────
  pA <- ggplot() +
    geom_raster(data = layers$sim, aes(x = i, y = j, fill = v)) +
    scale_fill_gradientn(
      colours = c("#F8F8F8","#A8DBC2","#F2DC78","#E08838","#7E1F1F"),
      values  = scales::rescale(c(0, sim_q_lo, (sim_q_lo + sim_q_hi)/2, sim_q_hi, 1)),
      limits = c(0, 1), oob = scales::squish, name = "Similarity"
    ) +
    ggnewscale::new_scale_fill() +
    geom_raster(data = layers$z, aes(x = i, y = j, fill = v)) +
    scale_fill_gradient2(
      low = "#2C5AA0", mid = "#FAFAFA", high = "#B22222", midpoint = 0,
      limits = c(-Z_CLIP, Z_CLIP), oob = scales::squish,
      name = "Z (local)"
    ) +
    # D17 boundary marker: vertical line on the diagonal at the original bp
    geom_point(
      data = data.table(x = bp_to_window(bdy_bp), y = bp_to_window(bdy_bp)),
      aes(x = x, y = y),
      shape = 4, size = 4, colour = "#1F1F1F", stroke = 1.0,
      inherit.aes = FALSE
    )

  # Add edge arcs first (so glyphs render on top)
  if (nrow(edges) > 0L) {
    # Need anchor_w for both endpoints
    node_map <- setNames(nodes$anchor_w, nodes$node_id)
    edges[, a_w := node_map[node_a]]
    edges[, b_w := node_map[node_b]]
    edges <- edges[!is.na(a_w) & !is.na(b_w) & a_w != b_w]
    if (nrow(edges) > 0L) {
      # Three relation classes — keep paired_flanks visible, fade weak_overlap
      ed_pf <- edges[inferred_relation == "paired_flanks"]
      ed_co <- edges[inferred_relation == "co_event"]
      ed_wk <- edges[inferred_relation == "weak_overlap"]
      mk_xy <- function(d) {
        # arc x,y: smaller window index, bigger; the arc curves through upper triangle
        d[, x  := pmin(a_w, b_w)]
        d[, y  := pmin(a_w, b_w)]
        d[, xe := pmax(a_w, b_w)]
        d[, ye := pmax(a_w, b_w)]
        d
      }
      if (nrow(ed_pf) > 0L) {
        ed_pf <- mk_xy(ed_pf)
        pA <- pA + geom_curve(
          data = ed_pf, aes(x = x, y = y, xend = xe, yend = ye),
          curvature = -0.30, colour = "#B22222", linewidth = 1.0,
          inherit.aes = FALSE
        )
      }
      if (nrow(ed_co) > 0L) {
        ed_co <- mk_xy(ed_co)
        pA <- pA + geom_curve(
          data = ed_co, aes(x = x, y = y, xend = xe, yend = ye),
          curvature = -0.30, colour = "#666666", linewidth = 0.5,
          linetype = "dashed", inherit.aes = FALSE
        )
      }
      if (nrow(ed_wk) > 0L) {
        ed_wk <- mk_xy(ed_wk)
        pA <- pA + geom_curve(
          data = ed_wk, aes(x = x, y = y, xend = xe, yend = ye),
          curvature = -0.30, colour = "#AAAAAA", linewidth = 0.3,
          linetype = "dotted", alpha = 0.5, inherit.aes = FALSE
        )
      }
    }
  }

  # Glyph sizing: sqrt(n_samples) scaled to a readable range
  n_max <- max(nodes$n_samples, na.rm = TRUE)
  nodes[, glyph_size := pmax(4, pmin(14, 2 + 1.5 * sqrt(n_samples)))]

  pA <- pA + ggnewscale::new_scale_fill()

  if (use_bands) {
    # Band-majority categorical fill
    pA <- pA + geom_point(
      data = nodes,
      aes(x = anchor_w, y = anchor_w, fill = band_id, size = glyph_size),
      shape = 21, colour = "#000000", stroke = 0.6, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = band_palette, name = "Band majority",
                      na.value = "#CCCCCC") +
    scale_size_identity()
  } else {
    # Continuous mean_feat fill (theta_pi or div). Gradient: blue (low) → white → red (high)
    feat_min <- min(nodes$mean_feat, na.rm = TRUE)
    feat_max <- max(nodes$mean_feat, na.rm = TRUE)
    if (!is.finite(feat_min) || !is.finite(feat_max) || feat_min == feat_max) {
      feat_mid <- feat_min; feat_min <- feat_min - 0.01; feat_max <- feat_max + 0.01
    } else feat_mid <- (feat_min + feat_max) / 2
    pA <- pA + geom_point(
      data = nodes,
      aes(x = anchor_w, y = anchor_w, fill = mean_feat, size = glyph_size),
      shape = 21, colour = "#000000", stroke = 0.6, inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2C5AA0", mid = "#FAFAFA", high = "#B22222",
      midpoint = feat_mid,
      limits = c(feat_min, feat_max), oob = scales::squish,
      name = if (LABEL == "ghsl") "mean GHSL div" else "mean θπ"
    ) +
    scale_size_identity()
  }

  # n_samples text inside each glyph
  pA <- pA + geom_text(
    data = nodes,
    aes(x = anchor_w, y = anchor_w, label = n_samples),
    colour = "#FFFFFF", fontface = "bold", size = 3.0,
    inherit.aes = FALSE
  )

  pA <- pA +
    coord_equal(xlim = c(w_lo - 0.5, w_hi + 0.5),
                ylim = c(w_lo - 0.5, w_hi + 0.5), expand = FALSE) +
    scale_x_continuous(name = "window", breaks = axes$breaks,
                       sec.axis = sec_axis(~ ., breaks = axes$breaks,
                                           labels = axes$labels, name = "Mb"),
                       expand = c(0, 0)) +
    scale_y_continuous(name = "window", breaks = axes$breaks,
                       sec.axis = sec_axis(~ ., breaks = axes$breaks,
                                           labels = axes$labels, name = "Mb"),
                       expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), legend.position = "right",
          legend.key.height = unit(0.5, "cm"),
          plot.title = element_text(size = 11, face = "bold")) +
    labs(title = sprintf("%s | %s | refinement: %s @ %.3f Mb (window ±%g kb)",
                         CHR, LABEL, bdy_idx, bdy_bp / 1e6, WINDOW_KB))

  # PANEL B: nodes table ─────────────────────────────────────────────────────
  nodes_disp <- nodes[, .(
    node_id,
    bp_Mb       = sprintf("%.4f", bp_median / 1e6),
    bp_CI       = sprintf("%.4f-%.4f", bp_ci_lo / 1e6, bp_ci_hi / 1e6),
    n_samples,
    direction   = direction_pattern,
    score_med   = round(score_med, 3),
    mean_feat   = round(mean_feat, 4),
    band_majority = if (use_bands)
      sprintf("%s (%.0f%%)", band_id, 100 * band_frac)
      else rep("(no bands)", .N)
  )]
  tt <- ttheme_minimal(
    core    = list(fg_params = list(fontsize = 8),
                   bg_params = list(fill = c("#FAFAFA", "#FFFFFF"))),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"))
  )
  pB <- gridExtra::tableGrob(nodes_disp, rows = NULL, theme = tt)

  # PANEL C: edges table ─────────────────────────────────────────────────────
  if (nrow(edges) > 0L) {
    edges_disp <- edges[, .(
      node_a, node_b,
      distance_kb = round(distance_bp / 1000, 1),
      n_shared = n_shared_samples,
      jaccard  = round(jaccard, 3),
      relation = inferred_relation
    )]
  } else {
    edges_disp <- data.table(note = "no edges (single node or all jaccards below threshold)")
  }
  pC <- gridExtra::tableGrob(edges_disp, rows = NULL, theme = tt)

  # Compose page with patchwork — heights 60/25/15
  pA_wrapped <- patchwork::wrap_elements(full = pA)
  pB_wrapped <- patchwork::wrap_elements(full = pB)
  pC_wrapped <- patchwork::wrap_elements(full = pC)

  page <- pA_wrapped / pB_wrapped / pC_wrapped +
    patchwork::plot_layout(heights = c(0.60, 0.25, 0.15))
  page
}

# ─────────────────────────────────────────────────────────────────────────────
# Render multi-page PDF
# ─────────────────────────────────────────────────────────────────────────────
out_pdf <- file.path(OUTDIR, sprintf("%s_%s.refinement_overlay.pdf", LABEL, CHR))
cat(LP, "writing ", out_pdf, "\n", sep = "")
pdf(out_pdf, width = PAGE_SIZE, height = PAGE_SIZE)

# Page 1
print(build_page1())

# Pages 2..N+1
if (nrow(summary_dt) > 0L) {
  setorder(summary_dt, bp_boundary)
  for (k in seq_len(nrow(summary_dt))) {
    row <- summary_dt[k]
    cat(LP, "page ", k + 1L, ": ", row$boundary_idx, " @ ",
        round(row$bp_boundary / 1e6, 3), " Mb\n", sep = "")
    page <- build_refinement_page(row)
    if (!is.null(page)) print(page)
  }
}

dev.off()
cat(LP, "done\n", sep = "")
