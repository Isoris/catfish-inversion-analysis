#!/usr/bin/env Rscript

# =============================================================================
# 07_plot_L2_localpca_zblocks.R
#
# Multi-page L1+L2 overlay PDF for one chromosome.
#
#   Page 1            whole chromosome (nn80 by default), L1 envelopes +
#                     L1 boundaries, same as step 05 page 1.
#   Pages 2..(N+1)    one zoomed page per L1 SEGMENT, drawn on nn40 by
#                     default (--sim_mat_l2). Shows the segment's L1 outline
#                     plus all L2 boundaries inside, as red squares.
#
# Pipeline position:
#   06 detect_L2  -> 07 PLOT_L2
#
# ── Default NN scales: 80 chrom-wide + 40 inside-segment ──────────────────
# Override via --nn N (chrom) and --nn_l2 N (inside).
#
# ── Inputs (canonical mode) ───────────────────────────────────────────────
#   --precomp_dir <dir>    auto-resolves precomp + sim_mat (nn80) + sim_mat_l2 (nn40)
#   --L1_dir      <dir>    auto-resolves L1 envelopes + L1 boundaries
#   --L2_dir      <dir>    auto-resolves L2 envelopes + L2 boundaries
#   --chr         <label>
#   --outdir      <dir>
#   [--nn         80]      chrom-wide page NN scale
#   [--nn_l2      40]      inside-segment NN scale
#
# ── Inputs (legacy) ───────────────────────────────────────────────────────
#   --precomp <path> --sim_mat <path> --sim_mat_l2 <path>
#   --l1_catalogue <path> --l1_boundaries <path>
#   --l2_envelopes <path> --l2_boundaries <path>
#
# ── Plot controls ─────────────────────────────────────────────────────────
#   --boundary_filter   stable | non_decay | all   (default: stable)
#   Per-level cosmetics (sizes/colours independently tunable for L1 vs L2):
#       --boundary_number_size_l1   3.5    (L2 default 1.75)
#       --boundary_label_size_l1    2.4    (L2 default 1.2)
#       --boundary_label_color_l1   #0D2A7A (deep navy)
#       --boundary_label_color_l2   #1B5E20 (dark green)
#       --envelope_label_size_l1    1.6    (L2 default 1.0)
#       --envelope_corner_size      1.75
#
# ── Output (in --outdir) ──────────────────────────────────────────────────
#   <chr>.L2_overlay.pdf    multi-page PDF (1 chrom-wide + N segment pages)
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: STEP_D17c_overlay_plot_L2_v8.R)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggnewscale)
  library(scales)
  library(ggrepel)
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ---- CLI --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

precomp_f      <- get_arg("--precomp")
sim_mat_f      <- get_arg("--sim_mat")
sim_mat_l2_f   <- get_arg("--sim_mat_l2", NA_character_)
l1_cat_f       <- get_arg("--l1_catalogue")
l1_bnd_f       <- get_arg("--l1_boundaries")
l2_bnd_f       <- get_arg("--l2_boundaries", NA_character_)
l2_env_f       <- get_arg("--l2_envelopes", NA_character_)
chr_label      <- get_arg("--chr", "chr")
outdir         <- get_arg("--outdir", ".")

# ---- Auto-resolution for L2 overlay -----------------------------------------
# --precomp_dir: resolves precomp + sim_mat (chrom-wide page nn80) +
#                sim_mat_l2 (inside-segment pages nn40).
# --L1_dir / --L2_dir: resolve catalogue + boundary TSVs.
# Tunable via --nn (default 80, chrom-wide) and --nn_l2 (default 40, inside).
precomp_dir    <- get_arg("--precomp_dir")
nn_chrom       <- as.integer(get_arg("--nn",    "80"))
nn_inside      <- as.integer(get_arg("--nn_l2", "40"))
L1_dir         <- get_arg("--L1_dir")
L2_dir         <- get_arg("--L2_dir")

if (is.na(precomp_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, paste0(chr_label, ".precomp.rds"))
  if (file.exists(cand)) { precomp_f <- cand; cat("[plot_L2] resolved --precomp: ", cand, "\n") }
}
if (is.na(sim_mat_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, "sim_mats",
                    paste0(chr_label, ".sim_mat_nn", nn_chrom, ".rds"))
  if (file.exists(cand)) {
    sim_mat_f <- cand
    cat("[plot_L2] resolved --sim_mat (nn=", nn_chrom, "): ", cand, "\n", sep = "")
  }
}
if (is.na(sim_mat_l2_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, "sim_mats",
                    paste0(chr_label, ".sim_mat_nn", nn_inside, ".rds"))
  if (file.exists(cand)) {
    sim_mat_l2_f <- cand
    cat("[plot_L2] resolved --sim_mat_l2 (nn=", nn_inside, "): ", cand, "\n", sep = "")
  }
}
if (is.na(sim_mat_l2_f)) sim_mat_l2_f <- sim_mat_f

if (is.na(l1_cat_f) && !is.na(L1_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(L1_dir, paste0(chr_label, ".L1_envelopes.tsv"))
  if (file.exists(cand)) { l1_cat_f <- cand; cat("[plot_L2] resolved --l1_catalogue: ", cand, "\n") }
}
if (is.na(l1_bnd_f) && !is.na(L1_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(L1_dir, paste0(chr_label, ".L1_boundaries.tsv"))
  if (file.exists(cand)) { l1_bnd_f <- cand; cat("[plot_L2] resolved --l1_boundaries: ", cand, "\n") }
}
if (is.na(l2_env_f) && !is.na(L2_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(L2_dir, paste0(chr_label, ".L2_envelopes.tsv"))
  if (file.exists(cand)) { l2_env_f <- cand; cat("[plot_L2] resolved --l2_envelopes: ", cand, "\n") }
}
if (is.na(l2_bnd_f) && !is.na(L2_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(L2_dir, paste0(chr_label, ".L2_boundaries.tsv"))
  if (file.exists(cand)) { l2_bnd_f <- cand; cat("[plot_L2] resolved --l2_boundaries: ", cand, "\n") }
}

boundary_filter <- tolower(get_arg("--boundary_filter", "stable"))
page_size      <- as.numeric(get_arg("--page_size", "14"))
# Cosmetics — boundary text now styled PER LEVEL.
#   L1 boundary text   = dark blue, larger    (matches L1 v7 sizes)
#   L2 boundary text   = dark green, smaller  (room for crowded L2 pages)
# All four sizes and both colours are independently CLI-tunable.
boundary_number_size_l1 <- as.numeric(get_arg("--boundary_number_size_l1", "3.5"))
boundary_number_size_l2 <- as.numeric(get_arg("--boundary_number_size_l2", "1.75"))
boundary_label_size_l1  <- as.numeric(get_arg("--boundary_label_size_l1",  "2.4"))
boundary_label_size_l2  <- as.numeric(get_arg("--boundary_label_size_l2",  "1.2"))
# Dark blue ~ deep navy reads well on white halo (#0D2A7A). Matches the
# L1 envelope hue family but darker.
boundary_label_color_l1 <- get_arg("--boundary_label_color_l1", "#0D2A7A")
# Dark green from same hue as L2 envelope (#00E676) but darker.
boundary_label_color_l2 <- get_arg("--boundary_label_color_l2", "#1B5E20")
# Envelope title (the long candidate_id rendered above each rect by
# ggrepel). Per-level so L2 can be smaller on crowded per-segment pages.
# L1 default 1.6; L2 default 1.0 (~ 40% smaller).
envelope_label_size_l1 <- as.numeric(get_arg("--envelope_label_size_l1", "1.6"))
envelope_label_size_l2 <- as.numeric(get_arg("--envelope_label_size_l2", "1.0"))
# Bold number rendered inside upper-left of each envelope rectangle, so
# envelopes can be matched to their full label by index. Default 1.75
# (was 3.5 -- the bigger size dominated per-segment pages).
envelope_corner_size <- as.numeric(get_arg("--envelope_corner_size", "1.75"))
# Hover the boundary square by N windows ABOVE the diagonal so it doesn't
# touch the diagonal line. Visual-only. Default 1 (a one-window gap).
draw_W_hover         <- as.integer(get_arg("--draw_W_hover", "1"))
boundary_lw          <- as.numeric(get_arg("--boundary_lw",  "0.4"))
envelope_lw_l1       <- as.numeric(get_arg("--envelope_lw_l1", "0.6"))
envelope_lw_l2       <- as.numeric(get_arg("--envelope_lw_l2", "0.9"))
# L1 envelope outline colour: bright blue #0042FF (was deep navy #1F3A6E).
# L2 envelope outline defaults to electric green so L2 envelopes pop against
# the heatmap and are clearly distinguishable from L1 on combined pages.
envelope_color_l1    <- get_arg("--envelope_color_l1", "#0042FF")
envelope_color_l2    <- get_arg("--envelope_color_l2", "#00E676")
z_clip               <- as.numeric(get_arg("--z_clip", "5"))
# Contrast knobs for the heatmap colour scales:
#   --sim_q_lo / --sim_q_hi : quantile anchors for the similarity colormap.
#     Default 0.05 / 0.95 (matches L1 v7 q05/q95 spread). Lower q_hi -> more
#     contrast on uniformly-high segments.
#   --z_max_min : minimum value for |z| at the colour-scale extreme. The
#     scale collapses on uniform segments where actual max|z| is small;
#     this floor keeps the colour mapping useful. Default 2.5.
sim_q_lo  <- as.numeric(get_arg("--sim_q_lo", "0.05"))
sim_q_hi  <- as.numeric(get_arg("--sim_q_hi", "0.95"))
z_max_min <- as.numeric(get_arg("--z_max_min", "2.5"))

# Whether the L2 boundaries TSV already stores boundary_w at the visible
# center (the new default for STEP_D17_multipass_L2_v8.R with
# --anchor_at_gap_center yes). When TRUE, no visual shift is applied to
# the L2 squares. When FALSE, the legacy +ceiling((G+1)/2) shift is
# applied so the squares land on the visible patch. L1 always uses the
# legacy shift since the L1 multipass writes kernel-anchor coordinates.
l2_already_visible_center <- tolower(get_arg("--l2_already_visible_center", "yes"))
l2_already_visible_center <- l2_already_visible_center %in%
  c("yes", "y", "true", "t", "1", "on")

if (is.na(precomp_f) || !file.exists(precomp_f))
  stop("[D17L2c] --precomp is required and must exist")
if (is.na(sim_mat_f) || !file.exists(sim_mat_f))
  stop("[D17L2c] --sim_mat is required and must exist")
if (is.na(l1_cat_f) || !file.exists(l1_cat_f))
  stop("[D17L2c] --l1_catalogue is required")
if (is.na(l1_bnd_f) || !file.exists(l1_bnd_f))
  stop("[D17L2c] --l1_boundaries is required")
# --l2_boundaries is OPTIONAL. If not given, the per-segment pages just
# show the L1 segment outline with no L2 overlays -- effectively turning
# this script into a multi-page L1 viewer.
have_l2 <- !is.na(l2_bnd_f) && file.exists(l2_bnd_f)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[D17L2c] precomp:        ", precomp_f, "\n", sep = "")
cat("[D17L2c] sim_mat (L1):   ", sim_mat_f, "\n", sep = "")
cat("[D17L2c] sim_mat (L2):   ", sim_mat_l2_f, "\n", sep = "")
cat("[D17L2c] l1_catalogue:   ", l1_cat_f, "\n", sep = "")
cat("[D17L2c] l1_boundaries:  ", l1_bnd_f, "\n", sep = "")
cat("[D17L2c] l2_boundaries:  ",
    if (have_l2) l2_bnd_f else "(none -- multi-page L1-only mode)",
    "\n", sep = "")
cat("[D17L2c] l2_envelopes:   ", if (is.na(l2_env_f)) "(none)" else l2_env_f, "\n", sep = "")
cat("[D17L2c] boundary_filter:", boundary_filter, "\n", sep = "")
cat("[D17L2c] outdir:         ", outdir, "\n", sep = "")
cat("[D17L2c] envelope colours: L1=", envelope_color_l1,
    "  L2=", envelope_color_l2,
    "  (lw L1=", envelope_lw_l1, ", L2=", envelope_lw_l2, ")\n", sep = "")
cat("[D17L2c] boundary text L1: size_num=", boundary_number_size_l1,
    "  size_label=", boundary_label_size_l1,
    "  colour=", boundary_label_color_l1, "\n", sep = "")
cat("[D17L2c] boundary text L2: size_num=", boundary_number_size_l2,
    "  size_label=", boundary_label_size_l2,
    "  colour=", boundary_label_color_l2, "\n", sep = "")
cat("[D17L2c] envelope label:   size L1=", envelope_label_size_l1,
    "  size L2=", envelope_label_size_l2,
    "  corner=", envelope_corner_size,
    "  square hover=", draw_W_hover, "\n", sep = "")
cat("[D17L2c] contrast:       sim quantile anchors q",
    sim_q_lo, "..q", sim_q_hi,
    "  z_max_min=", z_max_min, "\n", sep = "")
cat("[D17L2c] l2_already_visible_center=", l2_already_visible_center,
    "  (L2 squares: ",
    if (l2_already_visible_center) "no extra shift, data already at visible center"
    else "legacy +ceiling((G+1)/2) shift applied",
    ")\n", sep = "")

# ---- Status palette (same as L1 overlay) ------------------------------------

status_palette <- c(
  "STABLE_BLUE" = "#C8102E",  # red
  "DECAYS"      = "#999999",  # gray
  "MARGINAL"    = "#E07B3F",  # orange
  "EDGE"        = "#7A4FBF"   # purple
)

filter_status_set <- switch(boundary_filter,
  "stable"    = c("STABLE_BLUE"),
  "non_decay" = c("STABLE_BLUE", "MARGINAL", "EDGE"),
  "all"       = c("STABLE_BLUE", "DECAYS", "MARGINAL", "EDGE"),
  c("STABLE_BLUE")
)

# ---- Loaders ----------------------------------------------------------------

load_precomp <- function(f) {
  pc <- readRDS(f)
  dt_pc <- if (!is.null(pc$dt)) {
    as.data.table(pc$dt)
  } else if (!is.null(pc$pc) && !is.null(pc$pc$dt)) {
    as.data.table(pc$pc$dt)
  } else {
    stop("[D17L2c] cannot find dt inside precomp")
  }
  if (!is.null(dt_pc$start)) {
    list(start_bp = dt_pc$start, end_bp = dt_pc$end)
  } else if (!is.null(dt_pc$start_bp)) {
    list(start_bp = dt_pc$start_bp, end_bp = dt_pc$end_bp)
  } else if (!is.null(dt_pc$window_start)) {
    list(start_bp = dt_pc$window_start, end_bp = dt_pc$window_end)
  } else if (!is.null(dt_pc$bp_start)) {
    list(start_bp = dt_pc$bp_start, end_bp = dt_pc$bp_end)
  } else {
    stop("[D17L2c] precomp missing window start/end columns")
  }
}

load_sim <- function(f) {
  sm_obj <- readRDS(f)
  if (is.matrix(sm_obj))                   return(sm_obj)
  if (!is.null(sm_obj$sim_mat))            return(sm_obj$sim_mat)
  if (is.list(sm_obj) && length(sm_obj))   return(sm_obj[[1]])
  stop("[D17L2c] sim_mat structure not recognized")
}

cat("[D17L2c] loading precomp\n")
pc_lst <- load_precomp(precomp_f)
window_start_bp <- pc_lst$start_bp
window_end_bp   <- pc_lst$end_bp
N <- length(window_start_bp)
cat("[D17L2c] N windows: ", N, "\n", sep = "")

cat("[D17L2c] loading sim_mat L1\n")
sim_mat_l1 <- load_sim(sim_mat_f)
if (nrow(sim_mat_l1) != N)
  stop("[D17L2c] sim_mat (L1) rows != N windows")

if (sim_mat_l2_f == sim_mat_f) {
  sim_mat_l2 <- sim_mat_l1
} else {
  cat("[D17L2c] loading sim_mat L2\n")
  sim_mat_l2 <- load_sim(sim_mat_l2_f)
  if (nrow(sim_mat_l2) != N)
    stop("[D17L2c] sim_mat (L2) rows != N windows")
}

cat("[D17L2c] loading L1 catalogue + boundaries\n")
l1_envs <- fread(l1_cat_f)
l1_bnd  <- fread(l1_bnd_f)
l1_bnd_keep <- if ("validation_status" %in% names(l1_bnd))
  l1_bnd[validation_status %in% filter_status_set] else l1_bnd
# Tag rows with their level so build_brect can apply the right shift
# convention per row (L1 = legacy kernel anchor + visual shift,
# L2 = visible center if l2_already_visible_center else legacy).
if (nrow(l1_bnd_keep) > 0L) l1_bnd_keep[, level := "L1"]

if (have_l2) {
  cat("[D17L2c] loading L2 boundaries\n")
  l2_bnd <- fread(l2_bnd_f)
  l2_bnd_keep <- if (nrow(l2_bnd) > 0L && "validation_status" %in% names(l2_bnd))
    l2_bnd[validation_status %in% filter_status_set] else l2_bnd[FALSE]
  if (nrow(l2_bnd_keep) > 0L) l2_bnd_keep[, level := "L2"]
} else {
  cat("[D17L2c] L2 boundaries not provided -- multi-page L1-only mode\n")
  l2_bnd      <- data.table()
  l2_bnd_keep <- data.table()
}

l2_envs <- if (have_l2 && !is.na(l2_env_f) && file.exists(l2_env_f))
  fread(l2_env_f) else NULL

# ---- Per-distance Z (local within an x-range) ------------------------------
# Computes a Z(i,j) matrix on a sub-matrix using the sub-matrix's own
# per-distance mean/sd. Returns the Z matrix clipped at +/- z_clip.
compute_local_z <- function(sm, z_clip = 5) {
  Nl <- nrow(sm)
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
    out[cbind(jj, ii)] <- z   # symmetric Z
  }
  out[out >  z_clip] <-  z_clip
  out[out < -z_clip] <- -z_clip
  out
}

# ---- Boundary square geometry helper ---------------------------------------
# Returns the brect data.table with sq_xmin..sq_ymax columns AND the
# anchor_w column. Inputs: a peaks data.table, plus min/max window in the
# plot view (for clipping). draw_W is the visual square edge.
build_brect <- function(bnd_dt, w_min, w_max, draw_W = 5L) {
  if (nrow(bnd_dt) == 0L) return(NULL)
  brect <- copy(bnd_dt)
  # Per-row visual shift: depends on the row's `level` (L1 vs L2) and
  # the global flag l2_already_visible_center.
  #   - L1 rows always get the legacy +ceiling((G+1)/2) shift, since the
  #     L1 multipass writes kernel-anchor boundary_w.
  #   - L2 rows: if l2_already_visible_center == TRUE (default), no shift,
  #     because the L2 multipass with --anchor_at_gap_center yes (default)
  #     already writes visible-center coordinates. If FALSE, apply the
  #     legacy +ceiling((G+1)/2) shift for backward compat with old TSVs.
  #   - Pages that pass a single-stream `bnd` without a `level` column
  #     (page 1 = L1 only, per-segment pages = L2 only) are tagged by
  #     the caller via attr(bnd, "level") OR a `level` column. Default
  #     fallback when neither is set: assume L1 (legacy convention).
  bG <- if ("boundary_offset" %in% names(brect))
    as.integer(brect$boundary_offset[1]) else 5L
  if (!is.finite(bG)) bG <- 5L
  legacy_shift <- as.integer(ceiling((bG + 1L) / 2L))
  if (!"level" %in% names(brect)) {
    inferred_level <- attr(bnd_dt, "level")
    if (is.null(inferred_level) || !nzchar(inferred_level))
      inferred_level <- "L1"
    brect[, level := inferred_level]
  }
  brect[, row_shift := fifelse(
    level == "L2" & isTRUE(l2_already_visible_center),
    0L,
    legacy_shift
  )]
  brect[, anchor_w := boundary_w + row_shift]
  # Y-hover lifts the square off the diagonal by `draw_W_hover` windows,
  # so the square's lower edge sits above the diagonal cell rather than
  # touching it. Visual-only; no effect on data semantics. Default 1.
  hover <- as.integer(draw_W_hover)
  brect[, sq_xmin := anchor_w - draw_W + 0.5]
  brect[, sq_xmax := anchor_w + 0.5]
  brect[, sq_ymin := anchor_w - 0.5 + hover]
  brect[, sq_ymax := anchor_w + draw_W - 0.5 + hover]
  brect[, draw_W_used := draw_W]
  brect <- brect[anchor_w >= w_min & anchor_w <= w_max]
  brect
}

# ---- Page builder -----------------------------------------------------------
# Builds one ggplot for a window range [w_lo, w_hi]. The heatmap below the
# diagonal is `sm[w_lo..w_hi, w_lo..w_hi]` similarity (lower triangle), the
# heatmap above is the local Z (upper triangle). Boundaries to overlay are
# in `bnd` (data.table with boundary_w in absolute coords).
#
# `envs_list` is a list of envelope layers to draw, each a named list:
#   list(envs = <data.table with start_w/end_w/candidate_id>,
#        color = <hex string>,
#        linewidth = <numeric>,
#        level = <string, e.g. "L1" or "L2">)
# Layers are drawn in order, so later layers paint on top.
build_page <- function(sm, w_lo, w_hi, bnd, envs_list, page_title) {
  Nl <- w_hi - w_lo + 1L
  sub <- sm[w_lo:w_hi, w_lo:w_hi]
  # local Z
  z_mat <- compute_local_z(sub, z_clip = z_clip)
  # build long-form: x = row index (varies fast), y = column index (varies slow).
  # Matches the L1 overlay convention: aes(x = i (row), y = j (col)).
  abs_idx <- w_lo:w_hi
  i_row <- rep(abs_idx, times = Nl)   # row, varies fastly (column-major)
  j_col <- rep(abs_idx, each  = Nl)   # col, varies slowly
  sim_v <- as.numeric(sub)
  z_v   <- as.numeric(z_mat)
  # In this convention the diagonal is i_row == j_col. We render the lower
  # triangle (j_col < i_row) as similarity and the upper triangle
  # (j_col > i_row) as local Z, which matches the L1 overlay's split.
  is_lower <- j_col < i_row
  is_upper <- j_col > i_row
  hm_dt <- data.table(x = i_row, y = j_col,
                      lower = ifelse(is_lower, sim_v, NA_real_),
                      upper = ifelse(is_upper, z_v,   NA_real_))
  # Drop fully-NA rows (diagonal) for speed
  hm_dt <- hm_dt[is.finite(lower) | is.finite(upper)]

  # Quantile-based color stretching for the similarity palette so the
  # heatmap reads the same way as the L1 overlay (q05..q98 anchors).
  sim_v_finite <- sim_v[is.finite(sim_v)]
  q_lo <- if (length(sim_v_finite) > 0L)
            as.numeric(quantile(sim_v_finite, sim_q_lo, na.rm = TRUE)) else 0.0
  q_hi <- if (length(sim_v_finite) > 0L)
            as.numeric(quantile(sim_v_finite, sim_q_hi, na.rm = TRUE)) else 1.0
  if (!is.finite(q_lo)) q_lo <- 0.0
  if (!is.finite(q_hi) || q_hi <= q_lo) q_hi <- min(1.0, q_lo + 0.1)
  z_max_local <- max(abs(hm_dt$upper), na.rm = TRUE)
  if (!is.finite(z_max_local) || z_max_local < 1e-3) z_max_local <- 1
  z_max_local <- min(z_max_local, z_clip)
  # Floor to keep the gradient2 mapping useful on uniform segments where
  # actual max|z| is small. Without this floor, on a flat page the scale
  # would map +/-0.3 to deep blue/red and look almost-binary.
  z_max_local <- max(z_max_local, z_max_min)

  p <- ggplot() +
    # Lower triangle: similarity (matches L1 overlay's pale->green->yellow->orange->red)
    geom_raster(data = hm_dt[is.finite(lower)],
                aes(x = x, y = y, fill = lower)) +
    scale_fill_gradientn(
      colours = c("#F8F8F8", "#A8DBC2", "#F2DC78", "#E08838", "#7E1F1F"),
      values  = scales::rescale(c(0, q_lo, (q_lo + q_hi) / 2, q_hi, 1)),
      limits  = c(0, 1), oob = scales::squish,
      name    = "Similarity\n(lower)"
    ) +
    ggnewscale::new_scale_fill() +
    # Upper triangle: Z (matches L1 overlay's deeper blue->white->deeper red)
    geom_raster(data = hm_dt[is.finite(upper)],
                aes(x = x, y = y, fill = upper)) +
    scale_fill_gradient2(
      low      = "#2C5AA0",
      mid      = "#FAFAFA",
      high     = "#B22222",
      midpoint = 0,
      limits   = c(-z_max_local, z_max_local),
      oob      = scales::squish,
      name     = "Z (diagonal)\n(upper)"
    ) +
    # Diagonal as a thin yellow line (drawn via segment for crispness)
    annotate("segment",
             x = w_lo - 0.5, y = w_lo - 0.5,
             xend = w_hi + 0.5, yend = w_hi + 0.5,
             colour = "#E8C547", linewidth = 0.4, alpha = 0.8) +
    coord_fixed(xlim = c(w_lo - 0.5, w_hi + 0.5),
                ylim = c(w_lo - 0.5, w_hi + 0.5),
                expand = FALSE) +
    labs(title = page_title, x = "window index", y = "window index") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          plot.title  = element_text(size = 12, face = "bold"))

  # Envelope outlines + labels — one rendering pass per envelope layer.
  # Drawing order: layer 1 first (typically L1), layer 2 on top (typically
  # L2). Each layer gets its own colour and a ggrepel label matching the
  # L1 v7 envelope-label styling (force=1.5, white halo, etc).
  if (!is.null(envs_list) && length(envs_list) > 0L) {
    for (lyr in envs_list) {
      e_dt <- lyr$envs
      if (is.null(e_dt) || nrow(e_dt) == 0L) next
      e_view <- e_dt[end_w >= w_lo & start_w <= w_hi]
      if (nrow(e_view) == 0L) next
      e_view <- copy(e_view)
      e_view[, draw_start := pmax(start_w, w_lo)]
      e_view[, draw_end   := pmin(end_w,   w_hi)]
      e_view[, lab_x := (draw_start + draw_end) / 2]
      # Place each label slightly above the envelope's upper edge, scaled
      # to the page span (same convention as L1 v7).
      e_view[, lab_y := draw_end + (Nl * 0.015)]
      # Build the label text. Prefer candidate_id if present, else fall
      # back to a synthesised "<level>_<index>" string.
      e_view[, label_text := {
        if ("candidate_id" %in% names(e_view))
          as.character(candidate_id)
        else
          sprintf("%s_%04d", lyr$level, seq_len(.N))
      }]
      # Corner-index marker: bold number at the upper-left INSIDE the
      # envelope rectangle. Reads as a quick legend marker so the user
      # can match the envelope outline to its label by number.
      e_view[, corner_text := {
        # Trailing digits from candidate_id, e.g. "C_gar_LG28_d17L2_0008_03"
        # -> "8_3" or "_0008_03" -> "8_3". Fallback: row index.
        cid <- if ("candidate_id" %in% names(e_view))
                 as.character(candidate_id) else
                 sprintf("%s_%04d", lyr$level, seq_len(.N))
        # Take last two underscore-separated number groups, strip leading zeros
        m <- regmatches(cid, gregexpr("[0-9]+", cid))
        sapply(m, function(nums) {
          if (length(nums) == 0L) return("")
          tail_n <- tail(nums, 2)
          paste(as.integer(tail_n), collapse = "_")
        })
      }]
      # Corner position: just inside the upper-left of the rect, with a
      # small offset proportional to the page span.
      e_view[, corner_x := draw_start + (Nl * 0.012)]
      e_view[, corner_y := draw_end   - (Nl * 0.012)]
      lyr_color <- lyr$color %||% "#1F3A6E"
      lyr_lw    <- lyr$linewidth %||% 0.6
      lyr_level <- lyr$level %||% "L1"
      lyr_label_size <- if (identical(lyr_level, "L2"))
        envelope_label_size_l2 else envelope_label_size_l1
      p <- p +
        geom_rect(
          data = e_view,
          aes(xmin = draw_start - 0.5, xmax = draw_end + 0.5,
              ymin = draw_start - 0.5, ymax = draw_end + 0.5),
          colour    = lyr_color,
          fill      = NA,
          linewidth = lyr_lw,
          linetype  = "solid",
          inherit.aes = FALSE
        ) +
        ggrepel::geom_text_repel(
          data = e_view,
          aes(x = lab_x, y = lab_y, label = label_text),
          colour    = lyr_color,
          size      = lyr_label_size,
          fontface  = "bold",
          bg.colour = "white",
          bg.r      = 0.10,
          segment.size  = lyr_lw,
          segment.alpha = 0.5,
          min.segment.length = 0,
          seed = 1,
          max.overlaps = Inf,
          force = 1.5,
          box.padding   = 0.30,
          point.padding = 0.10,
          inherit.aes = FALSE,
          show.legend = FALSE
        ) +
        # Corner-index marker (geom_text, stable position inside rect)
        geom_text(
          data = e_view,
          aes(x = corner_x, y = corner_y, label = corner_text),
          colour    = lyr_color,
          size      = envelope_corner_size,
          fontface  = "bold",
          hjust     = 0,
          vjust     = 1,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
  }

  # Boundary squares + labels
  brect <- build_brect(bnd, w_min = w_lo, w_max = w_hi)
  if (!is.null(brect) && nrow(brect) > 0L) {
    number_dx <- Nl * 0.005
    text_dx   <- Nl * 0.020
    has_perp   <- "right_frac_blue" %in% names(brect)
    has_grow   <- "grow_max_z"      %in% names(brect)
    brect[, label_text_b := {
      base <- sprintf("%s  %.2f Mb  s=%.1f",
                      boundary_idx, boundary_bp / 1e6, boundary_score)
      if (has_grow) base <- sprintf("%s  gmz=%+.2f", base, grow_max_z)
      if (has_perp) base <- sprintf("%s  rfb=%.2f lfb=%.2f rmz=%+.1f lmz=%+.1f",
                                    base, right_frac_blue, left_frac_blue,
                                    right_max_z, left_max_z)
      base
    }]
    brect[, number_text := {
      raw <- sub("^.*?b0*([0-9]+)$", "\\1", boundary_idx)
      ifelse(grepl("^[0-9]+$", raw), raw, as.character(seq_len(.N)))
    }]
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_rect(data = brect,
                aes(xmin = sq_xmin, xmax = sq_xmax,
                    ymin = sq_ymin, ymax = sq_ymax,
                    fill = validation_status),
                colour = NA, linewidth = 0, inherit.aes = FALSE) +
      scale_fill_manual(
        values = status_palette,
        breaks = c("STABLE_BLUE", "MARGINAL", "EDGE", "DECAYS"),
        name   = "Boundary status",
        drop   = TRUE
      )

    # Per-level boundary text: split brect by level, render with that
    # level's own size + colour. Skip a level if it has no rows.
    brect_l1 <- brect[level == "L1"]
    brect_l2 <- brect[level == "L2"]

    if (nrow(brect_l1) > 0L) {
      p <- p +
        ggrepel::geom_text_repel(
          data = brect_l1,
          aes(x = anchor_w, y = anchor_w, label = label_text_b),
          size      = boundary_label_size_l1,
          colour    = boundary_label_color_l1,
          fontface  = "plain",
          nudge_x   = text_dx,
          direction = "y",
          force     = 0.5,
          force_pull = 0,
          hjust     = 0,
          bg.colour = "white",
          bg.r      = 0.15,
          segment.size  = boundary_lw * 0.6,
          segment.alpha = 0.6,
          min.segment.length = 0,
          max.overlaps = Inf,
          seed = 2,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
    if (nrow(brect_l2) > 0L) {
      p <- p +
        ggrepel::geom_text_repel(
          data = brect_l2,
          aes(x = anchor_w, y = anchor_w, label = label_text_b),
          size      = boundary_label_size_l2,
          colour    = boundary_label_color_l2,
          fontface  = "plain",
          nudge_x   = text_dx,
          direction = "y",
          force     = 0.5,
          force_pull = 0,
          hjust     = 0,
          bg.colour = "white",
          bg.r      = 0.15,
          segment.size  = boundary_lw * 0.6,
          segment.alpha = 0.6,
          min.segment.length = 0,
          max.overlaps = Inf,
          seed = 2,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }

    p <- p + ggnewscale::new_scale_colour()

    # Boundary number: bold, ggrepel with white halo, per-level colour
    # and size. Positioned just to the right of the square's anchor via
    # nudge_x. Same level split as the label text above.
    if (nrow(brect_l1) > 0L) {
      p <- p +
        ggrepel::geom_text_repel(
          data = brect_l1,
          aes(x = anchor_w, y = anchor_w, label = number_text),
          size       = boundary_number_size_l1,
          colour     = boundary_label_color_l1,
          fontface   = "bold",
          nudge_x    = number_dx,
          force      = 0,
          force_pull = 0,
          hjust      = 0,
          bg.colour  = "white",
          bg.r       = 0.20,
          segment.size      = 0,
          min.segment.length = Inf,
          max.overlaps = Inf,
          seed = 3,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
    if (nrow(brect_l2) > 0L) {
      p <- p +
        ggrepel::geom_text_repel(
          data = brect_l2,
          aes(x = anchor_w, y = anchor_w, label = number_text),
          size       = boundary_number_size_l2,
          colour     = boundary_label_color_l2,
          fontface   = "bold",
          nudge_x    = number_dx,
          force      = 0,
          force_pull = 0,
          hjust      = 0,
          bg.colour  = "white",
          bg.r       = 0.20,
          segment.size      = 0,
          min.segment.length = Inf,
          max.overlaps = Inf,
          seed = 3,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
  }
  p
}

# ---- Build all pages --------------------------------------------------------

n_l1 <- nrow(l1_envs)
n_l2_envs <- if (!is.null(l2_envs)) nrow(l2_envs) else 0L
n_l2_bnd  <- nrow(l2_bnd_keep)
cat(sprintf("[D17L2c] L1 segments to render as zoom pages: %d\n", n_l1))

outpath <- file.path(outdir, paste0(chr_label, ".L2_overlay.pdf"))
cat("[D17L2c] writing PDF: ", outpath, "\n", sep = "")
pdf(outpath, width = page_size, height = page_size)

# Helper: stack two boundary tables (L1 and L2) into one rendering set.
# L1 and L2 share the same column schema for the validator-status fields,
# so we only need the columns that build_brect uses downstream.
combine_bnd <- function(b1, b2) {
  if (is.null(b1) || nrow(b1) == 0L) return(b2)
  if (is.null(b2) || nrow(b2) == 0L) return(b1)
  rbindlist(list(b1, b2), fill = TRUE)
}

# Layer 1 = L1 envelopes (dark blue, thinner). Layer 2 = L2 envelopes
# (electric green, thicker, drawn on top).
l1_layer <- list(envs = l1_envs,
                 color = envelope_color_l1, linewidth = envelope_lw_l1,
                 level = "L1")
make_l2_layer <- function(envs_dt) {
  list(envs = envs_dt,
       color = envelope_color_l2, linewidth = envelope_lw_l2,
       level = "L2")
}

# Page 1: whole chromosome at nn80 (or whatever --sim_mat is). L1 only --
# the L2 sim_mat may have different smoothing, and L2 envelopes are only
# meaningful on L2 sim_mat. Keep page 1 clean for the L1 reference view.
title1 <- sprintf("%s | whole chromosome (L1 sim_mat) | L1 envelopes (%d) | L1 STABLE_BLUE boundaries (%d)",
                  chr_label, n_l1, nrow(l1_bnd_keep))
cat("[D17L2c]   page 1: whole chromosome (L1 sim_mat) -- L1 only\n")
p1 <- build_page(sim_mat_l1, 1L, N, l1_bnd_keep,
                 list(l1_layer), title1)
print(p1)

# Page 2: whole chromosome at the L2 sim_mat (e.g. nn40). Now shows
# BOTH L1 and L2 overlays so the user can see the L2 partition against the
# L2 heatmap and compare it to the L1 partition. If sim_mat_l2 == sim_mat
# then page 2 is essentially page 1 but with L2 added on top.
page_offset <- 1L
if (sim_mat_l2_f != sim_mat_f) {
  combined_bnd_p2 <- combine_bnd(l1_bnd_keep, l2_bnd_keep)
  title2 <- sprintf("%s | whole chromosome (L2 sim_mat) | L1 envs (%d) + L2 envs (%d) | L1 STABLE (%d) + L2 STABLE (%d)",
                    chr_label, n_l1, n_l2_envs,
                    nrow(l1_bnd_keep), n_l2_bnd)
  cat("[D17L2c]   page 2: whole chromosome (L2 sim_mat) -- L1 + L2\n")
  envs_layers_p2 <- if (!is.null(l2_envs) && nrow(l2_envs) > 0L)
    list(l1_layer, make_l2_layer(l2_envs)) else list(l1_layer)
  p2 <- build_page(sim_mat_l2, 1L, N, combined_bnd_p2,
                   envs_layers_p2, title2)
  print(p2)
  page_offset <- 2L
}

# Pages (page_offset+1)..(page_offset+n_l1): per-L1-segment zoomed views.
# Each shows: the parent L1 envelope (deep blue, for context), the L2
# envelopes inside it (electric green), and the L2 boundaries inside it.
for (k in seq_len(n_l1)) {
  row <- l1_envs[k]
  s <- as.integer(row$start_w); e <- as.integer(row$end_w)
  pid <- as.character(row$candidate_id)
  l2_in <- if (nrow(l2_bnd_keep) > 0L)
    l2_bnd_keep[parent_l1_id == pid] else l2_bnd_keep
  l2_envs_in <- if (!is.null(l2_envs))
    l2_envs[parent_l1_id == pid] else NULL
  # Just the parent L1 envelope row, so context is preserved as a single
  # outer rectangle on the zoomed page.
  l1_parent_only <- l1_envs[k]
  title_k <- if (have_l2) {
    sprintf("%s | L1 segment %d/%d (%s) | %.2f-%.2f Mb | L2 envs (%d) | L2 boundaries (%d)",
            chr_label, k, n_l1, pid,
            window_start_bp[s] / 1e6, window_end_bp[e] / 1e6,
            if (!is.null(l2_envs_in)) nrow(l2_envs_in) else 0L,
            nrow(l2_in))
  } else {
    sprintf("%s | L1 segment %d/%d (%s) | %.2f-%.2f Mb (zoomed L1 view)",
            chr_label, k, n_l1, pid,
            window_start_bp[s] / 1e6, window_end_bp[e] / 1e6)
  }
  cat(sprintf("[D17L2c]   page %d: %s  [%d..%d]  N=%d\n",
              k + page_offset, pid, s, e, e - s + 1L))
  # Parent L1 envelope drawn in deep blue (context); L2 envelopes inside
  # drawn in electric green (the level we are tuning).
  envs_layers_k <- list(
    list(envs = l1_parent_only, color = envelope_color_l1,
         linewidth = envelope_lw_l1, level = "L1"),
    if (!is.null(l2_envs_in) && nrow(l2_envs_in) > 0L)
      make_l2_layer(l2_envs_in) else NULL
  )
  envs_layers_k <- envs_layers_k[!vapply(envs_layers_k, is.null, logical(1))]
  p_k <- build_page(sim_mat_l2, s, e, l2_in,
                    envs_layers_k, title_k)
  print(p_k)
}

dev.off()
cat("[D17L2c] done\n")
