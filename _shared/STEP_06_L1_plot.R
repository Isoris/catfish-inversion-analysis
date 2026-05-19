#!/usr/bin/env Rscript

# =============================================================================
# STEP_06_L1_plot.R  (path-agnostic; was local_PCA_MDS_z/STEP_ZO_I_plot_L1.R)
#
# Multi-page L1 overlay PDF for one chromosome.
#
# Facultative: short-circuits early when MAKE_L1_PLOTS is not truthy
# (accepts true|1|yes|y|on, case-insensitive; default = true). Set
# MAKE_L1_PLOTS=false in 00_inversion_config.sh to skip L1 plotting
# across all three pipelines without removing the step from launchers.
#
#   Page 1            whole-chromosome heatmap with L1 envelopes outlined
#                     and L1 boundaries marked. GLOBAL per-distance Z
#                     (computed once over the whole chromosome).
#   Pages 2..(N+1)    one zoomed page per L1 segment. LOCAL per-distance Z
#                     computed from the segment's own sub-matrix, so colours
#                     normalize to that segment's distribution.
#
# Heatmap encoding:
#   lower triangle = similarity   (5-color gradient #F8F8F8 -> #7E1F1F)
#   upper triangle = Z            (diverging #2C5AA0 -> #FAFAFA -> #B22222)
#
# Pipeline position:
#   04 detect_L1  -> 05 PLOT_L1
#
# ── Default NN scale: 80 ──────────────────────────────────────────────────
# Matches step 04. Override via --nn N.
#
# ── Inputs (canonical mode) ───────────────────────────────────────────────
#   --precomp_dir <dir>    auto-resolves precomp + sim_mat_nn{N}
#   --L1_dir      <dir>    auto-resolves L1_envelopes.tsv + L1_boundaries.tsv
#   --chr         <label>
#   --outdir      <dir>
#   [--nn         80]
#
# ── Inputs (legacy) ───────────────────────────────────────────────────────
#   --precomp     <path>   --sim_mat   <path>
#   --catalogue   <path>   --boundaries <path>   (the two L1 TSVs)
#
# ── Plot controls ─────────────────────────────────────────────────────────
#   --boundary_filter  stable | non_decay | all     (default: stable)
#   --toggle_L1        yes | no                     hide envelope outlines if no
#   --debug                                         yellow diag + green anchor px
#                                                   on whole-chr page only
#   --page_size        14    --outline_lw 0.6
#   --label_size       1.2   --z_clip     5
#
# ── Style (locked over many sessions) ─────────────────────────────────────
#   L1 envelope outline       deep blue #1F3A6E
#   L1 envelope label         deep blue + white halo (ggrepel)
#   Boundary square fill      status colour (STABLE_BLUE = red,
#                                            DECAYS = grey, etc.)
#   Boundary text label       status colour + halo, ggrepel direction='y'
#   Boundary number           dark red #7A0000 + white halo, no leader
#
# ── Output (in --outdir) ──────────────────────────────────────────────────
#   <chr>.L1_overlay.pdf    multi-page PDF (1 chrom-wide + N segment pages)
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: STEP_D17c_overlay_plot_L1_only_v7.R)
# =============================================================================

.is_true <- function(x) tolower(as.character(x)) %in% c("1", "true", "yes", "y", "on")
.MAKE_L1_PLOTS <- Sys.getenv("MAKE_L1_PLOTS", "true")
if (!.is_true(.MAKE_L1_PLOTS)) {
  message("[STEP_06_L1_plot] MAKE_L1_PLOTS='", .MAKE_L1_PLOTS,
          "' is not truthy — skipping L1 plot stage.")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(scales)
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ---- CLI --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

precomp_f    <- get_arg("--precomp")
sim_mat_f    <- get_arg("--sim_mat")
catalogue_f  <- get_arg("--catalogue")
boundaries_f <- get_arg("--boundaries")  # optional sidecar TSV
chr_label    <- get_arg("--chr", "chr")
outdir       <- get_arg("--outdir", ".")
nn_label     <- get_arg("--nn_label", "nn80")

# ---- Auto-resolution: --precomp_dir + --chr + --nn (default 80 for L1) ------
# Plus --L1_dir for catalogue/boundaries auto-resolution.
precomp_dir  <- get_arg("--precomp_dir")
nn_scale     <- as.integer(get_arg("--nn", "80"))
L1_dir       <- get_arg("--L1_dir")

if (is.na(precomp_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, paste0(chr_label, ".precomp.rds"))
  if (file.exists(cand)) {
    precomp_f <- cand
    cat("[plot_L1] resolved --precomp: ", cand, "\n")
  }
}
if (is.na(sim_mat_f) && !is.na(precomp_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(precomp_dir, "sim_mats",
                    paste0(chr_label, ".sim_mat_nn", nn_scale, ".rds"))
  if (file.exists(cand)) {
    sim_mat_f <- cand
    cat("[plot_L1] resolved --sim_mat (nn=", nn_scale, "): ", cand, "\n", sep = "")
  }
}
if (is.na(catalogue_f) && !is.na(L1_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(L1_dir, paste0(chr_label, ".L1_envelopes.tsv"))
  if (file.exists(cand)) {
    catalogue_f <- cand
    cat("[plot_L1] resolved --catalogue: ", cand, "\n")
  }
}
if (is.na(boundaries_f) && !is.na(L1_dir) && !is.na(chr_label) && chr_label != "chr") {
  cand <- file.path(L1_dir, paste0(chr_label, ".L1_boundaries.tsv"))
  if (file.exists(cand)) {
    boundaries_f <- cand
    cat("[plot_L1] resolved --boundaries: ", cand, "\n")
  }
}
page_size    <- as.numeric(get_arg("--page_size", "14"))
outline_lw   <- as.numeric(get_arg("--outline_lw", "0.6"))
label_size   <- as.numeric(get_arg("--label_size", "1.2"))
z_clip       <- as.numeric(get_arg("--z_clip", "5"))
boundary_lw    <- as.numeric(get_arg("--boundary_lw", "0.4"))

# Boundary status filter & colors. Available statuses:
#   STABLE_BLUE  -- real long-range separation, persists at all val_W
#   DECAYS       -- trigger fired but val_score collapses as W grows
#   MARGINAL     -- ambiguous
#   EDGE         -- too close to chr edge to validate
# Filter values:
#   "stable"  draws only STABLE_BLUE       (recommended default)
#   "all"     draws every status (each in its own color)
#   "non_decay" draws everything except DECAYS
boundary_filter <- get_arg("--boundary_filter", "stable")
boundary_filter <- tolower(boundary_filter)

# v7: --toggle_L1 yes/no (default yes). When "no", L1 envelope outlines
# and labels are not drawn.
toggle_L1_str   <- tolower(get_arg("--toggle_L1", "yes"))
toggle_L1       <- toggle_L1_str %in% c("yes", "y", "true", "t", "1")

# v7: --boundary_show_all_status yes/no (default no). When yes, override
# --boundary_filter and draw every status in its own color. Lets the user
# see what was killed and why.
show_all_str    <- tolower(get_arg("--boundary_show_all_status", "no"))
show_all_status <- show_all_str %in% c("yes", "y", "true", "t", "1")

# v7: --debug yes/no (default no). When yes, draws a yellow diagonal line
# and green anchor dots at each STABLE_BLUE boundary's (boundary_w,
# boundary_w) pixel. Lets the user verify visually that the anchor pixel
# is on the matrix diagonal.
debug_str       <- tolower(get_arg("--debug", "no"))
debug_mode      <- debug_str %in% c("yes", "y", "true", "t", "1")

# v7: --dedup_overlap yes/no (default yes). When yes, after the status
# filter is applied, pairs of surviving boundaries closer than W windows
# apart are de-duplicated by keeping the one with the larger boundary_score
# (greedy, by descending score).
dedup_str       <- tolower(get_arg("--dedup_overlap", "yes"))
dedup_overlap   <- dedup_str %in% c("yes", "y", "true", "t", "1")

# v7: --boundary_style fill (default) | outline. Filled squares stay
# visible even at small W; outlined matches the older v6 look.
boundary_style  <- tolower(get_arg("--boundary_style", "fill"))
if (!boundary_style %in% c("fill", "outline")) {
  cat("[D17L1c] WARN: --boundary_style must be 'fill' or 'outline'; got '",
      boundary_style, "'. Falling back to 'fill'.\n", sep = "")
  boundary_style <- "fill"
}

# v7: --boundary_number yes (default) | no. Draw the index number of
# each surviving boundary in the center of its square. The number is the
# trailing digits of boundary_idx, with leading zeros stripped (e.g.
# "C_gar_LG28_d17L1_b0076" -> "76").
bnum_str        <- tolower(get_arg("--boundary_number", "yes"))
boundary_number <- bnum_str %in% c("yes", "y", "true", "t", "1")

# v7: --boundary_number_size (default 1.4). geom_text size for the number
# inside each square.
boundary_number_size <- as.numeric(get_arg("--boundary_number_size", "3.5"))

status_palette <- c(
  STABLE_BLUE = "#C8102E",   # red — real boundary
  DECAYS      = "#999999",   # gray — fake / inside an inversion
  MARGINAL    = "#E07B3F",   # orange — ambiguous
  EDGE        = "#7A4FBF"    # purple — chromosome edge
)

if (is.na(precomp_f)   || !file.exists(precomp_f))   stop("[D17L1c] --precomp required")
if (is.na(catalogue_f) || !file.exists(catalogue_f)) stop("[D17L1c] --catalogue required")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[D17L1c] precomp:   ", precomp_f, "\n")
cat("[D17L1c] sim_mat:   ", sim_mat_f %||% "(auto)", "\n")
cat("[D17L1c] catalogue: ", catalogue_f, "\n")
cat("[D17L1c] boundaries:", boundaries_f %||% "(none)", "\n")
cat("[D17L1c] boundary_filter: ", boundary_filter, "\n", sep = "")
cat("[D17L1c] toggle_L1:        ", if (toggle_L1) "yes" else "no", "\n", sep = "")
cat("[D17L1c] show_all_status:  ", if (show_all_status) "yes" else "no", "\n", sep = "")
cat("[D17L1c] debug:            ", if (debug_mode) "yes" else "no", "\n", sep = "")
cat("[D17L1c] dedup_overlap:    ", if (dedup_overlap) "yes" else "no", "\n", sep = "")
cat("[D17L1c] boundary_style:   ", boundary_style, "\n", sep = "")
cat("[D17L1c] boundary_number:  ", if (boundary_number) "yes" else "no",
    "  size=", boundary_number_size, "\n", sep = "")
cat("[D17L1c] chr:       ", chr_label, "\n")
cat("[D17L1c] outdir:    ", outdir, "\n")
cat("[D17L1c] page_size: ", page_size, "in\n")

# ---- Load precomp -----------------------------------------------------------

cat("[D17L1c] loading precomp\n")
pc_obj <- readRDS(precomp_f)
if (!is.null(pc_obj$dt)) {
  pc <- pc_obj
} else if (!is.null(pc_obj$pc)) {
  pc <- pc_obj$pc
} else stop("[D17L1c] cannot find pc$dt in precomp file")
dt_pc <- as.data.table(pc$dt)
n_windows_total <- nrow(dt_pc)

# Coordinate columns
if (all(c("start", "end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start
  window_end_bp   <- dt_pc$end
} else if (all(c("start_bp", "end_bp") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start_bp
  window_end_bp   <- dt_pc$end_bp
} else if (all(c("window_start", "window_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$window_start
  window_end_bp   <- dt_pc$window_end
} else if (all(c("bp_start", "bp_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$bp_start
  window_end_bp   <- dt_pc$bp_end
} else {
  stop("[D17L1c] precomp dt missing usable coordinate columns. Found: ",
       paste(names(dt_pc), collapse = ", "))
}
cat("[D17L1c] N windows: ", n_windows_total, "\n")

# ---- Load sim_mat -----------------------------------------------------------

if (is.na(sim_mat_f)) {
  candidates_path <- c(
    file.path(dirname(precomp_f), "sim_mat_nn80.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn160.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn40.rds")
  )
  hit <- candidates_path[file.exists(candidates_path)]
  if (length(hit) == 0) stop("[D17L1c] could not auto-find sim_mat_nn*.rds")
  sim_mat_f <- hit[1]
}
cat("[D17L1c] loading sim_mat: ", sim_mat_f, "\n")
sm_obj <- readRDS(sim_mat_f)
if (is.matrix(sm_obj)) {
  sim_mat <- sm_obj
} else if (is.list(sm_obj) && !is.null(sm_obj$sim_mat)) {
  sim_mat <- sm_obj$sim_mat
} else if (is.list(sm_obj) && length(sm_obj) == 1) {
  sim_mat <- sm_obj[[1]]
} else stop("[D17L1c] sim_mat structure not recognized")
storage.mode(sim_mat) <- "double"

if (!isTRUE(nrow(sim_mat) == n_windows_total &&
            ncol(sim_mat) == n_windows_total)) {
  stop("[D17L1c] sim_mat dim ", nrow(sim_mat), "x", ncol(sim_mat),
       " does not match N windows ", n_windows_total)
}

# ---- Load catalogue ---------------------------------------------------------

cat("[D17L1c] loading catalogue\n")
cat_dt <- fread(catalogue_f)
cat("[D17L1c] catalogue rows: ", nrow(cat_dt), "\n", sep = "")

# Filter to envelope rows only
cat_dt <- cat_dt[status == "ENVELOPE"]
cat_dt <- cat_dt[!is.na(start_w) & !is.na(end_w)]
cat("[D17L1c] envelopes to draw: ", nrow(cat_dt), "\n", sep = "")

# ---- Load boundaries (optional) --------------------------------------------
bdt <- data.table()
if (!is.na(boundaries_f) && file.exists(boundaries_f)) {
  bdt <- fread(boundaries_f)
  cat("[D17L1c] boundaries to draw: ", nrow(bdt), "\n", sep = "")
} else if (!is.na(boundaries_f)) {
  cat("[D17L1c] boundaries file not found: ", boundaries_f,
      " (skipping)\n", sep = "")
}

# ---- Boundary preprocessing (filter + dedup + brect construction) ----------
# Done ONCE up front so all pages render against the same surviving set of
# boundaries. Per-page plotting just filters brect to the page's window
# range; no re-classification or re-dedup happens per page.
if (nrow(bdt) > 0L) {
  if ("validation_status" %in% names(bdt)) {
    n_before <- nrow(bdt)
    if (show_all_status) {
      keep_set <- c("STABLE_BLUE", "DECAYS", "MARGINAL", "EDGE")
      cat("[D17L1c] show_all_status=yes -- drawing every status\n")
    } else {
      keep_set <- switch(boundary_filter,
        "stable"    = c("STABLE_BLUE"),
        "non_decay" = c("STABLE_BLUE", "MARGINAL", "EDGE"),
        "all"       = c("STABLE_BLUE", "DECAYS", "MARGINAL", "EDGE"),
        c("STABLE_BLUE")
      )
    }
    bdt <- bdt[validation_status %in% keep_set]
    cat("[D17L1c] boundary_filter='", boundary_filter,
        "': ", n_before, " -> ", nrow(bdt), " boundaries\n", sep = "")
  } else {
    cat("[D17L1c] no validation_status column -- drawing all boundaries\n")
    bdt[, validation_status := "STABLE_BLUE"]
  }

  # De-overlap: greedy by boundary_score, drop peaks within W of a kept one.
  if (dedup_overlap && !show_all_status && nrow(bdt) > 1L) {
    n_before_dedup <- nrow(bdt)
    W_dedup <- if ("boundary_W" %in% names(bdt))
                 as.integer(bdt$boundary_W[1]) else 5L
    bdt_sorted <- bdt[order(-boundary_score)]
    kept_idx <- integer(0); kept_w <- integer(0)
    for (k in seq_len(nrow(bdt_sorted))) {
      bw_k <- bdt_sorted$boundary_w[k]
      if (length(kept_w) == 0L ||
          all(abs(kept_w - bw_k) >= W_dedup)) {
        kept_idx <- c(kept_idx, k)
        kept_w   <- c(kept_w, bw_k)
      } else {
        cat(sprintf(
          "[D17L1c]   dedup: drop %s (w=%d, s=%.2f) -- overlaps kept boundary at w=%d\n",
          bdt_sorted$boundary_idx[k], bw_k, bdt_sorted$boundary_score[k],
          kept_w[which.min(abs(kept_w - bw_k))]))
      }
    }
    bdt <- bdt_sorted[kept_idx][order(boundary_w)]
    cat("[D17L1c] dedup_overlap (W=", W_dedup, "): ",
        n_before_dedup, " -> ", nrow(bdt), " boundaries\n", sep = "")
  }
}

# ---- Boundary square geometry (brect) construction --------------------------
# Builds the brect data.table with anchor coords, square coords, and label
# columns. brect is reused across all pages; per-page rendering filters it
# to anchors within the page's window range.
brect <- data.table()
if (nrow(bdt) > 0L) {
  bW_trig <- if ("boundary_W" %in% names(bdt))
               as.integer(bdt$boundary_W[1]) else 5L
  bG <- if ("boundary_offset" %in% names(bdt))
          as.integer(bdt$boundary_offset[1]) else 5L

  brect <- copy(bdt)
  brect[, draw_W := bW_trig]
  shift_to_gap_center <- as.integer(ceiling((bG + 1L) / 2L))
  brect[, anchor_w := boundary_w + shift_to_gap_center]
  cat("[D17L1c] gap-center shift (boundary_offset=", bG, "): +",
      shift_to_gap_center, " px on diagonal\n", sep = "")

  brect[, sq_xmin := anchor_w - draw_W + 0.5]
  brect[, sq_xmax := anchor_w + 0.5]
  brect[, sq_ymin := anchor_w - 0.5]
  brect[, sq_ymax := anchor_w + draw_W - 0.5]
  # Clip to plot area
  brect <- brect[sq_xmax <= n_windows_total + 0.5 &
                 sq_ymax <= n_windows_total + 0.5 &
                 sq_xmin >= 0.5 & sq_ymin >= 0.5]

  has_perp   <- "right_frac_blue" %in% names(brect)
  has_inside <- "inside_L1"       %in% names(brect)
  has_grow   <- "grow_max_z"      %in% names(brect)
  brect[, label_text_b := {
    base <- sprintf("%s  %.2f Mb  s=%.1f",
                    boundary_idx, boundary_bp / 1e6, boundary_score)
    if (has_grow) base <- sprintf("%s  gmz=%+.2f", base, grow_max_z)
    if (has_perp) base <- sprintf("%s  rfb=%.2f lfb=%.2f rmz=%+.1f lmz=%+.1f",
                                  base, right_frac_blue, left_frac_blue,
                                  right_max_z, left_max_z)
    if (has_inside)  base <- ifelse(inside_L1,
                                    sprintf("%s  IN-L1", base),
                                    base)
    base
  }]
  brect[, status_color := unname(status_palette[validation_status])]
  brect[, number_text := {
    raw <- sub("^.*?b0*([0-9]+)$", "\\1", boundary_idx)
    ifelse(grepl("^[0-9]+$", raw), raw, as.character(seq_len(.N)))
  }]
}

# ---- Diagonal-distance Z ----------------------------------------------------

cat("[D17L1c] computing diagonal-distance Z\n")
N <- n_windows_total
z_mat <- matrix(0, nrow = N, ncol = N)
for (d in 0:(N - 1L)) {
  idx_i <- seq_len(N - d)
  idx_j <- idx_i + d
  idx_upper <- cbind(idx_i, idx_j)
  idx_lower <- cbind(idx_j, idx_i)
  vals_upper <- sim_mat[idx_upper]
  vals_lower <- sim_mat[idx_lower]
  vals <- if (d == 0L) vals_upper else c(vals_upper, vals_lower)
  mu_d <- mean(vals, na.rm = TRUE)
  sd_d <- sd(vals, na.rm = TRUE)
  if (!is.finite(sd_d) || sd_d < 1e-9) sd_d <- 1
  z_mat[idx_upper] <- (vals_upper - mu_d) / sd_d
  if (d > 0L) z_mat[idx_lower] <- (vals_lower - mu_d) / sd_d
  if (d > 0L && d %% 1000L == 0L) {
    cat(sprintf("[D17L1c]   diag d=%d / %d\n", d, N - 1L))
  }
}
zr <- range(z_mat, na.rm = TRUE)
cat(sprintf("[D17L1c] Z range: [%.2f, %.2f]; clipping at +/- %.1f\n",
            zr[1], zr[2], z_clip))
z_mat[z_mat >  z_clip] <-  z_clip
z_mat[z_mat < -z_clip] <- -z_clip

# ---- Per-segment local Z ----------------------------------------------------
# Used on the zoomed per-L1-segment pages so that the heatmap colors are
# normalized to the segment's OWN per-distance distribution. Without this,
# small segments would render mostly desaturated because the global per-
# distance mean/sd was computed across all 4302 windows.
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

# ---- Build heatmap long-form (whole genome) --------------------------------

cat("[D17L1c] building long-form heatmap (",
    format(N * N, big.mark = ","), " cells)\n", sep = "")

# Note: for very large N this can be RAM-heavy. With N=4302, ~18.5M rows.
ii <- rep(seq_len(N), times = N)
jj <- rep(seq_len(N), each = N)
lower_mask <- jj <= ii
upper_mask <- jj >  ii

hm <- data.table(
  i       = ii,
  j       = jj,
  sim_val = ifelse(lower_mask, as.numeric(sim_mat), NA_real_),
  z_val   = ifelse(upper_mask, as.numeric(z_mat),   NA_real_)
)

sim_layer <- hm[!is.na(sim_val)]
z_layer   <- hm[!is.na(z_val)]

q_lo <- as.numeric(quantile(sim_layer$sim_val, 0.05, na.rm = TRUE))
q_hi <- as.numeric(quantile(sim_layer$sim_val, 0.98, na.rm = TRUE))

z_max <- max(abs(z_layer$z_val), na.rm = TRUE)
z_max <- min(z_max, z_clip)
if (!is.finite(z_max) || z_max < 1e-3) z_max <- 1

# ---- Envelope rectangles + labels -------------------------------------------

if (nrow(cat_dt) > 0L) {
  cat_dt[, draw_start := pmax(start_w, 1L)]
  cat_dt[, draw_end   := pmin(end_w, n_windows_total)]
  cat_dt[, label_text := sprintf("%s  %.2f-%.2f Mb (nW=%d)",
                                 candidate_id,
                                 start_bp / 1e6, end_bp / 1e6,
                                 n_windows)]
}

# ---- Build heatmap layers helper -------------------------------------------
# Builds the lower-triangle similarity layer and upper-triangle Z layer for
# a window range [w_lo, w_hi]. Returns a list with sim_layer and z_layer.
# When use_local_z = TRUE (per-segment pages), Z is recomputed from the
# sub-matrix using compute_local_z so colors normalize to the segment's
# own per-distance distribution.
build_heatmap_layers <- function(w_lo, w_hi, use_local_z) {
  if (w_lo == 1L && w_hi == n_windows_total && !use_local_z) {
    # Fast path for the whole-chr page: reuse the global sim_layer/z_layer
    # if they already exist in scope; otherwise build them fresh.
    return(list(sim_layer = sim_layer, z_layer = z_layer))
  }
  Nl <- w_hi - w_lo + 1L
  abs_idx <- w_lo:w_hi
  sm_sub <- sim_mat[abs_idx, abs_idx]
  if (use_local_z) {
    z_sub <- compute_local_z(sm_sub, z_clip = z_clip)
  } else {
    # Use the global z_mat slice (same per-distance normalization as page 1)
    z_sub <- z_mat[abs_idx, abs_idx]
  }
  # Long-form: x = row index (varies fast), y = col index (varies slow)
  ii <- rep(abs_idx, times = Nl)
  jj <- rep(abs_idx, each = Nl)
  lower_mask <- jj <= ii
  upper_mask <- jj >  ii
  hm <- data.table(
    i       = ii,
    j       = jj,
    sim_val = ifelse(lower_mask, as.numeric(sm_sub), NA_real_),
    z_val   = ifelse(upper_mask, as.numeric(z_sub),  NA_real_)
  )
  list(
    sim_layer = hm[!is.na(sim_val)],
    z_layer   = hm[!is.na(z_val)]
  )
}

# ---- Per-page axis breaks --------------------------------------------------
build_axis <- function(w_lo, w_hi) {
  br <- pretty(c(w_lo, w_hi), n = 8L)
  br <- unique(as.integer(round(br)))
  br <- br[br >= w_lo & br <= w_hi]
  lab <- sprintf("%.2f",
    window_start_bp[pmin(pmax(br, 1L), n_windows_total)] / 1e6
  )
  list(breaks = br, labels = lab)
}

# ---- The page builder ------------------------------------------------------
# Renders one page: heatmap (sim lower / Z upper), L1 envelopes (filtered to
# overlap with the window range), L1 boundary squares (filtered to anchors
# in range), axes, theme, title. The debug overlay is drawn only on the
# whole-chr page (where it makes sense visually).
build_page <- function(w_lo, w_hi, page_title, use_local_z, draw_debug) {
  layers <- build_heatmap_layers(w_lo, w_hi, use_local_z)
  sl <- layers$sim_layer
  zl <- layers$z_layer
  axes <- build_axis(w_lo, w_hi)

  cat(sprintf("[D17L1c]   page: [%d..%d] N=%d  local_Z=%s\n",
              w_lo, w_hi, w_hi - w_lo + 1L,
              if (use_local_z) "yes" else "no"))

  p <- ggplot() +
    # Lower triangle: similarity
    geom_raster(
      data = sl,
      aes(x = i, y = j, fill = sim_val)
    ) +
    scale_fill_gradientn(
      colours = c("#F8F8F8", "#A8DBC2", "#F2DC78", "#E08838", "#7E1F1F"),
      values  = scales::rescale(c(0, q_lo, (q_lo + q_hi) / 2, q_hi, 1)),
      limits  = c(0, 1),
      oob     = scales::squish,
      name    = "Similarity\n(lower)"
    ) +
    ggnewscale::new_scale_fill() +
    # Upper triangle: Z
    geom_raster(
      data = zl,
      aes(x = i, y = j, fill = z_val)
    ) +
    scale_fill_gradient2(
      low      = "#2C5AA0",
      mid      = "#FAFAFA",
      high     = "#B22222",
      midpoint = 0,
      limits   = c(-z_max, z_max),
      oob      = scales::squish,
      name     = "Z (diagonal)\n(upper)"
    )

  # Envelope outlines + labels (filtered to overlap with the page's range)
  if (toggle_L1 && nrow(cat_dt) > 0L) {
    cat_view <- cat_dt[end_w >= w_lo & start_w <= w_hi]
    if (nrow(cat_view) > 0L) {
      cat_view <- copy(cat_view)
      cat_view[, draw_start := pmax(start_w, w_lo)]
      cat_view[, draw_end   := pmin(end_w,   w_hi)]
      cat_view[, lab_x := (draw_start + draw_end) / 2]
      # Place labels above each envelope's upper edge, scaled to the page span
      cat_view[, lab_y := draw_end + ((w_hi - w_lo + 1L) * 0.015)]
      p <- p +
        geom_rect(
          data = cat_view,
          aes(xmin = draw_start - 0.5, xmax = draw_end + 0.5,
              ymin = draw_start - 0.5, ymax = draw_end + 0.5),
          colour    = "#1F3A6E",
          fill      = NA,
          linewidth = outline_lw,
          linetype  = "solid",
          inherit.aes = FALSE
        ) +
        ggrepel::geom_text_repel(
          data = cat_view,
          aes(x = lab_x, y = lab_y, label = label_text),
          colour = "#1F3A6E",
          size = label_size,
          fontface = "bold",
          bg.colour = "white",
          bg.r = 0.10,
          segment.size = outline_lw,
          segment.alpha = 0.5,
          min.segment.length = 0,
          seed = 1,
          max.overlaps = Inf,
          force = 1.5,
          box.padding = 0.30,
          point.padding = 0.10,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
  }

  # Boundary squares (filtered to anchors within the page's window range)
  if (nrow(brect) > 0L) {
    brect_view <- brect[anchor_w >= w_lo & anchor_w <= w_hi]
    if (nrow(brect_view) > 0L) {
      # Page-local label nudges (so spacing scales with zoom)
      number_dx <- (w_hi - w_lo + 1L) * 0.005
      text_dx   <- (w_hi - w_lo + 1L) * 0.020

      if (boundary_style == "fill") {
        p <- p +
          ggnewscale::new_scale_fill() +
          geom_rect(
            data = brect_view,
            aes(xmin = sq_xmin, xmax = sq_xmax,
                ymin = sq_ymin, ymax = sq_ymax,
                fill = validation_status),
            colour      = NA,
            linewidth   = 0,
            inherit.aes = FALSE
          ) +
          scale_fill_manual(
            values = status_palette,
            breaks = c("STABLE_BLUE", "MARGINAL", "EDGE", "DECAYS"),
            name   = "Boundary status",
            drop   = TRUE
          ) +
          ggrepel::geom_text_repel(
            data = brect_view,
            aes(x = anchor_w, y = anchor_w, label = label_text_b,
                colour = validation_status),
            size      = label_size,
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
          ) +
          scale_colour_manual(
            values = status_palette,
            breaks = c("STABLE_BLUE", "MARGINAL", "EDGE", "DECAYS"),
            name   = "Boundary status",
            drop   = TRUE,
            guide  = "none"
          )
      } else {
        p <- p +
          geom_rect(
            data = brect_view,
            aes(xmin = sq_xmin, xmax = sq_xmax,
                ymin = sq_ymin, ymax = sq_ymax,
                colour = validation_status),
            fill        = NA,
            linewidth   = boundary_lw,
            linetype    = "solid",
            inherit.aes = FALSE
          ) +
          ggrepel::geom_text_repel(
            data = brect_view,
            aes(x = anchor_w, y = anchor_w, label = label_text_b,
                colour = validation_status),
            size      = label_size,
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
          ) +
          scale_colour_manual(
            values = status_palette,
            breaks = c("STABLE_BLUE", "MARGINAL", "EDGE", "DECAYS"),
            name   = "Boundary status",
            drop   = TRUE
          )
      }

      # Boundary number: ggrepel with white halo, bold dark red, positioned
      # just to the right of the square's anchor via nudge_x.
      if (boundary_number) {
        p <- p +
          ggnewscale::new_scale_colour() +
          ggrepel::geom_text_repel(
            data = brect_view,
            aes(x = anchor_w, y = anchor_w, label = number_text),
            size       = boundary_number_size,
            colour     = "#7A0000",
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
  }

  # Debug overlay: only on the whole-chr page
  if (draw_debug && debug_mode) {
    cat("[D17L1c] DEBUG MODE: drawing yellow diagonal + green anchor pixels\n")
    diag_seg <- data.table(x = w_lo, y = w_lo, xend = w_hi, yend = w_hi)
    p <- p +
      ggnewscale::new_scale_colour() +
      geom_segment(
        data = diag_seg,
        aes(x = x, y = y, xend = xend, yend = yend),
        colour    = "#FFD700",
        linewidth = 0.12,
        linetype  = "solid",
        inherit.aes = FALSE
      )
    if (nrow(brect) > 0L) {
      anchors <- brect[, .(ax = anchor_w, ay = anchor_w)]
      p <- p +
        geom_tile(
          data = anchors,
          aes(x = ax, y = ay),
          width  = 1, height = 1,
          fill   = "#00FF00",
          colour = "#003300",
          linewidth = 0.15,
          inherit.aes = FALSE
        )
    }
  }

  # Axes + theme + title
  p <- p +
    coord_equal(xlim = c(w_lo - 0.5, w_hi + 0.5),
                ylim = c(w_lo - 0.5, w_hi + 0.5),
                expand = FALSE) +
    scale_x_continuous(
      expand = c(0, 0),
      name = "window index",
      breaks = axes$breaks,
      sec.axis = sec_axis(
        ~ .,
        breaks = axes$breaks,
        labels = axes$labels,
        name = "Mb"
      )
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      name = "window index",
      breaks = axes$breaks,
      sec.axis = sec_axis(
        ~ .,
        breaks = axes$breaks,
        labels = axes$labels,
        name = "Mb"
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1.2, "cm"),
      plot.title = element_text(size = 12, face = "bold"),
      axis.line  = element_line(colour = "black", linewidth = 0.25),
      axis.ticks = element_line(colour = "black", linewidth = 0.25)
    ) +
    labs(title = page_title)

  p
}

# ---- Render multi-page PDF -------------------------------------------------

n_l1 <- nrow(cat_dt)
cat(sprintf("[D17L1c] multipage rendering: 1 whole-chr page + %d zoomed L1 segment pages\n",
            n_l1))

outname <- paste0(chr_label, ".L1_overlay",
                  if (debug_mode) "_DEBUG" else "",
                  ".pdf")
outpath <- file.path(outdir, outname)
cat("[D17L1c] writing PDF: ", outpath, "\n", sep = "")
pdf(outpath, width = page_size, height = page_size)

# Page 1: whole chromosome (uses global Z)
title1 <- sprintf(
  "%s | %s | whole chromosome | %s%s%s",
  chr_label, nn_label,
  if (toggle_L1) sprintf("L1 envelopes (%d)", n_l1) else "L1 envelopes hidden",
  if (nrow(brect) > 0L)
    sprintf(" | L1 boundaries (%d%s)", nrow(brect),
            if (show_all_status) ", all status" else "")
  else "",
  if (debug_mode) " | DEBUG (yellow=diagonal, green=anchor px)" else ""
)
cat("[D17L1c] page 1: whole chromosome\n")
p1 <- build_page(1L, n_windows_total, title1,
                 use_local_z = FALSE, draw_debug = TRUE)
print(p1)

# Pages 2..N+1: per-L1-segment zoomed views (local Z for sharper contrast)
if (n_l1 > 0L) {
  for (k in seq_len(n_l1)) {
    row <- cat_dt[k]
    s <- as.integer(row$start_w); e <- as.integer(row$end_w)
    pid <- as.character(row$candidate_id)
    title_k <- sprintf(
      "%s | %s | L1 segment %d/%d (%s) | %.2f-%.2f Mb (nW=%d)",
      chr_label, nn_label, k, n_l1, pid,
      window_start_bp[s] / 1e6, window_end_bp[e] / 1e6,
      e - s + 1L
    )
    cat(sprintf("[D17L1c] page %d: %s [%d..%d]\n", k + 1L, pid, s, e))
    p_k <- build_page(s, e, title_k,
                      use_local_z = TRUE, draw_debug = FALSE)
    print(p_k)
  }
}

dev.off()
cat("[D17L1c] done\n")
