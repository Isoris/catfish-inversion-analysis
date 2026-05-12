#!/usr/bin/env Rscript
# =============================================================================
# STEP_E_refine_boundaries.R  (v1, 2026-05-12)
# =============================================================================
# Per-D17-boundary refinement via two-sided per-sample CUSUM. Builds a
# small breakpoint graph (nodes = bp positions, edges = sample-set links)
# rather than collapsing everything into a single bp number, because real
# refinement windows often contain multiple events (paired flanks of one
# inversion, nested inversions, partial overlaps) and different karyotype
# groups (HET, HOM_INV, HOM_REF) show different *directions* at the same
# physical breakpoint, not different bp positions.
#
# PATH-AGNOSTIC. Works on both:
#   - θπ:   theta_native.<CHR>.<scale>.tsv.gz feature data
#   - GHSL: <CHR>.ghsl_matrices.rds feature data
# Loader chosen by --label {theta,ghsl}. Algorithm identical for both.
#
# DESIGN PRINCIPLE: lean toward NOT collapsing structure. Defaults are
# slightly lenient (more peaks, wider clustering, more edges kept) because
# downstream merging is easy and downstream un-merging is impossible.
#
# Pipeline position:
#   D17 L1 detect  ->  E_refine_boundaries  ->  atlas + manuscript
#
# Inputs:
#   --boundaries_tsv <path>   D17 L1_boundaries.tsv (uses STABLE_BLUE rows)
#   --precomp_dir    <path>   for window grid (start_bp, end_bp) per chrom
#   --feature_source <path>   theta TSV dir, OR ghsl_matrices dir
#   --label          theta|ghsl   selects loader and per-path defaults
#   --scale          <str>    feature scale to refine on
#                             theta: pestPG scale label (e.g.
#                                    "win10000.step2000" = dense)
#                             ghsl:  rolling key (e.g. "raw","s10","s50")
#   --outdir         <path>
#   [--chr <CHR>]              optional: refine only one chrom's boundaries
#   [--precomp_suffix <s>]     "precomp.rds" (theta) or "ghsl_precomp.rds"
#
# Tuning knobs (all have per-label defaults; override as needed):
#   --window_kb         refinement halfwidth around each D17 boundary
#                       (default 200 kb theta, 100 kb ghsl)
#   --noise_k           CUSUM peak threshold = median + k*mad (default 1.5)
#   --min_peak_sep_bp   min separation between peaks within one sample
#                       (default 15000 bp)
#   --cluster_eps_bp    1D hierarchical clustering cut for grouping
#                       per-sample peaks into nodes (default 10000 bp)
#   --edge_jaccard_min  minimum jaccard to keep an edge (default 0.30)
#   --min_node_samples  minimum samples contributing to a node (default 5)
#
# Outputs (per D17 STABLE_BLUE boundary):
#   <label>_<chr>_<boundary_idx>.nodes.tsv
#   <label>_<chr>_<boundary_idx>.edges.tsv
#   <label>_<chr>_<boundary_idx>.sample_peaks.tsv
# Plus one summary across all boundaries:
#   <label>_<chr>.refinement_summary.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# ── Args ────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

BOUNDARIES_TSV <- get_arg("--boundaries_tsv")
PRECOMP_DIR    <- get_arg("--precomp_dir")
FEATURE_SOURCE <- get_arg("--feature_source")
LABEL          <- get_arg("--label", "theta")
SCALE          <- get_arg("--scale")
OUTDIR         <- get_arg("--outdir", ".")
CHR_FILTER     <- get_arg("--chr")
PRECOMP_SUFFIX <- get_arg("--precomp_suffix",
                          if (LABEL == "ghsl") "ghsl_precomp.rds" else "precomp.rds")

# Per-label defaults that lean toward not collapsing
default_window_kb     <- if (LABEL == "ghsl") 100 else 200
default_scale         <- if (LABEL == "ghsl") "raw" else "win10000.step2000"
default_min_peak_sep  <- if (LABEL == "ghsl") 15000 else 10000   # tighter for ghsl raw 5kb
default_cluster_eps   <- if (LABEL == "ghsl") 10000 else 8000

WINDOW_KB        <- as.numeric(get_arg("--window_kb",       as.character(default_window_kb)))
NOISE_K          <- as.numeric(get_arg("--noise_k",         "1.5"))
MIN_PEAK_SEP_BP  <- as.integer(get_arg("--min_peak_sep_bp", as.character(default_min_peak_sep)))
CLUSTER_EPS_BP   <- as.integer(get_arg("--cluster_eps_bp",  as.character(default_cluster_eps)))
EDGE_JACCARD_MIN <- as.numeric(get_arg("--edge_jaccard_min","0.30"))
MIN_NODE_SAMP    <- as.integer(get_arg("--min_node_samples","5"))
if (is.na(SCALE)) SCALE <- default_scale

stopifnot(!is.na(BOUNDARIES_TSV), !is.na(PRECOMP_DIR), !is.na(FEATURE_SOURCE))
stopifnot(LABEL %in% c("theta", "ghsl"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

LP <- paste0("[", LABEL, ".E_refine] ")
cat(LP, "boundaries=", BOUNDARIES_TSV, "\n", sep = "")
cat(LP, "precomp_dir=", PRECOMP_DIR, " (suffix=", PRECOMP_SUFFIX, ")\n", sep = "")
cat(LP, "feature_source=", FEATURE_SOURCE, " (label=", LABEL,
    ", scale=", SCALE, ")\n", sep = "")
cat(LP, "window_kb=", WINDOW_KB, " noise_k=", NOISE_K,
    " min_peak_sep_bp=", MIN_PEAK_SEP_BP,
    " cluster_eps_bp=", CLUSTER_EPS_BP,
    " edge_jaccard_min=", EDGE_JACCARD_MIN,
    " min_node_samples=", MIN_NODE_SAMP, "\n", sep = "")

# ── Load D17 boundaries, restrict to STABLE_BLUE ────────────────────────────
bdt <- fread(BOUNDARIES_TSV)
if ("validation_status" %in% names(bdt)) {
  bdt <- bdt[validation_status == "STABLE_BLUE"]
}
if (!is.na(CHR_FILTER)) bdt <- bdt[chr == CHR_FILTER]
if (nrow(bdt) == 0L) {
  cat(LP, "no STABLE_BLUE boundaries to refine — exit\n", sep = ""); quit(save = "no")
}
cat(LP, "refining ", nrow(bdt), " boundaries\n", sep = "")

# ── Per-label feature loaders ────────────────────────────────────────────────
# Each loader returns:
#   list(feat_mat   = matrix [n_samp × n_win],
#        win_starts = integer [n_win] start_bp,
#        win_ends   = integer [n_win] end_bp,
#        win_mids   = integer [n_win],
#        sample_names = character [n_samp])
load_feature_chrom_theta <- function(chr, feature_dir, scale_label) {
  tsv <- file.path(feature_dir, sprintf("theta_native.%s.%s.tsv.gz", chr, scale_label))
  if (!file.exists(tsv)) stop("[loader theta] missing ", tsv)
  long_dt <- fread(tsv)
  long_dt <- long_dt[chrom == chr]
  if (nrow(long_dt) == 0L) stop("[loader theta] no rows for chrom ", chr)
  win_grid <- unique(long_dt[, .(window_idx, start_bp, end_bp)])
  setorder(win_grid, window_idx)
  win_grid[, mid_bp := as.integer((start_bp + end_bp) / 2L)]
  n_win <- nrow(win_grid)
  sample_names <- sort(unique(long_dt$sample))
  n_samp <- length(sample_names)
  samp_to_row <- setNames(seq_along(sample_names), sample_names)
  win_to_col  <- setNames(seq_len(n_win), as.character(win_grid$window_idx))
  feat_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                     dimnames = list(sample_names, NULL))
  rows <- samp_to_row[long_dt$sample]
  cols <- win_to_col[as.character(long_dt$window_idx)]
  good <- !is.na(rows) & !is.na(cols)
  feat_mat[cbind(rows[good], cols[good])] <- long_dt$theta_pi[good]
  list(feat_mat = feat_mat, win_starts = win_grid$start_bp,
       win_ends = win_grid$end_bp, win_mids = win_grid$mid_bp,
       sample_names = sample_names)
}

load_feature_chrom_ghsl <- function(chr, feature_dir, scale_label) {
  rds <- file.path(feature_dir, sprintf("%s.ghsl_matrices.rds", chr))
  if (!file.exists(rds)) stop("[loader ghsl] missing ", rds)
  gm <- readRDS(rds)
  feat_mat <- if (scale_label == "raw") {
    gm$div_mat
  } else {
    if (is.null(gm$rolling[[scale_label]]))
      stop("[loader ghsl] scale '", scale_label, "' not in ghsl_matrices.rds")
    gm$rolling[[scale_label]]
  }
  list(feat_mat   = feat_mat,
       win_starts = gm$window_info$start_bp,
       win_ends   = gm$window_info$end_bp,
       win_mids   = gm$window_info$mid_bp,
       sample_names = gm$sample_names)
}

load_feature_chrom <- function(chr) {
  if (LABEL == "theta") load_feature_chrom_theta(chr, FEATURE_SOURCE, SCALE)
  else                  load_feature_chrom_ghsl( chr, FEATURE_SOURCE, SCALE)
}

# =============================================================================
# Core algorithm
# =============================================================================

# Two-sided CUSUM over a numeric vector. Returns a data.table with one row
# per detected peak: $pos (index into x), $score (peak height), $sign (+1/-1).
# Reference is the median of the input (robust). Sample-level traces tend
# to be short (a few dozen windows) so the median is stable.
two_sided_cusum <- function(x, noise_k, min_peak_sep_idx) {
  ok <- is.finite(x)
  if (sum(ok) < 5L) return(data.table(pos = integer(0), score = numeric(0), sign = integer(0)))
  ref <- median(x[ok])
  # Centered residuals
  r <- x - ref
  r[!is.finite(r)] <- 0
  # Standard two-sided CUSUM
  s_up   <- numeric(length(x)); s_dn <- numeric(length(x))
  for (k in seq_along(x)) {
    if (k == 1L) {
      s_up[k] <- max(0, r[k]); s_dn[k] <- max(0, -r[k])
    } else {
      s_up[k] <- max(0, s_up[k - 1L] + r[k])
      s_dn[k] <- max(0, s_dn[k - 1L] - r[k])
    }
  }
  # Detect peaks in s_up and s_dn separately. A peak = local max above an
  # adaptive threshold = median + noise_k * mad over each CUSUM trace.
  detect_peaks <- function(s, sign_val) {
    if (all(s == 0)) return(NULL)
    s_pos <- s[s > 0]
    thr_med <- if (length(s_pos) >= 5L) median(s_pos) else 0
    thr_mad <- if (length(s_pos) >= 5L) mad(s_pos)    else 1
    if (!is.finite(thr_mad) || thr_mad < 1e-9) thr_mad <- 0.1 * (thr_med + 1e-9)
    threshold <- thr_med + noise_k * thr_mad
    # 1D peak finding: a peak is a local max with score above threshold and
    # min_peak_sep_idx separation to the next peak.
    n <- length(s); is_peak <- logical(n)
    for (i in seq_len(n)) {
      if (s[i] < threshold) next
      lo <- max(1L, i - 1L); hi <- min(n, i + 1L)
      if (s[i] >= s[lo] && s[i] >= s[hi]) is_peak[i] <- TRUE
    }
    peaks <- which(is_peak)
    if (length(peaks) == 0L) return(NULL)
    # Enforce min_peak_sep: keep highest-scoring peak in each conflict cluster
    keep <- rep(TRUE, length(peaks))
    if (length(peaks) > 1L) {
      for (i in seq_along(peaks)) {
        if (!keep[i]) next
        for (j in seq.int(i + 1L, length(peaks))) {
          if (j > length(peaks) || !keep[j]) next
          if (abs(peaks[j] - peaks[i]) < min_peak_sep_idx) {
            if (s[peaks[i]] >= s[peaks[j]]) keep[j] <- FALSE
            else { keep[i] <- FALSE; break }
          } else break
        }
      }
    }
    peaks <- peaks[keep]
    data.table(pos = peaks, score = s[peaks], sign = as.integer(sign_val))
  }
  up   <- detect_peaks(s_up, +1L)
  down <- detect_peaks(s_dn, -1L)
  if (is.null(up) && is.null(down)) return(data.table(pos = integer(0), score = numeric(0), sign = integer(0)))
  rbindlist(list(up, down), use.names = TRUE, fill = TRUE)
}

# Build the breakpoint graph from a per-sample peak table:
#   peaks_dt: columns (sample, pos_bp, score, sign)
# Clustering is 1D in bp space (sign is a node attribute, not a clustering
# dimension — same physical breakpoint manifests up in HET, down weakly in
# HOM_INV; they should cluster together).
build_breakpoint_graph <- function(peaks_dt, cluster_eps_bp, min_node_samp,
                                   edge_jaccard_min) {
  if (nrow(peaks_dt) == 0L) {
    return(list(nodes = data.table(), edges = data.table()))
  }
  # Hierarchical 1D clustering on bp positions. Use complete linkage so
  # nodes stay tight; cut at cluster_eps_bp.
  pos <- peaks_dt$pos_bp
  if (length(unique(pos)) == 1L) {
    cluster_id <- rep(1L, length(pos))
  } else {
    d <- dist(pos)
    hc <- hclust(d, method = "complete")
    cluster_id <- cutree(hc, h = cluster_eps_bp)
  }
  peaks_dt[, cluster_id := cluster_id]

  # Aggregate to node table
  nodes <- peaks_dt[, .(
    bp_median  = as.integer(median(pos_bp)),
    bp_ci_lo   = as.integer(quantile(pos_bp, 0.25)),
    bp_ci_hi   = as.integer(quantile(pos_bp, 0.75)),
    n_samples  = uniqueN(sample),
    n_up       = sum(sign ==  1L),
    n_down     = sum(sign == -1L),
    score_med  = median(score),
    samples_csv = paste(sort(unique(sample)), collapse = ",")
  ), by = cluster_id]
  setorder(nodes, bp_median)
  # Filter nodes with too few contributing samples
  nodes <- nodes[n_samples >= min_node_samp]
  if (nrow(nodes) == 0L) return(list(nodes = data.table(), edges = data.table()))
  nodes[, node_id := paste0("n", seq_len(.N))]
  nodes[, direction_pattern := sprintf("%d↑+%d↓", n_up, n_down)]

  # Build edges by sample-set overlap
  edges <- list()
  node_samples <- lapply(nodes$samples_csv, function(s) strsplit(s, ",")[[1]])
  for (a in seq_len(nrow(nodes) - 1L)) {
    for (b in seq.int(a + 1L, nrow(nodes))) {
      shared <- intersect(node_samples[[a]], node_samples[[b]])
      union_n <- length(union(node_samples[[a]], node_samples[[b]]))
      if (union_n == 0L) next
      jac <- length(shared) / union_n
      if (jac < edge_jaccard_min) next
      # Direction relation: how do the shared samples' signs at A vs B compare?
      # paired_flanks  = shared samples mostly up at A and down at B (or vice versa)
      # co_event       = shared samples mostly same direction
      # weak_overlap   = jaccard borderline, no clear direction pattern
      shared_at_a <- peaks_dt[sample %in% shared & cluster_id == nodes$cluster_id[a],
                              .(sign_a = sign[which.max(score)]), by = sample]
      shared_at_b <- peaks_dt[sample %in% shared & cluster_id == nodes$cluster_id[b],
                              .(sign_b = sign[which.max(score)]), by = sample]
      pair_sign <- merge(shared_at_a, shared_at_b, by = "sample")
      n_opposite <- sum(pair_sign$sign_a == -pair_sign$sign_b)
      n_same     <- sum(pair_sign$sign_a ==  pair_sign$sign_b)
      n_total    <- nrow(pair_sign)
      rel <- if (n_total < 3L) {
        "weak_overlap"
      } else if (n_opposite / n_total >= 0.66) {
        "paired_flanks"
      } else if (n_same / n_total >= 0.66) {
        "co_event"
      } else {
        "weak_overlap"
      }
      edges[[length(edges) + 1L]] <- data.table(
        node_a = nodes$node_id[a],
        node_b = nodes$node_id[b],
        distance_bp = nodes$bp_median[b] - nodes$bp_median[a],
        n_shared_samples = length(shared),
        jaccard = round(jac, 3),
        n_opposite_dir = n_opposite,
        n_same_dir = n_same,
        inferred_relation = rel
      )
    }
  }
  edges_dt <- if (length(edges) > 0L) rbindlist(edges) else data.table()
  list(nodes = nodes, edges = edges_dt)
}

# =============================================================================
# Main loop
# =============================================================================

# Group boundaries by chrom so we load each chrom's feature data once
bdt[, chr_grp := chr]
chroms <- unique(bdt$chr_grp)
summary_rows <- list()

for (cur_chr in chroms) {
  cat(LP, "loading feature data for ", cur_chr, "\n", sep = "")
  fd <- tryCatch(load_feature_chrom(cur_chr), error = function(e) {
    cat(LP, "  load failed: ", conditionMessage(e), "\n", sep = ""); NULL
  })
  if (is.null(fd)) next

  feat_mat <- fd$feat_mat
  win_starts <- fd$win_starts; win_ends <- fd$win_ends; win_mids <- fd$win_mids
  sample_names <- fd$sample_names
  n_win <- length(win_starts)
  median_win_size <- median(win_ends - win_starts + 1L, na.rm = TRUE)
  min_peak_sep_idx <- max(1L, as.integer(round(MIN_PEAK_SEP_BP / median_win_size)))
  cat(LP, "  feat_mat ", nrow(feat_mat), "x", ncol(feat_mat),
      ", median_win_size=", median_win_size,
      " bp -> min_peak_sep_idx=", min_peak_sep_idx, " windows\n", sep = "")

  chr_bdt <- bdt[chr == cur_chr]
  for (br in seq_len(nrow(chr_bdt))) {
    b <- chr_bdt[br]
    bdy_idx <- b$boundary_idx
    bdy_bp  <- as.integer(b$boundary_bp)
    win_lo_bp <- bdy_bp - WINDOW_KB * 1000L
    win_hi_bp <- bdy_bp + WINDOW_KB * 1000L

    # Restrict to refinement window
    in_win <- which(win_mids >= win_lo_bp & win_mids <= win_hi_bp)
    if (length(in_win) < 5L) {
      cat(LP, "  ", bdy_idx, ": only ", length(in_win),
          " windows in refinement range — skip\n", sep = ""); next
    }
    sub_feat <- feat_mat[, in_win, drop = FALSE]
    sub_mids <- win_mids[in_win]

    # Per-sample CUSUM
    sample_peaks <- list()
    for (s_idx in seq_len(nrow(sub_feat))) {
      x <- sub_feat[s_idx, ]
      peaks <- two_sided_cusum(x, NOISE_K, min_peak_sep_idx)
      if (nrow(peaks) > 0L) {
        peaks[, sample := sample_names[s_idx]]
        peaks[, pos_bp := sub_mids[pos]]
        sample_peaks[[length(sample_peaks) + 1L]] <- peaks
      }
    }
    if (length(sample_peaks) == 0L) {
      cat(LP, "  ", bdy_idx, ": no per-sample peaks found — skip\n", sep = ""); next
    }
    peaks_dt <- rbindlist(sample_peaks, use.names = TRUE, fill = TRUE)
    peaks_dt <- peaks_dt[, .(sample, pos_bp, score, sign)]

    # Build the breakpoint graph
    graph <- build_breakpoint_graph(peaks_dt, CLUSTER_EPS_BP,
                                    MIN_NODE_SAMP, EDGE_JACCARD_MIN)
    n_nodes <- nrow(graph$nodes); n_edges <- nrow(graph$edges)

    # Identify the primary paired-flanks edge if any (highest combined score,
    # n_shared_samples and inferred_relation == paired_flanks)
    primary_bp_lo <- NA_integer_; primary_bp_hi <- NA_integer_
    primary_n_samples <- NA_integer_
    if (n_edges > 0L) {
      pf <- graph$edges[inferred_relation == "paired_flanks"]
      if (nrow(pf) > 0L) {
        # Pick the edge with the most shared samples (most carriers)
        best <- pf[which.max(n_shared_samples)]
        node_a <- graph$nodes[node_id == best$node_a]
        node_b <- graph$nodes[node_id == best$node_b]
        primary_bp_lo <- min(node_a$bp_median, node_b$bp_median)
        primary_bp_hi <- max(node_a$bp_median, node_b$bp_median)
        primary_n_samples <- as.integer(best$n_shared_samples)
      }
    }

    # Write outputs (per-boundary)
    base <- file.path(OUTDIR, sprintf("%s_%s_%s", LABEL, cur_chr, bdy_idx))
    fwrite(graph$nodes, paste0(base, ".nodes.tsv"), sep = "\t")
    fwrite(graph$edges, paste0(base, ".edges.tsv"), sep = "\t")
    fwrite(peaks_dt,    paste0(base, ".sample_peaks.tsv"), sep = "\t")

    cat(LP, "  ", bdy_idx, " @ ", round(bdy_bp / 1e6, 2), " Mb: ",
        n_nodes, " nodes, ", n_edges, " edges",
        if (!is.na(primary_bp_lo))
          sprintf(", primary paired flanks: %.3f-%.3f Mb (n=%d)",
                  primary_bp_lo / 1e6, primary_bp_hi / 1e6, primary_n_samples)
        else "",
        "\n", sep = "")

    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      chr = cur_chr,
      boundary_idx = bdy_idx,
      w_boundary   = as.integer(b$boundary_w),
      bp_boundary  = bdy_bp,
      n_nodes      = n_nodes,
      n_edges      = n_edges,
      n_paired_flanks = if (n_edges > 0L)
                          graph$edges[inferred_relation == "paired_flanks", .N] else 0L,
      n_co_event      = if (n_edges > 0L)
                          graph$edges[inferred_relation == "co_event", .N] else 0L,
      primary_bp_lo   = primary_bp_lo,
      primary_bp_hi   = primary_bp_hi,
      primary_n_samples = primary_n_samples,
      primary_size_bp = if (!is.na(primary_bp_lo)) primary_bp_hi - primary_bp_lo else NA_integer_
    )
  }
  rm(fd, feat_mat); invisible(gc(verbose = FALSE))
}

if (length(summary_rows) > 0L) {
  summary_dt <- rbindlist(summary_rows, fill = TRUE)
  out_sum <- file.path(OUTDIR, sprintf("%s_refinement_summary.tsv", LABEL))
  fwrite(summary_dt, out_sum, sep = "\t")
  cat(LP, "summary: ", out_sum, " (", nrow(summary_dt), " boundaries)\n", sep = "")
}

cat(LP, "done\n", sep = "")
