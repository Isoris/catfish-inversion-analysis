#!/usr/bin/env Rscript

# =============================================================================
# STEP10_v2_mds_multimode.R  (v8.0)
#
# CONFIGURABLE MDS BACKEND with 5 modes for inversion discovery.
#
# PROBLEM: Pure per-chromosome MDS distorts if one regime occupies >50%
# of a chromosome's windows. Full global MDS is slow.
#
# SOLUTION: Configurable modes:
#   chromosome  — fastest, per-chr only (can distort)
#   global      — slowest, all windows together
#   chunked_2x  — focal chr + 2N sampled background windows
#   chunked_3x  — focal chr + 3N sampled background
#   chunked_4x  — focal chr + 4N sampled background
#
# CRITICAL RULE: Snake harvesting remains chromosome-local regardless of
# MDS mode. Only focal-chromosome MDS coordinates are emitted for downstream.
#
# MDS is a regime-geometry detector, NOT a normal-vs-inversion classifier.
# Do not encode biological labels from MDS position.
#
# Seeding (outlier detection) uses MDS z-scores.
# Growth (snake continuity) uses MDS coordinates + similarity matrices.
# These are distinct operations — seeds ≠ continuity.
#
# INPUTS:
#   <step09_rds_dir>  — directory with STEP09(b) .window_pca.rds files
#   <outprefix>       — output path prefix
#   --mds_mode        — chromosome|global|chunked_2x|chunked_3x|chunked_4x
#   [--npc 2] [--mds_dims 20] [--z_thresh 3]
#   [--gap_bp 500000] [--min_windows 3] [--seed 42]
#
# OUTPUTS:
#   <outprefix>.mds.rds                  — per_chr list (Snake 1 compatible)
#   <outprefix>.window_mds.tsv.gz        — all windows with MDS coordinates
#   <outprefix>.candidate_regions.tsv.gz — clustered outlier regions
#   <outprefix>.candidate_window_membership.tsv.gz — per-window candidate membership
#   <outprefix>.mds_mode_metadata.tsv    — mode, seed, K, per-chr details
#   <outprefix>.mds_background_<chr>.txt — sampled background IDs (chunked)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript STEP10_v2_mds_multimode.R <step09_rds_dir> <outprefix> ",
       "--mds_mode chromosome [--npc 2] [--mds_dims 20] ...")
}

rds_dir   <- args[1]
outprefix <- args[2]

# Defaults
MDS_MODE    <- "chromosome"
NPC         <- 4L  # v9.4: 2L → 4L to match 2a/A02/A03
MDS_DIMS    <- 20L
Z_THRESH    <- 3.0
GAP_BP      <- 500000
MIN_WINDOWS <- 3L
SEED        <- 42L

# Parse named args
i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--mds_mode" && i < length(args)) {
    MDS_MODE <- args[i + 1]; i <- i + 2L
  } else if (a == "--npc" && i < length(args)) {
    NPC <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--mds_dims" && i < length(args)) {
    MDS_DIMS <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--z_thresh" && i < length(args)) {
    Z_THRESH <- as.numeric(args[i + 1]); i <- i + 2L
  } else if (a == "--gap_bp" && i < length(args)) {
    GAP_BP <- as.numeric(args[i + 1]); i <- i + 2L
  } else if (a == "--min_windows" && i < length(args)) {
    MIN_WINDOWS <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--seed" && i < length(args)) {
    SEED <- as.integer(args[i + 1]); i <- i + 2L
  } else {
    i <- i + 1L
  }
}

valid_modes <- c("chromosome", "global", "chunked_2x", "chunked_3x", "chunked_4x")
if (!(MDS_MODE %in% valid_modes)) {
  stop("Invalid --mds_mode. Must be one of: ", paste(valid_modes, collapse = ", "))
}

CHUNK_K <- switch(MDS_MODE,
  chunked_2x = 2L, chunked_3x = 3L, chunked_4x = 4L, 0L
)

message("[STEP10v2] Multi-mode MDS (", MDS_MODE, ")")
message("[STEP10v2] nPC=", NPC, " mds_dims=", MDS_DIMS, " z_thresh=", Z_THRESH)
if (CHUNK_K > 0) message("[STEP10v2] Chunk multiplier: ", CHUNK_K, "x, seed=", SEED)

# =============================================================================
# LOSTRUCT DISTANCE (from v7.4 — unchanged)
# =============================================================================

dist_sq_from_pcs <- function(values1, vectors1, values2, vectors2) {
  Xt <- crossprod(vectors2, vectors1)
  cross <- sum(outer(values2, values1) * Xt^2)
  sum(values1^2) + sum(values2^2) - 2 * cross
}

pc_dist_from_step09 <- function(pca_df, sample_names, npc, normalize = "L1") {
  values <- as.matrix(pca_df[, paste0("lam_", seq_len(npc)), drop = FALSE])
  vec_cols <- unlist(lapply(seq_len(npc), function(pc) paste0("PC_", pc, "_", sample_names)))
  vectors <- as.matrix(pca_df[, vec_cols, drop = FALSE])
  n_samples <- length(sample_names)

  if (normalize == "L1") {
    rs <- rowSums(abs(values))
    rs[rs == 0] <- NA
    values <- values / rs
  }

  n <- nrow(values)
  out <- matrix(NA_real_, n, n)

  emat <- function(u) matrix(u, nrow = n_samples, ncol = npc)

  for (i in seq_len(n)) {
    vi <- values[i, ]
    Ui <- emat(vectors[i, ])
    out[i, i] <- 0
    if (i < n) {
      for (j in (i + 1L):n) {
        vj <- values[j, ]
        Uj <- emat(vectors[j, ])
        if (anyNA(vi) || anyNA(vj) || anyNA(Ui) || anyNA(Uj)) {
          d2 <- NA_real_
        } else {
          d2 <- dist_sq_from_pcs(vi, Ui, vj, Uj)
          if (!is.finite(d2) || d2 < 0) d2 <- 0
        }
        out[i, j] <- sqrt(d2)
        out[j, i] <- out[i, j]
      }
    }
  }
  out
}

# =============================================================================
# LOAD ALL STEP09 DATA
# =============================================================================

rds_files <- sort(list.files(rds_dir, pattern = "\\.window_pca\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .window_pca.rds files in: ", rds_dir)

chr_data <- list()
sample_names_ref <- NULL

for (f in rds_files) {
  obj <- readRDS(f)
  if (is.null(sample_names_ref)) {
    sample_names_ref <- obj$sample_names
  } else if (!identical(sample_names_ref, obj$sample_names)) {
    stop("Sample names differ across STEP09 files")
  }
  chr <- obj$chrom
  chr_data[[chr]] <- list(
    meta = as.data.table(obj$window_meta),
    pca  = as.data.table(obj$pca),
    chrom = chr
  )
}

# Assign global IDs — reuse STEP09b window_id if present, else create new
global_id <- 0L
for (chr in names(chr_data)) {
  n <- nrow(chr_data[[chr]]$meta)
  if ("window_id" %in% names(chr_data[[chr]]$meta)) {
    # STEP09b already assigned globally unique window_id — reuse it
    chr_data[[chr]]$meta[, global_window_id := window_id]
    chr_data[[chr]]$pca[, global_window_id := window_id]
  } else {
    # Legacy STEP09 (no window_id) — assign sequentially
    chr_data[[chr]]$meta[, global_window_id := global_id + seq_len(.N)]
    chr_data[[chr]]$pca[, global_window_id := global_id + seq_len(.N)]
  }
  global_id <- global_id + n
}

total_windows <- sum(vapply(chr_data, function(x) nrow(x$meta), integer(1)))
message("[STEP10v2] Loaded ", length(chr_data), " chromosomes, ", total_windows, " windows total")

# =============================================================================
# CLUSTER OUTLIERS (reused from v7.4)
# =============================================================================

cluster_outliers_bp <- function(dt, flag_col, gap_bp, min_windows) {
  idx <- which(dt[[flag_col]])
  if (length(idx) == 0) return(NULL)
  clusters <- list(); cur <- idx[1]
  if (length(idx) > 1) {
    for (ii in idx[-1]) {
      if ((dt$start_bp[ii] - dt$end_bp[max(cur)]) <= gap_bp) {
        cur <- c(cur, ii)
      } else {
        if (length(cur) >= min_windows) clusters[[length(clusters) + 1]] <- cur
        cur <- ii
      }
    }
  }
  if (length(cur) >= min_windows) clusters[[length(clusters) + 1]] <- cur
  clusters
}

# =============================================================================
# MDS CORE: process one focal chromosome
# =============================================================================

process_focal_chr <- function(focal_chr, chr_data, mode, chunk_k, seed) {
  cd <- chr_data[[focal_chr]]
  n_focal <- nrow(cd$meta)
  if (n_focal < 3) return(NULL)

  dt_focal <- merge(cd$meta, cd$pca, by = "global_window_id")
  background_ids <- integer(0)

  if (mode == "chromosome") {
    # ── CHROMOSOME MODE: focal only ──────────────────────────────────
    dt_combined <- dt_focal

  } else if (mode == "global") {
    # ── GLOBAL MODE: all windows (slow) ──────────────────────────────
    all_dts <- lapply(chr_data, function(x) merge(x$meta, x$pca, by = "global_window_id"))
    dt_combined <- rbindlist(all_dts, fill = TRUE)

  } else if (grepl("^chunked_", mode)) {
    # ── CHUNKED MODE: focal + K*N sampled background ─────────────────
    other_chrs <- setdiff(names(chr_data), focal_chr)
    if (length(other_chrs) == 0) {
      dt_combined <- dt_focal
    } else {
      n_background_target <- chunk_k * n_focal

      # Pool all non-focal windows
      bg_list <- lapply(other_chrs, function(c) {
        merge(chr_data[[c]]$meta, chr_data[[c]]$pca, by = "global_window_id")
      })
      bg_pool <- rbindlist(bg_list, fill = TRUE)

      # Sample with fixed per-chromosome seed (reproducible + different per chr)
      chr_seed <- seed + match(focal_chr, names(chr_data))
      set.seed(chr_seed)
      n_sample <- min(n_background_target, nrow(bg_pool))
      bg_idx <- sort(sample.int(nrow(bg_pool), n_sample))
      bg_sampled <- bg_pool[bg_idx]
      background_ids <- bg_sampled$global_window_id

      dt_combined <- rbindlist(list(dt_focal, bg_sampled), fill = TRUE)
    }
  }

  # ── Compute distance matrix ────────────────────────────────────────
  dmat <- pc_dist_from_step09(as.data.frame(dt_combined), sample_names_ref,
                               npc = NPC, normalize = "L1")

  # Keep rows with all finite distances
  keep <- which(apply(dmat, 1, function(x) all(is.finite(x))))
  if (length(keep) < 3) return(NULL)

  dt_keep <- dt_combined[keep]
  dmat_keep <- dmat[keep, keep, drop = FALSE]

  # ── MDS ────────────────────────────────────────────────────────────
  k_mds <- min(MDS_DIMS, nrow(dmat_keep) - 1L)
  mds <- tryCatch(
    cmdscale(as.dist(dmat_keep), k = k_mds, eig = TRUE),
    error = function(e) { message("[WARN] MDS failed for ", focal_chr); NULL }
  )
  if (is.null(mds)) return(NULL)

  mds_dt <- as.data.table(mds$points)
  setnames(mds_dt, paste0("MDS", seq_len(ncol(mds_dt))))
  mds_dt[, global_window_id := dt_keep$global_window_id]

  # ── CRITICAL: Extract only focal chromosome windows ────────────────
  # Snake harvesting MUST remain chromosome-local
  focal_mask <- dt_keep$chrom == focal_chr
  focal_gwids <- dt_keep$global_window_id[focal_mask]

  mds_focal <- mds_dt[global_window_id %in% focal_gwids]
  dt_focal_out <- dt_keep[chrom == focal_chr]
  dmat_focal_idx <- which(focal_mask)

  if (length(dmat_focal_idx) < 3) return(NULL)
  dmat_focal <- dmat_keep[dmat_focal_idx, dmat_focal_idx, drop = FALSE]

  # ── Z-scores on FOCAL windows only (per-chromosome normalization) ──
  # ROBUST: median/MAD instead of mean/SD to resist inversion-heavy chromosomes
  for (ax in seq_len(ncol(mds$points))) {
    coln <- paste0("MDS", ax)
    if (!(coln %in% names(mds_focal))) next
    vv <- mds_focal[[coln]]
    zcol <- paste0("MDS", ax, "_z")
    med <- median(vv, na.rm = TRUE)
    mad_val <- mad(vv, na.rm = TRUE)
    if (is.na(mad_val) || !is.finite(mad_val) || mad_val < 1e-10) {
      sdev <- sd(vv, na.rm = TRUE)
      if (is.na(sdev) || !is.finite(sdev) || sdev == 0) {
        mds_focal[[zcol]] <- 0
      } else {
        mds_focal[[zcol]] <- (vv - mean(vv, na.rm = TRUE)) / sdev
      }
    } else {
      mds_focal[[zcol]] <- (vv - med) / mad_val
    }
    mds_focal[[paste0("MDS", ax, "_outlier")]] <- abs(mds_focal[[zcol]]) >= Z_THRESH
  }

  # Merge
  setkey(mds_focal, global_window_id)
  setkey(dt_focal_out, global_window_id)
  out_chr <- merge(dt_focal_out, mds_focal, by = "global_window_id")

  # Extract MDS matrix for focal windows (for Snake 1)
  mds_mat_focal <- as.matrix(mds_focal[, grep("^MDS\\d+$", names(mds_focal), value = TRUE),
                                         with = FALSE])

  list(
    out_dt = out_chr,
    dmat = dmat_focal,
    mds = list(points = mds_mat_focal, eig = mds$eig),
    background_ids = background_ids,
    n_background = length(background_ids),
    n_focal = sum(focal_mask)
  )
}

# =============================================================================
# MAIN: PROCESS ALL CHROMOSOMES
# =============================================================================

all_chr_results <- list()
metadata_rows <- list()

for (chr in names(chr_data)) {
  message("\n[STEP10v2] ═══════ ", chr, " (", nrow(chr_data[[chr]]$meta), " windows, mode=",
          MDS_MODE, ") ═══════")

  t0 <- proc.time()
  result <- tryCatch(
    process_focal_chr(chr, chr_data, MDS_MODE, CHUNK_K, SEED),
    error = function(e) { message("[WARN] ", chr, ": ", e$message); NULL }
  )
  elapsed <- (proc.time() - t0)[3]

  if (is.null(result)) {
    message("[SKIP] ", chr)
    next
  }

  all_chr_results[[chr]] <- result[c("out_dt", "dmat", "mds")]

  # Save background IDs for chunked modes
  if (CHUNK_K > 0 && length(result$background_ids) > 0) {
    bg_file <- paste0(outprefix, ".mds_background_", chr, ".txt")
    writeLines(as.character(result$background_ids), bg_file)
  }

  # Metadata
  metadata_rows[[chr]] <- data.table(
    run_id           = paste0(MDS_MODE, "_", chr),
    mds_mode         = MDS_MODE,
    focal_chrom      = chr,
    n_focal_windows  = result$n_focal,
    n_background_windows = result$n_background,
    chunk_multiplier = CHUNK_K,
    random_seed      = SEED,
    mds_dims         = MDS_DIMS,
    z_thresh         = Z_THRESH,
    elapsed_sec      = round(elapsed, 1),
    timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  n_outlier <- sum(result$out_dt$MDS1_outlier, na.rm = TRUE)
  message("[STEP10v2] ", chr, ": ", nrow(result$out_dt), " focal windows, ",
          n_outlier, " outliers, ", round(elapsed, 1), "s")
}

if (length(all_chr_results) == 0) stop("No chromosomes produced results")

# =============================================================================
# MERGE + CANDIDATE REGIONS
# =============================================================================

out_dt <- rbindlist(lapply(all_chr_results, function(x) x$out_dt), fill = TRUE)
message("\n[STEP10v2] Merged: ", nrow(out_dt), " windows across ",
        length(all_chr_results), " chromosomes")

# Cluster outliers into candidates
mds_cols <- grep("^MDS\\d+$", names(out_dt), value = TRUE)
cand_list <- list()
membership_list <- list()
cand_id <- 0L

for (ax in seq_along(mds_cols)) {
  flag_col <- paste0("MDS", ax, "_outlier")
  if (!(flag_col %in% names(out_dt))) next

  for (chr in unique(out_dt$chrom)) {
    sub <- out_dt[chrom == chr][order(start_bp)]
    clusters <- cluster_outliers_bp(sub, flag_col, GAP_BP, MIN_WINDOWS)
    if (!is.null(clusters)) {
      for (cl in clusters) {
        cand_id <- cand_id + 1L
        xx <- sub[cl]
        cand_list[[length(cand_list) + 1]] <- data.table(
          candidate_id = cand_id, mds_axis = ax, chrom = chr,
          start_bp = min(xx$start_bp), end_bp = max(xx$end_bp),
          center_bp = as.integer((min(xx$start_bp) + max(xx$end_bp)) / 2),
          n_windows = nrow(xx),
          first_global_window_id = min(xx$global_window_id),
          last_global_window_id = max(xx$global_window_id)
        )
        # Membership table: one row per window in this candidate
        for (k in seq_len(nrow(xx))) {
          membership_list[[length(membership_list) + 1]] <- data.table(
            candidate_id = cand_id,
            global_window_id = xx$global_window_id[k],
            chrom = chr,
            start_bp = xx$start_bp[k],
            end_bp = xx$end_bp[k],
            mds_axis = ax
          )
        }
      }
    }
  }
}

cand_dt <- if (length(cand_list) > 0) rbindlist(cand_list) else {
  data.table(candidate_id = integer(), chrom = character())
}

membership_dt <- if (length(membership_list) > 0) rbindlist(membership_list) else {
  data.table(candidate_id = integer(), global_window_id = integer(), chrom = character())
}

message("[STEP10v2] Candidate regions: ", nrow(cand_dt),
        " (", nrow(membership_dt), " member windows)")

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

meta_dt <- if (length(metadata_rows) > 0) rbindlist(metadata_rows) else data.table()

f1 <- paste0(outprefix, ".window_mds.tsv.gz")
f2 <- paste0(outprefix, ".candidate_regions.tsv.gz")
f3 <- paste0(outprefix, ".mds.rds")
f4 <- paste0(outprefix, ".mds_mode_metadata.tsv")
f5 <- paste0(outprefix, ".candidate_window_membership.tsv.gz")

fwrite(out_dt, f1, sep = "\t")
fwrite(cand_dt, f2, sep = "\t")
saveRDS(list(dt = out_dt, candidate_regions = cand_dt,
             per_chr = all_chr_results, mds_mode = MDS_MODE), f3)
fwrite(meta_dt, f4, sep = "\t")
fwrite(membership_dt, f5, sep = "\t")

message("\n[DONE] STEP10 v2 multi-mode MDS complete (mode=", MDS_MODE, ")")
message("  ", f1)
message("  ", f2)
message("  ", f5)
message("  ", f3)
message("  ", f4)
