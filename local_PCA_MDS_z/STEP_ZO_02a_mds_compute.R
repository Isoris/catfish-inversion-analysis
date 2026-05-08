#!/usr/bin/env Rscript

# =============================================================================
# 02a_mds_compute.R
#
# STAGE 1 (compute) — per-focal-chromosome MDS on lostruct distances. One
# SLURM array task PER FOCAL CHROMOSOME. Each task loads ALL per-chr
# .window_pca.rds (needed for chunked background sampling), computes the
# lostruct distance matrix on focal+background windows, runs cmdscale, and
# extracts only the focal chromosome's coordinates. Heaviest step in the
# pipeline (~12 h walltime per chromosome).
#
# Pipeline position:
#   01c local_pca_merge  -> 02a MDS_COMPUTE  -> 02b mds_merge  -> 03 precompute
#
# ── Inputs ────────────────────────────────────────────────────────────────
#   --rds_dir    <dir>   directory of <chr>.window_pca.rds from step 01c
#   --outdir     <dir>   output root; this stage writes to <outdir>/tmp/
#   --outprefix  <s>     filename stem (default 'inversion_localpca')
#   --focal_chr  <name>  chromosome under analysis for this array task
#   [--mds_mode  chunked_2x]
#                        chromosome | global | chunked_2x | chunked_3x | chunked_4x
#                        chunked_Nx samples N x focal-window-count of background
#                        from non-focal chroms, EXCLUDING high-inv_likeness windows
#                        (so other chroms' inversions don't leak into background).
#   [--npc       4]      eigvec count per window (must match 01b!)
#   [--mds_dims  20]     MDS dimensions
#   [--z_thresh  3.0]    candidate z threshold (consumed downstream)
#
# ── Output (in <outdir>/tmp/) ─────────────────────────────────────────────
#   <focal_chr>.mds_perchr.rds   focal-chr MDS coords + per-axis z scores +
#                                background metadata + chunked-mode bg indices
#
# ── Codebase ──────────────────────────────────────────────────────────────
#   inversion-popgen-toolkit v8.5 / consolidated layout v1.0
#   (was: STEP10v2_stage1_perchr_mds.R v8.3-parallel)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# PARSE ARGS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

rds_dir    <- NULL
outdir     <- NULL
outprefix  <- "inversion_localpca"
FOCAL_CHR  <- NULL
MDS_MODE   <- "chunked_2x"
NPC        <- 4L           # v9.4: default raised from 2L → 4L to match
                           # 2a_local_pca/STEP_A02|A03 default. The lostruct
                           # distance uses top-NPC eigvals/vecs; running B01
                           # at NPC=2 against A02/A03 windows computed at
                           # NPC=4 silently used only half the available
                           # signal and produced different distances than
                           # the legacy STEP_B01 monolithic (which was 4).
MDS_DIMS   <- 20L
Z_THRESH   <- 3.0
SEED       <- 42L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--rds_dir" && i < length(args)) {
    rds_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--outdir" && i < length(args)) {
    outdir <- args[i + 1]; i <- i + 2L
  } else if (a == "--outprefix" && i < length(args)) {
    outprefix <- args[i + 1]; i <- i + 2L
  } else if (a == "--focal_chr" && i < length(args)) {
    FOCAL_CHR <- args[i + 1]; i <- i + 2L
  } else if (a == "--mds_mode" && i < length(args)) {
    MDS_MODE <- args[i + 1]; i <- i + 2L
  } else if (a == "--npc" && i < length(args)) {
    NPC <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--mds_dims" && i < length(args)) {
    MDS_DIMS <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--z_thresh" && i < length(args)) {
    Z_THRESH <- as.numeric(args[i + 1]); i <- i + 2L
  } else if (a == "--seed" && i < length(args)) {
    SEED <- as.integer(args[i + 1]); i <- i + 2L
  } else {
    i <- i + 1L
  }
}

if (is.null(rds_dir) || is.null(outdir) || is.null(FOCAL_CHR)) {
  stop("Usage: Rscript STEP10v2_stage1_perchr_mds.R --rds_dir <dir> --outdir <dir> --focal_chr <chr> ...")
}

valid_modes <- c("chromosome", "global", "chunked_2x", "chunked_3x", "chunked_4x")
if (!(MDS_MODE %in% valid_modes)) {
  stop("Invalid --mds_mode. Must be one of: ", paste(valid_modes, collapse = ", "))
}

CHUNK_K <- switch(MDS_MODE,
  chunked_2x = 2L,
  chunked_3x = 3L,
  chunked_4x = 4L,
  0L
)

tmpdir <- file.path(outdir, "tmp")
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SKIP-IF-DONE
# =============================================================================

rds_out <- file.path(tmpdir, paste0(FOCAL_CHR, ".mds_perchr.rds"))

if (file.exists(rds_out)) {
  message("[STEP10v2-S1] ", FOCAL_CHR, ": output already exists, skipping")
  message("  ", rds_out)
  quit(status = 0)
}

message("[STEP10v2-S1] ═══════ ", FOCAL_CHR, " (mode=", MDS_MODE, ") ═══════")
message("[STEP10v2-S1] nPC=", NPC, " mds_dims=", MDS_DIMS, " z_thresh=", Z_THRESH)
if (CHUNK_K > 0) message("[STEP10v2-S1] Chunk multiplier: ", CHUNK_K, "x, seed=", SEED)

# =============================================================================
# PROGRESS BAR
# =============================================================================

format_eta <- function(eta_sec) {
  if (length(eta_sec) == 0 || is.na(eta_sec) || !is.finite(eta_sec) || eta_sec < 0) {
    return("NA")
  }
  if (eta_sec >= 3600) {
    return(sprintf("%dh%02dm",
                   as.integer(floor(eta_sec / 3600)),
                   as.integer(floor((eta_sec %% 3600) / 60))))
  }
  if (eta_sec >= 60) {
    return(sprintf("%dm%02ds",
                   as.integer(floor(eta_sec / 60)),
                   as.integer(floor(eta_sec %% 60))))
  }
  sprintf("%ds", as.integer(round(eta_sec)))
}

make_progress <- function(total, label) {
  total <- as.double(total)
  t0 <- proc.time()[3]
  last_print <- -Inf

  function(current) {
    current <- as.double(current)
    now <- proc.time()[3]

    if ((now - last_print) < 2 && current < total) return(invisible())
    last_print <<- now

    elapsed <- as.double(now - t0)
    pct <- if (total > 0) 100 * current / total else NA_real_
    eta_sec <- if (current > 0 && is.finite(current) && is.finite(total)) {
      elapsed * (total - current) / current
    } else {
      NA_real_
    }
    eta_str <- format_eta(eta_sec)

    bar_w <- 25L
    frac <- if (is.finite(pct)) max(0, min(1, pct / 100)) else 0
    filled <- as.integer(floor(bar_w * frac))
    bar <- paste0("[", strrep("=", filled), strrep(" ", bar_w - filled), "]")

    message(sprintf("[dist] %s  %s  %.0f/%.0f  (%.1f%%)  ETA %s",
                    FOCAL_CHR, bar, current, total,
                    ifelse(is.finite(pct), pct, 0), eta_str))
  }
}

# =============================================================================
# LOSTRUCT DISTANCE
# =============================================================================

dist_sq_from_pcs <- function(values1, vectors1, values2, vectors2) {
  Xt <- crossprod(vectors2, vectors1)
  cross <- sum(outer(values2, values1) * Xt^2)
  sum(values1^2) + sum(values2^2) - 2 * cross
}

pc_dist_from_step09 <- function(pca_df, sample_names, npc, normalize = "L1",
                                progress_fn = NULL) {
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

  total_pairs <- as.double(n) * (as.double(n) - 1) / 2
  pair_count <- 0

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

        pair_count <- pair_count + 1
        if (!is.null(progress_fn)) progress_fn(pair_count)
      }
    }
  }

  out
}

# =============================================================================
# LOAD ALL STEP09b DATA
# =============================================================================

message("[STEP10v2-S1] Loading per-chr RDS files...")
t_load <- proc.time()[3]

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
    meta  = as.data.table(obj$window_meta),
    pca   = as.data.table(obj$pca),
    chrom = chr
  )
}

for (chr in names(chr_data)) {
  if ("window_id" %in% names(chr_data[[chr]]$meta)) {
    chr_data[[chr]]$meta[, global_window_id := window_id]
    chr_data[[chr]]$pca[,  global_window_id := window_id]
  } else {
    stop("Missing window_id in ", chr, " — run STEP09b Stage 2 first")
  }
}

message("[STEP10v2-S1] Loaded ", length(chr_data), " chromosomes in ",
        round(proc.time()[3] - t_load, 1), "s")

if (!(FOCAL_CHR %in% names(chr_data))) {
  stop("Focal chromosome ", FOCAL_CHR, " not found in RDS files")
}

# =============================================================================
# PROCESS FOCAL CHROMOSOME
# =============================================================================

t0 <- proc.time()[3]

cd <- chr_data[[FOCAL_CHR]]
n_focal <- nrow(cd$meta)

if (n_focal < 3) {
  message("[STEP10v2-S1] ", FOCAL_CHR, ": only ", n_focal, " windows, skipping")
  saveRDS(list(skip = TRUE, chrom = FOCAL_CHR, reason = "too_few_windows"), rds_out)
  quit(status = 0)
}

dt_focal <- merge(cd$meta, cd$pca, by = "global_window_id")
background_ids <- integer(0)

if (MDS_MODE == "chromosome") {
  dt_combined <- dt_focal

} else if (MDS_MODE == "global") {
  all_dts <- lapply(chr_data, function(x) merge(x$meta, x$pca, by = "global_window_id"))
  dt_combined <- rbindlist(all_dts, fill = TRUE)

} else if (grepl("^chunked_", MDS_MODE)) {
  other_chrs <- setdiff(names(chr_data), FOCAL_CHR)

  if (length(other_chrs) == 0) {
    dt_combined <- dt_focal
  } else {
    n_background_target <- CHUNK_K * n_focal
    bg_list <- lapply(other_chrs, function(c) {
      merge(chr_data[[c]]$meta, chr_data[[c]]$pca, by = "global_window_id")
    })
    bg_pool <- rbindlist(bg_list, fill = TRUE)

    # ── SMART BACKGROUND FILTER ──────────────────────────────────────
    # Exclude windows with high inversion-likeness from background pool.
    # Without this, chunked_2x can accidentally sample other chromosomes'
    # inversions as "background," making the focal inversions look normal.
    # Uses PVE1 (proportion variance explained by PC1) as a fast proxy:
    # windows where PC1 dominates (pve1 > 0.5) are likely inversion-like.
    if ("lam_1" %in% names(bg_pool) && "lam_2" %in% names(bg_pool)) {
      bg_pool[, bg_pve1 := lam_1 / (lam_1 + lam_2)]
      n_before <- nrow(bg_pool)
      bg_pool <- bg_pool[is.na(bg_pve1) | bg_pve1 < 0.50]
      n_removed <- n_before - nrow(bg_pool)
      if (n_removed > 0) {
        message("[STEP10v2-S1] Background filter: removed ", n_removed,
                " inv-like windows (pve1 >= 0.50) from ", n_before, " → ", nrow(bg_pool))
      }
      bg_pool[, bg_pve1 := NULL]
    }

    chr_seed <- SEED + match(FOCAL_CHR, names(chr_data))
    set.seed(chr_seed)

    n_sample <- min(n_background_target, nrow(bg_pool))
    bg_idx <- sort(sample.int(nrow(bg_pool), n_sample))
    bg_sampled <- bg_pool[bg_idx]
    background_ids <- bg_sampled$global_window_id

    dt_combined <- rbindlist(list(dt_focal, bg_sampled), fill = TRUE)
  }
}

message("[STEP10v2-S1] ", FOCAL_CHR, ": ", n_focal, " focal + ",
        length(background_ids), " background = ", nrow(dt_combined), " total windows")

# =============================================================================
# DISTANCE MATRIX
# =============================================================================

n_total <- nrow(dt_combined)
n_pairs <- as.double(n_total) * (as.double(n_total) - 1) / 2

message("[STEP10v2-S1] Computing distance matrix: ", n_total, " windows (",
        format(round(n_pairs), big.mark = ","), " pairs)")

pb_dist <- make_progress(n_pairs, "dist")

dmat <- pc_dist_from_step09(
  as.data.frame(dt_combined),
  sample_names_ref,
  npc = NPC,
  normalize = "L1",
  progress_fn = pb_dist
)
message("")

keep <- which(apply(dmat, 1, function(x) all(is.finite(x))))
if (length(keep) < 3) {
  message("[STEP10v2-S1] ", FOCAL_CHR, ": <3 finite rows after filtering, skipping")
  saveRDS(list(skip = TRUE, chrom = FOCAL_CHR, reason = "insufficient_finite"), rds_out)
  quit(status = 0)
}

dt_keep <- dt_combined[keep]
dmat_keep <- dmat[keep, keep, drop = FALSE]

# =============================================================================
# MDS
# =============================================================================

message("[STEP10v2-S1] Running cmdscale (k=", min(MDS_DIMS, nrow(dmat_keep) - 1L), ")...")
t_mds <- proc.time()[3]

k_mds <- min(MDS_DIMS, nrow(dmat_keep) - 1L)
mds <- tryCatch(
  cmdscale(as.dist(dmat_keep), k = k_mds, eig = TRUE),
  error = function(e) {
    message("[WARN] MDS failed: ", e$message)
    NULL
  }
)

if (is.null(mds)) {
  saveRDS(list(skip = TRUE, chrom = FOCAL_CHR, reason = "mds_failed"), rds_out)
  quit(status = 0)
}

message("[STEP10v2-S1] MDS done in ", round(proc.time()[3] - t_mds, 1), "s")

mds_dt <- as.data.table(mds$points)
setnames(mds_dt, paste0("MDS", seq_len(ncol(mds_dt))))
mds_dt[, global_window_id := dt_keep$global_window_id]

# =============================================================================
# FOCAL-ONLY EXTRACTION
# =============================================================================

focal_mask <- dt_keep$chrom == FOCAL_CHR
focal_gwids <- dt_keep$global_window_id[focal_mask]

mds_focal <- mds_dt[global_window_id %in% focal_gwids]
dt_focal_out <- dt_keep[chrom == FOCAL_CHR]
dmat_focal_idx <- which(focal_mask)

if (length(dmat_focal_idx) < 3) {
  saveRDS(list(skip = TRUE, chrom = FOCAL_CHR, reason = "too_few_focal_after_filter"), rds_out)
  quit(status = 0)
}

dmat_focal <- dmat_keep[dmat_focal_idx, dmat_focal_idx, drop = FALSE]

# =============================================================================
# Z-SCORES ON FOCAL WINDOWS ONLY
# =============================================================================

for (ax in seq_len(ncol(mds$points))) {
  coln <- paste0("MDS", ax)
  if (!(coln %in% names(mds_focal))) next

  vv <- mds_focal[[coln]]
  zcol <- paste0("MDS", ax, "_z")

  # ROBUST z-score: median / MAD instead of mean / SD.
  # With chunked_2x, if an inversion occupies >30% of the chromosome,
  # mean/SD normalization pulls toward the inversion signal and compresses
  # z for everything. median/MAD is resistant to this because the inversion
  # windows are outliers that don't move the median.
  med <- median(vv, na.rm = TRUE)
  mad_val <- mad(vv, na.rm = TRUE)

  if (is.na(mad_val) || !is.finite(mad_val) || mad_val < 1e-10) {
    # Fallback to SD if MAD is degenerate (e.g. >50% identical values)
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

setkey(mds_focal, global_window_id)
setkey(dt_focal_out, global_window_id)
out_chr <- merge(dt_focal_out, mds_focal, by = "global_window_id")

mds_cols <- grep("^MDS\\d+$", names(mds_focal), value = TRUE)
mds_mat_focal <- as.matrix(mds_focal[, ..mds_cols])

elapsed <- round(proc.time()[3] - t0, 1)
n_outlier <- if ("MDS1_outlier" %in% names(out_chr)) sum(out_chr$MDS1_outlier, na.rm = TRUE) else 0L

# =============================================================================
# WRITE PER-CHR RESULT
# =============================================================================

result <- list(
  out_dt         = out_chr,
  dmat           = dmat_focal,
  mds            = list(points = mds_mat_focal, eig = mds$eig),
  background_ids = background_ids,
  n_background   = length(background_ids),
  n_focal        = sum(focal_mask)
)
saveRDS(result, rds_out)

meta_row <- data.table(
  run_id               = paste0(MDS_MODE, "_", FOCAL_CHR),
  mds_mode             = MDS_MODE,
  focal_chrom          = FOCAL_CHR,
  n_focal_windows      = result$n_focal,
  n_background_windows = result$n_background,
  chunk_multiplier     = CHUNK_K,
  random_seed          = SEED,
  mds_dims             = MDS_DIMS,
  z_thresh             = Z_THRESH,
  elapsed_sec          = elapsed,
  timestamp            = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
fwrite(meta_row, file.path(tmpdir, paste0(FOCAL_CHR, ".metadata.tsv")), sep = "\t")

if (CHUNK_K > 0 && length(background_ids) > 0) {
  writeLines(as.character(background_ids),
             file.path(tmpdir, paste0(FOCAL_CHR, ".background_ids.txt")))
}

message("")
message("[DONE] STEP10v2 Stage 1 — ", FOCAL_CHR, ": ", nrow(out_chr),
        " focal windows, ", n_outlier, " outliers, ", elapsed, "s")
message("  ", rds_out)
