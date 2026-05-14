#!/usr/bin/env Rscript

# =============================================================================
# STEP_ZO_E_mds_compute.R
#
# Per-focal-chromosome MDS on lostruct distances. One SLURM array task per
# chromosome (~4 min/chrom on Lanta under chromosome mode; the chunked_2x
# mode that produced ~12 h walltimes has been removed). Each task:
#   1. loads its focal chrom's per-window PCA fingerprints,
#   2. computes the lostruct distance matrix `dmat` (focal x focal),
#   3. runs `cmdscale(k = MDS_DIMS)` on dmat,
#   4. computes robust per-axis z scores on the focal MDS coords,
#   5. writes the per-chrom result to <outdir>/mds_perchr/<chr>.mds_perchr.rds.
#
# Pipeline position:
#   ZO_C/D (per-window PCA + global IDs)  ->  ZO_E (this script)
#       ->  ZO_G (precomp + sim_mats, reads mds_perchr/ directly)
#       ->  ZO_H / ZO_J (L1 / L2 stripe detection on sim_mat)
#
# Inputs:
#   --rds_dir    <dir>   directory of <chr>.window_pca.rds from ZO_C/D
#   --outdir     <dir>   output root; this stage writes to <outdir>/mds_perchr/
#   --outprefix  <s>     filename stem (default 'inversion_localpca')
#   --focal_chr  <name>  chromosome under analysis for this array task
#   [--npc        4]     eigvec count per window (must match ZO_C!)
#   [--mds_dims   5]     MDS dimensions kept after cmdscale (see comment)
#   [--z_thresh   3.0]   per-axis Z #1 threshold (descriptive only)
#   [--mds_mode / --seed] legacy chunked-mode flags, silently ignored.
#
# Output (in <outdir>/mds_perchr/):
#   <focal_chr>.mds_perchr.rds   list($out_dt, $dmat, $mds, $n_focal)
#
# Backward compat: ZO_G falls back to <outdir>/tmp/ if mds_perchr/ is absent,
# so data from older runs continues to work without renaming.
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
NPC        <- 4L           # top-NPC eigvals/vecs per window for the lostruct
                           # distance kernel. Must match the NPC used by the
                           # per-window PCA upstream (ZO_C / STEP_A02|A03).
MDS_DIMS   <- 5L           # Top-5 MDS axes. Base sim_mat is built directly
                           # from dmat_focal in ZO_G (not from MDS coords),
                           # so higher dims contribute nothing to L1/L2.
                           # NN-smoothed sim_mats average MDS coords across
                           # neighbours, so K=5 captures top variance with
                           # margin. K=2 only if you got time to test the
                           # eigenvalue spectrum per chrom; here we don't.
Z_THRESH   <- 3.0          # legacy; Z #1 (per-axis MDS-Z) is descriptive
                           # only and not consumed by L1/L2. Kept for atlas
                           # Z-profile plots.

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
  } else if (a == "--npc" && i < length(args)) {
    NPC <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--mds_dims" && i < length(args)) {
    MDS_DIMS <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--z_thresh" && i < length(args)) {
    Z_THRESH <- as.numeric(args[i + 1]); i <- i + 2L
  } else if (a == "--mds_mode" || a == "--seed") {
    # Legacy chunked-mode flags removed 2026-05-13; consume value and ignore.
    i <- i + 2L
  } else {
    i <- i + 1L
  }
}

if (is.null(rds_dir) || is.null(outdir) || is.null(FOCAL_CHR)) {
  stop("Usage: Rscript STEP_ZO_E_mds_compute.R --rds_dir <dir> --outdir <dir> --focal_chr <chr> ...")
}

# Per-chrom MDS results land here. Used to be called "tmp/" back when ZO_F
# merged these into a single inversion_localpca.mds.rds and downstream
# consumed only the merged file. ZO_F is gone now — these per-chrom RDS
# files ARE the final ZO_E output, consumed directly by ZO_G. The new name
# reflects that. ZO_G falls back to "tmp/" if "mds_perchr/" is absent, so
# data from older runs (still in tmp/) loads without renaming.
mds_perchr_dir <- file.path(outdir, "mds_perchr")
dir.create(mds_perchr_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SKIP-IF-DONE
# =============================================================================

rds_out <- file.path(mds_perchr_dir, paste0(FOCAL_CHR, ".mds_perchr.rds"))

if (file.exists(rds_out)) {
  message("[ZO_E] ", FOCAL_CHR, ": output already exists, skipping")
  message("  ", rds_out)
  quit(status = 0)
}

message("[ZO_E] ", FOCAL_CHR, ": nPC=", NPC, " mds_dims=", MDS_DIMS, " z_thresh=", Z_THRESH)

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

pc_dist <- function(pca_df, sample_names, npc, normalize = "L1",
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

dt_focal    <- merge(cd$meta, cd$pca, by = "global_window_id")
dt_combined <- dt_focal

message("[ZO_E] ", FOCAL_CHR, ": ", n_focal, " focal windows")

# =============================================================================
# DISTANCE MATRIX
# =============================================================================

n_total <- nrow(dt_combined)
n_pairs <- as.double(n_total) * (as.double(n_total) - 1) / 2

message("[STEP10v2-S1] Computing distance matrix: ", n_total, " windows (",
        format(round(n_pairs), big.mark = ","), " pairs)")

pb_dist <- make_progress(n_pairs, "dist")

dmat <- pc_dist(
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
  # If an inversion occupies a large fraction of the chromosome, mean/SD
  # normalization pulls toward the inversion signal and compresses z for
  # everything. median/MAD is resistant because the inversion windows are
  # outliers that don't move the median.
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
  out_dt  = out_chr,
  dmat    = dmat_focal,
  mds     = list(points = mds_mat_focal, eig = mds$eig),
  n_focal = sum(focal_mask)
)
saveRDS(result, rds_out)

meta_row <- data.table(
  focal_chrom     = FOCAL_CHR,
  n_focal_windows = result$n_focal,
  mds_dims        = MDS_DIMS,
  z_thresh        = Z_THRESH,
  elapsed_sec     = elapsed,
  timestamp       = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
fwrite(meta_row, file.path(mds_perchr_dir, paste0(FOCAL_CHR, ".metadata.tsv")), sep = "\t")

message("")
message("[DONE] ZO_E — ", FOCAL_CHR, ": ", nrow(out_chr),
        " focal windows, ", n_outlier, " z-outliers (|z|>=", Z_THRESH, "), ",
        elapsed, "s")
message("  ", rds_out)
