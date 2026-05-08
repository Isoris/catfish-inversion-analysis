#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR_A_compute_theta_matrices.R
# =============================================================================
# Heavy precomputation step for the theta-pi local-PCA / MDS path.
#
# Reads per-sample ANGSD pestPG files for one chromosome and emits a
# long-format TSV with one row per (sample, window) carrying the pairwise
# nucleotide diversity θπ and the per-window callable-site count.
#
# This is the heavy half of an A/B pair: STEP_TR_B_classify_theta.R is the
# light classifier that consumes the TSV, computes per-window population
# metrics + robust |Z|, runs per-window local PCA on the samples × windows
# θπ matrix, calls L2 envelopes from contiguous high-|Z| runs, and writes
# the consolidated atlas JSON.
#
# pestPG column layout (verified against ANGSD thetaStat output)
# --------------------------------------------------------------
#   col  1: (idxStart,idxEnd)(posStart,posEnd)(winStart,winEnd)
#   col  2: Chr
#   col  3: WinCenter
#   col  4: tW
#   col  5: tP            ← per-window SUM of per-site θπ (NOT a density)
#   col  6..13: tF, tH, tL, Tajima, fuf, fud, fayh, zeng
#   col 14: nSites        ← per-window callable-site count
#
# This script computes per-site θπ as tP / nSites (when nSites > 0). See
# README_theta_pi_scaling.md and ANGSD GH issue #329 for the rationale —
# in our 226-sample 9× cohort, nSites varies sample-to-sample so dividing
# is required before any cross-sample / cross-window comparison.
#
# Window-grid policy
# ------------------
# THETA_GRID_MODE controls how θπ values are assigned to window indices:
#   "native"  (default): use the pestPG grid as-is. At win10000.step2000
#                        this yields ~16,500 windows per chromosome. The
#                        atlas reconciles this with the dosage grid via
#                        Int32Array lookup tables built at JSON load.
#   "dosage":            pestPG values are nearest-midpoint joined onto
#                        the dosage scrubber's variable-bp window grid.
#                        Fallback if the native-grid per-sample loadings
#                        prove too noisy on real data — one config line
#                        to switch.
#
# Inputs (configured in 00_theta_config.sh)
# -----------------------------------------
#   $PESTPG_DIR/$SAMPLE.$PESTPG_SCALE.pestPG   per-sample pestPG files
#   $SAMPLE_LIST                                one sample id per line
#
# In dosage mode also requires:
#   $DOSAGE_WIN_BED_DIR/$CHROM/windows.bed     from the dosage scrubber
#
# Output
# ------
#   native mode: $THETA_TSV_DIR/theta_native.<CHROM>.<SCALE>.tsv.gz
#   dosage mode: $THETA_TSV_DIR/theta_dgrid.<CHROM>.tsv.gz
#   columns:     sample  chrom  window_idx  start_bp  end_bp
#                theta_pi (= tP/nSites, per-site)  tP_sum (raw tP)  n_sites
#   window_idx is 0-based and references the chosen output grid.
#
# Usage
# -----
#   source 00_theta_config.sh
#   Rscript STEP_TR_A_compute_theta_matrices.R --chrom <CHROM>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

# Named-arg parser
CHROM <- NULL
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom" && i < length(args)) {
    CHROM <- args[i + 1]; i <- i + 2L
  } else {
    i <- i + 1L
  }
}

if (is.null(CHROM)) {
  stop("Usage: Rscript STEP_TR_A_compute_theta_matrices.R --chrom <CHROM>")
}

# ── Config (from 00_theta_config.sh) ──────────────────────────────────────
PESTPG_DIR         <- Sys.getenv("PESTPG_DIR",         unset = NA)
PESTPG_SCALE       <- Sys.getenv("PESTPG_SCALE",       unset = "win10000.step2000")
SAMPLE_LIST        <- Sys.getenv("SAMPLE_LIST",        unset = NA)
THETA_TSV_DIR      <- Sys.getenv("THETA_TSV_DIR",      unset = NA)
THETA_GRID_MODE    <- Sys.getenv("THETA_GRID_MODE",    unset = "native")
DOSAGE_WIN_BED_DIR <- Sys.getenv("DOSAGE_WIN_BED_DIR", unset = NA)

stopifnot(!is.na(PESTPG_DIR), !is.na(SAMPLE_LIST), !is.na(THETA_TSV_DIR))
stopifnot(THETA_GRID_MODE %in% c("native", "dosage"))

# Output naming reflects the mode so both fallback paths can coexist on disk
out_basename <- if (THETA_GRID_MODE == "native") {
  sprintf("theta_native.%s.%s.tsv.gz", CHROM, PESTPG_SCALE)
} else {
  sprintf("theta_dgrid.%s.tsv.gz", CHROM)
}
out_tsv <- file.path(THETA_TSV_DIR, out_basename)

samples <- readLines(SAMPLE_LIST)
samples <- samples[nchar(samples) > 0]
message("[STEP_TR_A] CHROM=", CHROM,
        "  n_samples=", length(samples),
        "  scale=", PESTPG_SCALE,
        "  mode=", THETA_GRID_MODE)
message("           output = ", out_tsv)

# ── Helper: load one sample's pestPG, restrict to CHROM ───────────────────
load_pestpg <- function(sample) {
  pestpg_file <- file.path(PESTPG_DIR,
                           sprintf("%s.%s.pestPG", sample, PESTPG_SCALE))
  if (!file.exists(pestpg_file)) {
    return(NULL)
  }
  # pestPG header: #(...) Chr WinCenter tW tP tF tH tL Tajima ... nSites
  pp <- tryCatch(
    fread(cmd = sprintf("zcat -f %s | grep -v '^#'", shQuote(pestpg_file)),
          header = FALSE, fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(pp) || nrow(pp) < 1 || ncol(pp) < 14) return(NULL)
  setnames(pp, 1:14,
           c("ix", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL",
             "Tajima", "fuf", "fud", "fayh", "zeng", "nSites"))
  pp <- pp[Chr == CHROM]
  if (nrow(pp) == 0) return(NULL)

  # ── θπ per-site normalization (ANGSD GH #329) ─────────────────────────
  # `tP` (col 5) is the per-window SUM of per-site θπ; `nSites` (col 14)
  # is the per-window callable-site count. Divide to get a per-site value
  # comparable across windows and samples. Coverage varies sample-to-
  # sample so without this normalization, deeper samples falsely look
  # more diverse in any window where they have more callable sites.
  # Edge windows (chromosome ends) have nSites = 0 → NA.
  # `tP_sum` and `n_sites` are preserved alongside so STEP_TR_B / the
  # atlas can apply their own min-nSites masks without re-reading pestPG.
  pp[, theta_pi_per_site := fifelse(
        as.integer(nSites) > 0L,
        as.numeric(tP) / as.numeric(nSites),
        NA_real_
      )]

  # Reconstruct start_bp / end_bp from WinCenter (pestPG only stores center;
  # use SCALE_LABEL parsing to recover window size).
  win_size <- as.integer(sub("^win([0-9]+)\\..*$", "\\1", PESTPG_SCALE))
  pp[, `:=`(
    start_bp = pmax(0L, as.integer(WinCenter) - win_size %/% 2L),
    end_bp   = as.integer(WinCenter) + win_size %/% 2L,
    mid_bp   = as.integer(WinCenter)
  )]
  pp[, .(start_bp, end_bp, mid_bp,
         theta_pi = theta_pi_per_site,    # per-site (tP / nSites) — see comment above
         tP_sum   = as.numeric(tP),       # raw window sum, preserved
         n_sites  = as.integer(nSites))]
}

# ── Main loop: one sample at a time, append to output ────────────────────
out_rows <- vector("list", length(samples))
n_ok <- 0L
n_skip <- 0L

if (THETA_GRID_MODE == "native") {

  # Native mode — θπ uses its own pestPG window grid. window_idx is 0-based
  # over the unique sorted set of WinCenter values across all samples on
  # this chromosome (they should all share the same grid since pestPG is
  # deterministic at fixed -win/-step).

  # First pass: collect canonical window list from a reference sample
  ref_sample <- NULL
  for (s in samples) {
    pp_ref <- load_pestpg(s)
    if (!is.null(pp_ref)) { ref_sample <- pp_ref; break }
  }
  if (is.null(ref_sample)) {
    stop("[STEP_TR_A] Could not load any sample's pestPG for ", CHROM)
  }
  setorder(ref_sample, mid_bp)
  windows_grid <- ref_sample[, .(start_bp, end_bp, mid_bp)]
  windows_grid[, window_idx := .I - 1L]
  setkey(windows_grid, mid_bp)
  message("[STEP_TR_A] θπ-native grid: ", nrow(windows_grid), " windows")

  # Second pass: assign each sample's pestPG values to window_idx by mid_bp
  for (i in seq_along(samples)) {
    s <- samples[i]
    pp <- load_pestpg(s)
    if (is.null(pp)) {
      warning("[STEP_TR_A] Missing/malformed pestPG for ", s, " — skip")
      n_skip <- n_skip + 1L; next
    }
    setkey(pp, mid_bp)
    # Exact match on mid_bp (windows are canonical across samples; if a
    # sample has missing entries for some windows, those become NA).
    joined <- windows_grid[pp, on = "mid_bp", nomatch = NA]
    # Drop rows where window_idx is NA (sample had a window not in canonical
    # grid — shouldn't happen at fixed -win/-step but defensive)
    joined <- joined[!is.na(window_idx)]
    joined[, sample := s]
    joined[, chrom  := CHROM]
    out_rows[[i]] <- joined[, .(sample, chrom, window_idx,
                                start_bp = i.start_bp,
                                end_bp   = i.end_bp,
                                theta_pi, tP_sum, n_sites)]
    n_ok <- n_ok + 1L
    if (i %% 50 == 0) {
      message("[STEP_TR_A] processed ", i, " / ", length(samples))
    }
  }

} else {

  # Dosage mode — fallback. θπ values nearest-midpoint joined to dosage
  # scrubber windows (the v3 design). Useful if v4's finer-scale loadings
  # turn out too noisy on real data.
  if (is.na(DOSAGE_WIN_BED_DIR)) {
    stop("[STEP_TR_A] dosage mode requires DOSAGE_WIN_BED_DIR set")
  }
  dosage_bed <- file.path(DOSAGE_WIN_BED_DIR, CHROM, "windows.bed")
  if (!file.exists(dosage_bed)) {
    stop("[STEP_TR_A] dosage mode: missing ", dosage_bed)
  }
  dosage <- fread(dosage_bed,
                  col.names = c("chrom", "start_bp", "end_bp",
                                "window_idx", "n_snps"))
  dosage <- dosage[chrom == CHROM]
  dosage[, mid_bp := as.integer((start_bp + end_bp) / 2L)]
  setkey(dosage, mid_bp)
  message("[STEP_TR_A] dosage grid: ", nrow(dosage), " windows")

  for (i in seq_along(samples)) {
    s <- samples[i]
    pp <- load_pestpg(s)
    if (is.null(pp)) {
      warning("[STEP_TR_A] Missing/malformed pestPG for ", s, " — skip")
      n_skip <- n_skip + 1L; next
    }
    pp_for_join <- pp[, .(mid_bp, theta_pi, tP_sum, n_sites)]
    setkey(pp_for_join, mid_bp)
    joined <- pp_for_join[dosage, roll = "nearest", on = "mid_bp"]
    joined[, sample := s]
    joined[, chrom  := CHROM]
    out_rows[[i]] <- joined[, .(sample, chrom, window_idx,
                                start_bp, end_bp,
                                theta_pi, tP_sum, n_sites)]
    n_ok <- n_ok + 1L
    if (i %% 50 == 0) {
      message("[STEP_TR_A] processed ", i, " / ", length(samples))
    }
  }
}

if (n_ok == 0) {
  stop("[STEP_TR_A] No samples produced output. Check PESTPG_DIR + sample naming.")
}

result <- rbindlist(out_rows, use.names = TRUE)
message("[STEP_TR_A] joined: n_ok=", n_ok,
        " n_skip=", n_skip,
        " total_rows=", nrow(result))

# ── Write output TSV ──────────────────────────────────────────────────────
fwrite(result, out_tsv, sep = "\t", compress = "gzip", na = "")
fi <- file.info(out_tsv)
message("[STEP_TR_A] Wrote ", out_tsv,
        " (", round(fi$size / 1024 / 1024, 2), " MB)")

message("[STEP_TR_A] DONE — chrom=", CHROM, " mode=", THETA_GRID_MODE)
