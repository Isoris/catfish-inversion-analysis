#!/usr/bin/env Rscript

# =============================================================================
# 08a_build_sample_metadata.R
#
# Build the canonical sample-metadata TSV by reconciling THREE INDEPENDENT
# identity layers, ONCE, genome-wide. The output is consumed by step 08b
# (atlas JSON exporter) for every chromosome.
#
# The three layers (kept conceptually separate; reconciled here):
#
#   1. BAMLIST          → maps Ind0..IndN -> CGAxxx (sequencing identity)
#                         one CGA ID per line, in the SAME order as the
#                         per-window precomp PC_1_Ind* columns.
#
#   2. NGSRELATE PAIRS  → builds the family graph via union-find on the
#                         (id1, id2, theta) relatedness table, using the
#                         theta cutoff (default 0.177, Manichaikul 1st-deg).
#                         Each connected component gets a family_id; isolated
#                         samples get family_id = -1.
#
#   3. ANCESTRY         → NGSadmix K=8 cluster assignment (or other ancestry
#                         tag), CGA -> ancestry_label. Optional; samples
#                         without a record get "unknown".
#
# The reconciliation logic used to live inside the JSON exporter
# (08b_export_atlas_json_localpca_zblocks.R). Centralizing it here:
#   - lets you sanity-check the merged TSV before paying the JSON build cost
#   - means the exporter takes ONE flag (--sample_metadata) instead of three
#   - makes the merged identity reusable across paths 1/2/3 (this same TSV
#     feeds the atlas JSON for localpca_zblocks, localpca_thetapi, and
#     localpca_GHSL — all three discovery paths use the same 226 samples)
#
# Codebase:    inversion-popgen-toolkit v8.5 / consolidated layout v1.0
# Upstream:    bamlist, ngsRelate output, NGSadmix output (or ancestry tsv)
# Downstream:  08b_export_atlas_json_localpca_zblocks.R (and its peers in
#              paths 2 / 3 once they're reorganized).
#
# ── Usage ──
#
#   Rscript 08a_build_sample_metadata.R \
#     --bamlist     <one CGA per line, in PC_1_Ind* order> \
#     --pairs       <3-col: id1 id2 theta from ngsRelate> \
#     [--theta_cutoff 0.177] \
#     [--ancestry   <tsv: cga<TAB>ancestry_label, one row per sample>] \
#     --out         <out.tsv>
#
# ── Output TSV columns ──
#   ind         : Ind0, Ind1, ..., IndN (matches PC_1_Ind* column suffix)
#   cga         : CGAxxx (real sample ID after bamlist remap)
#   family_id   : integer; samples in same family share id; isolates = -1
#   ancestry    : ancestry label (default "unknown")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

BAMLIST    <- get_arg("--bamlist")
PAIRS      <- get_arg("--pairs")
THETA_CUT  <- as.numeric(get_arg("--theta_cutoff", "0.177"))
ANCESTRY   <- get_arg("--ancestry")
OUT        <- get_arg("--out")

if (is.na(BAMLIST) || !file.exists(BAMLIST))
  stop("[08a] --bamlist is required and must exist")
if (is.na(OUT))
  stop("[08a] --out is required")

# ---- Layer 1: bamlist -> ind/cga --------------------------------------------
message("[08a] Loading bamlist: ", BAMLIST)
bam_ids <- readLines(BAMLIST)
bam_ids <- bam_ids[nzchar(trimws(bam_ids))]

# Strip path + .bam/.cram/.sam extensions to get clean CGA labels
clean_cga <- vapply(bam_ids, function(x) {
  base <- basename(x)
  sub("\\.(bam|cram|sam)$", "", base, ignore.case = TRUE)
}, character(1))

n_samp <- length(clean_cga)
sample_meta <- data.table(
  ind = paste0("Ind", seq_len(n_samp) - 1L),
  cga = clean_cga,
  family_id = -1L,
  ancestry  = "unknown"
)
message("[08a] bamlist samples: ", n_samp)

# ---- Layer 2: ngsRelate pairs -> family_id (union-find) ---------------------
build_components_from_pairs <- function(pairs_path, theta_thr, all_ids) {
  message("[08a] Loading pairs: ", pairs_path, "  theta_cutoff=", theta_thr)
  first_line <- readLines(pairs_path, n = 1L)
  parts <- strsplit(first_line, "[\t ]+")[[1]]
  has_header <- !suppressWarnings(is.finite(as.numeric(parts[length(parts)])))
  pairs_dt <- fread(pairs_path, header = has_header,
                    col.names = c("id1", "id2", "theta"))
  pairs_dt[, theta := as.numeric(theta)]
  qs <- quantile(pairs_dt$theta, c(0.5, 0.9, 0.95, 0.99, 1.0), na.rm = TRUE)
  message(sprintf("[08a]   theta dist: median=%.4f p90=%.4f p95=%.4f p99=%.4f max=%.4f",
                  qs[1], qs[2], qs[3], qs[4], qs[5]))

  # Filter to edges above cutoff
  edges <- pairs_dt[theta >= theta_thr & id1 %in% all_ids & id2 %in% all_ids]
  message("[08a]   edges above cutoff: ", nrow(edges), " / ", nrow(pairs_dt))

  # Union-find
  id_to_idx <- setNames(seq_along(all_ids), all_ids)
  parent <- seq_along(all_ids)
  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  union <- function(a, b) {
    ra <- find(a); rb <- find(b)
    if (ra != rb) parent[ra] <<- rb
  }
  for (i in seq_len(nrow(edges))) {
    a <- id_to_idx[edges$id1[i]]
    b <- id_to_idx[edges$id2[i]]
    if (!is.na(a) && !is.na(b)) union(a, b)
  }
  roots <- vapply(seq_along(all_ids), find, integer(1))

  # Assign family_id: components of size >= 2 get sequential ids; isolates = -1
  comp_size <- table(roots)
  big_roots <- as.integer(names(comp_size)[comp_size >= 2])
  fam_lookup <- setNames(seq_along(big_roots), as.character(big_roots))
  fam <- ifelse(as.character(roots) %in% names(fam_lookup),
                fam_lookup[as.character(roots)], -1L)
  message("[08a]   families (size >=2): ", length(big_roots),
          "   isolates: ", sum(fam == -1L))
  data.table(cga = all_ids, family_id = as.integer(fam))
}

if (!is.na(PAIRS) && file.exists(PAIRS)) {
  fam_dt <- build_components_from_pairs(PAIRS, THETA_CUT, sample_meta$cga)
  sample_meta[, family_id := NULL]
  sample_meta <- merge(sample_meta, fam_dt, by = "cga", all.x = TRUE, sort = FALSE)
  sample_meta[is.na(family_id), family_id := -1L]
} else {
  message("[08a] no --pairs given; family_id = -1 for all samples")
}

# ---- Layer 3: ancestry ------------------------------------------------------
if (!is.na(ANCESTRY) && file.exists(ANCESTRY)) {
  message("[08a] Loading ancestry: ", ANCESTRY)
  anc_dt <- fread(ANCESTRY, header = TRUE)
  setnames(anc_dt, tolower(names(anc_dt)))
  if (!"cga" %in% names(anc_dt) || !"ancestry" %in% names(anc_dt))
    stop("[08a] ancestry TSV must have columns: cga, ancestry")
  sample_meta[, ancestry := NULL]
  sample_meta <- merge(sample_meta, anc_dt[, .(cga, ancestry)],
                       by = "cga", all.x = TRUE, sort = FALSE)
  sample_meta[is.na(ancestry), ancestry := "unknown"]
  message("[08a]   ancestry assigned: ",
          sum(sample_meta$ancestry != "unknown"), " / ", n_samp)
} else {
  message("[08a] no --ancestry given; ancestry = 'unknown' for all samples")
}

# Restore original ind ordering (PC_1_Ind* alignment is critical)
ord <- match(paste0("Ind", seq_len(n_samp) - 1L), sample_meta$ind)
sample_meta <- sample_meta[ord]
stopifnot(identical(sample_meta$ind, paste0("Ind", seq_len(n_samp) - 1L)))

# Final column order
sample_meta <- sample_meta[, .(ind, cga, family_id, ancestry)]

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(sample_meta, OUT, sep = "\t")
message("[08a] Wrote: ", OUT, "   (", n_samp, " rows)")
message("[08a] Preview:")
print(head(sample_meta, 5))
