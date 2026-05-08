#!/usr/bin/env python3
"""
STEP08_beagle_to_dosage_by_chr.py  (v2 — GL triplet archive support)

Read an ANGSD .beagle.gz file, compute expected genotype dosage per sample,
and split output by chromosome.

For each site:
  dosage = P(AB) + 2 * P(BB)

Outputs per chromosome:
  <outdir>/<chrom>.sites.tsv.gz
  <outdir>/<chrom>.dosage.tsv.gz

With --gl-triplet:
  <outdir>/<chrom>.gl_triplet.tsv.gz
  Format: marker, sample1_AA, sample1_AB, sample1_BB, sample2_AA, ...

The GL triplet preserves the full genotype likelihood vector, which is
irreversibly collapsed into a single scalar by the dosage computation.
This matters for within-stripe subgroup analysis where confident vs
uncertain genotypes must be distinguished.
"""

import argparse
import gzip
import os
import re
import sys
from typing import Dict, List, Optional, Tuple, TextIO


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beagle", required=True, help="Input ANGSD .beagle.gz file")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--digits", type=int, default=6, help="Decimal places")
    ap.add_argument("--prefix", default="", help="Optional output prefix")
    ap.add_argument("--gl-triplet", action="store_true",
                    help="Also output GL triplet archive (P(AA), P(AB), P(BB) per sample)")
    return ap.parse_args()


def infer_chr_pos(marker: str) -> Tuple[str, int]:
    marker2 = marker.replace(":", "_")
    parts = marker2.split("_")
    if len(parts) < 2:
        raise ValueError(f"Could not parse marker '{marker}'")

    pos_idx = None
    for i in range(len(parts) - 1, -1, -1):
        if re.fullmatch(r"\d+", parts[i]):
            pos_idx = i
            break

    if pos_idx is None:
        raise ValueError(f"Could not find numeric position in marker '{marker}'")

    pos = int(parts[pos_idx])
    chrom_parts = parts[:pos_idx]
    if not chrom_parts:
        raise ValueError(f"Could not infer chromosome from marker '{marker}'")

    chrom = "_".join(chrom_parts)
    return chrom, pos


def open_gzip_text_writer(path: str) -> TextIO:
    return gzip.open(path, "wt")


class ChromWriters:
    """Manages per-chromosome output file handles for dosage + optional GL triplet."""

    def __init__(self, outdir: str, prefix: str, samples: List[str],
                 write_gl: bool = False):
        self.outdir = outdir
        self.prefix = prefix
        self.samples = samples
        self.write_gl = write_gl
        self.handles: Dict[str, Tuple[TextIO, TextIO, Optional[TextIO]]] = {}

    def _base(self, chrom: str) -> str:
        return f"{self.prefix}{chrom}" if self.prefix else chrom

    def get(self, chrom: str) -> Tuple[TextIO, TextIO, Optional[TextIO]]:
        if chrom not in self.handles:
            base = self._base(chrom)
            sites_path = os.path.join(self.outdir, f"{base}.sites.tsv.gz")
            dosage_path = os.path.join(self.outdir, f"{base}.dosage.tsv.gz")

            sites_fh = open_gzip_text_writer(sites_path)
            dosage_fh = open_gzip_text_writer(dosage_path)

            sites_fh.write("marker\tchrom\tpos\tallele1\tallele2\n")
            dosage_fh.write("marker\t" + "\t".join(self.samples) + "\n")

            gl_fh: Optional[TextIO] = None
            if self.write_gl:
                gl_path = os.path.join(self.outdir, f"{base}.gl_triplet.tsv.gz")
                gl_fh = open_gzip_text_writer(gl_path)
                gl_cols = []
                for s in self.samples:
                    gl_cols.extend([f"{s}_AA", f"{s}_AB", f"{s}_BB"])
                gl_fh.write("marker\t" + "\t".join(gl_cols) + "\n")

            self.handles[chrom] = (sites_fh, dosage_fh, gl_fh)
        return self.handles[chrom]

    def close_all(self) -> None:
        for sites_fh, dosage_fh, gl_fh in self.handles.values():
            sites_fh.close()
            dosage_fh.close()
            if gl_fh is not None:
                gl_fh.close()


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    n_sites = 0
    chrom_counts: Dict[str, int] = {}

    with gzip.open(args.beagle, "rt") as fh:
        header = fh.readline().rstrip("\n").split()
        if len(header) < 6:
            raise RuntimeError("BEAGLE header too short")

        if (len(header) - 3) % 3 != 0:
            raise RuntimeError("Header genotype-likelihood columns after first 3 fields are not divisible by 3")

        raw_sample_cols = header[3:]
        n_samples = len(raw_sample_cols) // 3

        samples: List[str] = []
        for i in range(n_samples):
            s0 = raw_sample_cols[i * 3]
            s1 = raw_sample_cols[i * 3 + 1]
            s2 = raw_sample_cols[i * 3 + 2]
            if not (s0 == s1 == s2):
                raise RuntimeError(f"Expected repeated sample names in groups of 3, got: {s0}, {s1}, {s2}")
            samples.append(s0)

        writers = ChromWriters(args.outdir, args.prefix, samples,
                               write_gl=args.gl_triplet)
        fmt = "{:." + str(args.digits) + "f}"

        try:
            for line_num, line in enumerate(fh, start=2):
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split()

                expected_len = 3 + 3 * n_samples
                if len(fields) != expected_len:
                    raise RuntimeError(f"Line {line_num}: expected {expected_len} columns, found {len(fields)}")

                marker = fields[0]
                allele1 = fields[1]
                allele2 = fields[2]

                chrom, pos = infer_chr_pos(marker)
                sites_fh, dosage_fh, gl_fh = writers.get(chrom)

                sites_fh.write(f"{marker}\t{chrom}\t{pos}\t{allele1}\t{allele2}\n")

                dosages: List[str] = [marker]
                gl_vals: List[str] = [marker] if gl_fh else []
                gl = fields[3:]

                for i in range(n_samples):
                    p_aa = float(gl[i * 3])
                    p_ab = float(gl[i * 3 + 1])
                    p_bb = float(gl[i * 3 + 2])
                    dosage = p_ab + 2.0 * p_bb
                    dosages.append(fmt.format(dosage))

                    if gl_fh:
                        gl_vals.extend([
                            fmt.format(p_aa),
                            fmt.format(p_ab),
                            fmt.format(p_bb),
                        ])

                dosage_fh.write("\t".join(dosages) + "\n")
                if gl_fh:
                    gl_fh.write("\t".join(gl_vals) + "\n")

                n_sites += 1
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

        finally:
            writers.close_all()

    sys.stderr.write(f"[INFO] Samples: {n_samples}\n")
    sys.stderr.write(f"[INFO] Sites processed: {n_sites}\n")
    if args.gl_triplet:
        sys.stderr.write(f"[INFO] GL triplet archive: ENABLED\n")
    for chrom in sorted(chrom_counts):
        sys.stderr.write(f"[INFO] {chrom}: {chrom_counts[chrom]} sites\n")


if __name__ == "__main__":
    main()
