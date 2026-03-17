#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Concatenate per-OG trimmed alignments (FASTA) into a supermatrix.

Assumptions:
- Each *.trimmed.aln is a FASTA alignment.
- Headers are species names (first token after '>').
- Ideally each OG contains the same species set; if some species missing,
  gaps of alignment length will be inserted (safe fallback).

Outputs:
- concatenated FASTA (one record per species, in ref order)
- partitions file (start-end per OG, 1-based inclusive)

Example:
  python3 concat_trimmed_alignments.py \
    --aln_dir ../trimmed \
    --pattern "*.trimmed.aln" \
    --ref species_193.txt \
    --out_fasta supermatrix.fasta \
    --out_partitions partitions.txt
Optional:
  --og_list og_order.txt   # one filename (or OG id) per line
"""

import argparse
import glob
import os
from typing import Dict, List, Tuple


def load_ref(ref_path: str) -> List[str]:
    ref = []
    with open(ref_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s:
                ref.append(s)
    return ref


def read_fasta_alignment(fp: str) -> Dict[str, str]:
    """
    Read FASTA alignment, return dict: species -> sequence (concatenated, no spaces).
    """
    seqs: Dict[str, List[str]] = {}
    current = None
    with open(fp, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
                if current not in seqs:
                    seqs[current] = []
                continue
            if current is None:
                continue
            # keep alignment chars including '-' but remove spaces/tabs
            s = line.replace(" ", "").replace("\t", "").replace("\r", "")
            if s:
                seqs[current].append(s)
    out = {k: "".join(v) for k, v in seqs.items()}
    return out


def load_og_list(path: str) -> List[str]:
    items = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                items.append(s)
    return items


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln_dir", required=True, help="Directory with per-OG trimmed alignments")
    ap.add_argument("--pattern", default="*.trimmed.aln", help='Glob pattern (default "*.trimmed.aln")')
    ap.add_argument("--ref", required=True, help="species_193.txt (one species per line; output order)")
    ap.add_argument("--out_fasta", default="supermatrix.fasta", help="Output concatenated FASTA")
    ap.add_argument("--out_partitions", default="partitions.txt", help="Output partitions file")
    ap.add_argument("--og_list", default="", help="Optional OG order list (filenames or OG IDs)")
    ap.add_argument("--partition_prefix", default="gene", help="Prefix in partitions file (e.g. gene or LG)")
    args = ap.parse_args()

    ref = load_ref(args.ref)
    if not ref:
        raise SystemExit("[ERROR] reference species list is empty")

    # choose file order
    all_files = sorted(glob.glob(os.path.join(args.aln_dir, args.pattern)))
    if not all_files:
        raise SystemExit(f"[ERROR] no files matched: {os.path.join(args.aln_dir, args.pattern)}")

    if args.og_list:
        wanted = load_og_list(args.og_list)
        # build mapping from basename and OG id to full path
        by_base = {os.path.basename(p): p for p in all_files}
        by_stem = {os.path.splitext(os.path.basename(p))[0]: p for p in all_files}
        ordered_files = []
        missing = []
        for x in wanted:
            if x in by_base:
                ordered_files.append(by_base[x])
            elif x in by_stem:
                ordered_files.append(by_stem[x])
            else:
                missing.append(x)
        if missing:
            raise SystemExit(f"[ERROR] og_list contains items not found in aln_dir: {missing[:10]} ...")
    else:
        ordered_files = all_files

    # init concatenated sequences
    concat = {sp: [] for sp in ref}

    # partitions coordinates (1-based inclusive)
    pos = 1
    parts: List[Tuple[str, int, int]] = []

    for i, fp in enumerate(ordered_files, start=1):
        og_name = os.path.splitext(os.path.basename(fp))[0]  # e.g. OG0001234.trimmed
        aln = read_fasta_alignment(fp)

        # infer alignment length (within-file lengths should already be consistent)
        lengths = {len(seq) for seq in aln.values()}
        if not lengths:
            raise SystemExit(f"[ERROR] empty alignment: {fp}")
        if len(lengths) != 1:
            raise SystemExit(f"[ERROR] inconsistent lengths inside alignment (should be trimmed already): {fp}")
        L = next(iter(lengths))

        # append per species; if missing, fill gaps
        gap = "-" * L
        for sp in ref:
            concat[sp].append(aln.get(sp, gap))

        start = pos
        end = pos + L - 1
        parts.append((og_name, start, end))
        pos = end + 1

    # write concatenated FASTA
    with open(args.out_fasta, "w", encoding="utf-8") as out:
        for sp in ref:
            seq = "".join(concat[sp])
            out.write(f">{sp}\n")
            # wrap 60 chars per line
            for j in range(0, len(seq), 60):
                out.write(seq[j:j+60] + "\n")

    # write partitions (IQ-TREE style: name = start-end)
    with open(args.out_partitions, "w", encoding="utf-8") as out:
        for og_name, start, end in parts:
            out.write(f"{args.partition_prefix},{og_name} = {start}-{end}\n")

    total_len = pos - 1
    print(f"[DONE] Concatenated {len(ordered_files)} alignments.")
    print(f"[DONE] Wrote: {args.out_fasta}")
    print(f"[DONE] Wrote: {args.out_partitions}")
    print(f"[INFO] Species: {len(ref)}; Total length: {total_len}")

if __name__ == "__main__":
    main()

