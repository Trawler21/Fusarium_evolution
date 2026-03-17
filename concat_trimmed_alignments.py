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
from collections import OrderedDict


def read_fasta(path):
    seqs = OrderedDict()
    cur = None
    chunks = []
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur is not None:
                    seqs[cur] = ''.join(chunks)
                cur = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if cur is not None:
        seqs[cur] = ''.join(chunks)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True)
    ap.add_argument('--pattern', default='*.trimmed.aln')
    ap.add_argument('--out_fasta', default='concatenated_alignment.fasta')
    ap.add_argument('--out_partition', default='partition.tsv')
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.dir, args.pattern)))
    if not files:
        raise SystemExit(f'[ERROR] no files match {os.path.join(args.dir, args.pattern)}')

    species = OrderedDict()
    partitions = []
    start = 1

    for fp in files:
        og = os.path.basename(fp)
        seqs = read_fasta(fp)
        if not seqs:
            raise SystemExit(f'[ERROR] empty alignment: {fp}')

        lengths = {len(v) for v in seqs.values()}
        if len(lengths) != 1:
            raise SystemExit(f'[ERROR] inconsistent sequence lengths in {fp}')
        seg_len = next(iter(lengths))

        end = start + seg_len - 1
        partitions.append((og, start, end, seg_len))

        for sp in species:
            if sp not in seqs:
                species[sp].append('-' * seg_len)

        for sp, seq in seqs.items():
            if sp not in species:
                species[sp] = ['-' * (start - 1)]
            species[sp].append(seq)

        start = end + 1

    with open(args.out_fasta, 'w', encoding='utf-8') as out:
        for sp, parts in species.items():
            out.write(f'>{sp}\n')
            out.write(''.join(parts) + '\n')

    with open(args.out_partition, 'w', encoding='utf-8') as out:
        out.write('locus\tstart\tend\tlength\n')
        for og, s, e, l in partitions:
            out.write(f'{og}\t{s}\t{e}\t{l}\n')

    print(f'[DONE] concatenated taxa: {len(species)}')
    print(f'[DONE] loci: {len(partitions)}')
    print(f'[DONE] alignment: {args.out_fasta}')
    print(f'[DONE] partition: {args.out_partition}')


if __name__ == '__main__':
    main()
