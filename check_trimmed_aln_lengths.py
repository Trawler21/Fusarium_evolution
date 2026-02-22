#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Check whether all sequences inside each *.trimmed.aln (FASTA) file
have identical lengths (within-file), after concatenating multi-line sequences.

Usage:
  python3 check_trimmed_aln_lengths.py \
    --dir /path/to/trimmed_alignments \
    --pattern "*.trimmed.aln" \
    --report report.tsv

Exit code:
  0: all files OK
  2: at least one file has inconsistent lengths or other issues
"""

import argparse
import glob
import os
import sys
from collections import Counter
from typing import Dict, List, Tuple


def read_fasta_lengths(fp: str) -> Tuple[List[int], int, int]:
    """
    Return:
      lengths: list of sequence lengths (after concatenation)
      n_empty: number of records with empty sequence
      n_records: number of fasta records
    """
    lengths: List[int] = []
    n_empty = 0
    n_records = 0

    seq_chunks: List[str] = []
    in_record = False

    def flush_record():
        nonlocal seq_chunks, n_empty, n_records, in_record
        if not in_record:
            return
        n_records += 1
        seq = "".join(seq_chunks)
        # Remove whitespace characters inside sequence lines
        seq = seq.replace(" ", "").replace("\t", "").replace("\r", "").replace("\n", "")
        if len(seq) == 0:
            n_empty += 1
        lengths.append(len(seq))
        seq_chunks = []
        in_record = False

    with open(fp, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line:
                continue
            line = line.rstrip("\n")
            if line.startswith(">"):
                # finish previous record
                flush_record()
                in_record = True
                continue
            # sequence line
            if not in_record:
                # tolerate leading junk before first header; ignore
                continue
            if line.strip() == "":
                continue
            seq_chunks.append(line)

    flush_record()
    return lengths, n_empty, n_records


def summarize_counter(c: Counter) -> str:
    return " ".join(f"{k}({v})" for k, v in sorted(c.items(), key=lambda x: x[0]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="Directory containing *.trimmed.aln files")
    ap.add_argument("--pattern", default="*.trimmed.aln", help='Glob pattern (default: "*.trimmed.aln")')
    ap.add_argument("--report", default="", help="Optional TSV report path")
    ap.add_argument("--verbose", action="store_true", help="Print OK files too")
    args = ap.parse_args()

    target_glob = os.path.join(args.dir, args.pattern)
    files = sorted(glob.glob(target_glob))
    if not files:
        print(f"[ERROR] No files matched: {target_glob}", file=sys.stderr)
        return 2

    bad = 0
    rows = []
    for fp in files:
        lengths, n_empty, n_records = read_fasta_lengths(fp)
        base = os.path.basename(fp)

        if n_records == 0:
            bad += 1
            msg = "[BAD] no FASTA records"
            print(f"[BAD] {base}\t{msg}")
            rows.append((base, "BAD", "no_records", "", 0, n_empty))
            continue

        c = Counter(lengths)
        if n_empty > 0:
            bad += 1
            msg = f"[BAD] empty_seq_records={n_empty}; lengths: {summarize_counter(c)}"
            print(f"[BAD] {base}\t{msg}")
            rows.append((base, "BAD", "empty_seq", summarize_counter(c), n_records, n_empty))
            continue

        if len(c) > 1:
            bad += 1
            msg = f"[BAD] inconsistent_lengths; lengths: {summarize_counter(c)}"
            print(f"[BAD] {base}\t{msg}")
            rows.append((base, "BAD", "inconsistent", summarize_counter(c), n_records, 0))
        else:
            only_len = next(iter(c.keys()))
            if args.verbose:
                print(f"[OK]  {base}\tlen={only_len}\tn={n_records}")
            rows.append((base, "OK", "ok", str(only_len), n_records, 0))

    if args.report:
        with open(args.report, "w", encoding="utf-8") as out:
            out.write("file\tstatus\treason\tlength_info\tn_records\tn_empty\n")
            for r in rows:
                out.write("\t".join(map(str, r)) + "\n")

    if bad == 0:
        print(f"[SUMMARY] OK: {len(files)} files. All within-file sequence lengths are consistent.")
        return 0
    else:
        print(f"[SUMMARY] BAD: {bad} / {len(files)} files. See messages above.")
        return 2


if __name__ == "__main__":
    sys.exit(main())
