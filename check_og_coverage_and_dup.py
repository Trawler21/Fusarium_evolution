#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import re
from collections import Counter

re_suffix_digits = re.compile(r"^(.*?)(\d+)$")

def base_species(name: str) -> str:
    name = (name or "").strip()
    if not name:
        return ""
    name = name.split()[0]
    m = re_suffix_digits.match(name)
    if m:
        return m.group(1)
    return name

def load_ref(ref_path):
    ref_list = []
    with open(ref_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s:
                ref_list.append(s)
    return ref_list, set(ref_list)

def parse_headers(faa_path):
    headers = []
    with open(faa_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip()
                if h:
                    headers.append(h.split()[0])
    return headers

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="species_193.txt (one species per line)")
    ap.add_argument("--dir", required=True, help="Directory with OG .faa files")
    ap.add_argument("--pattern", default="*.faa", help='Glob pattern (default "*.faa")')
    ap.add_argument("--report", default="og_coverage_dup_report.tsv", help="Output TSV report")
    ap.add_argument("--write_lists", action="store_true",
                    help="Write per-OG missing/dup/extra lists to a folder")
    args = ap.parse_args()

    ref_list, ref_set = load_ref(args.ref)
    files = sorted(glob.glob(os.path.join(args.dir, args.pattern)))
    if not files:
        raise SystemExit(f"[ERROR] No files matched: {os.path.join(args.dir, args.pattern)}")

    outdir = None
    if args.write_lists:
        outdir = os.path.splitext(args.report)[0] + "_lists"
        os.makedirs(outdir, exist_ok=True)

    full_ok = 0
    any_dup = 0
    any_extra = 0

    with open(args.report, "w", encoding="utf-8") as out:
        out.write(
            "OG\tseq_n\tn_unique_base_species\tmissing_n\tdup_base_species_n\textra_n\t"
            "missing_species\tdup_base_species\textra_base_species\n"
        )

        for fp in files:
            og = os.path.basename(fp)
            headers = parse_headers(fp)
            seq_n = len(headers)

            base_sps = [base_species(h) for h in headers if base_species(h)]
            c = Counter(base_sps)
            uniq = set(c.keys())

            missing = [sp for sp in ref_list if sp not in uniq]
            dup = sorted([sp for sp, k in c.items() if k > 1])
            extra = sorted([sp for sp in uniq if sp not in ref_set])

            if len(missing) == 0 and len(dup) == 0 and len(extra) == 0:
                full_ok += 1
            if dup:
                any_dup += 1
            if extra:
                any_extra += 1

            out.write(
                f"{og}\t{seq_n}\t{len(uniq)}\t{len(missing)}\t{len(dup)}\t{len(extra)}\t"
                f"{','.join(missing)}\t{','.join(dup)}\t{','.join(extra)}\n"
            )

            if outdir and (missing or dup or extra):
                prefix = os.path.join(outdir, og)
                if missing:
                    with open(prefix + ".missing.txt", "w", encoding="utf-8") as f:
                        f.write("\n".join(missing) + "\n")
                if dup:
                    with open(prefix + ".dup.txt", "w", encoding="utf-8") as f:
                        for sp in dup:
                            f.write(f"{sp}\tcount={c[sp]}\n")
                if extra:
                    with open(prefix + ".extra.txt", "w", encoding="utf-8") as f:
                        f.write("\n".join(extra) + "\n")

    print(f"[DONE] Report written: {args.report}")
    print(f"[SUMMARY] Perfect OG (no missing/dup/extra): {full_ok} / {len(files)}")
    print(f"[SUMMARY] OG with dup species: {any_dup} / {len(files)}")
    print(f"[SUMMARY] OG with extra species: {any_extra} / {len(files)}")
    print("[NOTE] dup/extra are computed on base species (e.g. F_xxx2 -> F_xxx).")

if __name__ == "__main__":
    main()
