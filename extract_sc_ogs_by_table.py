#!/usr/bin/env python3
# pip install biopython pandas

import argparse
import os
import pandas as pd
from collections import defaultdict
from Bio import SeqIO


def load_species_fasta_map(proteomes_dir: str):
    """
    Assume each species has its own fasta file, e.g. F_asiaticum.faa
    """
    m = {}
    for fn in os.listdir(proteomes_dir):
        if not fn.endswith((".fa", ".faa", ".fasta", ".fas")):
            continue
        sp = fn.split(".")[0]   # F_asiaticum from F_asiaticum.faa
        m[sp] = os.path.join(proteomes_dir, fn)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table_long", required=True,
                    help="TSV with columns: Orthogroup,Species,Gene")
    ap.add_argument("--proteomes_dir", required=True,
                    help="Directory containing per-species protein fasta")
    ap.add_argument("--out_dir", default="over50%singlecopygene",
                    help="Output directory for OG fasta files")
    ap.add_argument("--min_coverage", type=float, default=0.5,
                    help="Keep OG only if present species / total species >= this")
    ap.add_argument("--total_species", type=int, default=193,
                    help="Total species number used to compute coverage")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.table_long, sep="\t", dtype=str).fillna("")
    required_cols = {"Orthogroup", "Species", "Gene"}
    if not required_cols.issubset(df.columns):
        raise SystemExit(f"[ERROR] table_long must contain columns: {required_cols}")

    sp2fa = load_species_fasta_map(args.proteomes_dir)

    written_ogs = 0
    skipped_by_cov = 0
    missing_gene_records = 0
    dup_species_events = 0

    for og, sub in df.groupby("Orthogroup", sort=True):
        sp_gene = sub.groupby("Species")["Gene"].apply(list).to_dict()

        present_species = [sp for sp in sp_gene if sp in sp2fa]
        cov = len(present_species) / args.total_species
        if cov < args.min_coverage:
            skipped_by_cov += 1
            continue

        records_out = []
        species_counter = defaultdict(int)

        for sp, genes in sp_gene.items():
            fa = sp2fa.get(sp)
            if not fa:
                continue

            for gid in genes:
                found = None
                for rec in SeqIO.parse(fa, "fasta"):
                    if rec.id == gid:
                        found = rec
                        break

                if found is None:
                    missing_gene_records += 1
                    continue

                species_counter[sp] += 1
                idx = species_counter[sp]

                if idx == 1:
                    header = sp
                else:
                    header = f"{sp}{idx}"
                    dup_species_events += 1

                outrec = found[:]
                outrec.id = header
                outrec.description = ""
                records_out.append(outrec)

        if records_out:
            out_fp = os.path.join(args.out_dir, f"{og}.faa")
            SeqIO.write(records_out, out_fp, "fasta")
            written_ogs += 1

    print(f"[DONE] Wrote {written_ogs} OG fasta files to: {args.out_dir}")
    print(f"[INFO] Skipped by coverage (<{args.min_coverage}): {skipped_by_cov}")
    print(f"[INFO] Missing gene records: {missing_gene_records}")
    print(f"[INFO] Duplicate species entries renamed (sp2, sp3, ...): {dup_species_events}")


if __name__ == "__main__":
    main()
