#!/usr/bin/env python3
# pip install pandas

import argparse
import pandas as pd

def split_genes(cell: str):
    cell = (cell or "").strip()
    if not cell:
        return []
    return [x.strip() for x in cell.split(",") if x.strip()]

def load_sc_ogs(sc_file: str):
    sc = set()
    with open(sc_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.lower().startswith("orthogroup"):
                continue
            og = line.split()[0]
            if og.startswith("OG"):
                sc.add(og)
    return sc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--og_tsv", required=True, help="Orthogroups.tsv")
    ap.add_argument("--sc_file", required=True, help="Orthogroups_SingleCopyOrthologues.txt")
    ap.add_argument("--out_wide", default="sc_og_species_gene_wide.tsv",
                    help="Output wide table: OG x species -> gene")
    ap.add_argument("--out_long", default="sc_og_species_gene_long.tsv",
                    help="Output long table: OG, species, gene")
    ap.add_argument("--expect_species_n", type=int, default=193,
                    help="Expected number of species columns (optional check)")
    args = ap.parse_args()

    sc_ogs = load_sc_ogs(args.sc_file)
    print(f"[INFO] Loaded single-copy OGs: {len(sc_ogs)}")

    df = pd.read_csv(args.og_tsv, sep="\t", dtype=str).fillna("")
    species_cols = list(df.columns[1:])
    if args.expect_species_n and len(species_cols) != args.expect_species_n:
        print(f"[WARN] Species columns in Orthogroups.tsv = {len(species_cols)}, "
              f"expected {args.expect_species_n}. (This may be OK if your run has different species count.)")

    df_sc = df[df.iloc[:, 0].isin(sc_ogs)].copy()
    print(f"[INFO] Rows kept in Orthogroups.tsv that are in sc_ogs: {df_sc.shape[0]}")

    out_wide = []
    out_long = []

    for _, row in df_sc.iterrows():
        og = row.iloc[0].strip()
        wide_row = {"Orthogroup": og}
        for sp in species_cols:
            genes = split_genes(row[sp])
            if len(genes) == 0:
                wide_row[sp] = ""
            elif len(genes) == 1:
                wide_row[sp] = genes[0]
                out_long.append((og, sp, genes[0]))
            else:           
                wide_row[sp] = "MULTI:" + ",".join(genes)
                for g in genes:
                    out_long.append((og, sp, g))
        out_wide.append(wide_row)

    wide_df = pd.DataFrame(out_wide)
    wide_df.to_csv(args.out_wide, sep="\t", index=False)

    long_df = pd.DataFrame(out_long, columns=["Orthogroup", "Species", "Gene"])
    long_df.to_csv(args.out_long, sep="\t", index=False)
   
    multi_cells = wide_df[species_cols].astype(str).applymap(lambda x: x.startswith("MULTI:")).to_numpy().sum()
    print(f"[INFO] MULTI cells (species has >1 gene in this OG): {int(multi_cells)}")


if __name__ == "__main__":
    main()
