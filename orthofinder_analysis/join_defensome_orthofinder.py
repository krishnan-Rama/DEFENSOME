#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict

def load_orthogroups_tsv(path: str):
    """
    Returns dict[(species, gene_id)] = orthogroup
    Orthogroups.tsv format: Orthogroup <tab> Species1 <tab> Species2 ...
    Cells contain comma-separated gene IDs.
    """
    mapping = {}
    with open(path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        species_cols = header[1:]
        for row in r:
            if not row:
                continue
            og = row[0]
            for sp, cell in zip(species_cols, row[1:]):
                cell = (cell or "").strip()
                if not cell:
                    continue
                genes = [g.strip() for g in cell.split(",") if g.strip()]
                for g in genes:
                    mapping[(sp, g)] = og
    return mapping

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--defensome_calls", required=True)
    ap.add_argument("--orthogroups_tsv", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--id_mode", choices=["plain", "species_pipe_protein"], default="plain",
                    help="How defensome protein IDs relate to OrthoFinder IDs.")
    args = ap.parse_args()

    og_map = load_orthogroups_tsv(args.orthogroups_tsv)

    total = 0
    missing = 0

    with open(args.defensome_calls, newline="") as fin, open(args.out_tsv, "w", newline="") as fout:
        r = csv.DictReader(fin, delimiter="\t")
        fieldnames = r.fieldnames + ["orthogroup"]
        w = csv.DictWriter(fout, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()

        for row in r:
            total += 1
            sp = row["species"]
            prot = row["protein"]

            if args.id_mode == "plain":
                key_gene = prot
            else:
                # OrthoFinder sometimes uses Species|Protein style IDs
                key_gene = f"{sp}|{prot}"

            og = og_map.get((sp, key_gene), "")
            if not og:
                missing += 1

            row["orthogroup"] = og
            w.writerow(row)

    print(f"[OK] Wrote: {args.out_tsv}")
    print(f"[INFO] Total defensome calls: {total}")
    print(f"[INFO] Missing OG assignments: {missing} ({missing/total:.3%})")
    if missing > 0:
        print("[WARN] If missing is high, your OrthoFinder gene IDs likely differ. Try --id_mode species_pipe_protein.")

if __name__ == "__main__":
    main()

