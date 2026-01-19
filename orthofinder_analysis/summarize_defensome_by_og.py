#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls_with_og", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--min_presence_frac", type=float, default=0.8)
    ap.add_argument("--expanded_max_copy", type=int, default=5)
    args = ap.parse_args()

    fam_og_sp_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    species_set = set()

    with open(args.calls_with_og, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            sp = row["species"]
            fam = row["family"]
            og = row.get("orthogroup", "")
            species_set.add(sp)
            if not og:
                continue
            fam_og_sp_counts[fam][og][sp] += 1

    species = sorted(species_set)
    nsp = len(species)
    if nsp == 0:
        raise SystemExit("[ERROR] No species detected")

    stats_out = args.out_prefix + ".family_orthogroup_stats.tsv"
    labels_out = args.out_prefix + ".family_orthogroup_labels.tsv"

    with open(stats_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family","orthogroup","presence_species","presence_frac","median_copy","max_copy","mean_copy","var_copy"])
        for fam in sorted(fam_og_sp_counts):
            for og in sorted(fam_og_sp_counts[fam]):
                copies = []
                present = 0
                for sp in species:
                    c = fam_og_sp_counts[fam][og].get(sp, 0)
                    copies.append(c)
                    if c > 0:
                        present += 1

                presence_frac = present / nsp
                copies_sorted = sorted(copies)
                mid = nsp // 2
                median = copies_sorted[mid] if nsp % 2 == 1 else 0.5*(copies_sorted[mid-1] + copies_sorted[mid])
                mean = sum(copies) / nsp
                var = sum((c - mean)**2 for c in copies) / nsp
                mx = max(copies)

                w.writerow([fam, og, present, f"{presence_frac:.3f}", f"{median:.3f}", mx, f"{mean:.3f}", f"{var:.3f}"])

    with open(stats_out, newline="") as fin, open(labels_out, "w", newline="") as fout:
        r = csv.DictReader(fin, delimiter="\t")
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["family","orthogroup","label","rule"])
        for row in r:
            presence_frac = float(row["presence_frac"])
            max_copy = int(row["max_copy"])
            median_copy = float(row["median_copy"])

            label = "other"
            rule = ""

            if presence_frac >= args.min_presence_frac and median_copy <= 1.1 and max_copy <= 2:
                label = "conserved"
                rule = f"presence_frac>={args.min_presence_frac} & median~1 & max<=2"
            elif max_copy >= args.expanded_max_copy:
                label = "expanded"
                rule = f"max_copy>={args.expanded_max_copy}"

            w.writerow([row["family"], row["orthogroup"], label, rule])

    print(f"[OK] Wrote:\n  {stats_out}\n  {labels_out}")

if __name__ == "__main__":
    main()

