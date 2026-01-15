#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict

def read_fasta(path: str) -> dict[str, str]:
    seqs = {}
    sid = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf)
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if sid is not None:
            seqs[sid] = "".join(buf)
    return seqs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls_tsv", required=True)
    ap.add_argument("--nr_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    want = defaultdict(lambda: defaultdict(set))
    families = set()

    with open(args.calls_tsv) as f:
        header = f.readline()
        for line in f:
            species, prot, fam, _ = line.rstrip("\n").split("\t", 3)
            want[species][fam].add(prot)
            families.add(fam)

    os.makedirs(args.out_dir, exist_ok=True)

    out_handles = {}
    for fam in sorted(families):
        out_handles[fam] = open(os.path.join(args.out_dir, f"{fam}.fa"), "w")

    for fn in sorted(os.listdir(args.nr_dir)):
        if not fn.endswith(".nr.fa"):
            continue
        species = fn[:-len(".nr.fa")]
        fasta = read_fasta(os.path.join(args.nr_dir, fn))

        for fam, prots in want.get(species, {}).items():
            oh = out_handles[fam]
            for p in sorted(prots):
                if p not in fasta:
                    continue
                oh.write(f">{species}|{p}\n")
                seq = fasta[p]
                for i in range(0, len(seq), 60):
                    oh.write(seq[i:i+60] + "\n")

    for h in out_handles.values():
        h.close()

    print(f"[OK] Wrote family FASTAs to: {args.out_dir}")

if __name__ == "__main__":
    main()

