#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict

def parse_fasta_lengths(fa_path: str) -> dict[str, int]:
    lengths = {}
    seq_id = None
    seq = []
    with open(fa_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if seq_id is not None:
                    lengths[seq_id] = len("".join(seq))
                seq_id = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if seq_id is not None:
            lengths[seq_id] = len("".join(seq))
    return lengths

def load_map_tsv(path: str) -> dict[str, set[str]]:
    fam_to_ids = defaultdict(set)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            fam, acc = parts[0], parts[1]
            fam_to_ids[fam].add(acc)
    return fam_to_ids

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ipr_dir", required=True)
    ap.add_argument("--nr_dir", required=True)
    ap.add_argument("--pfam_map", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--min_dom_cov", type=float, default=0.35)
    ap.add_argument("--min_len_cyp", type=int, default=350)
    ap.add_argument("--min_len_abc_nbd", type=int, default=150)
    args = ap.parse_args()

    pfam_map = load_map_tsv(args.pfam_map)

    # Load fasta lengths for length sanity checks
    fasta_len = {}
    for fn in sorted(os.listdir(args.nr_dir)):
        if not fn.endswith(".nr.fa"):
            continue
        species = fn[:-len(".nr.fa")]
        fasta_len[species] = parse_fasta_lengths(os.path.join(args.nr_dir, fn))

    # Collect calls[(species, prot)][fam] = list(evidence strings)
    calls = defaultdict(lambda: defaultdict(list))
    qc = {}

    ipr_files = sorted([p for p in os.listdir(args.ipr_dir) if p.endswith(".ipr.tsv")])
    if not ipr_files:
        raise SystemExit(f"[ERROR] No *.ipr.tsv found in {args.ipr_dir}")

    for fn in ipr_files:
        species = fn[:-len(".ipr.tsv")]
        path = os.path.join(args.ipr_dir, fn)

        n_rows = 0
        prots_seen = set()
        assigned_hits = 0

        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                # InterProScan TSV columns (typical):
                # 0 prot, 2 length, 3 analysis, 4 signature_accession, 6 start, 7 end, 11 interpro_accession ...
                if len(parts) < 8:
                    continue

                prot = parts[0]
                prots_seen.add(prot)
                n_rows += 1

                try:
                    ipr_len = int(parts[2])
                except Exception:
                    ipr_len = None

                analysis = parts[3] if len(parts) > 3 else ""
                sig_acc = parts[4] if len(parts) > 4 else ""

                start = None
                end = None
                if len(parts) > 7:
                    try:
                        start = int(parts[6])
                        end = int(parts[7])
                    except Exception:
                        start = None
                        end = None

                # Domain coverage filter
                dom_cov = None
                if start is not None and end is not None and ipr_len:
                    if ipr_len > 0:
                        dom_cov = (end - start + 1) / float(ipr_len)
                        if dom_cov < args.min_dom_cov:
                            continue

                # Assign families by PFAM signature accession
                for fam, ids in pfam_map.items():
                    if sig_acc in ids:
                        prot_len = fasta_len.get(species, {}).get(prot, ipr_len or 0)

                        # sanity filters
                        if fam == "CYP" and prot_len < args.min_len_cyp:
                            continue
                        if fam == "ABC_NBD" and prot_len < args.min_len_abc_nbd:
                            continue

                        if dom_cov is not None:
                            ev = f"{analysis}:{sig_acc}:{start}-{end}:cov={dom_cov:.3f}"
                        else:
                            ev = f"{analysis}:{sig_acc}:{start}-{end}"
                        calls[(species, prot)][fam].append(ev)
                        assigned_hits += 1

        qc[species] = {
            "species": species,
            "ipr_rows": n_rows,
            "proteins_in_ipr": len(prots_seen),
            "proteins_in_fasta": len(fasta_len.get(species, {})),
            "assigned_hits": assigned_hits,
        }

    # Post-processing: composite rules
    prot_calls = []
    for (species, prot), fams in calls.items():
        famset = set(fams.keys())

        # GST strict: require both N and C
        if "GST_N" in famset and "GST_C" in famset:
            fams["GST_strict"] = fams["GST_N"] + fams["GST_C"]

        # ABC strict: require NBD + TMD
        if "ABC_NBD" in famset and "ABC_TMD" in famset:
            fams["ABC_strict"] = fams["ABC_NBD"] + fams["ABC_TMD"]

        for fam, evs in fams.items():
            prot_calls.append((species, prot, fam, ";".join(evs)))

    # species x family counts (unique proteins)
    species_family_to_prots = defaultdict(set)
    for species, prot, fam, _ in prot_calls:
        species_family_to_prots[(species, fam)].add(prot)

    species_list = sorted({s for s, _ in species_family_to_prots.keys()} | set(qc.keys()))
    fam_list = sorted({f for _, f in species_family_to_prots.keys()})

    out_dir = os.path.dirname(args.out_prefix)
    os.makedirs(out_dir, exist_ok=True)

    prot_calls_path = args.out_prefix + ".protein_calls.tsv"
    species_counts_path = args.out_prefix + ".species_family_counts.tsv"
    qc_path = args.out_prefix + ".qc.tsv"

    with open(prot_calls_path, "w") as out:
        out.write("species\tprotein\tfamily\tevidence\n")
        for species, prot, fam, ev in sorted(prot_calls):
            out.write(f"{species}\t{prot}\t{fam}\t{ev}\n")

    with open(species_counts_path, "w") as out:
        out.write("species\t" + "\t".join(fam_list) + "\n")
        for s in species_list:
            row = [s]
            for fam in fam_list:
                row.append(str(len(species_family_to_prots.get((s, fam), set()))))
            out.write("\t".join(row) + "\n")

    with open(qc_path, "w") as out:
        out.write("species\tipr_rows\tproteins_in_ipr\tproteins_in_fasta\tassigned_hits\n")
        for s in sorted(qc.keys()):
            d = qc[s]
            out.write(f"{d['species']}\t{d['ipr_rows']}\t{d['proteins_in_ipr']}\t{d['proteins_in_fasta']}\t{d['assigned_hits']}\n")

    print(f"[OK] Wrote:\n  {prot_calls_path}\n  {species_counts_path}\n  {qc_path}")

if __name__ == "__main__":
    main()

