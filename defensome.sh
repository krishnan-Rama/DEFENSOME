#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# defensome_pipeline.sh
#
# End-to-end comparative "defensome" pipeline for insect proteomes:
#  0) Organize inputs
#  1) Reduce redundancy (CD-HIT) (optional but strongly recommended)
#  2) Annotate domains (InterProScan) (primary)
#  3) (Optional) Targeted HMMER scan if you provide defensome.hmm
#  4) Build:
#     - protein->family assignments with evidence
#     - species x family count matrix
#     - QC statistics per proteome
#  5) Extract per-family FASTAs
#  6) Align (MAFFT) + Tree (IQ-TREE) per family (optional but included)
#
# REQUIREMENTS (modules you said you have):
#   interproscan, cd-hit, python, mafft, iq-tree
#
# USAGE:
#   bash defensome_pipeline.sh \
#     /mnt/ecotox/GROUP-smbpk/c23048124/leps/peps \
#     /mnt/ecotox/GROUP-smbpk/c23048124/leps/analysis/defensome_run
#
# OPTIONAL ENV VARS:
#   THREADS=16
#   CDHIT=1|0               (default 1)
#   IPR=1|0                 (default 1)
#   HMMER=1|0               (default 0; requires DEFENSOME_HMM)
#   DEFENSOME_HMM=/path/to/defensome.hmm
#   TREES=1|0               (default 1)
#   MIN_LEN_CYP=350         (default 350)
#   MIN_LEN_ABC_NBD=150     (default 150)
#   MIN_DOM_COV=0.35        (default 0.35)  # domain coverage threshold in parser
#
# NOTES:
# - The PFAM/IPR mapping below is intentionally conservative and incomplete.
#   You MUST edit mapping tables for your exact "defensome" definition.
# - This pipeline does NOT do dN/dS or selection tests (needs CDS/orthology).
###############################################################################

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <proteomes_dir> <out_dir>"
  exit 1
fi

PEPS_DIR="$(readlink -f "$1")"
OUT_DIR="$(readlink -f "$2")"

THREADS="${THREADS:-16}"
CDHIT="${CDHIT:-1}"
IPR="${IPR:-1}"
HMMER="${HMMER:-0}"
TREES="${TREES:-1}"

MIN_LEN_CYP="${MIN_LEN_CYP:-350}"
MIN_LEN_ABC_NBD="${MIN_LEN_ABC_NBD:-150}"
MIN_DOM_COV="${MIN_DOM_COV:-0.35}"

mkdir -p "$OUT_DIR"/{00_input,00_nr,01_interpro,02_hmmer,03_tables,04_fastas,05_alignments,06_trees,workflow,logs,tmp}

# Copy (or symlink) inputs
echo "[INFO] Collecting input proteomes from: $PEPS_DIR"
find "$PEPS_DIR" -maxdepth 1 -type f -name "*.fa" -print0 | sort -z | while IFS= read -r -d '' f; do
  ln -sf "$f" "$OUT_DIR/00_input/$(basename "$f")"
done

# Helper: module load safely
ml() {
  local m="$1"
  if command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    module load "$m" >/dev/null 2>&1 || true
  fi
}

# Load likely modules (silently continue if module name differs)
ml "cd-hit/V4.8.1"
ml "cd-hit/V4.6.8"
ml "interproscan/5.53-87.0"
ml "python/3.10.5"
ml "python/3.11.6-63oqiza"
ml "python/3.12.1-yidmbj6"
ml "mafft/7.505-beg3hnq"
ml "iq-tree/2.2.2.7-2n5zbpa"
ml "HMMER/3.4-gompi-2023a"
ml "hmmer/3.1b2"

###############################################################################
# Step 0: basic sanity checks
###############################################################################
echo "[INFO] Checking tools..."
need() { command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Missing executable: $1"; exit 1; }; }

if [[ "$CDHIT" == "1" ]]; then
  need cd-hit
fi

if [[ "$IPR" == "1" ]]; then
  need interproscan.sh
fi

need python3 || need python

if [[ "$TREES" == "1" ]]; then
  # optional, only if present
  command -v mafft >/dev/null 2>&1 || echo "[WARN] mafft not found; alignments/trees may fail."
  command -v iqtree2 >/dev/null 2>&1 || echo "[WARN] iqtree2 not found; trees may fail."
fi

if [[ "$HMMER" == "1" ]]; then
  need hmmscan
  : "${DEFENSOME_HMM:?Set DEFENSOME_HMM=/path/to/defensome.hmm when HMMER=1}"
  [[ -s "$DEFENSOME_HMM" ]] || { echo "[ERROR] DEFENSOME_HMM not found or empty: $DEFENSOME_HMM"; exit 1; }
fi

###############################################################################
# Step 1: redundancy reduction (CD-HIT)
###############################################################################
if [[ "$CDHIT" == "1" ]]; then
  echo "[INFO] Running CD-HIT to reduce redundancy (0.99 identity)."
  for fa in "$OUT_DIR"/00_input/*.fa; do
    base="$(basename "$fa" .fa)"
    outfa="$OUT_DIR/00_nr/${base}.nr.fa"
    if [[ -s "$outfa" ]]; then
      echo "[INFO] CD-HIT output exists, skipping: $outfa"
      continue
    fi
    cd-hit -i "$fa" -o "$outfa" -c 0.99 -n 5 -d 0 -M 0 -T "$THREADS" \
      >"$OUT_DIR/logs/${base}.cdhit.log" 2>&1
  done
else
  echo "[INFO] CD-HIT disabled; using raw proteomes."
  for fa in "$OUT_DIR"/00_input/*.fa; do
    base="$(basename "$fa" .fa)"
    ln -sf "$fa" "$OUT_DIR/00_nr/${base}.nr.fa"
  done
fi

###############################################################################
# Step 2: InterProScan per proteome
###############################################################################
if [[ "$IPR" == "1" ]]; then
  echo "[INFO] Running InterProScan (TSV + GFF3)."
  for fa in "$OUT_DIR"/00_nr/*.nr.fa; do
    base="$(basename "$fa" .nr.fa)"
    outtsv="$OUT_DIR/01_interpro/${base}.ipr.tsv"
    if [[ -s "$outtsv" ]]; then
      echo "[INFO] InterProScan TSV exists, skipping: $outtsv"
      continue
    fi
    # -f tsv,gff3 can sometimes generate multiple outputs; we request TSV as primary.
    interproscan.sh \
      -i "$fa" \
      -f tsv \
      -dp \
      -goterms \
      -iprlookup \
      -pa \
      -cpu "$THREADS" \
      -o "$outtsv" \
      >"$OUT_DIR/logs/${base}.interpro.log" 2>&1
  done
else
  echo "[WARN] InterProScan disabled. Parser expects InterProScan TSVs unless you implement HMMER parsing."
fi

###############################################################################
# Step 3 (optional): Targeted HMMER scan (requires DEFENSOME_HMM)
###############################################################################
if [[ "$HMMER" == "1" ]]; then
  echo "[INFO] Running targeted HMMER scans using: $DEFENSOME_HMM"
  # hmmpress only needed for HMM databases; safe to run always
  if [[ ! -s "${DEFENSOME_HMM}.h3p" ]]; then
    hmmpress "$DEFENSOME_HMM" >"$OUT_DIR/logs/hmmpress.log" 2>&1 || true
  fi
  for fa in "$OUT_DIR"/00_nr/*.nr.fa; do
    base="$(basename "$fa" .nr.fa)"
    domtbl="$OUT_DIR/02_hmmer/${base}.domtblout"
    if [[ -s "$domtbl" ]]; then
      echo "[INFO] HMMER domtblout exists, skipping: $domtbl"
      continue
    fi
    hmmscan --cpu "$THREADS" --domtblout "$domtbl" \
      "$DEFENSOME_HMM" "$fa" \
      >"$OUT_DIR/logs/${base}.hmmscan.txt" 2>&1
  done
fi

###############################################################################
# Step 4: Build mapping tables + parse InterProScan TSV to defensome calls
###############################################################################
echo "[INFO] Writing defensome mapping tables (edit these!)."

cat > "$OUT_DIR/workflow/defensome_map_pfam.tsv" << 'EOF'
#family	pfam_id	notes
CYP	PF00067	Cytochrome P450
GST_N	PF02798	GST N-term
GST_C	PF00043	GST C-term
UGT	PF00201	UDP-glucuronosyltransferase
SULT	PF00685	Sulfotransferase
ABC_NBD	PF00005	ABC transporter ATP-binding
MFS	PF07690	Major facilitator superfamily (optional)
AKR	PF00248	Aldo-keto reductase
Catalase	PF00199	Catalase
SOD_CuZn	PF00080	SOD (Cu/Zn)
SOD_MnFe	PF00081	SOD (Mn/Fe)
Thioredoxin	PF00085	Thioredoxin
Peroxiredoxin	PF00578	Peroxiredoxin
Carboxylesterase	PF00135	Esterase/lipase (broad; needs curation)
Epoxide_hydrolase	PF00561	Epoxide hydrolase (alpha/beta hydrolase fold; broad)
FMO	PF00743	Flavin-containing monooxygenase (FMO-like)
EOF

# Optional InterPro IDs if you prefer them (you can extend this)
cat > "$OUT_DIR/workflow/defensome_map_interpro.tsv" << 'EOF'
#family	interpro_id	notes
CYP	IPR001128	Cytochrome P450
GST	IPR004045	Glutathione S-transferase (general)
UGT	IPR002213	UDP-glucuronosyl/UDP-glucosyltransferase
SULT	IPR000863	Sulfotransferase
ABC_NBD	IPR003439	ABC transporter-like
EOF

echo "[INFO] Parsing InterProScan TSV into defensome assignments + count matrix."

cat > "$OUT_DIR/workflow/parse_interpro_defensome.py" << 'PY'
import argparse
import os
import re
from collections import defaultdict

def parse_fasta_lengths(fa_path):
    lengths = {}
    with open(fa_path) as f:
        seq_id = None
        seq = []
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

def load_map_tsv(path):
    fam_to_ids = defaultdict(set)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fam, acc, *_ = line.split("\t")
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

    # InterProScan TSV columns (typical):
    # protein_accession, md5, length, analysis, signature_accession, signature_desc,
    # start, end, score, status, date, interpro_accession, interpro_desc, GO, pathways
    #
    # We will use: protein_accession, length, signature_accession, start, end, interpro_accession
    ipr_files = sorted([p for p in os.listdir(args.ipr_dir) if p.endswith(".ipr.tsv")])

    # Outputs
    prot_calls_path = args.out_prefix + ".protein_calls.tsv"
    species_counts_path = args.out_prefix + ".species_family_counts.tsv"
    qc_path = args.out_prefix + ".qc.tsv"

    # Collect per species/protein best evidence
    # calls[(species, protein)][family] = list of evidences
    calls = defaultdict(lambda: defaultdict(list))
    qc = {}

    # Preload fasta lengths to validate protein lengths
    fasta_len = {}
    for fn in sorted([p for p in os.listdir(args.nr_dir) if p.endswith(".nr.fa")]):
        species = fn.replace(".nr.fa", "")
        fasta_len[species] = parse_fasta_lengths(os.path.join(args.nr_dir, fn))

    for fn in ipr_files:
        species = fn.replace(".ipr.tsv", "")
        path = os.path.join(args.ipr_dir, fn)
        n_rows = 0
        n_prots_seen = set()
        n_hits = 0
        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                prot = parts[0]
                try:
                    ipr_len = int(parts[2])
                except Exception:
                    ipr_len = None
                analysis = parts[3]
                sig_acc = parts[4] if len(parts) > 4 else ""
                start = int(parts[6]) if len(parts) > 6 and parts[6].isdigit() else None
                end = int(parts[7]) if len(parts) > 7 and parts[7].isdigit() else None
                ipr_acc = parts[11] if len(parts) > 11 else ""

                n_rows += 1
                n_prots_seen.add(prot)

                # Domain coverage filter
                if start is None or end is None or ipr_len is None or ipr_len == 0:
                    dom_cov = None
                else:
                    dom_cov = (end - start + 1) / float(ipr_len)

                if dom_cov is not None and dom_cov < args.min_dom_cov:
                    continue

                # Assign families based on PFAM signature accession (preferred here)
                for fam, ids in pfam_map.items():
                    if sig_acc in ids:
                        # additional length sanity filters for known problematic families
                        prot_len = fasta_len.get(species, {}).get(prot, ipr_len or 0)

                        if fam == "CYP" and prot_len < args.min_len_cyp:
                            continue
                        if fam == "ABC_NBD" and prot_len < args.min_len_abc_nbd:
                            continue

                        ev = f"{analysis}:{sig_acc}:{start}-{end}:cov={dom_cov:.3f}" if dom_cov is not None else f"{analysis}:{sig_acc}:{start}-{end}"
                        calls[(species, prot)][fam].append(ev)
                        n_hits += 1

        # QC
        total_prots = len(fasta_len.get(species, {}))
        qc[species] = {
            "species": species,
            "ipr_rows": n_rows,
            "proteins_in_ipr": len(n_prots_seen),
            "proteins_in_fasta": total_prots,
            "assigned_hits": n_hits,
        }

    # Post-process for composite definitions:
    # - GST: require both GST_N and GST_C for "GST_strict"
    # - Keep also raw GST_N/GST_C counts for transparency
    prot_calls = []
    for (species, prot), fams in calls.items():
        famset = set(fams.keys())
        if "GST_N" in famset and "GST_C" in famset:
            fams["GST_strict"] = fams.get("GST_N", []) + fams.get("GST_C", [])
        # Optionally: CYP is direct (PF00067)
        for fam, evs in fams.items():
            prot_calls.append((species, prot, fam, ";".join(evs)))

    # Count matrix: species x family, counting unique proteins assigned to family
    species_family_to_prots = defaultdict(set)
    for species, prot, fam, _ in prot_calls:
        species_family_to_prots[(species, fam)].add(prot)

    species_list = sorted({s for s, _ in species_family_to_prots.keys()} | set(qc.keys()))
    fam_list = sorted({f for _, f in species_family_to_prots.keys()})

    # Write protein calls
    with open(prot_calls_path, "w") as out:
        out.write("species\tprotein\tfamily\tevidence\n")
        for species, prot, fam, ev in sorted(prot_calls):
            out.write(f"{species}\t{prot}\t{fam}\t{ev}\n")

    # Write count matrix
    with open(species_counts_path, "w") as out:
        out.write("species\t" + "\t".join(fam_list) + "\n")
        for s in species_list:
            row = [s]
            for fam in fam_list:
                row.append(str(len(species_family_to_prots.get((s, fam), set()))))
            out.write("\t".join(row) + "\n")

    # Write QC
    with open(qc_path, "w") as out:
        out.write("species\tipr_rows\tproteins_in_ipr\tproteins_in_fasta\tassigned_hits\n")
        for s in sorted(qc.keys()):
            d = qc[s]
            out.write(f"{d['species']}\t{d['ipr_rows']}\t{d['proteins_in_ipr']}\t{d['proteins_in_fasta']}\t{d['assigned_hits']}\n")

    print(f"[OK] Wrote:\n  {prot_calls_path}\n  {species_counts_path}\n  {qc_path}")

if __name__ == "__main__":
    main()
PY

PYTHON_BIN="$(command -v python3 || command -v python)"

"$PYTHON_BIN" "$OUT_DIR/workflow/parse_interpro_defensome.py" \
  --ipr_dir "$OUT_DIR/01_interpro" \
  --nr_dir "$OUT_DIR/00_nr" \
  --pfam_map "$OUT_DIR/workflow/defensome_map_pfam.tsv" \
  --out_prefix "$OUT_DIR/03_tables/defensome" \
  --min_dom_cov "$MIN_DOM_COV" \
  --min_len_cyp "$MIN_LEN_CYP" \
  --min_len_abc_nbd "$MIN_LEN_ABC_NBD" \
  >"$OUT_DIR/logs/parse_interpro_defensome.log" 2>&1

###############################################################################
# Step 5: Extract per-family FASTAs from NR proteomes (based on protein_calls.tsv)
###############################################################################
echo "[INFO] Extracting per-family FASTAs."

cat > "$OUT_DIR/workflow/extract_family_fastas.py" << 'PY'
import argparse
import os
from collections import defaultdict

def read_fasta(path):
    seqs = {}
    with open(path) as f:
        sid = None
        buf = []
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

    # species->family->set(prot)
    want = defaultdict(lambda: defaultdict(set))
    families = set()

    with open(args.calls_tsv) as f:
        header = f.readline()
        for line in f:
            species, prot, fam, _ = line.rstrip("\n").split("\t", 3)
            want[species][fam].add(prot)
            families.add(fam)

    os.makedirs(args.out_dir, exist_ok=True)

    # Initialize output handles per family
    out_handles = {}
    for fam in sorted(families):
        out_handles[fam] = open(os.path.join(args.out_dir, f"{fam}.fa"), "w")

    for fn in sorted([p for p in os.listdir(args.nr_dir) if p.endswith(".nr.fa")]):
        species = fn.replace(".nr.fa", "")
        fasta = read_fasta(os.path.join(args.nr_dir, fn))
        for fam, prots in want.get(species, {}).items():
            oh = out_handles[fam]
            for p in sorted(prots):
                if p in fasta:
                    oh.write(f">{species}|{p}\n")
                    seq = fasta[p]
                    for i in range(0, len(seq), 60):
                        oh.write(seq[i:i+60] + "\n")

    for h in out_handles.values():
        h.close()

    print(f"[OK] Wrote family FASTAs to: {args.out_dir}")

if __name__ == "__main__":
    main()
PY

"$PYTHON_BIN" "$OUT_DIR/workflow/extract_family_fastas.py" \
  --calls_tsv "$OUT_DIR/03_tables/defensome.protein_calls.tsv" \
  --nr_dir "$OUT_DIR/00_nr" \
  --out_dir "$OUT_DIR/04_fastas" \
  >"$OUT_DIR/logs/extract_family_fastas.log" 2>&1

###############################################################################
# Step 6: Align + tree per family (optional)
###############################################################################
if [[ "$TREES" == "1" ]]; then
  echo "[INFO] Aligning and building trees per family (MAFFT + IQ-TREE if available)."
  mkdir -p "$OUT_DIR/05_alignments" "$OUT_DIR/06_trees"

  for fa in "$OUT_DIR"/04_fastas/*.fa; do
    fam="$(basename "$fa" .fa)"
    # Skip tiny families
    nseq="$(grep -c '^>' "$fa" || true)"
    if [[ "$nseq" -lt 4 ]]; then
      echo "[INFO] Skipping $fam (nseq=$nseq < 4)."
      continue
    fi

    aln="$OUT_DIR/05_alignments/${fam}.aln.fa"
    tree_prefix="$OUT_DIR/06_trees/${fam}"

    if command -v mafft >/dev/null 2>&1; then
      if [[ ! -s "$aln" ]]; then
        mafft --thread "$THREADS" --auto "$fa" > "$aln" 2>"$OUT_DIR/logs/${fam}.mafft.log"
      fi
    else
      echo "[WARN] mafft not found; cannot align $fam"
      continue
    fi

    if command -v iqtree2 >/dev/null 2>&1; then
      if [[ ! -s "${tree_prefix}.treefile" ]]; then
        iqtree2 -s "$aln" -m MFP -B 1000 -T "$THREADS" --prefix "$tree_prefix" \
          >"$OUT_DIR/logs/${fam}.iqtree.log" 2>&1 || echo "[WARN] IQ-TREE failed for $fam (see log)."
      fi
    else
      echo "[WARN] iqtree2 not found; skipping tree for $fam"
    fi
  done
fi

###############################################################################
# Step 7: Summary report stub
###############################################################################
cat > "$OUT_DIR/03_tables/README_next_steps.txt" << 'EOF'
Key outputs:
- 03_tables/defensome.species_family_counts.tsv
- 03_tables/defensome.protein_calls.tsv
- 03_tables/defensome.qc.tsv
- 04_fastas/*.fa (per-family sequences)
- 05_alignments/*.aln.fa (per-family alignments)
- 06_trees/* (per-family IQ-TREE outputs)

Critical next steps you MUST do to make this biologically defensible:
1) Edit workflow/defensome_map_pfam.tsv to include all families you mean by "defensome".
   The provided table is conservative and incomplete.
2) Add stricter per-family rules where necessary:
   - ABC transporters: require NBD + TMD domains (current parser only ensures NBD).
   - Carboxylesterases and epoxide hydrolases: PF families used here are broad and will overcall.
3) Validate CYP calls:
   - length filter is on; but you should motif-check (heme-binding region) and/or use CYP reference placement
     if you need clan/family classification.
4) Consider isoform inflation: CD-HIT at 0.99 is a pragmatic step, not true longest-isoform per locus.
EOF

echo "[DONE] Pipeline finished. See $OUT_DIR/03_tables/defensome.species_family_counts.tsv"

