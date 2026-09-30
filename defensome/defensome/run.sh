#!/usr/bin/env bash
# bash run.sh <proteome_dir> <Pfam-A.hmm> <outdir> [threads] [metadata.tsv] [group_column]
set -euo pipefail
usage() { sed -n '2p' "$0"; exit 1; }
[[ $# -ge 3 ]] || usage
PROT=$1; PFAM=$2; OUT=$3; THREADS=${4:-8}; META=${5:-}; GROUP=${6:-}
HERE=$(cd "$(dirname "$0")" && pwd)

# Catch unset shell variables: "$PROJ/proteomes" with PROJ unset becomes
# "/proteomes", which fails much later with a confusing message.
# if-statements, not `[[ ]] && cmd`: under `set -e` the latter exits the script
# when the condition is false.
for v in PROT PFAM; do
  if [[ "${!v}" == /* && ! -e "${!v}" ]]; then
    echo "ERROR: $v does not exist: ${!v}"
    if [[ "${!v}" =~ ^/[^/]+$ ]]; then echo "  (looks like an unset shell variable)"; fi
    exit 1
  fi
done
[[ -d "$PROT" ]] || { echo "ERROR: not a directory: $PROT"; exit 1; }
[[ -f "$PFAM" ]] || { echo "ERROR: not a file: $PFAM"; exit 1; }
command -v hmmsearch >/dev/null || { echo "ERROR: hmmsearch not found. module load HMMER/3.4-gompi-2023b"; exit 1; }
python3 -c 'import pandas, numpy' 2>/dev/null || { echo "ERROR: pandas/numpy missing. module load SciPy-bundle"; exit 1; }

python3 "$HERE/defensome.py" scan     --proteomes "$PROT" --pfam "$PFAM" --out "$OUT" --threads "$THREADS"
python3 "$HERE/defensome.py" annotate --proteomes "$PROT" --out "$OUT"
if [[ -n "$META" && -n "$GROUP" ]]; then
  python3 "$HERE/defensome.py" report --out "$OUT" --metadata "$META" --group-by "$GROUP"
else
  python3 "$HERE/defensome.py" report --out "$OUT"
fi
echo "Done. Results in $OUT"
