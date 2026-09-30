#!/bin/bash
#SBATCH --job-name=def_scan
set -uo pipefail
source "$DEF_HOME/slurm/config.sh"
load_scan

# submit.sh wrote scan_list.txt, one sample_id per line, before this array
# existed. Task N reads line N, so the index-to-sample mapping is fixed at
# submission and auditable afterwards.
SID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DEF_OUT/scan_list.txt")
[[ -n "$SID" ]] || { echo "no sample at index $SLURM_ARRAY_TASK_ID"; exit 1; }
FAA="$DEF_OUT/proteomes/$SID.faa"
[[ -s "$FAA" ]] || { echo "ERROR: collapsed proteome missing or empty: $FAA"; exit 1; }

echo "=== $SID  task $SLURM_ARRAY_TASK_ID  $(hostname)  $(date) ==="
echo "python3: $(command -v python3)   hmmsearch: $(command -v hmmsearch)"
DB="${DEF_DB}/defensome.hmm"
[[ -s "$DB" ]] || { echo "ERROR: $DB missing; the setup job did not finish"; exit 1; }
ARCH=()
if [[ "${DEF_FULL_ARCH:-0}" == "1" && -n "${DEF_PFAM:-}" ]]; then ARCH=(--full-pfam "$DEF_PFAM"); fi
python3 "$DEF_HOME/defensome.py" scan --proteomes "$DEF_OUT/proteomes" \
        --pfam "$DB" --out "$DEF_OUT" \
        --threads "${SLURM_CPUS_PER_TASK:-8}" --species "$SID" "${ARCH[@]}"
rc=$?
echo "=== $SID done rc=$rc $(date) ==="
exit $rc
