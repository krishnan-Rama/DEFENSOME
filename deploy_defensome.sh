#!/usr/bin/env bash
set -euo pipefail

# deploy_defensome.sh (rerun end-to-end, safely)
# Chain: 1-prep_nr.slurm -> 2-ipr_array.slurm (array) -> 3-postprocess.slurm
#
# Key fixes vs your original:
# - Ensures required workflow scripts exist (parse/extract) before submitting anything
# - Regenerates defensome_map_pfam.for_parser.tsv from your master mapping (tabbed)
# - Avoids "wrong latest submit_array log" by writing array jobid to a dedicated file
# - Optional clean-run mode to wipe old outputs (set CLEAN=1)
#
# Usage:
#   bash deploy_defensome.sh
# Optional:
#   CLEAN=1 bash deploy_defensome.sh

OUT_DIR=/mnt/ecotox/GROUP-smbpk/c23048124/leps/analysis/defensome_run
SCRIPTS_DIR=/mnt/ecotox/GROUP-smbpk/c23048124/leps/analysis/workflow_scripts

CLEAN="${CLEAN:-0}"

mkdir -p "$OUT_DIR/logs" "$OUT_DIR/workflow"

# ----------------------------
# Preconditions (hard fail)
# ----------------------------
need_file() { [[ -s "$1" ]] || { echo "[ERROR] Missing/empty: $1" >&2; exit 1; }; }

need_file "1-prep_nr.slurm"
need_file "2-ipr_array.slurm"
need_file "3-postprocess.slurm"

# You must have a master mapping in tab-delimited format:
need_file "$OUT_DIR/workflow/defensome_map_pfam.tsv"

# You must have the python scripts somewhere stable:
need_file "$SCRIPTS_DIR/parse_interpro_defensome.py"
need_file "$SCRIPTS_DIR/extract_family_fastas.py"

# Optionally wipe old outputs
if [[ "$CLEAN" == "1" ]]; then
  echo "[WARN] CLEAN=1: removing prior outputs under $OUT_DIR/{00_nr,01_interpro,03_tables,04_fastas,05_alignments,06_trees,tmp}"
  rm -rf "$OUT_DIR"/{00_nr,01_interpro,03_tables,04_fastas,05_alignments,06_trees,tmp}
  mkdir -p "$OUT_DIR"/{00_nr,01_interpro,03_tables,04_fastas,05_alignments,06_trees,tmp}
fi

# ----------------------------
# Build parser-compatible map
# ----------------------------
# Master file columns: category, family, pfam_id, notes, refs
# Parser needs: family, pfam_id
awk 'BEGIN{FS=OFS="\t"}
     /^#/ {next}
     $3 ~ /^PF[0-9]{5}$/ {print $2,$3}
' "$OUT_DIR/workflow/defensome_map_pfam.tsv" > "$OUT_DIR/workflow/defensome_map_pfam.for_parser.tsv"

need_file "$OUT_DIR/workflow/defensome_map_pfam.for_parser.tsv"
echo "[INFO] Wrote: $OUT_DIR/workflow/defensome_map_pfam.for_parser.tsv ($(wc -l < "$OUT_DIR/workflow/defensome_map_pfam.for_parser.tsv") rows)"

# ----------------------------
# Submit prep job
# ----------------------------
prep_jobid=$(sbatch --parsable 1-prep_nr.slurm)
echo "Submitted prep job: $prep_jobid"

# ----------------------------
# Submit array-wrapper: submits InterProScan array after prep completes
# Writes array jobid to a dedicated file so we don't parse "latest log" incorrectly.
# ----------------------------
array_jobid_file="$OUT_DIR/logs/array_jobid.txt"
rm -f "$array_jobid_file"

submit_array_wrapper=$(mktemp)
cat > "$submit_array_wrapper" <<EOF
#!/usr/bin/env bash
set -euo pipefail
OUT_DIR=$OUT_DIR
LIST="\$OUT_DIR/workflow/nr_fastas.list"
JOBID_FILE="$array_jobid_file"

if [[ ! -s "\$LIST" ]]; then
  echo "[ERROR] Missing list: \$LIST" >&2
  exit 1
fi

N=\$(wc -l < "\$LIST")
echo "Array size = \$N"

jid=\$(sbatch --parsable --array=1-"\$N" 2-ipr_array.slurm)
echo "\$jid" | tee "\$JOBID_FILE"
EOF

chmod +x "$submit_array_wrapper"

array_wrapper_jobid=$(sbatch --parsable --dependency=afterok:"$prep_jobid" \
  --job-name=submit_ipr_array \
  --output="$OUT_DIR/logs/submit_array_%j.out" \
  --error="$OUT_DIR/logs/submit_array_%j.err" \
  "$submit_array_wrapper")

echo "Submitted array-wrapper job: $array_wrapper_jobid"

# ----------------------------
# Submit post-wrapper: submits postprocess after array completes
# Uses array_jobid_file written by the array-wrapper.
# ----------------------------
submit_post_wrapper=$(mktemp)
cat > "$submit_post_wrapper" <<EOF
#!/usr/bin/env bash
set -euo pipefail
OUT_DIR=$OUT_DIR
JOBID_FILE="$array_jobid_file"

if [[ ! -s "\$JOBID_FILE" ]]; then
  echo "[ERROR] Missing array jobid file: \$JOBID_FILE" >&2
  exit 1
fi

array_jobid=\$(tr -d '[:space:]' < "\$JOBID_FILE")
if [[ -z "\$array_jobid" ]]; then
  echo "[ERROR] Empty array jobid in \$JOBID_FILE" >&2
  exit 1
fi

echo "Parsed array jobid: \$array_jobid"

post_jobid=\$(sbatch --parsable --dependency=afterok:"\$array_jobid" 3-postprocess.slurm)
echo "Submitted postprocess job: \$post_jobid"
EOF

chmod +x "$submit_post_wrapper"

post_wrapper_jobid=$(sbatch --parsable --dependency=afterok:"$array_wrapper_jobid" \
  --job-name=submit_post \
  --output="$OUT_DIR/logs/submit_post_%j.out" \
  --error="$OUT_DIR/logs/submit_post_%j.err" \
  "$submit_post_wrapper")

echo "Submitted post-wrapper job: $post_wrapper_jobid"

echo "Deployment complete."
echo "Track jobs with: squeue -u $USER"

