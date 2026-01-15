#!/bin/bash
#SBATCH --job-name=def_ipr
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --output=logs/ipr_%A_%a.out
#SBATCH --error=logs/ipr_%A_%a.err

set -euo pipefail

module load interproscan/5.53-87.0

OUT_DIR=/mnt/ecotox/GROUP-smbpk/c23048124/leps/analysis/defensome_run
LIST="$OUT_DIR/workflow/nr_fastas.list"

# You must submit this with the correct --array size (see below)
FA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
BASE=$(basename "$FA" .nr.fa)
OUTTSV="$OUT_DIR/01_interpro/${BASE}.ipr.tsv"

[[ -s "$OUTTSV" ]] && { echo "Exists, skip: $OUTTSV"; exit 0; }

# node-local scratch (big I/O win)
SCRATCH="${TMPDIR:-/tmp}/ipr_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SCRATCH"

# Pfam-only is the key speedup for your Pfam-based defensome mapping
interproscan.sh \
  -i "$FA" \
  -f tsv \
  -appl Pfam \
  -cpu "$SLURM_CPUS_PER_TASK" \
  -T "$SCRATCH" \
  -o "$OUTTSV"
