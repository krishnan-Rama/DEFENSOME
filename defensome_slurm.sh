#!/bin/bash
#SBATCH --job-name=ipr
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --array=1-$(wc -l < proteomes.list)
#SBATCH --output=logs/ipr_%A_%a.out
#SBATCH --error=logs/ipr_%A_%a.err

module load interproscan/5.53-87.0

FA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" proteomes.list)
BASE=$(basename "$FA" .fa)

OUTDIR=defensome_run/01_interpro
mkdir -p "$OUTDIR" logs

# Use local scratch if available (big speedup on many clusters)
SCRATCH=${TMPDIR:-/tmp}
mkdir -p "$SCRATCH/ipr_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

interproscan.sh \
  -i "$FA" \
  -f tsv \
  -cpu "$SLURM_CPUS_PER_TASK" \
  -o "$OUTDIR/${BASE}.ipr.tsv" \
  -T "$SCRATCH/ipr_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

