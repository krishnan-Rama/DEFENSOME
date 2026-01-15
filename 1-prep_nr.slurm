#!/bin/bash
#SBATCH --job-name=def_prep
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --output=logs/prep_%j.out
#SBATCH --error=logs/prep_%j.err

set -euo pipefail

module load cd-hit/V4.8.1 || module load cd-hit/V4.6.8

PEPS_DIR=/mnt/ecotox/GROUP-smbpk/c23048124/leps/peps
OUT_DIR=/mnt/ecotox/GROUP-smbpk/c23048124/leps/analysis/defensome_run
THREADS=${SLURM_CPUS_PER_TASK}

mkdir -p "$OUT_DIR"/{00_input,00_nr,logs,workflow,01_interpro}

# link inputs
find "$PEPS_DIR" -maxdepth 1 -type f -name "*.fa" -print0 | sort -z | while IFS= read -r -d '' f; do
  ln -sf "$f" "$OUT_DIR/00_input/$(basename "$f")"
done

# cd-hit per proteome
for fa in "$OUT_DIR"/00_input/*.fa; do
  base="$(basename "$fa" .fa)"
  outfa="$OUT_DIR/00_nr/${base}.nr.fa"
  [[ -s "$outfa" ]] && continue
  cd-hit -i "$fa" -o "$outfa" -c 0.99 -n 5 -d 0 -M 0 -T "$THREADS" \
    >"$OUT_DIR/logs/${base}.cdhit.log" 2>&1
done

# list for array
ls -1 "$OUT_DIR"/00_nr/*.nr.fa > "$OUT_DIR/workflow/nr_fastas.list"
echo "Wrote: $OUT_DIR/workflow/nr_fastas.list"
