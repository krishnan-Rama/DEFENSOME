#!/bin/bash
#SBATCH --job-name=def_post
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --output=defensome_run/logs/post_%j.out
#SBATCH --error=defensome_run/logs/post_%j.err

set -euo pipefail

module load python/3.10.5 || true
module load mafft/7.505-beg3hnq || true
module load iq-tree/2.2.2.7-2n5zbpa || true

OUT_DIR=/mnt/ecotox/GROUP-smbpk/c23048124/leps/analysis/defensome_run
SCRIPTS=/mnt/ecotox/GROUP-smbpk/c23048124/leps/analysis/workflow_scripts
THREADS=${SLURM_CPUS_PER_TASK}
PYTHON_BIN="$(command -v python3 || command -v python)"

mkdir -p "$OUT_DIR"/{03_tables,04_fastas,05_alignments,06_trees,logs,tmp}

# Preconditions
[[ -s "$SCRIPTS/parse_interpro_defensome.py" ]] || { echo "[ERROR] Missing $SCRIPTS/parse_interpro_defensome.py" >&2; exit 1; }
[[ -s "$SCRIPTS/extract_family_fastas.py" ]] || { echo "[ERROR] Missing $SCRIPTS/extract_family_fastas.py" >&2; exit 1; }
[[ -s "$OUT_DIR/workflow/defensome_map_pfam.for_parser.tsv" ]] || { echo "[ERROR] Missing $OUT_DIR/workflow/defensome_map_pfam.for_parser.tsv" >&2; exit 1; }

# Atomic outputs to avoid header-only truncation
OUTP="$OUT_DIR/03_tables/defensome"
TMPP="$OUT_DIR/tmp/defensome_${SLURM_JOB_ID}"

echo "[INFO] Parsing InterProScan -> defensome calls"
$PYTHON_BIN "$SCRIPTS/parse_interpro_defensome.py" \
  --ipr_dir "$OUT_DIR/01_interpro" \
  --nr_dir "$OUT_DIR/00_nr" \
  --pfam_map "$OUT_DIR/workflow/defensome_map_pfam.for_parser.tsv" \
  --out_prefix "$TMPP" \
  >"$OUT_DIR/logs/parse_interpro_defensome.${SLURM_JOB_ID}.log" 2>&1

# Move into place atomically
mv -f "${TMPP}.protein_calls.tsv" "${OUTP}.protein_calls.tsv"
mv -f "${TMPP}.species_family_counts.tsv" "${OUTP}.species_family_counts.tsv"
mv -f "${TMPP}.qc.tsv" "${OUTP}.qc.tsv"

echo "[INFO] Extracting family FASTAs"
$PYTHON_BIN "$SCRIPTS/extract_family_fastas.py" \
  --calls_tsv "$OUTP.protein_calls.tsv" \
  --nr_dir "$OUT_DIR/00_nr" \
  --out_dir "$OUT_DIR/04_fastas" \
  >"$OUT_DIR/logs/extract_family_fastas.${SLURM_JOB_ID}.log" 2>&1

echo "[INFO] Aligning and building trees (nseq>=4)"
for fa in "$OUT_DIR"/04_fastas/*.fa; do
  fam="$(basename "$fa" .fa)"
  nseq="$(grep -c '^>' "$fa" || true)"
  [[ "$nseq" -lt 4 ]] && continue

  aln="$OUT_DIR/05_alignments/${fam}.aln.fa"
  tree_prefix="$OUT_DIR/06_trees/${fam}"

  [[ -s "$aln" ]] || mafft --thread "$THREADS" --auto "$fa" > "$aln" 2>"$OUT_DIR/logs/${fam}.mafft.log"
  [[ -s "${tree_prefix}.treefile" ]] || iqtree2 -s "$aln" -m MFP -B 1000 -T "$THREADS" --prefix "$tree_prefix" \
    >"$OUT_DIR/logs/${fam}.iqtree.log" 2>&1 || true
done

echo "[DONE] Postprocess complete."
