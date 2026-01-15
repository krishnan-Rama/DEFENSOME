#!/usr/bin/env bash
set -euo pipefail

OUT=defensome_run

module load python/3.10.5 || true

echo "[INFO] Parsing InterPro defensome calls"

python3 $OUT/workflow/parse_interpro_defensome.py \
  --ipr_dir $OUT/01_interpro \
  --nr_dir $OUT/00_nr \
  --pfam_map $OUT/workflow/defensome_map_pfam.for_parser.tsv \
  --out_prefix $OUT/03_tables/defensome

echo "[INFO] Extracting family FASTAs"

python3 $OUT/workflow/extract_family_fastas.py \
  --calls_tsv $OUT/03_tables/defensome.protein_calls.tsv \
  --nr_dir $OUT/00_nr \
  --out_dir $OUT/04_fastas

echo "[DONE]"
