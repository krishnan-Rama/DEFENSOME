#!/bin/bash
#SBATCH --job-name=def_analysis
set -uo pipefail
source "$DEF_HOME/slurm/config.sh"
load_analysis

D="$DEF_HOME/defensome.py"
T="${SLURM_CPUS_PER_TASK:-16}"
P="$DEF_OUT/proteomes"
FAILED=()

# must: a link in the core chain; a failure stops the job.
must() {
  echo; echo "########## $1  $(date +%H:%M:%S) ##########"
  if "${@:2}"; then echo "---------- $1 OK"; else
    rc=$?; echo "---------- $1 FAILED (rc=$rc)"
    echo "STOPPING: every later step depends on '$1'."; exit "$rc"; fi
}
# step: an independent leaf; a failure is recorded and the job continues.
step() {
  echo; echo "########## $1  $(date +%H:%M:%S) ##########"
  if "${@:2}"; then echo "---------- $1 OK"; else
    echo "---------- $1 FAILED (rc=$?)"; FAILED+=("$1"); fi
}

# DEF_STEPS (colon-separated, from submit.sh --steps) limits which steps run
want() { [[ -z "${DEF_STEPS:-}" || ":${DEF_STEPS}:" == *":$1:"* ]]; }

echo "=== analysis $(date) on $(hostname) ==="
echo "steps: ${DEF_STEPS:-all}"
python3 "$D" doctor

want annotate && must annotate python3 "$D" annotate --proteomes "$P" --out "$DEF_OUT"
if want short-orfs && python3 -c 'import csv,sys; r=list(csv.DictReader(open(sys.argv[1]),delimiter="\t")); sys.exit(0 if any((x.get("nucleotide") or "").strip() for x in r) else 1)' "$DEF_SHEET"; then
  step short-orfs python3 "$D" short-orfs --samplesheet "$DEF_SHEET" --out "$DEF_OUT"
fi
want qc       && step qc       python3 "$D" qc       --out "$DEF_OUT"
want extract  && step extract  python3 "$D" extract  --out "$DEF_OUT" --proteomes "$P"
want domains  && step domains  python3 "$D" domains  --out "$DEF_OUT"

if want cyp; then
  if [[ -n "${DEF_CYP_REFS:-}" && -f "${DEF_CYP_REFS}" ]]; then
    [[ -s "$DEF_OUT/fasta/CYP.faa" ]] || step extract python3 "$D" extract --out "$DEF_OUT" --proteomes "$P"
    step cyp python3 "$D" cyp --out "$DEF_OUT" --refs "$DEF_CYP_REFS" --family CYP --threads "$T"
  else
    echo; echo "########## cyp ##########"
    echo "WARNING: no CYP clan reference, so CYP clans are NOT assigned."
    echo "         The dashboard clan ring and Clans tab will be empty."
    FAILED+=("cyp(no reference)")
  fi
fi

if want trees && python3 "$D" doctor --require trees >/dev/null 2>&1; then
  step trees python3 "$D" trees --out "$DEF_OUT" --threads "$T" \
        --families "${DEF_TREE_FAMILIES:-CYP,CCE,GST,UGT,ABC_full}"
elif want trees; then
  echo; echo "skipping trees: MAFFT/FastTree not found (see doctor output above)"
fi

# publication PDFs of the CYP gene tree, coloured by clan and by method
if want genetree && [[ -s "$DEF_OUT/trees/CYP.tre" ]]; then
  [[ -s "$DEF_OUT/cyp/CYP_clan_calls.tsv" ]] && \
    step genetree_clan python3 "$D" genetree --out "$DEF_OUT" --family CYP --by clan
  step genetree_method python3 "$D" genetree --out "$DEF_OUT" --family CYP \
        --samplesheet "$DEF_SHEET"
fi

# The method comparison. --complete-only strips Trinity fragments; raw counts
# are correct for a within-species paired design.
want compare && step compare python3 "$D" compare --samplesheet "$DEF_SHEET" --out "$DEF_OUT" \
      --reference "${DEF_REFERENCE:-DToL}" --complete-only \
      ${DEF_EXCLUDE:+--exclude "$DEF_EXCLUDE"}

if ! want report; then :
elif [[ -n "${DEF_META:-}" && -n "${DEF_GROUP:-}" ]]; then
  step report python3 "$D" report --out "$DEF_OUT" --metadata "$DEF_META" --group-by "$DEF_GROUP"
else
  step report python3 "$D" report --out "$DEF_OUT"
fi

if want tree && [[ -n "${DEF_SPECIES_TREE:-}" && -f "${DEF_SPECIES_TREE}" ]]; then
  step tree python3 "$D" tree --out "$DEF_OUT" --tree "$DEF_SPECIES_TREE" \
        ${DEF_META:+--metadata "$DEF_META"} ${DEF_GROUP:+--group-by "$DEF_GROUP"}
fi

want dashboard && step dashboard python3 "$D" dashboard --out "$DEF_OUT" \
      --samplesheet "$DEF_SHEET" \
      ${DEF_META:+--metadata "$DEF_META"} ${DEF_GROUP:+--group-by "$DEF_GROUP"} \
      ${DEF_SPECIES_TREE:+--species-tree "$DEF_SPECIES_TREE"}

echo; echo "########## SUMMARY $(date) ##########"
[[ ${#FAILED[@]} -eq 0 ]] && echo "all steps OK" || echo "FAILED: ${FAILED[*]}"
echo "results:    $DEF_OUT"
echo "comparison: $DEF_OUT/comparison/"
echo "dashboard:  $DEF_OUT/dashboard.html"
[[ ${#FAILED[@]} -eq 0 ]]
