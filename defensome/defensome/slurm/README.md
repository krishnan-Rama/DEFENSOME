# SLURM

| file | what it is |
|---|---|
| `config.sh` | **edit this**: partitions, account, resources, module loading |
| `submit.sh` | the one command; preflight, resume logic, dependency chain |
| `00_collapse.sh` | isoform collapse, one job |
| `01_scan.sh` | hmmsearch, array job, one task per proteome |
| `02_analysis.sh` | annotate through dashboard, one job |

```bash
bash slurm/submit.sh --samplesheet results/samplesheet.tsv \
     --pfam /db/Pfam-A.hmm --out results --drop-utrorf --exclude A_lych
```

Add `--dry-run` first to see the preflight and the scan list without
submitting. Rerun the same command to resume.
