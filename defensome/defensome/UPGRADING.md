# Upgrading from an earlier copy

Earlier versions needed `dashboard.html`, `dashboard.js` and
`defensome_map.tsv` next to `defensome.py`. If you extracted a tarball and also
have loose copies of those files in a parent directory, the loose ones shadow
the new ones and you get errors like:

```
ERROR: template not found: .../tool/dashboard.html
```

From 1.2.0 everything is embedded in `defensome.py`, so the clean fix is:

```bash
cd /path/to/tool
python3 defensome.py --version        # confirm which file you are running
```

If that prints a path you did not expect, or a version below 1.2.0, delete the
stale loose copies and keep only the extracted directory:

```bash
rm -f dashboard.html dashboard.js defensome_map.tsv defensome.py \
      01_scan.sh 02_report.sh 03_deep.sh submit.sh test_dashboard.sh
tar xzf defensome.tar.gz
cd defensome
python3 defensome.py --version
```

Nothing in `results/` is affected. Everything downstream of `hmmsearch` is
recomputed from the cached `domtblout` files in seconds.

---

# Troubleshooting an implausible gene count

`inspect` flags `IMPLAUSIBLE_GENE_COUNT` when a genome annotation still has more
than about 35,000 genes after collapsing. An insect genome should land roughly
between 10,000 and 30,000.

The usual cause is an unfiltered ab initio AUGUSTUS set. Its IDs look like
`AUGUSTUSALYP00000000001.t1`: the `.t1` suffix only groups isoforms when the
part before it is a *gene* ID (`g123.t1`, `g123.t2`). Here it is a unique
protein serial, so nothing groups and every prediction counts as a gene,
including many short spurious ORFs.

Diagnose it:

```bash
F=peps/A_lych_pep_DToL.fasta
grep -c '^>' $F                                   # total entries
head -3 $F                                        # any gene: tag at all?
awk '/^>/{if(l)print l;l=0;next}{l+=length}END{print l}' $F \
  | sort -n | awk '{a[NR]=$1}END{print "median length",a[int(NR/2)]; \
    n=0; for(i in a) if(a[i]<100) n++; print n, "proteins under 100 aa"}'
```

A large spike under 100 aa confirms spurious ORFs. Options, in order of
preference: obtain the filtered annotation from the DToL portal or Ensembl;
restrict to proteins with a BUSCO or Pfam hit; or leave the sample out of the
paired comparison with `compare --exclude A_lych`. Do not simply raise the
threshold. The family filters (`min_len`, `min_cov`) keep most spurious ORFs out
of the defensome counts, but they cannot fix a gene total inflated threefold,
and that total is what any normalisation divides by.
