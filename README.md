# Pipeline to retrieve multiple species defensome ancestry. (Not in production, yet)

# Lepidoptera Defensome Comparative Genomics Pipeline

## Overview

This pipeline performs systematic identification, quantification, and phylogenetic analysis of detoxification gene families ("defensome") across 42 Lepidoptera species plus outgroups. The defensome comprises genes involved in xenobiotic metabolism and oxidative stress response, which are central to insect adaptation to host plant toxins and synthetic insecticides.

---

## Biological Rationale

### The Defensome Concept

Insects encounter diverse toxic compounds from host plants (alkaloids, glucosinolates, furanocoumarins) and anthropogenic sources (insecticides). The "defensome" refers to the integrated network of genes that detect, metabolize, and excrete these xenobiotics:

| Phase | Function | Gene Families |
|-------|----------|---------------|
| **Phase 0** | Transcriptional regulation of detox response | Nuclear receptors (HR96), bHLH-PAS TFs, Keap1/CncC pathway |
| **Phase I** | Functionalization (oxidation, reduction, hydrolysis) | CYP450s, Carboxylesterases (CCE), Flavin monooxygenases (FMO), Epoxide hydrolases (EH) |
| **Phase II** | Conjugation (increase solubility) | GSTs, UGTs, SULTs |
| **Phase III** | Excretion (efflux transport) | ABC transporters |
| **Oxidative Stress** | ROS detoxification | SOD, Catalase, Peroxiredoxins |
| **Redox Support** | Maintain cellular redox balance | Thioredoxins |

### Why Lepidoptera?

Lepidoptera (butterflies and moths) represent an ideal system for studying defensome evolution:

1. Extreme dietary diversity (monophagous to polyphagous)
2. Well-documented host plant associations with known toxins
3. Major agricultural pests with insecticide resistance
4. Growing genomic resources from Darwin Tree of Life and similar projects

---

## Species Dataset

**42 species total:**

- **40 Lepidoptera** spanning major families:
  - Noctuidae (owlet moths): 12 species
  - Geometridae (geometer moths): 4 species
  - Nymphalidae (brush-footed butterflies): 5 species
  - Pieridae (whites and sulphurs): 4 species
  - Lycaenidae (blues): 4 species
  - Hesperiidae (skippers): 2 species
  - Sphingidae (hawkmoths): 1 species (*Manduca sexta*)
  - Plutellidae: 1 species (*Plutella xylostella*)
  - Erebidae: 1 species

- **2 Outgroups:**
  - *Drosophila melanogaster* (Diptera) - well-annotated reference
  - *Bombus terrestris* (Hymenoptera) - represents Holometabola

**Proteome sources:** Predicted protein sequences from genome annotations (primarily Darwin Tree of Life assemblies)

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        MODULE 1                                  │
│              Per-species Pfam annotation                         │
│         (SLURM array job: 42 parallel tasks)                     │
├─────────────────────────────────────────────────────────────────┤
│  Input: Species proteome (.fa)                                   │
│                    ↓                                             │
│  hmmscan vs Pfam-A.hmm (--domE 1e-5)                            │
│                    ↓                                             │
│  Filter: i-evalue ≤ 1e-5, HMM coverage ≥ 35%, query cov ≥ 35%   │
│                    ↓                                             │
│  Map Pfam hits → Defensome families (with co-occurrence rules)  │
│                    ↓                                             │
│  Extract CYP candidates (PF00067, length ≥ 250 aa)              │
│                    ↓                                             │
│  Output: .defensome_gene_calls.tsv, .cyp_candidates.faa         │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                        MODULE 2                                  │
│                Cross-species aggregation                         │
├─────────────────────────────────────────────────────────────────┤
│  Aggregate defensome counts across all species                   │
│                    ↓                                             │
│  Build presence/absence matrix                                   │
│                    ↓                                             │
│  Merge all CYP candidates → all_lep_cyps.faa                    │
│                    ↓                                             │
│  Assign CYP clan membership via HMM (FlyBase Dmel references)   │
│     - Build per-clan HMMs from D. melanogaster CYPs             │
│     - hmmscan all Lep CYPs against clan HMMs                    │
│     - Best-hit assignment with QC (bits delta threshold)        │
│                    ↓                                             │
│  Generate summary statistics and figures                         │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                        MODULE 3                                  │
│               Per-clan CYP phylogenetics                         │
├─────────────────────────────────────────────────────────────────┤
│  For each CYP clan (CYP2, CYP3, CYP4, MITO):                    │
│     - Extract assigned sequences                                 │
│     - MAFFT alignment (--auto)                                  │
│     - IQ-TREE2 ML tree (ModelFinder, 1000 UFBoot)               │
│                    ↓                                             │
│  Output: Per-clan alignments and phylogenies                     │
└─────────────────────────────────────────────────────────────────┘
```

---

## Detailed Methods

### Module 1: Pfam Annotation and Defensome Mapping

**Tool:** HMMER 3.x hmmscan

**Database:** Pfam-A.hmm (current release)

**Parameters:**
- Domain E-value threshold: 1e-5 (`--domE`, `--incdomE`)
- Post-filtering i-Evalue: ≤ 1e-5
- Minimum HMM coverage: 35%
- Minimum query coverage: 35%

**Defensome Mapping Rules:**

The mapping file (`defensome_map_pfam.tsv`) defines 23 Pfam domains across 17 gene families. Key logic:

| Family | Pfam ID(s) | Special Rules |
|--------|------------|---------------|
| CYP | PF00067 | Length filter ≥ 250 aa for candidate extraction |
| CCE | PF00135 | None |
| GST | PF02798 + PF00043 | Both GST_N and GST_C required for "strict" GST call |
| ABC transporter | PF00005 + PF00664 | Both NBD and TMD required for "strict" call |
| Nuclear receptor | PF00104 + PF00105 | Both LBD and DBD required (co-occurrence) |
| bHLH-PAS | PF00010 + (PF00989 or PF13426) | bHLH domain requires PAS domain |
| Keap1 | PF00651, PF07707, PF01344 | Individual domains counted (BTB, BACK, Kelch) |

**Rationale for thresholds:**

- **35% coverage filters** balance sensitivity vs. specificity for divergent insect sequences
- **Co-occurrence rules** reduce false positives from partial domain matches (e.g., standalone bHLH domains that are not bHLH-PAS transcription factors)
- **250 aa CYP filter** excludes truncated predictions and pseudogenes

### Module 2: Aggregation and CYP Clan Assignment

**CYP Clan Classification:**

Cytochrome P450s are classified into four major metazoan clans based on deep evolutionary relationships:

| Clan | Typical Functions | Drosophila Examples |
|------|-------------------|---------------------|
| CYP2 | Ecdysteroid metabolism, xenobiotic | Cyp18a1, Cyp306a1 |
| CYP3 | Xenobiotic metabolism (major expansion) | Cyp6, Cyp9, Cyp28 families |
| CYP4 | Fatty acid metabolism, pheromones, xenobiotic | Cyp4 family |
| MITO | Mitochondrial, ecdysteroid biosynthesis | Cyp12, Cyp314a1 |

**Method:**

1. Download FlyBase *D. melanogaster* CYP sequences and gene group annotations
2. Split Dmel CYPs by clan membership
3. Build per-clan profile HMMs (MAFFT alignment → hmmbuild)
4. Scan all Lepidoptera CYPs against the four clan HMMs
5. Assign to best-scoring clan
6. QC flag: sequences where best_bits - second_bits < threshold (ambiguous)

**Outputs:**
- `defensome_summary.tsv`: Copy numbers per species per family
- `defensome_presence_absence.tsv`: Binary matrix
- `cyp_clan_calls.hmm.tsv`: Per-sequence clan assignments

### Module 3: Phylogenetic Analysis

**Alignment:** MAFFT v7.505 with `--auto` (selects L-INS-i for <200 seqs, FFT-NS-2 for larger)

**Tree inference:** IQ-TREE2 v2.2.2
- Model selection: ModelFinder (`-m MFP`)
- Branch support: 1000 ultrafast bootstrap replicates (`-B 1000`)
- Output: Newick tree files per clan

---

## Quality Control

### Proteome QC Metrics

For each input proteome, we computed:

| Metric | Good | Suspect | Failed |
|--------|------|---------|--------|
| Pfam hits per protein | > 4.0 | 1.0 - 4.0 | < 1.0 |
| Defensome genes | > 300 | 100 - 300 | < 100 |
| CYP candidates | > 80 | 40 - 80 | < 40 |

**Issue identified:** Initial runs showed 16 species with < 0.1 hits/protein due to SLURM job truncation (12h timeout). Resolution: Extended walltime to 24h, added completion marker check (`# [ok]`) to hmmscan output validation.

### CYP Clan Assignment QC

Sequences flagged as ambiguous if:
- `best_bits - second_bits < min_bits_delta` (default: 50 bits)
- Indicates uncertain clan membership, possibly novel/divergent CYPs

Flagged sequences can be included or excluded from downstream phylogenetics via `--include_ambiguous`.

---

## Output Files

### Per-species (Module 1)
```
results_python/
├── 01_pfam/{species}.domtblout          # Raw hmmscan output
├── 02_parsed/{species}.pfam_hits.filtered.tsv   # Filtered Pfam hits
├── 03_defensome/
│   ├── {species}.defensome_gene_calls.tsv       # Gene-level calls
│   ├── {species}.defensome_family_counts.tsv    # Family counts
│   └── {species}.cyp_candidates.faa             # CYP sequences
```

### Aggregated (Module 2)
```
results_python/
├── 04_cyp/
│   ├── all_lep_cyps.faa                 # Merged CYP sequences
│   └── clan_calls/cyp_clan_calls.hmm.tsv  # Clan assignments
├── 05_figures/
│   ├── defensome_summary.tsv            # Full copy number table
│   ├── defensome_presence_absence.tsv   # Binary matrix
│   └── defensome_core_shell_cloud.tsv   # Family prevalence classification
```

### Phylogenetics (Module 3)
```
results_python/
├── 04_cyp/tree_by_clan/
│   ├── {CLAN}/{CLAN}.aln.faa            # Multiple sequence alignment
│   ├── {CLAN}/{CLAN}.aln.faa.treefile   # ML tree (Newick)
│   └── {CLAN}/{CLAN}.aln.faa.iqtree     # IQ-TREE report
```

---

## Computational Resources

| Module | Jobs | CPUs/job | RAM/job | Time/job |
|--------|------|----------|---------|----------|
| Module 1 | 42 (array) | 8 | 8 GB | 2-8 hours |
| Module 2 | 1 | 8 | 16 GB | 30-60 min |
| Module 3 | 1 | 8 | 16 GB | 4-12 hours |

**Total:** ~200-400 CPU-hours depending on proteome sizes

---

## Software Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| HMMER | 3.1b2+ | Pfam annotation, clan HMM building |
| MAFFT | 7.505 | Multiple sequence alignment |
| IQ-TREE2 | 2.2.2 | Maximum likelihood phylogenetics |
| Python | 3.10+ | Pipeline orchestration |
| pandas | 1.5+ | Data manipulation |
| matplotlib | 3.5+ | Visualization |

---

## Key Design Decisions

1. **Pfam-based annotation** rather than BLAST: Provides standardized, interpretable domain assignments with established family definitions.

2. **Conservative filtering**: 35% coverage thresholds and co-occurrence rules reduce false positives at the cost of some sensitivity. Appropriate for comparative analysis where consistency across species matters more than exhaustive detection within species.

3. **FlyBase-anchored CYP clans**: Using *D. melanogaster* as reference ensures reproducible clan assignments with well-characterized functional annotations.

4. **Per-clan phylogenetics**: Separate trees for each CYP clan enable:
   - Manageable tree sizes (hundreds rather than thousands of tips)
   - Biologically meaningful groupings
   - Detection of lineage-specific expansions within clans

5. **Modular pipeline**: Three independent modules with checkpointing enable:
   - Parallel execution across HPC nodes
   - Recovery from failures without full re-run
   - Flexible re-analysis (e.g., change tree parameters without redoing Pfam annotation)

---

## Limitations and Caveats

1. **Annotation quality dependence**: Results are only as good as input proteomes. Species with fragmented assemblies or incomplete gene prediction will undercount defensome genes.

2. **Pfam coverage**: Some defensome-relevant genes may lack Pfam domains or have divergent domains that fail HMM thresholds.

3. **Orthology vs. paralogy**: Copy number reflects total genes per family, not 1:1 orthologs. Gene family expansions in one lineage inflate counts.

4. **CYP clan assignment**: Based on sequence similarity to Dmel references. Highly divergent or novel insect-specific CYP clades may be mis-assigned or flagged as ambiguous.

5. **No expression data**: Presence does not imply function. Defensome genes may be tissue-specific, life-stage-specific, or induced only under xenobiotic exposure.

---

## Citation

If using this pipeline, cite:

- HMMER: Eddy SR (2011) PLOS Comput Biol 7:e1002195
- Pfam: Mistry J et al. (2021) Nucleic Acids Res 49:D412-D419
- MAFFT: Katoh K & Standley DM (2013) Mol Biol Evol 30:772-780
- IQ-TREE2: Minh BQ et al. (2020) Mol Biol Evol 37:1530-1534
- FlyBase: Gramates LS et al. (2022) Genetics 221:iyac035

---

*Report generated: January 2026*
*Pipeline version: 1.0*
