# Wastewater Hybrid-Capture Viromics

A reproducible workflow for **viral detection, genome reconstruction, phylogenetic placement, and diagnostic oligonucleotide mismatch surveillance** from hybrid-capture wastewater metagenomic sequencing.

This project uses a publicly available hybrid-capture wastewater sequencing library from **PRJEB87273** and demonstrates an end-to-end computational workflow for GastroCap-enriched viromics.

---

## Overview

The workflow integrates:

- putative PCR duplicate removal
- read quality control and trimming
- human-read removal
- de novo metagenomic assembly
- viral contig identification with geNomad
- target-virus identification using nucleotide similarity searches
- candidate-contig read-back validation
- competitive reference mapping
- genome breadth and depth estimation
- normalized sequencing support
- reference-guided consensus reconstruction
- genotype and phylogenetic placement
- in-silico diagnostic primer/probe mismatch surveillance

```text
FASTQ
  │
  ├── Putative PCR deduplication
  │
  ├── QC / trimming
  │
  ├── Human-read removal
  │
  ├── De novo assembly
  │
  ├── Viral classification
  │
  ├── Target-virus identification
  │
  ├── Candidate read-back validation
  │
  ├── Competitive reference mapping
  │
  ├── Genome breadth / depth
  │
  ├── Normalized read support
  │
  ├── Consensus reconstruction
  │
  ├── Genotype / phylogenetic placement
  │
  └── Diagnostic primer/probe mismatch surveillance
```

---

## Dataset

Primary hybrid-capture library:

| Field | Value |
|---|---|
| Project | `PRJEB87273` |
| Run | `ERR14788990` |
| BioSample | `SAMEA117873011` |
| Location | Copenhagen, Denmark |
| Collection date | 2018-03-01 |
| Sequencing | Paired-end Illumina MiSeq |
| Raw read pairs | 199,495 |

The raw sequencing data are publicly available and are **not redistributed in this repository**.

---

## Workflow

### 1. Putative PCR duplicate removal

Raw paired-end reads were processed with **CD-HIT-DUP** before trimming.

A 150-nt paired-read prefix criterion was used to identify sequence-level duplicate pairs.

Because the library does not contain UMIs, these reads are described as **putative PCR duplicates**, rather than experimentally verified molecular duplicates.

### 2. Quality control

Reads were processed with **fastp** using quality trimming and filtering criteria including:

- Phred quality threshold: 20
- sliding-window trimming
- maximum low-quality base fraction: 40%
- maximum ambiguous bases: 5
- minimum retained read length: 50 bp

### 3. Human-read removal

Quality-controlled reads were mapped against **GRCh38** using Bowtie2 in `--very-sensitive` mode.

Only paired reads for which both mates remained unmapped were retained for downstream analysis.

### 4. De novo assembly

Nonhuman reads were assembled using **MEGAHIT** with multiple k-mer sizes.

Contigs ≥1 kb were retained for viral discovery.

### 5. Viral discovery

Contigs were screened with **geNomad**.

Candidate viral contigs belonging to GastroCap-targeted viral families were retained for target-focused analysis.

### 6. Target identification

Candidate viral contigs were compared against public nucleotide reference sequences.

Target identification integrated evidence from:

- geNomad viral classification
- nucleotide sequence similarity
- contig-level read-back support

### 7. Competitive reference mapping

A curated target reference panel was constructed from complete or near-complete public viral sequences.

Capture reads were competitively mapped against the reference panel.

Reference-level metrics included:

- mapped reads
- proper-pair fragments
- breadth ≥1×
- breadth ≥5×
- breadth ≥10×
- mean depth
- MAPQ ≥20 breadth and depth

### 8. Consensus reconstruction

Targets selected for genome reconstruction were remapped independently to their corresponding references to reduce competition between related references.

Variants were called with **bcftools**.

Consensus bases required:

- mapping quality ≥20
- base quality ≥20
- depth ≥10×

Positions below 10× depth were masked as `N`.

### 9. Genotype and phylogenetic placement

Classical human astrovirus typing focused on the **ORF2 capsid region**.

Sample and public reference sequences were aligned with **MAFFT**, and exploratory phylogenetic placement was performed using **FastTree** under a GTR+Gamma model.

### 10. Diagnostic oligonucleotide mismatch surveillance

Published human astrovirus RT-qPCR primer/probe sequences were compared in silico against reconstructed viral sequences.

The analysis distinguishes:

- internal primer mismatches
- primer 3′-terminal mismatches
- probe mismatches
- uncallable consensus positions

This analysis evaluates **sequence compatibility only** and does not directly establish experimental PCR sensitivity or assay failure.

---

## Key results

### Preprocessing and assembly

Starting from **199,495 raw read pairs**:

| Stage | Read pairs |
|---|---:|
| Raw | 199,495 |
| After putative PCR deduplication | 198,528 |
| After QC | 188,573 |
| After human-read removal | 188,573 |

Assembly produced:

- **1,597 contigs ≥500 bp**
- **1,183,881 bp** total assembly length
- **42,702 bp** maximum contig length
- **673 bp** N50
- **139 contigs ≥1 kb**

Among the ≥1-kb contigs:

- **54** were classified as viral by geNomad
- **18** were assigned to strict GastroCap target families
- **1** additional contig was classified as an unresolved Picornavirales candidate

---

## Reference-guided viral recovery

Several human enteric viral targets showed strong reference-guided support.

| Target | Reference | Callable sequence | Mean depth | Interpretation |
|---|---|---:|---:|---|
| HAstV-1 population A | `LC694990.1` | 98.92% | 397.9× | Near-complete consensus |
| HAstV-1 population B | `LC694997.1` | 47.84% | 21.6× | Partial divergent population |
| HAstV-3 | `MN444721.1` | 99.96% | 255.6× | Near-complete consensus |
| Rotavirus A segment 3 | `OR890193.1` | 83.64% | 32.4× | Partial high-confidence segment consensus |

---

## Multiple HAstV-1 sequence populations

Independent assembly, competitive mapping, consensus reconstruction, and ORF2 phylogenetic analysis supported the presence of **at least two genetically distinct HAstV-1 sequence populations**.

Two assembled contigs showed reciprocal reference specificity:

| Contig | HAstV1_A identity | HAstV1_B identity |
|---|---:|---:|
| `k141_843` | 99.08% | 88.97% |
| `k141_1995` | 90.56% | 99.64% |

The reconstructed HAstV-1 ORF2 sequences shared only:

- **1,997 comparable nucleotide sites**
- **1,802 matches**
- **90.24% nucleotide identity**

Population A was dominant and nearly completely reconstructed.

Population B showed substantially lower read support and was only partially reconstructed, but phylogenetic analysis placed it within a distinct HAstV-1 public-sequence cluster.

---

## HAstV-3 reconstruction

The HAstV-3 consensus showed:

- reference length: **6,790 bp**
- callable sequence: **99.96%**
- masked positions: **3 bp**
- mean depth: **255.6×**
- high-confidence reference-relative SNPs: **15**

The reconstructed ORF2 sequence showed **99.75% nucleotide identity** to `MN444721.1`.

ORF2 phylogenetic placement supported classification as **HAstV-3**.

---

## Rotavirus A segment 3

Rotavirus A segment 3 showed:

- reference length: **2,591 bp**
- breadth ≥1×: **97.61%**
- breadth ≥10×: **83.64%**
- mean depth: **32.42×**
- callable sequence: **83.64%**
- masked positions: **424 bp**

This sequence is therefore treated as a **partial high-confidence segment consensus**, rather than a near-complete segment reconstruction.

---

## Diagnostic oligonucleotide compatibility

Two published classical human astrovirus RT-qPCR assay designs were evaluated in silico.

Summary:

| Sample | Assay | Best forward primer mismatches | Reverse mismatches | Probe mismatches | Primer 3′ mismatches |
|---|---|---:|---:|---:|---:|
| HAstV1_A | BCCDC classical HAstV | 1 | 0 | 1 | 0 |
| HAstV1_A | Classic HAstV Japan | 1 | 0 | 0 | 0 |
| HAstV1_B | BCCDC classical HAstV | 1 | 0 | 0 | 0 |
| HAstV1_B | Classic HAstV Japan | 0 | 0 | 0 | 0 |
| HAstV3 | BCCDC classical HAstV | 1 | 0 | 2 | 0 |
| HAstV3 | Classic HAstV Japan | 0 | 1 | 0 | 0 |

No primer 3′-terminal mismatches were observed in the evaluated sequences.

Notably, the divergent HAstV1_B population remained fully compatible with the forward, reverse, and probe oligonucleotides of one evaluated classical HAstV assay despite substantial ORF2 divergence from HAstV1_A.

The BCCDC HAstV probe contained two mismatches against the HAstV-3 consensus.

These results illustrate that overall viral sequence divergence does not necessarily predict divergence at diagnostic oligonucleotide binding sites.

---

## Repository structure

```text
.
├── README.md
├── LICENSE
├── .gitignore
│
├── scripts/
│   ├── 03_dedup_raw.sh
│   ├── 04_host_removal.sh
│   ├── 05_assembly.sh
│   ├── 06_target_readback_mapping.sh
│   ├── 07_fetch_target_references.sh
│   ├── 08_reference_mapping.sh
│   ├── 09_consensus_target.sh
│   └── 10_primer_probe_mismatch_v2.py
│
├── docs/
│   ├── 01_workflow.md
│   ├── 02_methods.md
│   └── 03_results.md
│
├── refs/
│   └── diagnostic_assays/
│       └── astrovirus_public_assays.tsv
│
├── results_summary/
│   ├── preprocessing_summary.tsv
│   ├── reference_mapping_summary.tsv
│   ├── consensus_summary.tsv
│   ├── hastv_population_summary.tsv
│   └── diagnostic_mismatch_summary.tsv
│
└── environment/
    └── software_versions.txt
```

Large sequencing files, BAM files, databases, reference indices, and intermediate analysis outputs are intentionally excluded from the repository.

---

## Main software

The workflow uses:

- fastp
- CD-HIT-DUP
- Bowtie2
- samtools
- bcftools
- MEGAHIT
- geNomad
- NCBI BLAST+
- MAFFT
- FastTree
- Python 3

Exact software versions used for this analysis are recorded in:

```text
environment/software_versions.txt
```

---

## Important methodological notes

### Putative PCR duplicates

Sequence-level deduplication was performed using a 150-nt paired-read prefix criterion.

Because the sequencing library does not contain UMIs, duplicate calls should not be interpreted as experimentally confirmed PCR duplicates.

### Hybrid-capture abundance

Hybrid capture changes the relative representation of viral targets.

Therefore, normalized mapped-fragment counts are interpreted as **normalized sequencing support**, not as absolute viral abundance in the original wastewater sample.

### Consensus interpretation

Consensus genomes represent the dominant nucleotide state among sufficiently supported reads.

They should not be interpreted as evidence that only one viral strain was present in the sample.

### Phylogenetic interpretation

FastTree analyses were used for exploratory genotype and within-genotype sequence placement.

FastTree internal-node values are interpreted as **local support values**, not conventional bootstrap replicates.

### Diagnostic mismatch interpretation

Primer/probe mismatch analysis evaluates sequence compatibility only.

A mismatch does not by itself demonstrate reduced assay sensitivity, increased Ct values, or diagnostic failure.

### Capture probes

The GastroCap capture-bait sequences are proprietary and were not evaluated.

Diagnostic oligonucleotide analysis was restricted to publicly available primer/probe sequences.

---

## Reproducibility

Core scripts are organized in approximate workflow order under `scripts/`.

Additional documentation is available in:

- [`docs/01_workflow.md`](docs/01_workflow.md)
- [`docs/02_methods.md`](docs/02_methods.md)
- [`docs/03_results.md`](docs/03_results.md)

Compact analysis outputs are provided under:

- [`results_summary/`](results_summary/)

Raw sequencing data can be retrieved independently from the corresponding public sequencing archive using the accession information above.

---

## Scope

This repository is intended as a focused demonstration of **hybrid-capture wastewater viromics and sequence-resolved viral surveillance**.

It is not intended to serve as a validated clinical diagnostic pipeline.

---

## Author

**Qi Huang**

Environmental biotechnology · wastewater viromics · microbial genomics · multi-omics · bioinformatics
