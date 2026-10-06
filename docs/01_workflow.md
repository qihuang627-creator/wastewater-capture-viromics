# Workflow

## 1. Input

Hybrid-capture paired-end wastewater metagenomic sequencing data.

Primary demonstration sample:

`ERR14788990`

## 2. Putative PCR duplicate removal

CD-HIT-DUP was applied before read trimming using a 150-nt paired-read prefix criterion.

Because the library contains no UMI information, duplicate calls are interpreted as putative PCR duplicates.

## 3. Quality control

Reads were processed with fastp using:

- sliding-window quality trimming
- minimum Phred score 20
- maximum 40% low-quality bases
- maximum 5 ambiguous bases
- minimum retained length 50 bp

## 4. Human-read removal

Reads were aligned against GRCh38 with Bowtie2 `--very-sensitive`.

Only read pairs in which both mates were unmapped were retained.

## 5. De novo assembly

Nonhuman reads were assembled with MEGAHIT.

Contigs ≥1 kb were retained for downstream viral discovery.

## 6. Viral identification

geNomad was used for viral contig classification.

Candidate viral contigs assigned to GastroCap target families were retained for target-focused analysis.

## 7. Target confirmation

Candidate contigs were compared against NCBI nucleotide reference sequences.

Evidence from:

- viral classification
- nucleotide similarity
- contig read-back mapping

was integrated for assembly-supported target detection.

## 8. Competitive reference mapping

Capture reads were competitively mapped to a curated viral reference panel.

Reference-level metrics included:

- mapped read support
- proper-pair fragment support
- breadth ≥1×
- breadth ≥5×
- breadth ≥10×
- mean depth
- MAPQ ≥20 breadth

## 9. Consensus reconstruction

High-quality targets were remapped independently to avoid competition between closely related references.

Variants were called with bcftools.

Consensus bases required:

- mapping quality ≥20
- base quality ≥20
- depth ≥10×

Lower-depth positions were masked as `N`.

## 10. Genotype and phylogenetic placement

Classical HAstV ORF2 sequences were aligned with MAFFT.

FastTree GTR+Gamma trees were used for exploratory genotype and within-genotype phylogenetic placement.

FastTree internal node values are interpreted as local support values rather than conventional bootstrap replicates.

## 11. Diagnostic mismatch surveillance

Publicly available HAstV RT-qPCR primers and probes were compared against reconstructed viral sequences.

Mismatch interpretation distinguishes:

- internal primer mismatches
- primer 3'-terminal mismatches
- probe mismatches
- uncallable sequence positions

No proprietary capture-bait sequences were analyzed.
