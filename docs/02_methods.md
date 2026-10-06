# Methods

## Dataset

Hybrid-capture wastewater metagenomic sequencing data were obtained from public project PRJEB87273.

Run ERR14788990 was selected as the primary capture library.

The library contained 199,495 paired-end read pairs.

## Read preprocessing

Putative PCR duplicates were collapsed with CD-HIT-DUP using a 150-nt paired-read prefix criterion.

Read quality control was performed with fastp.

Human-associated reads were screened against GRCh38 using Bowtie2 in very-sensitive mode, and only pairs with both mates unmapped were retained.

## Assembly and viral discovery

Nonhuman reads were assembled with MEGAHIT using multiple k-mer sizes.

Contigs ≥1 kb were screened with geNomad.

Target-family viral contigs were identified by combining geNomad taxonomy with nucleotide similarity searches against NCBI nucleotide databases.

## Read-back support

Final nonhuman reads were mapped competitively against candidate viral contigs.

Assembly-supported detection required:

1. geNomad classification to a targeted viral family;
2. ≥95% contig breadth at ≥1×; and
3. mean contig depth ≥10×.

## Reference-guided analysis

A curated reference panel was constructed using complete or near-complete public viral sequences.

Reads were competitively mapped against the panel using Bowtie2.

Genome breadth and depth were calculated using samtools.

Normalized sequencing support was expressed as proper-pair fragments per million input read pairs.

Because hybrid capture alters relative target representation, this metric is not interpreted as absolute viral abundance.

## Consensus generation

Targets selected for genome reconstruction were remapped independently to the corresponding reference sequence.

Variants were called using bcftools mpileup and bcftools call with haploid ploidy.

Variant filtering required QUAL ≥20 and depth ≥10.

Consensus positions with depth <10 were masked as N.

## Astrovirus typing

Classical HAstV genotype placement was evaluated using ORF2.

Sample sequences and public HAstV reference sequences were aligned using MAFFT.

Exploratory phylogenetic trees were generated using FastTree under a GTR+Gamma model.

## Diagnostic oligonucleotide surveillance

Published HAstV RT-qPCR primer/probe sequences were compared in silico with reconstructed HAstV sequences.

IUPAC degeneracy was explicitly supported.

Unknown consensus positions were recorded as uncallable rather than counted as nucleotide mismatches.

Primer 3'-terminal mismatches were tracked separately because they may be more consequential for amplification efficiency.

Sequence mismatch analysis does not by itself establish experimental assay failure.
