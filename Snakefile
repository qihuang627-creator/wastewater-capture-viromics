import csv
import os

configfile: "config/multisample/workflow.yaml"

SAMPLE_SHEET = config["samples"]

with open(SAMPLE_SHEET) as f:
    reader = csv.DictReader(f, delimiter="\t")
    SAMPLES = [row["sample"] for row in reader]

NT_DB = config["paths"]["nt_viruses_db"]
MAX_TARGET_SEQS = config["blast"]["max_target_seqs"]

PLOT_PYTHON = os.environ.get("PLOT_PYTHON", "python3")

BLAST_OUTFMT = (
    "6 qseqid sacc pident length qlen "
    "qstart qend qcovs evalue bitscore "
    "staxids sscinames stitle"
)


# ============================================================
# Final portfolio products
# ============================================================

rule all:
    input:
        "results_summary/multisample/"
        "all96_annotation_reconciled.tsv",

        "results_summary/multisample/"
        "all96_evidence_master.tsv",

        "results_summary/multisample/"
        "biological_group_evidence_matrix.tsv",

        "results_summary/multisample/"
        "portfolio_enteric_group_matrix.tsv",

        "figures/multisample/"
        "portfolio_enteric_evidence_heatmap.png",

        "figures/multisample/"
        "portfolio_enteric_evidence_heatmap.pdf",


# ============================================================
# PER-SAMPLE STAGE
# ============================================================

rule preprocess:
    input:
        r1="raw/{sample}_1.fastq.gz",
        r2="raw/{sample}_2.fastq.gz",
    output:
        clean_r1=(
            "qc/per_sample/clean/"
            "{sample}_R1.clean.fastq.gz"
        ),
        clean_r2=(
            "qc/per_sample/clean/"
            "{sample}_R2.clean.fastq.gz"
        ),
        nonhuman_r1=(
            "results/per_sample/host_removal/"
            "{sample}_R1.nonhuman.fastq.gz"
        ),
        nonhuman_r2=(
            "results/per_sample/host_removal/"
            "{sample}_R2.nonhuman.fastq.gz"
        ),
    threads:
        config["threads"]["preprocess"]
    shell:
        r"""
        THREADS={threads} \
        bash scripts/multisample/02_preprocess_one.sh \
        {wildcards.sample}
        """


rule assembly:
    input:
        r1=(
            "results/per_sample/host_removal/"
            "{sample}_R1.nonhuman.fastq.gz"
        ),
        r2=(
            "results/per_sample/host_removal/"
            "{sample}_R2.nonhuman.fastq.gz"
        ),
    output:
        contigs=(
            "results/per_sample/assembly/"
            "{sample}_megahit/final.contigs.fa"
        ),
        contigs_1kb=(
            "results/per_sample/viral_screen/"
            "{sample}_contigs_1kb.fa"
        ),
    threads:
        config["threads"]["assembly"]
    shell:
        r"""
        THREADS={threads} \
        bash scripts/multisample/04_assembly_one.sh \
        {wildcards.sample}
        """


rule genomad:
    input:
        contigs=(
            "results/per_sample/viral_screen/"
            "{sample}_contigs_1kb.fa"
        ),
    output:
        summary=(
            "results/per_sample/genomad/{sample}/"
            "{sample}_contigs_1kb_summary/"
            "{sample}_contigs_1kb_virus_summary.tsv"
        ),
    threads:
        config["threads"]["genomad"]
    shell:
        r"""
        THREADS={threads} \
        bash scripts/multisample/06_genomad_one.sh \
        {wildcards.sample}
        """


rule strict_targets:
    input:
        summary=(
            "results/per_sample/genomad/{sample}/"
            "{sample}_contigs_1kb_summary/"
            "{sample}_contigs_1kb_virus_summary.tsv"
        ),
        contigs=(
            "results/per_sample/viral_screen/"
            "{sample}_contigs_1kb.fa"
        ),
    output:
        fasta=(
            "results/multisample/target_candidates/"
            "{sample}/{sample}_strict_target_candidates.fa"
        ),
        metadata=(
            "results/multisample/target_candidates/"
            "{sample}/{sample}_strict_target_candidates.tsv"
        ),
    shell:
        r"""
        python3 \
        scripts/multisample/08_extract_strict_targets.py \
        {wildcards.sample}
        """


rule target_readback:
    input:
        r1=(
            "results/per_sample/host_removal/"
            "{sample}_R1.nonhuman.fastq.gz"
        ),
        r2=(
            "results/per_sample/host_removal/"
            "{sample}_R2.nonhuman.fastq.gz"
        ),
        fasta=(
            "results/multisample/target_candidates/"
            "{sample}/{sample}_strict_target_candidates.fa"
        ),
        metadata=(
            "results/multisample/target_candidates/"
            "{sample}/{sample}_strict_target_candidates.tsv"
        ),
    output:
        bam=(
            "results/multisample/target_readback/"
            "{sample}/{sample}_target_readback.bam"
        ),
        bai=(
            "results/multisample/target_readback/"
            "{sample}/{sample}_target_readback.bam.bai"
        ),
        idxstats=(
            "results/multisample/target_readback/"
            "{sample}/{sample}_idxstats.tsv"
        ),
        depth=(
            "results/multisample/target_readback/"
            "{sample}/{sample}_depth.tsv"
        ),
        support=(
            "results/multisample/target_readback/"
            "{sample}/{sample}_target_support.tsv"
        ),
    threads:
        config["threads"]["readback"]
    shell:
        r"""
        THREADS={threads} \
        bash scripts/multisample/09_target_readback_one.sh \
        {wildcards.sample}
        """


# ============================================================
# MULTISAMPLE AGGREGATION
# ============================================================

rule summarize_target_readback:
    input:
        expand(
            "results/multisample/target_readback/"
            "{sample}/{sample}_target_support.tsv",
            sample=SAMPLES,
        ),
    output:
        supported=(
            "results_summary/multisample/"
            "supported_target_contigs.tsv"
        ),
        sample_summary=(
            "results_summary/multisample/"
            "target_readback_summary_8samples.tsv"
        ),
        family_all=(
            "results_summary/multisample/"
            "target_family_contig_matrix.tsv"
        ),
        family_supported=(
            "results_summary/multisample/"
            "supported_target_family_matrix.tsv"
        ),
    shell:
        r"""
        python3 \
        scripts/multisample/10_summarize_target_readback.py
        """


rule prepare_supported_pool:
    input:
        "results_summary/multisample/"
        "supported_target_contigs.tsv",
    output:
        pool_dir=directory(
            "results/multisample/clustering/input"
        ),
        manifest=(
            "results_summary/multisample/"
            "supported_contig_manifest.tsv"
        ),
    shell:
        r"""
        python3 \
        scripts/multisample/11_prepare_supported_pool.py
        """


rule cluster_supported_targets:
    input:
        pool_dir=(
            "results/multisample/clustering/input"
        ),
        manifest=(
            "results_summary/multisample/"
            "supported_contig_manifest.tsv"
        ),
    output:
        cluster_dir=directory(
            "results/multisample/clustering/clusters"
        ),
    threads:
        config["threads"]["clustering"]
    shell:
        r"""
        THREADS={threads} \
        bash scripts/multisample/12_cluster_supported_targets.sh
        """


rule summarize_clusters:
    input:
        cluster_dir=(
            "results/multisample/clustering/clusters"
        ),
        manifest=(
            "results_summary/multisample/"
            "supported_contig_manifest.tsv"
        ),
    output:
        "results_summary/multisample/"
        "supported_sequence_clusters.tsv",
    shell:
        r"""
        python3 \
        scripts/multisample/13_summarize_clusters.py
        """


rule prepare_all96_representatives:
    input:
        clusters=(
            "results_summary/multisample/"
            "supported_sequence_clusters.tsv"
        ),
        cluster_dir=(
            "results/multisample/clustering/clusters"
        ),
    output:
        "results/multisample/clustering/"
        "all96_cluster_representatives.fa",
    shell:
        r"""
        python3 \
        scripts/multisample/21_prepare_all_cluster_representatives.py
        """


# ============================================================
# LOCAL nt_viruses ANNOTATION
# ============================================================

rule megablast_all96:
    input:
        fasta=(
            "results/multisample/clustering/"
            "all96_cluster_representatives.fa"
        ),
    output:
        "results/multisample/blast/"
        "local_nt_viruses_all96/all96.megablast.tsv",
    params:
        db=NT_DB,
        max_targets=MAX_TARGET_SEQS,
    threads:
        config["threads"]["blast"]
    shell:
        r"""
        mkdir -p \
          results/multisample/blast/local_nt_viruses_all96

        blastn \
          -task megablast \
          -query {input.fasta} \
          -db {params.db} \
          -out {output} \
          -outfmt '{BLAST_OUTFMT}' \
          -max_target_seqs {params.max_targets} \
          -num_threads {threads}
        """


rule summarize_megablast:
    input:
        blast=(
            "results/multisample/blast/"
            "local_nt_viruses_all96/all96.megablast.tsv"
        ),
        clusters=(
            "results_summary/multisample/"
            "supported_sequence_clusters.tsv"
        ),
    output:
        "results_summary/multisample/"
        "all96_ntviruses_annotation.tsv",
    shell:
        r"""
        python3 \
        scripts/multisample/22_summarize_all96_megablast.py
        """


rule extract_nonstrong_queries:
    input:
        annotation=(
            "results_summary/multisample/"
            "all96_ntviruses_annotation.tsv"
        ),
        fasta=(
            "results/multisample/clustering/"
            "all96_cluster_representatives.fa"
        ),
    output:
        "results/multisample/blast/"
        "local_nt_viruses_all96/nonstrong10.fa",
    shell:
        r"""
        python3 \
        scripts/multisample/22b_extract_nonstrong_queries.py
        """


rule sensitive_blastn:
    input:
        fasta=(
            "results/multisample/blast/"
            "local_nt_viruses_all96/nonstrong10.fa"
        ),
    output:
        "results/multisample/blast/"
        "local_nt_viruses_all96/nonstrong10.blastn.tsv",
    params:
        db=NT_DB,
        max_targets=MAX_TARGET_SEQS,
    threads:
        config["threads"]["blast"]
    shell:
        r"""
        if grep -q '^>' {input.fasta}; then

            blastn \
              -task blastn \
              -query {input.fasta} \
              -db {params.db} \
              -out {output} \
              -outfmt '{BLAST_OUTFMT}' \
              -max_target_seqs {params.max_targets} \
              -num_threads {threads}

        else
            : > {output}
        fi
        """


rule merge_all96_annotation:
    input:
        mega=(
            "results_summary/multisample/"
            "all96_ntviruses_annotation.tsv"
        ),
        sensitive=(
            "results/multisample/blast/"
            "local_nt_viruses_all96/nonstrong10.blastn.tsv"
        ),
    output:
        "results_summary/multisample/"
        "all96_annotation_master.tsv",
    shell:
        r"""
        python3 \
        scripts/multisample/23_merge_all96_annotation.py
        """


# ============================================================
# ALL96 REFERENCE + COMPETITIVE MAPPING
# ============================================================

rule build_all96_reference:
    input:
        "results/multisample/clustering/"
        "all96_cluster_representatives.fa",
    output:
        fasta=(
            "refs/multisample/all96/"
            "all96_cluster_representatives.fa"
        ),
        i1="refs/multisample/all96/bowtie2/all96.1.bt2",
        i2="refs/multisample/all96/bowtie2/all96.2.bt2",
        i3="refs/multisample/all96/bowtie2/all96.3.bt2",
        i4="refs/multisample/all96/bowtie2/all96.4.bt2",
        ir1=(
            "refs/multisample/all96/bowtie2/"
            "all96.rev.1.bt2"
        ),
        ir2=(
            "refs/multisample/all96/bowtie2/"
            "all96.rev.2.bt2"
        ),
    shell:
        r"""
        mkdir -p refs/multisample/all96/bowtie2

        cp {input} {output.fasta}

        rm -f \
          refs/multisample/all96/bowtie2/all96*.bt2

        bowtie2-build \
          {output.fasta} \
          refs/multisample/all96/bowtie2/all96
        """


rule all96_mapping:
    input:
        r1=(
            "results/per_sample/host_removal/"
            "{sample}_R1.nonhuman.fastq.gz"
        ),
        r2=(
            "results/per_sample/host_removal/"
            "{sample}_R2.nonhuman.fastq.gz"
        ),
        ref=(
            "refs/multisample/all96/"
            "all96_cluster_representatives.fa"
        ),
        i1="refs/multisample/all96/bowtie2/all96.1.bt2",
        i2="refs/multisample/all96/bowtie2/all96.2.bt2",
        i3="refs/multisample/all96/bowtie2/all96.3.bt2",
        i4="refs/multisample/all96/bowtie2/all96.4.bt2",
        ir1=(
            "refs/multisample/all96/bowtie2/"
            "all96.rev.1.bt2"
        ),
        ir2=(
            "refs/multisample/all96/bowtie2/"
            "all96.rev.2.bt2"
        ),
    output:
        bam=(
            "results/multisample/all96_mapping/{sample}/"
            "{sample}_all96.sorted.bam"
        ),
        bai=(
            "results/multisample/all96_mapping/{sample}/"
            "{sample}_all96.sorted.bam.bai"
        ),
        depth=(
            "results/multisample/all96_mapping/{sample}/"
            "{sample}_all96.mapq20.depth.tsv"
        ),
        fragments=(
            "results/multisample/all96_mapping/{sample}/"
            "{sample}_all96.mapq20.fragments.tsv"
        ),
        support=(
            "results/multisample/all96_mapping/{sample}/"
            "{sample}_all96_support.tsv"
        ),
    threads:
        config["threads"]["all96_mapping"]
    shell:
        r"""
        THREADS={threads} \
        bash scripts/multisample/24_all96_mapping_one.sh \
        {wildcards.sample}
        """


# ============================================================
# INTEGRATED EVIDENCE
# ============================================================

rule integrate_all96_evidence:
    input:
        clusters=(
            "results_summary/multisample/"
            "supported_sequence_clusters.tsv"
        ),
        annotation=(
            "results_summary/multisample/"
            "all96_annotation_master.tsv"
        ),
        mapping=expand(
            "results/multisample/all96_mapping/"
            "{sample}/{sample}_all96_support.tsv",
            sample=SAMPLES,
        ),
    output:
        master=(
            "results_summary/multisample/"
            "all96_evidence_master.tsv"
        ),
        matrix=(
            "results_summary/multisample/"
            "all96_evidence_matrix.tsv"
        ),
        sample_summary=(
            "results_summary/multisample/"
            "all96_evidence_sample_summary.tsv"
        ),
        recovered=(
            "results_summary/multisample/"
            "all96_mapping_recovered.tsv"
        ),
    shell:
        r"""
        python3 \
        scripts/multisample/26_integrate_all96_evidence.py
        """


rule reconcile_taxonomy:
    input:
        "results_summary/multisample/"
        "all96_annotation_master.tsv",
    output:
        "results_summary/multisample/"
        "all96_annotation_reconciled.tsv",
    shell:
        r"""
        python3 \
        scripts/multisample/27_reconcile_all96_taxonomy.py
        """


rule biological_grouping:
    input:
        "results_summary/multisample/"
        "all96_annotation_reconciled.tsv",
    output:
        "results_summary/multisample/"
        "all96_biological_grouping.tsv",
    shell:
        r"""
        python3 \
        scripts/multisample/28_build_biological_grouping.py
        """


rule group_level_evidence:
    input:
        grouping=(
            "results_summary/multisample/"
            "all96_biological_grouping.tsv"
        ),
        evidence=(
            "results_summary/multisample/"
            "all96_evidence_master.tsv"
        ),
    output:
        long=(
            "results_summary/multisample/"
            "biological_group_evidence.tsv"
        ),
        matrix=(
            "results_summary/multisample/"
            "biological_group_evidence_matrix.tsv"
        ),
        summary=(
            "results_summary/multisample/"
            "biological_group_sample_summary.tsv"
        ),
    shell:
        r"""
        python3 \
        scripts/multisample/29_build_group_level_evidence.py
        """


rule portfolio_matrix:
    input:
        "results_summary/multisample/"
        "biological_group_evidence_matrix.tsv",
    output:
        "results_summary/multisample/"
        "portfolio_enteric_group_matrix.tsv",
    shell:
        r"""
        python3 \
        scripts/multisample/31_build_portfolio_matrix.py
        """


rule portfolio_heatmap:
    input:
        "results_summary/multisample/"
        "portfolio_enteric_group_matrix.tsv",
    output:
        png=(
            "figures/multisample/"
            "portfolio_enteric_evidence_heatmap.png"
        ),
        pdf=(
            "figures/multisample/"
            "portfolio_enteric_evidence_heatmap.pdf"
        ),
    params:
        python=PLOT_PYTHON,
    shell:
        r"""
        MPLBACKEND=Agg {params.python} \
        scripts/multisample/32_plot_portfolio_heatmap.py
        """


# ============================================================
# Optional full biological-group heatmap
# Not required by rule all.
# ============================================================

rule full_biological_heatmap:
    input:
        "results_summary/multisample/"
        "biological_group_evidence_matrix.tsv",
    output:
        png=(
            "figures/multisample/"
            "full_biological_group_evidence_heatmap.png"
        ),
        pdf=(
            "figures/multisample/"
            "full_biological_group_evidence_heatmap.pdf"
        ),
    params:
        python=PLOT_PYTHON,
    shell:
        r"""
        MPLBACKEND=Agg {params.python} \
        scripts/multisample/30_plot_biological_group_heatmap.py
        """
