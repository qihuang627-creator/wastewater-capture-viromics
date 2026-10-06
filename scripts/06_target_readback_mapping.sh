#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
THREADS="${THREADS:-8}"

FA="results/target_candidates/${RUN}_target_candidates.fa"

R1="results/host_removal_dedup/${RUN}_R1.nonhuman.fastq.gz"
R2="results/host_removal_dedup/${RUN}_R2.nonhuman.fastq.gz"

OUTDIR="results/target_readback/${RUN}"
INDEX="${OUTDIR}/target_candidates"
BAM="${OUTDIR}/${RUN}_target_candidates.sorted.bam"

mkdir -p "${OUTDIR}" logs

echo "=== Building Bowtie2 index ==="

bowtie2-build \
  "${FA}" \
  "${INDEX}" \
  > "logs/${RUN}_target_index.log" 2>&1

echo "=== Competitive read-back mapping ==="

bowtie2 \
  --very-sensitive \
  -x "${INDEX}" \
  -1 "${R1}" \
  -2 "${R2}" \
  -p "${THREADS}" \
  2> "logs/${RUN}_target_readback_bowtie2.log" \
| samtools view -@ "${THREADS}" -b - \
| samtools sort -@ "${THREADS}" -o "${BAM}" -

samtools index "${BAM}"

samtools flagstat \
  -@ "${THREADS}" \
  "${BAM}" \
  > "${OUTDIR}/${RUN}_target_readback.flagstat.txt"

samtools idxstats \
  "${BAM}" \
  > "${OUTDIR}/${RUN}_target_readback.idxstats.tsv"

echo "=== Calculating depth ==="

samtools depth \
  -aa \
  -Q 10 \
  "${BAM}" \
  > "${OUTDIR}/${RUN}_target_readback.depth.tsv"

echo "=== Done ==="
