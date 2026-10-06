#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
THREADS="${THREADS:-8}"

DB="${HOME}/QI/wastewater-shotgun/databases/human_GRCh38/GRCh38"

R1="qc/dedup_clean/${RUN}_R1.clean.fastq.gz"
R2="qc/dedup_clean/${RUN}_R2.clean.fastq.gz"

OUTDIR="results/host_removal_dedup"
mkdir -p "${OUTDIR}" logs

echo "=== GRCh38 screening: ${RUN} ==="

bowtie2 \
  --very-sensitive \
  -x "${DB}" \
  -1 "${R1}" \
  -2 "${R2}" \
  -p "${THREADS}" \
  2> "logs/${RUN}_dedup_GRCh38.log" \
| samtools view \
    -@ "${THREADS}" \
    -b \
    -o "${OUTDIR}/${RUN}_vs_GRCh38.bam" -

samtools flagstat \
  -@ "${THREADS}" \
  "${OUTDIR}/${RUN}_vs_GRCh38.bam" \
  > "${OUTDIR}/${RUN}_GRCh38.flagstat.txt"

samtools view \
  -@ "${THREADS}" \
  -b \
  -f 12 \
  -F 2304 \
  "${OUTDIR}/${RUN}_vs_GRCh38.bam" \
| samtools sort \
    -n \
    -@ "${THREADS}" \
    -o "${OUTDIR}/${RUN}_both_unmapped.namesort.bam" -

samtools fastq \
  -@ "${THREADS}" \
  -n \
  -1 "${OUTDIR}/${RUN}_R1.nonhuman.fastq.gz" \
  -2 "${OUTDIR}/${RUN}_R2.nonhuman.fastq.gz" \
  -0 /dev/null \
  -s /dev/null \
  "${OUTDIR}/${RUN}_both_unmapped.namesort.bam"

echo "=== Completed ==="
cat "${OUTDIR}/${RUN}_GRCh38.flagstat.txt"
