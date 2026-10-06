#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
THREADS="${THREADS:-4}"

mkdir -p qc/raw logs

echo "=== Raw FastQC: ${RUN} ==="

fastqc \
  -t "${THREADS}" \
  -o qc/raw \
  "raw/${RUN}_1.fastq.gz" \
  "raw/${RUN}_2.fastq.gz" \
  > "logs/${RUN}_raw_fastqc.log" 2>&1

echo "=== FastQC completed ==="
ls -lh qc/raw/${RUN}_*
