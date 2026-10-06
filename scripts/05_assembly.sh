#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
THREADS="${THREADS:-8}"

R1="results/host_removal_dedup/${RUN}_R1.nonhuman.fastq.gz"
R2="results/host_removal_dedup/${RUN}_R2.nonhuman.fastq.gz"

OUTDIR="results/assembly/${RUN}_megahit"
LOG="logs/${RUN}_megahit.log"

mkdir -p results/assembly logs

rm -rf "${OUTDIR}"

conda run -n metagenome megahit \
  -1 "${R1}" \
  -2 "${R2}" \
  -o "${OUTDIR}" \
  --min-contig-len 500 \
  --k-list 21,29,39,59,79,99,119,141 \
  -t "${THREADS}" \
  > "${LOG}" 2>&1

echo "=== MEGAHIT completed ==="
ls -lh "${OUTDIR}/final.contigs.fa"
