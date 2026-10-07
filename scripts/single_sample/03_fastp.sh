#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
THREADS="${THREADS:-4}"

mkdir -p qc/per_sample/clean logs

R1="results/per_sample/dedup/u150/${RUN}_R1.dedup.fastq.gz"
R2="results/per_sample/dedup/u150/${RUN}_R2.dedup.fastq.gz"

[[ -f "${R1}" ]] || {
    echo "Missing deduplicated R1: ${R1}" >&2
    exit 1
}

[[ -f "${R2}" ]] || {
    echo "Missing deduplicated R2: ${R2}" >&2
    exit 1
}

conda run -n wwshotgun fastp \
  -i "${R1}" \
  -I "${R2}" \
  -o "qc/per_sample/clean/${RUN}_R1.clean.fastq.gz" \
  -O "qc/per_sample/clean/${RUN}_R2.clean.fastq.gz" \
  --thread "${THREADS}" \
  --detect_adapter_for_pe \
  --cut_right \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --length_required 50 \
  --json "qc/per_sample/clean/${RUN}.fastp.json" \
  --html "qc/per_sample/clean/${RUN}.fastp.html" \
  > "logs/${RUN}_fastp.log" 2>&1

echo "=== fastp completed ==="
ls -lh qc/per_sample/clean/${RUN}*
