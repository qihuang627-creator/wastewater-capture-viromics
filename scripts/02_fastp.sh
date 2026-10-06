#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
THREADS="${THREADS:-4}"

mkdir -p qc/clean logs

conda run -n wwshotgun fastp \
  -i "raw/${RUN}_1.fastq.gz" \
  -I "raw/${RUN}_2.fastq.gz" \
  -o "qc/clean/${RUN}_1.clean.fastq.gz" \
  -O "qc/clean/${RUN}_2.clean.fastq.gz" \
  --thread "${THREADS}" \
  --detect_adapter_for_pe \
  --cut_right \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --length_required 50 \
  --json "qc/clean/${RUN}_fastp.json" \
  --html "qc/clean/${RUN}_fastp.html" \
  > "logs/${RUN}_fastp.log" 2>&1

echo "=== fastp completed ==="
ls -lh qc/clean/${RUN}*
