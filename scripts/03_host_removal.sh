#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
THREADS="${THREADS:-8}"

DB="${HOME}/QI/wastewater-shotgun/databases/human_GRCh38/GRCh38"

R1="qc/clean/${RUN}_1.clean.fastq.gz"
R2="qc/clean/${RUN}_2.clean.fastq.gz"

OUTDIR="results/host_removal"
mkdir -p "${OUTDIR}" logs

# Sanity checks
[[ -f "${R1}" ]] || { echo "Missing ${R1}" >&2; exit 1; }
[[ -f "${R2}" ]] || { echo "Missing ${R2}" >&2; exit 1; }

if ! ls "${DB}"*.bt2* >/dev/null 2>&1; then
    echo "Bowtie2 GRCh38 index not found: ${DB}" >&2
    exit 1
fi

echo "=== Mapping ${RUN} against GRCh38 ==="

bowtie2 \
  --very-sensitive \
  -x "${DB}" \
  -1 "${R1}" \
  -2 "${R2}" \
  -p "${THREADS}" \
  2> "logs/${RUN}_bowtie2_GRCh38.log" \
| samtools view \
    -@ "${THREADS}" \
    -b \
    -o "${OUTDIR}/${RUN}_vs_GRCh38.bam" -

echo "=== Flagstat ==="

samtools flagstat \
  -@ "${THREADS}" \
  "${OUTDIR}/${RUN}_vs_GRCh38.bam" \
  > "${OUTDIR}/${RUN}_GRCh38.flagstat.txt"

# Keep only proper read pairs for which BOTH mates are unmapped.
# 0x4 = read unmapped
# 0x8 = mate unmapped
# 0x100/0x800 excluded = secondary/supplementary
echo "=== Extracting non-human pairs ==="

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

echo "=== Host removal completed ==="
cat "${OUTDIR}/${RUN}_GRCh38.flagstat.txt"

echo
echo "=== Output FASTQ ==="
ls -lh \
  "${OUTDIR}/${RUN}_R1.nonhuman.fastq.gz" \
  "${OUTDIR}/${RUN}_R2.nonhuman.fastq.gz"
