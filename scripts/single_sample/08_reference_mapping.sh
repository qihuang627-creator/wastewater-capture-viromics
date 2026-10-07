#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
THREADS="${THREADS:-8}"

REF="refs/single_sample/target_panel/target_references.fa"
R1="results/per_sample/host_removal/${RUN}_R1.nonhuman.fastq.gz"
R2="results/per_sample/host_removal/${RUN}_R2.nonhuman.fastq.gz"

OUTDIR="results/single_sample/reference_mapping/${RUN}"
INDEX="${OUTDIR}/target_panel"
BAM="${OUTDIR}/${RUN}_target_panel.sorted.bam"

mkdir -p "${OUTDIR}" logs

echo "=== Building Bowtie2 reference panel ==="

bowtie2-build \
  "${REF}" \
  "${INDEX}" \
  > "logs/${RUN}_reference_index.log" 2>&1

echo "=== Competitive mapping ==="

bowtie2 \
  --very-sensitive \
  --no-mixed \
  --no-discordant \
  -x "${INDEX}" \
  -1 "${R1}" \
  -2 "${R2}" \
  -p "${THREADS}" \
  2> "logs/${RUN}_reference_mapping.log" \
| samtools view -@ "${THREADS}" -b - \
| samtools sort -@ "${THREADS}" -o "${BAM}" -

samtools index "${BAM}"

samtools flagstat \
  -@ "${THREADS}" \
  "${BAM}" \
  > "${OUTDIR}/${RUN}_target_panel.flagstat.txt"

samtools idxstats \
  "${BAM}" \
  > "${OUTDIR}/${RUN}_target_panel.idxstats.tsv"

echo "=== Depth: all primary mappings ==="

samtools depth \
  -aa \
  -q 20 \
  -Q 0 \
  "${BAM}" \
  > "${OUTDIR}/${RUN}_target_panel.depth.tsv"

echo "=== Depth: MAPQ >=20 ==="

samtools depth \
  -aa \
  -q 20 \
  -Q 20 \
  "${BAM}" \
  > "${OUTDIR}/${RUN}_target_panel.MQ20.depth.tsv"

echo "=== Properly paired fragments ==="

samtools view \
  -@ "${THREADS}" \
  -f 2 \
  -F 2308 \
  "${BAM}" \
| awk '{
    n[$3]++
}
END {
    for (r in n)
        printf "%s\t%.0f\n", r, n[r]/2
}' \
> "${OUTDIR}/${RUN}_proper_pair_fragments.tsv"

echo "=== Reference mapping complete ==="
