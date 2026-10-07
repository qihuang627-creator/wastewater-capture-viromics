#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
PREFIX="${2:-150}"

CDHIT_DUP="${CDHIT_DUP:-cd-hit-dup}"

R1="raw/${RUN}_1.fastq.gz"
R2="raw/${RUN}_2.fastq.gz"

OUTDIR="results/per_sample/dedup/u${PREFIX}"
TMPDIR="${OUTDIR}/tmp"
LOG="logs/${RUN}_cdhitdup_u${PREFIX}.log"

mkdir -p "${OUTDIR}" "${TMPDIR}" logs

if ! command -v "${CDHIT_DUP}" >/dev/null 2>&1; then
    echo "ERROR: cd-hit-dup not found: ${CDHIT_DUP}" >&2
    exit 1
fi

[[ -f "${R1}" ]] || {
    echo "ERROR: missing ${R1}" >&2
    exit 1
}

[[ -f "${R2}" ]] || {
    echo "ERROR: missing ${R2}" >&2
    exit 1
}

echo "=== Decompressing input FASTQ ==="

gzip -dc "${R1}" > "${TMPDIR}/${RUN}_R1.fastq"
gzip -dc "${R2}" > "${TMPDIR}/${RUN}_R2.fastq"

echo "=== Running cd-hit-dup: prefix=${PREFIX} nt ==="

"${CDHIT_DUP}" \
    -i  "${TMPDIR}/${RUN}_R1.fastq" \
    -i2 "${TMPDIR}/${RUN}_R2.fastq" \
    -o  "${TMPDIR}/${RUN}_R1.dedup.fastq" \
    -o2 "${TMPDIR}/${RUN}_R2.dedup.fastq" \
    -u "${PREFIX}" \
    > "${LOG}" 2>&1

echo "=== Compressing deduplicated FASTQ ==="

gzip -c "${TMPDIR}/${RUN}_R1.dedup.fastq" \
    > "${OUTDIR}/${RUN}_R1.dedup.fastq.gz"

gzip -c "${TMPDIR}/${RUN}_R2.dedup.fastq" \
    > "${OUTDIR}/${RUN}_R2.dedup.fastq.gz"

echo "=== Calculating duplicate statistics ==="

INPUT_R1=$(awk 'END{print NR/4}' "${TMPDIR}/${RUN}_R1.fastq")
INPUT_R2=$(awk 'END{print NR/4}' "${TMPDIR}/${RUN}_R2.fastq")

OUTPUT_R1=$(awk 'END{print NR/4}' "${TMPDIR}/${RUN}_R1.dedup.fastq")
OUTPUT_R2=$(awk 'END{print NR/4}' "${TMPDIR}/${RUN}_R2.dedup.fastq")

if [[ "${INPUT_R1}" != "${INPUT_R2}" ]]; then
    echo "ERROR: input R1/R2 pair counts differ" >&2
    exit 1
fi

if [[ "${OUTPUT_R1}" != "${OUTPUT_R2}" ]]; then
    echo "ERROR: deduplicated R1/R2 pair counts differ" >&2
    exit 1
fi

python3 - "${INPUT_R1}" "${OUTPUT_R1}" "${PREFIX}" \
    "${OUTDIR}/${RUN}_dedup_stats.tsv" <<'PY'
import sys

input_pairs = int(sys.argv[1])
unique_pairs = int(sys.argv[2])
prefix = int(sys.argv[3])
outfile = sys.argv[4]

duplicate_pairs = input_pairs - unique_pairs
duplicate_fraction = duplicate_pairs / input_pairs

with open(outfile, "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"input_pairs\t{input_pairs}\n")
    f.write(f"unique_pairs\t{unique_pairs}\n")
    f.write(f"duplicate_pairs\t{duplicate_pairs}\n")
    f.write(f"duplicate_fraction\t{duplicate_fraction:.6f}\n")
    f.write(f"duplicate_percent\t{duplicate_fraction*100:.2f}\n")
    f.write(f"prefix_length\t{prefix}\n")

print(f"Input pairs:       {input_pairs:,}")
print(f"Unique pairs:      {unique_pairs:,}")
print(f"Duplicate pairs:   {duplicate_pairs:,}")
print(f"Duplicate percent: {duplicate_fraction*100:.2f}%")
print(f"Prefix length:     {prefix} nt")
PY

echo
echo "=== Output files ==="
ls -lh \
    "${OUTDIR}/${RUN}_R1.dedup.fastq.gz" \
    "${OUTDIR}/${RUN}_R2.dedup.fastq.gz" \
    "${OUTDIR}/${RUN}_dedup_stats.tsv"

rm -rf "${TMPDIR}"

echo
echo "=== cd-hit-dup completed ==="
