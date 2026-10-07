#!/usr/bin/env bash
set -euo pipefail

RUN="${1:?Usage: $0 RUN}"
THREADS="${THREADS:-8}"

IMAGE="community.wave.seqera.io/library/genomad:1.12.0--27836e6e665e84b5"
DB="${GENOMAD_DB:-databases/genomad/genomad_db}"

INPUT="results/per_sample/viral_screen/${RUN}_contigs_1kb.fa"
OUTDIR="results/per_sample/genomad/${RUN}"
LOG="logs/multisample/${RUN}_genomad.log"

if [[ ! -d "${DB}" ]]; then
    echo "ERROR: geNomad database directory not found: ${DB}" >&2
    exit 1
fi

DB="$(cd "${DB}" && pwd)"

mkdir -p "${OUTDIR}" logs/multisample

SUMMARY="${OUTDIR}/${RUN}_contigs_1kb_summary/${RUN}_contigs_1kb_virus_summary.tsv"

if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: missing input ${INPUT}" >&2
    exit 1
fi

echo
echo "=========================================="
echo "geNomad: ${RUN}"
echo "=========================================="

if [[ -s "${SUMMARY}" ]]; then
    echo "Existing geNomad summary found — skipping"
    echo "${SUMMARY}"
    exit 0
fi

docker run --rm \
  -v "${PWD}:/work" \
  -v "${DB}:/genomad_db:ro" \
  -w /work \
  "${IMAGE}" \
  genomad end-to-end \
  "/work/${INPUT}" \
  "/work/${OUTDIR}" \
  /genomad_db \
  --threads "${THREADS}" \
  > "${LOG}" 2>&1

if [[ ! -s "${SUMMARY}" ]]; then
    echo "ERROR: geNomad finished but summary not found:" >&2
    echo "${SUMMARY}" >&2
    exit 1
fi

N=$(awk 'NR>1 {n++} END{print n+0}' "${SUMMARY}")

echo "Viral contigs: ${N}"
echo "Completed: ${RUN}"
