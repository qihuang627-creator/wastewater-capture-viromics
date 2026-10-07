#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"

CDHIT_EST="${CDHIT_EST:-cd-hit-est}"

INDIR="results/multisample/clustering/input"
OUTDIR="results/multisample/clustering/clusters"

if ! command -v "${CDHIT_EST}" >/dev/null 2>&1; then
    echo "ERROR: cd-hit-est not found: ${CDHIT_EST}" >&2
    exit 1
fi

if [[ ! -d "${INDIR}" ]]; then
    echo "ERROR: input directory not found: ${INDIR}" >&2
    echo "Run scripts/multisample/11_prepare_supported_pool.py first." >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

for FA in "${INDIR}"/*viridae.fa
do
    [[ -e "${FA}" ]] || continue

    FAMILY=$(basename "${FA}" .fa)
    OUT="${OUTDIR}/${FAMILY}_ani95"

    echo
    echo "=========================================="
    echo "Clustering ${FAMILY}"
    echo "=========================================="

    N_INPUT=$(grep -c '^>' "${FA}")

    "${CDHIT_EST}" \
      -i "${FA}" \
      -o "${OUT}.fa" \
      -c 0.95 \
      -n 10 \
      -G 0 \
      -aS 0.80 \
      -g 1 \
      -r 1 \
      -d 0 \
      -T "${THREADS}" \
      -M 4000 \
      > "${OUT}.log" 2>&1

    N_CLUSTER=$(grep -c '^>' "${OUT}.fa")

    echo "${FAMILY}: ${N_INPUT} contigs -> ${N_CLUSTER} clusters"
done

echo
echo "ALL FAMILY CLUSTERING COMPLETE"
