#!/usr/bin/env bash
set -euo pipefail

RUN="${1:?Usage: $0 RUN}"
THREADS="${THREADS:-8}"

MEGAHIT="/media/desk16/iy19302/miniconda3/envs/metagenome/bin/megahit"

R1="results/per_sample/host_removal/${RUN}_R1.nonhuman.fastq.gz"
R2="results/per_sample/host_removal/${RUN}_R2.nonhuman.fastq.gz"

OUTDIR="results/per_sample/assembly/${RUN}_megahit"
SCREEN="results/per_sample/viral_screen"
LOGDIR="logs/multisample"

mkdir -p "${SCREEN}" "${LOGDIR}"

FINAL="${OUTDIR}/final.contigs.fa"
FA1K="${SCREEN}/${RUN}_contigs_1kb.fa"

if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
    echo "ERROR: missing nonhuman FASTQ for ${RUN}" >&2
    exit 1
fi

echo
echo "=========================================="
echo "ASSEMBLY: ${RUN}"
echo "=========================================="

if [[ -s "${FINAL}" ]]; then

    echo "Existing MEGAHIT assembly found — skipping"

else

    rm -rf "${OUTDIR}"

    "${MEGAHIT}" \
      -1 "${R1}" \
      -2 "${R2}" \
      -o "${OUTDIR}" \
      --min-contig-len 500 \
      --k-list 21,29,39,59,79,99,119,141 \
      -t "${THREADS}" \
      > "${LOGDIR}/${RUN}_megahit.log" 2>&1
fi

# Extract contigs >=1 kb
awk '
    /^>/ {
        if (seq != "" && length(seq) >= 1000) {
            print header
            print seq
        }
        header=$0
        seq=""
        next
    }

    {
        seq = seq $0
    }

    END {
        if (seq != "" && length(seq) >= 1000) {
            print header
            print seq
        }
    }
' "${FINAL}" > "${FA1K}"

echo "Assembly complete: ${RUN}"
