#!/usr/bin/env bash
set -euo pipefail

RUNLIST="${1:-config/multisample/run_lists/new_runs.txt}"
OUTDIR="${2:-raw}"

mkdir -p "${OUTDIR}" logs/multisample

while read -r RUN
do
    [[ -z "${RUN}" ]] && continue

    echo
    echo "=========================================="
    echo "Downloading ${RUN}"
    echo "=========================================="

    META=$(mktemp)

    curl -fsSL --retry 5 --retry-delay 3 \
      "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv" \
      -o "${META}"

    DATA=$(tail -n +2 "${META}")

    if [[ -z "${DATA}" ]]; then
        echo "ERROR: ENA returned no FASTQ metadata for ${RUN}" >&2
        rm -f "${META}"
        exit 1
    fi

    FASTQ_FTP=$(echo "${DATA}" | cut -f2)
    FASTQ_MD5=$(echo "${DATA}" | cut -f3)

    IFS=';' read -r URL1 URL2 <<< "${FASTQ_FTP}"
    IFS=';' read -r MD51 MD52 <<< "${FASTQ_MD5}"

    if [[ -z "${URL1:-}" || -z "${URL2:-}" ]]; then
        echo "ERROR: ${RUN} does not appear to contain paired FASTQ files" >&2
        rm -f "${META}"
        exit 1
    fi

    URL1="https://${URL1}"
    URL2="https://${URL2}"

    OUT1="${OUTDIR}/${RUN}_1.fastq.gz"
    OUT2="${OUTDIR}/${RUN}_2.fastq.gz"

    echo "R1: ${URL1}"
    echo "R2: ${URL2}"

    curl -fL \
      -C - \
      --retry 10 \
      --retry-delay 5 \
      --retry-all-errors \
      --connect-timeout 30 \
      "${URL1}" \
      -o "${OUT1}"

    curl -fL \
      -C - \
      --retry 10 \
      --retry-delay 5 \
      --retry-all-errors \
      --connect-timeout 30 \
      "${URL2}" \
      -o "${OUT2}"

    echo
    echo "Checking gzip integrity..."

    gzip -t "${OUT1}"
    gzip -t "${OUT2}"

    echo "Checking MD5..."

    GOT1=$(md5sum "${OUT1}" | awk '{print $1}')
    GOT2=$(md5sum "${OUT2}" | awk '{print $1}')

    if [[ "${GOT1}" != "${MD51}" ]]; then
        echo "ERROR: ${RUN} R1 MD5 mismatch"
        echo "Expected: ${MD51}"
        echo "Observed: ${GOT1}"
        exit 1
    fi

    if [[ "${GOT2}" != "${MD52}" ]]; then
        echo "ERROR: ${RUN} R2 MD5 mismatch"
        echo "Expected: ${MD52}"
        echo "Observed: ${GOT2}"
        exit 1
    fi

    PAIRS1=$(zcat "${OUT1}" | awk 'END{print NR/4}')
    PAIRS2=$(zcat "${OUT2}" | awk 'END{print NR/4}')

    if [[ "${PAIRS1}" != "${PAIRS2}" ]]; then
        echo "ERROR: R1/R2 read counts differ for ${RUN}" >&2
        exit 1
    fi

    echo "${RUN}: ${PAIRS1} read pairs"
    echo "${RUN}: DOWNLOAD_OK"

    rm -f "${META}"

done < "${RUNLIST}"

echo
echo "All requested FASTQs downloaded and verified."
