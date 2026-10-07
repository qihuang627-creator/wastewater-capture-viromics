#!/usr/bin/env bash
set -euo pipefail

MANIFEST="refs/single_sample/target_panel/reference_manifest.tsv"
OUTDIR="refs/single_sample/target_panel"
RAWDIR="${OUTDIR}/raw"
COMBINED="${OUTDIR}/target_references.fa"
SUMMARY="${OUTDIR}/reference_sequences.tsv"

mkdir -p "${RAWDIR}"

# Start fresh
: > "${COMBINED}"
echo -e "reference_id\taccession\ttarget\tlength_bp" > "${SUMMARY}"

tail -n +2 "${MANIFEST}" |
while IFS=$'\t' read -r REFID ACC TARGET
do
    echo "=== Fetching ${REFID} (${ACC}) ==="

    RAW="${RAWDIR}/${ACC}.fa"

    curl -fsSL \
      --retry 4 \
      --retry-delay 2 \
      --connect-timeout 20 \
      --max-time 120 \
      -G "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
      --data-urlencode "db=nuccore" \
      --data-urlencode "id=${ACC}" \
      --data-urlencode "rettype=fasta" \
      --data-urlencode "retmode=text" \
      -o "${RAW}"

    if ! grep -q '^>' "${RAW}"; then
        echo "ERROR: no FASTA record returned for ${ACC}" >&2
        exit 1
    fi

    NREC=$(grep -c '^>' "${RAW}")
    if [[ "${NREC}" -ne 1 ]]; then
        echo "ERROR: expected 1 FASTA record for ${ACC}, got ${NREC}" >&2
        exit 1
    fi

    LEN=$(
      awk '
        /^>/ {next}
        {gsub(/[[:space:]]/, ""); n += length($0)}
        END {print n+0}
      ' "${RAW}"
    )

    if [[ "${LEN}" -le 0 ]]; then
        echo "ERROR: zero-length sequence for ${ACC}" >&2
        exit 1
    fi

    # Use clean, unique IDs for downstream competitive mapping
    awk \
      -v ref="${REFID}" \
      -v acc="${ACC}" \
      -v target="${TARGET}" '
        /^>/ {
            print ">" ref "|" acc "|" target
            next
        }
        {print}
      ' "${RAW}" >> "${COMBINED}"

    echo -e "${REFID}\t${ACC}\t${TARGET}\t${LEN}" >> "${SUMMARY}"

    echo "Downloaded: ${LEN} bp"

    # Stay below NCBI unauthenticated request-rate limits
    sleep 0.5
done

echo
echo "=== Reference panel complete ==="
echo "Sequences: $(grep -c '^>' "${COMBINED}")"

if command -v column >/dev/null 2>&1; then
    column -t "${SUMMARY}"
else
    cat "${SUMMARY}"
fi
