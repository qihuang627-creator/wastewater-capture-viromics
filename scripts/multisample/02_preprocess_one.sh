#!/usr/bin/env bash
set -euo pipefail

RUN="${1:?Usage: $0 RUN}"
THREADS="${THREADS:-8}"

RAW1="raw/${RUN}_1.fastq.gz"
RAW2="raw/${RUN}_2.fastq.gz"

CDHIT="${HOME}/bin/cd-hit-dup"
FASTP="/media/desk16/iy19302/miniconda3/envs/metagenome/bin/fastp"
HUMAN_INDEX="${HOME}/QI/wastewater-shotgun/databases/human_GRCh38/GRCh38"

DEDUPDIR="results/per_sample/dedup/u150"
QCDIR="qc/per_sample/clean"
HOSTDIR="results/per_sample/host_removal"
REPORTDIR="results_summary/multisample/preprocessing"
LOGDIR="logs/multisample"

mkdir -p \
  "${DEDUPDIR}" \
  "${QCDIR}" \
  "${HOSTDIR}" \
  "${REPORTDIR}" \
  "${LOGDIR}"

for F in "${RAW1}" "${RAW2}"
do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: missing ${F}" >&2
        exit 1
    fi
done

if [[ ! -x "${CDHIT}" ]]; then
    echo "ERROR: cd-hit-dup not found at ${CDHIT}" >&2
    exit 1
fi

# -----------------------------
# Helpers
# -----------------------------

count_pairs () {
    local f="$1"
    gzip -dc "$f" | awk 'END {printf "%.0f\n", NR/4}'
}

RAW_PAIRS=$(count_pairs "${RAW1}")

echo
echo "=========================================="
echo "Sample: ${RUN}"
echo "Raw pairs: ${RAW_PAIRS}"
echo "=========================================="

# -----------------------------
# 1. Putative PCR deduplication
# -----------------------------

D1="${DEDUPDIR}/${RUN}_R1.dedup.fastq.gz"
D2="${DEDUPDIR}/${RUN}_R2.dedup.fastq.gz"

if [[ -s "${D1}" && -s "${D2}" ]]; then

    echo "[1/3] Dedup already exists — skipping"

else

    echo "[1/3] CD-HIT-DUP u150"

    TMP=$(mktemp -d)

    cleanup () {
        rm -rf "${TMP}"
    }

    trap cleanup EXIT

    gzip -dc "${RAW1}" > "${TMP}/R1.fastq"
    gzip -dc "${RAW2}" > "${TMP}/R2.fastq"

    "${CDHIT}" \
      -i "${TMP}/R1.fastq" \
      -i2 "${TMP}/R2.fastq" \
      -o "${TMP}/R1.dedup.fastq" \
      -o2 "${TMP}/R2.dedup.fastq" \
      -u 150 \
      > "${LOGDIR}/${RUN}_cdhitdup.log" 2>&1

    gzip -c "${TMP}/R1.dedup.fastq" > "${D1}"
    gzip -c "${TMP}/R2.dedup.fastq" > "${D2}"

    gzip -t "${D1}"
    gzip -t "${D2}"

    cleanup
    trap - EXIT
fi

DEDUP_PAIRS=$(count_pairs "${D1}")

# -----------------------------
# 2. fastp
# -----------------------------

C1="${QCDIR}/${RUN}_R1.clean.fastq.gz"
C2="${QCDIR}/${RUN}_R2.clean.fastq.gz"

if [[ -s "${C1}" && -s "${C2}" ]]; then

    echo "[2/3] fastp output already exists — skipping"

else

    echo "[2/3] fastp"

    "${FASTP}" \
      -i "${D1}" \
      -I "${D2}" \
      -o "${C1}" \
      -O "${C2}" \
      --thread 4 \
      --detect_adapter_for_pe \
      --cut_right \
      --cut_window_size 4 \
      --cut_mean_quality 20 \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 40 \
      --n_base_limit 5 \
      --length_required 50 \
      --json "${QCDIR}/${RUN}.fastp.json" \
      --html "${QCDIR}/${RUN}.fastp.html" \
      > "${LOGDIR}/${RUN}_fastp.log" 2>&1
fi

QC_PAIRS=$(count_pairs "${C1}")

# -----------------------------
# 3. Human removal
# -----------------------------

NH1="${HOSTDIR}/${RUN}_R1.nonhuman.fastq.gz"
NH2="${HOSTDIR}/${RUN}_R2.nonhuman.fastq.gz"

if [[ -s "${NH1}" && -s "${NH2}" ]]; then

    echo "[3/3] Nonhuman FASTQ already exists — skipping"

else

    echo "[3/3] GRCh38 removal"

    TMP=$(mktemp -d)

    cleanup () {
        rm -rf "${TMP}"
    }

    trap cleanup EXIT

    bowtie2 \
      --very-sensitive \
      -x "${HUMAN_INDEX}" \
      -1 "${C1}" \
      -2 "${C2}" \
      -p "${THREADS}" \
      2> "${LOGDIR}/${RUN}_human_bowtie2.log" \
    | samtools view \
        -@ "${THREADS}" \
        -b \
        -o "${TMP}/human_screen.bam" -

    samtools view \
      -@ "${THREADS}" \
      -b \
      -f 12 \
      -F 2304 \
      "${TMP}/human_screen.bam" \
    | samtools sort \
        -@ "${THREADS}" \
        -n \
        -o "${TMP}/both_unmapped.namesort.bam" -

    samtools fastq \
      -@ "${THREADS}" \
      -n \
      -1 "${NH1}" \
      -2 "${NH2}" \
      -0 /dev/null \
      -s /dev/null \
      "${TMP}/both_unmapped.namesort.bam" \
      2> "${LOGDIR}/${RUN}_human_fastq.log"

    gzip -t "${NH1}"
    gzip -t "${NH2}"

    cleanup
    trap - EXIT
fi

NONHUMAN_PAIRS=$(count_pairs "${NH1}")

# -----------------------------
# Summary
# -----------------------------

HOST_RATE=$(
    grep -oE '[0-9.]+% overall alignment rate' \
      "${LOGDIR}/${RUN}_human_bowtie2.log" \
    | tail -1 \
    | awk '{print $1}' \
    | tr -d '%' \
    || true
)

HOST_RATE="${HOST_RATE:-NA}"

DEDUP_REMOVED=$(
    awk -v a="${RAW_PAIRS}" -v b="${DEDUP_PAIRS}" \
      'BEGIN {printf "%.3f", 100*(a-b)/a}'
)

QC_RETAINED=$(
    awk -v a="${DEDUP_PAIRS}" -v b="${QC_PAIRS}" \
      'BEGIN {printf "%.3f", 100*b/a}'
)

FINAL_RETAINED=$(
    awk -v a="${RAW_PAIRS}" -v b="${NONHUMAN_PAIRS}" \
      'BEGIN {printf "%.3f", 100*b/a}'
)

SUMMARY="${REPORTDIR}/${RUN}_preprocessing.tsv"

cat > "${SUMMARY}" <<EOF2
sample	raw_pairs	dedup_pairs	qc_pairs	nonhuman_pairs	dedup_removed_pct	qc_retained_pct	human_alignment_rate_pct	final_retained_pct
${RUN}	${RAW_PAIRS}	${DEDUP_PAIRS}	${QC_PAIRS}	${NONHUMAN_PAIRS}	${DEDUP_REMOVED}	${QC_RETAINED}	${HOST_RATE}	${FINAL_RETAINED}
EOF2

echo
echo "=========================================="
echo "Completed: ${RUN}"
echo "=========================================="
column -t "${SUMMARY}"
