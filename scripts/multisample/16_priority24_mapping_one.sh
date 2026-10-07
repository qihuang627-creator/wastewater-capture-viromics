#!/usr/bin/env bash
set -euo pipefail

RUN="${1:?Usage: $0 RUN}"
THREADS="${THREADS:-8}"

R1="results/per_sample/host_removal/${RUN}_R1.nonhuman.fastq.gz"
R2="results/per_sample/host_removal/${RUN}_R2.nonhuman.fastq.gz"

INDEX="refs/multisample/priority_clusters/bowtie2/priority24"
REF="refs/multisample/priority_clusters/priority24_cluster_representatives.fa"

OUTDIR="results/multisample/priority_mapping/${RUN}"
LOGDIR="logs/multisample"

mkdir -p "${OUTDIR}" "${LOGDIR}"

BAM="${OUTDIR}/${RUN}_priority24.sorted.bam"
DEPTH="${OUTDIR}/${RUN}_priority24.mapq20.depth.tsv"
SUMMARY="${OUTDIR}/${RUN}_priority24_support.tsv"

if [[ -s "${SUMMARY}" ]]; then
    echo "Existing summary found — skipping ${RUN}"
    exit 0
fi

if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
    echo "ERROR: missing input FASTQ for ${RUN}" >&2
    exit 1
fi

echo
echo "=========================================="
echo "PRIORITY24 MAPPING: ${RUN}"
echo "=========================================="

# Count input fragments (paired reads)
INPUT_FRAGMENTS=$(
    gzip -dc "${R1}" \
    | awk 'END {printf "%.0f\n", NR/4}'
)

echo "Input fragments: ${INPUT_FRAGMENTS}"

# Competitive mapping
bowtie2 \
  --very-sensitive \
  --no-mixed \
  --no-discordant \
  -x "${INDEX}" \
  -1 "${R1}" \
  -2 "${R2}" \
  -p "${THREADS}" \
  2> "${LOGDIR}/${RUN}_priority24_bowtie2.log" \
| samtools sort \
    -@ "${THREADS}" \
    -o "${BAM}" -

samtools index "${BAM}"

# MAPQ >=20 depth
samtools depth \
  -aa \
  -Q 20 \
  "${BAM}" \
  > "${DEPTH}"

python3 - \
  "${RUN}" \
  "${INPUT_FRAGMENTS}" \
  "${REF}" \
  "${BAM}" \
  "${DEPTH}" \
  "${SUMMARY}" <<'PY'
import sys
import subprocess
from collections import defaultdict

run, input_fragments, ref_fn, bam_fn, depth_fn, out_fn = sys.argv[1:]

input_fragments = int(input_fragments)

# -----------------------------
# Reference lengths
# -----------------------------
lengths = {}
name = None
seq = []

with open(ref_fn) as f:
    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if name is not None:
                lengths[name] = len("".join(seq))

            name = line[1:].split()[0]
            seq = []

        else:
            seq.append(line)

if name is not None:
    lengths[name] = len("".join(seq))

# -----------------------------
# MAPQ20 depth
# -----------------------------
covered1 = defaultdict(int)
covered10 = defaultdict(int)
depth_sum = defaultdict(int)

with open(depth_fn) as f:
    for line in f:
        chrom, pos, depth = line.rstrip().split("\t")
        d = int(depth)

        depth_sum[chrom] += d

        if d >= 1:
            covered1[chrom] += 1

        if d >= 10:
            covered10[chrom] += 1

# -----------------------------
# Properly paired MAPQ20 fragments
#
# -f 66:
#   paired + properly paired + read1
#
# -F 2304:
#   exclude secondary + supplementary
#
# Counting read1 means one count per fragment.
# -----------------------------
def fragment_count(ref):
    result = subprocess.run(
        [
            "samtools", "view",
            "-c",
            "-q", "20",
            "-f", "66",
            "-F", "2304",
            bam_fn,
            ref,
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    return int(result.stdout.strip())

# -----------------------------
# Output
# -----------------------------
with open(out_fn, "w") as out:

    out.write(
        "sample\tcluster_id\tlength\t"
        "mapq20_fragments\tFPM\t"
        "breadth1_mapq20_pct\t"
        "breadth10_mapq20_pct\t"
        "mean_depth_mapq20\n"
    )

    for ref, length in lengths.items():

        fragments = fragment_count(ref)

        fpm = (
            fragments / input_fragments * 1_000_000
            if input_fragments else 0
        )

        breadth1 = (
            covered1[ref] / length * 100
            if length else 0
        )

        breadth10 = (
            covered10[ref] / length * 100
            if length else 0
        )

        mean_depth = (
            depth_sum[ref] / length
            if length else 0
        )

        out.write(
            f"{run}\t"
            f"{ref}\t"
            f"{length}\t"
            f"{fragments}\t"
            f"{fpm:.2f}\t"
            f"{breadth1:.2f}\t"
            f"{breadth10:.2f}\t"
            f"{mean_depth:.2f}\n"
        )

print(f"Completed: {run}")
PY

echo
column -t "${SUMMARY}"
