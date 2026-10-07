#!/usr/bin/env bash
set -euo pipefail

RUN="${1:?Usage: $0 RUN}"
THREADS="${THREADS:-8}"

R1="results/per_sample/host_removal/${RUN}_R1.nonhuman.fastq.gz"
R2="results/per_sample/host_removal/${RUN}_R2.nonhuman.fastq.gz"

INDEX="refs/multisample/all96/bowtie2/all96"
REF="refs/multisample/all96/all96_cluster_representatives.fa"

OUTDIR="results/multisample/all96_mapping/${RUN}"
LOGDIR="logs/multisample"

mkdir -p "${OUTDIR}" "${LOGDIR}"

BAM="${OUTDIR}/${RUN}_all96.sorted.bam"
DEPTH="${OUTDIR}/${RUN}_all96.mapq20.depth.tsv"
FRAGS="${OUTDIR}/${RUN}_all96.mapq20.fragments.tsv"
SUMMARY="${OUTDIR}/${RUN}_all96_support.tsv"

if [[ -s "${SUMMARY}"    && -s "${DEPTH}"    && -s "${FRAGS}"    && -s "${BAM}"    && -s "${BAM}.bai" ]]; then
    echo "Existing complete all96 mapping outputs found — skipping ${RUN}"
    exit 0
fi

INPUT_FRAGMENTS=$(
    gzip -dc "${R1}" |
    awk 'END {printf "%.0f\n", NR/4}'
)

echo "=== ${RUN} ==="
echo "Input fragments: ${INPUT_FRAGMENTS}"

bowtie2 \
  --very-sensitive \
  --no-mixed \
  --no-discordant \
  -x "${INDEX}" \
  -1 "${R1}" \
  -2 "${R2}" \
  -p "${THREADS}" \
  2> "${LOGDIR}/${RUN}_all96_bowtie2.log" \
| samtools sort \
    -@ "${THREADS}" \
    -o "${BAM}" -

samtools index "${BAM}"

samtools depth \
  -aa \
  -Q 20 \
  "${BAM}" \
  > "${DEPTH}"

# one count per properly paired MAPQ20 fragment
samtools view \
  -q 20 \
  -f 66 \
  -F 2304 \
  "${BAM}" \
| awk '{count[$3]++}
       END {
         for (r in count)
           print r "\t" count[r]
       }' \
> "${FRAGS}"

python3 - \
  "${RUN}" \
  "${INPUT_FRAGMENTS}" \
  "${REF}" \
  "${DEPTH}" \
  "${FRAGS}" \
  "${SUMMARY}" <<'PY'
import sys
from collections import defaultdict

run, nfrag, ref_fn, depth_fn, frag_fn, out_fn = sys.argv[1:]
nfrag = int(nfrag)

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

covered1 = defaultdict(int)
covered10 = defaultdict(int)
depth_sum = defaultdict(int)

with open(depth_fn) as f:
    for line in f:
        ref, pos, depth = line.rstrip().split("\t")
        d = int(depth)

        depth_sum[ref] += d

        if d >= 1:
            covered1[ref] += 1

        if d >= 10:
            covered10[ref] += 1

fragments = defaultdict(int)

with open(frag_fn) as f:
    for line in f:
        ref, n = line.rstrip().split("\t")
        fragments[ref] = int(n)

with open(out_fn, "w") as out:

    out.write(
        "sample\tcluster_id\tlength\t"
        "mapq20_fragments\tFPM\t"
        "breadth1_mapq20_pct\t"
        "breadth10_mapq20_pct\t"
        "mean_depth_mapq20\n"
    )

    for ref, length in lengths.items():

        nf = fragments[ref]

        fpm = (
            nf / nfrag * 1_000_000
            if nfrag else 0
        )

        b1 = covered1[ref] / length * 100
        b10 = covered10[ref] / length * 100
        md = depth_sum[ref] / length

        out.write(
            f"{run}\t{ref}\t{length}\t"
            f"{nf}\t{fpm:.2f}\t"
            f"{b1:.2f}\t{b10:.2f}\t{md:.2f}\n"
        )

print(f"Completed: {run}")
PY
