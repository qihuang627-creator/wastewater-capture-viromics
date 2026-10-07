#!/usr/bin/env bash
set -euo pipefail

RUN="${1:?Usage: $0 RUN}"
THREADS="${THREADS:-8}"

R1="results/per_sample/host_removal/${RUN}_R1.nonhuman.fastq.gz"
R2="results/per_sample/host_removal/${RUN}_R2.nonhuman.fastq.gz"

FA="results/multisample/target_candidates/${RUN}/${RUN}_strict_target_candidates.fa"
META="results/multisample/target_candidates/${RUN}/${RUN}_strict_target_candidates.tsv"

OUT="results/multisample/target_readback/${RUN}"
IDX="${OUT}/index/candidates"

mkdir -p "${OUT}/index" logs/multisample

BAM="${OUT}/${RUN}_target_readback.bam"
IDXSTATS="${OUT}/${RUN}_idxstats.tsv"
DEPTH="${OUT}/${RUN}_depth.tsv"
SUMMARY="${OUT}/${RUN}_target_support.tsv"

if [[ -s "${SUMMARY}" \
   && -s "${DEPTH}" \
   && -s "${IDXSTATS}" \
   && -s "${BAM}" \
   && -s "${BAM}.bai" ]]; then
    echo "Existing complete target read-back outputs found — skipping ${RUN}"
    exit 0
fi

echo
echo "=========================================="
echo "TARGET READ-BACK: ${RUN}"
echo "=========================================="

bowtie2-build \
  "${FA}" \
  "${IDX}" \
  > "logs/multisample/${RUN}_target_bowtie2build.log" 2>&1

bowtie2 \
  --very-sensitive \
  -x "${IDX}" \
  -1 "${R1}" \
  -2 "${R2}" \
  -p "${THREADS}" \
  2> "logs/multisample/${RUN}_target_readback_bowtie2.log" \
| samtools sort \
    -@ "${THREADS}" \
    -o "${BAM}" -

samtools index "${BAM}"

samtools idxstats \
  "${BAM}" \
  > "${IDXSTATS}"

samtools depth \
  -aa \
  -Q 10 \
  "${BAM}" \
  > "${DEPTH}"

python3 - \
  "${RUN}" \
  "${META}" \
  "${IDXSTATS}" \
  "${DEPTH}" \
  "${SUMMARY}" <<'PY'
import csv
import sys
from collections import defaultdict

run, meta_fn, idx_fn, depth_fn, out_fn = sys.argv[1:]

meta = []

with open(meta_fn) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        meta.append(r)

mapped = {}

with open(idx_fn) as f:
    for line in f:
        name, length, nmap, nunmap = line.rstrip().split("\t")
        mapped[name] = int(nmap)

covered = defaultdict(int)
depth_sum = defaultdict(int)

with open(depth_fn) as f:
    for line in f:
        chrom, pos, depth = line.rstrip().split("\t")
        d = int(depth)

        depth_sum[chrom] += d

        if d >= 1:
            covered[chrom] += 1

rows = []
npass = 0

for r in meta:

    name = r["seq_name"]
    length = int(r["length"])

    breadth = (
        100.0 * covered[name] / length
        if length else 0.0
    )

    mean_depth = (
        depth_sum[name] / length
        if length else 0.0
    )

    passed = (
        breadth >= 95.0
        and mean_depth >= 10.0
    )

    if passed:
        npass += 1

    rows.append({
        "sample": run,
        "contig": name,
        "family": r["family"],
        "length": length,
        "mapped_reads": mapped.get(name, 0),
        "breadth1_pct": f"{breadth:.2f}",
        "mean_depth": f"{mean_depth:.2f}",
        "assembly_supported_pass": "YES" if passed else "NO",
    })

with open(out_fn, "w") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "sample",
            "contig",
            "family",
            "length",
            "mapped_reads",
            "breadth1_pct",
            "mean_depth",
            "assembly_supported_pass",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)

print(
    f"{run}: "
    f"{npass}/{len(rows)} strict target contigs "
    f"pass breadth>=95% and mean_depth>=10x"
)
PY

echo "Completed: ${RUN}"
