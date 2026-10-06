#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-ERR14788990}"
REFID="${2:?Usage: $0 RUN REFID}"
THREADS="${THREADS:-8}"
MIN_DP="${MIN_DP:-10}"

PANEL="refs/target_panel/target_references.fa"
MANIFEST="refs/target_panel/reference_manifest.tsv"

COMP_BAM="results/reference_mapping/${RUN}/${RUN}_target_panel.sorted.bam"

OUTDIR="results/consensus/${RUN}/${REFID}"
mkdir -p "${OUTDIR}" logs

# --------------------------------------------------
# 1. Resolve accession / target and extract reference
# --------------------------------------------------

ACC=$(awk -F'\t' -v r="${REFID}" 'NR>1 && $1==r {print $2}' "${MANIFEST}")
TARGET=$(awk -F'\t' -v r="${REFID}" 'NR>1 && $1==r {print $3}' "${MANIFEST}")

if [[ -z "${ACC}" ]]; then
    echo "ERROR: ${REFID} not found in manifest" >&2
    exit 1
fi

REF="${OUTDIR}/${REFID}.fa"

python3 - "${PANEL}" "${REFID}" "${REF}" <<'PY'
import sys

panel, target, out = sys.argv[1:]

keep = False
found = False

with open(panel) as fi, open(out, "w") as fo:
    for line in fi:
        if line.startswith(">"):
            name = line[1:].split("|")[0]
            keep = (name == target)
            if keep:
                fo.write(f">{target}\n")
                found = True
        elif keep:
            fo.write(line)

if not found:
    raise SystemExit(f"Reference {target} not found")
PY

samtools faidx "${REF}"

# --------------------------------------------------
# 2. Select reads uniquely assigned in panel mapping
# --------------------------------------------------

PANEL_REF="${REFID}|${ACC}|${TARGET}"
NAMES="${OUTDIR}/${REFID}.MQ20.readnames.txt"

samtools view \
  -q 20 \
  -F 2308 \
  "${COMP_BAM}" \
  "${PANEL_REF}" \
| cut -f1 \
| sort -u \
> "${NAMES}"

echo "Unique read names selected: $(wc -l < "${NAMES}")"

# Pull both mates belonging to selected read names
samtools view \
  -b \
  -N "${NAMES}" \
  "${COMP_BAM}" \
| samtools sort \
    -n \
    -@ "${THREADS}" \
    -o "${OUTDIR}/${REFID}.selected.namesort.bam" -

samtools fastq \
  -@ "${THREADS}" \
  -n \
  -1 "${OUTDIR}/${REFID}.R1.fastq.gz" \
  -2 "${OUTDIR}/${REFID}.R2.fastq.gz" \
  -0 /dev/null \
  -s /dev/null \
  "${OUTDIR}/${REFID}.selected.namesort.bam" \
  2> "logs/${RUN}_${REFID}_selected_fastq.log"

# --------------------------------------------------
# 3. Remap selected reads to the single reference
# --------------------------------------------------

bowtie2-build \
  "${REF}" \
  "${OUTDIR}/${REFID}_index" \
  > "logs/${RUN}_${REFID}_index.log" 2>&1

bowtie2 \
  --very-sensitive \
  -x "${OUTDIR}/${REFID}_index" \
  -1 "${OUTDIR}/${REFID}.R1.fastq.gz" \
  -2 "${OUTDIR}/${REFID}.R2.fastq.gz" \
  -p "${THREADS}" \
  2> "logs/${RUN}_${REFID}_remap.log" \
| samtools view -@ "${THREADS}" -b - \
| samtools sort \
    -@ "${THREADS}" \
    -o "${OUTDIR}/${REFID}.remap.sorted.bam" -

samtools index "${OUTDIR}/${REFID}.remap.sorted.bam"

# --------------------------------------------------
# 4. High-quality depth
# --------------------------------------------------

samtools depth \
  -aa \
  -q 20 \
  -Q 20 \
  "${OUTDIR}/${REFID}.remap.sorted.bam" \
  > "${OUTDIR}/${REFID}.MQ20.depth.tsv"

# --------------------------------------------------
# 5. Haploid variant calling
# --------------------------------------------------

bcftools mpileup \
  -d 100000 \
  -Ou \
  -f "${REF}" \
  -q 20 \
  -Q 20 \
  -a FORMAT/DP,FORMAT/AD \
  "${OUTDIR}/${REFID}.remap.sorted.bam" \
| bcftools call \
    --ploidy 1 \
    -mv \
    -Oz \
    -o "${OUTDIR}/${REFID}.raw.vcf.gz"

bcftools index -f "${OUTDIR}/${REFID}.raw.vcf.gz"

bcftools filter \
  -i 'QUAL>=20 && FORMAT/DP>=10' \
  -Oz \
  -o "${OUTDIR}/${REFID}.filtered.vcf.gz" \
  "${OUTDIR}/${REFID}.raw.vcf.gz"

bcftools index -f "${OUTDIR}/${REFID}.filtered.vcf.gz"

# --------------------------------------------------
# 6. Apply variants to reference
# --------------------------------------------------

bcftools consensus \
  -f "${REF}" \
  "${OUTDIR}/${REFID}.filtered.vcf.gz" \
  > "${OUTDIR}/${REFID}.consensus.unmasked.fa"

# --------------------------------------------------
# 7. Mask positions with DP < MIN_DP and summarize
# --------------------------------------------------

python3 - \
  "${OUTDIR}/${REFID}.consensus.unmasked.fa" \
  "${OUTDIR}/${REFID}.MQ20.depth.tsv" \
  "${OUTDIR}/${REFID}.consensus.fa" \
  "${MIN_DP}" \
  "${OUTDIR}/${REFID}.filtered.vcf.gz" <<'PY'
import gzip
import statistics
import subprocess
import sys

fa, depth_file, out, min_dp, vcf = sys.argv[1:]
min_dp = int(min_dp)

seq = []
with open(fa) as f:
    for line in f:
        if not line.startswith(">"):
            seq.append(line.strip())

seq = list("".join(seq))

depths = []
with open(depth_file) as f:
    for line in f:
        _, _, d = line.rstrip().split("\t")
        depths.append(int(d))

if len(depths) != len(seq):
    raise SystemExit(
        f"Depth length {len(depths)} != reference length {len(seq)}"
    )

for i, d in enumerate(depths):
    if d < min_dp:
        seq[i] = "N"

with open(out, "w") as f:
    f.write(">consensus\n")
    s = "".join(seq)
    for i in range(0, len(s), 80):
        f.write(s[i:i+80] + "\n")

L = len(depths)

def breadth(x):
    return 100 * sum(d >= x for d in depths) / L

variant_count = int(
    subprocess.check_output(
        ["bcftools", "view", "-H", vcf],
        text=True
    ).count("\n")
)

print(f"Reference length: {L:,}")
print(f"Breadth >=1x:  {breadth(1):.2f}%")
print(f"Breadth >=5x:  {breadth(5):.2f}%")
print(f"Breadth >=10x: {breadth(10):.2f}%")
print(f"Breadth >=20x: {breadth(20):.2f}%")
print(f"Mean depth:     {statistics.mean(depths):.2f}x")
print(f"Median depth:   {statistics.median(depths):.2f}x")
print(f"Masked N:       {sum(x == 'N' for x in seq):,}")
print(f"Callable pct:   {100*sum(x != 'N' for x in seq)/L:.2f}%")
print(f"PASS variants:  {variant_count}")
PY

echo
echo "=== Completed ${REFID} ==="
