#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

support_file = Path(
    "results_summary/multisample/supported_target_contigs.tsv"
)

outdir = Path("results/multisample/clustering/input")
outdir.mkdir(parents=True, exist_ok=True)

manifest_out = Path(
    "results_summary/multisample/"
    "supported_contig_manifest.tsv"
)

def read_fasta(fn):
    seqs = {}
    name = None
    parts = []

    with fn.open() as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:
                    seqs[name] = "".join(parts)

                name = line[1:].split()[0]
                parts = []

            else:
                parts.append(line)

    if name is not None:
        seqs[name] = "".join(parts)

    return seqs


with support_file.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

by_sample = defaultdict(list)

for r in rows:
    by_sample[r["sample"]].append(r)

family_records = defaultdict(list)
manifest_rows = []

for sample, sample_rows in by_sample.items():

    fa = Path(
        f"results/multisample/target_candidates/{sample}/"
        f"{sample}_strict_target_candidates.fa"
    )

    if not fa.exists():
        raise SystemExit(f"Missing FASTA: {fa}")

    seqs = read_fasta(fa)

    for r in sample_rows:

        contig = r["contig"]
        family = r["family"]

        if contig not in seqs:
            raise SystemExit(
                f"{sample} {contig} missing from {fa}"
            )

        seq = seqs[contig]

        unique_id = f"{sample}|{contig}"

        record = {
            "unique_id": unique_id,
            "sample": sample,
            "contig": contig,
            "family": family,
            "length": len(seq),
            "breadth1_pct": r["breadth1_pct"],
            "mean_depth_mapq10": r["mean_depth_mapq10"],
            "mapped_reads": r["mapped_reads"],
            "sequence": seq,
        }

        family_records[family].append(record)
        manifest_rows.append(record)


# Sort each family by descending contig length.
# CD-HIT therefore preferentially uses a long sequence as representative.
for family, records in family_records.items():

    records.sort(
        key=lambda x: x["length"],
        reverse=True
    )

    outfile = outdir / f"{family}.fa"

    with outfile.open("w") as f:

        for r in records:

            f.write(f">{r['unique_id']}\n")

            seq = r["sequence"]

            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


# Combined FASTA
all_records = sorted(
    manifest_rows,
    key=lambda x: x["length"],
    reverse=True
)

pooled = outdir / "all_supported_targets.fa"

with pooled.open("w") as f:

    for r in all_records:

        f.write(f">{r['unique_id']}\n")

        seq = r["sequence"]

        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")


# Manifest
fields = [
    "unique_id",
    "sample",
    "contig",
    "family",
    "length",
    "breadth1_pct",
    "mean_depth_mapq10",
    "mapped_reads",
]

with manifest_out.open("w") as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )

    w.writeheader()
    w.writerows(manifest_rows)


print(f"Total supported sequences: {len(manifest_rows)}")

for family in sorted(family_records):
    print(
        f"{family}\t"
        f"{len(family_records[family])}"
    )
