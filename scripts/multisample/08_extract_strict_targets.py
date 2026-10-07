#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

run = sys.argv[1]

families = [
    "Adenoviridae",
    "Astroviridae",
    "Caliciviridae",
    "Hepeviridae",
    "Parvoviridae",
    "Picornaviridae",
    "Sedoreoviridae",
    "Spinareoviridae",
]

summary = Path(
    f"results/per_sample/genomad/{run}/"
    f"{run}_contigs_1kb_summary/"
    f"{run}_contigs_1kb_virus_summary.tsv"
)

fasta = Path(
    f"results/per_sample/viral_screen/{run}_contigs_1kb.fa"
)

outdir = Path(
    f"results/multisample/target_candidates/{run}"
)
outdir.mkdir(parents=True, exist_ok=True)

outfa = outdir / f"{run}_strict_target_candidates.fa"
outtsv = outdir / f"{run}_strict_target_candidates.tsv"

selected = {}
rows = []

with summary.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for r in reader:
        tax = r["taxonomy"] or ""

        matches = [
            fam for fam in families
            if fam in tax
        ]

        if not matches:
            continue

        fam = matches[0]
        seq = r["seq_name"]

        selected[seq] = fam

        rows.append({
            "seq_name": seq,
            "family": fam,
            "length": r["length"],
            "taxonomy": tax,
        })

seqs = {}
header = None
parts = []

with fasta.open() as f:
    for line in f:

        if line.startswith(">"):

            if header is not None:
                seqs[header] = "".join(parts)

            header = line[1:].strip().split()[0]
            parts = []

        else:
            parts.append(line.strip())

    if header is not None:
        seqs[header] = "".join(parts)

with outfa.open("w") as f:
    for r in rows:

        name = r["seq_name"]

        if name not in seqs:
            raise SystemExit(
                f"ERROR: {name} missing from {fasta}"
            )

        f.write(f">{name}\n")

        seq = seqs[name]

        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

with outtsv.open("w") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "seq_name",
            "family",
            "length",
            "taxonomy",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)

print(f"{run}\t{len(rows)} strict target contigs")
