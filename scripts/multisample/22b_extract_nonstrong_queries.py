#!/usr/bin/env python3

import csv
from pathlib import Path

ANNOT = Path(
    "results_summary/multisample/"
    "all96_ntviruses_annotation.tsv"
)

FASTA = Path(
    "results/multisample/clustering/"
    "all96_cluster_representatives.fa"
)

OUT = Path(
    "results/multisample/blast/"
    "local_nt_viruses_all96/"
    "nonstrong10.fa"
)

OUT.parent.mkdir(parents=True, exist_ok=True)


def read_fasta(fn):
    seqs = {}
    name = None
    buf = []

    with fn.open() as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:
                    seqs[name] = "".join(buf)

                name = line[1:].split()[0]
                buf = []

            else:
                buf.append(line)

    if name is not None:
        seqs[name] = "".join(buf)

    return seqs


with ANNOT.open() as f:
    rows = list(
        csv.DictReader(f, delimiter="\t")
    )

wanted = [
    r["cluster_id"]
    for r in rows
    if r["match_quality"] != "strong"
]

seqs = read_fasta(FASTA)

missing = [
    x for x in wanted
    if x not in seqs
]

if missing:
    raise SystemExit(
        "Missing representative sequences: "
        + ", ".join(missing)
    )

with OUT.open("w") as f:

    for cid in wanted:

        f.write(f">{cid}\n")

        seq = seqs[cid]

        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")


print(
    f"Non-strong queries: {len(wanted)}"
)

print("Saved:", OUT)
