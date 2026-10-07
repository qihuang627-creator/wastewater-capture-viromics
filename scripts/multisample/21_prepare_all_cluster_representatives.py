#!/usr/bin/env python3

import csv
from pathlib import Path

clusters_file = Path(
    "results_summary/multisample/"
    "supported_sequence_clusters.tsv"
)

search_root = Path(
    "results/multisample/clustering"
)

outfile = search_root / "all96_cluster_representatives.fa"

with clusters_file.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

print("Clusters in table:", len(rows))

rep_to_cluster = {
    r["representative"]: r["cluster_id"]
    for r in rows
}

wanted = set(rep_to_cluster)

fasta_files = []

for pattern in ("*.fa", "*.fasta", "*.fna"):
    fasta_files.extend(search_root.rglob(pattern))

fasta_files = sorted(set(fasta_files))

seqs = {}

def save_record(name, seq):
    if name is None:
        return

    if name not in wanted:
        return

    seq = "".join(seq).upper()

    if name in seqs and seqs[name] != seq:
        raise SystemExit(
            f"Conflicting sequences found for {name}"
        )

    seqs[name] = seq


for fn in fasta_files:

    if fn == outfile:
        continue

    name = None
    buf = []

    with fn.open() as f:

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                save_record(name, buf)

                name = line[1:].split()[0]
                buf = []

            else:
                buf.append(line)

        save_record(name, buf)

missing = sorted(wanted - set(seqs))

if missing:
    print("\nMissing representatives:")
    for x in missing:
        print(x)

    raise SystemExit(
        f"\nERROR: found {len(seqs)}/{len(wanted)} representatives"
    )

with outfile.open("w") as out:

    for r in rows:
        cid = r["cluster_id"]
        rep = r["representative"]

        out.write(f">{cid}\n")

        seq = seqs[rep]

        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

print(f"Representatives found: {len(seqs)}")
print(f"Saved: {outfile}")
