#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

cluster_dir = Path(
    "results/multisample/clustering/clusters"
)

manifest_file = Path(
    "results_summary/multisample/"
    "supported_contig_manifest.tsv"
)

out = Path(
    "results_summary/multisample/"
    "supported_sequence_clusters.tsv"
)

if not manifest_file.exists():
    raise SystemExit(
        f"Missing manifest: {manifest_file}"
    )

with manifest_file.open() as f:
    manifest = {
        r["unique_id"]: r
        for r in csv.DictReader(f, delimiter="\t")
    }


def parse_clstr(fn):

    clusters = []
    current = None

    with fn.open() as f:

        for line in f:

            line = line.strip()

            if line.startswith(">Cluster"):

                if current is not None:
                    clusters.append(current)

                current = {
                    "cluster": line.split()[-1],
                    "members": [],
                    "representative": None,
                }

                continue

            if ">" not in line:
                continue

            seq_id = (
                line.split(">", 1)[1]
                .split("...", 1)[0]
            )

            current["members"].append(seq_id)

            if line.endswith("*"):
                current["representative"] = seq_id

    if current is not None:
        clusters.append(current)

    return clusters


rows = []

for clstr in sorted(
    cluster_dir.glob("*_ani95.fa.clstr")
):

    family = clstr.name.replace(
        "_ani95.fa.clstr", ""
    )

    clusters = parse_clstr(clstr)

    for c in clusters:

        rep = c["representative"]

        if rep is None:
            raise RuntimeError(
                f"No representative in "
                f"{clstr}, cluster {c['cluster']}"
            )

        members = c["members"]

        for member in members:
            if member not in manifest:
                raise RuntimeError(
                    f"{member} missing from manifest"
                )

        samples = sorted({
            manifest[x]["sample"]
            for x in members
        })

        rep_meta = manifest[rep]

        rows.append({
            "family": family,
            "cluster_id":
                f"{family}_C{int(c['cluster']) + 1:03d}",
            "representative": rep,
            "representative_length":
                rep_meta["length"],
            "cluster_size":
                len(members),
            "sample_count":
                len(samples),
            "samples":
                ",".join(samples),
            "members":
                ",".join(members),
        })


rows.sort(
    key=lambda r: (
        r["family"],
        -int(r["sample_count"]),
        -int(r["cluster_size"]),
        -int(r["representative_length"]),
    )
)

fields = [
    "family",
    "cluster_id",
    "representative",
    "representative_length",
    "cluster_size",
    "sample_count",
    "samples",
    "members",
]

with out.open("w") as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(rows)


family_clusters = defaultdict(int)
family_cross = defaultdict(int)
family_contigs = defaultdict(int)

for x in manifest.values():
    family_contigs[x["family"]] += 1

for r in rows:

    family_clusters[r["family"]] += 1

    if int(r["sample_count"]) >= 2:
        family_cross[r["family"]] += 1


print(f"Total sequence clusters: {len(rows)}")
print()

print(
    "family\t"
    "supported_contigs\t"
    "sequence_clusters\t"
    "cross_sample_clusters"
)

for family in sorted(family_clusters):

    print(
        f"{family}\t"
        f"{family_contigs[family]}\t"
        f"{family_clusters[family]}\t"
        f"{family_cross[family]}"
    )
