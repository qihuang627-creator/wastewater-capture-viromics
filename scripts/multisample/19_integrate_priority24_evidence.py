#!/usr/bin/env python3

import csv
import re
from pathlib import Path

ann_file = Path(
    "results_summary/multisample/"
    "priority_cluster_annotation.tsv"
)

map_file = Path(
    "results_summary/multisample/"
    "priority24_mapping_long.tsv"
)

out_long = Path(
    "results_summary/multisample/"
    "priority24_evidence_integrated.tsv"
)

out_matrix = Path(
    "results_summary/multisample/"
    "priority24_evidence_matrix.tsv"
)

sample_order = [
    "ERR14788990",
    "ERR14788988",
    "ERR14788991",
    "ERR14788992",
    "ERR14788995",
    "ERR14789044",
    "ERR14788864",
    "ERR14788954",
]

# ---------------------------------------
# Load cluster annotations and recover
# samples contributing assembled contigs
# ---------------------------------------
annotations = {}

with ann_file.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):

        # Robustly recover ERR accessions regardless of separator
        assembly_samples = set(
            re.findall(r"ERR\d+", r["samples"])
        )

        r["assembly_samples_set"] = assembly_samples
        annotations[r["cluster_id"]] = r


def classify(assembly_member, breadth):

    # Highest-confidence evidence:
    # sample itself produced a read-back-supported contig
    if assembly_member:
        return "ASSEMBLY"

    # Mapping-only descriptive tiers
    if breadth >= 80:
        return "MAPPING"

    if breadth >= 50:
        return "PARTIAL"

    if breadth >= 10:
        return "LOCAL"

    return "NONE"


rows = []

with map_file.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):

        cid = r["cluster_id"]
        sample = r["sample"]
        breadth = float(r["breadth1_mapq20_pct"])

        a = annotations[cid]

        assembly_member = (
            sample in a["assembly_samples_set"]
        )

        evidence = classify(
            assembly_member,
            breadth
        )

        rows.append({
            "sample": sample,
            "family": r["family"],
            "cluster_id": cid,
            "biological_label": r["biological_label"],
            "virus_group": r["virus_group"],
            "assembly_member":
                "YES" if assembly_member else "NO",
            "evidence_class": evidence,
            "mapq20_fragments": r["mapq20_fragments"],
            "FPM": r["FPM"],
            "breadth1_mapq20_pct":
                r["breadth1_mapq20_pct"],
            "breadth10_mapq20_pct":
                r["breadth10_mapq20_pct"],
            "mean_depth_mapq20":
                r["mean_depth_mapq20"],
        })

# ---------------------------------------
# Long table
# ---------------------------------------
fields = [
    "sample",
    "family",
    "cluster_id",
    "biological_label",
    "virus_group",
    "assembly_member",
    "evidence_class",
    "mapq20_fragments",
    "FPM",
    "breadth1_mapq20_pct",
    "breadth10_mapq20_pct",
    "mean_depth_mapq20",
]

with out_long.open("w") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(rows)

# ---------------------------------------
# Matrix
# ---------------------------------------
cluster_order = []

for r in rows:
    if r["cluster_id"] not in cluster_order:
        cluster_order.append(r["cluster_id"])

lookup = {
    (r["cluster_id"], r["sample"]): r
    for r in rows
}

with out_matrix.open("w") as f:

    f.write(
        "\t".join([
            "cluster_id",
            "biological_label",
            *sample_order
        ]) + "\n"
    )

    for cid in cluster_order:

        label = annotations[cid]["biological_label"]

        vals = [
            lookup[(cid, s)]["evidence_class"]
            for s in sample_order
        ]

        f.write(
            "\t".join(
                [cid, label, *vals]
            ) + "\n"
        )

print(f"Rows: {len(rows)}")
print(f"Saved: {out_long}")
print(f"Saved: {out_matrix}")

# ---------------------------------------
# Summary
# ---------------------------------------
counts = {}

for r in rows:
    key = r["evidence_class"]
    counts[key] = counts.get(key, 0) + 1

print()
print("Evidence class counts:")

for k in [
    "ASSEMBLY",
    "MAPPING",
    "PARTIAL",
    "LOCAL",
    "NONE"
]:
    print(f"{k}: {counts.get(k, 0)}")
