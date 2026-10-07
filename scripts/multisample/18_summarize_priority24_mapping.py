#!/usr/bin/env python3

import csv
from pathlib import Path

mapping_root = Path(
    "results/multisample/priority_mapping"
)

annotation_file = Path(
    "results_summary/multisample/"
    "priority_cluster_annotation.tsv"
)

out_long = Path(
    "results_summary/multisample/"
    "priority24_mapping_long.tsv"
)

out_breadth = Path(
    "results_summary/multisample/"
    "priority24_breadth1_matrix.tsv"
)

out_fpm = Path(
    "results_summary/multisample/"
    "priority24_FPM_matrix.tsv"
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

# ----------------------------
# Load annotation
# ----------------------------
annotation = {}

with annotation_file.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):
        annotation[r["cluster_id"]] = r

# ----------------------------
# Load mapping results
# ----------------------------
rows = []

for sample in sample_order:

    fn = (
        mapping_root
        / sample
        / f"{sample}_priority24_support.tsv"
    )

    if not fn.exists():
        raise SystemExit(f"Missing: {fn}")

    with fn.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):

            cid = r["cluster_id"]

            if cid not in annotation:
                raise SystemExit(
                    f"Missing annotation for {cid}"
                )

            a = annotation[cid]

            rows.append({
                "sample": sample,
                "family": a["family"],
                "cluster_id": cid,
                "biological_label": a["biological_label"],
                "virus_group": a["virus_group"],
                "genome_unit": a["genome_unit"],
                "representative_length": r["length"],
                "mapq20_fragments": r["mapq20_fragments"],
                "FPM": r["FPM"],
                "breadth1_mapq20_pct":
                    r["breadth1_mapq20_pct"],
                "breadth10_mapq20_pct":
                    r["breadth10_mapq20_pct"],
                "mean_depth_mapq20":
                    r["mean_depth_mapq20"],
            })

# ----------------------------
# Long-format table
# ----------------------------
fields = [
    "sample",
    "family",
    "cluster_id",
    "biological_label",
    "virus_group",
    "genome_unit",
    "representative_length",
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

# ----------------------------
# Cluster order
# ----------------------------
cluster_order = []

for r in rows:
    cid = r["cluster_id"]

    if cid not in cluster_order:
        cluster_order.append(cid)

lookup = {
    (r["cluster_id"], r["sample"]): r
    for r in rows
}

# ----------------------------
# Matrix writer
# ----------------------------
def write_matrix(path, metric):

    with path.open("w") as f:

        header = [
            "cluster_id",
            "biological_label",
            *sample_order,
        ]

        f.write("\t".join(header) + "\n")

        for cid in cluster_order:

            label = annotation[cid]["biological_label"]

            vals = [
                lookup[(cid, s)][metric]
                for s in sample_order
            ]

            f.write(
                "\t".join(
                    [cid, label, *vals]
                )
                + "\n"
            )

write_matrix(
    out_breadth,
    "breadth1_mapq20_pct",
)

write_matrix(
    out_fpm,
    "FPM",
)

print(f"Samples: {len(sample_order)}")
print(f"Clusters: {len(cluster_order)}")
print(f"Rows: {len(rows)}")
print(f"Saved: {out_long}")
print(f"Saved: {out_breadth}")
print(f"Saved: {out_fpm}")
