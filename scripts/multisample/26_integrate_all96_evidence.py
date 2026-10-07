#!/usr/bin/env python3

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

clusters_file = Path(
    "results_summary/multisample/"
    "supported_sequence_clusters.tsv"
)

annotation_file = Path(
    "results_summary/multisample/"
    "all96_annotation_master.tsv"
)

mapping_root = Path(
    "results/multisample/all96_mapping"
)

out_long = Path(
    "results_summary/multisample/"
    "all96_evidence_master.tsv"
)

out_matrix = Path(
    "results_summary/multisample/"
    "all96_evidence_matrix.tsv"
)

out_summary = Path(
    "results_summary/multisample/"
    "all96_evidence_sample_summary.tsv"
)

samples = [
    "ERR14788990",
    "ERR14788988",
    "ERR14788991",
    "ERR14788992",
    "ERR14788995",
    "ERR14789044",
    "ERR14788864",
    "ERR14788954",
]

# --------------------------------------------------
# Cluster membership from assembly-supported contigs
# --------------------------------------------------
cluster_meta = {}

with clusters_file.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):

        r["assembly_samples"] = set(
            re.findall(r"ERR\d+", r["samples"])
        )

        cluster_meta[r["cluster_id"]] = r

# --------------------------------------------------
# Annotation
# --------------------------------------------------
annotation = {}

with annotation_file.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):
        annotation[r["cluster_id"]] = r

if set(cluster_meta) != set(annotation):
    raise SystemExit(
        "ERROR: cluster sets differ between "
        "cluster table and annotation table"
    )

# --------------------------------------------------
# Evidence classes
#
# Descriptive workflow tiers only.
# Not diagnostic positivity thresholds.
# --------------------------------------------------
def classify(assembly_member, breadth):

    if assembly_member:
        return "ASSEMBLY"

    if breadth >= 80:
        return "HIGH_BREADTH_MAPPING"

    if breadth >= 50:
        return "PARTIAL_MAPPING"

    if breadth >= 10:
        return "LOCAL_MAPPING"

    return "NONE"

rows = []

for sample in samples:

    fn = (
        mapping_root
        / sample
        / f"{sample}_all96_support.tsv"
    )

    if not fn.exists():
        raise SystemExit(f"Missing mapping file: {fn}")

    with fn.open() as f:

        for m in csv.DictReader(f, delimiter="\t"):

            cid = m["cluster_id"]

            c = cluster_meta[cid]
            a = annotation[cid]

            assembly_member = (
                sample in c["assembly_samples"]
            )

            breadth = float(
                m["breadth1_mapq20_pct"]
            )

            evidence = classify(
                assembly_member,
                breadth
            )

            rows.append({
                "sample": sample,
                "family": c["family"],
                "cluster_id": cid,

                "representative":
                    c["representative"],

                "representative_length":
                    c["representative_length"],

                "cluster_size":
                    c["cluster_size"],

                "cluster_sample_count":
                    c["sample_count"],

                "cluster_samples":
                    c["samples"],

                "assembly_member":
                    "YES" if assembly_member else "NO",

                "evidence_class":
                    evidence,

                "best_nt_accession":
                    a["best_nt_accession"],

                "best_nt_name":
                    a["best_nt_name"],

                "query_coverage_pct":
                    a["query_coverage_pct"],

                "top_hsp_identity_pct":
                    a["top_hsp_identity_pct"],

                "match_quality":
                    a["match_quality"],

                "annotation_source":
                    a["annotation_source"],

                "mapq20_fragments":
                    m["mapq20_fragments"],

                "FPM":
                    m["FPM"],

                "breadth1_mapq20_pct":
                    m["breadth1_mapq20_pct"],

                "breadth10_mapq20_pct":
                    m["breadth10_mapq20_pct"],

                "mean_depth_mapq20":
                    m["mean_depth_mapq20"],
            })

# --------------------------------------------------
# Long master table
# --------------------------------------------------
fields = [
    "sample",
    "family",
    "cluster_id",
    "representative",
    "representative_length",
    "cluster_size",
    "cluster_sample_count",
    "cluster_samples",
    "assembly_member",
    "evidence_class",
    "best_nt_accession",
    "best_nt_name",
    "query_coverage_pct",
    "top_hsp_identity_pct",
    "match_quality",
    "annotation_source",
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

# --------------------------------------------------
# Evidence matrix
# --------------------------------------------------
cluster_order = list(cluster_meta)

lookup = {
    (r["cluster_id"], r["sample"]):
        r["evidence_class"]
    for r in rows
}

with out_matrix.open("w") as f:

    f.write(
        "\t".join([
            "family",
            "cluster_id",
            "best_nt_name",
            *samples
        ]) + "\n"
    )

    for cid in cluster_order:

        a = annotation[cid]

        vals = [
            lookup[(cid, s)]
            for s in samples
        ]

        f.write(
            "\t".join([
                cluster_meta[cid]["family"],
                cid,
                a["best_nt_name"],
                *vals
            ]) + "\n"
        )

# --------------------------------------------------
# Sample summary
# --------------------------------------------------
counts = defaultdict(Counter)

for r in rows:
    counts[r["sample"]][
        r["evidence_class"]
    ] += 1

summary_fields = [
    "sample",
    "ASSEMBLY",
    "HIGH_BREADTH_MAPPING",
    "PARTIAL_MAPPING",
    "LOCAL_MAPPING",
    "NONE",
]

with out_summary.open("w") as f:

    w = csv.DictWriter(
        f,
        fieldnames=summary_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()

    for sample in samples:

        w.writerow({
            "sample": sample,
            "ASSEMBLY":
                counts[sample]["ASSEMBLY"],

            "HIGH_BREADTH_MAPPING":
                counts[sample]["HIGH_BREADTH_MAPPING"],

            "PARTIAL_MAPPING":
                counts[sample]["PARTIAL_MAPPING"],

            "LOCAL_MAPPING":
                counts[sample]["LOCAL_MAPPING"],

            "NONE":
                counts[sample]["NONE"],
        })

print(f"Samples: {len(samples)}")
print(f"Clusters: {len(cluster_meta)}")
print(f"Rows: {len(rows)}")

print(f"Saved: {out_long}")
print(f"Saved: {out_matrix}")
print(f"Saved: {out_summary}")
