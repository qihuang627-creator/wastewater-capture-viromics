#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

clusters_file = Path(
    "results_summary/multisample/"
    "supported_sequence_clusters.tsv"
)

blast_file = Path(
    "results/multisample/blast/"
    "local_nt_viruses_all96/"
    "all96.megablast.tsv"
)

outfile = Path(
    "results_summary/multisample/"
    "all96_ntviruses_annotation.tsv"
)

fields = [
    "qseqid",
    "sacc",
    "pident",
    "length",
    "qlen",
    "qstart",
    "qend",
    "qcovs",
    "evalue",
    "bitscore",
    "staxids",
    "sscinames",
    "stitle",
]

# -----------------------------------
# Cluster metadata
# -----------------------------------
with clusters_file.open() as f:
    cluster_rows = list(
        csv.DictReader(f, delimiter="\t")
    )

cluster_meta = {
    r["cluster_id"]: r
    for r in cluster_rows
}

# -----------------------------------
# Keep BLAST's first-ranked subject
# for each query, then aggregate HSPs
# belonging to that subject.
# -----------------------------------
top_accession = {}
hits = defaultdict(list)

with blast_file.open() as f:

    reader = csv.DictReader(
        f,
        fieldnames=fields,
        delimiter="\t",
    )

    for r in reader:

        q = r["qseqid"]
        acc = r["sacc"]

        if q not in top_accession:
            top_accession[q] = acc

        if acc == top_accession[q]:
            hits[q].append(r)


def triage(qcov, ident):

    # Workflow triage only:
    # NOT species/genotype thresholds.

    if qcov >= 80 and ident >= 90:
        return "strong"

    if qcov >= 50 and ident >= 80:
        return "moderate"

    return "partial_or_divergent"


out_fields = [
    "family",
    "cluster_id",
    "representative",
    "representative_length",
    "cluster_size",
    "sample_count",
    "samples",
    "best_nt_accession",
    "best_nt_name",
    "query_coverage_pct",
    "top_hsp_identity_pct",
    "hsp_count",
    "total_bitscore",
    "match_quality",
    "best_nt_title",
]

out_rows = []

for c in cluster_rows:

    cid = c["cluster_id"]

    if cid not in hits:

        out_rows.append({
            "family": c["family"],
            "cluster_id": cid,
            "representative": c["representative"],
            "representative_length":
                c["representative_length"],
            "cluster_size": c["cluster_size"],
            "sample_count": c["sample_count"],
            "samples": c["samples"],
            "best_nt_accession": "",
            "best_nt_name": "",
            "query_coverage_pct": "",
            "top_hsp_identity_pct": "",
            "hsp_count": "0",
            "total_bitscore": "",
            "match_quality": "no_megablast_hit",
            "best_nt_title": "",
        })

        continue

    hsps = hits[cid]

    # qcovs is subject-level query coverage
    qcov = max(
        float(x["qcovs"])
        for x in hsps
    )

    # Highest-scoring HSP identity,
    # not whole-query identity.
    best_hsp = max(
        hsps,
        key=lambda x: float(x["bitscore"])
    )

    ident = float(best_hsp["pident"])

    total_bitscore = sum(
        float(x["bitscore"])
        for x in hsps
    )

    first = hsps[0]

    out_rows.append({
        "family": c["family"],
        "cluster_id": cid,
        "representative": c["representative"],
        "representative_length":
            c["representative_length"],
        "cluster_size": c["cluster_size"],
        "sample_count": c["sample_count"],
        "samples": c["samples"],
        "best_nt_accession": first["sacc"],
        "best_nt_name": first["sscinames"],
        "query_coverage_pct":
            f"{qcov:.1f}",
        "top_hsp_identity_pct":
            f"{ident:.3f}",
        "hsp_count": len(hsps),
        "total_bitscore":
            f"{total_bitscore:.1f}",
        "match_quality":
            triage(qcov, ident),
        "best_nt_title": first["stitle"],
    })

with outfile.open("w") as f:

    w = csv.DictWriter(
        f,
        fieldnames=out_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(out_rows)

print(f"Clusters: {len(out_rows)}")

counts = defaultdict(int)

for r in out_rows:
    counts[r["match_quality"]] += 1

for status in [
    "strong",
    "moderate",
    "partial_or_divergent",
    "no_megablast_hit",
]:
    print(
        f"{status}: "
        f"{counts[status]}"
    )

print(f"Saved: {outfile}")
