#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

mega_file = Path(
    "results_summary/multisample/"
    "all96_ntviruses_annotation.tsv"
)

sens_file = Path(
    "results/multisample/blast/"
    "local_nt_viruses_all96/"
    "nonstrong10.blastn.tsv"
)

outfile = Path(
    "results_summary/multisample/"
    "all96_annotation_master.tsv"
)

blast_fields = [
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

# --------------------------------
# Load existing megablast table
# --------------------------------
with mega_file.open() as f:
    mega_rows = list(
        csv.DictReader(f, delimiter="\t")
    )

# --------------------------------
# Sensitive blastn:
# use BLAST first-ranked subject
# and summarize all its HSPs
# --------------------------------
top_acc = {}
sens_hits = defaultdict(list)

with sens_file.open() as f:

    reader = csv.DictReader(
        f,
        fieldnames=blast_fields,
        delimiter="\t"
    )

    for r in reader:

        q = r["qseqid"]

        if q not in top_acc:
            top_acc[q] = r["sacc"]

        if r["sacc"] == top_acc[q]:
            sens_hits[q].append(r)


def triage(qcov, ident):

    # Workflow triage only.
    # Not formal taxonomic thresholds.

    if qcov >= 80 and ident >= 90:
        return "strong"

    if qcov >= 50 and ident >= 80:
        return "moderate"

    return "partial_or_divergent"


master = []

for r in mega_rows:

    cid = r["cluster_id"]

    row = dict(r)

    if cid in sens_hits:

        hsps = sens_hits[cid]

        first = hsps[0]

        qcov = max(
            float(x["qcovs"])
            for x in hsps
        )

        best_hsp = max(
            hsps,
            key=lambda x: float(x["bitscore"])
        )

        ident = float(best_hsp["pident"])

        total_bitscore = sum(
            float(x["bitscore"])
            for x in hsps
        )

        row["best_nt_accession"] = first["sacc"]
        row["best_nt_name"] = first["sscinames"]

        row["query_coverage_pct"] = (
            f"{qcov:.1f}"
        )

        row["top_hsp_identity_pct"] = (
            f"{ident:.3f}"
        )

        row["hsp_count"] = str(len(hsps))

        row["total_bitscore"] = (
            f"{total_bitscore:.1f}"
        )

        row["match_quality"] = triage(
            qcov,
            ident
        )

        row["best_nt_title"] = first["stitle"]
        row["annotation_source"] = "sensitive_blastn"

    else:
        row["annotation_source"] = "megablast"

    master.append(row)


fields = list(mega_rows[0].keys()) + [
    "annotation_source"
]

with outfile.open("w") as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(master)


counts = defaultdict(int)

for r in master:
    counts[r["match_quality"]] += 1


print(f"Clusters: {len(master)}")

for x in [
    "strong",
    "moderate",
    "partial_or_divergent",
    "no_megablast_hit",
]:
    print(f"{x}: {counts[x]}")

print(
    "Sensitive replacements:",
    sum(
        r["annotation_source"] == "sensitive_blastn"
        for r in master
    )
)

print("Saved:", outfile)
