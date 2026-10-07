#!/usr/bin/env python3

import csv
from pathlib import Path

clusters_file = Path(
    "results_summary/multisample/"
    "supported_sequence_clusters.tsv"
)

hits_file = Path(
    "results_summary/multisample/"
    "shortlist24_ntviruses_bestsubjects.tsv"
)

shortlist_file = Path(
    "results_summary/multisample/"
    "blast_shortlist.tsv"
)

outfile = Path(
    "results_summary/multisample/"
    "priority_cluster_annotation.tsv"
)

# Biological labels are deliberately conservative.
labels = {

    "Adenoviridae_C001":
        ("Human adenovirus 41-like", "HAdV41", "genome"),

    "Adenoviridae_C002":
        ("Human adenovirus 41-like", "HAdV41", "genome"),

    "Astroviridae_C001":
        ("Human astrovirus 1-like", "HAstV1", "genome"),

    "Astroviridae_C002":
        ("Human astrovirus 1-like", "HAstV1", "genome"),

    "Astroviridae_C005":
        ("Human astrovirus 3-like", "HAstV3", "genome"),

    "Astroviridae_C006":
        ("Human astrovirus 5-like", "HAstV5", "genome"),

    "Astroviridae_C008":
        ("Human astrovirus 3-like", "HAstV3", "genome"),

    "Caliciviridae_C004":
        ("Sapovirus GI.2-like", "SaV_GI2", "genome"),

    "Caliciviridae_C009":
        ("Sapovirus GI.2-like", "SaV_GI2", "genome"),

    "Parvoviridae_C001":
        ("Human bocavirus 3-like", "HBoV3", "genome"),

    "Parvoviridae_C002":
        ("Rodent bocavirus-like", "Rodent_BoV", "genome"),

    "Parvoviridae_C004":
        ("Human bocavirus 2-like", "HBoV2", "genome"),

    "Parvoviridae_C006":
        ("Human bocavirus 2-like", "HBoV2", "genome"),

    "Picornaviridae_C001":
        (
            "Divergent Enterovirus B-like (B106 best hit)",
            "EVB_like",
            "genome"
        ),

    "Picornaviridae_C002":
        ("Enterovirus B84-like", "EVB84", "genome"),

    "Picornaviridae_C003":
        ("Aichivirus A1-like", "Aichi_A1", "genome"),

    "Picornaviridae_C004":
        (
            "Kobuvirus-like (Aichivirus-related)",
            "Kobuvirus_like",
            "genome"
        ),

    "Picornaviridae_C005":
        ("Coxsackievirus A19-like", "CVA19", "genome"),

    "Sedoreoviridae_C001":
        ("Rotavirus A segment 1 / VP1-like", "RVA_seg1", "segment1"),

    "Sedoreoviridae_C004":
        ("Rotavirus A segment 3 / VP3-like", "RVA_seg3", "segment3"),

    "Sedoreoviridae_C011":
        ("Rotavirus A segment 1 / VP1-like", "RVA_seg1", "segment1"),

    "Sedoreoviridae_C012":
        ("Rotavirus A segment 1 / VP1-like", "RVA_seg1", "segment1"),

    "Sedoreoviridae_C013":
        ("Rotavirus A segment 3 / VP3-like", "RVA_seg3", "segment3"),

    "Sedoreoviridae_C014":
        ("Rotavirus A segment 2 / VP2-like", "RVA_seg2", "segment2"),
}


def load(fn, key):
    with fn.open() as f:
        return {
            r[key]: r
            for r in csv.DictReader(f, delimiter="\t")
        }


clusters = load(clusters_file, "cluster_id")
hits = load(hits_file, "cluster_id")
shortlist = load(shortlist_file, "cluster_id")

if set(shortlist) != set(labels):
    missing_labels = sorted(set(shortlist) - set(labels))
    extra_labels = sorted(set(labels) - set(shortlist))

    raise SystemExit(
        f"Annotation mismatch\n"
        f"Missing labels: {missing_labels}\n"
        f"Extra labels: {extra_labels}"
    )

rows = []

for cid in shortlist:

    c = clusters[cid]
    h = hits[cid]
    biological_label, virus_group, genome_unit = labels[cid]

    rows.append({
        "priority": shortlist[cid]["priority"],
        "family": c["family"],
        "cluster_id": cid,
        "biological_label": biological_label,
        "virus_group": virus_group,
        "genome_unit": genome_unit,
        "representative": c["representative"],
        "representative_length": c["representative_length"],
        "cluster_size": c["cluster_size"],
        "sample_count": c["sample_count"],
        "samples": c["samples"],
        "best_nt_accession": h["accession"],
        "best_nt_name": h["scientific_name"],
        "query_coverage_pct": h["query_coverage_pct"],
        "max_local_identity_pct": h["max_pident_pct"],
        "hsp_count": h["hsp_count"],
        "best_nt_title": h["title"],
    })

rows.sort(
    key=lambda r: (
        -int(r["sample_count"]),
        r["family"],
        r["cluster_id"],
    )
)

fields = [
    "priority",
    "family",
    "cluster_id",
    "biological_label",
    "virus_group",
    "genome_unit",
    "representative",
    "representative_length",
    "cluster_size",
    "sample_count",
    "samples",
    "best_nt_accession",
    "best_nt_name",
    "query_coverage_pct",
    "max_local_identity_pct",
    "hsp_count",
    "best_nt_title",
]

with outfile.open("w") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(rows)

print(f"Saved: {outfile}")
print(f"Annotated clusters: {len(rows)}")
