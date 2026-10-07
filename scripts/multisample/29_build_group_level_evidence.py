#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict
from pathlib import Path

group_file = Path(
    "results_summary/multisample/"
    "all96_biological_grouping.tsv"
)

evidence_file = Path(
    "results_summary/multisample/"
    "all96_evidence_master.tsv"
)

out_long = Path(
    "results_summary/multisample/"
    "biological_group_evidence.tsv"
)

out_matrix = Path(
    "results_summary/multisample/"
    "biological_group_evidence_matrix.tsv"
)

out_summary = Path(
    "results_summary/multisample/"
    "biological_group_sample_summary.tsv"
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

family_order = {
    "Adenoviridae": 0,
    "Astroviridae": 1,
    "Caliciviridae": 2,
    "Parvoviridae": 3,
    "Picornaviridae": 4,
    "Sedoreoviridae": 5,
}

rank = {
    "NONE": 0,
    "LOCAL_MAPPING": 1,
    "PARTIAL_MAPPING": 2,
    "HIGH_BREADTH_MAPPING": 3,
    "ASSEMBLY": 4,
}

# --------------------------------------------------
# Human-readable group labels
# --------------------------------------------------
pretty = {
    "HAdV41": "Human adenovirus 41-like",
    "HAdV_D_like": "Human mastadenovirus D-like",

    "HAstV1_A": "HAstV-1-like population A",
    "HAstV1_B": "HAstV-1-like population B",
    "HAstV1": "HAstV-1-like (population unresolved)",
    "HAstV2": "HAstV-2-like",
    "HAstV3": "HAstV-3-like",
    "HAstV4": "HAstV-4-like",
    "HAstV5": "HAstV-5-like",
    "HAstV8": "HAstV-8-like",
    "HAstV_unresolved": "Human astrovirus-like, unresolved type",
    "MLB1": "Astrovirus MLB1-like",
    "MLB2": "Astrovirus MLB2-like",
    "VA5": "Astrovirus VA5-like",
    "Mamastrovirus1_like": "Mamastrovirus 1-like",
    "Mamastrovirus2_like": "Mamastrovirus 2-like",
    "Rodent_astrovirus_like": "Rodent astrovirus-like (reference match)",

    "NoV_GI": "Norovirus GI-like",
    "NoV_GII": "Norovirus GII-like",
    "SaV_GI": "Sapovirus GI-like",
    "SaV_GII": "Sapovirus GII-like",
    "SaV_GV": "Sapovirus GV-like",

    "HBoV2": "Human bocavirus 2-like",
    "HBoV3": "Human bocavirus 3-like",
    "Canine_BoV_like": "Canine bocavirus-like (reference match)",
    "Porcine_BoV_like": "Porcine bocavirus-like (reference match)",
    "Rodent_BoV_like": "Rodent bocavirus-like (reference match)",
    "Rat_parvovirus_like": "Rat parvovirus-like (reference match)",

    "Aichivirus_A1_like": "Aichivirus A1-like",
    "Aichivirus_A_like": "Divergent Aichivirus A-like",
    "Enterovirus_A_like": "Divergent Enterovirus A-like",
    "Enterovirus_B_like": "Divergent Enterovirus B-like",
    "Enterovirus_unresolved": "Enterovirus-like, unresolved",
    "Kobuvirus_like": "Kobuvirus-like",
    "Picornavirus_unresolved": "Picornavirus-like, unresolved",
    "Salivirus_A_like": "Salivirus A-like",
    "Coxsackievirus_A13_like": "Coxsackievirus A13-like",
    "Coxsackievirus_A19_like": "Coxsackievirus A19-like",
    "Echovirus_E18_like": "Echovirus 18-like",
    "Echovirus_E9_like": "Echovirus 9-like",
}

def display_label(group):
    if group in pretty:
        return pretty[group]

    if group.startswith("RVA_segment_"):
        n = group.split("_")[-1]
        return f"Rotavirus A | segment {n}"

    if group.startswith("RVC_segment_"):
        n = group.split("_")[-1]
        return f"Rotavirus C | segment {n}"

    if group.startswith("Rotavirus_alpha_like_segment_"):
        n = group.split("_")[-1]
        return (
            "Rotavirus alphagastroenteritidis-like "
            f"| segment {n}"
        )

    return group.replace("_", " ")

# --------------------------------------------------
# Load grouping
# --------------------------------------------------
grouping = {}

with group_file.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):
        grouping[r["cluster_id"]] = r

# --------------------------------------------------
# Merge cluster evidence with grouping
# --------------------------------------------------
merged = []

with evidence_file.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):

        cid = r["cluster_id"]

        if cid not in grouping:
            raise SystemExit(
                f"Missing biological grouping: {cid}"
            )

        g = grouping[cid]

        x = dict(r)

        x["reconciled_family"] = g["reconciled_family"]
        x["biological_group"] = g["biological_group"]
        x["biological_label"] = g["biological_label"]

        merged.append(x)

# --------------------------------------------------
# Aggregate by sample × biological group
# --------------------------------------------------
bucket = defaultdict(list)

for r in merged:
    key = (
        r["sample"],
        r["reconciled_family"],
        r["biological_group"],
    )
    bucket[key].append(r)

out_rows = []

for (sample, family, group), rs in bucket.items():

    # Best cluster = highest evidence class first,
    # then strongest coverage/depth.
    best = max(
        rs,
        key=lambda r: (
            rank[r["evidence_class"]],
            float(r["breadth1_mapq20_pct"]),
            float(r["breadth10_mapq20_pct"]),
            float(r["mean_depth_mapq20"]),
            float(r["FPM"]),
        )
    )

    counts = Counter(
        r["evidence_class"]
        for r in rs
    )

    max_rank = max(
        rank[r["evidence_class"]]
        for r in rs
    )

    group_class = next(
        k for k, v in rank.items()
        if v == max_rank
    )

    assembly_ids = [
        r["cluster_id"]
        for r in rs
        if r["evidence_class"] == "ASSEMBLY"
    ]

    high_ids = [
        r["cluster_id"]
        for r in rs
        if r["evidence_class"] == "HIGH_BREADTH_MAPPING"
    ]

    out_rows.append({
        "sample": sample,
        "family": family,
        "biological_group": group,
        "display_label": display_label(group),
        "group_evidence_class": group_class,

        "clusters_in_group": len(rs),

        "assembly_cluster_count":
            counts["ASSEMBLY"],

        "high_breadth_mapping_cluster_count":
            counts["HIGH_BREADTH_MAPPING"],

        "partial_mapping_cluster_count":
            counts["PARTIAL_MAPPING"],

        "local_mapping_cluster_count":
            counts["LOCAL_MAPPING"],

        "assembly_cluster_ids":
            ",".join(assembly_ids),

        "high_breadth_mapping_cluster_ids":
            ",".join(high_ids),

        "best_supported_cluster":
            best["cluster_id"],

        "max_breadth1_mapq20_pct":
            f'{max(float(r["breadth1_mapq20_pct"]) for r in rs):.2f}',

        "max_breadth10_mapq20_pct":
            f'{max(float(r["breadth10_mapq20_pct"]) for r in rs):.2f}',

        "max_mean_depth_mapq20":
            f'{max(float(r["mean_depth_mapq20"]) for r in rs):.2f}',

        # Peak support only; DO NOT sum FPM across clusters.
        "max_FPM":
            f'{max(float(r["FPM"]) for r in rs):.2f}',
    })

# --------------------------------------------------
# Stable group order
# --------------------------------------------------
groups = sorted(
    {
        (r["family"], r["biological_group"])
        for r in out_rows
    },
    key=lambda x: (
        family_order.get(x[0], 999),
        display_label(x[1])
    )
)

sample_index = {
    s: i for i, s in enumerate(samples)
}

out_rows.sort(
    key=lambda r: (
        family_order.get(r["family"], 999),
        display_label(r["biological_group"]),
        sample_index[r["sample"]],
    )
)

# --------------------------------------------------
# Long-format table
# --------------------------------------------------
fields = list(out_rows[0])

with out_long.open("w") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(out_rows)

# --------------------------------------------------
# Evidence matrix
# --------------------------------------------------
lookup = {
    (r["family"], r["biological_group"], r["sample"]):
        r["group_evidence_class"]
    for r in out_rows
}

with out_matrix.open("w") as f:

    f.write(
        "\t".join([
            "family",
            "biological_group",
            "display_label",
            *samples
        ]) + "\n"
    )

    for family, group in groups:

        vals = [
            lookup[(family, group, sample)]
            for sample in samples
        ]

        f.write(
            "\t".join([
                family,
                group,
                display_label(group),
                *vals
            ]) + "\n"
        )

# --------------------------------------------------
# Per-sample group summary
# --------------------------------------------------
summary = defaultdict(Counter)

for r in out_rows:
    summary[r["sample"]][
        r["group_evidence_class"]
    ] += 1

with out_summary.open("w") as f:

    cols = [
        "sample",
        "ASSEMBLY",
        "HIGH_BREADTH_MAPPING",
        "PARTIAL_MAPPING",
        "LOCAL_MAPPING",
        "NONE",
    ]

    w = csv.DictWriter(
        f,
        fieldnames=cols,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()

    for sample in samples:
        w.writerow({
            "sample": sample,
            "ASSEMBLY": summary[sample]["ASSEMBLY"],
            "HIGH_BREADTH_MAPPING":
                summary[sample]["HIGH_BREADTH_MAPPING"],
            "PARTIAL_MAPPING":
                summary[sample]["PARTIAL_MAPPING"],
            "LOCAL_MAPPING":
                summary[sample]["LOCAL_MAPPING"],
            "NONE":
                summary[sample]["NONE"],
        })

print("Samples:", len(samples))
print("Biological groups:", len(groups))
print("Rows:", len(out_rows))

print("Saved:", out_long)
print("Saved:", out_matrix)
print("Saved:", out_summary)
