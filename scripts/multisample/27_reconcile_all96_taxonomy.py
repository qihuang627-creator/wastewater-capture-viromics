#!/usr/bin/env python3

import csv
from collections import Counter
from pathlib import Path

inp = Path(
    "results_summary/multisample/"
    "all96_annotation_master.tsv"
)

out = Path(
    "results_summary/multisample/"
    "all96_annotation_reconciled.tsv"
)

# Strong nucleotide-level conflicts with geNomad family calls.
# Keep original family; add reconciled family separately.
corrections = {

    "Picornaviridae_C019": {
        "family": "Caliciviridae",
        "group": "Sapovirus GV.1",
        "reason":
            "nt_viruses: Sapovirus GV.1; "
            "99% query coverage; 97.047% top-HSP identity",
    },

    "Picornaviridae_C020": {
        "family": "Caliciviridae",
        "group": "Norovirus GII",
        "reason":
            "nt_viruses: Norovirus GII; "
            "94% query coverage; 96.832% top-HSP identity",
    },

    "Picornaviridae_C021": {
        "family": "Caliciviridae",
        "group": "Norovirus GII",
        "reason":
            "nt_viruses: Norovirus GII; "
            "100% query coverage; 98.805% top-HSP identity",
    },
}

with inp.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

output = []

for r in rows:

    cid = r["cluster_id"]

    original_family = r["family"]

    if cid in corrections:

        c = corrections[cid]

        reconciled_family = c["family"]
        conflict = "YES"
        reason = c["reason"]
        group_candidate = c["group"]

    else:

        reconciled_family = original_family
        conflict = "NO"
        reason = ""
        group_candidate = r["best_nt_name"]

    x = dict(r)

    # Do not silently overwrite the original taxonomy.
    x["original_family"] = original_family
    x["reconciled_family"] = reconciled_family
    x["taxonomy_conflict"] = conflict
    x["reconciliation_reason"] = reason
    x["biological_group_candidate"] = group_candidate

    output.append(x)

fields = list(rows[0].keys()) + [
    "original_family",
    "reconciled_family",
    "taxonomy_conflict",
    "reconciliation_reason",
    "biological_group_candidate",
]

with out.open("w") as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(output)


print("Clusters:", len(output))
print(
    "Taxonomy conflicts corrected:",
    sum(x["taxonomy_conflict"] == "YES" for x in output)
)

print("\nOriginal family counts:")
for k, v in sorted(
    Counter(x["original_family"] for x in output).items()
):
    print(f"{k}\t{v}")

print("\nReconciled family counts:")
for k, v in sorted(
    Counter(x["reconciled_family"] for x in output).items()
):
    print(f"{k}\t{v}")

print("\nSaved:", out)
