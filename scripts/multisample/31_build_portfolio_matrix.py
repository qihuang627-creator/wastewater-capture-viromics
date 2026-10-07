#!/usr/bin/env python3

import csv
from pathlib import Path

src = Path(
    "results_summary/multisample/"
    "biological_group_evidence_matrix.tsv"
)

dst = Path(
    "results_summary/multisample/"
    "portfolio_enteric_group_matrix.tsv"
)

# Keep full results elsewhere, but remove reference-match /
# unresolved groups that distract from the main portfolio figure.
exclude = {
    "Rodent_astrovirus_like",
    "Canine_BoV_like",
    "Porcine_BoV_like",
    "Rodent_BoV_like",
    "Rat_parvovirus_like",
    "Picornavirus_unresolved",
    "Enterovirus_unresolved",
    "Rotavirus_alpha_like_segment_1",
    "Rotavirus_alpha_like_segment_2",
}

with src.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

fields = list(rows[0])

kept = []

for r in rows:

    if r["biological_group"] in exclude:
        continue

    # Avoid implying a third HAstV-1 population.
    if r["biological_group"] == "HAstV1":
        r["display_label"] = (
            "HAstV-1-like | population assignment unresolved"
        )

    elif r["biological_group"] == "HAstV1_A":
        r["display_label"] = (
            "HAstV-1-like | population A"
        )

    elif r["biological_group"] == "HAstV1_B":
        r["display_label"] = (
            "HAstV-1-like | population B"
        )

    kept.append(r)

with dst.open("w") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(kept)

print("Full groups:", len(rows))
print("Portfolio groups:", len(kept))
print("Excluded:", len(rows) - len(kept))
print("Saved:", dst)
