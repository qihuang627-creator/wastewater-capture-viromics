#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

inp = Path(
    "results_summary/multisample/"
    "biological_group_evidence_matrix.tsv"
)

outdir = Path("figures/multisample")
outdir.mkdir(parents=True, exist_ok=True)

png = outdir / "full_biological_group_evidence_heatmap.png"
pdf = outdir / "full_biological_group_evidence_heatmap.pdf"

rows = list(
    csv.DictReader(
        inp.open(),
        delimiter="\t"
    )
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

score = {
    "NONE": 0,
    "LOCAL_MAPPING": 1,
    "PARTIAL_MAPPING": 2,
    "HIGH_BREADTH_MAPPING": 3,
    "ASSEMBLY": 4,
}

matrix = [
    [
        score[r[s]]
        for s in samples
    ]
    for r in rows
]

labels = [
    f'{r["family"]} | {r["display_label"]}'
    for r in rows
]

fig, ax = plt.subplots(
    figsize=(13, 19)
)

cmap = ListedColormap([
    "#F2F2F2",  # none
    "#C7D5DC",  # local mapping
    "#9CCBC4",  # partial mapping
    "#4C9F9B",  # high-breadth mapping
    "#2F5D7C",  # assembly-supported
])

im = ax.imshow(
    matrix,
    aspect="auto",
    interpolation="nearest",
    cmap=cmap,
    vmin=-0.5,
    vmax=4.5,
)

ax.set_xticks(range(len(samples)))
ax.set_xticklabels(
    samples,
    rotation=45,
    ha="right",
)

ax.set_yticks(range(len(labels)))
ax.set_yticklabels(
    labels,
    fontsize=8,
)

ax.set_xlabel("Wastewater capture library")
ax.set_ylabel("Biological virus group / genome segment")

ax.set_title(
    "Integrated assembly and competitive-mapping evidence "
    "across eight capture libraries"
)

legend = [
    Patch(
        facecolor="#2F5D7C",
        label="Assembly-supported"
    ),
    Patch(
        facecolor="#4C9F9B",
        label="High-breadth mapping (≥80%)"
    ),
    Patch(
        facecolor="#9CCBC4",
        label="Partial mapping (50–<80%)"
    ),
    Patch(
        facecolor="#C7D5DC",
        label="Local mapping (10–<50%)"
    ),
    Patch(
        facecolor="#F2F2F2",
        edgecolor="#D0D0D0",
        label="<10% breadth / none"
    ),
]

ax.legend(
    handles=legend,
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False,
    fontsize=9,
)

plt.tight_layout()

fig.savefig(
    png,
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    pdf,
    bbox_inches="tight"
)

print("Saved:", png)
print("Saved:", pdf)
