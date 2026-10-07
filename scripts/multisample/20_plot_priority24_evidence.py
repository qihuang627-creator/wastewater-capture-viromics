#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

inp = Path(
    "results_summary/multisample/"
    "priority24_evidence_integrated.tsv"
)

outdir = Path(
    "figures/multisample"
)
outdir.mkdir(parents=True, exist_ok=True)

png = outdir / "priority24_evidence_heatmap.png"
pdf = outdir / "priority24_evidence_heatmap.pdf"

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

rows = list(
    csv.DictReader(
        inp.open(),
        delimiter="\t"
    )
)

cluster_order = []
labels = {}

for r in rows:
    cid = r["cluster_id"]

    if cid not in cluster_order:
        cluster_order.append(cid)
        labels[cid] = (
            f'{r["biological_label"]}  [{cid}]'
        )

lookup = {
    (r["cluster_id"], r["sample"]):
        r["evidence_class"]
    for r in rows
}

classes = {
    "NONE": 0,
    "LOCAL": 1,
    "PARTIAL": 2,
    "MAPPING": 3,
    "ASSEMBLY": 4,
}

matrix = [
    [
        classes[
            lookup[(cid, sample)]
        ]
        for sample in samples
    ]
    for cid in cluster_order
]

# Deliberately simple categorical palette.
# Colors are assigned only for evidence-class readability.
cmap = ListedColormap([
    "#f2f2f2",
    "#d9d9d9",
    "#9ecae1",
    "#3182bd",
    "#08519c",
])

fig, ax = plt.subplots(
    figsize=(12, 13)
)

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

ax.set_yticks(range(len(cluster_order)))
ax.set_yticklabels([
    labels[cid]
    for cid in cluster_order
])

ax.set_xlabel("Wastewater capture library")
ax.set_ylabel("Priority viral sequence cluster")

ax.set_title(
    "Assembly and competitive-mapping evidence "
    "across eight capture libraries"
)

legend = [
    Patch(
        facecolor="#08519c",
        label="Assembly-supported"
    ),
    Patch(
        facecolor="#3182bd",
        label="Mapping breadth ≥80%"
    ),
    Patch(
        facecolor="#9ecae1",
        label="Mapping breadth 50–<80%"
    ),
    Patch(
        facecolor="#d9d9d9",
        label="Local mapping 10–<50%"
    ),
    Patch(
        facecolor="#f2f2f2",
        label="<10% breadth / none"
    ),
]

ax.legend(
    handles=legend,
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False,
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
