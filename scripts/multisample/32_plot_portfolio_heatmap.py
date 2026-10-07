#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

inp = Path(
    "results_summary/multisample/"
    "portfolio_enteric_group_matrix.tsv"
)

outdir = Path("figures/multisample")
outdir.mkdir(parents=True, exist_ok=True)

png = outdir / "portfolio_enteric_evidence_heatmap.png"
pdf = outdir / "portfolio_enteric_evidence_heatmap.pdf"

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

with inp.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

score = {
    "NONE": 0,
    "LOCAL_MAPPING": 1,
    "PARTIAL_MAPPING": 2,
    "HIGH_BREADTH_MAPPING": 3,
    "ASSEMBLY": 4,
}

matrix = [
    [score[r[s]] for s in samples]
    for r in rows
]

labels = [
    r["display_label"]
    for r in rows
]

families = [
    r["family"]
    for r in rows
]

cmap = ListedColormap([
    "#f2f2f2",
    "#d9d9d9",
    "#9ecae1",
    "#3182bd",
    "#08519c",
])

# Wider figure
fig, ax = plt.subplots(
    figsize=(16, 14)
)

ax.imshow(
    matrix,
    aspect="auto",
    interpolation="nearest",
    cmap=cmap,
    vmin=-0.5,
    vmax=4.5,
)

# X axis
ax.set_xticks(range(len(samples)))
ax.set_xticklabels(
    samples,
    rotation=45,
    ha="right",
    fontsize=10,
)

# Y axis
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(
    labels,
    fontsize=9,
)

ax.set_xlabel(
    "Wastewater hybrid-capture library",
    fontsize=12,
    labelpad=10,
)

ax.set_ylabel(
    "Virus group / genome segment",
    fontsize=12,
    labelpad=12,
)

ax.set_title(
    "Integrated viral evidence across eight "
    "hybrid-capture wastewater libraries",
    fontsize=16,
    pad=18,
)

# --------------------------------------------------
# Family separators and labels
# --------------------------------------------------
family_ranges = []
start = 0

for i in range(1, len(families) + 1):

    if (
        i == len(families)
        or families[i] != families[start]
    ):

        family_ranges.append(
            (families[start], start, i - 1)
        )

        if i < len(families):
            ax.axhline(
                i - 0.5,
                linewidth=0.9,
                color="black",
                alpha=0.35,
            )

        start = i

# Put family labels on RIGHT side
for family, first, last in family_ranges:

    midpoint = (first + last) / 2

    ax.text(
        1.02,
        midpoint,
        family,
        transform=ax.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=10,
        fontweight="bold",
    )

# --------------------------------------------------
# Legend below heatmap
# --------------------------------------------------
legend = [
    Patch(
        facecolor="#08519c",
        label="Assembly-supported"
    ),
    Patch(
        facecolor="#3182bd",
        label="High-breadth mapping (≥80%)"
    ),
    Patch(
        facecolor="#9ecae1",
        label="Partial mapping (50–<80%)"
    ),
    Patch(
        facecolor="#d9d9d9",
        label="Local mapping (10–<50%)"
    ),
    Patch(
        facecolor="#f2f2f2",
        label="<10% breadth / none"
    ),
]

ax.legend(
    handles=legend,
    bbox_to_anchor=(0.5, -0.14),
    loc="upper center",
    ncol=3,
    frameon=False,
    fontsize=10,
)

# Footnote
fig.text(
    0.5,
    0.02,
    (
        "Evidence tiers are descriptive workflow categories, "
        "not diagnostic positivity thresholds. "
        "Capture-derived read support is not interpreted as abundance."
    ),
    ha="center",
    fontsize=9,
)

# Allocate more width to heatmap body
fig.subplots_adjust(
    left=0.34,
    right=0.78,
    bottom=0.22,
    top=0.93,
)

fig.savefig(
    png,
    dpi=300,
    bbox_inches="tight",
)

fig.savefig(
    pdf,
    bbox_inches="tight",
)

print("Groups:", len(rows))
print("Saved:", png)
print("Saved:", pdf)
