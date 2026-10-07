#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

SAMPLE_SHEET = Path("config/multisample/samples.tsv")

OUTDIR = Path("results_summary/multisample")
OUTDIR.mkdir(parents=True, exist_ok=True)

SUPPORTED_OUT = OUTDIR / "supported_target_contigs.tsv"
SUMMARY_OUT = OUTDIR / "target_readback_summary_8samples.tsv"
ALL_FAMILY_OUT = OUTDIR / "target_family_contig_matrix.tsv"
SUPPORTED_FAMILY_OUT = OUTDIR / "supported_target_family_matrix.tsv"

TARGET_FAMILIES = [
    "Adenoviridae",
    "Astroviridae",
    "Caliciviridae",
    "Hepeviridae",
    "Parvoviridae",
    "Picornaviridae",
    "Sedoreoviridae",
    "Spinareoviridae",
]


def value(row, *names):
    for name in names:
        if name in row and row[name] != "":
            return row[name]
    return ""


with SAMPLE_SHEET.open() as f:
    samples = [
        r["sample"]
        for r in csv.DictReader(f, delimiter="\t")
    ]

supported_rows = []
summary_rows = []

all_family_counts = defaultdict(lambda: defaultdict(int))
supported_family_counts = defaultdict(lambda: defaultdict(int))

for sample in samples:

    fn = Path(
        f"results/multisample/target_readback/"
        f"{sample}/{sample}_target_support.tsv"
    )

    if not fn.exists():
        raise SystemExit(f"Missing support table: {fn}")

    with fn.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    passed = []

    for r in rows:

        family = value(r, "family")
        contig = value(r, "contig")

        if not family or not contig:
            raise SystemExit(
                f"Missing family/contig column in {fn}"
            )

        all_family_counts[family][sample] += 1

        if value(
            r,
            "assembly_supported_pass"
        ).upper() != "YES":
            continue

        supported_family_counts[family][sample] += 1
        passed.append(r)

        supported_rows.append({
            "sample": sample,
            "contig": contig,
            "family": family,
            "length": value(r, "length"),
            "mapped_reads": value(r, "mapped_reads"),
            "breadth1_pct": value(r, "breadth1_pct"),
            "mean_depth_mapq10": value(
                r,
                "mean_depth_mapq10",
                "mean_depth",
            ),
            "assembly_supported_pass": "YES",
        })

    strict_n = len(rows)
    supported_n = len(passed)

    pass_pct = (
        100.0 * supported_n / strict_n
        if strict_n
        else 0.0
    )

    summary_rows.append({
        "sample": sample,
        "strict_target_contigs": strict_n,
        "supported_target_contigs": supported_n,
        "pass_pct": f"{pass_pct:.1f}",
    })


# ------------------------------------------------------------
# Supported contigs
# ------------------------------------------------------------

with SUPPORTED_OUT.open("w") as f:

    fields = [
        "sample",
        "contig",
        "family",
        "length",
        "mapped_reads",
        "breadth1_pct",
        "mean_depth_mapq10",
        "assembly_supported_pass",
    ]

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(supported_rows)


# ------------------------------------------------------------
# Per-sample summary
# ------------------------------------------------------------

with SUMMARY_OUT.open("w") as f:

    fields = [
        "sample",
        "strict_target_contigs",
        "supported_target_contigs",
        "pass_pct",
    ]

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(summary_rows)


# ------------------------------------------------------------
# Family matrices
# ------------------------------------------------------------

observed = {
    r["family"]
    for r in supported_rows
}

all_observed = set(all_family_counts)

families = [
    x for x in TARGET_FAMILIES
    if x in observed or x in all_observed
]

for outfile, counts in [
    (ALL_FAMILY_OUT, all_family_counts),
    (SUPPORTED_FAMILY_OUT, supported_family_counts),
]:

    with outfile.open("w") as f:

        fields = ["family"] + samples

        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        w.writeheader()

        for family in families:

            row = {"family": family}

            for sample in samples:
                row[sample] = counts[family][sample]

            w.writerow(row)


print(
    f"Supported contigs: {len(supported_rows)}"
)

for r in summary_rows:
    print(
        r["sample"],
        r["strict_target_contigs"],
        r["supported_target_contigs"],
        r["pass_pct"],
        sep="\t",
    )

print("Saved:", SUPPORTED_OUT)
print("Saved:", SUMMARY_OUT)
print("Saved:", ALL_FAMILY_OUT)
print("Saved:", SUPPORTED_FAMILY_OUT)
