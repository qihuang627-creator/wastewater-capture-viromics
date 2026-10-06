#!/usr/bin/env python3

from pathlib import Path
import csv

ASSAYS = Path("refs/diagnostic_assays/astrovirus_public_assays.tsv")

SAMPLES = {
    "HAstV1_A":
        Path("results/consensus/ERR14788990/HAstV1_A/HAstV1_A.consensus.fa"),

    "HAstV1_B":
        Path("results/consensus/ERR14788990/HAstV1_B/HAstV1_B.consensus.fa"),

    "HAstV3":
        Path("results/consensus/ERR14788990/HAstV3/HAstV3.consensus.fa"),
}

OUT = Path(
    "results/diagnostic_mismatch/"
    "ERR14788990_HAstV_primer_probe_mismatch.tsv"
)

IUPAC = {
    "A": set("A"),
    "C": set("C"),
    "G": set("G"),
    "T": set("T"),
    "R": set("AG"),
    "Y": set("CT"),
    "S": set("GC"),
    "W": set("AT"),
    "K": set("GT"),
    "M": set("AC"),
    "B": set("CGT"),
    "D": set("AGT"),
    "H": set("ACT"),
    "V": set("ACG"),
    "N": set("ACGT"),
}

COMP = {
    "A":"T", "C":"G", "G":"C", "T":"A",
    "R":"Y", "Y":"R", "S":"S", "W":"W",
    "K":"M", "M":"K", "B":"V", "V":"B",
    "D":"H", "H":"D", "N":"N"
}

def revcomp(s):
    return "".join(COMP[x] for x in reversed(s.upper()))

def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

def compare(oligo, target):
    """
    oligo is kept in its original 5'->3' orientation.
    target is oriented to match that oligo.
    """
    mismatches = []
    uncallable = []

    for i, (o, b) in enumerate(zip(oligo, target)):
        # Sample N = unknown, not mismatch
        if b not in "ACGT":
            uncallable.append(i + 1)
            continue

        if b not in IUPAC[o]:
            mismatches.append(i + 1)

    # Last 5 nt of primer = important 3' end
    terminal_start = max(1, len(oligo) - 4)

    terminal_mm = [
        x for x in mismatches
        if x >= terminal_start
    ]

    return mismatches, uncallable, terminal_mm


def best_hit(genome, oligo, role):

    L = len(oligo)
    candidates = []

    for start in range(0, len(genome) - L + 1):

        genomic_window = genome[start:start+L]

        # F primer sequence corresponds to plus-strand sequence.
        if role == "F":
            oriented = genomic_window
            orientations = [("+", oriented)]

        # Reverse primer binds plus strand; compare original primer
        # against reverse-complemented genomic window.
        elif role == "R":
            oriented = revcomp(genomic_window)
            orientations = [("-", oriented)]

        # Probe orientation can vary between published assays.
        else:
            orientations = [
                ("+", genomic_window),
                ("-", revcomp(genomic_window)),
            ]

        for strand, target in orientations:

            mm, unknown, terminal = compare(
                oligo, target
            )

            # Prefer:
            # fewer mismatches,
            # fewer unknown positions,
            # fewer primer 3'-terminal mismatches
            score = (
                len(mm),
                len(unknown),
                len(terminal)
            )

            candidates.append(
                (
                    score,
                    start + 1,
                    start + L,
                    strand,
                    target,
                    mm,
                    unknown,
                    terminal,
                )
            )

    return min(candidates, key=lambda x: x[0])


assays = []

with open(ASSAYS) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        assays.append(r)


header = [
    "sample",
    "assay",
    "role",
    "oligo",
    "oligo_sequence",
    "strand",
    "genome_start",
    "genome_end",
    "mismatches",
    "mismatch_positions_in_oligo",
    "uncallable_positions",
    "primer_3prime_mismatches_last5",
]

rows = []

for sample, fn in SAMPLES.items():

    genome = read_fasta(fn)

    for a in assays:

        oligo = a["sequence"].upper()
        role = a["role"]

        (
            score,
            start,
            end,
            strand,
            target,
            mm,
            unknown,
            terminal,
        ) = best_hit(genome, oligo, role)

        rows.append([
            sample,
            a["assay"],
            role,
            a["name"],
            oligo,
            strand,
            start,
            end,
            len(mm),
            ",".join(map(str, mm)) if mm else "none",
            ",".join(map(str, unknown)) if unknown else "none",
            len(terminal) if role in ("F","R") else "NA",
        ])


with open(OUT, "w") as f:
    w = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )
    w.writerow(header)
    w.writerows(rows)


print("\t".join(header))

for r in rows:
    print("\t".join(map(str, r)))

print()
print("Saved:", OUT)
