#!/usr/bin/env python3

from pathlib import Path
import csv

ASSAYS = Path(
    "refs/single_sample/diagnostic_assays/astrovirus_public_assays.tsv"
)

SAMPLES = {
    "HAstV1_A":
        Path("results/single_sample/consensus/ERR14788990/HAstV1_A/"
             "HAstV1_A.consensus.fa"),

    "HAstV1_B":
        Path("results/single_sample/consensus/ERR14788990/HAstV1_B/"
             "HAstV1_B.consensus.fa"),

    "HAstV3":
        Path("results/single_sample/consensus/ERR14788990/HAstV3/"
             "HAstV3.consensus.fa"),
}

OUT = Path(
    "results/single_sample/diagnostic_mismatch/"
    "ERR14788990_HAstV_primer_probe_mismatch_v2.tsv"
)

# 1-based inclusive search windows.
# Both assays target the ORF1b/ORF2 junction.
WINDOWS = {
    "Classic_HAstV_Japan": (4100, 4700),
    "BCCDC_classic_HAstV": (4000, 4700),
}

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
    return "".join(
        COMP.get(x, "N")
        for x in reversed(s.upper())
    )


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def compare(oligo, target):

    mismatches = []
    unknown = []

    for i, (o, b) in enumerate(zip(oligo, target), 1):

        if b not in "ACGT":
            unknown.append(i)
            continue

        if b not in IUPAC[o]:
            mismatches.append(i)

    terminal_start = max(1, len(oligo) - 4)

    terminal_mm = [
        i for i in mismatches
        if i >= terminal_start
    ]

    return mismatches, unknown, terminal_mm


def best_hit(genome, oligo, role, lo, hi):

    L = len(oligo)

    lo0 = max(0, lo - 1)
    hi0 = min(len(genome), hi)

    candidates = []

    for start in range(lo0, hi0 - L + 1):

        genomic = genome[start:start+L]

        if role == "F":
            orientations = [("+", genomic)]

        elif role == "R":
            orientations = [("-", revcomp(genomic))]

        else:
            # Probe polarity differs between assays/sources;
            # evaluate both orientations.
            orientations = [
                ("+", genomic),
                ("-", revcomp(genomic)),
            ]

        for strand, target in orientations:

            mm, unknown, terminal = compare(
                oligo, target
            )

            # Unknown bases must NOT behave as perfect matches.
            #
            # 1 unknown is penalized more heavily than
            # a normal internal mismatch.
            penalty = (
                len(mm)
                + 3 * len(unknown)
                + 2 * len(terminal)
            )

            score = (
                penalty,
                len(unknown),
                len(terminal),
                len(mm),
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

    if not candidates:
        raise RuntimeError("No candidate windows")

    return min(candidates, key=lambda x: x[0])


with open(ASSAYS) as f:
    assays = list(csv.DictReader(f, delimiter="\t"))


header = [
    "sample",
    "assay",
    "role",
    "oligo",
    "oligo_sequence",
    "strand",
    "genome_start",
    "genome_end",
    "target_sequence",
    "mismatches",
    "mismatch_positions_in_oligo",
    "uncallable_positions",
    "primer_3prime_mismatches_last5",
]

rows = []

for sample, fasta in SAMPLES.items():

    genome = read_fasta(fasta)

    for a in assays:

        assay = a["assay"]
        role = a["role"]
        oligo = a["sequence"].upper()

        lo, hi = WINDOWS[assay]

        (
            score,
            start,
            end,
            strand,
            target,
            mm,
            unknown,
            terminal,
        ) = best_hit(
            genome,
            oligo,
            role,
            lo,
            hi,
        )

        rows.append([
            sample,
            assay,
            role,
            a["name"],
            oligo,
            strand,
            start,
            end,
            target,
            len(mm),
            ",".join(map(str, mm))
                if mm else "none",
            ",".join(map(str, unknown))
                if unknown else "none",
            (
                len(terminal)
                if role in ("F", "R")
                else "NA"
            ),
        ])


with open(OUT, "w") as f:
    w = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writerow(header)
    w.writerows(rows)


print("\t".join(header))

for r in rows:
    print("\t".join(map(str, r)))

print()
print("Saved:", OUT)
