#!/usr/bin/env python3

import csv
import re
from collections import Counter
from pathlib import Path

inp = Path(
    "results_summary/multisample/"
    "all96_annotation_reconciled.tsv"
)

out = Path(
    "results_summary/multisample/"
    "all96_biological_grouping.tsv"
)

rows = list(
    csv.DictReader(
        inp.open(),
        delimiter="\t"
    )
)

# --------------------------------------------------
# Cluster-specific refinements supported by the
# earlier anchor / priority analysis.
# --------------------------------------------------
cluster_override = {

    # Astrovirus populations already characterized
    "Astroviridae_C001":
        ("HAstV-1-like population A", "HAstV1_A"),

    "Astroviridae_C002":
        ("HAstV-1-like population B", "HAstV1_B"),

    "Astroviridae_C005":
        ("HAstV-3-like", "HAstV3"),

    "Astroviridae_C006":
        ("HAstV-5-like", "HAstV5"),

    "Astroviridae_C008":
        ("HAstV-3-like", "HAstV3"),
}


def rotavirus_segment(title):
    """
    Conservative segment inference.

    Prefer explicit 'segment N'.
    For Rotavirus A only, recognize standard gene/segment
    relationships when the title contains a protein/gene name.
    """

    t = title.lower()

    m = re.search(r"segment\s+([0-9]{1,2})", t)
    if m:
        return f"segment_{m.group(1)}"

    # Standard RVA genome-segment relationships
    gene_map = [
        (r"\bvp1\b",  "segment_1"),
        (r"\bvp2\b",  "segment_2"),
        (r"\bvp3\b",  "segment_3"),
        (r"\bvp4\b",  "segment_4"),
        (r"\bnsp1\b", "segment_5"),
        (r"\bvp6\b",  "segment_6"),
        (r"\bnsp3\b", "segment_7"),
        (r"\bnsp2\b", "segment_8"),
        (r"\bvp7\b",  "segment_9"),
        (r"\bnsp4\b", "segment_10"),
        (r"\bnsp5\b", "segment_11"),
        (r"\bnsp6\b", "segment_11"),
    ]

    for pattern, seg in gene_map:
        if re.search(pattern, t):
            return seg

    return "segment_unresolved"


def classify(r):

    cid = r["cluster_id"]
    fam = r["reconciled_family"]
    name = r["best_nt_name"]
    title = r["best_nt_title"]
    quality = r["match_quality"]

    # ----------------------------------------------
    # Cluster-specific validated refinements
    # ----------------------------------------------
    if cid in cluster_override:
        label, group = cluster_override[cid]

        return {
            "biological_label": label,
            "biological_group": group,
            "grouping_level": "cluster_refined",
            "segment": "",
            "review_flag": "NO",
        }

    # ----------------------------------------------
    # Adenoviridae
    # ----------------------------------------------
    if fam == "Adenoviridae":

        if name == "Human adenovirus 41":
            label = "Human adenovirus 41-like"
            group = "HAdV41"

        elif name == "Human mastadenovirus D":
            label = "Human mastadenovirus D-like"
            group = "HAdV_D_like"

        else:
            label = f"{name}-like"
            group = "Adenovirus_unresolved"

        return {
            "biological_label": label,
            "biological_group": group,
            "grouping_level": "virus_like",
            "segment": "",
            "review_flag": "NO",
        }

    # ----------------------------------------------
    # Astroviridae
    # ----------------------------------------------
    if fam == "Astroviridae":

        mappings = {
            "Human astrovirus 1":
                ("HAstV-1-like", "HAstV1"),

            "Human astrovirus 2":
                ("HAstV-2-like", "HAstV2"),

            "Human astrovirus 3":
                ("HAstV-3-like", "HAstV3"),

            "Human astrovirus 4":
                ("HAstV-4-like", "HAstV4"),

            "Human astrovirus 8":
                ("HAstV-8-like", "HAstV8"),

            "Astrovirus MLB1":
                ("Astrovirus MLB1-like", "MLB1"),

            "Astrovirus MLB2":
                ("Astrovirus MLB2-like", "MLB2"),

            "Astrovirus VA5":
                ("Astrovirus VA5-like", "VA5"),

            "Jingmen rodent astrovirus 2":
                (
                    "Jingmen rodent astrovirus 2-like",
                    "Rodent_astrovirus_like"
                ),

            "Mamastrovirus 1":
                (
                    "Mamastrovirus 1-like",
                    "Mamastrovirus1_like"
                ),

            "Mamastrovirus 2":
                (
                    "Mamastrovirus 2-like",
                    "Mamastrovirus2_like"
                ),

            "Human astrovirus":
                (
                    "Human astrovirus-like, unresolved type",
                    "HAstV_unresolved"
                ),
        }

        if name in mappings:
            label, group = mappings[name]
            flag = (
                "YES"
                if "unresolved" in group.lower()
                else "NO"
            )

        else:
            label = f"{name}-like"
            group = "Astrovirus_unresolved"
            flag = "YES"

        return {
            "biological_label": label,
            "biological_group": group,
            "grouping_level": "virus_like",
            "segment": "",
            "review_flag": flag,
        }

    # ----------------------------------------------
    # Caliciviridae
    #
    # Collapse to genogroup where appropriate rather
    # than overclaim exact genotype from best hit.
    # ----------------------------------------------
    if fam == "Caliciviridae":

        # Important: test GII before GI because
        # "Norovirus GII".startswith("Norovirus GI") is True.
        if name.startswith("Norovirus GII"):
            label = "Norovirus GII-like"
            group = "NoV_GII"

        elif name.startswith("Norovirus GI"):
            label = "Norovirus GI-like"
            group = "NoV_GI"

        elif name.startswith("Sapovirus GII"):
            label = "Sapovirus GII-like"
            group = "SaV_GII"

        elif (
            name.startswith("Sapovirus GI")
            or "Sapovirus Hu/GI.2" in name
        ):
            label = "Sapovirus GI-like"
            group = "SaV_GI"

        elif name.startswith("Sapovirus GV"):
            label = "Sapovirus GV-like"
            group = "SaV_GV"

        else:
            label = f"{name}-like"
            group = "Calicivirus_unresolved"

        return {
            "biological_label": label,
            "biological_group": group,
            "grouping_level": "genogroup_like",
            "segment": "",
            "review_flag":
                "YES" if group == "Calicivirus_unresolved" else "NO",
        }

    # ----------------------------------------------
    # Parvoviridae
    # ----------------------------------------------
    if fam == "Parvoviridae":

        mappings = {
            "Human bocavirus 2":
                ("Human bocavirus 2-like", "HBoV2"),

            "Human bocavirus 3":
                ("Human bocavirus 3-like", "HBoV3"),

            "Rodent bocavirus":
                ("Rodent bocavirus-like", "Rodent_BoV_like"),

            "Canine bocavirus":
                ("Canine bocavirus-like", "Canine_BoV_like"),

            "Porcine bocavirus":
                ("Porcine bocavirus-like", "Porcine_BoV_like"),

            "Rattus rat parvovirus":
                (
                    "Rattus rat parvovirus-like",
                    "Rat_parvovirus_like"
                ),
        }

        if name in mappings:
            label, group = mappings[name]
        else:
            label = f"{name}-like"
            group = "Parvovirus_unresolved"

        return {
            "biological_label": label,
            "biological_group": group,
            "grouping_level": "virus_like",
            "segment": "",
            "review_flag":
                "YES" if group == "Parvovirus_unresolved" else "NO",
        }

    # ----------------------------------------------
    # Picornaviridae
    #
    # Strong = retain best-hit virus-like label.
    # Moderate = deliberately broaden annotation.
    # ----------------------------------------------
    if fam == "Picornaviridae":

        if quality == "moderate":

            if "Enterovirus B" in name or "Echovirus" in name:
                label = (
                    f"Divergent Enterovirus B-like "
                    f"({name} best hit)"
                )
                group = "Enterovirus_B_like"

            elif (
                "Coxsackievirus B" in name
            ):
                label = (
                    f"Divergent Enterovirus B-like "
                    f"({name} best hit)"
                )
                group = "Enterovirus_B_like"

            elif (
                "Coxsackievirus A" in name
                or "Enterovirus A" in name
            ):
                label = (
                    f"Divergent Enterovirus A-like "
                    f"({name} best hit)"
                )
                group = "Enterovirus_A_like"

            elif name == "aichivirus A1":
                label = "Divergent Aichivirus A-like"
                group = "Aichivirus_A_like"

            elif name == "Human enterovirus":
                label = "Divergent human enterovirus-like"
                group = "Enterovirus_unresolved"

            else:
                label = f"Divergent {name}-like"
                group = "Picornavirus_unresolved"

            level = "broad_divergent_group"

        else:

            label = f"{name}-like"

            replacements = {
                "aichivirus A1": "Aichivirus_A1_like",
                "Kobuvirus sp.": "Kobuvirus_like",
                "Salivirus A": "Salivirus_A_like",
                "Picornaviridae sp.": "Picornavirus_unresolved",
            }

            if name in replacements:
                group = replacements[name]
            else:
                group = (
                    name.replace(" ", "_")
                    .replace("/", "_")
                    + "_like"
                )

            level = "best_hit_like"

        return {
            "biological_label": label,
            "biological_group": group,
            "grouping_level": level,
            "segment": "",
            "review_flag":
                "YES"
                if "unresolved" in group.lower()
                else "NO",
        }

    # ----------------------------------------------
    # Sedoreoviridae
    # ----------------------------------------------
    if fam == "Sedoreoviridae":

        seg = rotavirus_segment(title)

        # Cluster-specific refinement:
        # RVA methyltransferase/capping enzyme = VP3,
        # encoded by genome segment 3.
        if cid == "Sedoreoviridae_C018":
            seg = "segment_3"

        if (
            name == "Rotavirus A"
            or name == "Human rotavirus A"
            or name.startswith("Rotavirus A human/")
            or name == "Human rotavirus P[8]"
        ):
            species = "RVA"
            label = "Rotavirus A-like"

        elif name == "Human rotavirus C":
            species = "RVC"
            label = "Rotavirus C-like"

        elif name == "Rotavirus alphagastroenteritidis":
            species = "Rotavirus_alpha_like"
            label = "Rotavirus alphagastroenteritidis-like"

        else:
            species = "Rotavirus_unresolved"
            label = f"{name}-like"

        if seg != "segment_unresolved":
            group = f"{species}_{seg}"
            label = (
                f"{label} | "
                f"{seg.replace('_', ' ')}"
            )
            flag = "NO"

        else:
            group = f"{species}_segment_unresolved"
            flag = "YES"

        return {
            "biological_label": label,
            "biological_group": group,
            "grouping_level": "species_segment",
            "segment": seg,
            "review_flag": flag,
        }

    # ----------------------------------------------
    # Fallback
    # ----------------------------------------------
    return {
        "biological_label": f"{name}-like",
        "biological_group": "Unresolved",
        "grouping_level": "unresolved",
        "segment": "",
        "review_flag": "YES",
    }


out_rows = []

for r in rows:

    g = classify(r)

    x = {
        "cluster_id": r["cluster_id"],
        "original_family": r["original_family"],
        "reconciled_family": r["reconciled_family"],
        "best_nt_name": r["best_nt_name"],
        "best_nt_accession": r["best_nt_accession"],
        "match_quality": r["match_quality"],
        "query_coverage_pct": r["query_coverage_pct"],
        "top_hsp_identity_pct":
            r["top_hsp_identity_pct"],
        "biological_label":
            g["biological_label"],
        "biological_group":
            g["biological_group"],
        "segment":
            g["segment"],
        "grouping_level":
            g["grouping_level"],
        "review_flag":
            g["review_flag"],
        "representative":
            r["representative"],
        "representative_length":
            r["representative_length"],
        "sample_count":
            r["sample_count"],
        "samples":
            r["samples"],
        "best_nt_title":
            r["best_nt_title"],
    }

    out_rows.append(x)


fields = list(out_rows[0])

with out.open("w") as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(out_rows)


print("Clusters:", len(out_rows))
print(
    "Biological groups:",
    len(set(x["biological_group"] for x in out_rows))
)

print(
    "Review flagged:",
    sum(x["review_flag"] == "YES" for x in out_rows)
)

print("\nGroup counts:")

counts = Counter(
    (
        x["reconciled_family"],
        x["biological_group"]
    )
    for x in out_rows
)

for (fam, group), n in sorted(counts.items()):
    print(f"{fam}\t{group}\t{n}")

print("\nSaved:", out)
