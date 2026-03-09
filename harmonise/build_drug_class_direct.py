#!/usr/bin/env python3
"""
build_drug_class_direct.py
==========================
Generate harmonise/drug_class_direct.tsv from the curated CURATED list below.

This file is a manual override table: drugs that harmonise.py cannot resolve
to a canonical class by source-name lookup alone (because the drug name in
CARD/ResFinder/NCBI doesn't match the class lookup key, the drug is
CARD-specific, or the assignment requires judgement).

Edit CURATED here and re-run to regenerate the TSV.
"""

import csv
from pathlib import Path

ROOT = Path(__file__).parent.parent
OUT = ROOT / "harmonise" / "drug_class_direct.tsv"

# ---------------------------------------------------------------------------
# Curated overrides
# Each entry: canonical_drug, canonical_class, notes
#
# canonical_drug   – INN drug name (lowercase), matching drug_canonical.tsv
# canonical_class  – CARD ARO drug class name, matching class_mapping.tsv
# notes            – rationale / cross-source reconciliation
# ---------------------------------------------------------------------------
CURATED = [
    # ---- aminoglycosides ----
    dict(canonical_drug="hygromycin b",
         canonical_class="aminoglycoside antibiotic",
         notes="aminoglycoside; research/veterinary use; INN: hygromycin B"),

    # ---- fluoroquinolones not resolved by class-name lookup ----
    dict(canonical_drug="levofloxacin",
         canonical_class="fluoroquinolone antibiotic",
         notes="3rd-gen fluoroquinolone; INN: levofloxacin"),
    dict(canonical_drug="moxifloxacin",
         canonical_class="fluoroquinolone antibiotic",
         notes="4th-gen fluoroquinolone; INN: moxifloxacin"),
    dict(canonical_drug="ofloxacin",
         canonical_class="fluoroquinolone antibiotic",
         notes="2nd-gen fluoroquinolone; INN: ofloxacin"),
    dict(canonical_drug="enrofloxacin",
         canonical_class="fluoroquinolone antibiotic",
         notes="veterinary fluoroquinolone; INN: enrofloxacin"),

    # ---- rifamycins ----
    dict(canonical_drug="rifabutin",
         canonical_class="rifamycin antibiotic",
         notes="rifamycin class; antituberculosis; INN: rifabutin"),

    # ---- nitrofurans ----
    dict(canonical_drug="nitrofurantoin",
         canonical_class="nitrofuran antibiotic",
         notes="INN: nitrofurantoin"),

    # ---- thioamides ----
    dict(canonical_drug="prothionamide",
         canonical_class="thioamide antibiotic",
         notes="antituberculosis thioamide; INN: prothionamide"),

    # ---- nitroimidazoles ----
    dict(canonical_drug="delamanid",
         canonical_class="nitroimidazole antibiotic",
         notes="nitroimidazole-oxazole; antituberculosis; INN: delamanid"),

    # ---- riminophenazines ----
    dict(canonical_drug="clofazimine",
         canonical_class="riminophenazine antibiotic",
         notes="anti-leprosy; antituberculosis; INN: clofazimine"),

    # ---- sulphones ----
    dict(canonical_drug="dapsone",
         canonical_class="sulphone antibiotic",
         notes="anti-leprosy; INN: dapsone; UK spelling: dapsone (same)"),

    # ---- elfamycins ----
    dict(canonical_drug="kirromycin",
         canonical_class="elfamycin antibiotic",
         notes="EF-Tu inhibitor; INN: kirromycin"),
    dict(canonical_drug="ge2270a",
         canonical_class="elfamycin antibiotic",
         notes="EF-Tu inhibitor; research compound; also known as MDL 62879"),
    dict(canonical_drug="pulvomycin",
         canonical_class="elfamycin antibiotic",
         notes="EF-Tu inhibitor; also known as GE2270A-related"),
    dict(canonical_drug="enacyloxin iia",
         canonical_class="elfamycin antibiotic",
         notes="EF-Tu inhibitor; polyketide; INN: enacyloxin"),
    dict(canonical_drug="elfamycin",
         canonical_class="elfamycin antibiotic",
         notes="class name used as drug entry in CARD; maps to own class"),

    # ---- peptide antibiotics ----
    dict(canonical_drug="edeine",
         canonical_class="peptide antibiotic",
         notes="translation inhibitor; research compound"),
    dict(canonical_drug="lysocin (e)",
         canonical_class="peptide antibiotic",
         notes="menaquinone-targeting cyclic peptide; research compound"),

    # ---- pactamycin-like ----
    dict(canonical_drug="pactamycin",
         canonical_class="pactamycin-like antibiotic",
         notes="translation inhibitor; INN: pactamycin"),

    # ---- zoliflodacin-like ----
    dict(canonical_drug="zoliflodacin",
         canonical_class="zoliflodacin-like antibiotic",
         notes="novel GyrB inhibitor; INN: zoliflodacin"),

    # ---- biocides ----
    dict(canonical_drug="triclosan",
         canonical_class="quaternary ammonium compound",
         notes="biocide/antiseptic; enoyl-ACP reductase inhibitor; context=biocide"),

    # ---- thiosemicarbazone (grouped after biocides to match original ordering) ----
    dict(canonical_drug="perchlozone",
         canonical_class="thioamide antibiotic",
         notes="thiosemicarbazone class; antituberculosis; also known as terizidone precursor"),

    # ---- streptogramins ----
    dict(canonical_drug="ostreogrycin",
         canonical_class="streptogramin antibiotic",
         notes="streptogramin mixture (ostreogrycin A = streptogramin A component; "
               "ostreogrycin B = streptogramin B component); INN: ostreogrycin"),

    # ---- oxazolidinone class token (placed here to match original ordering) ----
    dict(canonical_drug="oxazolidinone",
         canonical_class="oxazolidinone antibiotic",
         notes="appears as drug-level token in NCBI when gene confers class-wide resistance"),

    dict(canonical_drug="pristinamycin",
         canonical_class="streptogramin antibiotic",
         notes="mixture of pristinamycin IA (streptogramin B) + IIA (streptogramin A); "
               "INN: pristinamycin"),

    # ---- pleuromutilins (placed here to match original ordering) ----
    dict(canonical_drug="retapamulin",
         canonical_class="pleuromutilin antibiotic",
         notes="topical pleuromutilin; INN: retapamulin"),

    dict(canonical_drug="streptogramin a",
         canonical_class="streptogramin A antibiotic",
         notes="appears as drug-level token in NCBI; INN: streptogramin A component class"),
    dict(canonical_drug="streptogramin b",
         canonical_class="streptogramin B antibiotic",
         notes="appears as drug-level token in NCBI; INN: streptogramin B component class"),

    # ---- oxazolidinones ----
    dict(canonical_drug="tedizolid",
         canonical_class="oxazolidinone antibiotic",
         notes="second-generation oxazolidinone; INN: tedizolid"),

    # ---- phenicols ----
    dict(canonical_drug="thiamphenicol",
         canonical_class="phenicol antibiotic",
         notes="phenicol class; INN: thiamphenicol; not in ResFinder but related to chloramphenicol"),

    dict(canonical_drug="vernamycin b",
         canonical_class="streptogramin B antibiotic",
         notes="also known as virginiamycin S; streptogramin B component"),
    dict(canonical_drug="virginiamycin",
         canonical_class="streptogramin antibiotic",
         notes="mixture of virginiamycin M (streptogramin A) + virginiamycin S "
               "(streptogramin B); veterinary use"),

    # ---- combination products ----
    dict(canonical_drug="ceftazidime-avibactam",
         canonical_class="beta-lactam antibiotic",
         notes="CARD CZA; cephalosporin + BLI combination; class reflects the beta-lactam component"),
]


def main() -> None:
    fieldnames = ["canonical_drug", "canonical_class", "notes"]
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(CURATED)
    print(f"Wrote {len(CURATED)} rows to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
