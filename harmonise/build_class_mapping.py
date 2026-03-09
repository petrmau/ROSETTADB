#!/usr/bin/env python3
"""
build_class_mapping.py
======================
Generate harmonise/class_mapping.tsv from curated class definitions
and live ARO accession lookups from sources/CARD/card.json.

The cross-source alias mappings (ResFinder ↔ NCBI ↔ CARD canonical names)
are curated knowledge stored in CURATED_CLASSES below.
ARO accessions are resolved from CARD JSON where the name matches exactly,
or from ARO_NAME_MAP where our canonical name differs from CARD's name,
and left blank for classes that have no CARD drug-class entry.

Re-run this script whenever CARD is updated to refresh ARO accessions.
"""

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
CARD_JSON = ROOT / "sources" / "CARD" / "card.json"
OUT = ROOT / "harmonise" / "class_mapping.tsv"

# ---------------------------------------------------------------------------
# Curated class definitions
# Fields: canonical_name, resfinder_alias, ncbi_alias,
#         card_abbrev, card_class_abbrev, category, notes
#
# canonical_name        – CARD ARO drug class name (or closest INN equivalent)
# resfinder_alias       – class label used in ResFinder antibiotic_classes.txt
# ncbi_alias            – single-token Class value from NCBI AMRFinderPlus
#                         (slash-delimited multi-class tokens are NOT aliases)
# card_abbrev           – CARD shortname abbreviation for a representative drug
# card_class_abbrev     – CARD shortname abbreviation for the class itself
# category              – antibiotic | biocide | inhibitor | resistance_mechanism
#                         | veterinary_ionophore | non_therapeutic
# notes                 – free-text; cross-source reconciliation decisions
# ---------------------------------------------------------------------------
CURATED_CLASSES = [
    # ---- aminoglycosides ----
    dict(canonical_name="aminoglycoside antibiotic",
         resfinder_alias="Aminoglycoside", ncbi_alias="AMINOGLYCOSIDE",
         card_abbrev="AMK", card_class_abbrev="AMG", category="antibiotic",
         notes="AMG is CARD class abbrev; individual drugs: AMK=amikacin etc."),

    # ---- aminocoumarins ----
    dict(canonical_name="aminocoumarin antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="AMU", category="antibiotic",
         notes="e.g. novobiocin"),

    # ---- beta-lactams (parent + subclasses) ----
    dict(canonical_name="beta-lactam antibiotic",
         resfinder_alias="Beta-lactam", ncbi_alias="BETA-LACTAM",
         card_abbrev="", card_class_abbrev="BLA", category="antibiotic",
         notes="parent class; includes penicillins, cephalosporins, carbapenems, monobactams"),
    dict(canonical_name="carbapenem",
         resfinder_alias="", ncbi_alias="CARBAPENEM",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="subclass of beta-lactam"),
    dict(canonical_name="cephalosporin",
         resfinder_alias="", ncbi_alias="CEPHALOSPORIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="subclass of beta-lactam"),
    dict(canonical_name="penicillin antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="subclass of beta-lactam; CARD name: penicillin beta-lactam"),
    dict(canonical_name="monobactam",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="subclass of beta-lactam; e.g. aztreonam"),

    # ---- cycloserine ----
    dict(canonical_name="cycloserine-like antibiotic",
         resfinder_alias="Analog of d-alanine", ncbi_alias="",
         card_abbrev="DCS", card_class_abbrev="", category="antibiotic",
         notes="INN: cycloserine"),

    # ---- diarylquinolines ----
    dict(canonical_name="diarylquinoline antibiotic",
         resfinder_alias="Diarylquinoline", ncbi_alias="",
         card_abbrev="BDQ", card_class_abbrev="", category="antibiotic",
         notes="e.g. bedaquiline (antituberculosis)"),

    # ---- quinolones / fluoroquinolones ----
    dict(canonical_name="fluoroquinolone antibiotic",
         resfinder_alias="Quinolone", ncbi_alias="QUINOLONE",
         card_abbrev="", card_class_abbrev="FLO", category="antibiotic",
         notes='ResFinder "Quinolone" covers both quinolones and fluoroquinolones; FLO is CARD class abbrev'),
    dict(canonical_name="quinolone antibiotic",
         resfinder_alias="Quinolone", ncbi_alias="QUINOLONE",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="true (non-fluorinated) quinolones e.g. nalidixic acid; ResFinder lumps with fluoroquinolones"),

    # ---- folate-pathway ----
    dict(canonical_name="diaminopyrimidine antibiotic",
         resfinder_alias="Folate pathway antagonist", ncbi_alias="TRIMETHOPRIM",
         card_abbrev="TMP", card_class_abbrev="", category="antibiotic",
         notes="trimethoprim class"),
    dict(canonical_name="sulfonamide antibiotic",
         resfinder_alias="Folate pathway antagonist", ncbi_alias="SULFONAMIDE",
         card_abbrev="", card_class_abbrev="SLF", category="antibiotic",
         notes="SLF is CARD class abbrev"),

    # ---- fosfomycin ----
    dict(canonical_name="phosphonic acid antibiotic",
         resfinder_alias="Fosfomycin", ncbi_alias="FOSFOMYCIN",
         card_abbrev="FOF", card_class_abbrev="", category="antibiotic",
         notes="INN: fosfomycin"),

    # ---- glycopeptides ----
    dict(canonical_name="glycopeptide antibiotic",
         resfinder_alias="Glycopeptide", ncbi_alias="GLYCOPEPTIDE",
         card_abbrev="VAN", card_class_abbrev="", category="antibiotic",
         notes="e.g. vancomycin, teicoplanin"),

    # ---- ionophores (veterinary) ----
    dict(canonical_name="ionophore with antibiotic activity",
         resfinder_alias="Ionophores", ncbi_alias="IONOPHORE",
         card_abbrev="", card_class_abbrev="", category="veterinary_ionophore",
         notes="e.g. narasin, salinomycin, maduramicin; veterinary use"),

    # ---- isoniazid ----
    dict(canonical_name="isonicotinic acid hydrazide",
         resfinder_alias="Isonicotinic acid hydrazide", ncbi_alias="",
         card_abbrev="INH", card_class_abbrev="", category="antibiotic",
         notes="INN: isoniazid; antituberculosis; CARD name: isoniazid-like antibiotic"),

    # ---- lincosamides ----
    dict(canonical_name="lincosamide antibiotic",
         resfinder_alias="Lincosamide", ncbi_alias="LINCOSAMIDE",
         card_abbrev="CLI", card_class_abbrev="", category="antibiotic",
         notes="e.g. clindamycin, lincomycin"),

    # ---- macrolides ----
    dict(canonical_name="macrolide antibiotic",
         resfinder_alias="Macrolide", ncbi_alias="MACROLIDE",
         card_abbrev="", card_class_abbrev="MAC", category="antibiotic",
         notes="MAC is CARD class abbrev"),

    # ---- nitroimidazoles ----
    dict(canonical_name="nitroimidazole antibiotic",
         resfinder_alias="Nitroimidazole", ncbi_alias="NITROIMIDAZOLE",
         card_abbrev="MTZ", card_class_abbrev="", category="antibiotic",
         notes="e.g. metronidazole, delamanid (antituberculosis)"),

    # ---- oxazolidinones ----
    dict(canonical_name="oxazolidinone antibiotic",
         resfinder_alias="Oxazolidinone", ncbi_alias="OXAZOLIDINONE",
         card_abbrev="OXZ", card_class_abbrev="", category="antibiotic",
         notes="e.g. linezolid"),

    # ---- phenicols ----
    dict(canonical_name="phenicol antibiotic",
         resfinder_alias="Amphenicol", ncbi_alias="PHENICOL",
         card_abbrev="CHL", card_class_abbrev="", category="antibiotic",
         notes='ResFinder uses "Amphenicol"; INN class = phenicol'),

    # ---- pleuromutilins ----
    dict(canonical_name="pleuromutilin antibiotic",
         resfinder_alias="Pleuromutilin", ncbi_alias="PLEUROMUTILIN",
         card_abbrev="PLM", card_class_abbrev="", category="antibiotic",
         notes="e.g. tiamulin, retapamulin"),

    # ---- peptides / polymyxins ----
    dict(canonical_name="peptide antibiotic",
         resfinder_alias="Polymyxin", ncbi_alias="COLISTIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes='ResFinder uses class "Polymyxin"; NCBI names the drug "COLISTIN" as class'),

    # ---- mupirocin ----
    dict(canonical_name="mupirocin-like antibiotic",
         resfinder_alias="Pseudomonic acid", ncbi_alias="MUPIROCIN",
         card_abbrev="MUP", card_class_abbrev="", category="antibiotic",
         notes='ResFinder "Pseudomonic acid" = mupirocin class'),

    # ---- rifamycins ----
    dict(canonical_name="rifamycin antibiotic",
         resfinder_alias="Rifamycin", ncbi_alias="RIFAMYCIN",
         card_abbrev="RIF", card_class_abbrev="", category="antibiotic",
         notes='ResFinder lists class as "Rifamycin" (same as drug name); INN drug = rifampicin (UK) / rifampin (US)'),

    # ---- salicylic acid / PAS ----
    dict(canonical_name="salicylic acid antibiotic",
         resfinder_alias="Salicylic acid - anti-folate", ncbi_alias="",
         card_abbrev="PAS", card_class_abbrev="", category="antibiotic",
         notes="INN: para-aminosalicylic acid; antituberculosis"),

    # ---- fusidic acid ----
    dict(canonical_name="fusidane antibiotic",
         resfinder_alias="Steroid antibacterial", ncbi_alias="FUSIDIC ACID",
         card_abbrev="FA", card_class_abbrev="", category="antibiotic",
         notes='ResFinder "Steroid antibacterial"; INN: fusidic acid'),

    # ---- streptogramins ----
    dict(canonical_name="streptogramin A antibiotic",
         resfinder_alias="Streptogramin A", ncbi_alias="STREPTOGRAMIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic", notes=""),
    dict(canonical_name="streptogramin B antibiotic",
         resfinder_alias="Streptogramin B", ncbi_alias="STREPTOGRAMIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic", notes=""),
    dict(canonical_name="streptogramin antibiotic",
         resfinder_alias="", ncbi_alias="STREPTOGRAMIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="parent class when A/B not specified"),

    # ---- tetracyclines ----
    dict(canonical_name="tetracycline antibiotic",
         resfinder_alias="Tetracycline", ncbi_alias="TETRACYCLINE",
         card_abbrev="TET", card_class_abbrev="", category="antibiotic", notes=""),

    # ---- thioamides ----
    dict(canonical_name="thioamide antibiotic",
         resfinder_alias="Thioamide", ncbi_alias="",
         card_abbrev="ETO", card_class_abbrev="", category="antibiotic",
         notes="INN: ethionamide; also prothionamide; antituberculosis"),

    # ---- pyrazinamide ----
    dict(canonical_name="pyrazinamide",
         resfinder_alias="Synthetic derivative of nicotinamide", ncbi_alias="",
         card_abbrev="PZA", card_class_abbrev="", category="antibiotic",
         notes="INN: pyrazinamide; antituberculosis; structurally unique"),

    # ---- ethambutol ----
    dict(canonical_name="ethambutol",
         resfinder_alias="Unspecified", ncbi_alias="",
         card_abbrev="EMB", card_class_abbrev="", category="antibiotic",
         notes='INN: ethambutol; antituberculosis; listed as "Unspecified" in ResFinder'),

    # ---- NCBI-only classes (no ResFinder equivalent) ----
    dict(canonical_name="bacitracin",
         resfinder_alias="", ncbi_alias="BACITRACIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="peptide antibiotic; no ResFinder class entry"),
    dict(canonical_name="lipopeptide antibiotic",
         resfinder_alias="", ncbi_alias="LIPOPEPTIDE",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. daptomycin"),
    dict(canonical_name="streptothricin antibiotic",
         resfinder_alias="", ncbi_alias="STREPTOTHRICIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic", notes=""),
    dict(canonical_name="tuberactinomycin",
         resfinder_alias="", ncbi_alias="TUBERACTINOMYCIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. viomycin, capreomycin; antituberculosis"),
    dict(canonical_name="nitrofuran antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="NIT", card_class_abbrev="", category="antibiotic",
         notes="e.g. nitrofurantoin"),
    dict(canonical_name="tetracenomycin antibiotic",
         resfinder_alias="", ncbi_alias="TETRACENOMYCIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic", notes=""),
    dict(canonical_name="orthosomycin antibiotic",
         resfinder_alias="", ncbi_alias="AVILAMYCIN",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. avilamycin; orthosomycin class"),
    dict(canonical_name="thiostrepton",
         resfinder_alias="", ncbi_alias="THIOSTREPTON",
         card_abbrev="", card_class_abbrev="", category="antibiotic", notes=""),
    dict(canonical_name="bleomycin",
         resfinder_alias="", ncbi_alias="BLEOMYCIN",
         card_abbrev="", card_class_abbrev="", category="non_therapeutic",
         notes="cytotoxic glycopeptide; not a clinical antibiotic; AMR genes real"),

    # ---- biocides / disinfectants ----
    dict(canonical_name="quaternary ammonium compound",
         resfinder_alias="", ncbi_alias="QUATERNARY AMMONIUM",
         card_abbrev="", card_class_abbrev="", category="biocide",
         notes="disinfectant/antiseptic class; not an antibiotic"),
    dict(canonical_name="disinfecting agents and antiseptics",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="biocide",
         notes="CARD broad biocide category; covers triclosan, QACs, etc."),

    # ---- resistance mechanism pseudo-class ----
    dict(canonical_name="efflux",
         resfinder_alias="", ncbi_alias="EFFLUX",
         card_abbrev="", card_class_abbrev="", category="resistance_mechanism",
         notes="NCBI uses as pseudo-class; actually a resistance mechanism — exclude from drug class joins"),

    # ---- inhibitors ----
    dict(canonical_name="beta-lactamase inhibitor",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="AVI", card_class_abbrev="", category="inhibitor",
         notes="e.g. avibactam, clavulanic acid, tazobactam; not antibiotics themselves"),

    # ---- CARD-only / less common classes ----
    dict(canonical_name="riminophenazine antibiotic",
         resfinder_alias="Iminophenazine", ncbi_alias="",
         card_abbrev="CFZ", card_class_abbrev="", category="antibiotic",
         notes="e.g. clofazimine (anti-leprosy, antituberculosis)"),
    dict(canonical_name="sulphone antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="DAO", card_class_abbrev="", category="antibiotic",
         notes="e.g. dapsone (anti-leprosy); note: sulfone (INN) = sulphone (UK); CARD name: sulfone antibiotic"),
    dict(canonical_name="elfamycin antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="ELF", card_class_abbrev="", category="antibiotic",
         notes="e.g. kirromycin, GE2270A, pulvomycin, enacyloxin IIa; EF-Tu inhibitors"),
    dict(canonical_name="pactamycin-like antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="PAC", card_class_abbrev="", category="antibiotic",
         notes="e.g. pactamycin; translation inhibitor"),
    dict(canonical_name="macrocyclic antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. fidaxomicin; CARD has this class"),
    dict(canonical_name="zoliflodacin-like antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="ZOL", card_class_abbrev="", category="antibiotic",
         notes="e.g. zoliflodacin; DNA gyrase inhibitor, novel class"),
    dict(canonical_name="triazaacenaphthylene antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="triclosan target class"),
    dict(canonical_name="nucleoside antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. capreomycin, tuberactinomycins; some overlap with tuberactinomycin class"),
    dict(canonical_name="glycylcycline",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. tigecycline; subclass of tetracyclines in CARD"),
    dict(canonical_name="polyamine antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="CARD class; polyamine-based antibiotics"),
    dict(canonical_name="thiosemicarbazone antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. thiacetazone, perchlozone; antituberculosis"),
    dict(canonical_name="pyrazine antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. pyrazinamide; antituberculosis"),
    dict(canonical_name="moenomycin antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. moenomycin A; phosphoglycolipid; PBP glycosyltransferase inhibitor"),
    dict(canonical_name="bicyclomycin-like antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. bicyclomycin (bicozamycin); Rho factor inhibitor"),
    dict(canonical_name="nybomycin-like antibiotic",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. nybomycin; quinolone-related natural product"),
    dict(canonical_name="antibacterial free fatty acids",
         resfinder_alias="", ncbi_alias="",
         card_abbrev="", card_class_abbrev="", category="antibiotic",
         notes="e.g. MCFA-based antibacterials; CARD-specific class"),
]

# ---------------------------------------------------------------------------
# Canonical name → CARD drug class name (where they differ)
# ---------------------------------------------------------------------------
ARO_NAME_MAP = {
    "penicillin antibiotic":      "penicillin beta-lactam",
    "isonicotinic acid hydrazide": "isoniazid-like antibiotic",
    "sulphone antibiotic":        "sulfone antibiotic",
    "riminophenazine antibiotic": "riminophenazine antibiotic",  # may not exist in CARD
}

# ---------------------------------------------------------------------------
# Load CARD JSON and build drug-class name → ARO accession lookup
# ---------------------------------------------------------------------------

def load_card_drug_classes(card_json: Path) -> dict[str, str]:
    """Return {drug_class_name: 'ARO:XXXXXXX'} from CARD JSON."""
    data = json.loads(card_json.read_text())
    classes: dict[str, str] = {}
    for v in data.values():
        if not isinstance(v, dict):
            continue
        for cat in v.get("ARO_category", {}).values():
            if isinstance(cat, dict) and cat.get("category_aro_class_name") == "Drug Class":
                acc = "ARO:" + cat["category_aro_accession"]
                name = cat["category_aro_name"]
                classes[name] = acc
    return classes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    card_classes = load_card_drug_classes(CARD_JSON)

    fieldnames = [
        "canonical_name", "aro_accession",
        "resfinder_alias", "ncbi_alias",
        "card_abbrev", "card_class_abbrev",
        "category", "notes",
    ]

    rows = []
    missing_aro: list[str] = []

    for entry in CURATED_CLASSES:
        cname = entry["canonical_name"]
        card_lookup_name = ARO_NAME_MAP.get(cname, cname)
        aro = card_classes.get(card_lookup_name, "")
        if not aro:
            missing_aro.append(f"  {cname!r} (looked up as {card_lookup_name!r})")
        rows.append({**entry, "aro_accession": aro})

    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUT.relative_to(ROOT)}")
    if missing_aro:
        print(f"\n{len(missing_aro)} class(es) not found in CARD JSON "
              f"(ARO accession left blank):", file=sys.stderr)
        for m in missing_aro:
            print(m, file=sys.stderr)


if __name__ == "__main__":
    main()
