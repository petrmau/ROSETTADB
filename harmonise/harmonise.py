#!/usr/bin/env python3
"""
harmonise.py
============
Parse all three AMR data sources (ResFinder, CARD, NCBI AMRFinderPlus),
apply harmonisation rules, and emit:

  harmonise/drug_canonical.tsv   – one row per canonical drug
  harmonise/drug_alias.tsv       – spelling variants, brand names, source aliases
  harmonise/drug_class_member.tsv – drug → canonical drug class (many-to-many)
  harmonise/gene_drug_link.tsv   – gene → canonical drug + class (from NCBI)

Canonical drug names follow INN (International Nonproprietary Name).
Canonical class names follow CARD ARO drug class names.
"""

import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import parse_aro_obo

ROOT = Path(__file__).parent.parent

# ---------------------------------------------------------------------------
# INN normalisation helpers
# ---------------------------------------------------------------------------

# UK → INN corrections (non-exhaustive; covers known AMR variants)
UK_TO_INN = {
    "amoxycillin": "amoxicillin",
    "sulphamethoxazole": "sulfamethoxazole",
    "cephuroxime": "cefuroxime",
    "cephamycin": "cephamycin",   # same
    "rifampicin": "rifampicin",   # INN is rifampicin; US uses rifampin
    "rifampin": "rifampicin",     # US spelling → INN
    "para-aminosalicyclic acid": "para-aminosalicylic acid",  # typo in ResFinder
}

# Explicit INN overrides for names that appear differently across sources
SOURCE_TO_INN = {
    # CARD shortname → INN
    "hyrgomycin b": "hygromycin b",   # typo in CARD ("hyrgomycin")
    "kasugamicin": "kasugamycin",     # CARD typo
    # NCBI subclass tokens → INN
    "rifampin": "rifampicin",
    "methicillin": "meticillin",      # INN is meticillin; keep methicillin as alias
}

def normalise_name(name: str) -> str:
    """Lowercase, strip, collapse whitespace, apply UK→INN and typo fixes."""
    name = name.strip().lower()
    name = re.sub(r"\s+", " ", name)
    name = UK_TO_INN.get(name, name)
    name = SOURCE_TO_INN.get(name, name)
    return name


def split_ncbi_combo(value: str) -> list[str]:
    """Split NCBI slash-delimited multi-drug/class strings into individual tokens."""
    return [normalise_name(t) for t in value.split("/") if t.strip()]


# ---------------------------------------------------------------------------
# Load class mapping (our curated TSV)
# ---------------------------------------------------------------------------

def load_class_mapping() -> dict[str, dict]:
    """
    Returns:
      mapping: dict keyed by canonical_name
      resfinder_lookup: resfinder_alias (lower) → canonical_name
      ncbi_lookup: ncbi_alias token (lower) → canonical_name
      card_class_abbrevs: set of CARD abbreviations that represent classes (not drugs)
    """
    mapping = {}
    resfinder_lookup: dict[str, str] = {}
    ncbi_lookup: dict[str, str] = {}
    card_class_abbrevs: set[str] = set()

    with open(ROOT / "harmonise/class_mapping.tsv") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cn = row["canonical_name"].strip()
            mapping[cn] = row

            rf = row["resfinder_alias"].strip().lower()
            if rf:
                resfinder_lookup[rf] = cn

            for token in row["ncbi_alias"].strip().lower().split("/"):
                token = token.strip()
                if token:
                    ncbi_lookup[token] = cn

            cca = row.get("card_class_abbrev", "").strip()
            if cca:
                card_class_abbrevs.add(cca)

    return mapping, resfinder_lookup, ncbi_lookup, card_class_abbrevs


def load_drug_class_direct() -> list[dict]:
    """Load manually curated drug → class direct mappings."""
    rows = []
    path = ROOT / "harmonise/drug_class_direct.tsv"
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Parse ResFinder
# ---------------------------------------------------------------------------

def parse_resfinder(resfinder_lookup: dict[str, str]) -> tuple[list, list]:
    """
    Returns:
      drugs: list of (canonical_drug_name, source_name, resfinder_class_raw)
      aliases: list of (canonical_drug_name, alias, alias_type, source)
    """
    drugs = []
    aliases = []
    path = ROOT / "sources/resfinder_db/antibiotic_classes.txt"

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            class_raw = parts[0].strip()

            # skip header and UNDER_DEVELOPMENT lines
            if class_raw in ("Class", "UNDER_DEVELOPMENT"):
                continue

            # individual drug names start from index 2 (index 1 = "Unknown X")
            drug_names = parts[2:] if len(parts) > 2 else []
            for raw in drug_names:
                raw = raw.strip()
                if not raw:
                    continue

                # Handle combination products: store as-is + record components
                if "+" in raw:
                    canonical = normalise_name(raw)
                    components = [normalise_name(c) for c in raw.split("+")]
                    drugs.append({
                        "canonical_name": canonical,
                        "source_name": raw,
                        "resfinder_class": class_raw,
                        "is_combination": True,
                        "components": components,
                        "context": "clinical",
                    })
                    aliases.append((canonical, raw, "source_name", "resfinder"))
                    # Also emit each component as a standalone drug entry so it
                    # appears in drug_canonical (satisfies FK in drug_class_member)
                    for comp in components:
                        ctx = _context_flag(comp)
                        drugs.append({
                            "canonical_name": comp,
                            "source_name": comp,
                            "resfinder_class": class_raw,
                            "is_combination": False,
                            "components": [],
                            "context": ctx,
                        })
                else:
                    canonical = normalise_name(raw)
                    ctx = _context_flag(canonical)
                    drugs.append({
                        "canonical_name": canonical,
                        "source_name": raw,
                        "resfinder_class": class_raw,
                        "is_combination": False,
                        "components": [],
                        "context": ctx,
                    })
                    if raw.lower() != canonical:
                        aliases.append((canonical, raw, "source_name", "resfinder"))

    return drugs, aliases


# ---------------------------------------------------------------------------
# Parse CARD shortname_antibiotics.tsv
# ---------------------------------------------------------------------------

def parse_card_shortnames(card_class_abbrevs: set[str]) -> tuple[list, list]:
    """
    Returns drugs and aliases from CARD shortname TSV.
    Skips entries whose abbreviation is a class-level abbreviation (e.g. AMG, BLA, MAC)
    — those represent drug classes, not individual molecules.
    """
    drugs = []
    aliases = []
    seen_abbrevs = set()

    # Additional entries that are class names used as molecule names in CARD
    CLASS_MOLECULE_NAMES = {
        "aminoglycosides", "beta-lactams", "fluoroquinolones", "macrolides",
        "sulfonamides", "aminocoumarin", "multiple antibiotics",
    }

    path = ROOT / "sources/CARD/shortname_antibiotics.tsv"
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            abbrev = row["AAC Abbreviation"].strip()
            molecule = row["Molecule"].strip()

            # deduplicate (CAP appears twice)
            if abbrev in seen_abbrevs:
                continue
            seen_abbrevs.add(abbrev)

            # Skip class-level abbreviations
            if abbrev in card_class_abbrevs:
                continue

            canonical = normalise_name(molecule)

            # Skip if the molecule name is itself a class name
            if canonical in CLASS_MOLECULE_NAMES:
                continue

            # Handle & -separated combinations (e.g. "Ethambutol & Capreomycin")
            is_combo = " & " in molecule
            components = [normalise_name(c) for c in molecule.split(" & ")] if is_combo else []

            ctx = "clinical" if is_combo else _context_flag(canonical)
            drugs.append({
                "canonical_name": canonical,
                "card_abbrev": abbrev,
                "source_name": molecule,
                "context": ctx,
                "is_combination": is_combo,
                "components": components,
            })
            if molecule.lower() != canonical:
                aliases.append((canonical, molecule, "source_name", "card"))
            aliases.append((canonical, abbrev, "card_abbreviation", "card"))

    return drugs, aliases


# ---------------------------------------------------------------------------
# Parse NCBI AMRFinderPlus JSONL
# ---------------------------------------------------------------------------

def parse_ncbi() -> tuple[list, list, list]:
    """
    Returns:
      drugs: individual drug/class tokens extracted from subclass field
      aliases: alias rows
      gene_links: list of dicts for gene_drug_link.tsv
    """
    drugs = []
    aliases = []
    gene_links = []

    path = ROOT / "sources/amr_finder_plus/ncbi_dataset/data/data_report.jsonl"
    with open(path) as f:
        for line in f:
            record = json.loads(line)
            gene_name = record.get("geneFamily", "")
            accession = (
                (record.get("refseqProtein") or {}).get("accessionVersion", "")
                or (record.get("refseqNucleotide") or {}).get("accessionVersion", "")
            )
            class_raw = record.get("class", "")
            subclass_raw = record.get("subclass", "")
            element_type = record.get("type", "")

            # Split multi-class and multi-drug slash combos
            class_tokens = split_ncbi_combo(class_raw) if class_raw else []
            subclass_tokens = split_ncbi_combo(subclass_raw) if subclass_raw else []

            for token in subclass_tokens:
                if not token:
                    continue
                ctx = _context_flag(token)
                drugs.append({
                    "canonical_name": token,
                    "source_name": token,
                    "ncbi_subclass_raw": subclass_raw,
                    "context": ctx,
                })

            # gene → drug/class links
            for cls in class_tokens:
                for sub in subclass_tokens if subclass_tokens else [""]:
                    gene_links.append({
                        "gene_name": gene_name,
                        "accession": accession,
                        "element_type": element_type,
                        "ncbi_class_raw": class_raw,
                        "ncbi_subclass_raw": subclass_raw,
                        "canonical_class_token": cls,
                        "canonical_drug_token": sub,
                    })

    return drugs, aliases, gene_links


# ---------------------------------------------------------------------------
# Context flags for known non-standard entries
# ---------------------------------------------------------------------------

RESEARCH_TOOLS = {"g418"}
VETERINARY_DRUGS = {"narasin", "salinomycin", "maduramicin", "tylosin", "tiamulin",
                    "avilamycin", "thiostrepton", "virginiamycin", "virginiamycin m",
                    "virginiamycin s", "pristinamycin ia", "pristinamycin iia",
                    "ostreogrycin", "carbomycin"}
NON_THERAPEUTIC = {"bleomycin"}
INHIBITORS = {"clavulanic acid", "tazobactam", "avibactam", "sulbactam",
              "taniborbactam", "xeruborbactam", "relebactam"}
ANTITUBERCULOSIS = {"isoniazid", "pyrazinamide", "ethambutol", "rifampicin",
                    "ethionamide", "prothionamide", "para-aminosalicylic acid",
                    "delamanid", "bedaquiline", "capreomycin", "viomycin",
                    "clofazimine", "cycloserine", "d-cycloserine"}
BIOCIDES = {"triclosan", "quaternary ammonium"}

def _context_flag(name: str) -> str:
    if name in RESEARCH_TOOLS:
        return "research_tool"
    if name in NON_THERAPEUTIC:
        return "non_therapeutic"
    if name in INHIBITORS:
        return "inhibitor"
    if name in BIOCIDES:
        return "biocide"
    if name in VETERINARY_DRUGS:
        return "veterinary"
    if name in ANTITUBERCULOSIS:
        return "antituberculosis"
    return "clinical"


# ---------------------------------------------------------------------------
# Merge and deduplicate drug canonical table
# ---------------------------------------------------------------------------

def merge_drugs(*drug_lists) -> list[dict]:
    """Merge drug entries from multiple sources; deduplicate by canonical_name."""
    seen: dict[str, dict] = {}
    for drug_list in drug_lists:
        for d in drug_list:
            cn = d["canonical_name"]
            if cn not in seen:
                seen[cn] = {
                    "canonical_name": cn,
                    "is_combination": d.get("is_combination", False),
                    "components": d.get("components", []),
                    "context": d.get("context", "clinical"),
                    "card_abbrev": d.get("card_abbrev", ""),
                    "atc_code": "",         # to be filled by ATC lookup step
                    "inchikey": "",         # to be filled by PubChem lookup step
                    "pubchem_cid": "",
                    "chebi_id": "",
                    "sources": set(),
                }
            entry = seen[cn]
            # track which sources mention this drug
            src = d.get("card_abbrev") and "card" or \
                  d.get("ncbi_subclass_raw") and "ncbi" or \
                  d.get("resfinder_class") and "resfinder" or \
                  d.get("source") or "unknown"
            entry["sources"].add(src)
            # prefer richer context flag (non-clinical beats clinical)
            if entry["context"] == "clinical" and d.get("context", "clinical") != "clinical":
                entry["context"] = d["context"]
            if not entry["card_abbrev"] and d.get("card_abbrev"):
                entry["card_abbrev"] = d["card_abbrev"]

    # serialise sources set
    result = []
    for entry in seen.values():
        entry["sources"] = "|".join(sorted(entry["sources"]))
        entry["components"] = "|".join(entry["components"])
        result.append(entry)

    return sorted(result, key=lambda x: x["canonical_name"])


# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------

def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows):>5} rows → {path.relative_to(ROOT)}")


def write_alias_tsv(path: Path, aliases: list[tuple]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["canonical_name", "alias", "alias_type", "source"])
        for row in sorted(set(aliases)):
            writer.writerow(row)
    print(f"  Wrote {len(set(aliases)):>5} rows → {path.relative_to(ROOT)}")


# ---------------------------------------------------------------------------
# Build drug_class_member table
# ---------------------------------------------------------------------------

def build_class_membership(
    resfinder_drugs: list[dict],
    class_mapping: dict,
    resfinder_lookup: dict[str, str],
    ncbi_gene_links: list[dict],
    ncbi_lookup: dict[str, str],
    direct_mappings: list[dict],
) -> list[dict]:
    rows = []

    # From ResFinder: drug → class via resfinder_lookup
    for d in resfinder_drugs:
        rf_class = d.get("resfinder_class", "").lower()
        canonical_class = resfinder_lookup.get(rf_class)
        if canonical_class:
            cm = class_mapping.get(canonical_class, {})
            rows.append({
                "canonical_drug": d["canonical_name"],
                "canonical_class": canonical_class,
                "aro_accession": cm.get("aro_accession", ""),
                "category": cm.get("category", ""),
                "source": "resfinder",
            })
        # Also map combination product components
        for component in [c for c in d.get("components", []) if c]:
            canonical_class_c = resfinder_lookup.get(rf_class)
            if canonical_class_c:
                rows.append({
                    "canonical_drug": component,
                    "canonical_class": canonical_class_c,
                    "aro_accession": class_mapping.get(canonical_class_c, {}).get("aro_accession", ""),
                    "category": class_mapping.get(canonical_class_c, {}).get("category", ""),
                    "source": "resfinder_combination_component",
                })

    # From NCBI gene links: subclass drug token → class token
    # Only infer drug→class membership when the gene has a SINGLE class
    # (no slash in class field), to avoid false cross-pairings like
    # chloramphenicol → aminoglycoside antibiotic from multi-class entries.
    for link in ncbi_gene_links:
        cls_token = link["canonical_class_token"]
        drug_token = link["canonical_drug_token"]
        if not drug_token or not cls_token:
            continue
        # Skip if the original class field contained slashes (multi-class gene)
        if "/" in link["ncbi_class_raw"]:
            continue
        canonical_class = ncbi_lookup.get(cls_token)
        if canonical_class:
            cm = class_mapping.get(canonical_class, {})
            rows.append({
                "canonical_drug": drug_token,
                "canonical_class": canonical_class,
                "aro_accession": cm.get("aro_accession", ""),
                "category": cm.get("category", ""),
                "source": "ncbi",
            })

    # From manually curated direct drug→class mappings
    for dm in direct_mappings:
        drug = dm["canonical_drug"].strip()
        cls = dm["canonical_class"].strip()
        if not drug or not cls:
            continue
        cm = class_mapping.get(cls, {})
        rows.append({
            "canonical_drug": drug,
            "canonical_class": cls,
            "aro_accession": cm.get("aro_accession", ""),
            "category": cm.get("category", ""),
            "source": "curated",
        })

    # Deduplicate
    seen = set()
    deduped = []
    for r in rows:
        key = (r["canonical_drug"], r["canonical_class"])
        if key not in seen:
            seen.add(key)
            deduped.append(r)

    return sorted(deduped, key=lambda x: (x["canonical_class"], x["canonical_drug"]))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading class mapping …")
    class_mapping, resfinder_lookup, ncbi_lookup, card_class_abbrevs = load_class_mapping()

    print("Loading direct drug→class mappings …")
    direct_mappings = load_drug_class_direct()

    print("Parsing ResFinder …")
    rf_drugs, rf_aliases = parse_resfinder(resfinder_lookup)

    print("Parsing CARD shortnames …")
    card_drugs, card_aliases = parse_card_shortnames(card_class_abbrevs)

    print("Parsing NCBI AMRFinderPlus …")
    ncbi_drugs, ncbi_aliases, ncbi_gene_links = parse_ncbi()

    aro_drug_members, aro_gene_links = parse_aro_obo.parse()

    # Convert ARO drug→class entries into drug dicts for merge_drugs()
    aro_drugs = [
        {
            "canonical_name": r["canonical_drug"],
            "source_name":    r["canonical_drug"],
            "context":        _context_flag(r["canonical_drug"]),
            "is_combination": False,
            "components":     [],
            "source":         "CARD.obo",
        }
        for r in aro_drug_members
    ]

    print("Merging drug canonical table …")
    all_drugs = merge_drugs(rf_drugs, card_drugs, ncbi_drugs, aro_drugs)

    print("Building class membership …")
    membership = build_class_membership(
        rf_drugs, class_mapping, resfinder_lookup,
        ncbi_gene_links, ncbi_lookup, direct_mappings
    )

    # Append ARO-derived drug→class memberships (converted to standard format)
    seen_membership = {(r["canonical_drug"], r["canonical_class"]) for r in membership}
    for r in aro_drug_members:
        key = (r["canonical_drug"], r["canonical_class"])
        if key not in seen_membership:
            seen_membership.add(key)
            cm = class_mapping.get(r["canonical_class"], {})
            membership.append({
                "canonical_drug": r["canonical_drug"],
                "canonical_class": r["canonical_class"],
                "aro_accession":  r["class_aro_accession"],
                "category":       cm.get("category", ""),
                "source":         "aro_obo",
            })
    membership.sort(key=lambda x: (x["canonical_class"], x["canonical_drug"]))

    print("Building alias table …")
    all_aliases = rf_aliases + card_aliases + ncbi_aliases

    out = ROOT / "harmonise"

    print("\nWriting outputs:")
    write_tsv(
        out / "drug_canonical.tsv",
        all_drugs,
        ["canonical_name", "is_combination", "components", "context",
         "card_abbrev", "atc_code", "inchikey", "pubchem_cid", "chebi_id", "sources"],
    )

    write_alias_tsv(out / "drug_alias.tsv", all_aliases)

    write_tsv(
        out / "drug_class_member.tsv",
        membership,
        ["canonical_drug", "canonical_class", "aro_accession", "category", "source"],
    )

    write_tsv(
        out / "gene_drug_link.tsv",
        ncbi_gene_links,
        ["gene_name", "accession", "element_type",
         "ncbi_class_raw", "ncbi_subclass_raw",
         "canonical_class_token", "canonical_drug_token"],
    )

    write_tsv(
        out / "aro_gene_class.tsv",
        aro_gene_links,
        ["aro_accession", "gene_name", "canonical_class",
         "class_aro_accession", "source"],
    )

    # Summary stats
    aro_new_drugs = len({r["canonical_drug"] for r in aro_drug_members})
    aro_new_links = sum(1 for r in membership if r.get("source") == "aro_obo")
    print(f"\nSummary:")
    print(f"  Canonical drug classes : {len(class_mapping)}")
    print(f"  Canonical drugs        : {len(all_drugs)}")
    print(f"    found in ARO OBO     : {aro_new_drugs}")
    print(f"  Drug aliases           : {len(set(all_aliases))}")
    print(f"  Drug→class links       : {len(membership)}")
    print(f"    of which ARO-sourced : {aro_new_links}")
    print(f"  Gene→drug links (NCBI) : {len(ncbi_gene_links)}")
    print(f"  Gene→class (ARO OBO)   : {len(aro_gene_links)}")

    # Warn about drugs with no class assignment
    classed = {r["canonical_drug"] for r in membership}
    unclassed = [d["canonical_name"] for d in all_drugs
                 if d["canonical_name"] not in classed
                 and d["context"] not in ("inhibitor", "research_tool")
                 and not d["is_combination"]]
    if unclassed:
        print(f"\n  WARN: {len(unclassed)} drugs have no class assignment (review needed):")
        for u in sorted(unclassed)[:30]:
            print(f"    - {u}")
        if len(unclassed) > 30:
            print(f"    … and {len(unclassed)-30} more")


if __name__ == "__main__":
    main()
