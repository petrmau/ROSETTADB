#!/usr/bin/env python3
"""
parse_aro_obo.py
================
Parse the CARD Antibiotic Resistance Ontology (ARO) OBO file to extract:

1. Drug → drug class is_a hierarchy  →  aro_drug_class_member.tsv
   One row per (drug, canonical_class) pair where the drug term's is_a
   ancestry includes a known canonical class from class_mapping.tsv.

2. Gene → drug class confers_resistance_to_drug_class edges  →  aro_gene_class.tsv
   One row per (gene ARO accession, canonical_class) pair.

Both outputs are keyed against canonical class names in class_mapping.tsv.
"""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).parent.parent
OBO_PATH = ROOT / "sources/CARD/aro.obo"


# ---------------------------------------------------------------------------
# Load canonical class lookups from class_mapping.tsv
# ---------------------------------------------------------------------------

def load_canonical_classes(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """
    Returns:
      aro_to_canonical : ARO accession  → canonical_name
      name_to_canonical: lowercase name → canonical_name
    """
    aro_to_canonical: dict[str, str] = {}
    name_to_canonical: dict[str, str] = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cn = row["canonical_name"].strip()
            aro = row["aro_accession"].strip()
            if aro:
                aro_to_canonical[aro] = cn
            name_to_canonical[cn.lower()] = cn
    return aro_to_canonical, name_to_canonical


# ---------------------------------------------------------------------------
# OBO parser
# ---------------------------------------------------------------------------

def parse_obo(path: Path) -> dict[str, dict]:
    """
    Parse OBO file into a dict: ARO ID → term dict with keys:
      name        : str
      is_a        : list[str]   (ARO IDs)
      relationships: dict[str, list[str]]  (rel_type → list of ARO IDs)
      obsolete    : bool
    """
    terms: dict[str, dict] = {}
    current: dict | None = None

    with open(path) as f:
        for raw in f:
            line = raw.rstrip()

            if line == "[Term]":
                current = {"is_a": [], "relationships": {}, "obsolete": False}
                continue

            if line == "" or line == "[Typedef]":
                if current is not None and "id" in current:
                    terms[current["id"]] = current
                current = None
                continue

            if current is None:
                continue

            if line.startswith("id: "):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("is_a: "):
                m = re.match(r"is_a: (ARO:\d+)", line)
                if m:
                    current["is_a"].append(m.group(1))
            elif line.startswith("relationship: "):
                m = re.match(r"relationship: (\S+) (ARO:\d+)", line)
                if m:
                    rel, target = m.group(1), m.group(2)
                    current["relationships"].setdefault(rel, []).append(target)
            elif line == "is_obsolete: true":
                current["obsolete"] = True

    # Flush last term if file doesn't end with blank line
    if current is not None and "id" in current:
        terms[current["id"]] = current

    return terms


# ---------------------------------------------------------------------------
# Resolve canonical class via is_a ancestry
# ---------------------------------------------------------------------------

def resolve_class(
    aro_id: str,
    terms: dict,
    aro_to_canonical: dict[str, str],
    name_to_canonical: dict[str, str],
    _depth: int = 0,
    _visited: frozenset = frozenset(),
) -> str | None:
    """
    Walk up the is_a hierarchy from aro_id to find the nearest ancestor
    that maps to a known canonical class.  Returns the canonical name or None.
    """
    if _depth > 8 or aro_id in _visited:
        return None

    if aro_id in aro_to_canonical:
        return aro_to_canonical[aro_id]

    term = terms.get(aro_id)
    if term is None:
        return None

    if term.get("name", "").lower() in name_to_canonical:
        return name_to_canonical[term["name"].lower()]

    visited = _visited | {aro_id}
    for parent in term.get("is_a", []):
        result = resolve_class(parent, terms, aro_to_canonical,
                               name_to_canonical, _depth + 1, visited)
        if result:
            return result

    return None


# ---------------------------------------------------------------------------
# Extract drug → drug class membership
# ---------------------------------------------------------------------------

def extract_drug_class_members(
    terms: dict,
    aro_to_canonical: dict[str, str],
    name_to_canonical: dict[str, str],
) -> list[dict]:
    """
    For each OBO term that is NOT itself a known drug class, look for a
    canonical class in its is_a ancestry.  Emit one row per
    (drug_name, canonical_class) pair.
    """
    class_aro_ids = set(aro_to_canonical)
    class_names_lower = set(name_to_canonical)
    # Reverse lookup: canonical_class → its ARO accession
    canonical_to_aro = {v: k for k, v in aro_to_canonical.items()}

    rows = []
    for aro_id, term in terms.items():
        if term.get("obsolete") or not term.get("is_a"):
            continue
        name = term.get("name", "").strip()
        if not name:
            continue
        name_lower = name.lower()

        # Skip terms that are themselves drug classes
        if aro_id in class_aro_ids or name_lower in class_names_lower:
            continue

        # Find nearest canonical class in ancestry
        canonical_class = None
        for parent_id in term["is_a"]:
            canonical_class = resolve_class(
                parent_id, terms, aro_to_canonical, name_to_canonical
            )
            if canonical_class:
                break

        if canonical_class:
            rows.append({
                "canonical_drug":      name_lower,
                "drug_aro_accession":  aro_id,
                "canonical_class":     canonical_class,
                "class_aro_accession": canonical_to_aro.get(canonical_class, ""),
                "source":              "aro_obo",
            })

    return sorted(rows, key=lambda x: (x["canonical_class"], x["canonical_drug"]))


# ---------------------------------------------------------------------------
# Extract gene → drug class (confers_resistance_to_drug_class)
# ---------------------------------------------------------------------------

def extract_gene_class_links(
    terms: dict,
    aro_to_canonical: dict[str, str],
    name_to_canonical: dict[str, str],
) -> list[dict]:
    """
    For each term with confers_resistance_to_drug_class relationships,
    map the target ARO ID to a canonical class and emit a row.
    """
    rows = []
    for aro_id, term in terms.items():
        if term.get("obsolete"):
            continue
        name = term.get("name", "").strip()
        if not name:
            continue

        targets = term.get("relationships", {}).get(
            "confers_resistance_to_drug_class", []
        )
        if not targets:
            continue

        for class_aro in targets:
            canonical_class = aro_to_canonical.get(class_aro) or resolve_class(
                class_aro, terms, aro_to_canonical, name_to_canonical
            )
            if canonical_class:
                rows.append({
                    "aro_accession":      aro_id,
                    "gene_name":          name,
                    "canonical_class":    canonical_class,
                    "class_aro_accession": class_aro,
                    "source":             "aro_obo",
                })

    return sorted(rows, key=lambda x: (x["gene_name"], x["canonical_class"]))


# ---------------------------------------------------------------------------
# Write TSV helper
# ---------------------------------------------------------------------------

def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows):>5} rows → {path.relative_to(ROOT)}")


# ---------------------------------------------------------------------------
# Public entry point (called from harmonise.py)
# ---------------------------------------------------------------------------

def parse(obo_path: Path = OBO_PATH) -> tuple[list[dict], list[dict]]:
    """
    Returns:
      drug_class_members : list of drug→class dicts (for drug_class_member.tsv)
      gene_class_links   : list of gene→class dicts (for aro_gene_class.tsv)
    """
    class_mapping_path = ROOT / "harmonise/class_mapping.tsv"
    aro_to_canonical, name_to_canonical = load_canonical_classes(class_mapping_path)

    print(f"Parsing ARO OBO ({obo_path.name}) …")
    terms = parse_obo(obo_path)
    print(f"  Loaded {len(terms)} terms")

    drug_members = extract_drug_class_members(terms, aro_to_canonical, name_to_canonical)
    gene_links   = extract_gene_class_links(terms, aro_to_canonical, name_to_canonical)

    return drug_members, gene_links


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    out = ROOT / "harmonise"

    drug_members, gene_links = parse()

    write_tsv(
        out / "aro_drug_class_member.tsv",
        drug_members,
        ["canonical_drug", "drug_aro_accession", "canonical_class",
         "class_aro_accession", "source"],
    )
    write_tsv(
        out / "aro_gene_class.tsv",
        gene_links,
        ["aro_accession", "gene_name", "canonical_class",
         "class_aro_accession", "source"],
    )

    # Stats
    unique_drugs   = len({r["canonical_drug"]   for r in drug_members})
    unique_classes = len({r["canonical_class"]  for r in drug_members})
    unique_genes   = len({r["aro_accession"]    for r in gene_links})
    print(f"\nSummary:")
    print(f"  Drug→class pairs : {len(drug_members)}  ({unique_drugs} drugs, {unique_classes} classes)")
    print(f"  Gene→class pairs : {len(gene_links)}  ({unique_genes} genes)")
