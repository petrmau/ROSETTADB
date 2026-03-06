#!/usr/bin/env python3
"""
parse_card_aro.py
=================
Parse CARD aro_index.tsv to build a gene → canonical drug class table,
cross-referencing canonical class names from class_mapping.tsv.

Output: harmonise/card_gene_class.tsv
  aro_accession      – ARO accession of the AMR gene model
  gene_name          – ARO Name of the gene/model
  card_short_name    – CARD Short Name (15-char abbreviation)
  gene_family        – AMR Gene Family
  resistance_mechanism
  model_type         – protein homolog | variant | rRNA | …  (from aro_index CARD Short Name column absent; use card.json)
  drug_class_raw     – raw Drug Class string from aro_index.tsv
  canonical_class    – mapped canonical class name (FK → amr.drug_class)

One row per (aro_accession, canonical_class) pair.  Genes conferring
resistance to multiple classes produce multiple rows (not slash-combined).
"""

import csv
import json
from pathlib import Path

ROOT = Path(__file__).parent.parent

# ── CARD-specific class name → canonical class alias table ───────────────────
# CARD uses some names that differ from our canonical (ARO-based) names.
# Keys are lowercase CARD drug class strings; values are canonical names
# (must exist in class_mapping.tsv) or None to skip.
CARD_CLASS_ALIASES: dict[str, str | None] = {
    # Spelling / phrasing variants
    "penicillin beta-lactam":               "penicillin antibiotic",
    "sulfone antibiotic":                   "sulphone antibiotic",
    "isoniazid-like antibiotic":            "isonicotinic acid hydrazide",
    # Subclass → parent class (acceptable lossy mapping)
    "glycylcycline":                        "glycylcycline",   # kept as own class (see class_mapping)
    # CARD broad category → nearest canonical
    "disinfecting agents and antiseptics":  "disinfecting agents and antiseptics",
    # Skip – no meaningful canonical mapping
    "antibiotic without defined classification": None,
}


def load_canonical_classes() -> set[str]:
    path = ROOT / "harmonise/class_mapping.tsv"
    with open(path) as f:
        return {r["canonical_name"].strip()
                for r in csv.DictReader(f, delimiter="\t")
                if r.get("canonical_name", "").strip()}


def map_class(raw: str, canonical_set: set[str]) -> str | None:
    """Map a CARD drug class string to a canonical class name, or None to skip."""
    raw = raw.strip()
    lower = raw.lower()

    # Check alias table first
    if lower in CARD_CLASS_ALIASES:
        return CARD_CLASS_ALIASES[lower]

    # Direct match (case-insensitive)
    if lower in {c.lower() for c in canonical_set}:
        # Return the properly-cased canonical name
        for c in canonical_set:
            if c.lower() == lower:
                return c

    return None   # unmappable → skip


def parse(out_path: Path):
    canonical_set = load_canonical_classes()

    aro_index = ROOT / "sources/CARD/aro_index.tsv"
    rows_out = []
    unmapped: dict[str, int] = {}

    with open(aro_index) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            aro_acc   = row["ARO Accession"].strip()
            gene_name = row["ARO Name"].strip()
            short     = row["CARD Short Name"].strip()
            family    = row["AMR Gene Family"].strip()
            mechanism = row["Resistance Mechanism"].strip()
            drug_class_raw = row["Drug Class"].strip()

            # Drug Class is semicolon-delimited
            for cls_raw in drug_class_raw.split(";"):
                cls_raw = cls_raw.strip()
                if not cls_raw:
                    continue

                canonical = map_class(cls_raw, canonical_set)
                if canonical is None:
                    unmapped[cls_raw] = unmapped.get(cls_raw, 0) + 1
                    continue

                rows_out.append({
                    "aro_accession":        aro_acc,
                    "gene_name":            gene_name,
                    "card_short_name":      short,
                    "gene_family":          family,
                    "resistance_mechanism": mechanism,
                    "drug_class_raw":       cls_raw,
                    "canonical_class":      canonical,
                })

    # Write TSV
    fieldnames = [
        "aro_accession", "gene_name", "card_short_name",
        "gene_family", "resistance_mechanism",
        "drug_class_raw", "canonical_class",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"  Wrote {len(rows_out):>5} rows → {out_path.relative_to(ROOT)}")

    # Report unmapped classes
    if unmapped:
        skipped_rows = sum(unmapped.values())
        print(f"  Skipped {skipped_rows} rows across {len(unmapped)} unmapped CARD classes:")
        for cls, n in sorted(unmapped.items(), key=lambda x: -x[1]):
            print(f"    {n:>4}  {cls}")

    # Summary
    unique_genes   = len({r["aro_accession"] for r in rows_out})
    unique_classes = len({r["canonical_class"] for r in rows_out})
    print(f"  Unique genes: {unique_genes}, canonical classes covered: {unique_classes}")
    return rows_out


if __name__ == "__main__":
    out = ROOT / "harmonise/card_gene_class.tsv"
    print("Parsing CARD ARO index …")
    parse(out)
