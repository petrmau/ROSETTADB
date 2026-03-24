#!/usr/bin/env python3
"""
enrich_atc_groups.py
====================
Fill missing atc_group1 / atc_group2 values in drug_canonical.tsv by
deriving them directly from the ATC code prefix, using a lookup table built
from the AMR-R antimicrobials.txt reference.

ATC hierarchy used:
  atc_group1  →  4-char prefix  (e.g. J01D → "Other beta-lactam antibacterials")
  atc_group2  →  5-char prefix  (e.g. J01DB → "First-generation cephalosporins")

The lookup is learned from rows in antimicrobials.txt that already have both
an ATC code and group labels, so no external API calls or internet access are
required.

A drug with multiple ATC codes (comma-separated) uses the first code whose
prefix appears in the lookup.  J-prefixed codes are tried before others,
mirroring the ATC priority convention in enrich_amr_r.py.

Usage:
    python harmonise/enrich_atc_groups.py [--dry-run]

    --dry-run   Print what would change without writing any files.
"""

import argparse
import csv
from pathlib import Path

ROOT      = Path(__file__).parent.parent
DRUG_TSV  = ROOT / "harmonise/drug_canonical.tsv"
AMR_TSV   = ROOT / "harmonise/antimicrobials.txt"

# ATC priority order (J first, then Q, then others) — same as enrich_amr_r.py
ATC_PRIORITY = ["J", "Q", "P", "D", "A", "L", "B", "C", "G", "H", "M", "N", "R", "S", "V"]


def _atc_rank(code: str) -> int:
    letter = code[0].upper() if code else "Z"
    try:
        return ATC_PRIORITY.index(letter)
    except ValueError:
        return len(ATC_PRIORITY)


def _split_atcs(field: str) -> list[str]:
    """Return individual ATC codes from a comma-separated field, sorted by priority."""
    if not field or field.strip().upper() == "NA":
        return []
    codes = [c.strip() for c in field.split(",") if c.strip() and c.strip().upper() != "NA"]
    return sorted(codes, key=_atc_rank)


def build_prefix_maps(amr_tsv: Path) -> tuple[dict[str, str], dict[str, str]]:
    """
    Scan antimicrobials.txt and build:
      g1_map: 4-char ATC prefix → atc_group1 label
      g2_map: 5-char ATC prefix → atc_group2 label

    Only rows that already have both an ATC code and group labels contribute.
    First occurrence wins (INN rows appear before synonym rows in the file).
    """
    g1_map: dict[str, str] = {}
    g2_map: dict[str, str] = {}

    with open(amr_tsv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            atc_raw = row.get("atc", "").strip().strip('"')
            g1      = row.get("atc_group1", "").strip().strip('"')
            g2      = row.get("atc_group2", "").strip().strip('"')

            if not atc_raw or atc_raw.upper() == "NA" or not g1:
                continue

            for atc in _split_atcs(atc_raw):
                if len(atc) >= 4:
                    g1_map.setdefault(atc[:4], g1)
                if len(atc) >= 5 and g2:
                    g2_map.setdefault(atc[:5], g2)

    return g1_map, g2_map


def lookup_groups(
    atc_field: str,
    g1_map: dict[str, str],
    g2_map: dict[str, str],
) -> tuple[str, str]:
    """
    Given a comma-separated ATC field, return (group1, group2) by scanning
    codes in priority order and taking the first prefix hit for each level.
    Returns ("", "") when nothing matches.
    """
    codes = _split_atcs(atc_field)
    g1 = next((g1_map[c[:4]] for c in codes if len(c) >= 4 and c[:4] in g1_map), "")
    g2 = next((g2_map[c[:5]] for c in codes if len(c) >= 5 and c[:5] in g2_map), "")
    return g1, g2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would change without writing files.")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Build lookup from AMR-R reference
    # ------------------------------------------------------------------
    print(f"Building ATC prefix maps from {AMR_TSV.relative_to(ROOT)} …")
    g1_map, g2_map = build_prefix_maps(AMR_TSV)
    print(f"  {len(g1_map)} group1 prefixes (4-char), "
          f"{len(g2_map)} group2 prefixes (5-char)")

    # ------------------------------------------------------------------
    # Load drug_canonical.tsv
    # ------------------------------------------------------------------
    with open(DRUG_TSV, newline="", encoding="utf-8") as f:
        drugs = list(csv.DictReader(f, delimiter="\t"))

    if not drugs:
        print("drug_canonical.tsv is empty — nothing to do.")
        return

    fieldnames = list(drugs[0].keys())
    for col in ("atc_group1", "atc_group2"):
        if col not in fieldnames:
            fieldnames.append(col)
            for row in drugs:
                row.setdefault(col, "")

    # ------------------------------------------------------------------
    # Fill missing groups
    # ------------------------------------------------------------------
    filled_g1 = filled_g2 = skipped_no_atc = skipped_already = 0
    g2_still_missing: list[str] = []

    for row in drugs:
        atc = row.get("atc_code", "").strip()
        g1  = row.get("atc_group1", "").strip()
        g2  = row.get("atc_group2", "").strip()

        if not atc:
            skipped_no_atc += 1
            continue

        if g1 and g2:
            skipped_already += 1
            continue

        new_g1, new_g2 = lookup_groups(atc, g1_map, g2_map)

        if not g1 and new_g1:
            if args.dry_run:
                print(f"  [dry] {row['canonical_name']}: atc_group1 ← {new_g1!r}")
            else:
                row["atc_group1"] = new_g1
            filled_g1 += 1

        if not g2 and new_g2:
            if args.dry_run:
                print(f"  [dry] {row['canonical_name']}: atc_group2 ← {new_g2!r}")
            else:
                row["atc_group2"] = new_g2
            filled_g2 += 1

        if not g2 and not new_g2:
            g2_still_missing.append(f"{row['canonical_name']} ({atc})")

    # ------------------------------------------------------------------
    # Write (unless dry-run)
    # ------------------------------------------------------------------
    print(f"\nResults:")
    print(f"  Skipped (no atc_code)   : {skipped_no_atc}")
    print(f"  Skipped (already set)   : {skipped_already}")
    print(f"  Filled  atc_group1      : {filled_g1}")
    print(f"  Filled  atc_group2      : {filled_g2}")

    if g2_still_missing:
        print(f"  atc_group2 still empty  : {len(g2_still_missing)}")
        for entry in g2_still_missing:
            print(f"    {entry}")

    if args.dry_run:
        print("\n[dry-run] No files written.")
        return

    with open(DRUG_TSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(drugs)

    print(f"\nWritten → {DRUG_TSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
