#!/usr/bin/env python3
"""
enrich_atc_groups.py
====================
Two-pass offline enrichment of drug_canonical.tsv:

Pass 1 — ATC code fill (name lookup)
  For drugs that have no atc_code, look up the canonical_name
  (case-insensitive) in the local WHO ATC/ATCvet code tables
  (atc_codes_all.tsv and, if present, atcvet_codes_all.tsv).
  J-category codes are preferred over other categories when a name
  appears in multiple trees (same priority as enrich_amr_r.py /
  enrich.py).  Only fills; never overwrites an existing atc_code.

Pass 2 — ATC group fill (prefix lookup)
  For drugs that now have an atc_code but still lack atc_group1 /
  atc_group2, derives the group labels from the ATC code prefix using
  a lookup table built from the AMR-R antimicrobials.txt rows that
  already carry group labels:
    4-char prefix (e.g. J01D)  → atc_group1
    5-char prefix (e.g. J01DB) → atc_group2
  No network calls; no external dependencies beyond the local TSV files.

Usage:
    python harmonise/enrich_atc_groups.py [--dry-run]

    --dry-run   Print what would change without writing any files.
"""

import argparse
import csv
from pathlib import Path

ROOT         = Path(__file__).parent.parent
DRUG_TSV     = ROOT / "harmonise/drug_canonical.tsv"
AMR_TSV      = ROOT / "harmonise/antimicrobials.txt"
ATC_TSV      = ROOT / "harmonise/atc_codes_all.tsv"
ATCVET_TSV   = ROOT / "harmonise/atcvet_codes_all.tsv"

# ATC priority order — J first, Q (ATCvet) immediately after, then others.
# Same convention used by enrich_amr_r.py and enrich.py.
ATC_PRIORITY = ["J", "Q", "P", "D", "A", "L", "B", "C", "G", "H", "M", "N", "R", "S", "V"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _atc_rank(code: str) -> int:
    letter = code[0].upper() if code else "Z"
    try:
        return ATC_PRIORITY.index(letter)
    except ValueError:
        return len(ATC_PRIORITY)


def _is_pharmacological(code: str) -> bool:
    """Return True only for codes whose hierarchy reliably encodes pharmacological class.

    The WHO ATC system uses J (human) and QJ (ATCvet) for systemic
    antimicrobials; in both, each level of the code encodes a progressively
    finer pharmacological grouping.

    Route/indication-based ATCvet categories such as QG51 (intramammary
    preparations) or QS (sensory organs) mix drugs from completely different
    pharmacological classes under the same prefix — using their prefixes for
    group-label lookup produces wrong assignments (e.g. amoxicillin's QG51AA03
    would propagate "Beta-lactam antibacterials, penicillins" to cefquinome's
    QG51AA07 and gentamicin's QG51AA04).

    Restricting to J / QJ codes avoids all such cross-contamination.
    """
    u = code.upper()
    return u.startswith("J") or u.startswith("QJ")


def _atc_prefix_lengths(code: str) -> tuple[int, int]:
    """Return (group1_prefix_len, group2_prefix_len) for an ATC code.

    ATCvet codes start with 'Q', which shifts every level boundary by one
    character relative to human ATC codes:

      Human   J01CA04  → level-3 prefix = 4 chars (J01C)
                          level-4 prefix = 5 chars (J01CA)
      ATCvet  QJ01CA04 → level-3 prefix = 5 chars (QJ01C)
                          level-4 prefix = 6 chars (QJ01CA)

    Using a flat 4/5-char cut for ATCvet codes produces level-2 and
    level-3 prefixes respectively — one level too shallow — causing broad
    class labels (e.g. "Penicillins with extended spectrum") to be
    incorrectly assigned to unrelated subclasses (e.g. fourth-generation
    cephalosporins).
    """
    if code.upper().startswith("Q"):
        return 5, 6
    return 4, 5



def _split_atcs(field: str) -> list[str]:
    """Return individual ATC codes from a comma-separated field, sorted by priority."""
    if not field or field.strip().upper() == "NA":
        return []
    codes = [c.strip() for c in field.split(",") if c.strip() and c.strip().upper() != "NA"]
    return sorted(codes, key=_atc_rank)


# ---------------------------------------------------------------------------
# Pass 1 — build name → ATC code lookup from local WHO tables
# ---------------------------------------------------------------------------

def build_name_atc_map(*tsv_paths: Path) -> dict[str, str]:
    """
    Build a lowercase-name → best ATC code mapping from one or more ATC TSVs.

    When the same name appears in multiple categories the highest-priority
    category (J > Q > P > …) wins.  ATCvet (Q-prefix) codes from
    atcvet_codes_all.tsv are included automatically if the file exists.
    """
    # name → list of codes (we pick best after collecting all)
    name_codes: dict[str, list[str]] = {}

    for path in tsv_paths:
        if not path.exists():
            continue
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                code = row.get("atc_code", "").strip()
                name = row.get("name", "").strip().lower()
                if not code or not name:
                    continue
                name_codes.setdefault(name, []).append(code)

    # For each name keep the highest-priority code
    return {
        name: sorted(codes, key=_atc_rank)[0]
        for name, codes in name_codes.items()
    }


# ---------------------------------------------------------------------------
# Pass 2 — build ATC prefix → group label lookup from AMR-R reference
# ---------------------------------------------------------------------------

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
                if not _is_pharmacological(atc):
                    continue   # skip route-based categories (QG51, QS, …)
                l1, l2 = _atc_prefix_lengths(atc)
                if len(atc) >= l1:
                    g1_map.setdefault(atc[:l1], g1)
                if len(atc) >= l2 and g2:
                    g2_map.setdefault(atc[:l2], g2)

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
    codes = [c for c in _split_atcs(atc_field) if _is_pharmacological(c)]
    g1 = next(
        (g1_map[c[:_atc_prefix_lengths(c)[0]]]
         for c in codes
         if len(c) >= _atc_prefix_lengths(c)[0]
         and c[:_atc_prefix_lengths(c)[0]] in g1_map),
        "",
    )
    g2 = next(
        (g2_map[c[:_atc_prefix_lengths(c)[1]]]
         for c in codes
         if len(c) >= _atc_prefix_lengths(c)[1]
         and c[:_atc_prefix_lengths(c)[1]] in g2_map),
        "",
    )
    return g1, g2


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would change without writing files.")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Pass 1 setup — name → ATC code lookup
    # ------------------------------------------------------------------
    atc_sources = [p for p in (ATC_TSV, ATCVET_TSV) if p.exists()]
    print("Building name→ATC map from: %s" %
          ", ".join(p.relative_to(ROOT).as_posix() for p in atc_sources))
    name_atc_map = build_name_atc_map(*atc_sources)
    print(f"  {len(name_atc_map)} unique drug names indexed")

    # ------------------------------------------------------------------
    # Pass 2 setup — prefix → group label lookup
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
    for col in ("atc_code", "atc_group1", "atc_group2"):
        if col not in fieldnames:
            fieldnames.append(col)
            for row in drugs:
                row.setdefault(col, "")

    # ------------------------------------------------------------------
    # Pass 1 — fill missing atc_code from name lookup
    # ------------------------------------------------------------------
    print("\n--- Pass 1: fill atc_code from local ATC name lookup ---")
    filled_atc = 0

    for row in drugs:
        if row.get("atc_code", "").strip():
            continue  # already has a code
        name = row["canonical_name"].strip().lower()
        code = name_atc_map.get(name)
        if code:
            if args.dry_run:
                print(f"  [dry] {row['canonical_name']}: atc_code ← {code!r}")
            else:
                row["atc_code"] = code
            filled_atc += 1

    # ------------------------------------------------------------------
    # Pass 2 — fill missing atc_group1 / atc_group2 from prefix
    # ------------------------------------------------------------------
    print("\n--- Pass 2: fill atc_group1/2 from ATC code prefix ---")
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
    # Summary
    # ------------------------------------------------------------------
    print(f"\nResults:")
    print(f"  Pass 1 — filled atc_code    : {filled_atc}")
    print(f"  Pass 2 — skipped (no atc)   : {skipped_no_atc}")
    print(f"  Pass 2 — skipped (complete) : {skipped_already}")
    print(f"  Pass 2 — filled atc_group1  : {filled_g1}")
    print(f"  Pass 2 — filled atc_group2  : {filled_g2}")

    if g2_still_missing:
        print(f"  atc_group2 still empty      : {len(g2_still_missing)}")
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
