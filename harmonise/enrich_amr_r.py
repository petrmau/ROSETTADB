#!/usr/bin/env python3
"""
enrich_amr_r.py
===============
Pre-enrich drug_canonical.tsv and drug_alias.tsv using the msberends/AMR
R-package antimicrobials reference dataset.

Source:
  https://raw.githubusercontent.com/msberends/AMR/main/data-raw/datasets/antimicrobials.txt

For each canonical drug in drug_canonical.tsv the script:
  - Fills in pubchem_cid  from the `cid` column (if blank; not overwritten)
  - Updates atc_code      from the `atc` column — stores ALL codes, sorted:
                          J first, then Q immediately after J, then others.
                          Always written when the reference has data, so a
                          prior single-code value becomes the full list.
  - Updates loinc_codes   from the `loinc` column (always written when present)
  - Updates atc_group1    from the `atc_group1` column (ATC level-2 group,
                          e.g. "Aminoglycoside antibacterials")
  - Updates atc_group2    from the `atc_group2` column (ATC level-3 group,
                          e.g. "Other aminoglycosides")
  - Adds synonyms and abbreviations to drug_alias.tsv

Combination drugs (e.g. "amoxicillin+clavulanic acid") are matched against
the reference by normalising both separators (/, +, &) and comparing the
sorted frozenset of component names — order-independent.

Running this before enrich.py reduces the number of live ChEBI/PubChem API
calls needed.

Usage:
  python harmonise/enrich_amr_r.py [--download]

  --download   Re-download antimicrobials.txt even if a local copy exists.
               Without this flag, the cached local copy is used.
"""

import argparse
import csv
import re
import urllib.request
from pathlib import Path

ROOT       = Path(__file__).parent.parent
DRUG_TSV   = ROOT / "harmonise/drug_canonical.tsv"
ALIAS_TSV  = ROOT / "harmonise/drug_alias.tsv"
AMR_TSV    = ROOT / "harmonise/antimicrobials.txt"
AMR_URL    = ("https://raw.githubusercontent.com/msberends/AMR"
              "/main/data-raw/datasets/antimicrobials.txt")

# ATC category priority.
# Q (veterinary) comes immediately after J (human antiinfectives) because
# most AMR drugs have a matching QJ... code and we want them adjacent.
# All other categories follow standard ATC order.
ATC_PRIORITY = ["J", "Q", "P", "D", "A", "L", "B", "C", "G", "H", "M", "N", "R", "S", "V"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalise(s: str) -> str:
    """Lowercase, strip, collapse whitespace."""
    return re.sub(r"\s+", " ", s.strip().lower())


def sort_atcs(atc_field: str) -> str:
    """
    Return all ATC codes from a comma-separated field, sorted by ATC_PRIORITY
    (J first, then Q, then others).  Duplicate and NA values are dropped.
    Returns an empty string when nothing useful is present.
    """
    if not atc_field or atc_field.strip().upper() == "NA":
        return ""
    codes_raw = [c.strip() for c in atc_field.split(",")
                 if c.strip() and c.strip().upper() != "NA"]
    # deduplicate, preserving first occurrence before sorting
    seen: set[str] = set()
    codes: list[str] = []
    for c in codes_raw:
        if c not in seen:
            seen.add(c)
            codes.append(c)
    if not codes:
        return ""

    def _rank(code: str) -> int:
        letter = code[0].upper() if code else "Z"
        try:
            return ATC_PRIORITY.index(letter)
        except ValueError:
            return len(ATC_PRIORITY)

    return ",".join(sorted(codes, key=_rank))


def clean_loinc(loinc_field: str) -> str:
    """
    Normalise the loinc field: strip NA tokens and return a
    comma-separated string of unique codes, or empty string.
    """
    if not loinc_field or loinc_field.strip().upper() == "NA":
        return ""
    seen: set[str] = set()
    unique: list[str] = []
    for c in loinc_field.split(","):
        c = c.strip()
        if c and c.upper() != "NA" and c not in seen:
            seen.add(c)
            unique.append(c)
    return ",".join(unique)


def combo_key(name: str) -> frozenset[str]:
    """
    Split a combination drug name on any of /, +, & and return a frozenset
    of normalised component names.  Returns an empty frozenset for singles.
    """
    parts = re.split(r"[/+&]", name)
    return frozenset(normalise(p) for p in parts if p.strip())


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def download_amr_txt() -> None:
    print(f"Downloading {AMR_URL} …")
    with urllib.request.urlopen(AMR_URL, timeout=30) as r:
        content = r.read()
    AMR_TSV.write_bytes(content)
    print(f"  Saved → {AMR_TSV.relative_to(ROOT)}")


# ---------------------------------------------------------------------------
# Load and index the AMR reference
# ---------------------------------------------------------------------------

def load_amr_reference() -> list[dict]:
    rows = []
    with open(AMR_TSV, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
    return rows


def build_name_index(amr_rows: list[dict]) -> dict[str, dict]:
    """
    Build a lookup: normalised name/alias → amr_row.
    Indexes `name` first (highest priority), then `synonyms`, then
    `abbreviations`.  First match wins so the INN name takes priority.
    """
    index: dict[str, dict] = {}

    # Pass 1: canonical INN names
    for row in amr_rows:
        key = normalise(row.get("name", ""))
        if key:
            index.setdefault(key, row)

    # Pass 2: synonyms
    for row in amr_rows:
        raw = row.get("synonyms", "")
        if not raw or raw.strip().upper() == "NA":
            continue
        for alias in raw.split(","):
            key = normalise(alias)
            if key:
                index.setdefault(key, row)

    # Pass 3: abbreviations
    for row in amr_rows:
        raw = row.get("abbreviations", "")
        if not raw or raw.strip().upper() == "NA":
            continue
        for alias in raw.split(","):
            key = normalise(alias)
            if key:
                index.setdefault(key, row)

    return index


def build_combo_index(amr_rows: list[dict]) -> dict[frozenset, dict]:
    """
    Build a lookup: frozenset of normalised component names → amr_row,
    restricted to entries whose name contains a combination separator (/ + &).
    Used to match our canonical combo names (using + or &) against the
    reference (which uses /).
    """
    index: dict[frozenset, dict] = {}
    for row in amr_rows:
        name = row.get("name", "")
        if not re.search(r"[/+&]", name):
            continue
        key = combo_key(name)
        if len(key) >= 2:
            index.setdefault(key, row)
    return index


# ---------------------------------------------------------------------------
# Apply enrichment from a matched reference row
# ---------------------------------------------------------------------------

def _clean_group(raw: str) -> str:
    """Return stripped group string, or empty string if absent/NA."""
    s = raw.strip() if raw else ""
    return "" if s.upper() == "NA" else s


def apply_hit(row: dict, hit: dict, cn: str,
              new_aliases: list[tuple]) -> tuple[bool, bool, bool, bool, bool]:
    """
    Apply data from a matched AMR reference row to a drug_canonical row.
    Returns (cid_filled, atc_updated, loinc_updated, g1_updated, g2_updated).
    """
    cid_filled = atc_updated = loinc_updated = g1_updated = g2_updated = False

    # PubChem CID — only if currently blank (existing value is preserved)
    if not row.get("pubchem_cid"):
        cid = hit.get("cid", "").strip()
        if cid and cid.upper() != "NA":
            row["pubchem_cid"] = cid
            cid_filled = True

    # ATC codes — always update to get the full sorted list
    atc = sort_atcs(hit.get("atc", ""))
    if atc and atc != row.get("atc_code", ""):
        row["atc_code"] = atc
        atc_updated = True

    # LOINC codes — always update when present
    loinc = clean_loinc(hit.get("loinc", ""))
    if loinc and loinc != row.get("loinc_codes", ""):
        row["loinc_codes"] = loinc
        loinc_updated = True

    # ATC groups — always update when present
    g1 = _clean_group(hit.get("atc_group1", ""))
    if g1 and g1 != row.get("atc_group1", ""):
        row["atc_group1"] = g1
        g1_updated = True

    g2 = _clean_group(hit.get("atc_group2", ""))
    if g2 and g2 != row.get("atc_group2", ""):
        row["atc_group2"] = g2
        g2_updated = True

    # Aliases — only from single-drug entries (combo aliases are not useful)
    for field, alias_type in (
        ("name",          "source_name"),
        ("synonyms",      "synonym"),
        ("abbreviations", "abbreviation"),
    ):
        raw = hit.get(field, "")
        if not raw or raw.strip().upper() == "NA":
            continue
        values = raw.split(",") if field != "name" else [raw]
        for v in values:
            v = v.strip()
            if v and normalise(v) != normalise(cn):
                new_aliases.append((cn, v, alias_type, "amr_r"))

    return cid_filled, atc_updated, loinc_updated, g1_updated, g2_updated


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Pre-enrich drug_canonical.tsv and drug_alias.tsv "
                    "from the msberends/AMR reference dataset"
    )
    ap.add_argument("--download", action="store_true",
                    help="Re-download antimicrobials.txt even if a local copy exists")
    args = ap.parse_args()

    if args.download or not AMR_TSV.exists():
        download_amr_txt()
    else:
        print(f"Using cached {AMR_TSV.relative_to(ROOT)}  "
              f"(pass --download to refresh)")

    print("Loading AMR reference …")
    amr_rows    = load_amr_reference()
    name_index  = build_name_index(amr_rows)
    combo_index = build_combo_index(amr_rows)
    print(f"  {len(amr_rows)} entries, {len(name_index)} name/alias keys, "
          f"{len(combo_index)} combination keys")

    # ------------------------------------------------------------------
    # Enrich drug_canonical.tsv
    # ------------------------------------------------------------------
    with open(DRUG_TSV, newline="", encoding="utf-8") as f:
        drugs = list(csv.DictReader(f, delimiter="\t"))

    if not drugs:
        print("drug_canonical.tsv is empty — nothing to do.")
        return

    # Ensure new columns exist (safe to re-run against an older TSV)
    fieldnames = list(drugs[0].keys())
    for col in ("loinc_codes", "atc_group1", "atc_group2"):
        if col not in fieldnames:
            fieldnames.append(col)
            for row in drugs:
                row.setdefault(col, "")

    new_aliases: list[tuple] = []
    matched_single = matched_combo = 0
    filled_cid = 0
    updated_atc = updated_loinc = updated_g1 = updated_g2 = 0

    for row in drugs:
        cn  = row["canonical_name"]
        hit = name_index.get(normalise(cn))

        if hit is None and row.get("is_combination", "").strip() == "True":
            # Try combo matching by component set
            key = combo_key(cn)
            if len(key) >= 2:
                hit = combo_index.get(key)
            if hit is not None:
                matched_combo += 1
        elif hit is not None:
            matched_single += 1

        if hit is None:
            continue

        cid_f, atc_u, loinc_u, g1_u, g2_u = apply_hit(row, hit, cn, new_aliases)
        if cid_f:   filled_cid   += 1
        if atc_u:   updated_atc  += 1
        if loinc_u: updated_loinc += 1
        if g1_u:    updated_g1   += 1
        if g2_u:    updated_g2   += 1

    total_matched = matched_single + matched_combo
    print(f"\ndrug_canonical.tsv: {total_matched}/{len(drugs)} drugs matched "
          f"({matched_single} single, {matched_combo} combination)")
    print(f"  Filled  pubchem_cid : {filled_cid}")
    print(f"  Updated atc_code    : {updated_atc}")
    print(f"  Updated loinc_codes : {updated_loinc}")
    print(f"  Updated atc_group1  : {updated_g1}")
    print(f"  Updated atc_group2  : {updated_g2}")

    with open(DRUG_TSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(drugs)

    # ------------------------------------------------------------------
    # Enrich drug_alias.tsv
    # ------------------------------------------------------------------
    existing_aliases: set[tuple] = set()
    with open(ALIAS_TSV, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            existing_aliases.add((
                row["canonical_name"], row["alias"],
                row["alias_type"],     row["source"],
            ))

    new_unique = [a for a in new_aliases if a not in existing_aliases]
    all_aliases = sorted(existing_aliases | set(new_aliases))

    with open(ALIAS_TSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["canonical_name", "alias", "alias_type", "source"])
        writer.writerows(all_aliases)

    print(f"\ndrug_alias.tsv: added {len(new_unique)} new aliases (source=amr_r)")
    print(f"  Total aliases now : {len(all_aliases)}")


if __name__ == "__main__":
    main()
