#!/usr/bin/env python3
"""
enrich_amr_r.py
===============
Pre-enrich drug_canonical.tsv and drug_alias.tsv using the msberends/AMR
R-package antimicrobials reference dataset.

Source:
  https://raw.githubusercontent.com/msberends/AMR/main/data-raw/datasets/antimicrobials.txt

For each canonical drug in drug_canonical.tsv the script:
  - Fills in pubchem_cid  from the `cid`  column (if blank)
  - Fills in atc_code     from the `atc`  column (J-category wins; same
                          priority as enrich.py)
  - Fills in loinc_codes  from the `loinc` column (comma-separated LOINC
                          susceptibility test identifiers)
  - Adds synonyms and abbreviations to drug_alias.tsv

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

# Same ATC category priority as enrich.py
ATC_PRIORITY = ["J", "P", "D", "A", "L", "B", "C", "G", "H", "M", "N", "R", "S", "V"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalise(s: str) -> str:
    """Lowercase, strip, collapse whitespace."""
    return re.sub(r"\s+", " ", s.strip().lower())


def best_atc(atc_field: str) -> str:
    """
    Pick the best ATC code from a comma-separated list.
    J (antiinfectives) wins; otherwise follows ATC_PRIORITY order.
    Returns empty string if the field is absent or NA.
    """
    if not atc_field or atc_field.strip().upper() == "NA":
        return ""
    codes = [c.strip() for c in atc_field.split(",")
             if c.strip() and c.strip().upper() != "NA"]
    if not codes:
        return ""
    if len(codes) == 1:
        return codes[0]

    def _rank(code: str) -> int:
        letter = code[0].upper() if code else "Z"
        try:
            return ATC_PRIORITY.index(letter)
        except ValueError:
            return len(ATC_PRIORITY)

    return sorted(codes, key=_rank)[0]


def clean_loinc(loinc_field: str) -> str:
    """
    Normalise the loinc field: strip NA tokens and return a
    comma-separated string of unique codes, or empty string.
    """
    if not loinc_field or loinc_field.strip().upper() == "NA":
        return ""
    codes = [c.strip() for c in loinc_field.split(",")
             if c.strip() and c.strip().upper() != "NA"]
    # deduplicate while preserving order
    seen: set[str] = set()
    unique = []
    for c in codes:
        if c not in seen:
            seen.add(c)
            unique.append(c)
    return ",".join(unique)


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
    Indexes the `name` field first, then `synonyms`, then `abbreviations`
    (first match wins so the canonical INN name takes priority).
    """
    index: dict[str, dict] = {}

    # Pass 1: canonical INN names (highest priority)
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
    amr_rows = load_amr_reference()
    name_index = build_name_index(amr_rows)
    print(f"  {len(amr_rows)} entries indexed under {len(name_index)} name/alias keys")

    # ------------------------------------------------------------------
    # Enrich drug_canonical.tsv
    # ------------------------------------------------------------------
    with open(DRUG_TSV, newline="", encoding="utf-8") as f:
        drugs = list(csv.DictReader(f, delimiter="\t"))

    if not drugs:
        print("drug_canonical.tsv is empty — nothing to do.")
        return

    # Ensure loinc_codes column exists in the in-memory rows
    fieldnames = list(drugs[0].keys())
    if "loinc_codes" not in fieldnames:
        fieldnames.append("loinc_codes")
        for row in drugs:
            row.setdefault("loinc_codes", "")

    new_aliases: list[tuple] = []
    filled_cid = filled_atc = filled_loinc = 0
    matched = 0

    for row in drugs:
        cn = row["canonical_name"]
        hit = name_index.get(normalise(cn))
        if hit is None:
            continue
        matched += 1

        # PubChem CID
        if not row.get("pubchem_cid"):
            cid = hit.get("cid", "").strip()
            if cid and cid.upper() != "NA":
                row["pubchem_cid"] = cid
                filled_cid += 1

        # ATC code
        if not row.get("atc_code"):
            atc = best_atc(hit.get("atc", ""))
            if atc:
                row["atc_code"] = atc
                filled_atc += 1

        # LOINC codes
        if not row.get("loinc_codes"):
            loinc = clean_loinc(hit.get("loinc", ""))
            if loinc:
                row["loinc_codes"] = loinc
                filled_loinc += 1

        # Collect aliases: INN name, synonyms, abbreviations
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

    print(f"\ndrug_canonical.tsv: {matched}/{len(drugs)} drugs matched")
    print(f"  Filled pubchem_cid : {filled_cid}")
    print(f"  Filled atc_code    : {filled_atc}")
    print(f"  Filled loinc_codes : {filled_loinc}")

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

    print(f"\ndrug_alias.tsv: added {len(new_unique)} new aliases "
          f"(source=amr_r)")
    print(f"  Total aliases now : {len(all_aliases)}")


if __name__ == "__main__":
    main()
