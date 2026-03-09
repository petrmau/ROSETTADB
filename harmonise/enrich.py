#!/usr/bin/env python3
"""
enrich.py
=========
Enrich harmonise/drug_canonical.tsv with cross-database identifiers.

Strategy per drug (ChEBI-first):
  1. Search ChEBI by INN name → chebi_id, inchikey
     Pick the highest-star exact-name match (3★ = fully curated).
  2a. If ChEBI hit: search PubChem by InChIKey → pubchem_cid, atc_code
      (InChIKey lookup is unambiguous; avoids name-matching noise)
  2b. If no ChEBI hit: search PubChem by name → pubchem_cid, inchikey, atc_code
  3.  If PubChem hit but still no chebi_id: try PubChem xrefs as last resort

Skipped automatically:
  - Combination drugs (is_combination = True) — no single molecule CID
  - Class-level placeholders (context = resistance_mechanism)

Cache: harmonise/.enrich_cache.json
  Keys are drug canonical names; values store all four identifiers.
  Entries already in the cache are never re-fetched (run again to refresh
  a specific drug by deleting its key from the cache file).

Rate limits (no API key required):
  ChEBI  : ~3 req/s  → CHEBI_SLEEP = 0.35 s
  PubChem: ~5 req/s  → PC_SLEEP    = 0.22 s
enrich.py
=========
Enrich harmonise/drug_canonical.tsv with cross-database identifiers.

Strategy per drug (ChEBI-first):
  1. Search ChEBI by INN name → chebi_id, inchikey
     Pick the highest-star exact-name match (3★ = fully curated).
  2a. If ChEBI hit: search PubChem by InChIKey → pubchem_cid
      (InChIKey lookup is unambiguous; avoids name/salt ambiguity)
  2b. If no ChEBI hit: search PubChem by name → pubchem_cid, inchikey
  3.  ATC code: local WHO ATC table (harmonise/atc_codes_all.tsv) first;
      PubChem classification endpoint used only as fallback for misses.
  4.  If PubChem hit but still no chebi_id: try PubChem xrefs as last resort

Skipped automatically:
  - Combination drugs (is_combination = True) — no single molecule CID
  - Class-level placeholders (context = resistance_mechanism)

ATC disambiguation (when a name appears in multiple categories):
  J (antiinfectives) wins; otherwise earliest in ATC_PRIORITY is chosen.

Cache: harmonise/.enrich_cache.json
  Keys are drug canonical names; values store all four identifiers.
  Entries already in the cache are never re-fetched (run again to refresh
  a specific drug by deleting its key from the cache file).
  On each run, cached entries with an empty atc_code are back-filled from
  the local ATC table at zero network cost.

Rate limits (no API key required):
  ChEBI  : ~3 req/s  → CHEBI_SLEEP = 0.35 s
  PubChem: ~5 req/s  → PC_SLEEP    = 0.22 s
"""

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT           = Path(__file__).parent.parent
CACHE_PATH     = ROOT / "harmonise/.enrich_cache.json"
DRUG_TSV       = ROOT / "harmonise/drug_canonical.tsv"
ATC_TABLE_PATH = ROOT / "harmonise/atc_codes_all.tsv"

CHEBI_SLEEP = 0.35
PC_SLEEP    = 0.22
MAX_RETRIES = 4

# ATC category priority for disambiguation (most AMR-relevant first).
# J = antiinfectives, P = antiparasitics, D = dermatologicals (topical antifungals),
# A = alimentary (some antifungals/antiprotozoals), then remaining categories.
ATC_PRIORITY = ["J", "P", "D", "A", "L", "B", "C", "G", "H", "M", "N", "R", "S", "V"]

CHEBI_BASE = "https://www.ebi.ac.uk/webservices/chebi/2.0/api"
PC_BASE    = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# XML namespace used in all ChEBI web-service responses
CHEBI_NS = "https://www.ebi.ac.uk/webservices/chebi/2.0/api"


# ---------------------------------------------------------------------------
# Generic HTTP helper
# ---------------------------------------------------------------------------

def _get(url: str, retries: int = MAX_RETRIES) -> bytes | None:
    """Fetch URL with exponential-backoff retry. Returns None on failure."""
    delay = 2
    for _ in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code in (429, 503) or e.code >= 500:
                time.sleep(delay)
                delay *= 2
                continue
            return None
        except Exception:
            time.sleep(delay)
            delay *= 2
    return None


def _get_json(url: str) -> dict | None:
    raw = _get(url)
    if raw is None:
        return None
    try:
        return json.loads(raw)
    except Exception:
        return None


def _get_xml(url: str) -> ET.Element | None:
    raw = _get(url)
    if raw is None:
        return None
    try:
        return ET.fromstring(raw)
    except Exception:
        return None


def _xml_text(root: ET.Element, tag: str) -> str:
    """Find first element matching tag (wildcard namespace) and return its text."""
    el = root.find(f".//{{{CHEBI_NS}}}{tag}")
    if el is None:
        # fallback: no-namespace (some responses omit it)
        el = root.find(f".//{tag}")
    return (el.text or "").strip() if el is not None else ""


def _xml_all(root: ET.Element, tag: str) -> list[ET.Element]:
    els = root.findall(f".//{{{CHEBI_NS}}}{tag}")
    if not els:
        els = root.findall(f".//{tag}")
    return els


# ---------------------------------------------------------------------------
# ChEBI lookups
# ---------------------------------------------------------------------------

def chebi_search_by_name(name: str) -> dict:
    """
    Search ChEBI for an INN drug name.
    Returns the best-matching entity as {chebi_id, chebi_name, star, inchikey}
    or an empty dict if nothing found.

    Selection criteria (in priority order):
      1. Exact case-insensitive name match among 3-star entries
      2. Exact case-insensitive name match among any star rating
      3. Highest-star entry among all results (≥2 stars preferred)
    """
    enc = urllib.parse.quote(name)
    url = (f"{CHEBI_BASE}/getLiteEntity"
           f"?search={enc}&searchCategory=ALL_NAMES"
           f"&maximumResults=10&starsCategory=ALL")
    root = _get_xml(url)
    if root is None:
        return {}

    candidates = []
    for el in _xml_all(root, "ListElement"):
        chebi_id   = _xml_text(el, "chebiId")
        chebi_name = _xml_text(el, "chebiAsciiName")
        star_txt   = _xml_text(el, "entityStar")
        star       = int(star_txt) if star_txt.isdigit() else 0
        if chebi_id:
            candidates.append({
                "chebi_id":   chebi_id,
                "chebi_name": chebi_name,
                "star":       star,
            })

    if not candidates:
        return {}

    name_lower = name.lower().strip()

    # Priority 1: exact name match, 3-star
    for c in candidates:
        if c["chebi_name"].lower() == name_lower and c["star"] == 3:
            return c

    # Priority 2: exact name match, any star
    for c in candidates:
        if c["chebi_name"].lower() == name_lower:
            return c

    # Priority 3: highest-star result (accept ≥2 stars only to avoid noise)
    best = max(candidates, key=lambda c: c["star"])
    if best["star"] >= 2:
        return best

    return {}


def chebi_get_inchikey(chebi_id: str) -> str:
    """
    Fetch the InChIKey for a known ChEBI ID via getCompleteEntity.
    Returns empty string if unavailable.
    """
    enc = urllib.parse.quote(chebi_id)
    url = f"{CHEBI_BASE}/getCompleteEntity?chebiId={enc}"
    root = _get_xml(url)
    if root is None:
        return ""
    return _xml_text(root, "inchikey")


# ---------------------------------------------------------------------------
# PubChem lookups
# ---------------------------------------------------------------------------

def pc_cid_by_inchikey(inchikey: str) -> int | None:
    """Unambiguous CID lookup via InChIKey."""
    url = f"{PC_BASE}/compound/inchikey/{inchikey}/cids/JSON"
    data = _get_json(url)
    if data is None:
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


def pc_cid_by_name(name: str) -> int | None:
    """CID lookup via compound name. Less precise — used as fallback."""
    enc = urllib.parse.quote(name)
    # Try exact match first, then word match
    for suffix in ["", "?name_type=word"]:
        url = f"{PC_BASE}/compound/name/{enc}/cids/JSON{suffix}"
        data = _get_json(url)
        if data:
            cids = data.get("IdentifierList", {}).get("CID", [])
            if cids:
                return cids[0]
    return None


def pc_inchikey(cid: int) -> str:
    """Fetch InChIKey from PubChem for a given CID."""
    url = f"{PC_BASE}/compound/cid/{cid}/property/InChIKey/JSON"
    data = _get_json(url)
    if data is None:
        return ""
    props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
    return props.get("InChIKey", "")


def pc_chebi_xref(cid: int) -> str:
    """Return ChEBI ID from PubChem xrefs for a CID, or empty string."""
    url = f"{PC_BASE}/compound/cid/{cid}/xrefs/RegistryID/JSON"
    data = _get_json(url)
    if data is None:
        return ""
    ids = (data.get("InformationList", {})
               .get("Information", [{}])[0]
               .get("RegistryID", []))
    for rid in ids:
        if rid.upper().startswith("CHEBI:"):
            return rid
    return ""


def pc_atc(cid: int) -> str:
    """Return first WHO ATC code from PubChem classification, or empty string."""
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
           f"/data/compound/{cid}/JSON?heading=ATC+Code")
    data = _get_json(url)
    if data is None:
        return ""
    try:
        for section in data.get("Record", {}).get("Section", []):
            for sub in section.get("Section", []):
                if "ATC" in sub.get("TOCHeading", ""):
                    for info in sub.get("Information", []):
                        for val in info.get("Value", {}).get("StringWithMarkup", []):
                            s = val.get("String", "")
                            if s and len(s) >= 5:
                                return s.split()[0]
    except Exception:
        pass
    return ""


# ---------------------------------------------------------------------------
# Local WHO ATC table
# ---------------------------------------------------------------------------

def load_atc_table() -> dict[str, list[str]]:
    """
    Load harmonise/atc_codes_all.tsv into a lowercase-name → [atc_code, …] dict.
    Returns empty dict if the file is missing (enrich.py still works without it).
    """
    if not ATC_TABLE_PATH.exists():
        return {}
    table: dict[str, list[str]] = {}
    with open(ATC_TABLE_PATH, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            code = row.get("atc_code", "").strip()
            name = row.get("name", "").strip().lower()
            if code and name:
                table.setdefault(name, []).append(code)
    return table


def atc_from_table(name: str, table: dict[str, list[str]]) -> str:
    """
    Look up ATC code from the local WHO ATC table.
    If the same name maps to multiple codes (different categories), J wins;
    otherwise the earliest category in ATC_PRIORITY is chosen.
    Returns empty string on miss.
    """
    hits = table.get(name.lower(), [])
    if not hits:
        return ""
    if len(hits) == 1:
        return hits[0]

    def _rank(code: str) -> int:
        letter = code[0].upper() if code else "Z"
        try:
            return ATC_PRIORITY.index(letter)
        except ValueError:
            return len(ATC_PRIORITY)

    return sorted(hits, key=_rank)[0]


# ---------------------------------------------------------------------------
# Per-drug enrichment (ChEBI-first)
# ---------------------------------------------------------------------------

def enrich_drug(name: str, cache: dict, atc_table: dict | None = None) -> dict:
    """
    Return enrichment dict {chebi_id, inchikey, pubchem_cid, atc_code}.
    Consults cache first; if missing, runs ChEBI → PubChem pipeline.
    ATC code: local table lookup first, pc_atc() only as fallback.
    """
    if name in cache:
        return cache[name]

    result = {"chebi_id": "", "inchikey": "", "pubchem_cid": "", "atc_code": ""}

    # ── Step 1: ChEBI by name ──────────────────────────────────────────────
    chebi_hit = chebi_search_by_name(name)
    time.sleep(CHEBI_SLEEP)

    if chebi_hit:
        result["chebi_id"] = chebi_hit["chebi_id"]
        # Fetch InChIKey from ChEBI
        ik = chebi_get_inchikey(chebi_hit["chebi_id"])
        time.sleep(CHEBI_SLEEP)
        result["inchikey"] = ik

    # ── Step 2: PubChem ───────────────────────────────────────────────────
    cid = None

    if result["inchikey"]:
        # 2a: InChIKey → CID (unambiguous)
        cid = pc_cid_by_inchikey(result["inchikey"])
        time.sleep(PC_SLEEP)
    else:
        # 2b: name → CID (fallback)
        cid = pc_cid_by_name(name)
        time.sleep(PC_SLEEP)
        if cid and not result["inchikey"]:
            result["inchikey"] = pc_inchikey(cid)
            time.sleep(PC_SLEEP)

    if cid:
        result["pubchem_cid"] = str(cid)
        # ATC: local WHO table first (no network); PubChem as fallback
        if atc_table is not None:
            result["atc_code"] = atc_from_table(name, atc_table)
        if not result["atc_code"]:
            result["atc_code"] = pc_atc(cid)
            time.sleep(PC_SLEEP)

    # ── Step 3: ChEBI xref from PubChem (last resort) ─────────────────────
    if cid and not result["chebi_id"]:
        result["chebi_id"] = pc_chebi_xref(cid)
        time.sleep(PC_SLEEP)

    cache[name] = result
    return result


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def load_cache() -> dict:
    if CACHE_PATH.exists():
        try:
            return json.loads(CACHE_PATH.read_text())
        except Exception:
            return {}
    return {}


def save_cache(cache: dict):
    CACHE_PATH.write_text(json.dumps(cache, indent=2))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

SKIP_CONTEXTS = {"resistance_mechanism"}


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Enrich drug_canonical.tsv via ChEBI → PubChem")
    ap.add_argument("--refresh-missing-chebi", action="store_true",
                    help="Remove cache entries that have a PubChem CID but no ChEBI ID, "
                         "forcing a fresh ChEBI-first lookup for those drugs")
    args = ap.parse_args()

    with open(DRUG_TSV) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    fieldnames = list(rows[0].keys()) if rows else []
    cache = load_cache()
    print(f"Cache: {len(cache)} entries loaded from {CACHE_PATH.name}")

    # Load local WHO ATC table and back-fill any cached entries missing ATC
    atc_table = load_atc_table()
    if atc_table:
        print(f"ATC table: {sum(len(v) for v in atc_table.values())} entries "
              f"from {ATC_TABLE_PATH.name}")
        upgraded = 0
        for n, v in cache.items():
            if not v.get("atc_code"):
                code = atc_from_table(n, atc_table)
                if code:
                    v["atc_code"] = code
                    upgraded += 1
        if upgraded:
            print(f"  Back-filled ATC for {upgraded} cached entries (no network calls)")
    else:
        print(f"  Warning: {ATC_TABLE_PATH.name} not found — ATC will use PubChem only")

    if args.refresh_missing_chebi:
        to_drop = [n for n, v in cache.items()
                   if v.get("pubchem_cid") and not v.get("chebi_id")]
        for n in to_drop:
            del cache[n]
        print(f"  Dropped {len(to_drop)} entries missing chebi_id → will re-fetch via ChEBI")

    # Apply cached values to all rows first
    for row in rows:
        name = row["canonical_name"]
        if name in cache:
            data = cache[name]
            row["chebi_id"]    = data.get("chebi_id", "")
            row["inchikey"]    = data.get("inchikey", "")
            row["pubchem_cid"] = data.get("pubchem_cid", "")
            row["atc_code"]    = data.get("atc_code", "")

    # Drugs still needing a live lookup
    to_enrich = [
        r for r in rows
        if r["canonical_name"] not in cache
        and r.get("is_combination") != "True"
        and r.get("context", "") not in SKIP_CONTEXTS
    ]

    print(f"Enriching {len(to_enrich)} drugs (ChEBI → PubChem) …\n")

    for i, row in enumerate(to_enrich):
        name = row["canonical_name"]
        print(f"  [{i+1}/{len(to_enrich)}] {name} …", end=" ", flush=True)

        data = enrich_drug(name, cache, atc_table)

        row["chebi_id"]    = data.get("chebi_id", "")
        row["inchikey"]    = data.get("inchikey", "")
        row["pubchem_cid"] = data.get("pubchem_cid", "")
        row["atc_code"]    = data.get("atc_code", "")

        parts = []
        if data.get("chebi_id"):    parts.append(f"ChEBI={data['chebi_id']}")
        if data.get("pubchem_cid"): parts.append(f"CID={data['pubchem_cid']}")
        if data.get("inchikey"):    parts.append(f"IK={data['inchikey'][:14]}…")
        if data.get("atc_code"):    parts.append(f"ATC={data['atc_code']}")
        print(" | ".join(parts) if parts else "not found")

        if (i + 1) % 10 == 0:
            save_cache(cache)

    save_cache(cache)

    # Write enriched TSV in-place
    with open(DRUG_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    total       = len(rows)
    found_chebi = sum(1 for r in rows if r.get("chebi_id"))
    found_ik    = sum(1 for r in rows if r.get("inchikey"))
    found_cid   = sum(1 for r in rows if r.get("pubchem_cid"))
    found_atc   = sum(1 for r in rows if r.get("atc_code"))

    print(f"\nEnrichment complete ({total} drugs total):")
    print(f"  ChEBI ID    : {found_chebi}/{total}")
    print(f"  InChIKey    : {found_ik}/{total}")
    print(f"  PubChem CID : {found_cid}/{total}")
    print(f"  ATC code    : {found_atc}/{total}")

    missing = [
        r["canonical_name"] for r in rows
        if not r.get("chebi_id") and not r.get("pubchem_cid")
        and r.get("is_combination") != "True"
        and r.get("context", "") not in SKIP_CONTEXTS
    ]
    if missing:
        print(f"\n  No identifiers found for {len(missing)} drugs:")
        for m in sorted(missing):
            print(f"    - {m}")


if __name__ == "__main__":
    main()
