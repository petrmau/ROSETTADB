#!/usr/bin/env python3
"""
enrich_pubchem.py
=================
For each drug in harmonise/drug_canonical.tsv, query PubChem PUG REST to
retrieve:
  - pubchem_cid   (integer canonical CID)
  - inchikey      (standard InChIKey, 27 chars)
  - chebi_id      (if available via PubChem xrefs)
  - atc_code      (first ATC code from PubChem classification tree)

Results are written back into drug_canonical.tsv in-place.
A cache file (harmonise/.pubchem_cache.json) avoids redundant API calls on
re-runs.

Rate-limit: PubChem allows ~5 req/s without API key; we use 0.22 s sleep.
"""

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).parent.parent
CACHE_PATH = ROOT / "harmonise/.pubchem_cache.json"
DRUG_TSV = ROOT / "harmonise/drug_canonical.tsv"

SLEEP = 0.22   # seconds between requests
MAX_RETRIES = 4


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------

def _get(url: str, retries: int = MAX_RETRIES) -> bytes | None:
    delay = 2
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=15) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None          # not found — stop retrying
            if e.code == 429 or e.code >= 500:
                time.sleep(delay)
                delay *= 2
                continue
            return None
        except Exception:
            time.sleep(delay)
            delay *= 2
    return None


def pubchem_get_json(url: str) -> dict | None:
    raw = _get(url)
    if raw is None:
        return None
    try:
        return json.loads(raw)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# PubChem lookups
# ---------------------------------------------------------------------------

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def lookup_cid(name: str) -> int | None:
    """Return the first canonical CID for a compound name."""
    enc = urllib.parse.quote(name)
    url = f"{BASE}/compound/name/{enc}/cids/JSON?name_type=word"
    data = pubchem_get_json(url)
    if data is None:
        # try without name_type
        url2 = f"{BASE}/compound/name/{enc}/cids/JSON"
        data = pubchem_get_json(url2)
    if data is None:
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


def lookup_props(cid: int) -> dict:
    """Return InChIKey and IUPAC name for a CID."""
    url = f"{BASE}/compound/cid/{cid}/property/InChIKey,IUPACName/JSON"
    data = pubchem_get_json(url)
    if data is None:
        return {}
    props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
    return {
        "inchikey": props.get("InChIKey", ""),
        "iupac_name": props.get("IUPACName", ""),
    }


def lookup_chebi(cid: int) -> str:
    """Return ChEBI ID (e.g. 'CHEBI:28971') from PubChem xrefs, or empty string."""
    url = f"{BASE}/compound/cid/{cid}/xrefs/RegistryID/JSON"
    data = pubchem_get_json(url)
    if data is None:
        return ""
    ids = data.get("InformationList", {}).get("Information", [{}])[0].get("RegistryID", [])
    for rid in ids:
        if rid.upper().startswith("CHEBI:"):
            return rid
    return ""


def lookup_atc(cid: int) -> str:
    """
    Return first ATC code from PubChem's pharmacology classification for CID.
    PubChem stores ATC codes in the compound classification tree under
    the WHO ATC source.  We fetch the /classification endpoint.
    Falls back to empty string if none found.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading=ATC+Code"
    data = pubchem_get_json(url)
    if data is None:
        return ""
    # Navigate the PUG View tree looking for ATC Code sections
    try:
        record = data.get("Record", {})
        for section in record.get("Section", []):
            for sub in section.get("Section", []):
                if "ATC" in sub.get("TOCHeading", ""):
                    for info in sub.get("Information", []):
                        for val in info.get("Value", {}).get("StringWithMarkup", []):
                            s = val.get("String", "")
                            if s and len(s) >= 5:
                                return s.split()[0]  # first token is the code
    except Exception:
        pass
    return ""


# ---------------------------------------------------------------------------
# Per-drug enrichment (with cache)
# ---------------------------------------------------------------------------

def enrich_drug(name: str, cache: dict) -> dict:
    """
    Return enrichment dict for drug name.
    Checks cache first; falls back to live PubChem queries.
    """
    if name in cache:
        return cache[name]

    result = {"pubchem_cid": "", "inchikey": "", "chebi_id": "", "atc_code": ""}

    cid = lookup_cid(name)
    time.sleep(SLEEP)

    if cid is None:
        cache[name] = result
        return result

    result["pubchem_cid"] = str(cid)

    props = lookup_props(cid)
    time.sleep(SLEEP)
    result["inchikey"] = props.get("inchikey", "")

    chebi = lookup_chebi(cid)
    time.sleep(SLEEP)
    result["chebi_id"] = chebi

    atc = lookup_atc(cid)
    time.sleep(SLEEP)
    result["atc_code"] = atc

    cache[name] = result
    return result


# ---------------------------------------------------------------------------
# Main
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


def main():
    # Load existing table
    with open(DRUG_TSV) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    fieldnames = list(rows[0].keys()) if rows else []

    cache = load_cache()
    print(f"Cache has {len(cache)} entries")

    # Identify drugs that still need enrichment
    to_enrich = [
        r for r in rows
        if r["canonical_name"] not in cache   # not in cache at all → needs lookup
        and r.get("is_combination") != "True" # skip combinations (no single CID)
        and r.get("context") not in ("resistance_mechanism",)
    ]

    # Apply cache entries to all rows that have cache data (even if empty)
    for row in rows:
        name = row["canonical_name"]
        if name in cache:
            data = cache[name]
            row["pubchem_cid"] = data.get("pubchem_cid", "")
            row["inchikey"]    = data.get("inchikey", "")
            row["chebi_id"]    = data.get("chebi_id", "")
            row["atc_code"]    = data.get("atc_code", "")

    print(f"Enriching {len(to_enrich)} drugs via PubChem …")

    for i, row in enumerate(to_enrich):
        name = row["canonical_name"]
        print(f"  [{i+1}/{len(to_enrich)}] {name} …", end=" ", flush=True)
        data = enrich_drug(name, cache)

        # Write back into the row
        row["pubchem_cid"] = data.get("pubchem_cid", "")
        row["inchikey"]    = data.get("inchikey", "")
        row["chebi_id"]    = data.get("chebi_id", "")
        row["atc_code"]    = data.get("atc_code", "")

        status = f"CID={data['pubchem_cid'] or 'not found'}"
        if data.get("inchikey"):
            status += f"  IK={data['inchikey']}"
        if data.get("atc_code"):
            status += f"  ATC={data['atc_code']}"
        print(status)

        # Save cache every 10 lookups
        if (i + 1) % 10 == 0:
            save_cache(cache)

    save_cache(cache)

    # Write enriched TSV
    with open(DRUG_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    found_cid    = sum(1 for r in rows if r.get("pubchem_cid"))
    found_inchi  = sum(1 for r in rows if r.get("inchikey"))
    found_chebi  = sum(1 for r in rows if r.get("chebi_id"))
    found_atc    = sum(1 for r in rows if r.get("atc_code"))
    total        = len(rows)

    print(f"\nEnrichment complete ({total} drugs total):")
    print(f"  PubChem CID : {found_cid}/{total}")
    print(f"  InChIKey    : {found_inchi}/{total}")
    print(f"  ChEBI ID    : {found_chebi}/{total}")
    print(f"  ATC code    : {found_atc}/{total}")

    missing = [r["canonical_name"] for r in rows
               if not r.get("pubchem_cid")
               and r.get("is_combination") != "True"
               and r.get("context") not in ("resistance_mechanism",)]
    if missing:
        print(f"\n  Not found on PubChem ({len(missing)}):")
        for m in sorted(missing):
            print(f"    - {m}")


if __name__ == "__main__":
    main()
