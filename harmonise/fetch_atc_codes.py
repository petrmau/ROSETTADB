#!/usr/bin/env python3
"""
Fetch ATC codes and DDD information from https://atcddd.fhi.no/atc_ddd_index/
and save to a TSV file.

Usage:
    python fetch_atc_codes.py all                   # Fetch all 14 top-level categories
    python fetch_atc_codes.py A                     # Fetch a single top-level category
    python fetch_atc_codes.py A B C                 # Fetch multiple categories
    python fetch_atc_codes.py A10                   # Fetch a specific sub-level code
    python fetch_atc_codes.py A10BA02               # Fetch a specific leaf code
    python fetch_atc_codes.py all -o my_output.tsv  # Custom output file
    python fetch_atc_codes.py J --debug             # Debug mode: shows raw HTML parse info

Top-level categories:
    A  Alimentary Tract and Metabolism
    B  Blood and Blood Forming Organs
    C  Cardiovascular System
    D  Dermatologicals
    G  Genito Urinary System and Sex Hormones
    H  Systemic Hormonal Preparations
    J  Antiinfectives for Systemic Use
    L  Antineoplastic and Immunomodulating Agents
    M  Musculo-Skeletal System
    N  Nervous System
    P  Antiparasitic Products, Insecticides and Repellents
    R  Respiratory System
    S  Sensory Organs
    V  Various
"""

import argparse
import requests
from bs4 import BeautifulSoup
import csv
import time
import sys
import re

BASE_URL = "https://atcddd.fhi.no/atc_ddd_index/"

TOP_LEVEL_CODES = [
    "A", "B", "C", "D", "G", "H", "J", "L", "M", "N", "P", "R", "S", "V"
]

TOP_LEVEL_NAMES = {
    "A": "Alimentary Tract and Metabolism",
    "B": "Blood and Blood Forming Organs",
    "C": "Cardiovascular System",
    "D": "Dermatologicals",
    "G": "Genito Urinary System and Sex Hormones",
    "H": "Systemic Hormonal Preparations",
    "J": "Antiinfectives for Systemic Use",
    "L": "Antineoplastic and Immunomodulating Agents",
    "M": "Musculo-Skeletal System",
    "N": "Nervous System",
    "P": "Antiparasitic Products, Insecticides and Repellents",
    "R": "Respiratory System",
    "S": "Sensory Organs",
    "V": "Various",
}

HEADERS = {
    "User-Agent": "Mozilla/5.0 (compatible; ATC-Fetcher/1.0; research purposes)"
}

# Matches any ATC code pattern: A, A01, A01A, A01AA, A01AA01, etc.
ATC_CODE_RE = re.compile(r'^[A-Z]\d{0,2}[A-Z]{0,2}\d{0,2}$')

DEBUG = False


def fetch_page(code):
    """Fetch ATC index page for a given code."""
    url = f"{BASE_URL}?code={code}&showdescription=no"
    if DEBUG:
        print(f"    [debug] GET {url}")
    try:
        response = requests.get(url, headers=HEADERS, timeout=15)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        print(f"  ERROR fetching {code}: {e}", file=sys.stderr)
        return None


def extract_code_from_href(href):
    """Pull the ATC code out of a ?code=XXX href, or return None."""
    if "?code=" not in href:
        return None
    code = href.split("?code=")[1].split("&")[0].strip().upper()
    return code if ATC_CODE_RE.match(code) else None


def parse_atc_page(html, parent_code):
    """
    Parse an ATC index page.

    Strategy:
    1. Collect every href that contains ?code=  →  candidate sub-codes to recurse.
    2. Find the main data table (has th "ATC code" / "DDD" etc.) and parse its rows.
    3. If no table rows found, the page may be an intermediate grouping page that
       only lists sub-codes as links — those are already captured in step 1.

    Returns:
        sub_codes : list of child ATC codes to recurse into (deduped, excluding parent)
        rows      : list of data dicts for this level
    """
    soup = BeautifulSoup(html, "html.parser")

    # ------------------------------------------------------------------ #
    # 1. Gather all linked sub-codes (any href with ?code=)               #
    # ------------------------------------------------------------------ #
    seen_codes = set()
    sub_codes = []
    for a in soup.find_all("a", href=True):
        code = extract_code_from_href(a["href"])
        if code and code != parent_code and code not in seen_codes:
            # Only keep codes that are one level deeper than the parent
            # (avoids picking up breadcrumb links to ancestor levels)
            if code.upper().startswith(parent_code.upper()) and code != parent_code:
                seen_codes.add(code)
                sub_codes.append(code)

    if DEBUG:
        print(f"    [debug] sub_codes for {parent_code}: {sub_codes[:15]}")

    # ------------------------------------------------------------------ #
    # 2. Parse the data table                                             #
    # ------------------------------------------------------------------ #
    rows = []
    tables = soup.find_all("table")

    if DEBUG:
        print(f"    [debug] tables found: {len(tables)}")
        for i, t in enumerate(tables):
            ths = [th.get_text(strip=True) for th in t.find_all("th")]
            print(f"      table[{i}] headers={ths}  rows={len(t.find_all('tr'))}")

    for table in tables:
        ths = [th.get_text(strip=True).lower() for th in table.find_all("th")]
        joined = " ".join(ths)

        all_trs = table.find_all("tr")

        # If there are no <th> elements, check the first <tr> to see if it
        # looks like a header row made of <td> cells (the site does this).
        # We accept the table if either:
        #   (a) <th> cells mention ATC/DDD/name keywords, OR
        #   (b) there are no <th> at all (headerless table — trust the structure)
        has_th_keywords = any(kw in joined for kw in ("atc", "ddd", "adm", "name"))
        no_headers = len(ths) == 0

        if not has_th_keywords and not no_headers:
            continue

        for tr in all_trs:
            cells = tr.find_all("td")
            if not cells:
                continue

            cell_texts = [c.get_text(separator=" ", strip=True) for c in cells]

            # Pull ATC code: prefer the hyperlink in cell[0]
            code_link = cells[0].find("a", href=True)
            if code_link:
                atc_code = extract_code_from_href(code_link["href"]) or cell_texts[0]
            else:
                atc_code = cell_texts[0]

            # Skip rows where the first cell doesn't look like an ATC code
            if not atc_code or not ATC_CODE_RE.match(atc_code.strip().upper()):
                continue

            atc_code = atc_code.strip().upper()

            rows.append({
                "atc_code": atc_code,
                "name":     cell_texts[1] if len(cell_texts) > 1 else "",
                "ddd":      cell_texts[2] if len(cell_texts) > 2 else "",
                "unit":     cell_texts[3] if len(cell_texts) > 3 else "",
                "adm_r":    cell_texts[4] if len(cell_texts) > 4 else "",
                "note":     cell_texts[5] if len(cell_texts) > 5 else "",
            })

    if DEBUG:
        print(f"    [debug] rows parsed for {parent_code}: {len(rows)}")

    return sub_codes, rows


def crawl_atc(code, visited=None, depth=0):
    """Recursively crawl ATC codes starting from `code`. Returns list of row dicts."""
    if visited is None:
        visited = set()
    if code in visited:
        return []
    visited.add(code)

    indent = "  " * depth
    print(f"{indent}Fetching: {code}")

    html = fetch_page(code)
    if not html:
        return []

    sub_codes, rows = parse_atc_page(html, code)
    time.sleep(0.3)

    for sub_code in sub_codes:
        if sub_code not in visited:
            rows.extend(crawl_atc(sub_code, visited, depth + 1))

    return rows


def resolve_codes(targets):
    """
    Resolve user-supplied targets into ATC codes to crawl.
    Returns (codes_to_fetch, suggested_output_filename).
    """
    normalised = [t.strip().upper() for t in targets]

    if "ALL" in normalised:
        return TOP_LEVEL_CODES, "atc_codes_all.tsv"

    invalid = [c for c in normalised if not c or not ATC_CODE_RE.match(c) or c[0] not in TOP_LEVEL_CODES]
    if invalid:
        print(f"ERROR: Unknown or invalid ATC code(s): {', '.join(invalid)}", file=sys.stderr)
        print(f"Valid top-level codes: {', '.join(TOP_LEVEL_CODES)}", file=sys.stderr)
        print('Use "all" to fetch every category.', file=sys.stderr)
        sys.exit(1)

    filename = (
        f"atc_codes_{normalised[0]}.tsv"
        if len(normalised) == 1
        else "atc_codes_" + "_".join(normalised) + ".tsv"
    )
    return normalised, filename


def deduplicate(rows):
    seen = set()
    unique = []
    for row in rows:
        key = (row["atc_code"], row["name"])
        if key not in seen and row["atc_code"]:
            seen.add(key)
            unique.append(row)
    return unique


def write_tsv(rows, path):
    fieldnames = ["atc_code", "name", "ddd", "unit", "adm_r", "note"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    global DEBUG

    parser = argparse.ArgumentParser(
        description="Fetch ATC/DDD codes from atcddd.fhi.no and save to a TSV file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "codes",
        nargs="+",
        metavar="CODE",
        help='ATC code(s) to fetch, or "all". Examples: all  A  B C  A10  A10BA02',
    )
    parser.add_argument(
        "-o", "--output",
        metavar="FILE",
        default=None,
        help="Output TSV filename (default: auto-generated from code(s))",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print raw parse diagnostics to help troubleshoot empty results",
    )
    args = parser.parse_args()
    DEBUG = args.debug

    codes_to_fetch, default_filename = resolve_codes(args.codes)
    output_path = args.output or default_filename

    print("=" * 60)
    print("ATC Code Fetcher")
    print(f"Source : {BASE_URL}")
    if "ALL" in [c.upper() for c in args.codes]:
        print("Fetching: ALL categories")
    else:
        labels = [
            f"{c} — {TOP_LEVEL_NAMES[c]}" if c in TOP_LEVEL_NAMES else c
            for c in codes_to_fetch
        ]
        print(f"Fetching: {', '.join(labels)}")
    print(f"Output : {output_path}")
    print("=" * 60)

    all_rows = []
    visited = set()

    for code in codes_to_fetch:
        label = TOP_LEVEL_NAMES.get(code, code)
        print(f"\n--- {code}: {label} ---")
        rows = crawl_atc(code, visited, depth=0)
        all_rows.extend(rows)
        print(f"  Collected {len(rows)} entries under {code}")

    unique_rows = deduplicate(all_rows)

    print(f"\nTotal unique rows : {len(unique_rows)}")
    print(f"Writing to        : {output_path}")
    write_tsv(unique_rows, output_path)
    print("Done!")


if __name__ == "__main__":
    main()
