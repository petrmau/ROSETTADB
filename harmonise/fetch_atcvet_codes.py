#!/usr/bin/env python3
"""
Fetch ATCvet codes and DDD information from https://atcddd.fhi.no/atcvet/atcvet_index/
and save to a TSV file.

ATCvet codes mirror the human ATC system but carry a "Q" prefix at every level
(e.g. QJ01AA01 = tetracycline for veterinary systemic use).  The top-level
categories are two-letter codes (QA, QB, …, QV) rather than single letters.

Usage:
    python fetch_atcvet_codes.py all                    # Fetch all top-level categories
    python fetch_atcvet_codes.py QJ                     # Fetch a single category
    python fetch_atcvet_codes.py QJ QP                  # Fetch multiple categories
    python fetch_atcvet_codes.py QJ01                   # Fetch a specific sub-level code
    python fetch_atcvet_codes.py QJ01AA01               # Fetch a specific leaf code
    python fetch_atcvet_codes.py all -o my_output.tsv   # Custom output file
    python fetch_atcvet_codes.py QJ --debug             # Debug mode

Top-level categories:
    QA  Alimentary Tract and Metabolism
    QB  Blood and Blood Forming Organs
    QC  Cardiovascular System
    QD  Dermatologicals
    QG  Genito Urinary System and Sex Hormones
    QH  Systemic Hormonal Preparations
    QI  Immunologicals (veterinary-only; no human ATC equivalent)
    QJ  Antiinfectives for Systemic Use
    QL  Antineoplastic and Immunomodulating Agents
    QM  Musculo-Skeletal System
    QN  Nervous System
    QP  Antiparasitic Products, Insecticides and Repellents
    QR  Respiratory System
    QS  Sensory Organs
    QV  Various
"""

import argparse
import requests
from bs4 import BeautifulSoup
import csv
import time
import sys
import re

BASE_URL = "https://atcddd.fhi.no/atcvet/atcvet_index/"

TOP_LEVEL_CODES = [
    "QA", "QB", "QC", "QD", "QG", "QH", "QI",
    "QJ", "QL", "QM", "QN", "QP", "QR", "QS", "QV",
]

TOP_LEVEL_NAMES = {
    "QA": "Alimentary Tract and Metabolism",
    "QB": "Blood and Blood Forming Organs",
    "QC": "Cardiovascular System",
    "QD": "Dermatologicals",
    "QG": "Genito Urinary System and Sex Hormones",
    "QH": "Systemic Hormonal Preparations",
    "QI": "Immunologicals",
    "QJ": "Antiinfectives for Systemic Use",
    "QL": "Antineoplastic and Immunomodulating Agents",
    "QM": "Musculo-Skeletal System",
    "QN": "Nervous System",
    "QP": "Antiparasitic Products, Insecticides and Repellents",
    "QR": "Respiratory System",
    "QS": "Sensory Organs",
    "QV": "Various",
}

HEADERS = {
    "User-Agent": "Mozilla/5.0 (compatible; ATCvet-Fetcher/1.0; research purposes)"
}

# ATCvet codes: Q + human ATC pattern.
# Examples: QJ, QJ01, QJ01A, QJ01AA, QJ01AA01
ATCVET_CODE_RE = re.compile(r'^Q[A-Z]\d{0,2}[A-Z]{0,2}\d{0,2}$')

DEBUG = False
DUMP_HTML = None  # path prefix for dumping raw HTML, set via --dump-html


def fetch_page(code):
    """Fetch ATCvet index page for a given code."""
    url = f"{BASE_URL}?code={code}&showdescription=yes"
    if DEBUG:
        print(f"    [debug] GET {url}")
    try:
        response = requests.get(url, headers=HEADERS, timeout=15)
        response.raise_for_status()
        if DUMP_HTML:
            path = f"{DUMP_HTML}_{code}.html"
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(response.text)
            if DEBUG:
                print(f"    [debug] HTML dumped to {path}")
        return response.text
    except requests.RequestException as e:
        print(f"  ERROR fetching {code}: {e}", file=sys.stderr)
        return None


def extract_code_from_href(href):
    """Pull the ATCvet code out of a ?code=XXX href, or return None."""
    if "?code=" not in href:
        return None
    code = href.split("?code=")[1].split("&")[0].strip().upper()
    return code if ATCVET_CODE_RE.match(code) else None


def parse_atcvet_page(html, parent_code):
    """
    Parse an ATCvet index page.

    Strategy (mirrors fetch_atc_codes.py):
    1. Collect every href that contains ?code=  →  candidate sub-codes to recurse.
    2. Find the main data table and parse its rows.

    Returns:
        sub_codes : list of child ATCvet codes to recurse into
        rows      : list of data dicts for this level
    """
    soup = BeautifulSoup(html, "html.parser")

    # ------------------------------------------------------------------ #
    # 1. Gather all linked sub-codes                                      #
    # ------------------------------------------------------------------ #
    seen_codes = set()
    sub_codes = []
    for a in soup.find_all("a", href=True):
        code = extract_code_from_href(a["href"])
        if code and code != parent_code and code not in seen_codes:
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
        if not tables:
            # Show the tag types and classes present in the main content area
            content = soup.find("div", id="content") or soup.find("main") or soup.body
            if content:
                tag_summary = {}
                for tag in content.find_all(True):
                    key = f"<{tag.name} class={tag.get('class',[])}>"
                    tag_summary[key] = tag_summary.get(key, 0) + 1
                top = sorted(tag_summary.items(), key=lambda x: -x[1])[:20]
                print(f"    [debug] top tags in content: {top}")
                # Also show first 1000 chars of body text
                print(f"    [debug] body snippet:\n{content.get_text()[:500]}")

    for table in tables:
        ths = [th.get_text(strip=True).lower() for th in table.find_all("th")]
        joined = " ".join(ths)

        all_trs = table.find_all("tr")

        # Accept table if it has recognisable column headers OR no headers at all
        # (headerless tables trust structure; ATCvet uses "ATCvet code", "INN/common
        # name", "DDD", "U", "Adm.R" — all covered by the keywords below).
        has_th_keywords = any(
            kw in joined
            for kw in ("atc", "ddd", "adm", "name", "inn", "code", "route", "dose")
        )
        no_headers = len(ths) == 0

        if not has_th_keywords and not no_headers:
            continue

        for tr in all_trs:
            cells = tr.find_all("td")
            if not cells:
                continue

            cell_texts = [c.get_text(separator=" ", strip=True) for c in cells]

            # Try to find an ATCvet code in any cell (link first, then plain text).
            # ATCvet pages sometimes put the code in a column other than column 0.
            atcvet_code = None
            code_col = None
            for ci, cell in enumerate(cells):
                link = cell.find("a", href=True)
                candidate = None
                if link:
                    candidate = extract_code_from_href(link["href"])
                if not candidate:
                    candidate = cell_texts[ci].strip().upper()
                if candidate and ATCVET_CODE_RE.match(candidate):
                    atcvet_code = candidate
                    code_col = ci
                    break

            if not atcvet_code:
                continue

            # Build the row using positions relative to where the code was found
            def _get(i, col=code_col, texts=cell_texts):
                idx = col + i
                return texts[idx] if 0 <= idx < len(texts) else ""

            rows.append({
                "atc_code": atcvet_code,
                "name":     _get(1),
                "ddd":      _get(2),
                "unit":     _get(3),
                "adm_r":    _get(4),
                "note":     _get(5),
            })

    # ------------------------------------------------------------------ #
    # 3. Fallback: parse list-based layout (<li> / <p> / <div>)          #
    #    ATCvet pages may render data as lists rather than tables.        #
    # ------------------------------------------------------------------ #
    if not rows:
        # Walk every element that directly wraps an ATCvet code link.
        # Typical pattern: <li><a href="?code=QD06BB01">QD06BB01</a> name DDD u adm.r</li>
        # We collect text tokens from the parent element after stripping the code itself.
        DDD_UNITS = re.compile(r'^(mg|g|µg|ug|mmol|ml|MU|TU|U|IU|dose|doses?)$', re.IGNORECASE)
        ADM_ROUTES = re.compile(r'^(O|P|N|Inhal|V|TD|SL|R|SL|Impl|GI|loz|chew|TD|gum)$', re.IGNORECASE)
        DDD_VALUE = re.compile(r'^\d[\d.,]*$')

        seen_codes_fallback = set()
        for a_tag in soup.find_all("a", href=True):
            code = extract_code_from_href(a_tag["href"])
            if not code or not code.upper().startswith(parent_code.upper()):
                continue
            if code in seen_codes_fallback:
                continue
            seen_codes_fallback.add(code)

            # Climb to the nearest block-level or list-item parent
            container = a_tag.parent
            while container and container.name not in ("li", "p", "div", "td", "span", "body"):
                container = container.parent
            if not container:
                continue

            # Get the full text of the container, split into tokens
            full_text = container.get_text(separator=" ", strip=True)
            # Remove the code itself from the start of the text
            text_without_code = re.sub(r'^\s*' + re.escape(code) + r'\s*', '', full_text, flags=re.IGNORECASE).strip()
            tokens = text_without_code.split()
            if not tokens:
                # No name info — still record the code
                rows.append({"atc_code": code, "name": "", "ddd": "", "unit": "", "adm_r": "", "note": ""})
                continue

            # Heuristic: find a DDD value (number) surrounded by unit/route tokens
            name_tokens = []
            ddd = unit = adm_r = note = ""
            i = 0
            while i < len(tokens):
                tok = tokens[i]
                if DDD_VALUE.match(tok):
                    ddd = tok
                    # Next token might be unit
                    if i + 1 < len(tokens) and DDD_UNITS.match(tokens[i + 1]):
                        unit = tokens[i + 1]
                        i += 1
                    # Next token might be adm route
                    if i + 1 < len(tokens) and ADM_ROUTES.match(tokens[i + 1]):
                        adm_r = tokens[i + 1]
                        i += 1
                    # Remaining tokens become note
                    note = " ".join(tokens[i + 1:])
                    break
                else:
                    name_tokens.append(tok)
                i += 1

            name = " ".join(name_tokens)
            rows.append({"atc_code": code, "name": name, "ddd": ddd, "unit": unit, "adm_r": adm_r, "note": note})

        if DEBUG:
            print(f"    [debug] fallback list-parser rows: {len(rows)}")

    if DEBUG:
        print(f"    [debug] rows parsed for {parent_code}: {len(rows)}")

    return sub_codes, rows


def crawl_atcvet(code, visited=None, depth=0):
    """Recursively crawl ATCvet codes starting from `code`. Returns list of row dicts."""
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

    sub_codes, rows = parse_atcvet_page(html, code)
    time.sleep(0.3)

    for sub_code in sub_codes:
        if sub_code not in visited:
            rows.extend(crawl_atcvet(sub_code, visited, depth + 1))

    return rows


def resolve_codes(targets):
    """
    Resolve user-supplied targets into ATCvet codes to crawl.
    Returns (codes_to_fetch, suggested_output_filename).
    """
    normalised = [t.strip().upper() for t in targets]

    if "ALL" in normalised:
        return TOP_LEVEL_CODES, "atcvet_codes_all.tsv"

    invalid = [c for c in normalised if not c or not ATCVET_CODE_RE.match(c)]
    if invalid:
        print(f"ERROR: Unknown or invalid ATCvet code(s): {', '.join(invalid)}", file=sys.stderr)
        print(f"Valid top-level codes: {', '.join(TOP_LEVEL_CODES)}", file=sys.stderr)
        print('Use "all" to fetch every category.', file=sys.stderr)
        sys.exit(1)

    filename = (
        f"atcvet_codes_{normalised[0]}.tsv"
        if len(normalised) == 1
        else "atcvet_codes_" + "_".join(normalised) + ".tsv"
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
        description="Fetch ATCvet codes from atcddd.fhi.no/atcvet and save to a TSV file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "codes",
        nargs="+",
        metavar="CODE",
        help='ATCvet code(s) to fetch, or "all". Examples: all  QJ  QP QD  QJ01  QJ01AA01',
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
    parser.add_argument(
        "--dump-html",
        metavar="PREFIX",
        default=None,
        help="Save raw HTML for each fetched page to PREFIX_CODE.html (useful for debugging)",
    )
    args = parser.parse_args()
    DEBUG = args.debug
    global DUMP_HTML
    DUMP_HTML = args.dump_html

    codes_to_fetch, default_filename = resolve_codes(args.codes)
    output_path = args.output or default_filename

    print("=" * 60)
    print("ATCvet Code Fetcher")
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
        rows = crawl_atcvet(code, visited, depth=0)
        all_rows.extend(rows)
        print(f"  Collected {len(rows)} entries under {code}")

    unique_rows = deduplicate(all_rows)

    print(f"\nTotal unique rows : {len(unique_rows)}")
    print(f"Writing to        : {output_path}")
    write_tsv(unique_rows, output_path)
    print("Done!")


if __name__ == "__main__":
    main()
