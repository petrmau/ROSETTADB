#!/usr/bin/env python3
"""
export_protein_clusters.py
==========================
Export a TSV table with one row per protein cluster:

    protein_id  protein_sequence  cluster_id  representative_jrc  representative_sequence

Usage:
    python export_protein_clusters.py [--dsn <connstr>] [--output <file.tsv>]

Output goes to stdout if --output is not given.
"""

import argparse
import os
import sys

import psycopg2
import psycopg2.extras

QUERY = """
SELECT
    p.protein_id,
    p.protein_sequence,
    c.cluster_id,
    c.representative_jrc,
    s.sequence                  AS representative_sequence
FROM amr.protein p
JOIN amr.sequence s  ON s.jrc_id          = p.representative_jrc
JOIN amr.cluster  c  ON c.representative_jrc = p.representative_jrc
ORDER BY p.protein_id;
"""

COLUMNS = [
    "protein_id",
    "protein_sequence",
    "cluster_id",
    "representative_jrc",
    "representative_sequence",
]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--dsn",
        default=os.environ.get("ROSETTADB_DSN", ""),
        help="PostgreSQL connection string (default: $ROSETTADB_DSN)",
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output TSV file path (default: stdout)",
    )
    args = parser.parse_args()

    if not args.dsn:
        print("ERROR: provide --dsn or set $ROSETTADB_DSN", file=sys.stderr)
        sys.exit(1)

    conn = psycopg2.connect(args.dsn)
    try:
        with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
            print("Running query …", file=sys.stderr)
            cur.execute(QUERY)
            rows = cur.fetchall()
            print(f"Fetched {len(rows)} rows.", file=sys.stderr)
    finally:
        conn.close()

    out = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout
    try:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(str(row[col]) if row[col] is not None else "" for col in COLUMNS) + "\n")
    finally:
        if args.output:
            out.close()
            print(f"Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
