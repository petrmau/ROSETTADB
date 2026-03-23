#!/usr/bin/env python3
"""
export_gene_clusters.py
=======================
Export a TSV table with one row per gene sequence cluster:

    cluster_id  representative_jrc  representative_sequence  protein_id  protein_sequence
    gene_name_card  gene_name_ncbi  gene_name_resfinder

Usage:
    python export_gene_clusters.py [--dsn <connstr>] [--output <file.tsv>]

Output goes to stdout if --output is not given.
"""

import argparse
import os
import sys

import psycopg2
import psycopg2.extras

QUERY = """
SELECT
    c.cluster_id,
    c.representative_jrc,
    s.sequence                          AS representative_sequence,
    p.protein_id,
    p.protein_sequence,
    card.gene_name                      AS gene_name_card,
    ncbi.gene_name                      AS gene_name_ncbi,
    resfinder.gene_name                 AS gene_name_resfinder
FROM amr.cluster c
JOIN amr.sequence s  ON s.jrc_id          = c.representative_jrc
JOIN amr.protein  p  ON p.representative_jrc = c.representative_jrc
LEFT JOIN amr.sequence_metadata card
    ON card.jrc_id = c.representative_jrc AND card.source = 'CARD'
LEFT JOIN amr.sequence_metadata ncbi
    ON ncbi.jrc_id = c.representative_jrc AND ncbi.source = 'NCBI'
LEFT JOIN amr.sequence_metadata resfinder
    ON resfinder.jrc_id = c.representative_jrc AND resfinder.source = 'RESFINDER'
ORDER BY c.cluster_id;
"""

COLUMNS = [
    "cluster_id",
    "representative_jrc",
    "representative_sequence",
    "protein_id",
    "protein_sequence",
    "gene_name_card",
    "gene_name_ncbi",
    "gene_name_resfinder",
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
