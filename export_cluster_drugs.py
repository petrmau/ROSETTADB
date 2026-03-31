#!/usr/bin/env python3
"""
export_cluster_drugs.py
=======================
Export a TSV table with one row per (cluster × canonical drug) pair:

    cluster_id  canonical_drug  link_source

link_source records how the drug was linked to the cluster:
  NCBI        — via amr.sequence_drug (NCBI gene_name → gene_drug_link → drug)
  CARD        — via amr.sequence_drug_class + drug_class_member (CARD ARO path)
  RESFINDER   — via amr.sequence_drug_class + drug_class_member (ResFinder path)
  NCBI|CARD, NCBI|CARD|RESFINDER, … — multiple evidence paths (pipe-delimited)

The same drug may be reachable from multiple sequences within a cluster.
One row is emitted per (cluster_id, canonical_drug) pair; link_source is the
union of all evidence sources across every sequence in that cluster.

Usage:
    python export_cluster_drugs.py [--dsn <connstr>] [--output <file.tsv>]

Output goes to stdout if --output is not given.
"""

import argparse
import os
import sys

import psycopg2
import psycopg2.extras

# ---------------------------------------------------------------------------
# Query
# ---------------------------------------------------------------------------
# Aggregates evidence_sources across all sequences in a cluster, then merges
# the pipe-delimited tokens into a single sorted, deduplicated pipe string.
#
# Evidence path:
#   cluster → gene.jrc_id → sequence_drug.canonical_drug  (direct NCBI link)
#   cluster → gene.jrc_id → sequence_drug_class           (class path)
#                         → drug_class_member.canonical_drug
#
# Both paths write their evidence_sources into amr.sequence_drug already, so
# a single join to sequence_drug is sufficient — NCBI, CARD, and RESFINDER
# labels are already encoded there by ingest.py.
# ---------------------------------------------------------------------------

QUERY = """
WITH cluster_drug AS (
    SELECT
        g.cluster_id,
        sd.canonical_drug,
        sd.evidence_sources
    FROM amr.gene g
    JOIN amr.sequence_drug sd ON sd.jrc_id = g.jrc_id
),
aggregated AS (
    SELECT
        cluster_id,
        canonical_drug,
        -- collect all unique source tokens across every sequence in the cluster
        array_agg(DISTINCT source_token ORDER BY source_token) AS source_tokens
    FROM cluster_drug,
         -- unnest the pipe-delimited evidence_sources into individual tokens
         LATERAL unnest(string_to_array(evidence_sources, '|')) AS source_token
    GROUP BY cluster_id, canonical_drug
)
SELECT
    cluster_id,
    canonical_drug,
    array_to_string(source_tokens, '|') AS link_source
FROM aggregated
ORDER BY cluster_id, canonical_drug;
"""

COLUMNS = ["cluster_id", "canonical_drug", "link_source"]


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
            out.write(
                "\t".join(str(row[col]) if row[col] is not None else "" for col in COLUMNS)
                + "\n"
            )
    finally:
        if args.output:
            out.close()
            print(f"Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
