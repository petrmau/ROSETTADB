#!/usr/bin/env python3
"""
export_cluster_drugs.py
=======================
Export a TSV table with one row per (cluster × canonical drug) pair:

    cluster_id  canonical_drug  link_source
    atc_code  inchikey  pubchem_cid  chebi_id  sources  atc_group1  atc_group2

link_source is a pipe-delimited set of evidence sources:
  NCBI        — sequence_drug: NCBI gene_name → gene_drug_link → drug (direct)
  CARD        — sequence_drug_class + drug_class_member: CARD ARO path
  RESFINDER   — sequence_drug_class + drug_class_member: ResFinder path

Two evidence paths are combined:
  Path A (direct):  cluster → gene → sequence_drug → canonical_drug
  Path B (via class): cluster → gene → sequence_drug_class
                              → drug_class_member → canonical_drug

One row is emitted per (cluster_id, canonical_drug) pair; link_source is the
union of all evidence source tokens across every sequence in that cluster.

Options:
    --require-inchikey   Skip drugs that have no InChIKey in amr.drug.

Usage:
    python export_cluster_drugs.py [--dsn <connstr>] [--output <file.tsv>]
                                   [--require-inchikey]

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
# Two evidence paths are unioned before aggregation:
#
# Path A — direct drug link (NCBI only, via amr.sequence_drug):
#   cluster → gene.jrc_id → sequence_drug → canonical_drug
#   evidence_source token: taken directly from sequence_drug.evidence_sources
#
# Path B — via drug class (CARD + RESFINDER, via amr.sequence_drug_class):
#   cluster → gene.jrc_id → sequence_drug_class → drug_class_member → drug
#   evidence_source token: taken from sequence_drug_class.evidence_sources
#   (CARD and RESFINDER write there; NCBI class links also present but drugs
#    already covered by Path A)
#
# Both paths are deduplicated per (cluster_id, canonical_drug) and the source
# tokens are merged into a single sorted pipe-delimited string.
# ---------------------------------------------------------------------------

QUERY = """
WITH raw AS (

    -- Path A: direct sequence → drug (NCBI gene_drug_link path)
    SELECT
        g.cluster_id,
        sd.canonical_drug,
        sdt.source_token
    FROM amr.gene g
    JOIN amr.sequence_drug sd ON sd.jrc_id = g.jrc_id,
    LATERAL unnest(string_to_array(sd.evidence_sources, '|')) AS sdt(source_token)

    UNION

    -- Path B: sequence → drug class → drug (CARD + RESFINDER class path)
    SELECT
        g.cluster_id,
        dcm.canonical_drug,
        sdt.source_token
    FROM amr.gene g
    JOIN amr.sequence_drug_class sdc ON sdc.jrc_id = g.jrc_id
    JOIN amr.drug_class_member   dcm ON dcm.canonical_class = sdc.canonical_class,
    LATERAL unnest(string_to_array(sdc.evidence_sources, '|')) AS sdt(source_token)

),
aggregated AS (
    SELECT
        cluster_id,
        canonical_drug,
        array_to_string(
            array_agg(DISTINCT source_token ORDER BY source_token), '|'
        ) AS link_source
    FROM raw
    GROUP BY cluster_id, canonical_drug
)
SELECT
    a.cluster_id,
    a.canonical_drug,
    a.link_source,
    d.atc_code,
    d.inchikey,
    d.pubchem_cid,
    d.chebi_id,
    d.sources,
    d.atc_group1,
    d.atc_group2
FROM aggregated a
JOIN amr.drug d ON d.canonical_name = a.canonical_drug
{inchikey_filter}
ORDER BY a.cluster_id, a.canonical_drug;
"""

COLUMNS = [
    "cluster_id", "canonical_drug", "link_source",
    "atc_code", "inchikey", "pubchem_cid", "chebi_id",
    "sources", "atc_group1", "atc_group2",
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
    parser.add_argument(
        "--require-inchikey",
        action="store_true",
        help="Skip drugs that have no InChIKey in amr.drug.",
    )
    args = parser.parse_args()

    if not args.dsn:
        print("ERROR: provide --dsn or set $ROSETTADB_DSN", file=sys.stderr)
        sys.exit(1)

    inchikey_filter = "WHERE d.inchikey IS NOT NULL AND d.inchikey <> ''" \
                      if args.require_inchikey else ""
    query = QUERY.format(inchikey_filter=inchikey_filter)

    conn = psycopg2.connect(args.dsn)
    try:
        with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
            print("Running query …", file=sys.stderr)
            cur.execute(query)
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
