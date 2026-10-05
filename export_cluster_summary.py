#!/usr/bin/env python3
"""
export_cluster_summary.py
=========================
Export one TSV row per gene-cluster representative with all AMR annotation.

Columns
-------
cluster_id           JRCGRP_000042
representative_jrc   JRC3a7f8b1c2d
cluster_size         9
sources_present      CARD|NCBI|RESFINDER  (sorted, pipe-delimited)
is_coding            TRUE / FALSE  (TRUE iff sequence has an assigned protein_id)
protein_id           JRCPRO_000012 or NULL
gene_name_best       best name: CARD > NCBI > RESFINDER priority
aro_accession        ARO:3000078 or NULL
gene_family          from CARD aro_index or NCBI ReferenceGeneCatalog
resistance_mechanism from CARD aro_index or ResFinder phenotypes
gene_names_by_source CARD:TEM-1|NCBI:blaTEM-1|RESFINDER:TEM-1  (DB:NULL if absent)
drug_classes         beta-lactam antibiotic[CARD|CARD.aro]|...  (evidence in brackets)
drug_names           amoxicillin[CARD]|ampicillin[CARD|CARD.aro] (evidence in brackets)
drug_inchikeys       LSQZJLSUYDQPKJ-...|NULL  (positionally aligned with drug_names)
drug_atc_codes       J01CA04|J01CA01          (positionally aligned with drug_names)
drug_atc_groups      Penicillins|Penicillins  (positionally aligned with drug_names)

Positional alignment: position i across drug_names / drug_inchikeys /
drug_atc_codes / drug_atc_groups always refers to the same drug.  Missing
values are written as the literal string "NULL".

Usage:
    python export_cluster_summary.py --dsn "host=... dbname=..." --output summary.tsv

Output goes to stdout if --output is not given.
"""

import argparse
import os
import sys
from collections import defaultdict

import psycopg2
import psycopg2.extras


# ---------------------------------------------------------------------------
# SQL — one query pulls all per-cluster data; Python post-processes into rows
# ---------------------------------------------------------------------------

SQL_CLUSTERS = """
SELECT
    cl.cluster_id,
    cl.representative_jrc,
    -- cluster_size: number of gene FASTA records in the cluster
    (
        SELECT count(*)
        FROM amr.gene g_sz
        WHERE g_sz.cluster_id = cl.cluster_id
    ) AS cluster_size,
    -- sources_present: all distinct sources in the cluster
    (
        SELECT array_to_string(
            array_agg(DISTINCT g2.source ORDER BY g2.source), '|'
        )
        FROM amr.gene g2
        WHERE g2.cluster_id = cl.cluster_id
    ) AS sources_present,
    -- coding / protein_id from the representative sequence
    (seq.protein_id IS NOT NULL) AS is_coding,
    seq.protein_id AS protein_id,
    -- CARD metadata (from the cluster's CARD gene, if any)
    card_g.gene_name AS card_gene_name,
    card_g.product_name AS card_product_name,
    card_g.aro_accession AS aro_accession,
    card_g.gene_family AS gene_family_card,
    card_g.resistance_mechanism AS resistance_mechanism_card,
    -- NCBI metadata
    ncbi_g.gene_name AS ncbi_gene_name,
    ncbi_g.product_name AS ncbi_product_name,
    ncbi_g.gene_family AS gene_family_ncbi,
    -- ResFinder metadata
    rf_g.gene_name AS rf_gene_name
FROM amr.cluster cl
-- representative sequence row
JOIN amr.sequence seq ON seq.jrc_id = cl.representative_jrc
-- CARD gene for this representative (may be NULL)
LEFT JOIN LATERAL (
    SELECT g.gene_name, g.product_name, g.aro_accession,
           g.gene_family, g.resistance_mechanism
    FROM amr.gene g
    WHERE g.jrc_id = cl.representative_jrc AND g.source = 'CARD'
    LIMIT 1
) card_g ON TRUE
-- NCBI gene for this representative
LEFT JOIN LATERAL (
    SELECT g.gene_name, g.product_name, g.gene_family
    FROM amr.gene g
    WHERE g.jrc_id = cl.representative_jrc AND g.source = 'NCBI'
    LIMIT 1
) ncbi_g ON TRUE
-- ResFinder gene for this representative
LEFT JOIN LATERAL (
    SELECT g.gene_name
    FROM amr.gene g
    WHERE g.jrc_id = cl.representative_jrc AND g.source = 'RESFINDER'
    LIMIT 1
) rf_g ON TRUE
ORDER BY cl.cluster_id;
"""

SQL_DRUG_CLASSES = """
SELECT
    g.cluster_id,
    sdc.canonical_class,
    sdc.evidence_sources
FROM amr.gene g
JOIN amr.sequence_drug_class sdc ON sdc.jrc_id = g.jrc_id
WHERE g.jrc_id IN (
    SELECT representative_jrc FROM amr.cluster
)
ORDER BY g.cluster_id, sdc.canonical_class;
"""

SQL_DRUG_NAMES = """
SELECT
    g.cluster_id,
    sd.canonical_drug,
    sd.evidence_sources,
    d.inchikey,
    d.atc_code,
    d.atc_group1
FROM amr.gene g
JOIN amr.sequence_drug sd ON sd.jrc_id = g.jrc_id
LEFT JOIN amr.drug d ON d.canonical_name = sd.canonical_drug
WHERE g.jrc_id IN (
    SELECT representative_jrc FROM amr.cluster
)
ORDER BY g.cluster_id, sd.canonical_drug;
"""

COLUMNS = [
    "cluster_id",
    "representative_jrc",
    "cluster_size",
    "sources_present",
    "is_coding",
    "protein_id",
    "gene_name_best",
    "aro_accession",
    "gene_family",
    "resistance_mechanism",
    "gene_names_by_source",
    "drug_classes",
    "drug_names",
    "drug_inchikeys",
    "drug_atc_codes",
    "drug_atc_groups",
]


def _null(v):
    """Return 'NULL' string for None or empty, else the value as string."""
    if v is None or (isinstance(v, str) and v.strip() == ""):
        return "NULL"
    return str(v)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--dsn",
        default=os.environ.get(
            "ROSETTADB_DSN",
            "host=localhost dbname=rosettadb user=postgres password=postgres",
        ),
        help="PostgreSQL DSN (default: $ROSETTADB_DSN or localhost)",
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output TSV path (default: stdout)",
    )
    args = parser.parse_args()

    conn = psycopg2.connect(args.dsn)
    try:
        with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
            # --- drug classes per cluster representative ---
            print("Fetching drug classes …", file=sys.stderr)
            cur.execute(SQL_DRUG_CLASSES)
            # {cluster_id: [(canonical_class, evidence_sources), ...]}
            drug_classes_by_cluster: dict[str, list] = defaultdict(list)
            for row in cur.fetchall():
                drug_classes_by_cluster[row["cluster_id"]].append(
                    (row["canonical_class"], row["evidence_sources"])
                )

            # --- drug names + ATC / InChIKey per cluster representative ---
            print("Fetching drug names …", file=sys.stderr)
            cur.execute(SQL_DRUG_NAMES)
            # {cluster_id: [(canonical_drug, evidence_sources, inchikey, atc_code, atc_group1), ...]}
            drug_names_by_cluster: dict[str, list] = defaultdict(list)
            for row in cur.fetchall():
                drug_names_by_cluster[row["cluster_id"]].append((
                    row["canonical_drug"],
                    row["evidence_sources"],
                    row["inchikey"],
                    row["atc_code"],
                    row["atc_group1"],
                ))

            # --- cluster representative metadata ---
            print("Fetching cluster metadata …", file=sys.stderr)
            cur.execute(SQL_CLUSTERS)
            cluster_rows = cur.fetchall()
            print(f"Fetched {len(cluster_rows)} clusters.", file=sys.stderr)
    finally:
        conn.close()

    out = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout
    try:
        out.write("\t".join(COLUMNS) + "\n")
        for row in cluster_rows:
            cid = row["cluster_id"]

            # ── gene_name_best: CARD > NCBI > RESFINDER ──────────────────
            gene_name_best = (
                row["card_gene_name"]
                or row["ncbi_gene_name"]
                or row["rf_gene_name"]
                or "NULL"
            )

            # ── gene_family: CARD > NCBI ─────────────────────────────────
            gene_family = _null(row["gene_family_card"] or row["gene_family_ncbi"])

            # ── resistance_mechanism: CARD > ResFinder ───────────────────
            # (ResFinder phenotypes.txt has a 'resistance_mechanism' column
            # if loaded; CARD aro_index has AMR Gene Family / mechanism)
            resistance_mechanism = _null(row["resistance_mechanism_card"])

            # ── gene_names_by_source ──────────────────────────────────────
            parts = [
                f"CARD:{_null(row['card_gene_name'])}",
                f"NCBI:{_null(row['ncbi_gene_name'])}",
                f"RESFINDER:{_null(row['rf_gene_name'])}",
            ]
            gene_names_by_source = "|".join(parts)

            # ── drug classes (evidence in brackets) ───────────────────────
            dc_list = sorted(drug_classes_by_cluster.get(cid, []),
                             key=lambda x: x[0])
            drug_classes_str = "|".join(
                f"{cls}[{ev}]" for cls, ev in dc_list
            ) or "NULL"

            # ── drug names + positionally aligned columns ─────────────────
            dn_list = sorted(drug_names_by_cluster.get(cid, []),
                             key=lambda x: x[0])
            drug_names_str = "|".join(
                f"{drug}[{ev}]" for drug, ev, *_ in dn_list
            ) or "NULL"

            drug_inchikeys_str = "|".join(
                _null(inchikey) for _, _, inchikey, *_ in dn_list
            ) or "NULL"

            drug_atc_codes_str = "|".join(
                _null(atc) for _, _, _, atc, _ in dn_list
            ) or "NULL"

            drug_atc_groups_str = "|".join(
                _null(grp) for _, _, _, _, grp in dn_list
            ) or "NULL"

            values = [
                _null(cid),
                _null(row["representative_jrc"]),
                _null(row["cluster_size"]),
                _null(row["sources_present"]),
                "TRUE" if row["is_coding"] else "FALSE",
                _null(row["protein_id"]),
                gene_name_best,
                _null(row["aro_accession"]),
                gene_family,
                resistance_mechanism,
                gene_names_by_source,
                drug_classes_str,
                drug_names_str,
                drug_inchikeys_str,
                drug_atc_codes_str,
                drug_atc_groups_str,
            ]
            out.write("\t".join(values) + "\n")
    finally:
        if args.output:
            out.close()
            print(f"Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
