# ROSETTADB

PostgreSQL database for Antimicrobial Resistance (AMR) genes aggregated from
multiple public sources and clustered at 100% sequence identity.

## Sources

| Source | File | Sequences |
|---|---|---|
| CARD | `sources/CARD/nucleotide_fasta_protein_homolog_model.fasta` | ~6 052 |
| NCBI AMRFinderPlus | `sources/amr_finder_plus/ncbi_dataset/data/nucleotide.fna` | ~8 302 |
| ResFinder | `sources/resfinder_db/all.fsa` | TBD |

ResFinder records are enriched with drug class, resistance mechanism, and
phenotype data from `sources/resfinder_db/phenotypes.txt`.

## Sequence identifier

Every unique nucleotide sequence receives a stable JRC identifier:

```
JRC<first-10-chars-of-MD5(sequence)>
```

e.g. `JRC8214f7c788`

## Clustering

Sequences with 100% identity (same MD5) are grouped into a single cluster.
The **representative** is the member with the longest `original_header`.

Each cluster is assigned a stable string identifier of the form:

```
JRCGRP_000001
```

## Schema

```
amr.sequence           – one row per unique sequence (jrc_id PK, md5, raw sequence)
amr.cluster            – one row per cluster; cluster_id is JRCGRP_XXXXXX;
                         points to the representative jrc_id
amr.gene               – one row per FASTA record; links to sequence + cluster;
                         carries all source metadata (gene name, drug class, ARO, …)
amr.sequence_metadata  – canonical metadata per (jrc_id, source) after deduplication;
                         one row per unique sequence × source combination
```

### amr.sequence_metadata

This table provides a deduplicated, canonical view of metadata for each unique
sequence per source. It is populated/updated via upsert during ingestion and is
useful for cross-source comparisons without duplicates from the raw `amr.gene` table.

## Setup

### 1. Create the database

```bash
createdb rosettadb
```

### 2. Install Python dependencies

```bash
pip install psycopg2-binary biopython
```

### 3. Run the ingestion script

```bash
python ingest.py \
  --dsn "host=localhost dbname=rosettadb user=postgres password=postgres"
```

Optional overrides for non-default file paths:

```
--resfinder       PATH   (default: sources/resfinder_db/all.fsa)
--resfinder-pheno PATH   (default: sources/resfinder_db/phenotypes.txt)
--card            PATH   (default: sources/CARD/nucleotide_fasta_protein_homolog_model.fasta)
--ncbi            PATH   (default: sources/amr_finder_plus/ncbi_dataset/data/nucleotide.fna)
```

The script is idempotent (`ON CONFLICT DO NOTHING` / `ON CONFLICT DO UPDATE`) and safe to re-run.

## Useful queries

```sql
-- Total unique sequences
SELECT count(*) FROM amr.sequence;

-- Sequences per source
SELECT source, count(*) FROM amr.gene GROUP BY source;

-- Sequences shared between CARD and NCBI
SELECT s.jrc_id, count(DISTINCT g.source) AS n_sources
FROM amr.sequence s
JOIN amr.gene g USING (jrc_id)
GROUP BY s.jrc_id
HAVING count(DISTINCT g.source) > 1
LIMIT 10;

-- Cluster representative with all its aliases
SELECT c.cluster_id, c.representative_jrc, g.source, g.original_header
FROM amr.cluster c
JOIN amr.gene g ON g.jrc_id = c.representative_jrc
ORDER BY c.cluster_id;

-- Canonical metadata for a sequence across all sources
SELECT * FROM amr.sequence_metadata WHERE jrc_id = 'JRC8214f7c788';

-- Drug classes covered per source (deduplicated)
SELECT source, drug_class, count(*) AS n_sequences
FROM amr.sequence_metadata
WHERE drug_class IS NOT NULL
GROUP BY source, drug_class
ORDER BY source, n_sequences DESC;
```
