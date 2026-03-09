# ROSETTADB

PostgreSQL database for Antimicrobial Resistance (AMR) genes aggregated from
multiple public sources, clustered at 100% sequence identity, and linked to a
harmonised drug/antibiotic class vocabulary.

---

## Table of contents

1. [Sources](#sources)
2. [Schema](#schema)
3. [Sequence identifier](#sequence-identifier)
4. [Clustering](#clustering)
5. [Harmonised drug vocabulary](#harmonised-drug-vocabulary)
6. [Procedures](#procedures)
   - [First-time setup](#first-time-setup)
   - [Re-ingest sequences only](#re-ingest-sequences-only)
   - [Update drug vocabulary and re-create links](#update-drug-vocabulary-and-re-create-links)
   - [Full rebuild from scratch](#full-rebuild-from-scratch)
7. [ingest.py flags](#ingestpy-flags)
8. [Useful queries](#useful-queries)

---

## Sources

| Source | Sequences | Metadata |
|--------|-----------|----------|
| CARD (v4.0.1) | ~6 052 nucleotide sequences | ARO accession, gene family, drug class, resistance mechanism |
| NCBI AMRFinderPlus | ~8 302 sequences | Gene family, class, subclass, scope, type |
| ResFinder | ~3 212 sequences | Drug class, phenotype, PMID, required co-gene |

---

## Schema

### Core sequence tables

```
amr.sequence           – one row per unique nucleotide sequence
                         (jrc_id PK, MD5, raw sequence, length)
amr.cluster            – one row per 100%-identity cluster (JRCGRP_XXXXXX);
                         points to the representative jrc_id
amr.gene               – one row per FASTA record; FK → sequence + cluster;
                         carries all source metadata (gene name, drug class, ARO, …)
amr.sequence_metadata  – canonical metadata per (jrc_id × source) after dedup;
                         useful for cross-source comparisons without raw duplicates
```

### Harmonised drug vocabulary

```
amr.drug_class         – 64 canonical antibiotic/antimicrobial classes
                         (names follow CARD ARO; ARO accessions are stable keys)
amr.drug               – 169 canonical drugs (INN names); enriched with
                         PubChem CID, InChIKey, ATC code, ChEBI ID
amr.drug_alias         – 80 spelling variants, brand names, source-specific aliases
amr.drug_class_member  – 172 drug → class memberships (many-to-many);
                         provenance tracked (resfinder | ncbi | curated)
```

### Gene → drug class links (intermediate)

```
amr.gene_drug_link     – 10 388 NCBI AMRFinderPlus gene × drug/class pairings
                         (raw class and subclass tokens + normalised canonical tokens)
amr.card_gene_class    – 13 361 CARD ARO gene × canonical drug class links
                         (unique ARO models × drug classes)
```

### Sequence → drug class / drug (final bridge tables)

```
amr.sequence_drug_class – 20 713 rows; one row per unique sequence × canonical drug class
                          Evidence paths:
                            CARD:      gene.aro_accession → card_gene_class → drug_class
                            NCBI:      gene.amr_class (slash-split) → drug_class.ncbi_alias
                                       + gene.gene_name → gene_drug_link → drug_class
                            ResFinder: gene.drug_class text → drug_class.resfinder_alias
                          Coverage: 99.1 % of 9 222 unique sequences
amr.sequence_drug       –  3 376 rows; one row per unique sequence × canonical drug name
                          (NCBI gene_name → gene_drug_link → drug)
```

These are the primary output tables answering: *"for this non-redundant sequence,
what drug classes and antibiotics does it confer resistance to?"*

---

## Sequence identifier

Every unique nucleotide sequence receives a stable JRC identifier:

```
JRC<first-10-chars-of-MD5(sequence)>
```

e.g. `JRC8214f7c788`

---

## Clustering

Sequences sharing 100% identity (same MD5) are grouped into one cluster.
The **representative** is the member with the longest `original_header`.
Each cluster gets a stable identifier of the form `JRCGRP_000001`.

---

## Harmonised drug vocabulary

The `harmonise/` directory contains pipeline scripts and the generated TSVs
loaded into the database:

| File | Rows | Description |
|------|------|-------------|
| `class_mapping.tsv` | 64 | Canonical drug classes with ARO accessions, ResFinder/NCBI aliases |
| `drug_canonical.tsv` | 169 | Canonical drugs (INN) with PubChem/ATC/ChEBI identifiers |
| `drug_alias.tsv` | 80 | Source-specific aliases and abbreviations |
| `drug_class_member.tsv` | 172 | Drug → class links with evidence provenance |
| `gene_drug_link.tsv` | 10 388 | NCBI gene → drug/class pairings (raw + normalised) |
| `card_gene_class.tsv` | 13 361 | CARD ARO gene → canonical drug class links |

### Regenerate harmonised TSVs

Run these scripts whenever source data is updated. They write new TSVs to
`harmonise/` but do **not** touch the database — run `ingest.py --harmonise-only`
afterwards to reload.

```bash
# Rebuild drug/class tables from source data
python harmonise/harmonise.py

# Enrich drug_canonical.tsv with cross-database identifiers (ChEBI → PubChem)
python harmonise/enrich.py

# Rebuild CARD gene → class links from sources/CARD/aro_index.tsv
python harmonise/parse_card_aro.py
```

#### Enrichment strategy (`enrich.py`)

Identifiers are resolved in this order per drug:

| Step | API | Fields retrieved |
|------|-----|-----------------|
| 1 | ChEBI by INN name | `chebi_id`, `inchikey` — picks highest-star exact-name match |
| 2a | PubChem by **InChIKey** (if ChEBI hit) | `pubchem_cid`, `atc_code` — unambiguous |
| 2b | PubChem by **name** (fallback) | `pubchem_cid`, `inchikey`, `atc_code` |
| 3 | PubChem xrefs (last resort) | `chebi_id` if step 1 failed but CID was found |

ChEBI is used first because it is fully curated (3-star entries are manually reviewed),
and an InChIKey-based PubChem lookup avoids the salt/stereoisomer ambiguity of
name-based searches. ATC codes come exclusively from PubChem (ChEBI does not index them).

Cache: `harmonise/.enrich_cache.json` — all 169 current drugs are pre-cached; re-runs
make zero live API calls unless new drugs are added or cache entries are cleared.

To upgrade the 45 existing entries that have a PubChem CID but no ChEBI ID:

```bash
python harmonise/enrich.py --refresh-missing-chebi
```

---

## Procedures

### First-time setup

Complete steps to go from an empty PostgreSQL server to a fully populated
ROSETTADB.

**1. Install Python dependencies**

```bash
pip install psycopg2-binary
```

**2. Create the database**

```bash
createdb rosettadb
# or with explicit connection params:
psql -c "CREATE DATABASE rosettadb;" postgres
```

**3. Run the full ingestion**

```bash
python ingest.py \
  --dsn "host=localhost dbname=rosettadb user=postgres password=postgres"
```

This single command:
- Applies `schema.sql` (creates all tables and indexes)
- Parses the three FASTA sources, deduplicates at 100% identity, assigns
  `jrc_id` and `JRCGRP_` cluster IDs
- Inserts all gene records and canonical per-source metadata
- Loads all harmonised drug vocabulary tables from `harmonise/`
- Populates `sequence_drug_class` and `sequence_drug` from all evidence paths

The script is fully idempotent — safe to re-run (`ON CONFLICT DO UPDATE`).

---

### Re-ingest sequences only

Use this when source FASTA files change (new CARD release, updated NCBI
AMRFinderPlus, updated ResFinder) but the drug vocabulary has not changed.

```bash
python ingest.py \
  --dsn "host=localhost dbname=rosettadb user=postgres password=postgres" \
  --skip-harmonise
```

This re-parses FASTAs, upserts sequences/clusters/genes, and re-populates
the bridge tables (`sequence_drug_class`, `sequence_drug`) using the existing
vocabulary already in the DB. Drug vocabulary tables are left untouched.

> **Note:** existing `jrc_id` values for unchanged sequences are stable —
> they are derived from the sequence MD5 and will not change.

---

### Update drug vocabulary and re-create links

Use this after editing any `harmonise/` TSV (e.g. adding a new drug class,
correcting a mapping, adding an alias). Sequences do **not** need to be
re-ingested.

**Step 1 — regenerate TSVs** (only if the source scripts were modified):

```bash
python harmonise/harmonise.py      # drug classes + drugs + aliases + class membership
python harmonise/enrich_pubchem.py # PubChem/ATC/ChEBI identifiers
python harmonise/parse_card_aro.py # CARD gene → class links
```

Skip any script whose TSV you edited manually.

**Step 2 — reload into the database**:

```bash
python ingest.py \
  --dsn "host=localhost dbname=rosettadb user=postgres password=postgres" \
  --harmonise-only
```

This:
- Upserts all six harmonise TSVs into the drug vocabulary tables
- Truncates `sequence_drug_class` and `sequence_drug`
- Re-populates both bridge tables from scratch using the updated vocabulary

Runs in seconds — no FASTA parsing.

---

### Full rebuild from scratch

Use this to drop all data and reload everything cleanly (e.g. after a major
schema change).

```bash
# Drop and recreate the database
dropdb rosettadb
createdb rosettadb

# Full ingest
python ingest.py \
  --dsn "host=localhost dbname=rosettadb user=postgres password=postgres"
```

---

## ingest.py flags

| Flag | Default | Description |
|------|---------|-------------|
| `--dsn` | `host=localhost dbname=rosettadb user=postgres password=postgres` | PostgreSQL connection string |
| `--resfinder` | `sources/resfinder_db/all.fsa` | ResFinder FASTA path |
| `--card` | `sources/CARD/nucleotide_fasta_protein_homolog_model.fasta` | CARD FASTA path |
| `--ncbi` | `sources/amr_finder_plus/ncbi_dataset/data/nucleotide.fna` | NCBI FASTA path |
| `--schema` | `schema.sql` | DDL file path |
| `--skip-harmonise` | — | Skip loading harmonised drug vocabulary tables |
| `--skip-links` | — | Skip populating `sequence_drug_class` / `sequence_drug` tables |
| `--harmonise-only` | — | Skip FASTA parsing; only reload `harmonise/` TSVs and re-create sequence links (sequences must already be in DB) |

---

## Useful queries

```sql
-- Sequences per source
SELECT source, count(*) FROM amr.gene GROUP BY source;

-- Coverage summary: how many sequences have a drug class / drug assigned
SELECT
    count(DISTINCT jrc_id)                                        AS total_sequences,
    count(DISTINCT jrc_id) FILTER (WHERE jrc_id IN
        (SELECT jrc_id FROM amr.sequence_drug_class))             AS seqs_with_class,
    count(DISTINCT jrc_id) FILTER (WHERE jrc_id IN
        (SELECT jrc_id FROM amr.sequence_drug))                   AS seqs_with_drug
FROM amr.sequence;

-- All drug classes linked to a specific sequence (primary bridge table)
SELECT sdc.canonical_class, sdc.evidence_sources
FROM amr.sequence_drug_class sdc
WHERE sdc.jrc_id = 'JRC8214f7c788'
ORDER BY sdc.canonical_class;

-- All drug names linked to a specific sequence
SELECT sd.canonical_drug, sd.evidence_sources
FROM amr.sequence_drug sd
WHERE sd.jrc_id = 'JRC8214f7c788';

-- Full picture: sequence → gene name → drug class → individual drugs
SELECT DISTINCT
    sdc.jrc_id,
    g.gene_name,
    sdc.canonical_class,
    d.canonical_name  AS drug,
    d.atc_code,
    sdc.evidence_sources
FROM amr.sequence_drug_class sdc
JOIN amr.gene g              ON g.jrc_id = sdc.jrc_id
LEFT JOIN amr.drug_class_member dcm ON dcm.canonical_class = sdc.canonical_class
LEFT JOIN amr.drug d         ON d.canonical_name = dcm.canonical_drug
WHERE sdc.canonical_class = 'aminoglycoside antibiotic'
ORDER BY g.gene_name, d.canonical_name
LIMIT 20;

-- Sequences shared between two or more sources
SELECT s.jrc_id, count(DISTINCT g.source) AS n_sources,
       string_agg(DISTINCT g.source, '|' ORDER BY g.source) AS sources
FROM amr.sequence s
JOIN amr.gene g USING (jrc_id)
GROUP BY s.jrc_id
HAVING count(DISTINCT g.source) > 1
LIMIT 10;

-- All gene sequences conferring resistance to fluoroquinolones (via CARD)
SELECT cgc.gene_name, cgc.card_short_name, cgc.gene_family, cgc.resistance_mechanism
FROM amr.card_gene_class cgc
WHERE cgc.canonical_class = 'fluoroquinolone antibiotic'
ORDER BY cgc.gene_family, cgc.gene_name;

-- Drug class membership for a specific drug (with ARO and category)
SELECT dm.canonical_drug, dc.canonical_name AS drug_class,
       dc.aro_accession, dc.category, dm.evidence_source
FROM amr.drug_class_member dm
JOIN amr.drug_class dc ON dc.canonical_name = dm.canonical_class
WHERE dm.canonical_drug = 'ciprofloxacin';

-- All drugs in a class with their cross-database identifiers
SELECT d.canonical_name, d.atc_code, d.pubchem_cid, d.inchikey, d.chebi_id
FROM amr.drug d
JOIN amr.drug_class_member dm ON dm.canonical_drug = d.canonical_name
WHERE dm.canonical_class = 'aminoglycoside antibiotic'
ORDER BY d.canonical_name;

-- Drug classes covered per source (from raw gene metadata)
SELECT source, drug_class, count(*) AS n_sequences
FROM amr.sequence_metadata
WHERE drug_class IS NOT NULL
GROUP BY source, drug_class
ORDER BY source, n_sequences DESC;

-- Canonical metadata for a sequence across all sources
SELECT * FROM amr.sequence_metadata WHERE jrc_id = 'JRC8214f7c788';
```
