# ROSETTADB

PostgreSQL database for Antimicrobial Resistance (AMR) genes aggregated from
multiple public sources, clustered at 100% sequence identity, and linked to a
harmonised drug/antibiotic class vocabulary.

## Sources

| Source | Sequences | Metadata |
|--------|-----------|----------|
| CARD (v4.0.1) | ~6 052 nucleotide sequences | ARO accession, gene family, drug class, resistance mechanism |
| NCBI AMRFinderPlus | ~8 302 sequences | Gene family, class, subclass, scope, type |
| ResFinder | ~5 000+ sequences | Drug class, phenotype, PMID, required co-gene |

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
amr.drug_class         – 65 canonical antibiotic/antimicrobial classes
                         (names follow CARD ARO; ARO accessions are stable keys)
amr.drug               – 169 canonical drugs (INN names); enriched with
                         PubChem CID, InChIKey, ATC code, ChEBI ID
amr.drug_alias         – 80 spelling variants, brand names, source-specific aliases
amr.drug_class_member  – 172 drug → class memberships (many-to-many);
                         provenance tracked (resfinder | ncbi | curated)
```

### Gene → drug class links

```
amr.gene_drug_link     – 10 388 NCBI AMRFinderPlus gene × drug/class pairings
                         (raw class and subclass tokens + normalised canonical tokens)
amr.card_gene_class    – 13 377 CARD ARO gene × canonical drug class links
                         (6 440 unique gene models × 47 drug classes)
```

## Sequence identifier

Every unique nucleotide sequence receives a stable JRC identifier:

```
JRC<first-10-chars-of-MD5(sequence)>
```

e.g. `JRC8214f7c788`

## Clustering

Sequences sharing 100% identity (same MD5) are grouped into one cluster.
The **representative** is the member with the longest `original_header`.
Each cluster gets a stable identifier of the form `JRCGRP_000001`.

## Harmonised drug vocabulary

The `harmonise/` directory contains all pipeline scripts and generated TSVs:

| File | Description |
|------|-------------|
| `class_mapping.tsv` | 65 canonical drug classes with ARO accessions, ResFinder/NCBI aliases |
| `drug_canonical.tsv` | 169 canonical drugs (INN) with PubChem/ATC/ChEBI identifiers |
| `drug_alias.tsv` | 80 source-specific aliases and abbreviations |
| `drug_class_member.tsv` | 172 drug → class links with evidence provenance |
| `gene_drug_link.tsv` | 10 388 NCBI gene → drug/class pairings (raw + normalised) |
| `card_gene_class.tsv` | 13 377 CARD ARO gene → canonical drug class links |

### Regenerate harmonised TSVs

```bash
# Rebuild drug/class tables from source data
python harmonise/harmonise.py

# Re-enrich PubChem identifiers (uses .pubchem_cache.json; no live API calls if cache is warm)
python harmonise/enrich_pubchem.py

# Rebuild CARD gene → class links from aro_index.tsv
python harmonise/parse_card_aro.py
```

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

The script applies `schema.sql`, loads all three FASTA sources, and then loads
the harmonised drug vocabulary tables. It is fully idempotent (`ON CONFLICT DO UPDATE`).

Optional flags:

| Flag | Default | Description |
|------|---------|-------------|
| `--resfinder` | `sources/resfinder_db/all.fsa` | ResFinder FASTA path |
| `--card` | `sources/CARD/nucleotide_fasta_protein_homolog_model.fasta` | CARD FASTA path |
| `--ncbi` | `sources/amr_finder_plus/ncbi_dataset/data/nucleotide.fna` | NCBI FASTA path |
| `--schema` | `schema.sql` | DDL file path |
| `--skip-harmonise` | — | Skip loading harmonised drug vocabulary tables |

## Useful queries

```sql
-- Sequences per source
SELECT source, count(*) FROM amr.gene GROUP BY source;

-- Sequences shared between CARD and NCBI
SELECT s.jrc_id, count(DISTINCT g.source) AS n_sources
FROM amr.sequence s
JOIN amr.gene g USING (jrc_id)
GROUP BY s.jrc_id
HAVING count(DISTINCT g.source) > 1
LIMIT 10;

-- All genes conferring resistance to fluoroquinolones (via CARD ARO graph)
SELECT cgc.gene_name, cgc.card_short_name, cgc.gene_family, cgc.resistance_mechanism
FROM amr.card_gene_class cgc
WHERE cgc.canonical_class = 'fluoroquinolone antibiotic'
ORDER BY cgc.gene_family, cgc.gene_name;

-- Canonical drug class for a given drug name
SELECT dm.canonical_drug, dc.canonical_name AS drug_class, dc.aro_accession, dm.evidence_source
FROM amr.drug_class_member dm
JOIN amr.drug_class dc ON dc.canonical_name = dm.canonical_class
WHERE dm.canonical_drug = 'ciprofloxacin';

-- Drugs in a class with PubChem/ATC identifiers
SELECT d.canonical_name, d.atc_code, d.pubchem_cid, d.inchikey, d.chebi_id
FROM amr.drug d
JOIN amr.drug_class_member dm ON dm.canonical_drug = d.canonical_name
WHERE dm.canonical_class = 'aminoglycoside antibiotic'
ORDER BY d.canonical_name;

-- NCBI gene → drug links for a specific gene family
SELECT gdl.gene_name, gdl.ncbi_class_raw, gdl.ncbi_subclass_raw,
       gdl.canonical_class_token, gdl.canonical_drug_token
FROM amr.gene_drug_link gdl
WHERE gdl.gene_name LIKE 'aac%'
LIMIT 20;

-- Drug classes covered per source (from raw gene metadata)
SELECT source, drug_class, count(*) AS n_sequences
FROM amr.sequence_metadata
WHERE drug_class IS NOT NULL
GROUP BY source, drug_class
ORDER BY source, n_sequences DESC;

-- Canonical metadata for a sequence across all sources
SELECT * FROM amr.sequence_metadata WHERE jrc_id = 'JRC8214f7c788';
```
