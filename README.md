# ROSETTADB

PostgreSQL database for Antimicrobial Resistance (AMR) genes aggregated from
multiple public sources, clustered at 100% nucleotide sequence identity, with a
second layer of protein-level clustering, and linked to a harmonised
drug/antibiotic class vocabulary.

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
8. [protein_cluster.py flags](#protein_clusterpy-flags)
9. [Useful queries](#useful-queries)

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
amr.protein            – one row per unique translated protein sequence
                         (protein_id JRCPRO_XXXXXX, MD5, AA sequence, length,
                          representative_jrc)
amr.sequence           – one row per unique nucleotide sequence
                         (jrc_id PK, MD5, raw sequence, length,
                          protein_id FK → amr.protein; NULL for non-CDS)
amr.cluster            – one row per 100%-identity nucleotide cluster (JRCGRP_XXXXXX);
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
amr.drug               – 585 canonical drugs (INN names); enriched with
                         PubChem CID, InChIKey, pipe-delimited ATC leaf code(s),
                         ChEBI ID, LOINC susceptibility test codes, ATC group labels
                         (416 of these sourced from CARD ARO OBO)
amr.drug_alias         – 2 662 spelling variants, brand names, synonyms,
                         abbreviations, and source-specific aliases
amr.drug_class_member  – 621 drug → class memberships (many-to-many);
                         provenance tracked (resfinder | ncbi | curated | aro_obo)
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

### Nucleotide clusters

Sequences sharing 100% nucleotide identity (same MD5) are grouped into one
cluster. The **representative** is the member with the longest `original_header`.
Each cluster gets a stable identifier of the form `JRCGRP_000001`.

### Protein clusters

A second clustering layer groups nucleotide sequences that encode **exactly the
same protein**, despite synonymous codon differences or start/stop variation.
Run `protein_cluster.py` after `ingest.py` to populate `amr.protein` and
back-fill `amr.sequence.protein_id`.

Translation rules applied before comparing:

| Edge case | Rule |
|-----------|------|
| Alternative start codon (GTG, TTG, CTG, ATT, ATC, ATA) | Forced to Met (M) — bacteria load fMet regardless of the start codon spelling |
| Trailing stop codon (TAA, TAG, TGA) present | Stripped before hashing — so sequences with and without stop codon map to the same protein |
| Trailing stop codon absent | No action required |
| Length not divisible by 3, ambiguous nucleotides (N…), or internal stop | `protein_id = NULL` — sequence is non-CDS, truncated, or frameshifted |

The protein representative is the **longest nucleotide sequence** encoding that
protein (ties broken lexicographically by jrc_id). Each protein cluster gets a
stable identifier of the form `JRCPRO_000001`.

---

## Harmonised drug vocabulary

The `harmonise/` directory contains pipeline scripts and the generated TSVs
loaded into the database:

| File | Rows | Description |
|------|------|-------------|
| `class_mapping.tsv` | 64 | Canonical drug classes with ARO accessions, ResFinder/NCBI aliases — **generated by `build_class_mapping.py`** |
| `drug_class_direct.tsv` | 33 | Curated drug → class overrides for entries not resolvable by source-name lookup — **generated by `build_drug_class_direct.py`** |
| `drug_canonical.tsv` | 585 | Canonical drugs (INN) with PubChem/ATC/ChEBI/LOINC identifiers; `sources` column records which pipeline(s) introduced each drug (`card`, `ncbi`, `resfinder`, `CARD.obo`) |
| `drug_alias.tsv` | 2 662 | Source aliases, synonyms, and abbreviations (AMR R-package enriched) |
| `drug_class_member.tsv` | 621 | Drug → class links with evidence provenance (`resfinder \| ncbi \| curated \| aro_obo`) |
| `gene_drug_link.tsv` | 10 388 | NCBI gene → drug/class pairings (raw + normalised) |
| `card_gene_class.tsv` | 13 361 | CARD ARO gene → canonical drug class links (from `aro_index.tsv`) |
| `aro_drug_class_member.tsv` | 528 | Drug → class links extracted from `sources/CARD/aro.obo` via `is_a` ancestry — **generated by `parse_aro_obo.py`**; merged into `drug_class_member.tsv` by `harmonise.py` |
| `aro_gene_class.tsv` | 803 | Gene → class links from `confers_resistance_to_drug_class` edges in `aro.obo` — **generated by `parse_aro_obo.py`** |
| `atc_codes_all.tsv` | 5 680 | Full WHO ATC/DDD index (all 14 categories); used by `enrich.py` for ATC lookup |
| `antimicrobials.txt` | 505 | Cached copy of the [msberends/AMR](https://github.com/msberends/AMR) antimicrobials reference; downloaded by `enrich_amr_r.py --download` |

| Script | Description |
|--------|-------------|
| `build_class_mapping.py` | Generates `class_mapping.tsv`; looks up ARO accessions live from `sources/CARD/card.json`; re-run when CARD is updated |
| `build_drug_class_direct.py` | Generates `drug_class_direct.tsv` from the curated override list in the script; no live lookups |
| `harmonise.py` | Parses all three AMR sources + CARD ARO OBO and produces `drug_canonical.tsv`, `drug_alias.tsv`, `drug_class_member.tsv`, `gene_drug_link.tsv`, `aro_drug_class_member.tsv`, `aro_gene_class.tsv` |
| `parse_aro_obo.py` | Parses `sources/CARD/aro.obo` to extract (1) drug → drug class `is_a` hierarchy → `aro_drug_class_member.tsv`; (2) `confers_resistance_to_drug_class` edges → `aro_gene_class.tsv`; called automatically by `harmonise.py` |
| `enrich_amr_r.py` | Pre-enriches `drug_canonical.tsv` (PubChem CID, ATC code, LOINC codes) and expands `drug_alias.tsv` with synonyms/abbreviations from the AMR R-package reference; **zero API calls**; run before `enrich.py` |
| `enrich.py` | Enriches `drug_canonical.tsv` with InChIKey and ChEBI ID via ChEBI/PubChem APIs (fills only what `enrich_amr_r.py` could not) |
| `enrich_pubchem.py` | Low-level PubChem PUG REST helper used by `enrich.py`; `lookup_atc()` extracts the most-specific (level-5) ATC code per classification tree via regex, returning all trees pipe-delimited |
| `fetch_atc_codes.py` | Rebuilds `atc_codes_all.tsv` by scraping atcddd.fhi.no (~10 min) |
| `parse_card_aro.py` | Rebuilds `card_gene_class.tsv` from `sources/CARD/aro_index.tsv` |

### Regenerate harmonised TSVs

Run these scripts in order whenever source data is updated. They write new TSVs to
`harmonise/` but do **not** touch the database — run `ingest.py --harmonise-only`
afterwards to reload.

```bash
# Step 0 — Rebuild curated reference tables (re-run when CARD is updated or overrides change)
python harmonise/build_class_mapping.py    # → class_mapping.tsv (ARO accessions from CARD JSON)
python harmonise/build_drug_class_direct.py # → drug_class_direct.tsv (curated drug → class overrides)

# Step 1 — (Re)build the local WHO ATC code table (only needed when the ATC index is outdated)
python harmonise/fetch_atc_codes.py all
# → writes harmonise/atc_codes_all.tsv  (~5 680 rows, takes ~10 min; skip if file is current)

# Step 2 — Rebuild drug/class tables from source data
python harmonise/harmonise.py

# Step 2.5 — Pre-enrich from the AMR R-package reference (offline; zero API calls)
#            Fills pubchem_cid, atc_code, loinc_codes; expands drug_alias.tsv
#            Pass --download to refresh the cached antimicrobials.txt
python harmonise/enrich_amr_r.py

# Step 3 — Enrich drug_canonical.tsv with InChIKey and ChEBI ID via live APIs
#          (runs faster now because many CIDs/ATCs are already filled)
python harmonise/enrich.py

# Step 4 — Rebuild CARD gene → class links from sources/CARD/aro_index.tsv
python harmonise/parse_card_aro.py
```

#### Enrichment strategy (`enrich_amr_r.py` + `enrich.py`)

`enrich_amr_r.py` runs first and fills identifiers from the offline AMR R-package
reference (`antimicrobials.txt`) — no network calls:

| Field | Source column | Behaviour |
|-------|--------------|-----------|
| `pubchem_cid` | `cid` | Filled only if currently blank; not overwritten |
| `atc_code` | `atc` | All codes stored, sorted: J first, Q immediately after J, then other categories. Always updated when the reference has data (single-code values from a prior `enrich.py` run become full lists) |
| `loinc_codes` | `loinc` | Comma-separated susceptibility test identifiers; always updated when present |
| `atc_group1` | `atc_group1` | ATC level-2 group (e.g. `"Aminoglycoside antibacterials"`); always updated when present |
| `atc_group2` | `atc_group2` | ATC level-3 group (e.g. `"Other aminoglycosides"`); always updated when present |
| `drug_alias.tsv` | `name`, `synonyms`, `abbreviations` | Adds `synonym` and `abbreviation` rows with `source=amr_r` |

**Combination drug matching:** canonical names using `+` or ` & ` separators
(e.g. `piperacillin+tazobactam`) are matched against the reference (which uses
`/`) by comparing the sorted frozenset of normalised component names —
order-independent. This gives combination drugs their CID and full ATC code pair
(e.g. `J01CR05,QJ01CR05`).

`enrich.py` runs second and resolves the remaining identifiers via live APIs:

| Step | Source | Fields retrieved |
|------|--------|-----------------|
| 1 | ChEBI by INN name | `chebi_id`, `inchikey` — picks highest-star exact-name match |
| 2a | PubChem by **InChIKey** (if ChEBI hit) | `pubchem_cid` — unambiguous |
| 2b | PubChem by **name** (fallback) | `pubchem_cid`, `inchikey` |
| 3 | **Local ATC table** (`atc_codes_all.tsv`) | `atc_code` — no network call; J category wins when a name appears in multiple categories |
| 4 | PubChem classification (fallback) | `atc_code` — only if local table misses |
| 5 | PubChem xrefs (last resort) | `chebi_id` if step 1 failed but CID was found |

ChEBI is used first because it is fully curated (3-star entries are manually reviewed),
and an InChIKey-based PubChem lookup avoids the salt/stereoisomer ambiguity of
name-based searches. ATC codes are resolved from the local WHO ATC table (step 3)
rather than PubChem, which gives better coverage and is instantaneous.

**ATC disambiguation:** when a drug name appears in more than one ATC category, the
category is chosen by this priority order:
`J > P > D > A > L > B > C > G > H > M > N > R > S > V`
(J = antiinfectives for systemic use always wins for AMR drugs).

Cache: `harmonise/.enrich_cache.json` — all 169 current drugs are pre-cached; re-runs
make zero live API calls unless new drugs are added or cache entries are cleared.
On each run, cached entries with an empty `atc_code` are automatically back-filled from
the local ATC table at no network cost.

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

**4. Run protein clustering**

```bash
python protein_cluster.py \
  --dsn "host=localhost dbname=rosettadb user=postgres password=postgres"
```

This translates every nucleotide sequence, groups sequences encoding the same
protein, populates `amr.protein` (one row per unique protein, `JRCPRO_XXXXXX`),
and sets `amr.sequence.protein_id`. Also idempotent — only processes rows where
`protein_id IS NULL`.

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

After re-ingesting, run protein clustering to process any new sequences:

```bash
python protein_cluster.py \
  --dsn "host=localhost dbname=rosettadb user=postgres password=postgres"
```

---

### Update drug vocabulary and re-create links

Use this after editing any `harmonise/` TSV (e.g. adding a new drug class,
correcting a mapping, adding an alias). Sequences do **not** need to be
re-ingested.

**Step 1 — regenerate TSVs** (only if the source scripts were modified):

```bash
# Rebuild curated reference tables (if CARD updated or overrides changed)
python harmonise/build_class_mapping.py     # class_mapping.tsv
python harmonise/build_drug_class_direct.py # drug_class_direct.tsv

# Rebuild ATC table only if atcddd.fhi.no has been updated (slow — ~10 min)
python harmonise/fetch_atc_codes.py all

python harmonise/harmonise.py      # drug classes + drugs + aliases + class membership
python harmonise/enrich_amr_r.py   # offline pre-enrichment (CID/ATC/LOINC + aliases)
python harmonise/enrich.py         # ChEBI/InChIKey via live APIs
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

## protein_cluster.py flags

| Flag | Default | Description |
|------|---------|-------------|
| `--dsn` | `$ROSETTADB_DSN` | PostgreSQL connection string |
| `--dry-run` | — | Translate and group sequences but do not write to DB; useful for previewing how many proteins would be created |

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

-- Protein clustering summary
SELECT
    count(*)                                          AS total_proteins,
    sum(nucleotide_variants)                          AS sequences_with_protein,
    sum(nucleotide_variants) FILTER
        (WHERE nucleotide_variants > 1)               AS sequences_in_synonymous_group,
    max(nucleotide_variants)                          AS max_synonymous_variants
FROM (
    SELECT p.protein_id, count(s.jrc_id) AS nucleotide_variants
    FROM amr.protein p
    JOIN amr.sequence s USING (protein_id)
    GROUP BY p.protein_id
) sub;

-- All nucleotide sequences encoding the same protein as a given jrc_id
SELECT s.jrc_id, s.sequence_length, p.protein_id, p.protein_length
FROM amr.sequence s
JOIN amr.protein p ON p.protein_id = (
    SELECT protein_id FROM amr.sequence WHERE jrc_id = 'JRC8214f7c788'
)
WHERE s.protein_id = p.protein_id
ORDER BY s.sequence_length DESC;

-- Sequences with no protein assignment (non-CDS, truncated, or frameshifted)
SELECT count(*) FROM amr.sequence WHERE protein_id IS NULL;
```
