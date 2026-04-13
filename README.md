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
7. [Export scripts](#export-scripts)
8. [ingest.py flags](#ingestpy-flags)
9. [protein_cluster.py flags](#protein_clusterpy-flags)
10. [Useful queries](#useful-queries)

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
                         (416 of these sourced from CARD ARO OBO);
                         `context` field distinguishes specific drugs from
                         class-level tokens: clinical | veterinary | antituberculosis |
                         inhibitor | biocide | research_tool | non_therapeutic |
                         drug_class_name (20 entries, e.g. "aminoglycoside", "carbapenem")
amr.drug_alias         – 2 662 spelling variants, brand names, synonyms,
                         abbreviations, and source-specific aliases
amr.drug_class_member  – 621 drug → class memberships (many-to-many);
                         provenance tracked (resfinder | ncbi | curated | aro_obo)
```

### Gene → drug class / drug links (intermediate)

```
amr.gene_drug_link     – NCBI AMRFinderPlus gene × drug/class pairings
                         (raw class and subclass tokens + normalised canonical tokens)
amr.card_gene_class    – CARD ARO gene × canonical drug class links
                         (unique ARO models × drug classes; from aro_index.tsv)
amr.aro_gene_drug      – CARD ARO gene × canonical drug links derived from
                         confers_resistance_to_antibiotic edges in aro.obo
                         (high-specificity complement to card_gene_class:
                          encodes exact drug targets, e.g. TEM-1 → ampicillin)
```

### Sequence → drug class / drug (final bridge tables)

```
amr.sequence_drug_class – one row per unique sequence × canonical drug class
                          Evidence paths:
                            CARD:      gene.aro_accession → card_gene_class → drug_class
                            NCBI:      gene.amr_class (slash-split) → drug_class.ncbi_alias
                                       + gene.gene_name → gene_drug_link → drug_class
                            ResFinder: gene.drug_class text → drug_class.resfinder_alias
amr.sequence_drug       – one row per unique sequence × canonical drug name
                          Evidence paths:
                            NCBI:  gene.gene_name → gene_drug_link → drug
                            CARD:  gene.aro_accession → aro_gene_drug → drug
                                   (covers confers_resistance_to_antibiotic edges)
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
| `drug_canonical.tsv` | — | Canonical drugs (INN) with PubChem/ATC/ChEBI/LOINC identifiers; `sources` column records which pipeline(s) introduced each drug (`card`, `ncbi`, `resfinder`, `CARD.obo`); `context` column flags 20 class-level tokens as `drug_class_name` |
| `drug_alias.tsv` | — | Source aliases, synonyms, and abbreviations (AMR R-package enriched) |
| `drug_class_member.tsv` | — | Drug → class links with evidence provenance (`resfinder \| ncbi \| curated \| aro_obo`) |
| `gene_drug_link.tsv` | — | NCBI gene → drug/class pairings (raw + normalised) |
| `card_gene_class.tsv` | — | CARD ARO gene → canonical drug class links (from `aro_index.tsv`) |
| `aro_drug_class_member.tsv` | — | Drug → class links extracted from `sources/CARD/aro.obo` via `is_a` ancestry; merged into `drug_class_member.tsv` by `harmonise.py` — **auto-generated, do not run manually** |
| `aro_gene_class.tsv` | — | Gene → class links from `confers_resistance_to_drug_class` edges in `aro.obo` — **auto-generated, do not run manually** |
| `aro_gene_drug.tsv` | — | Gene → drug links from `confers_resistance_to_antibiotic` edges in `aro.obo`; loaded into `amr.aro_gene_drug` and used to populate `amr.sequence_drug` with CARD evidence — **auto-generated, do not run manually** |
| `atc_codes_all.tsv` | 5 680 | Full WHO ATC/DDD index (all 14 human categories); used by `enrich.py` for ATC lookup |
| `atcvet_codes_all.tsv` | ~3 500 | Full WHO ATCvet index (all 15 veterinary categories, Q-prefix codes); merged with `atc_codes_all.tsv` at load time by `enrich.py` |
| `antimicrobials.txt` | 505 | Cached copy of the [msberends/AMR](https://github.com/msberends/AMR) antimicrobials reference; downloaded by `enrich_amr_r.py --download`; also consulted by `harmonise.py` at startup to build a comprehensive synonym→INN map used by `normalise_name()` |

| Script | Description |
|--------|-------------|
| `build_class_mapping.py` | Generates `class_mapping.tsv`; looks up ARO accessions live from `sources/CARD/card.json`; re-run when CARD is updated |
| `build_drug_class_direct.py` | Generates `drug_class_direct.tsv` from the curated override list in the script; no live lookups |
| `harmonise.py` | Parses all three AMR sources + CARD ARO OBO and produces `drug_canonical.tsv`, `drug_alias.tsv`, `drug_class_member.tsv`, `gene_drug_link.tsv`, `aro_drug_class_member.tsv`, `aro_gene_class.tsv`, `aro_gene_drug.tsv`; loads `antimicrobials.txt` at startup to build a synonym→INN map (see *INN normalisation* below); assigns `context` flags via `_context_flag()` — edit the `CLASS_TERM_DRUGS` set here to manage class-level token entries |
| `parse_aro_obo.py` | Parses `sources/CARD/aro.obo` to extract: (1) drug → drug class `is_a` hierarchy → `aro_drug_class_member.tsv`; (2) `confers_resistance_to_drug_class` edges → `aro_gene_class.tsv`; (3) `confers_resistance_to_antibiotic` edges → `aro_gene_drug.tsv` (high-specificity gene→drug links, e.g. TEM-1 → ampicillin). **Called automatically by `harmonise.py` — never run directly.** |
| `enrich_amr_r.py` | Pre-enriches `drug_canonical.tsv` (PubChem CID, ATC code, LOINC codes) and expands `drug_alias.tsv` with synonyms/abbreviations from the AMR R-package reference; **zero API calls**; run before `enrich_atc_groups.py` and `enrich.py` |
| `enrich_atc_groups.py` | Fills any remaining `atc_group1` / `atc_group2` gaps in `drug_canonical.tsv` by deriving them from the ATC code prefix (4-char → group1, 5-char → group2) using a lookup built from the already-enriched rows in `antimicrobials.txt`; **zero API calls**; run after `enrich_amr_r.py` |
| `enrich.py` | Enriches `drug_canonical.tsv` with InChIKey and ChEBI ID via ChEBI/PubChem APIs (fills only what `enrich_amr_r.py` could not) |
| `enrich_pubchem.py` | Low-level PubChem PUG REST helper used by `enrich.py`; `lookup_atc()` extracts the most-specific (level-5) ATC code per classification tree via regex, returning all trees pipe-delimited |
| `fetch_atc_codes.py` | Rebuilds `atc_codes_all.tsv` by scraping atcddd.fhi.no (~10 min) |
| `fetch_atcvet_codes.py` | Rebuilds `atcvet_codes_all.tsv` by scraping atcddd.fhi.no/atcvet (~10 min); ATCvet codes carry a `Q` prefix (e.g. `QJ01AA01`); `enrich.py` merges both tables automatically |
| `parse_card_aro.py` | Rebuilds `card_gene_class.tsv` from `sources/CARD/aro_index.tsv` |

### INN normalisation

All drug names produced by the pipeline are resolved to their **International
Nonproprietary Name (INN)** via `normalise_name()`, which applies three layers
in order:

| Layer | Source | Example |
|-------|--------|---------|
| 1 — UK/regional spellings | `UK_TO_INN` dict (in-code) | `amoxycillin` → `amoxicillin` |
| 2 — Known source overrides | `SOURCE_TO_INN` dict (in-code) | `rifampin` → `rifampicin`, typos in CARD |
| 3 — Comprehensive synonym table | `antimicrobials.txt` (`name`, `synonyms`, `abbreviations` columns) | `methicillin` → `meticillin`, brand names, lab codes |

The synonym table is loaded from `antimicrobials.txt` at the start of each
`harmonise.py` run.  Synonyms that would map ambiguously to two different INNs
are silently excluded to avoid incorrect redirections.

This ensures that the same compound is never stored twice under different
spellings, regardless of which source introduced it — `meticillin` is always
preferred over `methicillin`, etc.

### Regenerate harmonised TSVs

Run these scripts in order whenever source data is updated. They write new TSVs to
`harmonise/` but do **not** touch the database — run `ingest.py --harmonise-only`
afterwards to reload.

```bash
# Step 0 — Rebuild curated reference tables (re-run when CARD is updated or overrides change)
python harmonise/build_class_mapping.py    # → class_mapping.tsv (ARO accessions from CARD JSON)
python harmonise/build_drug_class_direct.py # → drug_class_direct.tsv (curated drug → class overrides)

# Step 1 — (Re)build the local WHO ATC code tables (only needed when the indices are outdated)
python harmonise/fetch_atc_codes.py all
# → writes harmonise/atc_codes_all.tsv    (~5 680 rows, takes ~10 min; skip if file is current)
python harmonise/fetch_atcvet_codes.py all
# → writes harmonise/atcvet_codes_all.tsv (~3 500 rows, takes ~10 min; skip if file is current)

# Step 2 — Rebuild drug/class tables from source data
python harmonise/harmonise.py

# Step 2.5 — Pre-enrich from the AMR R-package reference (offline; zero API calls)
#            Fills pubchem_cid, atc_code, loinc_codes; expands drug_alias.tsv
#            Pass --download to refresh the cached antimicrobials.txt
python harmonise/enrich_amr_r.py

# Step 2.6 — Fill remaining atc_group1 / atc_group2 gaps from ATC code prefix (offline)
#            Derives group labels from 4/5-char ATC prefix using a lookup built from
#            already-labelled rows in antimicrobials.txt; run after enrich_amr_r.py
#            Use --dry-run to preview changes without writing
python harmonise/enrich_atc_groups.py

# Step 3 — Enrich drug_canonical.tsv with InChIKey and ChEBI ID via live APIs
#          (runs faster now because many CIDs/ATCs are already filled)
python harmonise/enrich.py

# Step 4 — Rebuild CARD gene → class links from sources/CARD/aro_index.tsv
python harmonise/parse_card_aro.py
```

#### Enrichment strategy (`enrich_amr_r.py` → `enrich_atc_groups.py` → `enrich.py`)

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

`enrich_atc_groups.py` runs second and fills any `atc_group1` / `atc_group2` values
that `enrich_amr_r.py` could not set (drugs present in the database but absent from
or unlabelled in `antimicrobials.txt`). It builds a prefix lookup from the rows that
already have group labels:

| Prefix length | Maps to | Example |
|---------------|---------|---------|
| 4 chars (`J01D`) | `atc_group1` | → "Other beta-lactam antibacterials" |
| 5 chars (`J01DB`) | `atc_group2` | → "First-generation cephalosporins" |

For drugs with multiple ATC codes the first code in J→Q→P→… priority order that
yields a hit is used.  No network calls; requires only the local `antimicrobials.txt`
cache.  Supports `--dry-run` to preview changes.  On current data this fills
**53 additional `atc_group1`** and **52 additional `atc_group2`** values.

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

# Rebuild ATC tables only if atcddd.fhi.no has been updated (slow — ~10 min each)
python harmonise/fetch_atc_codes.py all
python harmonise/fetch_atcvet_codes.py all

python harmonise/harmonise.py         # drug classes + drugs + aliases + class membership
python harmonise/enrich_amr_r.py      # offline pre-enrichment (CID/ATC/LOINC + aliases)
python harmonise/enrich_atc_groups.py # offline: fill remaining atc_group1/2 from ATC prefix
python harmonise/enrich.py            # ChEBI/InChIKey via live APIs
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

## Export scripts

| Script | Output | Description |
|--------|--------|-------------|
| `export_gene_clusters.py` | TSV | One row per nucleotide cluster: `cluster_id`, representative sequence, protein_id, gene names per source |
| `export_protein_clusters.py` | TSV | One row per protein cluster: `protein_id`, protein sequence, representative nucleotide sequence |
| `export_cluster_drugs.py` | TSV | One row per (cluster × canonical drug): `cluster_id`, `canonical_drug`, `link_source` (pipe-delimited evidence sources), plus drug metadata columns (`atc_code`, `inchikey`, `pubchem_cid`, `chebi_id`, `sources`, `atc_group1`, `atc_group2`) |

All three scripts accept `--dsn` (or `$ROSETTADB_DSN`) and `--output` (default: stdout).

`export_cluster_drugs.py` additional flags:

| Flag | Description |
|------|-------------|
| `--require-inchikey` | Skip drugs with no InChIKey |
| `--exclude-class-terms` | Skip entries with `context='drug_class_name'` — class-level tokens such as `"aminoglycoside"` or `"third-generation cephalosporin"` that are not specific drug entities (20 entries) |
| `--direct-links-only` | Only output Path A links (`amr.sequence_drug`): NCBI gene→drug and CARD `confers_resistance_to_antibiotic` edges. Excludes Path B (class fan-out via `sequence_drug_class → drug_class_member`), which can introduce broad class-level associations |

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

-- Direct drug links for a CARD gene via confers_resistance_to_antibiotic
SELECT agd.gene_name, agd.canonical_drug, agd.drug_aro_accession
FROM amr.aro_gene_drug agd
WHERE agd.gene_name = 'TEM-1'
ORDER BY agd.canonical_drug;

-- Sequences sharing the same jrc_id but carrying different gene names across sources
-- (indicates either a same-sequence annotation by two databases, or a CARD sequence
--  quality issue where the wrong protein was associated with a gene model)
SELECT jrc_id,
       count(DISTINCT gene_name) AS n_names,
       string_agg(DISTINCT gene_name, ', ') AS names,
       string_agg(DISTINCT source, ', ')    AS sources
FROM amr.sequence_metadata
GROUP BY jrc_id
HAVING count(DISTINCT gene_name) > 1
ORDER BY n_names DESC;
```
