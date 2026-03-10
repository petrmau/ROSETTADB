-- ROSETTADB: AMR Gene Database Schema
-- Sources: RESFINDER, CARD, NCBI AMRFinderPlus

CREATE SCHEMA IF NOT EXISTS amr;

-- ─────────────────────────────────────────────────────────
-- Core sequence table (one row per unique nucleotide sequence)
-- ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS amr.sequence (
    jrc_id          VARCHAR(13) PRIMARY KEY,  -- JRC + first 10 chars of md5
    sequence_md5    CHAR(32)    NOT NULL UNIQUE,
    sequence        TEXT        NOT NULL,
    sequence_length INTEGER     NOT NULL GENERATED ALWAYS AS (length(sequence)) STORED
);

-- ─────────────────────────────────────────────────────────
-- Cluster table (group of 100%-identical sequences)
-- ─────────────────────────────────────────────────────────
CREATE SEQUENCE IF NOT EXISTS amr.cluster_id_seq START 1;

CREATE TABLE IF NOT EXISTS amr.cluster (
    cluster_id          VARCHAR(15) PRIMARY KEY
                            DEFAULT 'JRCGRP_' || lpad(nextval('amr.cluster_id_seq')::text, 6, '0'),
    representative_jrc  VARCHAR(13) NOT NULL UNIQUE REFERENCES amr.sequence(jrc_id)
);

-- ─────────────────────────────────────────────────────────
-- Gene entry: one row per FASTA record (before dedup)
-- ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS amr.gene (
    gene_id         SERIAL PRIMARY KEY,
    jrc_id          VARCHAR(13) NOT NULL REFERENCES amr.sequence(jrc_id),
    cluster_id      VARCHAR(15) NOT NULL REFERENCES amr.cluster(cluster_id),
    source          VARCHAR(20) NOT NULL CHECK (source IN ('RESFINDER', 'CARD', 'NCBI')),
    original_header TEXT        NOT NULL,
    source_acc      TEXT,               -- primary accession parsed from header
    gene_name       TEXT,
    product_name    TEXT,
    drug_class      TEXT,
    resistance_mechanism TEXT,
    gene_family     TEXT,
    aro_accession   TEXT,               -- CARD-specific
    refseq_protein  TEXT,
    refseq_nucleotide TEXT,
    genbank_protein  TEXT,
    genbank_nucleotide TEXT,
    amr_class       TEXT,
    amr_subclass    TEXT,
    scope           TEXT,
    -- NCBI additional fields
    allele          TEXT,               -- NCBI allele designation (e.g. aac(3)-VIIIa)
    ncbi_type       TEXT,               -- NCBI type field (e.g. AMR)
    ncbi_subtype    TEXT,               -- NCBI subtype field
    -- ResFinder additional fields
    pmid            TEXT,               -- PubMed ID(s) for resistance reference
    notes           TEXT,               -- ResFinder notes / alternative names
    required_gene   TEXT,               -- ResFinder required co-gene
    -- CARD additional fields
    card_short_name TEXT,               -- CARD short name (15-char abbreviation)
    card_cvterm_id  TEXT,               -- CARD CVTERM ID
    card_model_id   TEXT,               -- CARD Model ID
    card_model_sequence_id TEXT,        -- CARD Model Sequence ID
    loaded_at       TIMESTAMPTZ NOT NULL DEFAULT now()
);

-- ─────────────────────────────────────────────────────────
-- Canonical metadata per unique sequence per source
-- (one row per jrc_id × source after deduplication)
-- ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS amr.sequence_metadata (
    jrc_id               VARCHAR(13) NOT NULL REFERENCES amr.sequence(jrc_id),
    source               VARCHAR(20) NOT NULL CHECK (source IN ('RESFINDER','CARD','NCBI')),
    gene_name            TEXT,
    product_name         TEXT,
    drug_class           TEXT,
    resistance_mechanism TEXT,
    gene_family          TEXT,
    aro_accession        TEXT,
    refseq_protein       TEXT,
    refseq_nucleotide    TEXT,
    genbank_protein      TEXT,
    genbank_nucleotide   TEXT,
    amr_class            TEXT,
    amr_subclass         TEXT,
    scope                TEXT,
    allele               TEXT,
    ncbi_type            TEXT,
    ncbi_subtype         TEXT,
    pmid                 TEXT,
    notes                TEXT,
    required_gene        TEXT,
    card_short_name      TEXT,
    card_cvterm_id       TEXT,
    card_model_id        TEXT,
    card_model_sequence_id TEXT,
    PRIMARY KEY (jrc_id, source)
);

-- ─────────────────────────────────────────────────────────
-- Indexes
-- ─────────────────────────────────────────────────────────
CREATE INDEX IF NOT EXISTS idx_gene_jrc      ON amr.gene(jrc_id);
CREATE INDEX IF NOT EXISTS idx_gene_source   ON amr.gene(source);
CREATE INDEX IF NOT EXISTS idx_gene_cluster  ON amr.gene(cluster_id);
CREATE INDEX IF NOT EXISTS idx_gene_name     ON amr.gene(gene_name);
CREATE INDEX IF NOT EXISTS idx_seq_md5       ON amr.sequence(sequence_md5);
CREATE INDEX IF NOT EXISTS idx_seqmeta_jrc   ON amr.sequence_metadata(jrc_id);
CREATE INDEX IF NOT EXISTS idx_seqmeta_src   ON amr.sequence_metadata(source);

-- ─────────────────────────────────────────────────────────
-- Harmonised drug vocabulary
-- ─────────────────────────────────────────────────────────

-- Canonical antibiotic/antimicrobial drug classes
-- (names follow CARD ARO; ARO accessions are stable cross-source keys)
CREATE TABLE IF NOT EXISTS amr.drug_class (
    canonical_name  TEXT PRIMARY KEY,    -- e.g. "fluoroquinolone antibiotic"
    aro_accession   TEXT,               -- e.g. "ARO:3000150"
    category        TEXT,               -- antibiotic | biocide | inhibitor | resistance_mechanism | non_therapeutic | veterinary_ionophore
    resfinder_alias TEXT,               -- raw class string used in ResFinder
    ncbi_alias      TEXT,               -- slash-delimited tokens used in NCBI class field
    card_abbrev     TEXT,               -- CARD drug-level abbreviation (if any)
    card_class_abbrev TEXT,             -- CARD class-level abbreviation (e.g. AMG, BLA, FLO)
    notes           TEXT
);

-- Canonical drug / antimicrobial molecule table (INN names)
CREATE TABLE IF NOT EXISTS amr.drug (
    canonical_name  TEXT PRIMARY KEY,   -- INN name, lower-case (e.g. "amikacin")
    is_combination  BOOLEAN NOT NULL DEFAULT FALSE,
    components      TEXT,               -- pipe-delimited component names for combinations
    context         TEXT NOT NULL DEFAULT 'clinical',
                                        -- clinical | antituberculosis | veterinary |
                                        -- inhibitor | biocide | research_tool | non_therapeutic
    card_abbrev     TEXT,               -- AAC abbreviation (e.g. "AMK")
    atc_code        TEXT,               -- WHO ATC code (e.g. "J01GB06")
    inchikey        TEXT,               -- standard InChIKey (27 chars)
    pubchem_cid     TEXT,               -- PubChem canonical CID
    chebi_id        TEXT,               -- ChEBI identifier (e.g. "CHEBI:2637")
    loinc_codes     TEXT,               -- comma-separated LOINC susceptibility test codes
    atc_group1      TEXT,               -- ATC level-2 group (e.g. "Aminoglycoside antibacterials")
    atc_group2      TEXT,               -- ATC level-3 group (e.g. "Other aminoglycosides")
    sources         TEXT                -- pipe-delimited: card|ncbi|resfinder
);

CREATE INDEX IF NOT EXISTS idx_drug_inchikey ON amr.drug(inchikey)
    WHERE inchikey IS NOT NULL AND inchikey <> '';
CREATE INDEX IF NOT EXISTS idx_drug_atc      ON amr.drug(atc_code)
    WHERE atc_code IS NOT NULL AND atc_code <> '';
CREATE INDEX IF NOT EXISTS idx_drug_pubchem  ON amr.drug(pubchem_cid)
    WHERE pubchem_cid IS NOT NULL AND pubchem_cid <> '';

-- Drug spelling variants, brand names, source-specific aliases
CREATE TABLE IF NOT EXISTS amr.drug_alias (
    canonical_name  TEXT NOT NULL REFERENCES amr.drug(canonical_name),
    alias           TEXT NOT NULL,
    alias_type      TEXT NOT NULL,  -- source_name | card_abbreviation
    source          TEXT NOT NULL,  -- resfinder | card | ncbi
    PRIMARY KEY (canonical_name, alias, source)
);

CREATE INDEX IF NOT EXISTS idx_drug_alias_alias ON amr.drug_alias(lower(alias));

-- Many-to-many: drug → drug class membership
-- A single drug can belong to multiple classes (e.g. pristinamycin is both
-- streptogramin A and streptogramin B; broad-spectrum agents span classes).
CREATE TABLE IF NOT EXISTS amr.drug_class_member (
    canonical_drug   TEXT NOT NULL REFERENCES amr.drug(canonical_name),
    canonical_class  TEXT NOT NULL REFERENCES amr.drug_class(canonical_name),
    aro_accession    TEXT,   -- class ARO accession (denormalised for fast lookup)
    category         TEXT,   -- class category (denormalised)
    evidence_source  TEXT NOT NULL DEFAULT 'curated',
                            -- resfinder | ncbi | curated | resfinder_combination_component
    PRIMARY KEY (canonical_drug, canonical_class)
);

CREATE INDEX IF NOT EXISTS idx_dcm_class ON amr.drug_class_member(canonical_class);

-- Gene → drug/class links derived from NCBI AMRFinderPlus metadata
-- Each row represents one gene × one class/drug pairing as reported by NCBI.
-- Retained at full cardinality (10k+ rows) for traceability.
CREATE TABLE IF NOT EXISTS amr.gene_drug_link (
    gene_name               TEXT NOT NULL,
    accession               TEXT,
    element_type            TEXT,
    ncbi_class_raw          TEXT,   -- original NCBI "class" field value
    ncbi_subclass_raw       TEXT,   -- original NCBI "subclass" field value
    canonical_class_token   TEXT,   -- normalised class token (may be empty)
    canonical_drug_token    TEXT    -- normalised drug token (may be empty)
);

CREATE INDEX IF NOT EXISTS idx_gdl_gene     ON amr.gene_drug_link(gene_name);
CREATE INDEX IF NOT EXISTS idx_gdl_drug     ON amr.gene_drug_link(canonical_drug_token)
    WHERE canonical_drug_token IS NOT NULL AND canonical_drug_token <> '';
CREATE INDEX IF NOT EXISTS idx_gdl_class    ON amr.gene_drug_link(canonical_class_token)
    WHERE canonical_class_token IS NOT NULL AND canonical_class_token <> '';

-- ─────────────────────────────────────────────────────────
-- CARD gene → canonical drug class links (from ARO graph)
-- One row per ARO gene model × canonical drug class.
-- Source: CARD aro_index.tsv, mapped through class_mapping.tsv.
-- ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS amr.card_gene_class (
    aro_accession        TEXT NOT NULL,  -- ARO accession of the AMR gene model
    gene_name            TEXT NOT NULL,  -- full ARO Name (e.g. CblA-1)
    card_short_name      TEXT,           -- CARD Short Name / abbreviation
    gene_family          TEXT,           -- AMR Gene Family label
    resistance_mechanism TEXT,           -- e.g. antibiotic inactivation
    drug_class_raw       TEXT NOT NULL,  -- raw CARD Drug Class string
    canonical_class      TEXT NOT NULL
        REFERENCES amr.drug_class(canonical_name),
    PRIMARY KEY (aro_accession, canonical_class)
);

CREATE INDEX IF NOT EXISTS idx_cgc_aro     ON amr.card_gene_class(aro_accession);
CREATE INDEX IF NOT EXISTS idx_cgc_class   ON amr.card_gene_class(canonical_class);
CREATE INDEX IF NOT EXISTS idx_cgc_family  ON amr.card_gene_class(gene_family)
    WHERE gene_family IS NOT NULL AND gene_family <> '';

-- ─────────────────────────────────────────────────────────
-- Sequence → harmonised drug class (non-redundant link table)
-- One row per unique sequence × canonical drug class.
-- Aggregates evidence from all three source paths:
--   CARD      : gene.aro_accession → card_gene_class → drug_class
--   NCBI      : gene.gene_name     → gene_drug_link  → drug_class
--   RESFINDER : gene.drug_class text → drug_class.resfinder_alias
-- ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS amr.sequence_drug_class (
    jrc_id           VARCHAR(13) NOT NULL REFERENCES amr.sequence(jrc_id),
    canonical_class  TEXT        NOT NULL REFERENCES amr.drug_class(canonical_name),
    evidence_sources TEXT        NOT NULL,  -- pipe-delimited: CARD|NCBI|RESFINDER
    PRIMARY KEY (jrc_id, canonical_class)
);

CREATE INDEX IF NOT EXISTS idx_sdc_class ON amr.sequence_drug_class(canonical_class);
CREATE INDEX IF NOT EXISTS idx_sdc_jrcid ON amr.sequence_drug_class(jrc_id);

-- ─────────────────────────────────────────────────────────
-- Sequence → harmonised drug name (non-redundant link table)
-- One row per unique sequence × canonical drug.
-- Evidence path:
--   NCBI  : gene.gene_name → gene_drug_link.canonical_drug_token → drug
-- ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS amr.sequence_drug (
    jrc_id           VARCHAR(13) NOT NULL REFERENCES amr.sequence(jrc_id),
    canonical_drug   TEXT        NOT NULL REFERENCES amr.drug(canonical_name),
    evidence_sources TEXT        NOT NULL,  -- pipe-delimited: NCBI|CARD|RESFINDER
    PRIMARY KEY (jrc_id, canonical_drug)
);

CREATE INDEX IF NOT EXISTS idx_sd_drug  ON amr.sequence_drug(canonical_drug);
CREATE INDEX IF NOT EXISTS idx_sd_jrcid ON amr.sequence_drug(jrc_id);
