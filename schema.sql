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
CREATE TABLE IF NOT EXISTS amr.cluster (
    cluster_id          SERIAL PRIMARY KEY,
    representative_jrc  VARCHAR(13) NOT NULL REFERENCES amr.sequence(jrc_id)
);

-- ─────────────────────────────────────────────────────────
-- Gene entry: one row per FASTA record (before dedup)
-- ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS amr.gene (
    gene_id         SERIAL PRIMARY KEY,
    jrc_id          VARCHAR(13) NOT NULL REFERENCES amr.sequence(jrc_id),
    cluster_id      INTEGER     NOT NULL REFERENCES amr.cluster(cluster_id),
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
    loaded_at       TIMESTAMPTZ NOT NULL DEFAULT now()
);

-- ─────────────────────────────────────────────────────────
-- Indexes
-- ─────────────────────────────────────────────────────────
CREATE INDEX IF NOT EXISTS idx_gene_jrc      ON amr.gene(jrc_id);
CREATE INDEX IF NOT EXISTS idx_gene_source   ON amr.gene(source);
CREATE INDEX IF NOT EXISTS idx_gene_cluster  ON amr.gene(cluster_id);
CREATE INDEX IF NOT EXISTS idx_gene_name     ON amr.gene(gene_name);
CREATE INDEX IF NOT EXISTS idx_seq_md5       ON amr.sequence(sequence_md5);
