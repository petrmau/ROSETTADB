-- =============================================================================
-- ROSETTADB — useful queries
-- =============================================================================


-- -----------------------------------------------------------------------------
-- Overview / coverage
-- -----------------------------------------------------------------------------

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

-- Drug classes covered per source (from raw gene metadata)
SELECT source, drug_class, count(*) AS n_sequences
FROM amr.sequence_metadata
WHERE drug_class IS NOT NULL
GROUP BY source, drug_class
ORDER BY source, n_sequences DESC;


-- -----------------------------------------------------------------------------
-- Sequence lookups
-- -----------------------------------------------------------------------------

-- All drug classes linked to a specific sequence (primary bridge table)
SELECT sdc.canonical_class, sdc.evidence_sources
FROM amr.sequence_drug_class sdc
WHERE sdc.jrc_id = 'JRC8214f7c788'
ORDER BY sdc.canonical_class;

-- All drug names linked to a specific sequence
SELECT sd.canonical_drug, sd.evidence_sources
FROM amr.sequence_drug sd
WHERE sd.jrc_id = 'JRC8214f7c788';

-- Canonical metadata for a sequence across all sources
SELECT * FROM amr.sequence_metadata WHERE jrc_id = 'JRC8214f7c788';

-- Sequences shared between two or more sources
SELECT s.jrc_id, count(DISTINCT g.source) AS n_sources,
       string_agg(DISTINCT g.source, '|' ORDER BY g.source) AS sources
FROM amr.sequence s
JOIN amr.gene g USING (jrc_id)
GROUP BY s.jrc_id
HAVING count(DISTINCT g.source) > 1
LIMIT 10;


-- -----------------------------------------------------------------------------
-- Drug class queries
-- -----------------------------------------------------------------------------

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

-- All gene sequences conferring resistance to fluoroquinolones (via CARD)
SELECT cgc.gene_name, cgc.card_short_name, cgc.gene_family, cgc.resistance_mechanism
FROM amr.card_gene_class cgc
WHERE cgc.canonical_class = 'fluoroquinolone antibiotic'
ORDER BY cgc.gene_family, cgc.gene_name;


-- -----------------------------------------------------------------------------
-- Drug / vocabulary queries
-- -----------------------------------------------------------------------------

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

-- All drugs in an ATC group with identifiers
SELECT d.canonical_name, d.atc_code, d.atc_group1, d.atc_group2,
       d.pubchem_cid, d.inchikey, d.chebi_id, d.loinc_codes
FROM amr.drug d
WHERE d.atc_group1 = 'Aminoglycoside antibacterials'
ORDER BY d.atc_group2, d.canonical_name;
