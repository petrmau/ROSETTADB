#!/usr/bin/env python3
"""
ROSETTADB ingestion script
Loads AMR gene sequences from RESFINDER, CARD and NCBI into PostgreSQL,
clusters at 100% sequence identity and selects representatives by
longest FASTA header.

Usage:
    python ingest.py [--dsn DSN] [--resfinder PATH] [--card PATH] [--ncbi PATH]

Defaults read from environment / constants below.
"""

import argparse
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import psycopg2
from psycopg2.extras import execute_values

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_DSN = "host=localhost dbname=rosettadb user=postgres password=postgres"

BASE = Path(__file__).parent / "sources"
DEFAULT_RESFINDER = BASE / "resfinder_db" / "all.fsa"
DEFAULT_CARD      = BASE / "CARD" / "nucleotide_fasta_protein_homolog_model.fasta"
DEFAULT_NCBI      = BASE / "amr_finder_plus" / "ncbi_dataset" / "data" / "nucleotide.fna"

NCBI_REPORT          = BASE / "amr_finder_plus" / "ncbi_dataset" / "data" / "data_report.jsonl"
CARD_ARO             = BASE / "CARD" / "aro_index.tsv"
RESFINDER_PHENOTYPES = BASE / "resfinder_db" / "phenotypes.txt"

HARMONISE = Path(__file__).parent / "harmonise"
DRUG_CLASS_TSV      = HARMONISE / "class_mapping.tsv"
DRUG_TSV            = HARMONISE / "drug_canonical.tsv"
DRUG_ALIAS_TSV      = HARMONISE / "drug_alias.tsv"
DRUG_CLASS_MEMBER_TSV = HARMONISE / "drug_class_member.tsv"
GENE_DRUG_LINK_TSV  = HARMONISE / "gene_drug_link.tsv"
CARD_GENE_CLASS_TSV = HARMONISE / "card_gene_class.tsv"


# ── Helpers ───────────────────────────────────────────────────────────────────

def jrc_id(seq: str) -> str:
    """Return JRC<first-10-chars-of-md5> identifier."""
    md5 = hashlib.md5(seq.encode()).hexdigest()
    return f"JRC{md5[:10]}"


def md5hex(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()


def parse_fasta(path: Path):
    """Yield (header, sequence) from a FASTA file. Handles multi-line seqs."""
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header = line[1:]
                parts = []
            else:
                parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts).upper()


# ── Metadata loaders ──────────────────────────────────────────────────────────

def load_ncbi_report(path: Path) -> dict:
    """Return dict keyed by RefSeq nucleotide accession."""
    meta = {}
    if not path.exists():
        return meta
    with open(path) as fh:
        for line in fh:
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            acc = (d.get("refseqNucleotide") or {}).get("accessionVersion", "")
            if acc:
                meta[acc.split(".")[0]] = d   # index without version
                meta[acc] = d                  # also with version
            # also index by genbank
            gb = (d.get("genbankNucleotide") or {}).get("accessionVersion", "")
            if gb:
                meta[gb.split(".")[0]] = d
                meta[gb] = d
    return meta


def load_card_aro(path: Path) -> dict:
    """Return dict keyed by DNA accession (e.g. AF002716.1)."""
    meta = {}
    if not path.exists():
        return meta
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if header is None:
                header = line.split("\t")
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                continue
            row = dict(zip(header, parts))
            dna = row.get("DNA Accession", "").strip()
            if dna:
                meta[dna] = row
                meta[dna.split(".")[0]] = row
    return meta


def load_resfinder_phenotypes(path: Path) -> dict:
    """
    Parse ResFinder phenotypes.txt.
    Returns dict keyed by Gene_accession (e.g. 'aac(6'')-Ib_2_M23634')
    with keys: drug_class, resistance_mechanism, phenotype.
    """
    meta = {}
    if not path.exists():
        return meta
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if header is None:
                header = line.split("\t")
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            row = dict(zip(header, parts))
            acc = row.get("Gene_accession no.", "").strip()
            if acc:
                meta[acc] = {
                    "drug_class":           row.get("Class", "").strip() or None,
                    "resistance_mechanism": row.get("Mechanism of resistance", "").strip() or None,
                    "phenotype":            row.get("Phenotype", "").strip() or None,
                    "pmid":                 row.get("PMID", "").strip() or None,
                    "notes":                row.get("Notes", "").strip() or None,
                    "required_gene":        row.get("Required_gene", "").strip() or None,
                }
    return meta


# ── Header parsers ────────────────────────────────────────────────────────────

_ACC_RE = re.compile(r"[A-Z]{1,2}_?\d{5,9}(\.\d+)?|[A-Z]{2}\d{6}(\.\d+)?")


def _first_acc(text: str) -> str:
    m = _ACC_RE.search(text)
    return m.group(0) if m else ""


def parse_card_header(header: str) -> dict:
    """
    Format: gb|GQ343019.1|+|132-1023|ARO:3002999|CblA-1 [organism]
    """
    parts = header.split("|")
    info: dict = {}
    if len(parts) >= 2:
        info["source_acc"] = parts[1].strip()
    if len(parts) >= 5:
        aro = parts[4].strip()
        if aro.startswith("ARO:"):
            info["aro_accession"] = aro
    if len(parts) >= 6:
        name_org = parts[5].strip()
        # gene name is before the space-bracket organism
        gene = re.split(r"\s+\[", name_org)[0].strip()
        info["gene_name"] = gene
    return info


def parse_ncbi_header(header: str) -> dict:
    """
    Format: NG_242157.1:101-637 Product name gene_name, description
    or:     NZ_CP012138.1:c1234-567 ...
    """
    info: dict = {}
    tok = header.split()
    if tok:
        acc_range = tok[0]
        acc = acc_range.split(":")[0]
        info["source_acc"] = acc
    return info


def parse_resfinder_header(header: str) -> dict:
    """
    ResFinder headers: <gene>_<allele_number>_<genbank_accession>
    e.g. ant(2'')-Ia_8_AY920928  -> allele=ant(2'')-Ia_8, source_acc=AY920928
         aac(6')-Ib_2_M23634     -> allele=aac(6')-Ib_2,  source_acc=M23634
    The accession is the last _-delimited token that looks like a GenBank accession.
    allele = all tokens before the accession joined with '_'.
    gene_name = first token (gene family, e.g. ant(2'')-Ia).
    genbank_nucleotide = source_acc (they are the same for ResFinder).
    """
    info: dict = {}
    parts = header.split("_")

    # Find the rightmost accession-like token
    acc_idx = None
    for i in range(len(parts) - 1, 0, -1):
        p = parts[i]
        # Match with or without version suffix (.1 etc.)
        if _ACC_RE.fullmatch(p) or _ACC_RE.fullmatch(p.split(".")[0] + ".1"):
            acc_idx = i
            break

    if acc_idx is not None:
        info["source_acc"] = parts[acc_idx]
        info["genbank_nucleotide"] = parts[acc_idx]
        info["allele"] = "_".join(parts[:acc_idx])
    else:
        # fallback: treat last token as accession
        info["source_acc"] = parts[-1]
        info["genbank_nucleotide"] = parts[-1]
        info["allele"] = "_".join(parts[:-1])

    if parts:
        info["gene_name"] = parts[0]
    return info


# ── Gene record builder ───────────────────────────────────────────────────────

def build_gene_records(source: str, fasta_path: Path,
                       ncbi_meta: dict, card_meta: dict,
                       resfinder_meta: dict | None = None) -> list[dict]:
    """
    Parse a FASTA file and return a list of dicts ready for DB insertion.
    Also returns the raw sequence keyed by jrc_id.
    """
    records = []
    for header, seq in parse_fasta(fasta_path):
        if not seq:
            continue
        seq_id   = jrc_id(seq)
        seq_md5  = md5hex(seq)
        rec = {
            "jrc_id":         seq_id,
            "sequence_md5":   seq_md5,
            "sequence":       seq,
            "source":         source,
            "original_header": header,
            "source_acc":     None,
            "gene_name":      None,
            "product_name":   None,
            "drug_class":     None,
            "resistance_mechanism": None,
            "gene_family":    None,
            "aro_accession":  None,
            "refseq_protein": None,
            "refseq_nucleotide": None,
            "genbank_protein":  None,
            "genbank_nucleotide": None,
            "amr_class":      None,
            "amr_subclass":   None,
            "scope":          None,
            # NCBI additional
            "allele":         None,
            "ncbi_type":      None,
            "ncbi_subtype":   None,
            # ResFinder additional
            "pmid":           None,
            "notes":          None,
            "required_gene":  None,
            # CARD additional
            "card_short_name": None,
            "card_cvterm_id":  None,
            "card_model_id":   None,
            "card_model_sequence_id": None,
        }

        if source == "CARD":
            parsed = parse_card_header(header)
            rec.update(parsed)
            acc = rec.get("source_acc", "")
            if acc in card_meta:
                m = card_meta[acc]
                rec["product_name"]          = m.get("Model Name") or m.get("ARO Name")
                rec["drug_class"]            = m.get("Drug Class")
                rec["resistance_mechanism"]  = m.get("Resistance Mechanism")
                rec["gene_family"]           = m.get("AMR Gene Family")
                rec["refseq_protein"]        = m.get("Protein Accession")
                rec["genbank_nucleotide"]    = m.get("DNA Accession")
                if not rec.get("aro_accession"):
                    rec["aro_accession"]     = m.get("ARO Accession")
                if not rec.get("gene_name"):
                    rec["gene_name"]         = m.get("ARO Name")
                rec["card_short_name"]       = m.get("CARD Short Name")
                rec["card_cvterm_id"]        = m.get("CVTERM ID")
                rec["card_model_id"]         = m.get("Model ID")
                rec["card_model_sequence_id"] = m.get("Model Sequence ID")

        elif source == "NCBI":
            parsed = parse_ncbi_header(header)
            rec.update(parsed)
            acc = rec.get("source_acc", "")
            m = ncbi_meta.get(acc) or ncbi_meta.get(acc.split(".")[0])
            if m:
                rec["product_name"]          = m.get("productName")
                rec["gene_family"]           = m.get("geneFamily")
                rec["allele"]                = m.get("allele")
                rec["amr_class"]             = m.get("class")
                rec["amr_subclass"]          = m.get("subclass")
                rec["ncbi_type"]             = m.get("type")
                rec["ncbi_subtype"]          = m.get("subtype")
                rec["scope"]                 = m.get("scope")
                rec["refseq_nucleotide"]     = (m.get("refseqNucleotide") or {}).get("accessionVersion")
                rec["refseq_protein"]        = (m.get("refseqProtein") or {}).get("accessionVersion")
                rec["genbank_nucleotide"]    = (m.get("genbankNucleotide") or {}).get("accessionVersion")
                rec["genbank_protein"]       = (m.get("genbankProtein") or {}).get("accessionVersion")
                # gene_name: allele first, then geneFamily
                if not rec.get("gene_name"):
                    rec["gene_name"] = rec.get("allele") or rec.get("gene_family")

        elif source == "RESFINDER":
            parsed = parse_resfinder_header(header)
            rec.update(parsed)
            if resfinder_meta:
                m = resfinder_meta.get(header)
                if m:
                    rec["drug_class"]           = m["drug_class"]
                    rec["resistance_mechanism"] = m["resistance_mechanism"]
                    rec["product_name"]         = m["phenotype"]
                    rec["pmid"]                 = m["pmid"]
                    rec["notes"]                = m["notes"]
                    rec["required_gene"]        = m["required_gene"]

        records.append(rec)
    return records


# ── Clustering ────────────────────────────────────────────────────────────────

def cluster_records(records: list[dict]) -> dict[str, dict]:
    """
    Group by jrc_id (== 100% sequence identity).
    For each cluster choose the representative with the LONGEST original_header.
    Returns dict jrc_id -> representative record.
    """
    groups: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        groups[r["jrc_id"]].append(r)

    representatives = {}
    for jid, members in groups.items():
        rep = max(members, key=lambda r: len(r["original_header"]))
        representatives[jid] = rep
    return representatives


# ── Database helpers ──────────────────────────────────────────────────────────

def upsert_sequences(cur, sequences: list[dict]):
    sql = """
        INSERT INTO amr.sequence (jrc_id, sequence_md5, sequence)
        VALUES %s
        ON CONFLICT (jrc_id) DO NOTHING
    """
    data = [(s["jrc_id"], s["sequence_md5"], s["sequence"]) for s in sequences]
    execute_values(cur, sql, data)


def upsert_clusters(cur, representatives: dict[str, dict]) -> dict[str, str]:
    """
    Insert one cluster per unique sequence (jrc_id).
    Returns mapping jrc_id -> cluster_id (e.g. 'JRCGRP_000001').
    """
    jrc_to_cluster = {}
    for jid in representatives:
        cur.execute("""
            INSERT INTO amr.cluster (representative_jrc)
            VALUES (%s)
            ON CONFLICT (representative_jrc) DO NOTHING
            RETURNING cluster_id
        """, (jid,))
        row = cur.fetchone()
        if row:
            jrc_to_cluster[jid] = row[0]
        else:
            cur.execute("SELECT cluster_id FROM amr.cluster WHERE representative_jrc = %s", (jid,))
            jrc_to_cluster[jid] = cur.fetchone()[0]
    return jrc_to_cluster


def insert_genes(cur, records: list[dict], jrc_to_cluster: dict[str, str]):
    cols = [
        "jrc_id", "cluster_id", "source", "original_header",
        "source_acc", "gene_name", "product_name", "drug_class",
        "resistance_mechanism", "gene_family", "aro_accession",
        "refseq_protein", "refseq_nucleotide", "genbank_protein",
        "genbank_nucleotide", "amr_class", "amr_subclass", "scope",
        "allele", "ncbi_type", "ncbi_subtype",
        "pmid", "notes", "required_gene",
        "card_short_name", "card_cvterm_id", "card_model_id", "card_model_sequence_id",
    ]
    sql = f"""
        INSERT INTO amr.gene ({", ".join(cols)})
        VALUES %s
        ON CONFLICT DO NOTHING
    """
    data = [
        tuple(
            jrc_to_cluster.get(r["jrc_id"]) if col == "cluster_id" else r.get(col)
            for col in cols
        )
        for r in records
        if r["jrc_id"] in jrc_to_cluster
    ]
    execute_values(cur, sql, data)


def upsert_sequence_metadata(cur, all_records: list[dict]):
    """
    For each unique (jrc_id, source) pair pick the record with the most
    metadata (longest header as tie-breaker) and upsert into
    amr.sequence_metadata.
    """
    # Best record per (jrc_id, source)
    best: dict[tuple, dict] = {}
    for r in all_records:
        key = (r["jrc_id"], r["source"])
        prev = best.get(key)
        if prev is None or len(r["original_header"]) > len(prev["original_header"]):
            best[key] = r

    cols = [
        "jrc_id", "source", "gene_name", "product_name", "drug_class",
        "resistance_mechanism", "gene_family", "aro_accession",
        "refseq_protein", "refseq_nucleotide", "genbank_protein",
        "genbank_nucleotide", "amr_class", "amr_subclass", "scope",
        "allele", "ncbi_type", "ncbi_subtype",
        "pmid", "notes", "required_gene",
        "card_short_name", "card_cvterm_id", "card_model_id", "card_model_sequence_id",
    ]
    sql = f"""
        INSERT INTO amr.sequence_metadata ({", ".join(cols)})
        VALUES %s
        ON CONFLICT (jrc_id, source) DO UPDATE SET
            {", ".join(f"{c} = EXCLUDED.{c}" for c in cols if c not in ("jrc_id", "source"))}
    """
    data = [tuple(r.get(c) for c in cols) for r in best.values()]
    execute_values(cur, sql, data)


# ── Harmonisation loaders ─────────────────────────────────────────────────────

def load_harmonise(cur):
    """Load all harmonise/ TSV files into the amr.drug_* tables."""
    import csv

    def _none(v):
        """Empty string / None → None for nullable columns."""
        v = (v or "").strip()
        return v or None

    # ── drug_class ──
    if DRUG_CLASS_TSV.exists():
        with open(DRUG_CLASS_TSV) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        data = [
            (
                r["canonical_name"].strip(),
                _none(( r.get("aro_accession") or "" ).strip()),
                _none(( r.get("category") or "" ).strip()),
                _none(( r.get("resfinder_alias") or "" ).strip()),
                _none(( r.get("ncbi_alias") or "" ).strip()),
                _none(( r.get("card_abbrev") or "" ).strip()),
                _none(( r.get("card_class_abbrev") or "" ).strip()),
                _none(( r.get("notes") or "" ).strip()),
            )
            for r in rows
            if ( r.get("canonical_name") or "" ).strip()
        ]
        execute_values(cur, """
            INSERT INTO amr.drug_class
                (canonical_name, aro_accession, category, resfinder_alias,
                 ncbi_alias, card_abbrev, card_class_abbrev, notes)
            VALUES %s
            ON CONFLICT (canonical_name) DO UPDATE SET
                aro_accession     = EXCLUDED.aro_accession,
                category          = EXCLUDED.category,
                resfinder_alias   = EXCLUDED.resfinder_alias,
                ncbi_alias        = EXCLUDED.ncbi_alias,
                card_abbrev       = EXCLUDED.card_abbrev,
                card_class_abbrev = EXCLUDED.card_class_abbrev,
                notes             = EXCLUDED.notes
        """, data)
        print(f"  drug_class: {len(data)} rows", file=sys.stderr)

    # ── drug ──
    if DRUG_TSV.exists():
        with open(DRUG_TSV) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        data = [
            (
                r["canonical_name"].strip(),
                r.get("is_combination", "False").strip() == "True",
                _none(( r.get("components") or "" ).strip()),
                r.get("context", "clinical").strip() or "clinical",
                _none(( r.get("card_abbrev") or "" ).strip()),
                _none(( r.get("atc_code") or "" ).strip()),
                _none(( r.get("inchikey") or "" ).strip()),
                _none(( r.get("pubchem_cid") or "" ).strip()),
                _none(( r.get("chebi_id") or "" ).strip()),
                _none(( r.get("sources") or "" ).strip()),
            )
            for r in rows
            if ( r.get("canonical_name") or "" ).strip()
        ]
        execute_values(cur, """
            INSERT INTO amr.drug
                (canonical_name, is_combination, components, context,
                 card_abbrev, atc_code, inchikey, pubchem_cid, chebi_id, sources)
            VALUES %s
            ON CONFLICT (canonical_name) DO UPDATE SET
                is_combination = EXCLUDED.is_combination,
                components     = EXCLUDED.components,
                context        = EXCLUDED.context,
                card_abbrev    = EXCLUDED.card_abbrev,
                atc_code       = EXCLUDED.atc_code,
                inchikey       = EXCLUDED.inchikey,
                pubchem_cid    = EXCLUDED.pubchem_cid,
                chebi_id       = EXCLUDED.chebi_id,
                sources        = EXCLUDED.sources
        """, data)
        print(f"  drug:       {len(data)} rows", file=sys.stderr)

    # ── drug_alias ──
    if DRUG_ALIAS_TSV.exists():
        with open(DRUG_ALIAS_TSV) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        data = [
            (
                r["canonical_name"].strip(),
                r["alias"].strip(),
                r.get("alias_type", "source_name").strip(),
                ( r.get("source") or "" ).strip(),
            )
            for r in rows
            if ( r.get("canonical_name") or "" ).strip() and ( r.get("alias") or "" ).strip()
        ]
        execute_values(cur, """
            INSERT INTO amr.drug_alias (canonical_name, alias, alias_type, source)
            VALUES %s
            ON CONFLICT (canonical_name, alias, source) DO NOTHING
        """, data)
        print(f"  drug_alias: {len(data)} rows", file=sys.stderr)

    # ── drug_class_member ──
    if DRUG_CLASS_MEMBER_TSV.exists():
        with open(DRUG_CLASS_MEMBER_TSV) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        data = [
            (
                r["canonical_drug"].strip(),
                r["canonical_class"].strip(),
                _none(( r.get("aro_accession") or "" ).strip()),
                _none(( r.get("category") or "" ).strip()),
                r.get("source", "curated").strip() or "curated",
            )
            for r in rows
            if ( r.get("canonical_drug") or "" ).strip() and ( r.get("canonical_class") or "" ).strip()
        ]
        execute_values(cur, """
            INSERT INTO amr.drug_class_member
                (canonical_drug, canonical_class, aro_accession, category, evidence_source)
            VALUES %s
            ON CONFLICT (canonical_drug, canonical_class) DO UPDATE SET
                aro_accession   = EXCLUDED.aro_accession,
                category        = EXCLUDED.category,
                evidence_source = EXCLUDED.evidence_source
        """, data)
        print(f"  drug_class_member: {len(data)} rows", file=sys.stderr)

    # ── gene_drug_link ──
    if GENE_DRUG_LINK_TSV.exists():
        with open(GENE_DRUG_LINK_TSV) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        data = [
            (
                ( r.get("gene_name") or "" ).strip(),
                _none(( r.get("accession") or "" ).strip()),
                _none(( r.get("element_type") or "" ).strip()),
                _none(( r.get("ncbi_class_raw") or "" ).strip()),
                _none(( r.get("ncbi_subclass_raw") or "" ).strip()),
                _none(( r.get("canonical_class_token") or "" ).strip()),
                _none(( r.get("canonical_drug_token") or "" ).strip()),
            )
            for r in rows
        ]
        execute_values(cur, """
            INSERT INTO amr.gene_drug_link
                (gene_name, accession, element_type,
                 ncbi_class_raw, ncbi_subclass_raw,
                 canonical_class_token, canonical_drug_token)
            VALUES %s
        """, data)
        print(f"  gene_drug_link: {len(data)} rows", file=sys.stderr)

    # ── card_gene_class ──
    if CARD_GENE_CLASS_TSV.exists():
        with open(CARD_GENE_CLASS_TSV) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        # Deduplicate on primary key (aro_accession, canonical_class)
        seen_cgc: set = set()
        data = []
        for r in rows:
            aro = r["aro_accession"].strip()
            cls = r["canonical_class"].strip()
            if not aro or not cls:
                continue
            pk = (aro, cls)
            if pk in seen_cgc:
                continue
            seen_cgc.add(pk)
            data.append((
                aro,
                r["gene_name"].strip(),
                _none(( r.get("card_short_name") or "" ).strip()),
                _none(( r.get("gene_family") or "" ).strip()),
                _none(( r.get("resistance_mechanism") or "" ).strip()),
                r["drug_class_raw"].strip(),
                cls,
            ))
        execute_values(cur, """
            INSERT INTO amr.card_gene_class
                (aro_accession, gene_name, card_short_name,
                 gene_family, resistance_mechanism,
                 drug_class_raw, canonical_class)
            VALUES %s
            ON CONFLICT (aro_accession, canonical_class) DO UPDATE SET
                gene_name            = EXCLUDED.gene_name,
                card_short_name      = EXCLUDED.card_short_name,
                gene_family          = EXCLUDED.gene_family,
                resistance_mechanism = EXCLUDED.resistance_mechanism,
                drug_class_raw       = EXCLUDED.drug_class_raw
        """, data)
        print(f"  card_gene_class: {len(data)} rows", file=sys.stderr)


# ── Sequence → drug class / drug population ───────────────────────────────────

def populate_sequence_links(cur):
    """
    Populate amr.sequence_drug_class and amr.sequence_drug from the three
    source-specific evidence paths.

    CARD path  : gene.aro_accession → card_gene_class → drug_class
    NCBI path  : gene.gene_name     → gene_drug_link  → drug_class / drug
    ResFinder  : gene.drug_class text → drug_class.resfinder_alias (case-insensitive)
    """

    # ── Step 1: collect raw (jrc_id, canonical_class, source) triples ──
    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS _sdc_raw (
            jrc_id          VARCHAR(13),
            canonical_class TEXT,
            evidence_source VARCHAR(20)
        ) ON COMMIT DROP
    """)

    # CARD path
    cur.execute("""
        INSERT INTO _sdc_raw (jrc_id, canonical_class, evidence_source)
        SELECT DISTINCT g.jrc_id, cgc.canonical_class, 'CARD'
        FROM amr.gene g
        JOIN amr.card_gene_class cgc ON g.aro_accession = cgc.aro_accession
        WHERE g.source = 'CARD'
          AND g.aro_accession IS NOT NULL
    """)
    cur.execute("SELECT count(*) FROM _sdc_raw WHERE evidence_source = 'CARD'")
    print(f"  CARD  → sequence_drug_class: {cur.fetchone()[0]} raw rows", file=sys.stderr)

    # NCBI path A – via gene_drug_link (gene_name match)
    cur.execute("""
        INSERT INTO _sdc_raw (jrc_id, canonical_class, evidence_source)
        SELECT DISTINCT g.jrc_id, dc.canonical_name, 'NCBI'
        FROM amr.gene g
        JOIN amr.gene_drug_link gdl ON g.gene_name = gdl.gene_name
        JOIN amr.drug_class dc      ON gdl.canonical_class_token = dc.canonical_name
        WHERE g.source = 'NCBI'
          AND gdl.canonical_class_token IS NOT NULL
          AND gdl.canonical_class_token <> ''
    """)

    # NCBI path B – via amr_class field directly (slash-delimited tokens → ncbi_alias)
    cur.execute("""
        INSERT INTO _sdc_raw (jrc_id, canonical_class, evidence_source)
        SELECT DISTINCT g.jrc_id, dc.canonical_name, 'NCBI'
        FROM amr.gene g
        CROSS JOIN LATERAL unnest(string_to_array(g.amr_class, '/')) AS t(token)
        JOIN amr.drug_class dc ON upper(trim(t.token)) = upper(dc.ncbi_alias)
        WHERE g.source = 'NCBI'
          AND g.amr_class IS NOT NULL
    """)
    cur.execute("SELECT count(*) FROM _sdc_raw WHERE evidence_source = 'NCBI'")
    print(f"  NCBI  → sequence_drug_class: {cur.fetchone()[0]} raw rows", file=sys.stderr)

    # ResFinder path – match on resfinder_alias (case-insensitive) or canonical_name
    cur.execute("""
        INSERT INTO _sdc_raw (jrc_id, canonical_class, evidence_source)
        SELECT DISTINCT g.jrc_id, dc.canonical_name, 'RESFINDER'
        FROM amr.gene g
        JOIN amr.drug_class dc
          ON lower(g.drug_class) = lower(dc.resfinder_alias)
          OR lower(g.drug_class) = lower(dc.canonical_name)
        WHERE g.source = 'RESFINDER'
          AND g.drug_class IS NOT NULL
    """)
    cur.execute("SELECT count(*) FROM _sdc_raw WHERE evidence_source = 'RESFINDER'")
    print(f"  RSFND → sequence_drug_class: {cur.fetchone()[0]} raw rows", file=sys.stderr)

    # ── Step 2: aggregate evidence sources and upsert ──
    cur.execute("""
        INSERT INTO amr.sequence_drug_class (jrc_id, canonical_class, evidence_sources)
        SELECT
            jrc_id,
            canonical_class,
            string_agg(DISTINCT evidence_source, '|' ORDER BY evidence_source) AS evidence_sources
        FROM _sdc_raw
        GROUP BY jrc_id, canonical_class
        ON CONFLICT (jrc_id, canonical_class) DO UPDATE
            SET evidence_sources = EXCLUDED.evidence_sources
    """)
    cur.execute("SELECT count(*) FROM amr.sequence_drug_class")
    print(f"  sequence_drug_class total: {cur.fetchone()[0]} rows", file=sys.stderr)

    # ── Step 3: sequence → drug (NCBI path only) ──
    cur.execute("""
        INSERT INTO amr.sequence_drug (jrc_id, canonical_drug, evidence_sources)
        SELECT DISTINCT g.jrc_id, d.canonical_name, 'NCBI'
        FROM amr.gene g
        JOIN amr.gene_drug_link gdl ON g.gene_name = gdl.gene_name
        JOIN amr.drug d             ON gdl.canonical_drug_token = d.canonical_name
        WHERE g.source = 'NCBI'
          AND gdl.canonical_drug_token IS NOT NULL
          AND gdl.canonical_drug_token <> ''
        ON CONFLICT (jrc_id, canonical_drug) DO UPDATE
            SET evidence_sources = EXCLUDED.evidence_sources
    """)
    cur.execute("SELECT count(*) FROM amr.sequence_drug")
    print(f"  sequence_drug total:        {cur.fetchone()[0]} rows", file=sys.stderr)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Load AMR genes into ROSETTADB")
    ap.add_argument("--dsn",       default=DEFAULT_DSN)
    ap.add_argument("--resfinder", type=Path, default=DEFAULT_RESFINDER)
    ap.add_argument("--card",      type=Path, default=DEFAULT_CARD)
    ap.add_argument("--ncbi",      type=Path, default=DEFAULT_NCBI)
    ap.add_argument("--schema",    type=Path, default=Path(__file__).parent / "schema.sql")
    ap.add_argument("--skip-harmonise", action="store_true",
                    help="Skip loading harmonise/ drug vocabulary tables")
    ap.add_argument("--skip-links", action="store_true",
                    help="Skip populating sequence_drug_class / sequence_drug link tables")
    ap.add_argument("--harmonise-only", action="store_true",
                    help="Skip FASTA parsing; only reload harmonise/ tables and re-create "
                         "sequence links (requires sequences already in DB)")
    args = ap.parse_args()

    # ── Connect to DB (needed for all paths) ──
    print(f"\nConnecting to: {args.dsn}", file=sys.stderr)
    conn = psycopg2.connect(args.dsn)
    conn.autocommit = False
    cur = conn.cursor()

    if args.harmonise_only:
        print("Applying schema (new tables only) …", file=sys.stderr)
        with open(args.schema) as fh:
            cur.execute(fh.read())
        print("\nReloading harmonised drug vocabulary …", file=sys.stderr)
        load_harmonise(cur)
        print("\nRepopulating sequence → drug class / drug links …", file=sys.stderr)
        # Truncate and re-fill so stale rows from old vocabulary are removed
        cur.execute("TRUNCATE amr.sequence_drug_class, amr.sequence_drug")
        populate_sequence_links(cur)
        conn.commit()
        cur.close()
        conn.close()
        print("\nDone (harmonise-only)!", file=sys.stderr)
        return

    # ── Load auxiliary metadata ──
    print("Loading NCBI metadata …", file=sys.stderr)
    ncbi_meta = load_ncbi_report(NCBI_REPORT)
    print(f"  {len(ncbi_meta)} accession entries", file=sys.stderr)

    print("Loading CARD ARO index …", file=sys.stderr)
    card_meta = load_card_aro(CARD_ARO)
    print(f"  {len(card_meta)} accession entries", file=sys.stderr)

    print("Loading ResFinder phenotypes …", file=sys.stderr)
    resfinder_meta = load_resfinder_phenotypes(RESFINDER_PHENOTYPES)
    print(f"  {len(resfinder_meta)} gene entries", file=sys.stderr)

    # ── Parse FASTA sources ──
    all_records: list[dict] = []

    sources_found = []
    for source, path in [
        ("RESFINDER", args.resfinder),
        ("CARD",      args.card),
        ("NCBI",      args.ncbi),
    ]:
        if not path.exists():
            print(f"  WARNING: {source} file not found: {path}", file=sys.stderr)
            continue
        print(f"Parsing {source}: {path} …", file=sys.stderr)
        recs = build_gene_records(source, path, ncbi_meta, card_meta,
                                  resfinder_meta if source == "RESFINDER" else None)
        print(f"  {len(recs)} sequences", file=sys.stderr)
        all_records.extend(recs)
        sources_found.append(source)

    if not all_records:
        print("No sequences found. Exiting.", file=sys.stderr)
        sys.exit(1)

    print(f"\nTotal raw records: {len(all_records)}", file=sys.stderr)

    # ── Cluster at 100% identity ──
    representatives = cluster_records(all_records)
    print(f"Unique clusters (100% id): {len(representatives)}", file=sys.stderr)

    # ── Connect to DB and load schema (normal full-ingest path) ──
    print("Applying schema …", file=sys.stderr)
    with open(args.schema) as fh:
        cur.execute(fh.read())

    # ── Insert sequences ──
    seq_list = [
        {"jrc_id": jid, "sequence_md5": md5hex(r["sequence"]), "sequence": r["sequence"]}
        for jid, r in representatives.items()
    ]
    print(f"Upserting {len(seq_list)} sequences …", file=sys.stderr)
    upsert_sequences(cur, seq_list)

    # ── Insert clusters ──
    print("Upserting clusters …", file=sys.stderr)
    jrc_to_cluster = upsert_clusters(cur, representatives)

    # ── Insert all gene records ──
    print(f"Inserting {len(all_records)} gene records …", file=sys.stderr)
    insert_genes(cur, all_records, jrc_to_cluster)

    # ── Upsert canonical metadata per (jrc_id, source) ──
    print("Upserting sequence metadata …", file=sys.stderr)
    upsert_sequence_metadata(cur, all_records)

    # ── Load harmonised drug vocabulary ──
    if not args.skip_harmonise:
        print("\nLoading harmonised drug vocabulary …", file=sys.stderr)
        load_harmonise(cur)

    # ── Populate sequence → drug class / drug link tables ──
    if not args.skip_links:
        print("\nPopulating sequence → drug class / drug links …", file=sys.stderr)
        populate_sequence_links(cur)

    conn.commit()
    cur.close()
    conn.close()

    print("\nDone!", file=sys.stderr)
    print(f"  Sources loaded  : {', '.join(sources_found)}", file=sys.stderr)
    print(f"  Total sequences : {len(seq_list)}", file=sys.stderr)
    print(f"  Total clusters  : {len(jrc_to_cluster)}", file=sys.stderr)
    print(f"  Total gene rows : {len(all_records)}", file=sys.stderr)
    if not args.skip_links:
        print("  sequence_drug_class and sequence_drug tables populated", file=sys.stderr)


if __name__ == "__main__":
    main()
