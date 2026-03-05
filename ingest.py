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
DEFAULT_RESFINDER = BASE / "sources" / "resfinder_db" / "all.fsa"
DEFAULT_CARD      = BASE / "CARD" / "nucleotide_fasta_protein_homolog_model.fasta"
DEFAULT_NCBI      = BASE / "amr_finder_plus" / "ncbi_dataset" / "data" / "nucleotide.fna"

NCBI_REPORT = BASE / "amr_finder_plus" / "ncbi_dataset" / "data" / "data_report.jsonl"
CARD_ARO    = BASE / "CARD" / "aro_index.tsv"


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
    ResFinder headers vary; common: gene_accession_something
    e.g. aac(3)-Ia_1_AF174129
    """
    info: dict = {}
    parts = header.split("_")
    if parts:
        info["gene_name"] = parts[0]
    # look for an accession-like token
    for p in parts[1:]:
        if _ACC_RE.fullmatch(p.split(".")[0] + ".1"):
            info["source_acc"] = p
            break
    if "source_acc" not in info and len(parts) >= 3:
        info["source_acc"] = parts[-1]
    return info


# ── Gene record builder ───────────────────────────────────────────────────────

def build_gene_records(source: str, fasta_path: Path,
                       ncbi_meta: dict, card_meta: dict) -> list[dict]:
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

        elif source == "NCBI":
            parsed = parse_ncbi_header(header)
            rec.update(parsed)
            acc = rec.get("source_acc", "")
            m = ncbi_meta.get(acc) or ncbi_meta.get(acc.split(".")[0])
            if m:
                rec["product_name"]          = m.get("productName")
                rec["gene_family"]           = m.get("geneFamily")
                rec["amr_class"]             = m.get("class")
                rec["amr_subclass"]          = m.get("subclass")
                rec["scope"]                 = m.get("scope")
                rec["refseq_nucleotide"]     = (m.get("refseqNucleotide") or {}).get("accessionVersion")
                rec["refseq_protein"]        = (m.get("refseqProtein") or {}).get("accessionVersion")
                rec["genbank_nucleotide"]    = (m.get("genbankNucleotide") or {}).get("accessionVersion")
                rec["genbank_protein"]       = (m.get("genbankProtein") or {}).get("accessionVersion")
                # gene_name: first token of productName or geneFamily
                if not rec.get("gene_name") and rec.get("gene_family"):
                    rec["gene_name"] = rec["gene_family"]

        elif source == "RESFINDER":
            parsed = parse_resfinder_header(header)
            rec.update(parsed)

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


def upsert_clusters(cur, representatives: dict[str, dict]) -> dict[str, int]:
    """
    Insert one cluster per unique sequence (jrc_id).
    Returns mapping jrc_id -> cluster_id.
    """
    jrc_to_cluster = {}
    for jid in representatives:
        cur.execute("""
            INSERT INTO amr.cluster (representative_jrc)
            VALUES (%s)
            ON CONFLICT DO NOTHING
            RETURNING cluster_id
        """, (jid,))
        row = cur.fetchone()
        if row:
            jrc_to_cluster[jid] = row[0]
        else:
            cur.execute("SELECT cluster_id FROM amr.cluster WHERE representative_jrc = %s", (jid,))
            jrc_to_cluster[jid] = cur.fetchone()[0]
    return jrc_to_cluster


def insert_genes(cur, records: list[dict], jrc_to_cluster: dict[str, int]):
    cols = [
        "jrc_id", "cluster_id", "source", "original_header",
        "source_acc", "gene_name", "product_name", "drug_class",
        "resistance_mechanism", "gene_family", "aro_accession",
        "refseq_protein", "refseq_nucleotide", "genbank_protein",
        "genbank_nucleotide", "amr_class", "amr_subclass", "scope",
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


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Load AMR genes into ROSETTADB")
    ap.add_argument("--dsn",       default=DEFAULT_DSN)
    ap.add_argument("--resfinder", type=Path, default=DEFAULT_RESFINDER)
    ap.add_argument("--card",      type=Path, default=DEFAULT_CARD)
    ap.add_argument("--ncbi",      type=Path, default=DEFAULT_NCBI)
    ap.add_argument("--schema",    type=Path, default=Path(__file__).parent / "schema.sql")
    args = ap.parse_args()

    # ── Load auxiliary metadata ──
    print("Loading NCBI metadata …", file=sys.stderr)
    ncbi_meta = load_ncbi_report(NCBI_REPORT)
    print(f"  {len(ncbi_meta)} accession entries", file=sys.stderr)

    print("Loading CARD ARO index …", file=sys.stderr)
    card_meta = load_card_aro(CARD_ARO)
    print(f"  {len(card_meta)} accession entries", file=sys.stderr)

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
        recs = build_gene_records(source, path, ncbi_meta, card_meta)
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

    # ── Connect to DB and load schema ──
    print(f"\nConnecting to: {args.dsn}", file=sys.stderr)
    conn = psycopg2.connect(args.dsn)
    conn.autocommit = False
    cur = conn.cursor()

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

    conn.commit()
    cur.close()
    conn.close()

    print("\nDone!", file=sys.stderr)
    print(f"  Sources loaded : {', '.join(sources_found)}", file=sys.stderr)
    print(f"  Total sequences: {len(seq_list)}", file=sys.stderr)
    print(f"  Total clusters : {len(jrc_to_cluster)}", file=sys.stderr)
    print(f"  Total gene rows: {len(all_records)}", file=sys.stderr)


if __name__ == "__main__":
    main()
