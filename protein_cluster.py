#!/usr/bin/env python3
"""
protein_cluster.py
==================
Post-processing step: translate every nucleotide sequence in amr.sequence,
group sequences that encode the same protein, and populate amr.protein +
set amr.sequence.protein_id.

Run AFTER ingest.py:
    python protein_cluster.py [--dsn <connstr>] [--dry-run]

Translation rules
-----------------
1. Strip trailing stop codon (TAA / TAG / TGA) if present — some FASTA
   entries include it, others do not.  We normalise before hashing so both
   forms map to the same protein.

2. Force the first codon to Met (M) regardless of its nucleotide spelling.
   Bacteria routinely use alternative start codons (GTG, TTG, CTG, ATT,
   ATC, ATA); the ribosome always loads fMet.

3. Require length % 3 == 0 after stop-stripping.  Sequences that fail
   this (partial / truncated CDSs) receive protein_id = NULL.

4. Ambiguous nucleotides (N, R, Y, …) or internal stop codons → NULL.
   These indicate assembly artefacts or non-CDS entries.

Output
------
- amr.protein  : one row per unique protein (protein_md5 is the dedup key)
- amr.sequence : protein_id column populated for translatable sequences

Idempotent: re-running updates only rows whose protein_id is currently NULL
or whose nucleotide sequence has changed (detected via sequence_md5 change,
which cannot happen without a new jrc_id — so in practice safe to re-run).
"""

import argparse
import hashlib
import os
import sys
from collections import defaultdict

import psycopg2
import psycopg2.extras

# ---------------------------------------------------------------------------
# Genetic code
# ---------------------------------------------------------------------------

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Codons that act as start codons in bacteria (ribosome loads fMet)
ALT_START_CODONS: frozenset[str] = frozenset({
    "GTG", "TTG", "CTG", "ATT", "ATC", "ATA",
})

STOP_CODONS: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})


# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------

def translate(seq: str) -> str | None:
    """
    Translate a nucleotide CDS to its single-letter amino acid sequence.

    Returns None (sequence gets protein_id = NULL) when:
    - Length is not divisible by 3 after stop-stripping
    - Any codon is absent from the standard table (ambiguous nucleotides)
    - An internal stop codon is encountered

    The returned protein string does NOT include the trailing '*'.
    """
    seq = seq.upper()

    # 1. Strip trailing stop codon if present
    if len(seq) >= 3 and seq[-3:] in STOP_CODONS:
        seq = seq[:-3]

    # 2. Must be divisible by 3
    if len(seq) == 0 or len(seq) % 3 != 0:
        return None

    aas: list[str] = []
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]

        # 3. Force alternative start to Met
        if i == 0 and codon in ALT_START_CODONS:
            aas.append("M")
            continue

        aa = CODON_TABLE.get(codon)
        if aa is None:
            # Ambiguous / degenerate nucleotide in this codon
            return None
        if aa == "*":
            # Internal stop codon — frameshift or non-CDS
            return None
        aas.append(aa)

    return "".join(aas)


def protein_md5(protein: str) -> str:
    return hashlib.md5(protein.encode()).hexdigest()


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

def connect(dsn: str) -> psycopg2.extensions.connection:
    return psycopg2.connect(dsn)


def fetch_sequences(cur) -> list[tuple[str, str]]:
    """Return (jrc_id, sequence) for all rows where protein_id IS NULL."""
    cur.execute("""
        SELECT jrc_id, sequence
        FROM amr.sequence
        WHERE protein_id IS NULL
        ORDER BY jrc_id
    """)
    return cur.fetchall()


def fetch_existing_proteins(cur) -> dict[str, str]:
    """Return protein_md5 → protein_id for all rows already in amr.protein."""
    cur.execute("SELECT protein_md5, protein_id FROM amr.protein")
    return {row[0]: row[1] for row in cur.fetchall()}


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def build_protein_groups(
    rows: list[tuple[str, str]],
) -> tuple[dict[str, str], dict[str, list[str]]]:
    """
    Translate all sequences and group jrc_ids by protein_md5.

    Returns:
      translations  : jrc_id → protein string  (only translatable entries)
      groups        : protein_md5 → [jrc_id, ...]
    """
    translations: dict[str, str] = {}
    groups: dict[str, list[str]] = defaultdict(list)
    skipped = 0

    for jrc_id, seq in rows:
        protein = translate(seq)
        if protein is None:
            skipped += 1
            continue
        md5 = protein_md5(protein)
        translations[jrc_id] = protein
        groups[md5].append(jrc_id)

    print(f"  Translatable       : {len(translations)}")
    print(f"  Skipped (non-CDS)  : {skipped}")
    print(f"  Unique proteins    : {len(groups)}")
    synonymous = sum(1 for v in groups.values() if len(v) > 1)
    print(f"  With synonymous nt : {synonymous}")

    return translations, groups


def choose_representative(jrc_ids: list[str], cur) -> str:
    """
    From a list of jrc_ids encoding the same protein, choose the one whose
    nucleotide sequence is longest (longest complete CDS wins; ties broken
    by jrc_id lexicographic order for determinism).
    """
    if len(jrc_ids) == 1:
        return jrc_ids[0]
    cur.execute(
        "SELECT jrc_id, sequence_length FROM amr.sequence WHERE jrc_id = ANY(%s)",
        (jrc_ids,),
    )
    rows = cur.fetchall()
    # longest sequence, then lexicographically smallest jrc_id for tie-breaking
    return max(rows, key=lambda r: (r[1], r[0]))[0]


def upsert_proteins(
    groups: dict[str, list[str]],
    translations: dict[str, str],
    existing: dict[str, str],
    cur,
    dry_run: bool,
) -> dict[str, str]:
    """
    Insert new proteins into amr.protein and return protein_md5 → protein_id
    for all groups (merged with existing).
    """
    md5_to_pid: dict[str, str] = dict(existing)

    new_proteins = [
        (md5, jrc_ids)
        for md5, jrc_ids in groups.items()
        if md5 not in existing
    ]

    if not new_proteins:
        print(f"  No new proteins to insert.")
        return md5_to_pid

    print(f"  Inserting {len(new_proteins)} new protein rows …")

    for md5, jrc_ids in new_proteins:
        protein_seq = translations[jrc_ids[0]]  # all identical; any member works
        rep_jrc = choose_representative(jrc_ids, cur)

        if dry_run:
            # Allocate a placeholder id so downstream code can still build
            # the protein_id mapping without touching the DB.
            md5_to_pid[md5] = f"JRCPRO_DRY{md5[:4]}"
            continue

        cur.execute(
            """
            INSERT INTO amr.protein (protein_md5, protein_sequence, representative_jrc)
            VALUES (%s, %s, %s)
            ON CONFLICT (protein_md5) DO NOTHING
            RETURNING protein_id
            """,
            (md5, protein_seq, rep_jrc),
        )
        row = cur.fetchone()
        if row:
            md5_to_pid[md5] = row[0]
        else:
            # Race condition or duplicate — fetch the existing id
            cur.execute(
                "SELECT protein_id FROM amr.protein WHERE protein_md5 = %s", (md5,)
            )
            md5_to_pid[md5] = cur.fetchone()[0]

    return md5_to_pid


def update_sequence_protein_ids(
    groups: dict[str, list[str]],
    md5_to_pid: dict[str, str],
    cur,
    dry_run: bool,
) -> int:
    """Set protein_id on amr.sequence for all translatable entries."""
    pairs: list[tuple[str, str]] = []  # (protein_id, jrc_id)
    for md5, jrc_ids in groups.items():
        pid = md5_to_pid[md5]
        for jrc_id in jrc_ids:
            pairs.append((pid, jrc_id))

    if not pairs:
        return 0

    if not dry_run:
        psycopg2.extras.execute_batch(
            cur,
            "UPDATE amr.sequence SET protein_id = %s WHERE jrc_id = %s",
            pairs,
            page_size=500,
        )

    return len(pairs)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dsn", default=os.environ.get("ROSETTADB_DSN", ""),
                        help="PostgreSQL connection string "
                             "(default: $ROSETTADB_DSN)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Translate and group sequences but do not write to DB")
    args = parser.parse_args()

    if not args.dsn and not args.dry_run:
        print("ERROR: provide --dsn or set $ROSETTADB_DSN", file=sys.stderr)
        sys.exit(1)

    print("Connecting to database …")
    if args.dry_run:
        # Still need a connection to read sequences; use a read-only approach
        # by rolling back at the end.
        conn = psycopg2.connect(args.dsn) if args.dsn else None
    else:
        conn = psycopg2.connect(args.dsn)

    if conn is None:
        print("ERROR: --dry-run without --dsn: cannot read sequences.", file=sys.stderr)
        sys.exit(1)

    try:
        with conn:
            with conn.cursor() as cur:
                print("Fetching untranslated sequences …")
                rows = fetch_sequences(cur)
                print(f"  Found {len(rows)} sequences with protein_id = NULL")

                if not rows:
                    print("Nothing to do.")
                    return

                print("\nTranslating …")
                translations, groups = build_protein_groups(rows)

                print("\nFetching existing protein records …")
                existing = fetch_existing_proteins(cur)
                print(f"  {len(existing)} proteins already in DB")

                print("\nUpserting proteins …")
                md5_to_pid = upsert_proteins(
                    groups, translations, existing, cur, args.dry_run
                )

                print("\nUpdating amr.sequence.protein_id …")
                n_updated = update_sequence_protein_ids(
                    groups, md5_to_pid, cur, args.dry_run
                )
                print(f"  Updated {n_updated} sequence rows")

                if args.dry_run:
                    conn.rollback()
                    print("\n[DRY RUN] No changes committed.")
                else:
                    print("\nCommitting …")
                    # conn context manager already commits on __exit__

        # Final summary
        if not args.dry_run:
            with conn.cursor() as cur:
                cur.execute("SELECT COUNT(*) FROM amr.protein")
                n_proteins = cur.fetchone()[0]
                cur.execute("SELECT COUNT(*) FROM amr.sequence WHERE protein_id IS NOT NULL")
                n_linked = cur.fetchone()[0]
                cur.execute("SELECT COUNT(*) FROM amr.sequence WHERE protein_id IS NULL")
                n_null = cur.fetchone()[0]

            print(f"\nSummary:")
            print(f"  Total proteins in DB            : {n_proteins}")
            print(f"  Sequences linked to a protein   : {n_linked}")
            print(f"  Sequences without protein (NULL): {n_null}")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
