#!/usr/bin/env python3
import os
import re
import sqlite3
import subprocess
import logging
import argparse
from typing import Tuple

def ensure_genome_tables(conn: sqlite3.Connection):
    cur = conn.cursor()
    cur.executescript("""
    PRAGMA foreign_keys=ON;

    CREATE TABLE IF NOT EXISTS assemblies (
        bacdive_id TEXT,
        accession TEXT,
        source TEXT,
        paired_asm TEXT,
        organism TEXT,
        taxid INTEGER,
        assembly_level TEXT,
        refseq_category TEXT,
        bioproject TEXT,
        biosample TEXT,
        isolate TEXT,
        release_date TEXT,
        ftp_path TEXT,
        chosen INTEGER DEFAULT 0,
        last_checked TEXT,
        notes TEXT,
        PRIMARY KEY (bacdive_id, accession)
    );

    CREATE TABLE IF NOT EXISTS strain_genome_links (
        bacdive_id TEXT PRIMARY KEY,
        accession TEXT
    );

    CREATE TABLE IF NOT EXISTS assembly_files (
        bacdive_id TEXT,
        accession TEXT,
        filetype TEXT,   -- fasta, gff, protein-fasta, etc.
        path TEXT,
        PRIMARY KEY (accession, filetype, path)
    );
    """)
    conn.commit()

def export_accessions_for_download(db_path: str, out_txt: str) -> int:
    """
    Export chosen assembly accessions (GCF_/GCA_) to a newline-delimited file.
    Priority: strain_genome_links (if table exists) → fall back to accessions table.
    """
    conn = sqlite3.connect(db_path)
    ensure_genome_tables(conn)
    cur = conn.cursor()

    # Is strain_genome_links populated?
    try:
        rows = cur.execute("""
          SELECT l.accession
          FROM strain_genome_links l
          WHERE l.accession LIKE 'GCF_%' OR l.accession LIKE 'GCA_%'
        """).fetchall()
    except sqlite3.OperationalError:
        rows = []

    accs = [r[0] for r in rows]

    if not accs:
        # Fallback: any accession directly from BacDive scrape (if it had any)
        logging.info("No chosen links yet; falling back to 'accessions' table.")
        try:
            rows = cur.execute("""
              SELECT DISTINCT acc_no
              FROM accessions
              WHERE acc_no LIKE 'GCF_%' OR acc_no LIKE 'GCA_%'
            """).fetchall()
        except sqlite3.OperationalError:
            rows = []
        accs = [r[0] for r in rows]

    conn.close()

    if not accs:
        logging.warning("No GCF_/GCA_ accessions found to export. Run the genome linker first.")
        return 0

    os.makedirs(os.path.dirname(out_txt) or ".", exist_ok=True)
    with open(out_txt, "w") as f:
        for a in accs:
            f.write(a + "\n")
    logging.info(f"Wrote {len(accs)} accessions → {out_txt}")
    return len(accs)

def download_assemblies_with_ngd(
    accessions_txt: str,
    out_dir: str = "genomes",
    section: str = "refseq",     # 'refseq' or 'genbank'
    formats: str = "fasta,gff",
    parallel: int = 4,
    dry_run: bool = False
) -> Tuple[int, str]:
    cmd = [
        "ncbi-genome-download",
        "--section", section,
        "--formats", formats,
        "--assembly-accessions", accessions_txt,
        "--output-folder", out_dir,
        "--parallel", str(parallel),
        "bacteria"
    ]
    if dry_run:
        cmd.insert(1, "--dry-run")

    logging.info("Running: " + " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        logging.error(res.stderr.strip() or res.stdout.strip())
    else:
        logging.info((res.stdout or "").strip()[:2000])
    return res.returncode, " ".join(cmd)

def ingest_download_manifest(db_path: str, out_dir: str) -> int:
    """
    Scan NGD output and record file paths into assembly_files.
    """
    conn = sqlite3.connect(db_path)
    ensure_genome_tables(conn)
    cur = conn.cursor()

    n = 0
    for root, _, files in os.walk(out_dir):
        for name in files:
            m = re.search(r'(GC[AF]_\d+\.\d+)', name)
            if not m:
                continue
            acc = m.group(1)

            if name.endswith((".fna.gz", "_genomic.fna.gz")):
                ft = "fasta"
            elif name.endswith(".gff.gz"):
                ft = "gff"
            elif name.endswith("_protein.faa.gz"):
                ft = "protein-fasta"
            else:
                continue

            # Try to map to a bacdive_id if we know it
            row = cur.execute("SELECT bacdive_id FROM assemblies WHERE accession=?", (acc,)).fetchone()
            bd = row[0] if row else None

            cur.execute("""
              INSERT OR IGNORE INTO assembly_files(bacdive_id, accession, filetype, path)
              VALUES (?,?,?,?)
            """, (bd, acc, ft, os.path.abspath(os.path.join(root, name))))
            n += cur.rowcount != 0

    conn.commit()
    conn.close()
    logging.info(f"Recorded {n} downloaded files into assembly_files")
    return n

def download_missing_genomes(
    db_path: str,
    out_dir: str = "genomes",
    section: str = "refseq",
    formats: str = "fasta,gff",
    parallel: int = 4,
    dry_run: bool = False
):
    os.makedirs(out_dir, exist_ok=True)
    acc_file = os.path.join(out_dir, "accessions.txt")
    n = export_accessions_for_download(db_path, acc_file)
    if n == 0:
        logging.warning("Nothing to download (no accessions).")
        return
    rc, _ = download_assemblies_with_ngd(
        acc_file, out_dir, section=section, formats=formats, parallel=parallel, dry_run=dry_run
    )
    if rc == 0 and not dry_run:
        ingest_download_manifest(db_path, out_dir)

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="bacdive_traits.sqlite")
    ap.add_argument("--out", default="genomes")
    ap.add_argument("--section", default="refseq", choices=["refseq","genbank"])
    ap.add_argument("--formats", default="fasta,gff")
    ap.add_argument("--parallel", default=4, type=int)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--log-level", default="INFO")
    return ap.parse_args()

def main():
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(levelname)s: %(message)s")
    download_missing_genomes(
        args.db, out_dir=args.out, section=args.section, formats=args.formats,
        parallel=args.parallel, dry_run=args.dry_run
    )

if __name__ == "__main__":
    main()

