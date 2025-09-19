#!/usr/bin/env python3
import os, re, sqlite3, subprocess, logging, argparse, time, datetime
import requests
from typing import Dict, List, Optional, Tuple

NCBI_DS = "https://api.ncbi.nlm.nih.gov/datasets/v2"

# -----------------------
# Utilities / logging
# -----------------------

def _http_get_json(url, params=None, email="you@example.org", tries=4, backoff=1.6):
    headers = {"User-Agent": f"syn_genome_linker/0.1 (+{email})"}
    for i in range(tries):
        r = requests.get(url, params=params or {}, headers=headers, timeout=30)
        if r.status_code == 200:
            return r.json()
        if r.status_code in (429,500,502,503,504):
            time.sleep(backoff**i)
            continue
        raise RuntimeError(f"NCBI {r.status_code}: {r.text[:300]}")
    raise RuntimeError("NCBI: too many retries")

def ensure_tables(conn: sqlite3.Connection):
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
        filetype TEXT,
        path TEXT,
        PRIMARY KEY (accession, filetype, path)
    );
    """)
    conn.commit()

# -----------------------
# Accession normalization + resolve
# -----------------------

def normalize_accession(raw: str) -> Optional[str]:
    """Turn sloppy tokens like 'GCA 000027305' into 'GCA_000027305' (no version yet)."""
    if not raw:
        return None
    s = str(raw).strip().upper()
    s = s.replace(" ", "_")
    # common typos: extra underscores
    s = re.sub(r"__+", "_", s)
    # pick only accession-like prefix
    m = re.search(r"(GC[AF]_\d+(\.\d+)?)", s)
    return m.group(1) if m else None

def resolve_canonical(acc_like: str, email="you@example.org") -> Optional[Dict]:
    """
    Resolve to a canonical assembly accession (with version) and metadata using NCBI Datasets.
    Works if acc_like is 'GCA_000027305' (no version) or full 'GCA_000027305.1'.
    Returns first report dict or None.
    """
    # try direct
    url = f"{NCBI_DS}/genome/assembly/{acc_like}/descriptor"
    try:
        js = _http_get_json(url, email=email)
        reps = js.get("reports") or []
        if reps:
            return reps[0]
    except RuntimeError:
        pass

    # try without version / with version variants
    base = acc_like.split(".")[0]
    for guess in (base, f"{base}.1", f"{base}.2", f"{base}.3"):
        try:
            js = _http_get_json(f"{NCBI_DS}/genome/assembly/{guess}/descriptor", email=email)
            reps = js.get("reports") or []
            if reps:
                return reps[0]
        except RuntimeError:
            continue
    return None

def parse_report(bacdive_id: Optional[str], rep: Dict) -> Dict:
    asm = rep.get("assembly", {})
    org = asm.get("organism", {}) or {}
    acc = asm.get("assembly_accession")
    src = "RefSeq" if (acc or "").startswith("GCF_") else ("GenBank" if (acc or "").startswith("GCA_") else "Unknown")
    return {
        "bacdive_id": bacdive_id,
        "accession": acc,
        "source": src,
        "paired_asm": (asm.get("paired_assembly") or {}).get("assembly_accession"),
        "organism": org.get("organism_name"),
        "taxid": org.get("tax_id"),
        "assembly_level": asm.get("assembly_level"),
        "refseq_category": asm.get("refseq_category"),
        "bioproject": (asm.get("bioproject_lineage") or [{}])[0].get("accession") if asm.get("bioproject_lineage") else None,
        "biosample": (asm.get("biosample") or {}).get("accession"),
        "isolate": (asm.get("isolates") or [{}])[0].get("isolate_name") if asm.get("isolates") else None,
        "release_date": asm.get("release_date"),
        "ftp_path": (asm.get("assembly_info") or {}).get("ftp_path"),
        "last_checked": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "notes": None,
    }

# -----------------------
# DB I/O
# -----------------------

def collect_accession_hints(conn: sqlite3.Connection) -> List[Tuple[str,str]]:
    """
    Return list of (bacdive_id, raw_acc) from the accessions table that look like GCA/GCF.
    """
    cur = conn.cursor()
    try:
        rows = cur.execute("""
          SELECT bacdive_id, acc_no FROM accessions
          WHERE acc_no LIKE 'GCF%' OR acc_no LIKE 'GCA%' OR acc_no LIKE 'gcF%' OR acc_no LIKE 'gcA%'
        """).fetchall()
        return [(str(bd) if bd is not None else None, str(a)) for bd,a in rows]
    except sqlite3.OperationalError:
        return []

def upsert_assembly(conn: sqlite3.Connection, meta: Dict):
    cur = conn.cursor()
    cur.execute("""
      INSERT INTO assemblies
        (bacdive_id, accession, source, paired_asm, organism, taxid, assembly_level,
         refseq_category, bioproject, biosample, isolate, release_date, ftp_path,
         chosen, last_checked, notes)
      VALUES
        (:bacdive_id, :accession, :source, :paired_asm, :organism, :taxid, :assembly_level,
         :refseq_category, :bioproject, :biosample, :isolate, :release_date, :ftp_path,
         0, :last_checked, :notes)
      ON CONFLICT(bacdive_id, accession) DO UPDATE SET
         source=excluded.source,
         paired_asm=excluded.paired_asm,
         organism=excluded.organism,
         taxid=excluded.taxid,
         assembly_level=excluded.assembly_level,
         refseq_category=excluded.refseq_category,
         bioproject=excluded.bioproject,
         biosample=excluded.biosample,
         isolate=excluded.isolate,
         release_date=excluded.release_date,
         ftp_path=excluded.ftp_path,
         last_checked=excluded.last_checked
    """, meta)
    conn.commit()

def choose_primary(conn: sqlite3.Connection, bacdive_id: str):
    """
    Simple heuristic: prefer RefSeq (GCF_) > GenBank (GCA_), then higher assembly_level, then newest date.
    """
    cur = conn.cursor()
    rows = cur.execute("""
      SELECT accession, source, assembly_level, release_date
      FROM assemblies WHERE bacdive_id=?""", (bacdive_id,)).fetchall()
    if not rows:
        return
    lvl = {"Complete Genome":4, "Chromosome":3, "Scaffold":2, "Contig":1}
    def score(r):
        acc, src, level, date = r
        s1 = 1 if acc.startswith("GCF_") else 0
        s2 = lvl.get(level, 0)
        try:
            d = datetime.date.fromisoformat((date or "1970-01-01")[:10]).toordinal()
        except Exception:
            d = 0
        return (s1, s2, d)
    best = max(rows, key=score)[0]
    cur.execute("UPDATE assemblies SET chosen=CASE WHEN accession=? THEN 1 ELSE 0 END WHERE bacdive_id=?", (best, bacdive_id))
    cur.execute("""
      INSERT INTO strain_genome_links(bacdive_id, accession)
      VALUES (?,?)
      ON CONFLICT(bacdive_id) DO UPDATE SET accession=excluded.accession
    """, (bacdive_id, best))
    conn.commit()

# -----------------------
# NGD integration
# -----------------------

def write_accession_lists(conn: sqlite3.Connection, out_dir: str) -> Tuple[Optional[str], Optional[str], int, int]:
    """
    Export chosen accessions into two files (refseq.txt for GCF_, genbank.txt for GCA_).
    If no chosen yet, fall back to *all* assemblies we resolved.
    """
    cur = conn.cursor()
    os.makedirs(out_dir, exist_ok=True)

    rows = cur.execute("""
      SELECT l.accession FROM strain_genome_links l WHERE l.accession LIKE 'GCF_%'
    """).fetchall()
    gcf = [r[0] for r in rows]
    if not gcf:
        rows = cur.execute("SELECT DISTINCT accession FROM assemblies WHERE accession LIKE 'GCF_%'").fetchall()
        gcf = [r[0] for r in rows]

    rows = cur.execute("""
      SELECT l.accession FROM strain_genome_links l WHERE l.accession LIKE 'GCA_%'
    """).fetchall()
    gca = [r[0] for r in rows]
    if not gca:
        rows = cur.execute("SELECT DISTINCT accession FROM assemblies WHERE accession LIKE 'GCA_%'").fetchall()
        gca = [r[0] for r in rows]

    path_refseq = os.path.join(out_dir, "refseq.txt") if gcf else None
    path_genbank = os.path.join(out_dir, "genbank.txt") if gca else None

    if gcf:
        with open(path_refseq, "w") as f:
            f.write("\n".join(gcf) + "\n")
        logging.info(f"Wrote {len(gcf)} GCF accessions → {path_refseq}")
    if gca:
        with open(path_genbank, "w") as f:
            f.write("\n".join(gca) + "\n")
        logging.info(f"Wrote {len(gca)} GCA accessions → {path_genbank}")

    return path_refseq, path_genbank, len(gcf), len(gca)

def run_ngd(section: str, acc_file: str, out_dir: str, formats="fasta,gff", parallel=4, dry_run=False) -> int:
    cmd = [
        "ncbi-genome-download",
        "--section", section,
        "--formats", formats,
        "--assembly-accessions", acc_file,
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
    return res.returncode

def ingest_download_manifest(db_path: str, out_dir: str) -> int:
    conn = sqlite3.connect(db_path)
    ensure_tables(conn)
    cur = conn.cursor()
    n = 0
    for root, _, files in os.walk(out_dir):
        for name in files:
            m = re.search(r"(GC[AF]_\d+\.\d+)", name)
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

# -----------------------
# End-to-end pipeline
# -----------------------

def link_and_download(db_path: str, out_dir: str, email="you@example.org",
                      formats="fasta,gff", parallel=4, dry_run=False):
    """
    Self-contained:
      1) Read GCA/GCF tokens from DB (accessions table)
      2) Normalize + resolve with NCBI (write to assemblies, pick primary)
      3) Export GCF/GCA lists separately
      4) Run NGD per section (refseq/genbank)
      5) Record file paths
    """
    os.makedirs(out_dir, exist_ok=True)
    conn = sqlite3.connect(db_path)
    ensure_tables(conn)

    hints = collect_accession_hints(conn)
    if not hints:
        logging.warning("No GCA/GCF-looking tokens in 'accessions' table.")
    else:
        logging.info(f"Found {len(hints)} raw accession tokens in DB")
    seen = set()
    resolved = 0

    for bd_id, raw in hints:
        norm = normalize_accession(raw)
        if not norm or norm in seen:
            continue
        seen.add(norm)
        rep = resolve_canonical(norm, email=email)
        if not rep:
            logging.warning(f"Could not resolve accession: {raw}")
            continue
        meta = parse_report(bd_id, rep)
        upsert_assembly(conn, meta)
        if bd_id:
            choose_primary(conn, bd_id)
        resolved += 1

    logging.info(f"Resolved {resolved} assemblies; choosing primaries done.")

    path_refseq, path_genbank, n_gcf, n_gca = write_accession_lists(conn, out_dir)

    if n_gcf == 0 and n_gca == 0:
        logging.warning("Nothing to download (no canonical accessions).")
        conn.close()
        return

    if path_refseq:
        rc = run_ngd("refseq", path_refseq, out_dir, formats=formats, parallel=parallel, dry_run=dry_run)
        if rc != 0:
            logging.error("RefSeq download step failed.")
    if path_genbank:
        rc = run_ngd("genbank", path_genbank, out_dir, formats=formats, parallel=parallel, dry_run=dry_run)
        if rc != 0:
            logging.error("GenBank download step failed.")

    if not dry_run:
        ingest_download_manifest(db_path, out_dir)
    conn.close()

# -----------------------
# CLI
# -----------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Self-contained linker/downloader: BacDive DB → canonical GCF/GCA → download via NGD."
    )
    ap.add_argument("--db", default="bacdive_traits.sqlite")
    ap.add_argument("--out", default="genomes")
    ap.add_argument("--formats", default="fasta,gff")
    ap.add_argument("--parallel", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--email", default="you@example.org", help="Contact email for NCBI User-Agent")
    ap.add_argument("--log-level", default="INFO")
    return ap.parse_args()

def main():
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(levelname)s: %(message)s")
    link_and_download(args.db, args.out, email=args.email,
                      formats=args.formats, parallel=args.parallel, dry_run=args.dry_run)

if __name__ == "__main__":
    main()

