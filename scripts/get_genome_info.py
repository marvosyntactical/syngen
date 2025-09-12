from typing import Callable, Dict, Iterable, List, Optional
import requests
import time
import datetime
import re
import sqlite3
import argparse


NCBI_DS = "https://api.ncbi.nlm.nih.gov/datasets/v2"

def parse_args():

    ap = argparse.ArgumentParser(description="Get phenotypical data from BacDive API")
    ap.add_argument("--email", type=str, default= "koss@cl.uni-heidelberg.de")

    return ap.parse_args()




def _http_get(url, params=None, email="you@example.org", tries=4, backoff=1.5):
    # NCBI likes a real user agent w/ contact. Be polite.
    headers = {"User-Agent": f"bacdive-linker/0.1 (+{email})"}
    for i in range(tries):
        r = requests.get(url, params=params or {}, headers=headers, timeout=30)
        if r.status_code == 200:
            return r.json()
        # Retry on 429/5xx
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep((backoff ** i) + 0.25)
            continue
        # Hard fail otherwise
        raise RuntimeError(f"NCBI API error {r.status_code}: {r.text[:300]}")
    raise RuntimeError("NCBI API: too many retries")

def _ensure_genome_tables(conn: sqlite3.Connection):
    cur = conn.cursor()
    cur.executescript("""
    PRAGMA foreign_keys=ON;

    CREATE TABLE IF NOT EXISTS assemblies (
        bacdive_id TEXT,
        accession TEXT,                 -- GCF_* or GCA_*
        source TEXT,                    -- 'RefSeq' or 'GenBank' (or 'Unknown')
        paired_asm TEXT,                -- GCA<->GCF counterpart if NCBI reports it
        organism TEXT,
        taxid INTEGER,
        assembly_level TEXT,            -- Complete Genome, Chromosome, Scaffold, Contig
        refseq_category TEXT,           -- reference genome, representative genome, etc.
        bioproject TEXT,
        biosample TEXT,
        isolate TEXT,
        release_date TEXT,
        ftp_path TEXT,
        chosen INTEGER DEFAULT 0,       -- 1 if this is the "primary" pick
        last_checked TEXT,
        notes TEXT,
        PRIMARY KEY (bacdive_id, accession),
        FOREIGN KEY (bacdive_id) REFERENCES strains(bacdive_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS strain_genome_links (
        bacdive_id TEXT PRIMARY KEY,
        accession TEXT,     -- the chosen one for quick joins
        FOREIGN KEY (bacdive_id, accession) REFERENCES assemblies(bacdive_id, accession)
    );

    CREATE INDEX IF NOT EXISTS idx_asm_accession ON assemblies(accession);
    """)
    conn.commit()

def _choose_best(asm_rows: List[Dict]) -> Optional[str]:
    """
    Pick a single 'best' assembly accession from candidates.
    Heuristic: Prefer RefSeq (GCF_) > GenBank (GCA_), then
    refseq_category: reference > representative > other,
    then higher assembly_level, then newest release_date.
    """
    if not asm_rows:
        return None

    level_rank = {
        "Complete Genome": 4,
        "Chromosome": 3,
        "Scaffold": 2,
        "Contig": 1,
        None: 0
    }
    cat_rank = {
        "reference genome": 3,
        "representative genome": 2,
        None: 0
    }

    def score(r):
        s_refseq = 1 if (r.get("accession","").startswith("GCF_")) else 0
        s_cat = cat_rank.get((r.get("refseq_category") or "").lower() or None, 0)
        s_lvl = level_rank.get(r.get("assembly_level"), 0)
        # newer better
        try:
            dt = datetime.date.fromisoformat((r.get("release_date") or "1970-01-01")[:10])
        except Exception:
            dt = datetime.date(1970,1,1)
        return (s_refseq, s_cat, s_lvl, dt.toordinal())

    return max(asm_rows, key=score)["accession"]

def _search_by_accession(acc: str, email: str) -> Optional[Dict]:
    # Direct descriptor fetch
    url = f"{NCBI_DS}/genome/assembly/{acc}/descriptor"
    try:
        js = _http_get(url, email=email)
    except RuntimeError:
        return None
    # The v2 response has a 'reports' list with 'assembly' objects.
    reps = (js.get("reports") or [])
    return reps[0] if reps else None

def _search_by_filter(biosample=None, bioproject=None, taxid=None, isolate=None, email="") -> List[Dict]:
    """
    Use descriptor search with filters. Returns list of 'report' objects.
    """
    url = f"{NCBI_DS}/genome/assembly/descriptors"
    params = {"page_size": 50}
    if biosample:
        params["filters.biosample"] = biosample
    if bioproject:
        params["filters.bioproject"] = bioproject
    if taxid:
        params["filters.taxonomy"] = str(taxid)
    # isolate filter isn’t first-class; try a text query as weak fallback
    if isolate and not (biosample or bioproject):
        params["query"] = isolate
    js = _http_get(url, params=params, email=email)
    return js.get("reports") or []

def _parse_report(bacdive_id: str, rep: Dict) -> Dict:
    asm = rep.get("assembly") or {}
    org = asm.get("organism") or {}
    acc = asm.get("assembly_accession")
    src = "RefSeq" if (acc or "").startswith("GCF_") else ("GenBank" if (acc or "").startswith("GCA_") else "Unknown")
    paired = asm.get("paired_assembly", {}).get("assembly_accession")
    meta = {
        "bacdive_id": bacdive_id,
        "accession": acc,
        "source": src,
        "paired_asm": paired,
        "organism": org.get("organism_name"),
        "taxid": (org.get("tax_id") or None),
        "assembly_level": asm.get("assembly_level"),
        "refseq_category": asm.get("refseq_category"),
        "bioproject": (asm.get("bioproject_lineage") or [{}])[0].get("accession") if asm.get("bioproject_lineage") else None,
        "biosample": asm.get("biosample", {}).get("accession"),
        "isolate": (asm.get("isolates") or [{}])[0].get("isolate_name") if asm.get("isolates") else None,
        "release_date": asm.get("release_date"),
        "ftp_path": (asm.get("assembly_info", {}) or {}).get("ftp_path"),
        "last_checked": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "notes": None,
    }
    return meta

def add_genome(conn: sqlite3.Connection, bacdive_id: str, email: str = "you@example.org") -> Tuple[int, Optional[str]]:
    """
    Populate genome info for one strain:
    - Try GCF_/GCA_ from 'accessions'
    - Else try BioSample, then BioProject
    - Else try TaxID + isolate/culture collection (weak)
    Inserts rows into 'assemblies' and picks a primary in 'strain_genome_links'.
    Returns (n_assemblies_added, chosen_accession_or_None).
    """
    _ensure_genome_tables(conn)
    cur = conn.cursor()

    # Pull hints from your DB
    s = cur.execute("""
        SELECT bacdive_id, ncbi_taxid, biosample, bioproject, isolate, culture_collection_ids
        FROM strains WHERE bacdive_id = ?
    """, (bacdive_id,)).fetchone()
    if not s:
        raise ValueError(f"Unknown bacdive_id: {bacdive_id}")

    _, taxid, biosample, bioproject, isolate, cc_ids = s
    if cc_ids:
        # for weak text query if all else fails
        isolate = isolate or cc_ids.split(";")[0]

    # 1) Direct by accession (best)
    acc_rows = cur.execute("""
      SELECT acc_no FROM accessions WHERE bacdive_id=? ORDER BY
        CASE WHEN acc_no LIKE 'GCF_%' THEN 0 WHEN acc_no LIKE 'GCA_%' THEN 1 ELSE 2 END
    """, (bacdive_id,)).fetchall()
    accs = [a[0] for a in acc_rows]

    reports: List[Dict] = []
    tried = set()

    def add_report_from_acc(acc):
        if acc in tried:
            return
        tried.add(acc)
        rep = _search_by_accession(acc, email=email)
        if rep:
            reports.append(rep)

    for acc in accs:
        if acc.startswith(("GCF_", "GCA_")):
            add_report_from_acc(acc)

    # 2) Search by BioSample / BioProject
    if not reports and biosample:
        reports = _search_by_filter(biosample=biosample, email=email)
    if not reports and bioproject:
        reports = _search_by_filter(bioproject=bioproject, email=email)

    # 3) Weak: taxid + isolate (or culture id text)
    if not reports and taxid:
        reports = _search_by_filter(taxid=taxid, isolate=isolate, email=email)

    # Normalize & upsert
    to_insert = []
    for rep in reports:
        meta = _parse_report(bacdive_id, rep)
        if not meta.get("accession"):
            continue
        to_insert.append(meta)

    n_added = 0
    for m in to_insert:
        cur.execute("""
          INSERT INTO assemblies(bacdive_id, accession, source, paired_asm, organism, taxid,
                                 assembly_level, refseq_category, bioproject, biosample, isolate,
                                 release_date, ftp_path, last_checked, notes)
          VALUES (:bacdive_id, :accession, :source, :paired_asm, :organism, :taxid,
                  :assembly_level, :refseq_category, :bioproject, :biosample, :isolate,
                  :release_date, :ftp_path, :last_checked, :notes)
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
        """, m)
        n_added += cur.rowcount != 0

    # Choose a primary assembly for this strain
    chosen = None
    if to_insert:
        chosen = _choose_best(to_insert)
        if chosen:
            cur.execute("""
              UPDATE assemblies SET chosen=CASE WHEN accession=? THEN 1 ELSE 0 END
              WHERE bacdive_id=?
            """, (chosen, bacdive_id))
            cur.execute("""
              INSERT INTO strain_genome_links(bacdive_id, accession)
              VALUES (?,?)
              ON CONFLICT(bacdive_id) DO UPDATE SET accession=excluded.accession
            """, (bacdive_id, chosen))

    conn.commit()
    return (len(to_insert), chosen)

def add_genomes_for_missing(db_path: str, email: str = "you@example.org"):
    """
    Find all strains with no primary assembly linked yet, attempt to add one.
    """
    conn = sqlite3.connect(db_path)
    _ensure_genome_tables(conn)
    cur = conn.cursor()

    missing = cur.execute("""
      SELECT s.bacdive_id
      FROM strains s
      LEFT JOIN strain_genome_links l ON s.bacdive_id = l.bacdive_id
      WHERE l.accession IS NULL
    """).fetchall()

    print(f"{len(missing)} strains missing genome link")
    for (bd_id,) in missing:
        try:
            n, chosen = add_genome(conn, bd_id, email=email)
            print(f"[{bd_id}] added {n} assemblies; chosen={chosen}")
        except Exception as e:
            print(f"[{bd_id}] ERROR: {e}")
    conn.close()

def main(args):

    # One strain by id:
    # conn = sqlite3.connect("bacdive_traits.sqlite")
    # n, chosen = add_genome(conn, bacdive_id="12345", email="you@example.org")
    # conn.close()
    # print("added:", n, "chosen:", chosen)

    # Or all missing:
    add_genomes_for_missing("bacdive_traits.sqlite", email=args.email)


if __name__ == "__main__":
    args = parse_args()
    main(args)

