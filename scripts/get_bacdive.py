#!/usr/bin/env python3
"""
get_phenotype_info.py

Pull phenotype-rich strain records from the BacDive API and write them into a
small SQLite database that you can later enrich with genome assemblies (RefSeq/ENA).

- Reads BacDive credentials from a two-line file: email, password
- Supports multiple search modes (taxonomy, explicit BacDive IDs, etc.)
- Logs the query, number of hits, and the first three hits for quick sanity checks
- Writes two tables: strains, accessions (same schema you already used)

Usage:
    python get_phenotype_info.py \
        --info-file .bacdive_info \
        --db bacdive_traits.sqlite \
        --search-mode taxonomy \
        --taxonomy "Escherichia coli" \
        --searchtype exact \
        --limit 50 \
        --log-level INFO

Tips (guaranteed-to-hit taxonomy terms):
  - "Escherichia coli" (species)
  - "Bacillus subtilis"
  - "Staphylococcus aureus"
  - "Pseudomonas aeruginosa"
  - "Halomonas" (genus; many entries)

If you want phenotype filtering later, add --require-halophilic --require-motile.
"""

import argparse
import itertools
import json
import logging
import sqlite3
from typing import Callable, Dict, Iterable, List, Optional

import bacdive


# -------------------------
# Helpers: general parsing
# -------------------------

def _as_list(x):
    if x is None:
        return []
    if isinstance(x, list):
        return x
    return [x]


def _get_first_int(text):
    """BacDive sometimes returns '562 (species)'; return the leading integer."""
    if not text:
        return None
    s = str(text).strip()
    num = []
    for ch in s:
        if ch.isdigit():
            num.append(ch)
        else:
            break
    return int("".join(num)) if num else None


def _extract_accessions(entry: dict) -> List[Dict[str, Optional[str]]]:
    """
    Return a list of dicts: each with acc_no, database, type, notes.
    Designed to capture INSDC assembly/sequence accessions (GCA_/GCF_/ENA/DDBJ/...).
    """
    accs = []
    seqinfo = entry.get("sequence information", {}) or {}
    genome_seq = _as_list(seqinfo.get("genome sequence"))
    other_seq = _as_list(seqinfo.get("sequence entries"))

    def push(obj, note=None):
        acc_no = (obj or {}).get("acc no") or (obj or {}).get("accession")
        if not acc_no:
            return
        db = (obj or {}).get("database") or (obj or {}).get("db")
        typ = (obj or {}).get("type") or (obj or {}).get("record type")
        accs.append({
            "acc_no": str(acc_no),
            "database": (str(db) if db else None),
            "type": (str(typ) if typ else None),
            "notes": (str(note) if note else None),
        })

    for g in genome_seq:
        push(g, note="genome sequence")
    for s in other_seq:
        push(s, note="sequence entry")

    # Fallback scan (defensive)
    for k, v in seqinfo.items():
        if isinstance(v, str) and any(tok in v.upper() for tok in ("GCA_", "GCF_", "ENA", "DDBJ", "AF", "NC_")):
            for token in v.replace(",", " ").split():
                if any(token.upper().startswith(p) for p in ("GCA_", "GCF_", "NC_", "NZ_")):
                    accs.append({"acc_no": token, "database": None, "type": None, "notes": f"from {k}"})
        elif isinstance(v, list):
            for it in v:
                if isinstance(it, str) and (it.upper().startswith("GCA_") or it.upper().startswith("GCF_")):
                    accs.append({"acc_no": it, "database": None, "type": None, "notes": f"from {k}"})

    # Deduplicate (preserve order)
    out, seen = [], set()
    for a in accs:
        key = (a["acc_no"], a.get("database"), a.get("type"))
        if key not in seen:
            seen.add(key)
            out.append(a)
    return out


def _extract_collections(entry: dict) -> List[str]:
    """Return culture collection identifiers like DSM 1234, ATCC 9999."""
    cc = []
    sn = entry.get("strain number", {}) or {}
    for k, v in sn.items():
        if not v:
            continue
        ks = k.lower()
        if "culture" in ks or "collection" in ks or "strain" in ks:
            if isinstance(v, list):
                cc.extend(str(x) for x in v if x)
            elif isinstance(v, str):
                cc.append(v)
    # Normalize + dedup
    out, seen = [], set()
    for s in cc:
        s = s.strip()
        if s and s not in seen:
            seen.add(s)
            out.append(s)
    return out


def _extract_biosample_bioproject(entry: dict):
    """Heuristic scan for BioSample / BioProject identifiers if present."""
    biosample = None
    bioproject = None
    seqinfo = entry.get("sequence information", {}) or {}
    flattened = []

    def flatten(x):
        if isinstance(x, dict):
            for k, v in x.items():
                flattened.append((str(k), v))
                flatten(v)
        elif isinstance(x, list):
            for y in x:
                flatten(y)
    flatten(seqinfo)

    def find_first(prefixes):
        for k, v in flattened:
            if isinstance(v, str) and any(p.lower() in k.lower() or v.upper().startswith(p) for p in prefixes):
                return v.strip()
        return None

    biosample = find_first(["biosample", "SAMN", "SAMEA", "SAMD"])
    bioproject = find_first(["bioproject", "PRJNA", "PRJEB", "PRJDB"])
    return biosample, bioproject


def _extract_normalized(entry: dict) -> Dict[str, Optional[str]]:
    """Make a flat, stable record from a raw BacDive JSON entry."""
    tax = entry.get("taxonomy", {}) or {}
    phys = entry.get("physiology", {}) or {}
    morph = entry.get("morphology", {}) or {}

    bd_id = entry.get("strain number", {}).get("BacDive-ID")
    bd_id = str(bd_id) if bd_id is not None else None

    ncbi_taxid = _get_first_int(tax.get("NCBI Taxonomy-ID"))

    organism_name = tax.get("name and status") or tax.get("name") or None
    organism_name = str(organism_name) if organism_name else None

    gram_stain = (morph.get("Gram stain") or None)
    cell_shape = (morph.get("cell shape") or None)
    motility = (morph.get("motility") or None)
    oxygen_tol = (phys.get("oxygen tolerance") or None)
    temp_opt = (phys.get("temperature optimum") or None)
    salinity_opt = (phys.get("NaCl concentration optimum") or None)

    accs = _extract_accessions(entry)
    biosample, bioproject = _extract_biosample_bioproject(entry)

    culture_ids = _extract_collections(entry)
    isolate = None
    for src in (tax, entry.get("strain number", {}) or {}):
        for k, v in (src or {}).items():
            if isinstance(v, str) and "isolate" in k.lower():
                isolate = v
                break

    return {
        "bacdive_id": bd_id,
        "organism_name": organism_name,
        "ncbi_taxid": ncbi_taxid,
        "gram_stain": gram_stain,
        "cell_shape": cell_shape,
        "motility": motility,
        "oxygen_tolerance": oxygen_tol,
        "temp_optimum": temp_opt,
        "salinity_optimum": salinity_opt,
        "biosample": biosample,
        "bioproject": bioproject,
        "culture_collection_ids": ";".join(culture_ids) if culture_ids else None,
        "isolate": isolate,
        "accessions": accs,          # join-critical; one-to-many goes to table 'accessions'
        "raw_json": json.dumps(entry),
    }


# -------------------------
# BacDive client & fetch
# -------------------------

def start_client(info_file: str) -> bacdive.BacdiveClient:
    with open(info_file, "r") as f:
        tokens = f.readlines()
    email = tokens[0].strip()
    pw = tokens[1].strip()
    client = bacdive.BacdiveClient(email, pw)
    return client


def fetch_bacdive_records(
    client: bacdive.BacdiveClient,
    *,
    search_kwargs: Optional[Dict] = None,
    bacdive_ids: Optional[Iterable[str]] = None,
    limit: Optional[int] = None,
    include_predictions: bool = False,
    phenotype_predicate: Optional[Callable[[dict], bool]] = None,
) -> List[Dict]:
    """
    Run a search or hydrate explicit BacDive IDs, retrieve entries, optionally filter,
    and return normalized dicts ready to persist.
    """
    if include_predictions and hasattr(client, "includePredictions"):
        client.includePredictions()

    # Handle 'searchtype' passthrough (exact / contains)
    if search_kwargs and "searchtype" in {k.lower() for k in search_kwargs.keys()}:
        for k in list(search_kwargs.keys()):
            if k.lower() == "searchtype":
                st = search_kwargs.pop(k)
                client.setSearchType(st)
                break

    # Run search
    if bacdive_ids is not None:
        bacdive_ids = list(bacdive_ids)
        logging.info(f"Searching by BacDive IDs (n={len(bacdive_ids)}): first={bacdive_ids[:3]}")
        client.search(id=bacdive_ids)
        query_repr = f"id={bacdive_ids[:5]}{'...' if len(bacdive_ids)>5 else ''}"
    elif search_kwargs:
        logging.info(f"Searching with kwargs: {search_kwargs}")
        client.search(**search_kwargs)
        query_repr = repr(search_kwargs)
    else:
        raise ValueError("Provide either search_kwargs or bacdive_ids.")

    # Retrieve stream and slice
    stream = client.retrieve()
    if limit is not None:
        stream = itertools.islice(stream, int(limit))

    # Materialize results
    raw_entries = list(stream)
    logging.info(f"Search completed. Query={query_repr} → raw hits: {len(raw_entries)}")

    # Show first three quickly for sanity
    preview = []
    for e in raw_entries[:3]:
        bd_id = (e.get("strain number", {}) or {}).get("BacDive-ID")
        name = (e.get("taxonomy", {}) or {}).get("name and status") or (e.get("taxonomy", {}) or {}).get("name")
        preview.append((bd_id, name))
    if preview:
        logging.info(f"First 3 hits: {preview}")
    else:
        logging.warning("No hits returned by BacDive. Consider using a more specific taxonomy or the presets below.")

    # Apply optional phenotype filter, normalize
    records = []
    for raw in raw_entries:
        if phenotype_predicate and not phenotype_predicate(raw):
            continue
        records.append(_extract_normalized(raw))

    logging.info(f"Normalized records retained: {len(records)}")
    return records


# -------------------------
# Persistence (SQLite)
# -------------------------

def write_sqlite(db_path: str, records: List[Dict]) -> None:
    """
    Creates/updates a SQLite DB with:
      - strains(bacdive_id PK, organism_name, ncbi_taxid, ... , raw_json)
      - accessions(bacdive_id FK, acc_no, database, type, notes)
    """
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.executescript("""
    PRAGMA foreign_keys=ON;

    CREATE TABLE IF NOT EXISTS strains (
        bacdive_id TEXT PRIMARY KEY,
        organism_name TEXT,
        ncbi_taxid INTEGER,
        gram_stain TEXT,
        cell_shape TEXT,
        motility TEXT,
        oxygen_tolerance TEXT,
        temp_optimum TEXT,
        salinity_optimum TEXT,
        biosample TEXT,
        bioproject TEXT,
        culture_collection_ids TEXT,
        isolate TEXT,
        raw_json TEXT
    );

    CREATE TABLE IF NOT EXISTS accessions (
        bacdive_id TEXT,
        acc_no TEXT,
        database TEXT,
        type TEXT,
        notes TEXT,
        PRIMARY KEY (bacdive_id, acc_no),
        FOREIGN KEY (bacdive_id) REFERENCES strains(bacdive_id) ON DELETE CASCADE
    );

    CREATE INDEX IF NOT EXISTS idx_accessions_acc_no ON accessions(acc_no);
    CREATE INDEX IF NOT EXISTS idx_strains_taxid ON strains(ncbi_taxid);
    """)

    if not records:
        logging.warning("No records to write. Database unchanged.")
        conn.close()
        return

    cur.executemany("""
    INSERT INTO strains (
        bacdive_id, organism_name, ncbi_taxid, gram_stain, cell_shape, motility,
        oxygen_tolerance, temp_optimum, salinity_optimum, biosample, bioproject,
        culture_collection_ids, isolate, raw_json
    ) VALUES (
        :bacdive_id, :organism_name, :ncbi_taxid, :gram_stain, :cell_shape, :motility,
        :oxygen_tolerance, :temp_optimum, :salinity_optimum, :biosample, :bioproject,
        :culture_collection_ids, :isolate, :raw_json
    )
    ON CONFLICT(bacdive_id) DO UPDATE SET
        organism_name=excluded.organism_name,
        ncbi_taxid=excluded.ncbi_taxid,
        gram_stain=excluded.gram_stain,
        cell_shape=excluded.cell_shape,
        motility=excluded.motility,
        oxygen_tolerance=excluded.oxygen_tolerance,
        temp_optimum=excluded.temp_optimum,
        salinity_optimum=excluded.salinity_optimum,
        biosample=excluded.biosample,
        bioproject=excluded.bioproject,
        culture_collection_ids=excluded.culture_collection_ids,
        isolate=excluded.isolate,
        raw_json=excluded.raw_json
    """, records)

    acc_rows = []
    for r in records:
        bd = r.get("bacdive_id")
        for a in r.get("accessions", []) or []:
            acc_rows.append({
                "bacdive_id": bd,
                "acc_no": a.get("acc_no"),
                "database": a.get("database"),
                "type": a.get("type"),
                "notes": a.get("notes"),
            })

    if acc_rows:
        cur.executemany("""
        INSERT INTO accessions (bacdive_id, acc_no, database, type, notes)
        VALUES (:bacdive_id, :acc_no, :database, :type, :notes)
        ON CONFLICT(bacdive_id, acc_no) DO UPDATE SET
            database=excluded.database,
            type=excluded.type,
            notes=excluded.notes
        """, acc_rows)
    else:
        logging.info("No accessions found in these records.")

    conn.commit()
    conn.close()
    logging.info(f"Wrote {len(records)} strains and {len(acc_rows)} accessions to {db_path}")


# -------------------------
# CLI + presets
# -------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Get phenotypical data from BacDive API into SQLite")
    ap.add_argument("--info-file", type=str, default=".bacdive_info",
                    help="File with BacDive email (line1) and password (line2)")
    ap.add_argument("--db", type=str, default="bacdive_traits.sqlite",
                    help="SQLite DB path to create/update")
    ap.add_argument("--limit", type=int, default=50, help="Max entries to retrieve")

    # Search mode
    ap.add_argument("--search-mode", choices=["taxonomy", "ids"], default="taxonomy",
                    help="Search by taxonomy (genus/species string) or a list of BacDive IDs")
    ap.add_argument("--taxonomy", type=str, default="Escherichia coli",
                    help="Taxonomy string when search-mode=taxonomy")
    ap.add_argument("--searchtype", choices=["exact", "contains"], default="exact",
                    help="BacDive search type: exact or contains")
    ap.add_argument("--ids", type=str, nargs="*", default=None,
                    help="Explicit BacDive IDs to hydrate when search-mode=ids")

    # Optional quick phenotype flags (example)
    ap.add_argument("--require-halophilic", action="store_true", help="Keep only halophilic")
    ap.add_argument("--require-motile", action="store_true", help="Keep only motile")

    ap.add_argument("--log-level", default="INFO", help="Logging level (DEBUG, INFO, WARNING, ERROR)")
    return ap.parse_args()


def make_predicate(args) -> Optional[Callable[[dict], bool]]:
    if not (args.require_halophilic or args.require_motile):
        return None

    def pred(raw):
        phys = (raw.get("physiology") or {})
        morph = (raw.get("morphology") or {})
        ok = True
        if args.require_halophilic:
            haloph = "haloph" in (phys.get("NaCl concentration optimum", "") or "").lower() \
                     or "haloph" in (phys.get("halophily", "") or "").lower()
            ok = ok and haloph
        if args.require_motile:
            motile = (morph.get("motility") or "").strip().lower() in {"motile", "yes", "positive"}
            ok = ok and motile
        return ok
    return pred


def main():
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(levelname)s: %(message)s")

    client = start_client(args.info_file)

    if args.search_mode == "taxonomy":
        # BacDive expects taxonomy as a tuple; contains/exact controlled via setSearchType
        search_kwargs = {"taxonomy": (args.taxonomy,), "searchtype": args.searchtype}
        logging.info(f"Running taxonomy search: '{args.taxonomy}' (searchtype={args.searchtype}, limit={args.limit})")
        records = fetch_bacdive_records(
            client,
            search_kwargs=search_kwargs,
            limit=args.limit,
            phenotype_predicate=make_predicate(args),
        )
    else:
        if not args.ids:
            raise SystemExit("search-mode=ids requires --ids BacDiveID1 BacDiveID2 ...")
        logging.info(f"Hydrating explicit BacDive IDs (n={len(args.ids)})")
        records = fetch_bacdive_records(
            client,
            bacdive_ids=args.ids,
            limit=args.limit,
            phenotype_predicate=make_predicate(args),
        )

    if not records:
        logging.warning("No records returned after normalization; database will not be updated.")
        return 0

    write_sqlite(args.db, records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

