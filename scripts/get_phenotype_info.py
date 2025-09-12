import argparse

import bacdive

import json
import sqlite3
import itertools
from typing import Callable, Dict, Iterable, List, Optional
# import requests
# import time
# import datetime
# import re


def start_client(info_file):

    with open(info_file, "r") as f:
        tokens = f.readlines()

    email = tokens[0][:-1]
    pw = tokens[1][:-1]

    client = bacdive.BacdiveClient(email, pw)

    return client


def parse_args():

    ap = argparse.ArgumentParser(description="Get phenotypical data from BacDive API")
    ap.add_argument("--info_file", type=str, default= ".bacdive_info", help="Info file containing BacDive email in the first line, and password in the second line")

    return ap.parse_args()


# helper functions to structure data from phenotype search

def _as_list(x):
    if x is None:
        return []
    if isinstance(x, list):
        return x
    return [x]


def _get_first_int(text):
    if not text:
        return None
    # BacDive sometimes gives "562 (species)" – keep leading digits
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
    Return a list of dicts: each with acc_no, db, type, notes.
    Designed to capture INSDC assembly or sequence accessions.
    """
    accs = []

    # Sequence/genome information
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

    # Sometimes BacDive tucks accessions in a flat list/dict; be defensive:
    for k, v in seqinfo.items():
        if isinstance(v, str) and any(tok in v.upper() for tok in ("GCA_", "GCF_", "ENA", "DDBJ", "AF", "NC_")):
            # naive token split, last-resort
            for token in v.replace(",", " ").split():
                if any(token.upper().startswith(p) for p in ("GCA_", "GCF_", "NC_", "NZ_")):
                    accs.append({"acc_no": token, "database": None, "type": None, "notes": f"from {k}"})
        elif isinstance(v, list):
            for it in v:
                if isinstance(it, str) and (it.upper().startswith("GCA_") or it.upper().startswith("GCF_")):
                    accs.append({"acc_no": it, "database": None, "type": None, "notes": f"from {k}"})

    # Deduplicate while preserving order
    seen = set()
    dedup = []
    for a in accs:
        key = (a["acc_no"], a.get("database"), a.get("type"))
        if key not in seen:
            seen.add(key)
            dedup.append(a)
    return dedup



def _extract_collections(entry: dict) -> List[str]:
    """Return culture collection identifiers like DSM 1234, ATCC 9999."""
    cc = []
    sn = entry.get("strain number", {}) or {}
    # BacDive often uses keys like "culture collection" or per-collection fields.
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
    """
    Heuristic scan for BioSample / BioProject identifiers if present in sequence info.
    """
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

    # Core IDs
    bd_id = entry.get("strain number", {}).get("BacDive-ID")
    bd_id = str(bd_id) if bd_id is not None else None

    ncbi_taxid = _get_first_int(tax.get("NCBI Taxonomy-ID"))

    # Name/status
    organism_name = tax.get("name and status") or tax.get("name") or None
    organism_name = str(organism_name) if organism_name else None

    # Phenotype highlights (extend as you like)
    gram_stain = (morph.get("Gram stain") or None)
    cell_shape = (morph.get("cell shape") or None)
    motility = (morph.get("motility") or None)
    oxygen_tol = (phys.get("oxygen tolerance") or None)
    temp_opt = (phys.get("temperature optimum") or None)
    salinity_opt = (phys.get("NaCl concentration optimum") or None)

    # Accessions + sample/project
    accs = _extract_accessions(entry)
    biosample, bioproject = _extract_biosample_bioproject(entry)

    # Culture collections & strain labels (good fuzzy join aids)
    culture_ids = _extract_collections(entry)
    isolate = None
    # Some entries expose isolate in taxonomy or strain number
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
        # Join-critical: keep *all* accessions (GCA_/GCF_/ENA/DDBJ/etc.) separately
        "accessions": accs,
        "raw_json": json.dumps(entry),
    }


def fetch_bacdive_records(
    client,
    *,
    # Option A: broad search (taxonomy, genome, etc.) via client.search(...)
    search_kwargs: Optional[Dict] = None,
    # Option B: hydrate explicit BacDive IDs (recommended after Advanced Search CSV)
    bacdive_ids: Optional[Iterable[str]] = None,
    limit: Optional[int] = None,
    include_predictions: bool = False,
    phenotype_predicate: Optional[Callable[[dict], bool]] = None,
) -> List[Dict]:
    """
    Connects to BacDive, runs a search or ID hydration, retrieves entries,
    applies optional client-side phenotype filtering, and returns normalized dicts.

    Either provide `search_kwargs` (e.g., {'taxonomy': ('Halomonas',), 'searchtype':'exact'})
    or provide `bacdive_ids` (from Advanced Search export). If both are given, IDs win.

    phenotype_predicate receives the *raw* entry and should return True/False.
    """
    if include_predictions and hasattr(client, "includePredictions"):
        client.includePredictions()

    # If the caller passed a 'searchtype' key, reflect it via setSearchType
    if search_kwargs and "searchtype" in {k.lower() for k in search_kwargs.keys()}:
        # normalize key access
        for k in list(search_kwargs.keys()):
            if k.lower() == "searchtype":
                client.setSearchType(search_kwargs.pop(k))
                break

    if bacdive_ids is not None:
        # Hydrate by explicit BacDive IDs
        client.search(id=list(bacdive_ids))
    elif search_kwargs:
        client.search(**search_kwargs)
    else:
        raise ValueError("Provide either search_kwargs or bacdive_ids.")

    # Retrieve stream and slice
    stream = client.retrieve()
    if limit is not None:
        stream = itertools.islice(stream, int(limit))

    records = []
    for raw in stream:
        if phenotype_predicate and not phenotype_predicate(raw):
            continue
        records.append(_extract_normalized(raw))
    return records



def write_sqlite(db_path: str, records: List[Dict]) -> None:
    """
    Creates/updates a SQLite DB with:
      - strains(bacdive_id PK, organism_name, ncbi_taxid, ... , raw_json)
      - accessions(bacdive_id FK, acc_no, database, type, notes)

    Idempotent upserts keyed by bacdive_id + acc_no.
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

    # Upsert strains
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

    # Upsert accessions
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

    conn.commit()
    conn.close()



def main(args):

    client = start_client(args.info_file)

    # 1) Pull 50 strains by genus, keep only halophilic & motile (client-side filter),
    #    then write to SQLite.
    def is_halophilic_and_motile(raw):
        phys = (raw.get("physiology") or {})
        morph = (raw.get("morphology") or {})
        haloph = "haloph" in (phys.get("NaCl concentration optimum","") or "").lower() \
                 or "haloph" in (phys.get("halophily","") or "").lower()
        motile = (morph.get("motility") or "").strip().lower() in {"motile","yes","positive"}
        return haloph and motile

    recs = fetch_bacdive_records(
        client,
        search_kwargs={"taxonomy": ("Escherichia",), "searchtype": "contains"},
        limit=50,
        # phenotype_predicate=is_halophilic_and_motile,
    )

    write_sqlite("bacdive_traits.sqlite", recs)

    return 0


if __name__ == "__main__":
    args = parse_args()
    main(args)
