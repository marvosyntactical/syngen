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
from pprint import pprint

import bacdive


# -------------------------
# Helpers: general parsing
# -------------------------


def _to_text(v):
    """Coerce BacDive values (str/int/float/list/dict/None) into a readable string or None."""
    if v is None:
        return None
    if isinstance(v, (str, int, float)):
        return str(v)

    if isinstance(v, list):
        parts = []
        for item in v:
            if isinstance(item, (str, int, float)):
                parts.append(str(item))
            elif isinstance(item, dict):
                # Prefer common human fields if present
                for key in ("value", "result", "activity", "text", "name", "label"):
                    if key in item and item[key] not in (None, ""):
                        parts.append(str(item[key]))
                        break
                else:
                    parts.append(json.dumps(item, ensure_ascii=False))
            else:
                parts.append(str(item))
        return "; ".join(p for p in parts if p)

    if isinstance(v, dict):
        # Prefer a representative value if present, else compact JSON
        for key in ("value", "result", "activity", "text", "name", "label"):
            if key in v and v[key] not in (None, ""):
                return str(v[key])
        return json.dumps(v, ensure_ascii=False)

    return str(v)

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

# ~~~~~~~~~~~~~ BACDIVE EXTRACTION ~~~~~~~~~~~~~~


def _get(d, *path, default=None):
    """Safe nested get: _get(x, 'A','B','C') -> x['A']['B']['C'] or default."""
    cur = d
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur



def _extract_bacdive_id(entry: dict) -> Optional[str]:
    """
    Robustly extract BacDive-ID from either:
      - New 'sectioned' schema: entry['General']['BacDive-ID']
      - Older schema: entry['strain number']['BacDive-ID']
      - Fallbacks: top-level variants
    """
    # New schema
    v = _get(entry, "General", "BacDive-ID")
    if v not in (None, ""):
        return str(v)

    # Older schema
    sn = entry.get("strain number") or entry.get("strainNumber") or {}
    for k in ("BacDive-ID", "BacDive ID", "BacDiveID", "bacdive_id", "BacDive Id"):
        if isinstance(sn, dict) and sn.get(k) not in (None, ""):
            return str(sn.get(k))

    # Rare fallbacks
    for k in ("BacDive-ID", "entry number", "entry_number"):
        v = entry.get(k)
        if v not in (None, ""):
            return str(v)

    return None


def _extract_ncbi_taxid(entry: dict) -> Optional[int]:
    """
    New schema: entry['General']['NCBI tax id']['NCBI tax id'] is an int.
    Older schema: entry['taxonomy']['NCBI Taxonomy-ID'] may be '562 (species)'.
    """
    v = _get(entry, "General", "NCBI tax id", "NCBI tax id")
    if isinstance(v, int):
        return v
    if v and str(v).isdigit():
        return int(v)

    old = _get(entry, "taxonomy", "NCBI Taxonomy-ID")
    if old:
        s = str(old).strip()
        num = []
        for ch in s:
            if ch.isdigit():
                num.append(ch)
            else:
                break
        if num:
            return int("".join(num))
    return None



def _extract_names(entry: dict) -> Dict[str, Optional[str]]:
    """
    New schema: use 'Name and taxonomic classification'
    Older schema: use 'taxonomy'
    """
    # New schema
    nm = entry.get("Name and taxonomic classification") or {}
    species = nm.get("species")
    full_name = nm.get("full scientific name")
    genus = nm.get("genus")
    strain = nm.get("strain designation")

    if species or full_name or genus:
        return {
            "organism_name": full_name or species or genus,
            "species": species,
            "genus": genus,
            "strain_designation": strain,
        }

    # Older schema
    tax = entry.get("taxonomy") or {}
    return {
        "organism_name": tax.get("name and status") or tax.get("name"),
        "species": tax.get("species"),
        "genus": tax.get("genus"),
        "strain_designation": None,
    }

def _extract_culture_collections(entry: dict) -> List[str]:
    """
    New schema: 'External links' -> 'culture collection no.' (string like 'DSM 4599')
    Older schema: under 'strain number' various collection fields.
    """
    out = []
    v = _get(entry, "External links", "culture collection no.")
    if isinstance(v, str) and v.strip():
        out.append(v.strip())

    sn = entry.get("strain number") or {}
    for k, val in (sn.items() if isinstance(sn, dict) else []):
        ks = str(k).lower()
        if any(tok in ks for tok in ("culture", "collection", "strain")):
            if isinstance(val, str) and val.strip():
                out.append(val.strip())
            elif isinstance(val, list):
                out.extend([str(x).strip() for x in val if x])

    # dedup, keep order
    seen, dedup = set(), []
    for s in out:
        if s and s not in seen:
            seen.add(s)
            dedup.append(s)
    return dedup

def _extract_accessions(entry: dict) -> List[Dict[str, Optional[str]]]:
    """
    Collect INSDC accessions (GCF_/GCA_/NZ_/NC_...) from the 'Sequence information'
    section in the new schema or the older 'sequence information'.
    Be defensive: scan strings for tokens too.
    """
    accs = []

    def push(acc_no, database=None, typ=None, note=None):
        if not acc_no:
            return
        accs.append({
            "acc_no": str(acc_no),
            "database": (str(database) if database else None),
            "type": (str(typ) if typ else None),
            "notes": (str(note) if note else None),
        })

    # New schema
    seqinfo = entry.get("Sequence information") or {}
    if isinstance(seqinfo, dict):
        # Common fields we’ve seen: 'genome sequence', 'sequence entries' (when present)
        for k, v in seqinfo.items():
            if isinstance(v, dict):
                # look for typical keys
                acc = v.get("acc no") or v.get("accession") or v.get("assembly accession")
                if acc:
                    push(acc, database=v.get("database"), typ=v.get("type"), note=k)
                # also scan stringy subfields for GCA/GCF/NC_
                for subv in v.values():
                    if isinstance(subv, str):
                        _scan_string_for_accessions(subv, accs, note=k)
            elif isinstance(v, list):
                for obj in v:
                    if isinstance(obj, dict):
                        acc = obj.get("acc no") or obj.get("accession") or obj.get("assembly accession")
                        if acc:
                            push(acc, database=obj.get("database"), typ=obj.get("type"), note=k)
                        else:
                            for subv in obj.values():
                                if isinstance(subv, str):
                                    _scan_string_for_accessions(subv, accs, note=k)
                    elif isinstance(obj, str):
                        _scan_string_for_accessions(obj, accs, note=k)
            elif isinstance(v, str):
                _scan_string_for_accessions(v, accs, note=k)

    # Older schema
    old = entry.get("sequence information") or {}
    if isinstance(old, dict):
        for k, v in old.items():
            if isinstance(v, dict):
                acc = v.get("acc no") or v.get("accession")
                if acc:
                    push(acc, database=v.get("database"), typ=v.get("type"), note=k)
            elif isinstance(v, list):
                for obj in v:
                    if isinstance(obj, dict):
                        acc = obj.get("acc no") or obj.get("accession")
                        if acc:
                            push(acc, database=obj.get("database"), typ=obj.get("type"), note=k)
                    elif isinstance(obj, str):
                        _scan_string_for_accessions(obj, accs, note=k)
            elif isinstance(v, str):
                _scan_string_for_accessions(v, accs, note=k)

    # Dedup while preserving order
    seen, dedup = set(), []
    for a in accs:
        key = (a["acc_no"], a.get("database"), a.get("type"))
        if key not in seen:
            seen.add(key)
            dedup.append(a)
    return dedup

def _scan_string_for_accessions(s: str, out: List[Dict[str, Optional[str]]], note=None):
    import re
    for m in re.findall(r'(GC[AF]_\d+\.\d+|N[CTZ]_\d+(\.\d+)?)', s):
        acc = m[0]
        out.append({"acc_no": acc, "database": None, "type": None, "notes": note})

def _extract_basic_phenotypes(entry: dict) -> Dict[str, Optional[str]]:
    """
    Very light phenotype extraction that works on both schemas:
      - Gram stain, motility, oxygen tolerance if present
      - Plus a few metabolic markers from the new 'Physiology and metabolism'
    """
    morph = entry.get("Morphology") or entry.get("morphology") or {}
    phys = entry.get("Physiology and metabolism") or entry.get("physiology") or {}

    gram = morph.get("Gram stain") or morph.get("Gram-stain")
    motility = morph.get("motility") or morph.get("Motility")
    oxygen = phys.get("oxygen tolerance") or phys.get("Oxygen tolerance")

    # Simple metabolic flags (new schema example shows API 20E & enzymes)
    indole = None
    vp = None
    beta_gal = None

    # API 20E list of dicts
    for rec in phys.get("API 20E", []) or []:
        if isinstance(rec, dict):
            # keys like 'IND' for indole, 'VP' for Voges-Proskauer, 'ONPG' ~ beta-gal
            indole = indole or rec.get("IND")
            vp = vp or rec.get("VP")
            beta_gal = beta_gal or rec.get("ONPG")

    # Enzymes list (look for beta-galactosidase)
    for rec in phys.get("enzymes", []) or []:
        if isinstance(rec, dict) and str(rec.get("value","")).lower() == "beta-galactosidase":
            beta_gal = beta_gal or rec.get("activity")

    return {
        "gram_stain": gram,
        "motility": motility,
        "oxygen_tolerance": oxygen,
        "indole_test": indole,
        "vp_test": vp,
        "beta_galactosidase": beta_gal,
    }

def _extract_normalized(entry: dict) -> Dict[str, Optional[str]]:
    """Normalize a BacDive entry (new or old schema) into our flat row dict."""
    bd_id = _extract_bacdive_id(entry)
    names = _extract_names(entry)
    taxid = _extract_ncbi_taxid(entry)
    accessions = _extract_accessions(entry)
    phen = _extract_basic_phenotypes(entry)
    culture_ids = _extract_culture_collections(entry)
    morph = entry.get("Morphology") or entry.get("morphology") or {}
    phys  = entry.get("Physiology and metabolism") or entry.get("physiology") or {}

    rec = {
        "bacdive_id": bd_id,
        "organism_name": names.get("organism_name"),
        "ncbi_taxid": taxid,

        # Coerced to scalars:
        "gram_stain": _to_text(morph.get("Gram stain") or morph.get("Gram-stain")),
        "cell_shape": _to_text(morph.get("cell shape") or morph.get("Cell shape")),
        "motility": _to_text(morph.get("motility") or morph.get("Motility")),
        "oxygen_tolerance": _to_text(phys.get("oxygen tolerance") or phys.get("Oxygen tolerance")),
        "temp_optimum": _to_text(phys.get("temperature optimum") or phys.get("Temperature optimum")),
        "salinity_optimum": _to_text(phys.get("NaCl concentration optimum") or phys.get("Salt requirement")),

        "biosample": None,
        "bioproject": None,
        "culture_collection_ids": ";".join(_extract_culture_collections(entry)) or None,
        "isolate": names.get("strain_designation"),

        "accessions": accessions,           # keep as list of dicts; handled separately for insert
        "raw_json": json.dumps(entry),
    }
    return rec


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
        cls_info = e['Name and taxonomic classification']
        # print(cls_info)
        # logging.debug(f"Entry: {cls_info['full scientific name'], cls_info['type strain']}")
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

    print(records)

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

    good, dropped = [], 0
    for r in records:
        if r.get("bacdive_id"):
            good.append(r)
        else:
            dropped += 1
    if dropped:
        logging.warning(f"Dropping {dropped} entries with missing BacDive-ID.")
    records = good

    acc_rows = []
    for r in records:
        print("Accessions for this record:", r["accessions"])
        print("Ncbi taxid:", r["ncbi_taxid"])
        print(r.keys())
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
    ap.add_argument("--taxonomy", type=str, default="Acetobacter aceti",
                    help="Taxonomy string when search-mode=taxonomy")
    ap.add_argument("--searchtype", choices=["exact", "contains"], default="contains",
                    help="BacDive search type: exact or contains")
    ap.add_argument("--ids", type=str, nargs="*", default=None,
                    help="Explicit BacDive IDs to hydrate when search-mode=ids")

    # Optional quick phenotype flags (example)
    ap.add_argument("--require-halophilic", action="store_true", help="Keep only halophilic")
    ap.add_argument("--require-motile", action="store_true", help="Keep only motile")

    ap.add_argument("--log-level", default="INFO", help="Logging level (DEBUG, INFO, WARNING, ERROR)")
    ap.add_argument("--dry-run", action="store_true", help="If set, don't write to sqlite3 DB.")
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
        search_kwargs = {"taxonomy": args.taxonomy} # , "searchtype": args.searchtype}
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

    if not args.dry_run:
        write_sqlite(args.db, records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

