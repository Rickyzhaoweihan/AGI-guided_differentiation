"""
catalog.py — Parse a CSV (or xlsx) dataset catalog and classify each accession.

Supported columns (detected by name, case-insensitive):
  - title / manuscript / manuscript/database
  - pmid
  - accession / database availability / data resource
  - url / link to data / url link to database
  - data type
  - tissue type
  - status / age

Handles multi-accession rows (e.g. "GSE158702, GSE114374" or "E-MTAB-9543, E-MTAB-9536").
"""

from __future__ import annotations
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ── accession type classification ────────────────────────────────────────────

DB_PATTERNS = [
    ("GEO",     r"\bGSE\d+\b",                       True),
    ("EBI",     r"\bE-[A-Z]+-\d+\b",                 True),
    ("HUBMAP",  r"\bHBM[A-Z0-9.]+\b",                True),
    ("S3",      r"s3://|s3\.console\.aws",            True),
    ("EGA",     r"\bEGA[SD]\d+\b",                    False),
    ("CNCB",    r"\bHRA\d+\b",                        False),
    ("DBGAP",   r"\bphs\d+\b",                        False),
    ("PANCDB",  r"hpap\.pmacs\.upenn\.edu|PANC\s*DB", False),
    ("DUOS",    r"\bDUOS-\d+\b",                      False),
]


def classify(raw: str) -> tuple[str, bool]:
    """Return (db_type, is_publicly_accessible) for a raw accession/URL string."""
    for db_type, pattern, public in DB_PATTERNS:
        if re.search(pattern, raw, re.IGNORECASE):
            return db_type, public
    return "UNKNOWN", False


def extract_accessions(raw: str) -> list[str]:
    """
    Split a cell that may contain multiple accessions separated by commas,
    semicolons, or newlines. Clean up trailing punctuation.
    """
    if not raw or (isinstance(raw, float)):
        return []
    raw = str(raw).strip()
    parts = re.split(r"[,;\n]+", raw)
    accessions = []
    for p in parts:
        p = p.strip().rstrip(".)").strip()
        if p:
            accessions.append(p)
    return accessions


# ── data classes ──────────────────────────────────────────────────────────────

@dataclass
class DatasetRecord:
    title: str
    pmid: Optional[str]
    accession: str          # single accession (one record per accession)
    url: Optional[str]
    db_type: str
    accessible: bool        # True = publicly downloadable
    tissue_type: Optional[str] = None
    data_type: Optional[str] = None
    status: Optional[str] = None   # embryo / fetal / adult / in vitro
    raw_accession_cell: str = ""   # original cell value before splitting


# ── column name normalisation ─────────────────────────────────────────────────

_COL_ALIASES = {
    "title":       ["title", "manuscript", "manuscript/database"],
    "pmid":        ["pmid"],
    "accession":   ["accession", "database availability", "data resource"],
    "url":         ["url", "link to data", "url link to database"],
    "data_type":   ["data type"],
    "tissue_type": ["tissue type"],
    "status":      ["status"],
}


def _find_col(columns: list[str], field: str) -> Optional[str]:
    """Return the actual column name matching a logical field."""
    aliases = _COL_ALIASES.get(field, [])
    for col in columns:
        if col.strip().lower() in aliases:
            return col
    return None


# ── main parser ───────────────────────────────────────────────────────────────

def load_catalog(path: str | Path) -> list[DatasetRecord]:
    """
    Load a CSV or xlsx catalog file and return one DatasetRecord per accession.
    """
    import pandas as pd

    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".csv":
        df = pd.read_csv(path, dtype=str)
    elif suffix in (".xlsx", ".xls"):
        df = pd.read_excel(path, dtype=str)
    else:
        raise ValueError(f"Unsupported file type: {suffix}. Use .csv or .xlsx")

    # Normalise column names for lookup
    cols = df.columns.tolist()
    col = {f: _find_col(cols, f) for f in _COL_ALIASES}

    records: list[DatasetRecord] = []

    for _, row in df.iterrows():
        def get(field):
            c = col.get(field)
            if c and c in row:
                v = row[c]
                return None if (isinstance(v, float) and str(v) == "nan") else str(v).strip()
            return None

        raw_acc = get("accession") or ""
        raw_url = get("url") or ""

        accessions = extract_accessions(raw_acc)
        urls = extract_accessions(raw_url)

        if not accessions:
            # Try to extract an accession from the URL itself
            accessions = extract_accessions(raw_url)
            if not accessions:
                continue

        for i, acc in enumerate(accessions):
            # Try to pair with a URL at the same index, else use the first URL
            url = urls[i] if i < len(urls) else (urls[0] if urls else None)

            # Classify using accession first, fall back to URL
            db_type, accessible = classify(acc)
            if db_type == "UNKNOWN" and url:
                db_type, accessible = classify(url)

            records.append(DatasetRecord(
                title=get("title") or "",
                pmid=get("pmid"),
                accession=acc,
                url=url,
                db_type=db_type,
                accessible=accessible,
                tissue_type=get("tissue_type"),
                data_type=get("data_type"),
                status=get("status"),
                raw_accession_cell=raw_acc,
            ))

    return records


def print_summary(records: list[DatasetRecord]) -> None:
    """Print a human-readable summary table."""
    public  = [r for r in records if r.accessible]
    private = [r for r in records if not r.accessible]

    print(f"\nTotal accessions parsed: {len(records)}")
    print(f"  Publicly accessible:  {len(public)}")
    print(f"  Controlled access:    {len(private)}\n")

    print(f"{'Accession':<25} {'DB':<8} {'Access':<10} Title")
    print("-" * 90)
    for r in records:
        flag = "public" if r.accessible else "CONTROLLED"
        title = (r.title[:45] + "...") if len(r.title) > 48 else r.title
        print(f"{r.accession:<25} {r.db_type:<8} {flag:<10} {title}")
