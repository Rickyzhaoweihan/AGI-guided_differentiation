"""
ebi.py — Download datasets from EBI ArrayExpress / BioStudies.

Tier 1: *.h5ad in study files  →  direct download
Tier 3: processed count matrix (TSV/CSV)  →  build AnnData → h5ad
"""

from __future__ import annotations
import gzip
import re
from pathlib import Path
from typing import Optional

import anndata as ad
import pandas as pd
import requests
import scipy.sparse

from .geo import DataUnavailableError, _download_file


BIOSTUDIES_API = "https://www.ebi.ac.uk/biostudies/api/v1"
# Direct FTP-over-HTTPS URL pattern (avoids a redirect that breaks requests)
# Pattern: E-MTAB-5061 → fire/E-MTAB-/061/E-MTAB-5061/Files/
BIOSTUDIES_FIRE = "https://ftp.ebi.ac.uk/biostudies/fire"


def _fire_url(accession: str, filename: str) -> str:
    """Build the direct ftp.ebi.ac.uk HTTPS URL for an ArrayExpress file."""
    # Extract the DB prefix and numeric part: E-MTAB-5061 → prefix=E-MTAB-, num=5061
    import re as _re
    m = _re.match(r"(E-[A-Z]+-)(\d+)", accession)
    if m:
        prefix, num = m.group(1), m.group(2)
        shard = num[-3:].zfill(3)   # last 3 digits, e.g. 5061 → 061
        return f"{BIOSTUDIES_FIRE}/{prefix}/{shard}/{accession}/Files/{filename}"
    # Fallback: go through the redirect
    return f"https://www.ebi.ac.uk/biostudies/files/{accession}/{filename}"


def _list_study_files(accession: str) -> list[dict]:
    """
    Query BioStudies REST API and return all files for the study.
    Parses the study JSON's nested section[].files[] structure.
    Each returned item: {name, url}
    """
    url = f"{BIOSTUDIES_API}/studies/{accession}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    files = []
    _collect_files(data, accession, files)
    return files


def _collect_files(node, accession: str, out: list):
    """Recursively collect file entries from BioStudies study JSON."""
    if isinstance(node, list):
        for item in node:
            _collect_files(item, accession, out)
    elif isinstance(node, dict):
        # A file node has 'path' and 'type': 'file'
        if node.get("type") == "file" and "path" in node:
            fname = node["path"]
            out.append({"name": fname, "url": _fire_url(accession, fname)})
        # Also handle nested 'files' lists
        if "files" in node and isinstance(node["files"], list):
            for f in node["files"]:
                if isinstance(f, dict) and f.get("type") == "file" and "path" in f:
                    fname = f["path"]
                    out.append({"name": fname, "url": _fire_url(accession, fname)})
        for v in node.values():
            if isinstance(v, (dict, list)):
                _collect_files(v, accession, out)


# ── Tier 1 ────────────────────────────────────────────────────────────────────

def _try_tier1(files: list[dict], out_dir: Path, accession: str) -> Optional[Path]:
    h5_files = [f for f in files if re.search(r"\.(h5ad|h5)$", f["name"], re.I)]
    if not h5_files:
        return None
    h5_files.sort(key=lambda f: (0 if f["name"].endswith(".h5ad") else 1))
    f = h5_files[0]
    dest = out_dir / f"{accession}.h5ad"
    _download_file(f["url"], dest)
    return dest


# ── Tier 3: TSV count matrix ──────────────────────────────────────────────────

def _try_tier3(files: list[dict], out_dir: Path, accession: str) -> Optional[Path]:
    tsv_files = [
        f for f in files
        if re.search(r"\.(tsv|csv|txt)(\.gz)?$", f["name"], re.I)
        and re.search(r"(count|expression|matrix|processed|rpkm|fpkm|tpm)", f["name"], re.I)
        and not re.search(r"\.(idf|sdrf)\.", f["name"], re.I)
    ]
    if not tsv_files:
        return None

    f = tsv_files[0]
    url = f["url"]
    suffix = ".gz" if f["name"].endswith(".gz") else ""
    tmp = out_dir / f"_tmp_{accession}{suffix}"
    if not tmp.exists():
        _download_file(url, tmp)
    else:
        print(f"  Reusing cached tmp file: {tmp.name} ({tmp.stat().st_size // 1024 // 1024} MB)")

    opener = gzip.open if suffix else open
    sep = "\t" if not f["name"].lower().endswith(".csv") else ","

    # Sniff first two lines to detect format:
    #   Standard (N header cols, N+1 data cols): index in col 0
    #   Two-ID rows (N header cols, N+2 data cols): col 0 = gene name, col 1 = transcript ID (drop)
    #   No index (N header cols, N data cols): set col 0 as index
    with opener(tmp, "rt") as fh:
        header_line = fh.readline().rstrip("\n")
        data_line   = fh.readline().rstrip("\n")
    header_len = len(header_line.split(sep))
    data_len   = len(data_line.split(sep))
    extra = data_len - header_len   # 0 = no index, 1 = standard, 2 = gene+transcript

    with opener(tmp, "rt") as fh:
        df = pd.read_csv(fh, sep=sep, header=0, index_col=False)

    if extra == 2:
        # Two gene-ID columns: use col 0 as gene name, drop col 1 (transcript ID)
        df = df.set_index(df.columns[0]).drop(columns=[df.columns[1]])
    elif extra >= 1:
        df = df.set_index(df.columns[0])
    else:
        df = df.set_index(df.columns[0])

    # Drop any remaining non-numeric columns
    df = df.select_dtypes(include="number")

    # Genes × cells → transpose to cells × genes
    if df.shape[0] > df.shape[1]:
        df = df.T

    adata = ad.AnnData(X=scipy.sparse.csr_matrix(df.values, dtype="float32"))
    adata.obs_names = df.index.tolist()
    adata.var_names = df.columns.tolist()

    dest = out_dir / f"{accession}.h5ad"
    adata.write_h5ad(dest)
    print(f"  Saved {adata.shape[0]} cells × {adata.shape[1]} genes → {dest}")
    tmp.unlink(missing_ok=True)
    return dest


# ── public entry point ────────────────────────────────────────────────────────

def download(accession: str, out_dir: str | Path) -> Path:
    """
    Download an EBI ArrayExpress / BioStudies dataset.
    Tries Tier 1 (h5ad) → Tier 3 (processed TSV).
    Raises DataUnavailableError if nothing suitable found.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dest = out_dir / f"{accession}.h5ad"
    if dest.exists():
        print(f"  Already downloaded: {dest}")
        return dest

    print(f"\n[EBI] {accession}")
    files = _list_study_files(accession)
    print(f"  Found {len(files)} file(s) in study")

    for f in files[:20]:   # show first 20
        print(f"    {f['name']}")
    if len(files) > 20:
        print(f"    ... and {len(files) - 20} more")

    result = _try_tier1(files, out_dir, accession)
    if result:
        print(f"  [Tier 1] h5ad downloaded: {result}")
        return result

    result = _try_tier3(files, out_dir, accession)
    if result:
        print(f"  [Tier 3] count matrix reconstructed: {result}")
        return result

    raise DataUnavailableError(
        f"{accession}: no h5ad or processed count matrix found. "
        f"Files may require manual download from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}"
    )
