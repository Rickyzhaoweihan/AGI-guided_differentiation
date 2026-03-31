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


def _list_study_files(accession: str) -> list[dict]:
    """
    Query BioStudies REST API and return all files for the study.
    Each item: {name, url, type}
    """
    url = f"{BIOSTUDIES_API}/studies/{accession}/files/tree"
    resp = requests.get(url, timeout=30)
    if not resp.ok:
        # Fallback: try the flat file list
        url = f"{BIOSTUDIES_API}/studies/{accession}/files"
        resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    files = []
    _walk(data, files)
    return files


def _walk(node, out: list):
    """Recursively walk BioStudies file tree nodes."""
    if isinstance(node, list):
        for item in node:
            _walk(item, out)
    elif isinstance(node, dict):
        if "name" in node and "url" in node:
            out.append({"name": node["name"], "url": node["url"]})
        for v in node.values():
            if isinstance(v, (dict, list)):
                _walk(v, out)


def _biostudies_download_url(accession: str, filename: str) -> str:
    """Build direct download URL for a BioStudies file."""
    return f"https://ftp.ebi.ac.uk/biostudies/arrayexpress/data/experiment/{accession[:7]}/{accession}/{filename}"


# ── Tier 1 ────────────────────────────────────────────────────────────────────

def _try_tier1(files: list[dict], out_dir: Path, accession: str) -> Optional[Path]:
    h5_files = [f for f in files if re.search(r"\.(h5ad|h5)$", f["name"], re.I)]
    if not h5_files:
        return None
    h5_files.sort(key=lambda f: (0 if f["name"].endswith(".h5ad") else 1))
    f = h5_files[0]
    dest = out_dir / f"{accession}.h5ad"
    url = f.get("url") or _biostudies_download_url(accession, f["name"])
    _download_file(url, dest)
    return dest


# ── Tier 3: TSV count matrix ──────────────────────────────────────────────────

def _try_tier3(files: list[dict], out_dir: Path, accession: str) -> Optional[Path]:
    tsv_files = [
        f for f in files
        if re.search(r"(count|expression|matrix|processed).*\.(tsv|csv|txt)(\.gz)?$", f["name"], re.I)
        and not re.search(r"(meta|sample|cell|barcode|annotation)", f["name"], re.I)
    ]
    if not tsv_files:
        return None

    f = tsv_files[0]
    url = f.get("url") or _biostudies_download_url(accession, f["name"])
    suffix = ".gz" if f["name"].endswith(".gz") else ""
    tmp = out_dir / f"_tmp_{accession}{suffix}"
    _download_file(url, tmp)

    opener = gzip.open if suffix else open
    sep = "\t" if not f["name"].lower().endswith(".csv") else ","
    with opener(tmp, "rt") as fh:
        df = pd.read_csv(fh, sep=sep, index_col=0)

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
