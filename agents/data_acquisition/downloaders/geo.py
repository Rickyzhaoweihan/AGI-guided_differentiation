"""
geo.py — Download datasets from NCBI GEO.

Tier 1: Supplementary *.h5ad / *.h5 file  →  direct download
Tier 2: 10x matrix.mtx.gz + barcodes + features  →  reconstruct AnnData
Tier 2b: *_count_matrix.txt.gz (tab-sep DGE)  →  reconstruct AnnData
Tier 4: FASTQ-only  →  skip, raise DataUnavailableError
"""

from __future__ import annotations
import gzip
import io
import re
import time
import urllib.request
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import requests
import scipy.io
import scipy.sparse


class DataUnavailableError(Exception):
    pass


# ── NCBI helpers ──────────────────────────────────────────────────────────────

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def _geo_ftp_dir(gse: str) -> str:
    """Return the FTP URL for a GSE accession's supplementary directory."""
    prefix = gse[:-(len(gse) - 3)] + "nnn"   # e.g. GSE158702 → GSE158nnn
    return f"{GEO_FTP}/{prefix}/{gse}/suppl/"


def _list_suppl_files(gse: str) -> list[dict]:
    """
    Use NCBI eutils to get the list of supplementary files for a GSE.
    Returns list of {name, url} dicts.
    """
    # NCBI rate limit: 3 req/sec without API key — add small delay
    time.sleep(0.4)

    # Search for the GEO series
    search_url = (
        f"{EUTILS}/esearch.fcgi?db=gds&term={gse}[Accession]&retmode=json"
    )
    resp = requests.get(search_url, timeout=30)
    resp.raise_for_status()
    ids = resp.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    # Fetch soft summary (much smaller than full SOFT)
    fetch_url = (
        f"{EUTILS}/esummary.fcgi?db=gds&id={ids[0]}&retmode=json"
    )
    resp = requests.get(fetch_url, timeout=30)
    resp.raise_for_status()
    result = resp.json().get("result", {})
    uid = ids[0]
    record = result.get(uid, {})

    files = []
    # FTPLink is a base path; supplementary files are listed in ftplink
    ftp_link = record.get("ftplink", "")
    if ftp_link:
        suppl_url = ftp_link.rstrip("/") + "/suppl/"
        # Convert ftp:// to https://
        suppl_url = suppl_url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
        try:
            r = requests.get(suppl_url, timeout=30)
            if r.ok:
                # Parse the directory listing — keep only bare filenames
                # (no slashes, not full URLs, not parent-dir links)
                names = re.findall(r'href="([^"/][^"]*)"', r.text)
                for name in names:
                    if name.startswith("http") or name.startswith("?"):
                        continue
                    name = name.rstrip("/")
                    if name:
                        files.append({"name": name, "url": suppl_url + name})
        except Exception:
            pass

    return files


def _download_file(url: str, dest: Path, chunk_size: int = 1024 * 1024) -> Path:
    """Stream-download url to dest, show progress. Falls back to wget on proxy errors."""
    import subprocess, shutil
    print(f"  Downloading {url.split('/')[-1]} ...")

    def _via_requests():
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            downloaded = 0
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded / total * 100
                        print(f"\r  {pct:.1f}%  ({downloaded // 1024 // 1024} MB)", end="", flush=True)
        print()

    def _via_wget():
        wget = shutil.which("wget")
        if not wget:
            return False
        result = subprocess.run(
            [wget, "-q", "--show-progress", "-O", str(dest), url],
            timeout=7200
        )
        print()
        return result.returncode == 0 and dest.exists() and dest.stat().st_size > 0

    try:
        _via_requests()
    except Exception as e:
        if dest.exists():
            dest.unlink()
        print(f"  requests failed ({type(e).__name__}), retrying with wget ...")
        if not _via_wget():
            raise

    return dest


# ── Tier 1: direct h5ad ───────────────────────────────────────────────────────

def _try_tier1(files: list[dict], out_dir: Path, gse: str) -> Optional[Path]:
    h5ad_files = [f for f in files if re.search(r"\.(h5ad|h5)$", f["name"], re.I)]
    if not h5ad_files:
        return None
    # Prefer .h5ad over .h5
    h5ad_files.sort(key=lambda f: (0 if f["name"].endswith(".h5ad") else 1))
    f = h5ad_files[0]
    dest = out_dir / f"{gse}.h5ad"
    _download_file(f["url"], dest)
    return dest


# ── Tier 2a: 10x matrix.mtx.gz ───────────────────────────────────────────────

def _try_tier2_10x(files: list[dict], out_dir: Path, gse: str) -> Optional[Path]:
    def find(pattern):
        matches = [f for f in files if re.search(pattern, f["name"], re.I)]
        return matches[0] if matches else None

    mtx  = find(r"matrix\.mtx\.gz$")
    bars = find(r"barcodes\.tsv\.gz$")
    feat = find(r"(features|genes)\.tsv\.gz$")

    if not (mtx and bars and feat):
        return None

    tmp = out_dir / "_tmp_10x"
    tmp.mkdir(exist_ok=True)

    mtx_path  = _download_file(mtx["url"],  tmp / "matrix.mtx.gz")
    bars_path = _download_file(bars["url"], tmp / "barcodes.tsv.gz")
    feat_path = _download_file(feat["url"], tmp / "features.tsv.gz")

    # Read
    with gzip.open(mtx_path)  as fh: matrix = scipy.io.mmread(fh).T.tocsr()
    with gzip.open(bars_path) as fh: barcodes = pd.read_csv(fh, header=None, sep="\t")[0].tolist()
    with gzip.open(feat_path) as fh:
        feat_df = pd.read_csv(fh, header=None, sep="\t")
        gene_ids   = feat_df[0].tolist()
        gene_names = feat_df[1].tolist() if feat_df.shape[1] > 1 else gene_ids

    adata = ad.AnnData(X=matrix)
    adata.obs_names = barcodes
    adata.var_names = gene_names
    adata.var["gene_id"] = gene_ids

    dest = out_dir / f"{gse}.h5ad"
    adata.write_h5ad(dest)
    print(f"  Saved {adata.shape[0]} cells × {adata.shape[1]} genes → {dest}")

    # Clean up tmp
    import shutil
    shutil.rmtree(tmp)
    return dest


# ── Tier 2b: tab-separated DGE matrix ────────────────────────────────────────

def _try_tier2_dge(files: list[dict], out_dir: Path, gse: str) -> Optional[Path]:
    dge_files = [f for f in files if re.search(r"(count|dge|expression|rpkm|fpkm|tpm).*\.(txt|tsv)(\.gz)?$", f["name"], re.I)]
    if not dge_files:
        return None

    f = dge_files[0]
    tmp_path = out_dir / "_tmp_dge.txt.gz"
    _download_file(f["url"], tmp_path)

    opener = gzip.open if f["name"].endswith(".gz") else open
    with opener(tmp_path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)

    # Genes × cells → transpose to cells × genes
    if df.shape[0] > df.shape[1]:
        df = df.T

    adata = ad.AnnData(X=scipy.sparse.csr_matrix(df.values))
    adata.obs_names = df.index.tolist()
    adata.var_names = df.columns.tolist()

    dest = out_dir / f"{gse}.h5ad"
    adata.write_h5ad(dest)
    print(f"  Saved {adata.shape[0]} cells × {adata.shape[1]} genes → {dest}")
    tmp_path.unlink(missing_ok=True)
    return dest


# ── Tier 2c: RAW.tar containing 10x sample directories ───────────────────────

def _try_tier2c_raw_tar(files: list[dict], out_dir: Path, gse: str) -> Optional[Path]:
    """
    Many GEO datasets package all samples in GSE*_RAW.tar.
    Each sample inside is typically a directory with matrix.mtx.gz / barcodes / features.
    Download the tar, extract the 10x files, concatenate samples → h5ad.
    """
    import tarfile, shutil
    tar_files = [f for f in files if re.search(r"_RAW\.tar$", f["name"], re.I)]
    if not tar_files:
        return None

    tar_f = tar_files[0]
    tmp = out_dir / f"_tmp_raw_{gse}"
    tmp.mkdir(exist_ok=True)
    tar_path = tmp / tar_f["name"]

    _download_file(tar_f["url"], tar_path)
    print(f"  Extracting {tar_f['name']} ...")

    with tarfile.open(tar_path) as tf:
        members = tf.getnames()
        print(f"  Tar contains {len(members)} entries")

        # Check if entries contain 10x-style files
        has_10x = any(re.search(r"matrix\.mtx\.gz|barcodes\.tsv\.gz", m, re.I) for m in members)
        if not has_10x:
            shutil.rmtree(tmp, ignore_errors=True)
            return None

        tf.extractall(tmp)

    # Find all matrix.mtx.gz files (one per sample)
    sample_adatas = []
    for mtx_path in sorted(tmp.rglob("matrix.mtx.gz")):
        sample_dir = mtx_path.parent
        bars_path = next(sample_dir.glob("barcodes.tsv.gz"), None)
        feat_path = next(sample_dir.glob("features.tsv.gz"), None) or \
                    next(sample_dir.glob("genes.tsv.gz"), None)
        if not (bars_path and feat_path):
            continue

        sample_name = sample_dir.name
        with gzip.open(mtx_path)  as fh: matrix = scipy.io.mmread(fh).T.tocsr()
        with gzip.open(bars_path) as fh: barcodes = pd.read_csv(fh, header=None, sep="\t")[0].tolist()
        with gzip.open(feat_path) as fh:
            feat_df = pd.read_csv(fh, header=None, sep="\t")
            gene_ids   = feat_df[0].tolist()
            gene_names = feat_df[1].tolist() if feat_df.shape[1] > 1 else gene_ids

        adata = ad.AnnData(X=matrix)
        adata.obs_names = [f"{sample_name}_{b}" for b in barcodes]
        adata.var_names = gene_names
        adata.var["gene_id"] = gene_ids
        adata.obs["sample"] = sample_name
        sample_adatas.append(adata)
        print(f"    Sample {sample_name}: {adata.shape}")

    if not sample_adatas:
        shutil.rmtree(tmp, ignore_errors=True)
        return None

    print(f"  Concatenating {len(sample_adatas)} samples ...")
    combined = ad.concat(sample_adatas, merge="same")
    dest = out_dir / f"{gse}.h5ad"
    combined.write_h5ad(dest)
    print(f"  Saved {combined.shape[0]} cells × {combined.shape[1]} genes → {dest}")
    shutil.rmtree(tmp, ignore_errors=True)
    return dest


# ── public entry point ────────────────────────────────────────────────────────

def download(accession: str, out_dir: str | Path, retries: int = 3) -> Path:
    """
    Download a GEO dataset to out_dir/{accession}.h5ad.
    Tries Tier 1 → Tier 2a → Tier 2b, raises DataUnavailableError if all fail.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dest = out_dir / f"{accession}.h5ad"
    if dest.exists():
        print(f"  Already downloaded: {dest}")
        return dest

    gse = accession.strip()
    print(f"\n[GEO] {gse}")

    for attempt in range(1, retries + 1):
        try:
            files = _list_suppl_files(gse)
            break
        except Exception as e:
            if attempt == retries:
                raise DataUnavailableError(f"Failed to list files for {gse}: {e}")
            time.sleep(2 ** attempt)

    if not files:
        raise DataUnavailableError(f"{gse}: no supplementary files found")

    print(f"  Found {len(files)} supplementary file(s):")
    for f in files:
        print(f"    {f['name']}")

    result = _try_tier1(files, out_dir, gse)
    if result:
        print(f"  [Tier 1] h5ad downloaded: {result}")
        return result

    result = _try_tier2_10x(files, out_dir, gse)
    if result:
        print(f"  [Tier 2a] 10x matrix reconstructed: {result}")
        return result

    result = _try_tier2c_raw_tar(files, out_dir, gse)
    if result:
        print(f"  [Tier 2c] RAW.tar 10x reconstructed: {result}")
        return result

    result = _try_tier2_dge(files, out_dir, gse)
    if result:
        print(f"  [Tier 2b] DGE matrix reconstructed: {result}")
        return result

    raise DataUnavailableError(
        f"{gse}: no h5ad, 10x matrix, or DGE matrix found in supplementary files. "
        f"May be FASTQ-only — download manually from SRA."
    )
