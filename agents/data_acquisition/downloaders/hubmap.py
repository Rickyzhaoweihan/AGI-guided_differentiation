"""
hubmap.py — Download datasets from the HuBMAP portal.

Resolves HBM* accession → dataset UUID → finds processed h5ad/h5 files
→ downloads via HuBMAP assets API.

Note: Some HuBMAP datasets require a token. Public datasets work without auth.
"""

from __future__ import annotations
import re
from pathlib import Path

import requests

from .geo import DataUnavailableError, _download_file


HUBMAP_SEARCH = "https://search.api.hubmapconsortium.org/v3/search"
HUBMAP_ASSETS = "https://assets.hubmapconsortium.org"


def _resolve_uuid(accession: str, token: str | None = None) -> str | None:
    """Resolve an HBM* display ID to a dataset UUID via HuBMAP Search API."""
    headers = {"Authorization": f"Bearer {token}"} if token else {}
    payload = {
        "query": {
            "bool": {
                "must": [
                    {"match": {"hubmap_id": accession}},
                    {"match": {"entity_type": "Dataset"}},
                ]
            }
        }
    }
    resp = requests.post(HUBMAP_SEARCH, json=payload, headers=headers, timeout=30)
    resp.raise_for_status()
    hits = resp.json().get("hits", {}).get("hits", [])
    if not hits:
        return None
    return hits[0]["_source"].get("uuid")


def _list_dataset_files(uuid: str, token: str | None = None) -> list[dict]:
    """List files in a HuBMAP dataset via the files API."""
    headers = {"Authorization": f"Bearer {token}"} if token else {}
    url = f"https://entity.api.hubmapconsortium.org/entities/{uuid}/files"
    resp = requests.get(url, headers=headers, timeout=30)
    if not resp.ok:
        return []
    return resp.json()


def download(accession: str, out_dir: str | Path, token: str | None = None) -> Path:
    """
    Download a HuBMAP dataset (RNA expression h5ad preferred).
    Pass token= for protected datasets.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dest = out_dir / f"{accession}.h5ad"
    if dest.exists():
        print(f"  Already downloaded: {dest}")
        return dest

    print(f"\n[HuBMAP] {accession}")

    uuid = _resolve_uuid(accession, token)
    if not uuid:
        raise DataUnavailableError(
            f"{accession}: could not resolve to a dataset UUID. "
            f"Check https://portal.hubmapconsortium.org/browse/collection/"
        )

    print(f"  UUID: {uuid}")
    files = _list_dataset_files(uuid, token)

    # Prefer RNA h5ad over ATAC; prefer processed over raw
    h5ad_files = [
        f for f in files
        if re.search(r"\.(h5ad|h5)$", f.get("rel_path", ""), re.I)
        and not re.search(r"atac|fragment", f.get("rel_path", ""), re.I)
    ]

    if not h5ad_files:
        raise DataUnavailableError(
            f"{accession} (UUID: {uuid}): no RNA h5ad found. "
            f"View dataset at https://portal.hubmapconsortium.org/browse/dataset/{uuid}"
        )

    # Prefer files with 'expr', 'rna', 'gene' in name; else take first
    h5ad_files.sort(key=lambda f: (
        0 if re.search(r"expr|rna|gene", f.get("rel_path", ""), re.I) else 1
    ))

    rel_path = h5ad_files[0]["rel_path"]
    url = f"{HUBMAP_ASSETS}/{uuid}/{rel_path}"
    headers = {"Authorization": f"Bearer {token}"} if token else {}

    print(f"  Downloading {rel_path} ...")
    _download_file(url, dest)
    return dest
