#!/usr/bin/env python3
"""
geo_search.py — Search GEO (NCBI GDS) and ArrayExpress (EBI BioStudies)
for scRNA-seq datasets, returning linked PMIDs and accession IDs.

Usage:
    python geo_search.py --query "human islet" --output geo_results.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import re
from pathlib import Path
from urllib.parse import urlencode, quote_plus
from urllib.request import urlopen, Request


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/api/v1"
EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
RATE_DELAY = 0.35


# ── GEO (GDS) search ─────────────────────────────────────────────────────────

def search_geo(query: str, max_results: int = 200) -> list[dict]:
    """
    Search NCBI GDS database for scRNA-seq datasets matching query.

    Returns list of dicts with keys:
        GEO_accession, Title, PMID, source
    """
    # Build query: user terms + single-cell filter
    sc_filter = '("single cell" OR "scRNA" OR "snRNA" OR "10x Genomics" OR "single nucleus")'
    organism_filter = '"Homo sapiens"[Organism]'
    full_query = f'{query} AND {sc_filter} AND {organism_filter}'

    print(f"[GEO] Searching GDS: {full_query[:80]}...")
    params = urlencode({
        "db": "gds",
        "term": full_query,
        "retmax": max_results,
        "retmode": "json",
    })
    url = f"{EUTILS}/esearch.fcgi?{params}"
    time.sleep(RATE_DELAY)
    with urlopen(url, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))

    ids = data.get("esearchresult", {}).get("idlist", [])
    total = data.get("esearchresult", {}).get("count", "0")
    print(f"  Found {total} datasets, fetching top {len(ids)}")

    if not ids:
        return []

    # Fetch summaries in batches of 100
    records = []
    for i in range(0, len(ids), 100):
        batch = ids[i:i + 100]
        time.sleep(RATE_DELAY)
        params2 = urlencode({
            "db": "gds",
            "id": ",".join(batch),
            "retmode": "json",
        })
        url2 = f"{EUTILS}/esummary.fcgi?{params2}"
        with urlopen(url2, timeout=30) as resp2:
            summary = json.loads(resp2.read().decode("utf-8"))

        for uid in summary.get("result", {}).get("uids", []):
            doc = summary["result"][uid]
            accession = doc.get("accession", "")
            # Only keep GSE (series) entries, skip GSM (samples) and GPL (platforms)
            if not accession.startswith("GSE"):
                continue
            pmids = doc.get("pubmedids", [])
            title = doc.get("title", "")
            for pmid in pmids:
                records.append({
                    "PMID": str(pmid),
                    "Title": title,
                    "GEO_accession": accession,
                    "source": "GEO",
                })
            # If no PMID linked, still record the dataset
            if not pmids:
                records.append({
                    "PMID": "",
                    "Title": title,
                    "GEO_accession": accession,
                    "source": "GEO",
                })

    print(f"  Extracted {len(records)} GEO records ({len(set(r['PMID'] for r in records if r['PMID']))} unique PMIDs)")
    return records


# ── ArrayExpress (BioStudies) search ──────────────────────────────────────────

def search_arrayexpress(query: str, max_results: int = 100) -> list[dict]:
    """
    Search EBI BioStudies (ArrayExpress collection) for scRNA-seq experiments.

    Returns list of dicts with keys:
        ArrayExpress_accession, Title, PMID, DOI, source
    """
    search_query = f"{query} single cell RNA-seq"
    print(f"[ArrayExpress] Searching BioStudies: {search_query[:80]}...")

    url = (
        f"{BIOSTUDIES}/search?"
        f"query={quote_plus(search_query)}"
        f"&pageSize={max_results}"
        f"&collection=ArrayExpress"
    )
    req = Request(url, headers={"Accept": "application/json"})
    try:
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        print(f"  BioStudies search error: {e}")
        return []

    hits = data.get("hits", [])
    total = data.get("totalHits", 0)
    print(f"  Found {total} experiments, processing top {len(hits)}")

    records = []
    for hit in hits:
        accession = hit.get("accession", "")
        title = hit.get("title", "")
        # Only keep E-MTAB / E-GEOD accessions
        if not (accession.startswith("E-MTAB") or accession.startswith("E-GEOD")):
            continue
        records.append({
            "PMID": "",
            "Title": title,
            "ArrayExpress_accession": accession,
            "source": "ArrayExpress",
        })

    # Resolve PMIDs for records by fetching study details and extracting DOIs
    print(f"  Resolving PMIDs for {len(records)} experiments...")
    resolved = 0
    for rec in records:
        pmid = _resolve_pmid_for_study(rec["ArrayExpress_accession"])
        if pmid:
            rec["PMID"] = pmid
            resolved += 1
        time.sleep(RATE_DELAY)

    print(f"  Resolved {resolved}/{len(records)} PMIDs")
    return records


def _resolve_pmid_for_study(accession: str) -> str:
    """Fetch study detail from BioStudies, extract DOI, resolve to PMID via EuropePMC."""
    try:
        url = f"{BIOSTUDIES}/studies/{accession}"
        req = Request(url, headers={"Accept": "application/json"})
        with urlopen(req, timeout=15) as resp:
            detail = json.loads(resp.read().decode("utf-8"))
    except Exception:
        return ""

    # Look for Publication subsections containing DOIs
    dois = []
    for section in detail.get("section", {}).get("subsections", []):
        if isinstance(section, dict) and section.get("type") == "Publication":
            attrs = {a["name"]: a.get("value", "") for a in section.get("attributes", [])}
            doi = attrs.get("DOI", "")
            if doi:
                dois.append(doi)

    # Resolve first DOI to PMID via EuropePMC
    for doi in dois:
        try:
            time.sleep(RATE_DELAY)
            epm_url = f"{EUROPEPMC}/search?query=DOI:{doi}&format=json&resultType=lite"
            with urlopen(epm_url, timeout=15) as resp:
                epm_data = json.loads(resp.read().decode("utf-8"))
            results = epm_data.get("resultList", {}).get("result", [])
            if results and results[0].get("pmid"):
                return str(results[0]["pmid"])
        except Exception:
            continue

    return ""


# ── Combined search ───────────────────────────────────────────────────────────

def search_data_repositories(
    query: str,
    geo_max: int = 200,
    ae_max: int = 100,
) -> list[dict]:
    """
    Search both GEO and ArrayExpress, merge results.

    Returns list of dicts. Each has at minimum:
        PMID, Title, source, GEO_accession (or ""), ArrayExpress_accession (or "")
    """
    geo_records = search_geo(query, geo_max)
    ae_records = search_arrayexpress(query, ae_max)

    # Normalize: ensure all records have both accession fields
    for r in geo_records:
        r.setdefault("ArrayExpress_accession", "")
    for r in ae_records:
        r.setdefault("GEO_accession", "")

    merged = geo_records + ae_records
    print(f"\n[Data repos] Total: {len(merged)} records "
          f"(GEO: {len(geo_records)}, ArrayExpress: {len(ae_records)})")
    return merged


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Search GEO + ArrayExpress for scRNA-seq datasets.")
    parser.add_argument("--query", required=True, help="Tissue/topic search query")
    parser.add_argument("--geo-max", type=int, default=200, help="Max GEO results")
    parser.add_argument("--ae-max", type=int, default=100, help="Max ArrayExpress results")
    parser.add_argument("--output", default="geo_results.csv", help="Output CSV path")
    args = parser.parse_args()

    records = search_data_repositories(args.query, args.geo_max, args.ae_max)

    if records:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        cols = ["PMID", "Title", "GEO_accession", "ArrayExpress_accession", "source"]
        with open(output, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(records)
        print(f"Wrote {len(records)} records -> {output}")
    else:
        print("No results found.")


if __name__ == "__main__":
    main()
