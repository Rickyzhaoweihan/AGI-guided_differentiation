# Literature Retrieval v2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a two-channel literature retrieval pipeline (PubMed + GEO/ArrayExpress) with LLM-based relevance filtering that achieves >40% precision and >80% recall for finding scRNA-seq dataset papers.

**Architecture:** Channel 1 searches PubMed with improved query expansion (atlas, organoid, embryo categories added, disease queries removed). Channel 2 searches GEO (NCBI GDS) and ArrayExpress (EBI BioStudies) directly for scRNA-seq datasets, resolving back to PMIDs. Results merge by PMID, then pass through a Claude API filter that classifies each abstract for relevance (understands developmental lineage, multi-organ atlases). Output is a CSV with relevance scores.

**Tech Stack:** Python 3.12, anthropic SDK 0.89.0, NCBI eutils, EBI BioStudies API, EuropePMC API. Conda env: `/nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env`

**Run all commands from:** `/nfs/turbo/umms-drjieliu/usr/rickyhan/AGI-guided_differentiation`

**Python prefix:** `conda run -p /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env python`

---

## File Map

```
agents/literature_mining/
├── query_expansion.py     ← NEW: refactored query expansion logic
├── geo_search.py          ← NEW: GEO (GDS) + ArrayExpress (BioStudies) search
├── llm_filter.py          ← NEW: Claude API relevance filter
├── mine_literature_v2.py  ← NEW: main orchestrator CLI
├── eval_retrieval.py      ← NEW: evaluation script against ground truth xlsx
├── pubmed_search.py       ← EXISTING: unchanged, reused
├── hirn_search.py         ← EXISTING: unchanged, reused
└── mine_literature.py     ← EXISTING: unchanged (v1, kept for reference)
```

---

### Task 1: query_expansion.py — Improved Query Expansion

**Files:**
- Create: `agents/literature_mining/query_expansion.py`

- [ ] **Step 1: Create `query_expansion.py` with `expand_query` function**

```python
#!/usr/bin/env python3
"""
query_expansion.py — Generate expanded PubMed queries for scRNA-seq dataset retrieval.

Given a user query (e.g. "human gut scRNA-seq"), produces multiple sub-queries
covering technology variants, tissue synonyms, organoid/iPSC, atlas, and
developmental categories. Disease-specific queries are NOT generated (user can
include disease terms in their base query if desired).
"""

from __future__ import annotations

import re


# Map common tissue keywords to synonyms for query expansion
TISSUE_SYNONYMS: dict[str, list[str]] = {
    "gut":       ["intestine", "intestinal", "colon", "duodenum", "ileum", "enterocyte", "gastrointestinal"],
    "intestin":  ["gut", "colon", "duodenum", "ileum", "enterocyte", "gastrointestinal", "bowel"],
    "colon":     ["intestine", "intestinal", "gut", "bowel"],
    "pancrea":   ["islet", "islets", "beta cell", "endocrine pancreas"],
    "islet":     ["pancreas", "pancreatic", "beta cell", "endocrine pancreas"],
    "heart":     ["cardiac", "cardiomyocyte", "myocardium"],
    "liver":     ["hepatocyte", "hepatic"],
    "lung":      ["pulmonary", "airway", "alveolar"],
    "kidney":    ["renal", "nephron"],
    "brain":     ["cerebral", "cortex", "neuron", "neural"],
}

# Regex to strip technology terms from the query to get the "base" tissue terms
_TECH_RE = re.compile(
    r"(?i)\b(scrna-?seq|snrna-?seq|single[- ]cell|single[- ]nucleus|"
    r"RNA[- ]seq(?:uencing)?|transcriptom\w*|scRNA|10x|chromium)\b"
)


def expand_query(query: str) -> list[str]:
    """
    Expand a user query into multiple PubMed sub-queries to maximize recall.

    Categories:
      1. Core query (as-is)
      2. Technology variants (scRNA-seq, snRNA-seq, single-cell transcriptome)
      3. Tissue synonyms (capped at top 3)
      4. Organoid / iPSC / stem cell differentiation
      5. Atlas / multi-organ
      6. Developmental (fetal, embryo, organogenesis)

    Returns deduplicated list of query strings.
    """
    queries: list[str] = [query]
    q_lower = query.lower()

    # Extract base tissue terms (remove technology keywords)
    base = _TECH_RE.sub("", query).strip()
    base = re.sub(r"\s+", " ", base).strip()

    # --- 1. Technology variants ---
    if base and len(base) > 3:
        queries.append(f"{base} single-cell RNA-seq")
        queries.append(f"{base} single-nucleus RNA-seq")
        queries.append(f"{base} single-cell transcriptome")

    # --- 2. Tissue synonyms (top 3) ---
    for key, syns in TISSUE_SYNONYMS.items():
        if key in q_lower:
            for syn in syns[:3]:
                if syn.lower() not in q_lower:
                    queries.append(f"human {syn} single-cell RNA-seq")
            break  # only match first tissue

    # --- 3. Organoid / iPSC ---
    if base and len(base) > 3:
        queries.append(f"{base} organoid single-cell RNA-seq")
        queries.append(f"human pluripotent stem cell {base} differentiation single-cell")
        queries.append(f"iPSC {base} single-cell RNA-seq")

    # --- 4. Atlas / multi-organ ---
    queries.append("human cell atlas single-cell transcriptome")
    queries.append("Tabula Sapiens single-cell")
    queries.append("human fetal gene expression atlas single-cell")
    if base and len(base) > 3:
        queries.append(f"{base} cell atlas single-cell")

    # --- 5. Developmental ---
    if base and len(base) > 3:
        queries.append(f"{base} fetal development single-cell RNA-seq")
        queries.append(f"human embryo organogenesis single-cell {base}")
        queries.append(f"{base} embryo single-cell transcriptome")

    # --- Deduplicate preserving order ---
    seen: set[str] = set()
    unique: list[str] = []
    for q in queries:
        q_norm = re.sub(r"\s+", " ", q.strip().lower())
        if q_norm not in seen:
            seen.add(q_norm)
            unique.append(q.strip())

    return unique
```

- [ ] **Step 2: Smoke-test query expansion**

Run:
```bash
conda run -p /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env python -c "
from agents.literature_mining.query_expansion import expand_query
queries = expand_query('human pancreatic islet scRNA-seq')
for i, q in enumerate(queries):
    print(f'{i+1}. {q}')
print(f'Total: {len(queries)} sub-queries')
"
```

Expected: ~15-18 sub-queries covering technology, synonym, organoid, atlas, and developmental axes. No disease queries.

- [ ] **Step 3: Commit**

```bash
git add agents/literature_mining/query_expansion.py
git commit -m "feat: add improved query expansion for literature retrieval v2"
```

---

### Task 2: geo_search.py — GEO + ArrayExpress Search

**Files:**
- Create: `agents/literature_mining/geo_search.py`

**API details discovered during research:**
- GEO: Use eutils `esearch` on `db=gds`, then `esummary` to get accession + PMIDs
- ArrayExpress: Use BioStudies API `https://www.ebi.ac.uk/biostudies/api/v1/search?query=...&collection=ArrayExpress`
- BioStudies returns DOIs not PMIDs — resolve via EuropePMC: `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:...`

- [ ] **Step 1: Create `geo_search.py`**

```python
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
```

- [ ] **Step 2: Smoke-test GEO search**

Run:
```bash
conda run -p /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env python agents/literature_mining/geo_search.py --query "human islet" --output /tmp/geo_test.csv
```

Expected: Some GDS results with GSE accessions and PMIDs. ArrayExpress results with E-MTAB accessions.

- [ ] **Step 3: Commit**

```bash
git add agents/literature_mining/geo_search.py
git commit -m "feat: add GEO + ArrayExpress search for literature retrieval v2"
```

---

### Task 3: llm_filter.py — Claude API Relevance Filter

**Files:**
- Create: `agents/literature_mining/llm_filter.py`

- [ ] **Step 1: Create `llm_filter.py`**

```python
#!/usr/bin/env python3
"""
llm_filter.py — Filter paper records for relevance using Claude API.

Sends batches of paper titles + abstracts to Claude, which classifies each
as relevant or not for the target tissue/topic. Understands multi-organ atlases,
developmental lineage, and organoid models.

Usage:
    # Typically called from mine_literature_v2.py, not directly.
    python llm_filter.py --input raw_results.csv --query "human gut scRNA-seq" \
        --api-key "sk-ant-..." --output filtered.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import anthropic


BATCH_SIZE = 10  # papers per API call
DEFAULT_MODEL = "claude-sonnet-4-20250514"

SYSTEM_PROMPT = "You are a bioinformatics expert specializing in single-cell genomics."

USER_PROMPT_TEMPLATE = """For each paper below, determine whether it reports or re-analyzes single-cell or single-nucleus RNA-seq data that includes {tissue_description}.

Consider:
- Multi-organ atlases that include this tissue (e.g., Tabula Sapiens, Human Cell Landscape, fetal gene expression atlas)
- Developmental biology: embryo/fetal datasets where this tissue lineage originates (e.g., endoderm for gut, mesoderm for heart)
- Organoid or iPSC/hPSC-derived models of this tissue
- Papers that generated original scRNA-seq/snRNA-seq data vs reviews/methods papers that don't contain data

Respond with a JSON array. For each paper, include:
- "pmid": the PMID string
- "relevant": "yes" or "no"
- "confidence": "high", "medium", or "low"
- "reasoning": one sentence explanation

Papers:
{papers_json}

Respond ONLY with the JSON array, no other text."""


def filter_papers(
    records: list[dict],
    tissue_description: str,
    api_key: str,
    model: str = DEFAULT_MODEL,
    verbose: bool = False,
) -> list[dict]:
    """
    Filter paper records using Claude API.

    Args:
        records: List of dicts with at least PMID, Title, Abstract keys.
        tissue_description: Description of the target tissue (e.g., "human pancreatic islet tissue").
        api_key: Anthropic API key.
        model: Claude model to use.
        verbose: Print progress details.

    Returns:
        The same records list, with three new keys added to each:
            LLM_relevant (yes/no), LLM_confidence (high/medium/low), LLM_reasoning (str)
    """
    client = anthropic.Anthropic(api_key=api_key)

    # Index records by PMID for fast lookup
    by_pmid: dict[str, dict] = {}
    for r in records:
        pmid = r.get("PMID", "")
        if pmid:
            by_pmid[pmid] = r

    # Process in batches
    papers_with_abstracts = [r for r in records if r.get("PMID") and r.get("Abstract")]
    papers_without_abstracts = [r for r in records if not r.get("Abstract")]

    total_batches = (len(papers_with_abstracts) + BATCH_SIZE - 1) // BATCH_SIZE
    print(f"[LLM Filter] {len(papers_with_abstracts)} papers with abstracts, "
          f"{len(papers_without_abstracts)} without abstracts")
    print(f"  Processing {total_batches} batches of {BATCH_SIZE}...")

    for batch_idx in range(0, len(papers_with_abstracts), BATCH_SIZE):
        batch = papers_with_abstracts[batch_idx:batch_idx + BATCH_SIZE]
        batch_num = batch_idx // BATCH_SIZE + 1

        # Build papers JSON for prompt
        papers_json = json.dumps([
            {
                "pmid": r["PMID"],
                "title": r.get("Title", ""),
                "abstract": r.get("Abstract", "")[:1500],  # cap abstract length
            }
            for r in batch
        ], indent=2)

        prompt = USER_PROMPT_TEMPLATE.format(
            tissue_description=tissue_description,
            papers_json=papers_json,
        )

        if verbose:
            print(f"  Batch {batch_num}/{total_batches}: {len(batch)} papers")

        try:
            response = client.messages.create(
                model=model,
                max_tokens=2048,
                system=SYSTEM_PROMPT,
                messages=[{"role": "user", "content": prompt}],
            )
            response_text = response.content[0].text.strip()

            # Parse JSON response — handle markdown code fences
            if response_text.startswith("```"):
                response_text = response_text.split("\n", 1)[1]
                response_text = response_text.rsplit("```", 1)[0]
            response_text = response_text.strip()

            results = json.loads(response_text)

            for item in results:
                pmid = str(item.get("pmid", ""))
                if pmid in by_pmid:
                    by_pmid[pmid]["LLM_relevant"] = item.get("relevant", "no")
                    by_pmid[pmid]["LLM_confidence"] = item.get("confidence", "low")
                    by_pmid[pmid]["LLM_reasoning"] = item.get("reasoning", "")

        except json.JSONDecodeError as e:
            print(f"  WARNING: Batch {batch_num} JSON parse error: {e}")
            if verbose:
                print(f"  Response was: {response_text[:200]}")
            # Mark batch as uncertain
            for r in batch:
                r.setdefault("LLM_relevant", "unknown")
                r.setdefault("LLM_confidence", "low")
                r.setdefault("LLM_reasoning", "JSON parse error")

        except anthropic.APIError as e:
            print(f"  WARNING: Batch {batch_num} API error: {e}")
            for r in batch:
                r.setdefault("LLM_relevant", "unknown")
                r.setdefault("LLM_confidence", "low")
                r.setdefault("LLM_reasoning", f"API error: {e}")

    # Papers without abstracts: mark as unknown (can't classify)
    for r in papers_without_abstracts:
        r.setdefault("LLM_relevant", "unknown")
        r.setdefault("LLM_confidence", "low")
        r.setdefault("LLM_reasoning", "no abstract available")

    # Ensure all records have the LLM fields
    for r in records:
        r.setdefault("LLM_relevant", "unknown")
        r.setdefault("LLM_confidence", "low")
        r.setdefault("LLM_reasoning", "not processed")

    # Summary
    yes_count = sum(1 for r in records if r.get("LLM_relevant") == "yes")
    no_count = sum(1 for r in records if r.get("LLM_relevant") == "no")
    unk_count = sum(1 for r in records if r.get("LLM_relevant") == "unknown")
    print(f"  Results: {yes_count} relevant, {no_count} not relevant, {unk_count} unknown")

    return records


# ── CLI (for standalone testing) ──────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Filter papers for relevance using Claude API.")
    parser.add_argument("--input", required=True, help="Input CSV with PMID, Title, Abstract columns")
    parser.add_argument("--query", required=True, help="Tissue description for relevance check")
    parser.add_argument("--api-key", required=True, help="Anthropic API key")
    parser.add_argument("--model", default=DEFAULT_MODEL, help=f"Claude model (default: {DEFAULT_MODEL})")
    parser.add_argument("--output", default="filtered.csv", help="Output CSV path")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    # Read input CSV
    with open(args.input, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        records = list(reader)

    print(f"Loaded {len(records)} records from {args.input}")
    records = filter_papers(records, args.query, args.api_key, args.model, args.verbose)

    # Write output
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    cols = list(records[0].keys()) if records else []
    with open(output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)
    print(f"Wrote {len(records)} records -> {output}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Test LLM filter with a small batch (3 papers)**

Run:
```bash
conda run -p /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env python -c "
from agents.literature_mining.llm_filter import filter_papers
import sys

# 3 test papers: 1 clearly relevant, 1 clearly not, 1 edge case (atlas)
test_records = [
    {'PMID': '27667667', 'Title': 'Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes', 'Abstract': 'Pancreatic islet cells are critical for glucose homeostasis. We performed single-cell RNA sequencing of human islets.'},
    {'PMID': '99999999', 'Title': 'Review of Machine Learning in Genomics', 'Abstract': 'This review covers applications of machine learning to genomic data analysis.'},
    {'PMID': '35549404', 'Title': 'The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans', 'Abstract': 'We created a multiple-organ single-cell transcriptomic atlas across 24 tissues.'},
]

# Replace YOUR_KEY with actual API key for testing
api_key = sys.argv[1] if len(sys.argv) > 1 else 'test'
if api_key == 'test':
    print('Skipping live test - pass API key as argument')
else:
    results = filter_papers(test_records, 'human pancreatic islet tissue', api_key, verbose=True)
    for r in results:
        print(f'  PMID {r[\"PMID\"]}: relevant={r[\"LLM_relevant\"]} confidence={r[\"LLM_confidence\"]} reason={r[\"LLM_reasoning\"]}')
"
```

Expected: PMID 27667667 → yes/high, 99999999 → no/high, 35549404 → yes/medium (atlas includes pancreas).

- [ ] **Step 3: Commit**

```bash
git add agents/literature_mining/llm_filter.py
git commit -m "feat: add Claude API relevance filter for literature retrieval v2"
```

---

### Task 4: mine_literature_v2.py — Main Orchestrator

**Files:**
- Create: `agents/literature_mining/mine_literature_v2.py`

- [ ] **Step 1: Create `mine_literature_v2.py`**

```python
#!/usr/bin/env python3
"""
mine_literature_v2.py — Two-channel literature retrieval with LLM filtering.

Channel 1: PubMed broad search (improved query expansion)
Channel 2: GEO + ArrayExpress dataset search
Filter:    Claude API relevance classification

Usage:
    python mine_literature_v2.py \
        --query "human pancreatic islet scRNA-seq" \
        --api-key "sk-ant-..." \
        --output results.csv

    python mine_literature_v2.py \
        --query "human gut intestinal development scRNA-seq" \
        --api-key "sk-ant-..." \
        --output results.csv \
        --verbose
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from pubmed_search import search_pubmed, efetch
from query_expansion import expand_query
from geo_search import search_data_repositories
from llm_filter import filter_papers

# Optional HIRN import
try:
    from hirn_search import search_hirn
    HIRN_AVAILABLE = True
except ImportError:
    HIRN_AVAILABLE = False

ISLET_KEYWORDS = re.compile(
    r"islet|pancrea|beta.cell|alpha.cell|delta.cell|insulin|glucagon|langerhans|endocrine.pancrea",
    re.I,
)

OUTPUT_COLUMNS = [
    "PMID", "Title", "Abstract", "Year", "Journal", "Authors", "DOI",
    "source", "GEO_accession", "ArrayExpress_accession",
    "LLM_relevant", "LLM_confidence", "LLM_reasoning",
]


def mine_v2(
    query: str,
    api_key: str,
    max_results: int = 50,
    model: str = "claude-sonnet-4-20250514",
    use_hirn: bool | str = "auto",
    no_geo: bool = False,
    no_filter: bool = False,
    verbose: bool = False,
) -> list[dict]:
    """
    Two-channel retrieval with LLM filtering.

    Args:
        query: User search query (e.g., "human pancreatic islet scRNA-seq").
        api_key: Anthropic API key for LLM filter.
        max_results: Max PubMed results per sub-query.
        model: Claude model for LLM filter.
        use_hirn: "auto", True, or False for HIRN search.
        no_geo: Skip GEO/ArrayExpress channel.
        no_filter: Skip LLM filter (return raw merged results).
        verbose: Print detailed progress.

    Returns:
        List of paper records (dicts).
    """

    # ── Channel 1: PubMed broad retrieval ──
    print("\n" + "=" * 60)
    print("CHANNEL 1: PubMed Retrieval")
    print("=" * 60)

    queries = expand_query(query)
    print(f"\n[Query expansion] {len(queries)} sub-queries:")
    for i, q in enumerate(queries):
        print(f"  {i + 1}. {q}")
    print()

    seen_pmids: set[str] = set()
    records: list[dict] = []

    for q in queries:
        pm_records = search_pubmed(q, max_results)
        new = 0
        for r in pm_records:
            pmid = str(r.get("PMID", ""))
            if pmid and pmid not in seen_pmids:
                r["source"] = "PubMed"
                r.setdefault("GEO_accession", "")
                r.setdefault("ArrayExpress_accession", "")
                records.append(r)
                seen_pmids.add(pmid)
                new += 1
        if verbose:
            print(f"  -> {new} new unique papers")

    print(f"\n[PubMed] Total unique papers: {len(records)}")

    # ── HIRN (if relevant) ──
    run_hirn = (
        HIRN_AVAILABLE
        and (use_hirn is True or (use_hirn == "auto" and ISLET_KEYWORDS.search(query)))
    )
    if run_hirn:
        print("\n[HIRN] Searching...")
        hirn_records = search_hirn(queries[0], max_results=30)
        hirn_new = 0
        for r in hirn_records:
            pmid = str(r.get("PMID", ""))
            if pmid and pmid not in seen_pmids:
                r.setdefault("GEO_accession", "")
                r.setdefault("ArrayExpress_accession", "")
                records.append(r)
                seen_pmids.add(pmid)
                hirn_new += 1
        print(f"  HIRN added {hirn_new} new papers")

    # ── Channel 2: GEO + ArrayExpress ──
    if not no_geo:
        print("\n" + "=" * 60)
        print("CHANNEL 2: GEO + ArrayExpress")
        print("=" * 60)

        # Extract base tissue terms for data repository search
        tech_re = re.compile(
            r"(?i)\b(scrna-?seq|snrna-?seq|single[- ]cell|single[- ]nucleus|"
            r"RNA[- ]seq(?:uencing)?|transcriptom\w*|scRNA|10x|chromium)\b"
        )
        base_query = tech_re.sub("", query).strip()
        base_query = re.sub(r"\s+", " ", base_query).strip()
        if not base_query:
            base_query = query

        geo_records = search_data_repositories(base_query)

        # Merge: for records with PMIDs we already have, add accessions
        # For new PMIDs, fetch metadata from PubMed and add
        new_pmids: list[str] = []
        for gr in geo_records:
            pmid = gr.get("PMID", "")
            geo_acc = gr.get("GEO_accession", "")
            ae_acc = gr.get("ArrayExpress_accession", "")

            if pmid and pmid in seen_pmids:
                # Update existing record with accession info
                for r in records:
                    if r["PMID"] == pmid:
                        if geo_acc:
                            existing = r.get("GEO_accession", "")
                            if geo_acc not in existing:
                                r["GEO_accession"] = f"{existing},{geo_acc}".strip(",")
                        if ae_acc:
                            existing = r.get("ArrayExpress_accession", "")
                            if ae_acc not in existing:
                                r["ArrayExpress_accession"] = f"{existing},{ae_acc}".strip(",")
                        if "GEO" not in r.get("source", "") and geo_acc:
                            r["source"] += ",GEO"
                        if "ArrayExpress" not in r.get("source", "") and ae_acc:
                            r["source"] += ",ArrayExpress"
                        break
            elif pmid and pmid not in seen_pmids:
                new_pmids.append(pmid)
                seen_pmids.add(pmid)

        # Fetch PubMed metadata for new PMIDs from data repos
        if new_pmids:
            print(f"\n[PubMed] Fetching metadata for {len(new_pmids)} new PMIDs from data repos...")
            pm_meta = efetch(list(set(new_pmids)))
            pm_by_pmid = {r["PMID"]: r for r in pm_meta}

            for gr in geo_records:
                pmid = gr.get("PMID", "")
                if pmid in pm_by_pmid and not any(r["PMID"] == pmid for r in records):
                    r = pm_by_pmid[pmid]
                    r["source"] = gr.get("source", "GEO")
                    r["GEO_accession"] = gr.get("GEO_accession", "")
                    r["ArrayExpress_accession"] = gr.get("ArrayExpress_accession", "")
                    records.append(r)

    print(f"\n[Merged] {len(records)} unique papers total")

    # ── LLM Filter ──
    if not no_filter:
        print("\n" + "=" * 60)
        print("LLM RELEVANCE FILTER")
        print("=" * 60)
        records = filter_papers(records, query, api_key, model, verbose)

    return records


def write_csv(records: list[dict], output: str | Path) -> Path:
    """Write records to CSV."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)
    print(f"\nWrote {len(records)} records -> {output}")
    return output


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Two-channel literature retrieval with LLM filtering.",
    )
    parser.add_argument("--query", required=True, help="Search query (e.g., 'human pancreatic islet scRNA-seq')")
    parser.add_argument("--api-key", required=True, help="Anthropic API key")
    parser.add_argument("--output", default="literature_results_v2.csv", help="Output CSV path")
    parser.add_argument("--max", type=int, default=50, help="Max PubMed results per sub-query (default: 50)")
    parser.add_argument("--model", default="claude-sonnet-4-20250514", help="Claude model for LLM filter")
    parser.add_argument("--no-geo", action="store_true", help="Skip GEO/ArrayExpress channel")
    parser.add_argument("--no-filter", action="store_true", help="Skip LLM filter (raw merged results)")
    parser.add_argument("--use-hirn", action="store_true", help="Force HIRN search")
    parser.add_argument("--verbose", action="store_true", help="Print detailed progress")
    args = parser.parse_args()

    use_hirn = True if args.use_hirn else "auto"
    records = mine_v2(
        query=args.query,
        api_key=args.api_key,
        max_results=args.max,
        model=args.model,
        use_hirn=use_hirn,
        no_geo=args.no_geo,
        no_filter=args.no_filter,
        verbose=args.verbose,
    )

    if records:
        write_csv(records, args.output)

        # Print summary
        if not args.no_filter:
            relevant = [r for r in records if r.get("LLM_relevant") == "yes"]
            print(f"\n{'=' * 60}")
            print(f"SUMMARY: {len(relevant)} relevant papers out of {len(records)} total")
            print(f"{'=' * 60}")
            for r in relevant:
                conf = r.get("LLM_confidence", "?")
                print(f"  [{conf}] PMID {r['PMID']}: {r.get('Title', '?')[:80]}")
    else:
        print("No results found.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run a quick end-to-end test with `--no-filter` (no API key needed)**

Run:
```bash
cd agents/literature_mining && conda run -p /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env python mine_literature_v2.py \
    --query "human pancreatic islet scRNA-seq" \
    --api-key "dummy" \
    --no-filter \
    --output /tmp/test_v2_raw.csv \
    --verbose
```

Expected: Should retrieve papers from PubMed + GEO + ArrayExpress channels and write merged CSV.

- [ ] **Step 3: Commit**

```bash
git add agents/literature_mining/mine_literature_v2.py
git commit -m "feat: add mine_literature_v2 orchestrator with two-channel retrieval"
```

---

### Task 5: eval_retrieval.py — Evaluation Against Ground Truth

**Files:**
- Create: `agents/literature_mining/eval_retrieval.py`

- [ ] **Step 1: Create `eval_retrieval.py`**

```python
#!/usr/bin/env python3
"""
eval_retrieval.py — Evaluate retrieval pipeline against ground truth xlsx tables.

Compares the PMID set retrieved by mine_literature_v2 against curated tables
to compute precision, recall, and F1.

Usage:
    python eval_retrieval.py \
        --results results.csv \
        --ground-truth ../../Databasehumanislets.xlsx \
        --filter-relevant
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd


def load_ground_truth(xlsx_path: str) -> set[str]:
    """Load PMIDs from ground truth xlsx. Searches all sheets for a 'PMID' column."""
    xl = pd.ExcelFile(xlsx_path)
    pmids: set[str] = set()
    for sheet in xl.sheet_names:
        df = pd.read_excel(xlsx_path, sheet_name=sheet)
        if "PMID" in df.columns:
            for val in df["PMID"].dropna():
                pmids.add(str(int(float(val))))
    return pmids


def load_results(csv_path: str, filter_relevant: bool = False) -> set[str]:
    """Load PMIDs from results CSV. Optionally filter to LLM_relevant=yes only."""
    pmids: set[str] = set()
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if filter_relevant and row.get("LLM_relevant") != "yes":
                continue
            pmid = row.get("PMID", "").strip()
            if pmid:
                pmids.add(pmid)
    return pmids


def evaluate(retrieved: set[str], ground_truth: set[str]) -> dict:
    """Compute precision, recall, F1."""
    tp = retrieved & ground_truth
    fp = retrieved - ground_truth
    fn = ground_truth - retrieved

    precision = len(tp) / len(retrieved) if retrieved else 0.0
    recall = len(tp) / len(ground_truth) if ground_truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    return {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "tp": sorted(tp),
        "fp_count": len(fp),
        "fn": sorted(fn),
        "retrieved_count": len(retrieved),
        "gt_count": len(ground_truth),
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate retrieval against ground truth.")
    parser.add_argument("--results", required=True, help="Results CSV from mine_literature_v2")
    parser.add_argument("--ground-truth", required=True, help="Ground truth xlsx")
    parser.add_argument("--filter-relevant", action="store_true",
                        help="Only count LLM_relevant=yes papers as retrieved")
    args = parser.parse_args()

    gt = load_ground_truth(args.ground_truth)
    retrieved = load_results(args.results, args.filter_relevant)

    print(f"Ground truth: {len(gt)} PMIDs from {args.ground_truth}")
    print(f"Retrieved:    {len(retrieved)} PMIDs from {args.results}"
          + (" (filtered to LLM_relevant=yes)" if args.filter_relevant else ""))

    metrics = evaluate(retrieved, gt)

    print(f"\nPrecision: {metrics['precision']:.1%} ({len(metrics['tp'])}/{metrics['retrieved_count']})")
    print(f"Recall:    {metrics['recall']:.1%} ({len(metrics['tp'])}/{metrics['gt_count']})")
    print(f"F1:        {metrics['f1']:.1%}")

    print(f"\nTrue positives ({len(metrics['tp'])}):")
    for pmid in metrics["tp"]:
        print(f"  {pmid}")

    print(f"\nFalse negatives ({len(metrics['fn'])}) — missed papers:")
    for pmid in metrics["fn"]:
        print(f"  {pmid}")

    print(f"\nFalse positives: {metrics['fp_count']} papers not in ground truth")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add agents/literature_mining/eval_retrieval.py
git commit -m "feat: add evaluation script for literature retrieval v2"
```

---

### Task 6: End-to-End Integration Test + Evaluation

This task requires an actual API key and makes live API calls.

- [ ] **Step 1: Run full pipeline on islet query**

```bash
cd agents/literature_mining && conda run -p /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env python mine_literature_v2.py \
    --query "human pancreatic islet scRNA-seq" \
    --api-key "$ANTHROPIC_API_KEY" \
    --output /tmp/islet_v2.csv \
    --verbose
```

- [ ] **Step 2: Evaluate islet results**

```bash
cd agents/literature_mining && conda run -p /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env python eval_retrieval.py \
    --results /tmp/islet_v2.csv \
    --ground-truth ../../Databasehumanislets.xlsx \
    --filter-relevant
```

Expected: Precision >40%, Recall >80%.

- [ ] **Step 3: Run full pipeline on gut query**

```bash
cd agents/literature_mining && conda run -p /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env python mine_literature_v2.py \
    --query "human gut intestinal development scRNA-seq" \
    --api-key "$ANTHROPIC_API_KEY" \
    --output /tmp/gut_v2.csv \
    --verbose
```

- [ ] **Step 4: Evaluate gut results**

```bash
cd agents/literature_mining && conda run -p /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env python eval_retrieval.py \
    --results /tmp/gut_v2.csv \
    --ground-truth "../../ML for gut development.xlsx" \
    --filter-relevant
```

Expected: Precision >40%, Recall >60% (gut is harder due to embryo/stem cell papers).

- [ ] **Step 5: If targets not met, iterate on query expansion or LLM prompt**

Review false negatives — if they're papers with no tissue terms AND no scRNA-seq terms in abstract, they may need the full-text retrieval step (out of scope for v2). Document remaining gaps.

- [ ] **Step 6: Final commit**

```bash
git add -A agents/literature_mining/
git commit -m "feat: literature retrieval v2 — two-channel + LLM filter pipeline"
```
