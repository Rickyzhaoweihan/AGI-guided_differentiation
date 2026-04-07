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
