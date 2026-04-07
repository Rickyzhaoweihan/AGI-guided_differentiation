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
