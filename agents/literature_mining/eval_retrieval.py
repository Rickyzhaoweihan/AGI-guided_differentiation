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
