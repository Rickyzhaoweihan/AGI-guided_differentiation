#!/usr/bin/env python3
"""
acquire_datasets.py — Download scRNA-seq datasets from a CSV (or xlsx) catalog.

Usage:
    python acquire_datasets.py --catalog <file.csv> --output-dir data/downloaded/
    python acquire_datasets.py --catalog <file.xlsx> --output-dir data/downloaded/ --dry-run
    python acquire_datasets.py --catalog <file.csv> --accession GSE158702   # single dataset

The catalog file must have columns for:
    - title / manuscript name
    - accession / data resource  (GSE*, E-MTAB-*, HBM*, ...)
    - url / link to data

Supported database types and their download method:
    GEO      (GSE*)        → NCBI FTP supplementary files (Tier 1 h5ad, Tier 2 10x matrix)
    EBI      (E-MTAB-*)    → BioStudies API (Tier 1 h5ad, Tier 3 TSV matrix)
    HuBMAP   (HBM*)        → HuBMAP Search + Assets API
    EGA      (EGAD*/EGAS*) → SKIP — controlled access, requires DAC application
    CNCB     (HRA*)        → SKIP — controlled access, requires CNCB registration
    dbGaP    (phs*)        → SKIP — controlled access, requires dbGaP approval
    PANC DB               → SKIP — requires account at hpap.pmacs.upenn.edu
    DUOS                  → SKIP — requires account
"""

from __future__ import annotations
import argparse
import sys
import traceback
from pathlib import Path

from catalog import load_catalog, print_summary, DatasetRecord
from downloaders.geo import download as geo_download, DataUnavailableError
from downloaders.ebi import download as ebi_download
from downloaders.hubmap import download as hubmap_download


# ── result tracking ───────────────────────────────────────────────────────────

class Result:
    def __init__(self, record: DatasetRecord):
        self.record = record
        self.status = "pending"   # pending | downloaded | skipped | failed
        self.output_path: Path | None = None
        self.message: str = ""


# ── dispatcher ────────────────────────────────────────────────────────────────

def acquire_one(record: DatasetRecord, out_dir: Path, hubmap_token: str | None = None) -> Result:
    res = Result(record)

    if not record.accessible:
        res.status = "skipped"
        res.message = f"Controlled access ({record.db_type}) — requires manual application"
        return res

    try:
        if record.db_type == "GEO":
            path = geo_download(record.accession, out_dir)
        elif record.db_type == "EBI":
            path = ebi_download(record.accession, out_dir)
        elif record.db_type == "HUBMAP":
            path = hubmap_download(record.accession, out_dir, token=hubmap_token)
        elif record.db_type == "S3":
            res.status = "skipped"
            res.message = "AWS S3 — install boto3 and download manually: " + (record.url or "")
            return res
        else:
            res.status = "skipped"
            res.message = f"No downloader implemented for {record.db_type}"
            return res

        res.status = "downloaded"
        res.output_path = path

    except DataUnavailableError as e:
        res.status = "failed"
        res.message = str(e)
    except Exception as e:
        res.status = "failed"
        res.message = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"

    return res


# ── summary printer ───────────────────────────────────────────────────────────

ICONS = {"downloaded": "✓", "skipped": "—", "failed": "✗", "pending": "?"}
COLORS = {
    "downloaded": "\033[32m",   # green
    "skipped":    "\033[33m",   # yellow
    "failed":     "\033[31m",   # red
}
RESET = "\033[0m"


def print_results(results: list[Result]) -> None:
    print("\n" + "=" * 90)
    print(f"{'Accession':<25} {'DB':<8} {'Status':<12} {'Output / Message'}")
    print("=" * 90)
    for r in results:
        icon = ICONS.get(r.status, "?")
        color = COLORS.get(r.status, "")
        out = str(r.output_path) if r.output_path else r.message[:55]
        print(f"{color}{r.record.accession:<25} {r.record.db_type:<8} {icon} {r.status:<10} {out}{RESET}")
    print("=" * 90)

    downloaded = [r for r in results if r.status == "downloaded"]
    skipped    = [r for r in results if r.status == "skipped"]
    failed     = [r for r in results if r.status == "failed"]

    print(f"\nSummary: {len(downloaded)} downloaded, {len(skipped)} skipped, {len(failed)} failed")

    if failed:
        print("\nFailed datasets:")
        for r in failed:
            print(f"  {r.record.accession}: {r.message}")

    if skipped:
        controlled = [r for r in skipped if "Controlled" in r.message]
        if controlled:
            print(f"\nControlled-access datasets (require manual application):")
            for r in controlled:
                print(f"  {r.record.accession} ({r.record.db_type}): {r.record.url or 'no URL'}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Download scRNA-seq datasets from a CSV/xlsx catalog.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--catalog",    required=True, help="Path to .csv or .xlsx catalog file")
    parser.add_argument("--output-dir", default="data/downloaded", help="Directory to save h5ad files (default: data/downloaded)")
    parser.add_argument("--accession",  help="Download only this accession (e.g. GSE158702)")
    parser.add_argument("--dry-run",    action="store_true", help="Parse catalog and show what would be downloaded, without downloading")
    parser.add_argument("--hubmap-token", default=None, help="HuBMAP API token for protected datasets")
    parser.add_argument("--skip-existing", action="store_true", default=True, help="Skip accessions already downloaded (default: True)")
    args = parser.parse_args()

    # Load catalog
    print(f"Loading catalog: {args.catalog}")
    records = load_catalog(args.catalog)
    print_summary(records)

    # Filter to single accession if requested
    if args.accession:
        records = [r for r in records if r.accession == args.accession]
        if not records:
            print(f"ERROR: accession '{args.accession}' not found in catalog.")
            sys.exit(1)

    out_dir = Path(args.output_dir)

    if args.dry_run:
        print("\n[DRY RUN] Would attempt to download:")
        for r in records:
            if r.accessible:
                dest = out_dir / f"{r.accession}.h5ad"
                exists = " (already exists)" if dest.exists() else ""
                print(f"  {r.accession:<25} {r.db_type:<8} → {dest}{exists}")
            else:
                print(f"  {r.accession:<25} {r.db_type:<8} SKIP — {r.accessible=}")
        return

    # Download
    results = []
    for i, record in enumerate(records):
        print(f"\n[{i+1}/{len(records)}] {record.accession} — {record.title[:60]}")

        # Skip if already exists
        if args.skip_existing:
            dest = out_dir / f"{record.accession}.h5ad"
            if dest.exists():
                res = Result(record)
                res.status = "downloaded"
                res.output_path = dest
                res.message = "already exists"
                results.append(res)
                print(f"  Already exists: {dest}")
                continue

        res = acquire_one(record, out_dir, hubmap_token=args.hubmap_token)
        results.append(res)

    print_results(results)


if __name__ == "__main__":
    main()
