#!/usr/bin/env python3
"""
test_pilot.py — Pilot test for the dataset acquisition skill.

Tests the full pipeline on two xlsx catalog files:
  1. Catalog parsing (both files)
  2. GEO file listing (no download) for 3 GEO accessions
  3. EBI file listing (no download) for 2 EBI accessions
  4. Actual download: 1 small GEO + 1 small EBI dataset
  5. Validate the downloaded h5ad files

Run:
    conda activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
    python agents/data_acquisition/test_pilot.py
"""

import re
import sys
import time
from pathlib import Path

# Add parent so we can import catalog / downloaders from here
sys.path.insert(0, str(Path(__file__).parent))

from catalog import load_catalog, print_summary
from downloaders.geo import _list_suppl_files, download as geo_download, DataUnavailableError
from downloaders.ebi import _list_study_files, download as ebi_download

ISLETS_XLSX = Path("Databasehumanislets.xlsx")
GUT_XLSX    = Path("ML for gut development.xlsx")
OUT_DIR     = Path("data/downloaded_pilot")

PASS = "\033[32m✓ PASS\033[0m"
FAIL = "\033[31m✗ FAIL\033[0m"
SKIP = "\033[33m— SKIP\033[0m"
SEP  = "-" * 70


def section(title: str):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


# ─────────────────────────────────────────────────────────────────────────────
# TEST 1: Catalog parsing
# ─────────────────────────────────────────────────────────────────────────────

def test_catalog_parsing():
    section("TEST 1 — Catalog parsing")
    results = {}

    for label, path in [("Islets", ISLETS_XLSX), ("Gut", GUT_XLSX)]:
        print(f"\n[{label}] {path}")
        try:
            records = load_catalog(path)
            public    = [r for r in records if r.accessible]
            private   = [r for r in records if not r.accessible]
            db_types  = {}
            for r in records:
                db_types[r.db_type] = db_types.get(r.db_type, 0) + 1

            print(f"  Total accessions : {len(records)}")
            print(f"  Public           : {len(public)}")
            print(f"  Controlled       : {len(private)}")
            print(f"  DB breakdown     : {db_types}")

            assert len(records) > 0, "No records parsed"
            assert len(public) > 0,  "No public records found"
            print(f"  {PASS}")
            results[label] = True
        except Exception as e:
            print(f"  {FAIL}: {e}")
            results[label] = False

    return results


# ─────────────────────────────────────────────────────────────────────────────
# TEST 2: GEO — list supplementary files (no download)
# ─────────────────────────────────────────────────────────────────────────────

# Pick 3 accessions of different ages/sizes from the islets catalog
GEO_PROBE = [
    ("GSE81608",  "2016 Smart-seq2, ~1700 cells — oldest, likely small"),
    ("GSE221237", "2023 10x, atlas — recent"),
    ("GSE158702", "2021 10x fetal gut — known to have h5ad"),
]

def test_geo_file_listing():
    section("TEST 2 — GEO supplementary file listing (no download)")
    results = {}

    for gse, note in GEO_PROBE:
        print(f"\n[GEO] {gse}  ({note})")
        try:
            files = _list_suppl_files(gse)
            if not files:
                print(f"  (no supplementary files found)")
                results[gse] = "empty"
            else:
                for f in files:
                    print(f"  {f['name']}")
                # Classify what's there
                has_h5ad    = any(f["name"].endswith((".h5ad", ".h5")) for f in files)
                has_10x     = any("matrix.mtx" in f["name"] for f in files)
                has_raw_tar = any(re.search(r"_RAW\.tar$", f["name"], re.I) for f in files)
                has_dge     = any(re.search(r"(count|dge|rpkm|fpkm|tpm)", f["name"], re.I) for f in files)
                tier = ("Tier1(h5ad)"    if has_h5ad else
                        "Tier2a(10x)"    if has_10x else
                        "Tier2c(RAW.tar)" if has_raw_tar else
                        "Tier2b(DGE)"    if has_dge else
                        "Tier4(FASTQ-only)")
                print(f"  → Would use: {tier}")
                results[gse] = tier
            print(f"  {PASS}")
        except Exception as e:
            print(f"  {FAIL}: {e}")
            results[gse] = "error"
        time.sleep(0.5)   # be polite to NCBI

    return results


# ─────────────────────────────────────────────────────────────────────────────
# TEST 3: EBI — list study files (no download)
# ─────────────────────────────────────────────────────────────────────────────

EBI_PROBE = [
    ("E-MTAB-5061", "2016 Smart-seq2 islets — oldest, likely small"),
    ("E-MTAB-8901", "2020 10x fetal gut"),
]

def test_ebi_file_listing():
    section("TEST 3 — EBI BioStudies file listing (no download)")
    results = {}

    for acc, note in EBI_PROBE:
        print(f"\n[EBI] {acc}  ({note})")
        try:
            files = _list_study_files(acc)
            print(f"  Total files: {len(files)}")
            # Show first 15
            for f in files[:15]:
                print(f"  {f['name']}")
            if len(files) > 15:
                print(f"  ... and {len(files)-15} more")

            has_h5ad = any(f["name"].endswith((".h5ad", ".h5")) for f in files)
            has_tsv  = any(
                ("count" in f["name"].lower() or "matrix" in f["name"].lower() or "expression" in f["name"].lower())
                and f["name"].endswith((".tsv", ".csv", ".txt", ".tsv.gz", ".txt.gz"))
                for f in files
            )
            tier = "Tier1(h5ad)" if has_h5ad else ("Tier3(TSV)" if has_tsv else "no processed data")
            print(f"  → Would use: {tier}")
            results[acc] = tier
            print(f"  {PASS}")
        except Exception as e:
            print(f"  {FAIL}: {e}")
            results[acc] = "error"

    return results


# ─────────────────────────────────────────────────────────────────────────────
# TEST 4: Actual download — 1 GEO + 1 EBI (small datasets)
# ─────────────────────────────────────────────────────────────────────────────

def test_actual_download(geo_tier_results: dict):
    section("TEST 4 — Actual download (small datasets)")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}

    # Pick the GEO accession that has a Tier1/Tier2 hit
    geo_candidates = [(gse, tier) for gse, tier in geo_tier_results.items()
                      if "Tier" in str(tier) and "FASTQ" not in str(tier)]
    # Prefer h5ad > 10x > DGE > RAW.tar (RAW.tar can be very large)
    tier_order = {"Tier1(h5ad)": 0, "Tier2a(10x)": 1, "Tier2b(DGE)": 2, "Tier2c(RAW.tar)": 3}
    geo_candidates.sort(key=lambda x: tier_order.get(x[1], 99))

    if not geo_candidates:
        print(f"\n[GEO] {SKIP} — no downloadable format found in probed accessions")
        results["geo"] = "skipped"
    else:
        gse, tier = geo_candidates[0]
        print(f"\n[GEO] Downloading {gse} ({tier}) ...")
        try:
            path = geo_download(gse, OUT_DIR)
            _validate_h5ad(gse, path)
            results["geo"] = "ok"
        except DataUnavailableError as e:
            print(f"  {FAIL}: {e}")
            results["geo"] = "failed"
        except Exception as e:
            import traceback
            print(f"  {FAIL}: {e}")
            traceback.print_exc()
            results["geo"] = "failed"

    # EBI: try E-MTAB-5061 (small Smart-seq2, 2209 cells)
    ebi_acc = "E-MTAB-5061"
    print(f"\n[EBI] Downloading {ebi_acc} ...")
    try:
        path = ebi_download(ebi_acc, OUT_DIR)
        _validate_h5ad(ebi_acc, path)
        results["ebi"] = "ok"
    except DataUnavailableError as e:
        print(f"  {FAIL}: {e}")
        results["ebi"] = "failed"
    except Exception as e:
        import traceback
        print(f"  {FAIL}: {e}")
        traceback.print_exc()
        results["ebi"] = "failed"

    return results


def _validate_h5ad(accession: str, path: Path):
    import anndata as ad
    print(f"  Validating {path.name} ...")
    adata = ad.read_h5ad(path)
    print(f"  Shape   : {adata.shape[0]} cells × {adata.shape[1]} genes")
    print(f"  obs cols: {adata.obs.columns.tolist()[:8]}")
    print(f"  var cols: {adata.var.columns.tolist()[:8]}")
    print(f"  X dtype : {adata.X.dtype}")
    size_mb = path.stat().st_size / 1024 / 1024
    print(f"  File    : {size_mb:.1f} MB")
    assert adata.shape[0] > 0, "0 cells"
    assert adata.shape[1] > 0, "0 genes"
    print(f"  {PASS}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("\nPILOT TEST — Dataset Acquisition Skill")
    print(f"Input files: {ISLETS_XLSX}, {GUT_XLSX}")
    print(f"Output dir : {OUT_DIR}")

    t0 = time.time()

    r1 = test_catalog_parsing()
    r2 = test_geo_file_listing()
    r3 = test_ebi_file_listing()
    r4 = test_actual_download(r2)

    elapsed = time.time() - t0

    section("SUMMARY")
    print(f"\nCatalog parsing   : {'OK' if all(r1.values()) else 'ISSUES'}")
    print(f"GEO file listing  : {r2}")
    print(f"EBI file listing  : {r3}")
    print(f"Downloads         : {r4}")
    print(f"\nElapsed: {elapsed:.1f}s")

    any_fail = (
        not all(r1.values()) or
        any(v == "error" for v in r2.values()) or
        any(v == "error" for v in r3.values()) or
        any(v == "failed" for v in r4.values())
    )
    if any_fail:
        print("\n\033[31mSome tests had issues — see above.\033[0m")
        sys.exit(1)
    else:
        print("\n\033[32mAll tests passed.\033[0m")


if __name__ == "__main__":
    main()
