#!/usr/bin/env python3
"""
Test for the GSEA skill.

Uses the existing perturbation predictions from the insilico perturbation skill:
    insilico_perturbation_pipeline/data/test_perturbation_predictions.h5ad

Tests:
  - log2FC computation (CPU, fast)
  - Volcano plot generation (CPU, fast)
  - GSEA pre-ranked run on one gene (CPU, requires gseapy + internet for gene sets)

Usage (from repo root):
    conda activate /nfs/turbo/umms-drjieliu/usr/zheyyz/miniforge-pypy3/envs/state_env
    python GSEA_skill/test_gsea.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

REPO  = Path(__file__).parent.parent
PRED  = REPO / "insilico_perturbation_pipeline/data/test_perturbation_predictions.h5ad"
OUT   = REPO / "GSEA_skill/test_results"

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"
SKIP = "\033[93mSKIP\033[0m"

results = []

def check(name, condition, detail=""):
    status = PASS if condition else FAIL
    print(f"  [{status}] {name}" + (f" — {detail}" if detail else ""))
    results.append((name, condition))
    return condition


# ── Step 0: Preflight ─────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 0: Preflight")
print("=" * 60)

check("predictions file exists", PRED.exists(), str(PRED))

try:
    import gseapy
    gsea_available = True
    check("gseapy installed", True, gseapy.__version__)
except ImportError:
    gsea_available = False
    print(f"  [{SKIP}] gseapy not installed — GSEA step will be skipped")
    print("           Run: pip install gseapy")

if not PRED.exists():
    print(f"\nNo predictions file found. Run the perturbation test first:")
    print("  python insilico_perturbation_pipeline/test_pipeline.py")
    sys.exit(1)

# ── Step 1: Load predictions ──────────────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 1: Load predictions")
print("=" * 60)

import anndata as ad
import numpy as np

pred = ad.read_h5ad(str(PRED))
print(f"  Shape: {pred.shape}")

check("predictions loaded",       pred.n_obs > 0,   f"{pred.n_obs:,} cells × {pred.n_vars:,} genes")
check("target_gene column",       "target_gene"      in pred.obs.columns)
check("non-targeting control",    "non-targeting"    in pred.obs["target_gene"].values)
check("cell_type column",         "cell_type"        in pred.obs.columns)

perturbations = [g for g in pred.obs["target_gene"].unique() if g != "non-targeting"]
check("perturbation genes present", len(perturbations) > 0, f"{len(perturbations)} perturbations: {perturbations[:4]}...")

# ── Step 2: log2FC computation ────────────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 2: log2FC computation")
print("=" * 60)

from run_gsea import compute_log2fc

test_gene = perturbations[0]
print(f"  Testing on: {test_gene}")

ranked = compute_log2fc(pred, test_gene)

check("log2fc computed",              ranked is not None)
check("has all 18080 genes",          ranked is not None and len(ranked) == pred.n_vars,
      f"{len(ranked) if ranked is not None else 0} genes")
check("log2fc values are finite",     ranked is not None and np.isfinite(ranked["log2fc"]).all())
check("genes ranked descending",      ranked is not None and ranked["log2fc"].iloc[0] >= ranked["log2fc"].iloc[-1])
check("has upregulated genes",        ranked is not None and (ranked["log2fc"] > 0).any())
check("has downregulated genes",      ranked is not None and (ranked["log2fc"] < 0).any())

if ranked is not None:
    top5 = ranked.head(5).index.tolist()
    bot5 = ranked.tail(5).index.tolist()
    print(f"  Top upregulated:   {top5}")
    print(f"  Top downregulated: {bot5}")

# ── Step 3: Per-cell-type log2FC ──────────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 3: Per-cell-type log2FC")
print("=" * 60)

cell_types = sorted(pred.obs["cell_type"].unique())
print(f"  Cell types: {cell_types}")

ct = cell_types[0]
ranked_ct = compute_log2fc(pred, test_gene, cell_type=ct)
check(f"log2fc for cell type '{ct}'", ranked_ct is not None,
      f"{len(ranked_ct) if ranked_ct is not None else 0} genes")

# ── Step 4: Volcano plot ──────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 4: Volcano plot")
print("=" * 60)

OUT.mkdir(parents=True, exist_ok=True)
volcano_path = OUT / f"{test_gene}_volcano.png"

try:
    from run_gsea import plot_volcano
    plot_volcano(ranked, test_gene, volcano_path, min_abs_lfc=0.1)
    check("volcano plot saved", volcano_path.exists(), str(volcano_path))
except Exception as e:
    print(f"  [{FAIL}] volcano plot failed: {e}")
    results.append(("volcano plot", False))

# ── Step 5: Save log2FC CSV ───────────────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 5: Save log2FC CSV")
print("=" * 60)

import pandas as pd

csv_path = OUT / f"{test_gene}_log2fc.csv"
try:
    ranked.to_csv(csv_path)
    loaded = pd.read_csv(csv_path, index_col=0)
    check("log2fc CSV saved",    csv_path.exists(), str(csv_path))
    check("CSV has correct cols", all(c in loaded.columns for c in ["log2fc", "mean_pert", "mean_ctrl"]))
    check("CSV row count",        len(loaded) == pred.n_vars, f"{len(loaded)} rows")
except Exception as e:
    print(f"  [{FAIL}] CSV save failed: {e}")
    results.append(("log2fc CSV", False))

# ── Step 6: GSEA (requires gseapy + internet) ────────────────────────────────

print("\n" + "=" * 60)
print("STEP 6: GSEA pre-ranked")
print("=" * 60)

if not gsea_available:
    print(f"  [{SKIP}] gseapy not installed")
else:
    try:
        from run_gsea import run_gsea_preranked
        gsea_out = OUT / test_gene
        gsea_out.mkdir(parents=True, exist_ok=True)

        print(f"  Running GSEA for {test_gene} against MSigDB_Hallmark_2020...")
        gsea_res = run_gsea_preranked(ranked, "MSigDB_Hallmark_2020", gsea_out, test_gene)

        check("GSEA returned results",    gsea_res is not None and len(gsea_res) > 0,
              f"{len(gsea_res)} pathways")
        check("GSEA has NES column",      "NES" in gsea_res.columns)
        check("GSEA has FDR column",      "FDR q-val" in gsea_res.columns)

        sig = gsea_res[gsea_res["FDR q-val"] < 0.25]
        print(f"  Significant pathways (FDR < 0.25): {len(sig)}")
        if len(sig) > 0:
            top = sig.sort_values("NES", ascending=False).iloc[0]
            print(f"  Top pathway: {top.name}  NES={top['NES']:.2f}  FDR={top['FDR q-val']:.3f}")

    except Exception as e:
        print(f"  [{FAIL}] GSEA failed: {e}")
        import traceback; traceback.print_exc()
        results.append(("gsea preranked", False))

# ── Summary ───────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)

passed = sum(1 for _, ok in results if ok)
total  = len(results)
for name, ok in results:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

print(f"\n{passed}/{total} checks passed")
if passed == total:
    print("All checks passed!")
else:
    print("Some checks failed — see output above.")
    sys.exit(1)
