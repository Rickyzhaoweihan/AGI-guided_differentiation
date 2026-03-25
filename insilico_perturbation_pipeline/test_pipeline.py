#!/usr/bin/env python3
"""
End-to-end test for the in-silico perturbation pipeline.

Runs all three steps using the bundled example data:
  - example_fetal_heart.h5ad  (fetal heart scRNA-seq)
  - example_cardiac_tfs.csv   (20 cardiac transcription factors)

Step 1 (preprocess) and Step 3 (CSV export) run on CPU.
Step 2 (inference) requires a GPU — skipped if none detected.

Usage (from repo root):
    conda activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
    python insilico_perturbation_pipeline/test_pipeline.py
"""

import sys
import subprocess
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────

PIPELINE = Path(__file__).parent
SCRIPTS  = PIPELINE / "scripts"
DATA     = PIPELINE / "data"

INPUT_H5AD   = DATA / "example_fetal_heart.h5ad"
GENE_LIST    = DATA / "example_cardiac_tfs.csv"
TEMPLATE     = DATA / "test_inference_template.h5ad"
PREDICTIONS  = DATA / "test_perturbation_predictions.h5ad"
CSV_DIR      = DATA / "test_perturbation_csvs"
CHECKPOINT   = PIPELINE / "models" / "checkpoints" / "step=step=18000-val_loss=val_loss=1.7692.ckpt"
MODEL_DIR    = PIPELINE / "models"

# Add scripts to path so we can import directly
sys.path.insert(0, str(SCRIPTS))

# ── Helpers ──────────────────────────────────────────────────────────────────

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"
SKIP = "\033[93mSKIP\033[0m"

results = []

def check(name, condition, detail=""):
    status = PASS if condition else FAIL
    msg = f"  [{status}] {name}"
    if detail:
        msg += f" — {detail}"
    print(msg)
    results.append((name, condition))
    return condition


def has_gpu():
    try:
        result = subprocess.run(["nvidia-smi"], capture_output=True, timeout=10)
        return result.returncode == 0
    except Exception:
        return False


# ── Step 0: Preflight checks ─────────────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 0: Preflight checks")
print("=" * 60)

check("example_fetal_heart.h5ad exists", INPUT_H5AD.exists(), str(INPUT_H5AD))
check("example_cardiac_tfs.csv exists",  GENE_LIST.exists(),  str(GENE_LIST))
check("checkpoint exists",               CHECKPOINT.exists(), str(CHECKPOINT))

gpu_available = has_gpu()
print(f"  [{'  OK' if gpu_available else SKIP}] GPU detected: {gpu_available}")

# ── Step 1: Preprocess ───────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 1: Preprocess (CPU)")
print("=" * 60)

try:
    from preprocess_for_inference import (
        load_gene_list,
        load_target_genes,
        stratified_sample,
        align_genes_to_target,
        create_inference_template,
    )
    import anndata as ad

    # 1a. Load gene list
    genes = load_gene_list(str(GENE_LIST))
    check("load_gene_list", len(genes) == 20, f"{len(genes)} genes loaded")

    # 1b. Load model target genes (18,080)
    target_genes = load_target_genes()
    check("load_target_genes", len(target_genes) == 18080, f"{len(target_genes)} target genes")

    # 1c. Load h5ad
    print(f"\nLoading {INPUT_H5AD.name}...")
    adata = ad.read_h5ad(str(INPUT_H5AD))
    print(f"  Shape: {adata.shape}")
    check("load h5ad", adata.n_obs > 0, f"{adata.n_obs} cells, {adata.n_vars} genes")
    cell_type_col = "celltype"  # actual column name in example_fetal_heart.h5ad
    check("celltype column present", cell_type_col in adata.obs.columns,
          str(adata.obs.columns.tolist()))

    # 1d. Stratified sample
    sampled = stratified_sample(adata, n_cells=500, cell_type_col=cell_type_col, random_state=42)
    check("stratified_sample", sampled.n_obs <= 500, f"{sampled.n_obs} cells sampled")

    # 1e. Align genes
    X_aligned = align_genes_to_target(sampled, target_genes)
    check("align_genes_to_target", X_aligned.shape == (sampled.n_obs, 18080),
          f"shape={X_aligned.shape}")

    # 1f. Create template (use 3 TFs only to keep it fast for testing)
    test_genes = genes[:3]
    template = create_inference_template(
        sampled, X_aligned, target_genes,
        perturbation_targets=test_genes,
        cell_type_col=cell_type_col,
    )
    expected_rows = sampled.n_obs * (len(test_genes) + 1)
    check("create_inference_template", template.n_obs == expected_rows,
          f"shape={template.shape}")
    check("non-targeting control present",
          "non-targeting" in template.obs["target_gene"].values)
    check("all perturbation genes present",
          all(g in template.obs["target_gene"].values for g in test_genes))

    # 1g. Write full template for inference step
    print(f"\nBuilding full template ({len(genes)} TFs, 500 cells)...")
    full_template = create_inference_template(
        sampled, X_aligned, target_genes,
        perturbation_targets=genes,
        cell_type_col=cell_type_col,
    )
    TEMPLATE.parent.mkdir(exist_ok=True)
    full_template.write_h5ad(str(TEMPLATE))
    check("write inference template", TEMPLATE.exists(), str(TEMPLATE))

except Exception as e:
    print(f"\n  [{FAIL}] Preprocessing crashed: {e}")
    import traceback; traceback.print_exc()
    results.append(("preprocess", False))

# ── Step 2: Inference (GPU required) ─────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 2: Inference (GPU)")
print("=" * 60)

if not gpu_available:
    print(f"  [{SKIP}] No GPU detected — skipping inference step")
    print("         Submit test_example.sh via sbatch to run on a GPU node")
elif not TEMPLATE.exists():
    print(f"  [{SKIP}] Template not found (preprocessing failed) — skipping")
else:
    try:
        result = subprocess.run(
            [
                sys.executable, str(SCRIPTS / "run_insilico_perturbation.py"),
                "--input",      str(TEMPLATE),
                "--output",     str(PREDICTIONS),
                "--model-dir",  str(MODEL_DIR),
                "--checkpoint", str(CHECKPOINT),
                "--seed",       "42",
            ],
            capture_output=False,
            timeout=3600,
        )
        check("inference exit code 0", result.returncode == 0,
              f"exit code: {result.returncode}")
        check("predictions file created", PREDICTIONS.exists(), str(PREDICTIONS))

        if PREDICTIONS.exists():
            pred = ad.read_h5ad(str(PREDICTIONS))
            check("predictions shape sane", pred.n_obs > 0,
                  f"shape={pred.shape}")
            check("target_gene column present",
                  "target_gene" in pred.obs.columns)
    except subprocess.TimeoutExpired:
        print(f"  [{FAIL}] Inference timed out (>1 hour)")
        results.append(("inference", False))
    except Exception as e:
        print(f"  [{FAIL}] Inference crashed: {e}")
        results.append(("inference", False))

# ── Step 3: CSV export ───────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("STEP 3: CSV export (CPU)")
print("=" * 60)

if not PREDICTIONS.exists():
    print(f"  [{SKIP}] No predictions file — skipping CSV export")
    print("         (Run inference step first)")
else:
    try:
        from extract_to_csv import extract_perturbations_to_csv

        extract_perturbations_to_csv(str(PREDICTIONS), str(CSV_DIR))

        csv_files = list(CSV_DIR.glob("*.csv"))
        check("CSV files created", len(csv_files) > 0, f"{len(csv_files)} files")
        check("control_baseline.csv present",
              (CSV_DIR / "control_baseline.csv").exists())

        if csv_files:
            import pandas as pd
            sample_csv = csv_files[0]
            df = pd.read_csv(sample_csv, index_col=0, nrows=5)
            check("CSV has gene columns", df.shape[1] > 100,
                  f"{df.shape[1]} columns in {sample_csv.name}")

    except Exception as e:
        print(f"\n  [{FAIL}] CSV export crashed: {e}")
        import traceback; traceback.print_exc()
        results.append(("csv_export", False))

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
