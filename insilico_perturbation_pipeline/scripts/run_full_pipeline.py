#!/usr/bin/env python3
"""Config-driven perturbation pipeline: prepare all datasets × cell types, generate SLURM jobs.

Reads a YAML config specifying datasets, cell types, filtering params, and model paths.
For each (dataset, cell_type): subsets h5ad, filters DEGs, batches genes, generates SLURM scripts.

Usage:
    # Prepare all datasets + generate SLURM scripts
    python run_full_pipeline.py prepare --config configs/collab_filtered_v1.yaml

    # Submit SLURM jobs with rate limiting
    python run_full_pipeline.py submit --config configs/collab_filtered_v1.yaml

    # Prepare a single dataset
    python run_full_pipeline.py prepare --config configs/collab_filtered_v1.yaml --dataset adult_gut

    # Dry run (show what would be done without executing)
    python run_full_pipeline.py prepare --config configs/collab_filtered_v1.yaml --dry-run
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

try:
    import yaml
except ImportError:
    print("PyYAML required: pip install pyyaml")
    sys.exit(1)


SCRIPT_DIR = Path(__file__).parent.resolve()


def load_config(config_path):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    # Resolve checkpoint relative to model_dir
    model = cfg["model"]
    if not os.path.isabs(model["checkpoint"]):
        model["checkpoint"] = os.path.join(model["model_dir"], model["checkpoint"])

    return cfg


def safe_name(cell_type):
    return cell_type.replace(" ", "_").replace("/", "_").lower()


# ── Prepare ──────────────────────────────────────────────────────────────────


def prepare(cfg, dataset_filter=None, dry_run=False):
    """Run prepare_celltype_run.py for each (dataset, cell_type), then generate SLURM jobs."""
    data_dir = cfg["data_dir"]
    output_dir = cfg["output_dir"]
    defaults = cfg["defaults"]
    prepare_script = SCRIPT_DIR / "prepare_celltype_run.py"
    preprocess_script = SCRIPT_DIR / "preprocess_for_inference.py"

    if not prepare_script.exists():
        print(f"ERROR: {prepare_script} not found")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "slurm_jobs"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "slurm_logs"), exist_ok=True)

    datasets = cfg["datasets"]
    if dataset_filter:
        datasets = {k: v for k, v in datasets.items() if k in dataset_filter}

    total_ct = sum(len(d["cell_types"]) for d in datasets.values())
    print(f"Pipeline: {len(datasets)} datasets, {total_ct} cell types")
    print(f"Output: {output_dir}")
    print()

    # Phase 1: prepare all cell types
    for ds_name, ds_cfg in datasets.items():
        h5ad = os.path.join(data_dir, ds_cfg["h5ad"])
        deg_csv = os.path.join(data_dir, ds_cfg["deg_csv"])
        ct_col = ds_cfg.get("cell_type_col", defaults["cell_type_col"])
        n_cells = ds_cfg.get("n_cells", defaults["n_cells"])
        batch_size = ds_cfg.get("batch_size", defaults["batch_size"])
        min_log2fc = ds_cfg.get("min_log2fc", defaults["min_log2fc"])
        max_padj = ds_cfg.get("max_padj", defaults["max_padj"])
        seed = ds_cfg.get("seed", defaults["seed"])

        print(f">>> {ds_name.upper()} <<<")
        print(f"    h5ad: {h5ad}")
        print(f"    DEG:  {deg_csv}")

        for ct in ds_cfg["cell_types"]:
            ct_safe = safe_name(ct)
            ct_dir = os.path.join(output_dir, ds_name, ct_safe)

            cmd = [
                sys.executable, str(prepare_script),
                "--h5ad", h5ad,
                "--deg-csv", deg_csv,
                "--cell-type", ct,
                "--cell-type-col", ct_col,
                "--output-dir", ct_dir,
                "--min-log2fc", str(min_log2fc),
                "--max-padj", str(max_padj),
                "--batch-size", str(batch_size),
                "--seed", str(seed),
            ]
            # Pass --gene-names-csv if bundled gene_names.csv exists
            gene_names = SCRIPT_DIR.parent / "models" / "gene_names.csv"
            if gene_names.exists():
                cmd.extend(["--gene-names-csv", str(gene_names)])

            print(f"  {ds_name}/{ct_safe} ...", end=" ", flush=True)
            if dry_run:
                print("[dry run]")
            else:
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode == 0:
                    # Extract gene count from output
                    for line in result.stdout.splitlines():
                        if "Final:" in line:
                            print(line.strip())
                            break
                    else:
                        print("done")
                else:
                    print(f"FAILED (exit {result.returncode})")
                    for line in result.stderr.splitlines()[-3:]:
                        print(f"    {line}")
        print()

    # Phase 2: generate SLURM job scripts
    print("Generating SLURM job scripts...")
    job_count = generate_slurm_jobs(cfg, output_dir, str(preprocess_script), dry_run)
    print(f"Generated {job_count} SLURM job scripts in {output_dir}/slurm_jobs/")


def generate_slurm_jobs(cfg, output_dir, preprocess_script, dry_run=False):
    """Walk output dirs, generate one SLURM script per batch."""
    model = cfg["model"]
    slurm_cfg = cfg["slurm"]
    defaults = cfg["defaults"]
    slurm_dir = os.path.join(output_dir, "slurm_jobs")
    job_count = 0

    for ds_name in cfg["datasets"]:
        ds_dir = os.path.join(output_dir, ds_name)
        if not os.path.isdir(ds_dir):
            continue

        for ct_name in sorted(os.listdir(ds_dir)):
            ct_dir = os.path.join(ds_dir, ct_name)
            batches_dir = os.path.join(ct_dir, "batches")
            if not os.path.isdir(batches_dir):
                continue

            subset_h5ads = [f for f in os.listdir(ct_dir) if f.endswith("_subset.h5ad")]
            if not subset_h5ads:
                continue
            subset_h5ad = os.path.join(ct_dir, subset_h5ads[0])

            ct_col = cfg["datasets"][ds_name].get("cell_type_col", defaults["cell_type_col"])
            n_cells = cfg["datasets"][ds_name].get("n_cells", defaults["n_cells"])
            seed = cfg["datasets"][ds_name].get("seed", defaults["seed"])

            for gene_file in sorted(os.listdir(batches_dir)):
                if not gene_file.startswith("gene_list_batch") or not gene_file.endswith(".txt"):
                    continue

                batch_name = gene_file.replace("gene_list_", "").replace(".txt", "")
                gene_list = os.path.join(batches_dir, gene_file)
                template_out = os.path.join(batches_dir, f"template_{batch_name}.h5ad")
                pred_out = os.path.join(batches_dir, f"predictions_{batch_name}.h5ad")
                n_genes = sum(1 for _ in open(gene_list))
                job_name = f"{ds_name}_{ct_name}_{batch_name}"

                script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={slurm_cfg['partition']}
#SBATCH --time={slurm_cfg['time']}
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem={slurm_cfg['mem']}
#SBATCH --output={output_dir}/slurm_logs/{job_name}_%j.out
#SBATCH --error={output_dir}/slurm_logs/{job_name}_%j.err
#SBATCH --mail-type=FAIL

echo "Job: {job_name} | ${{SLURM_JOB_ID}} | $(date)"
echo "Dataset: {ds_name} | Cell type: {ct_name} | Batch: {batch_name} | Genes: {n_genes}"

source activate {slurm_cfg['conda_env']}
cd {slurm_cfg['state_repo']}

echo "Step 1: Preprocessing..."
python {preprocess_script} \\
    --input {subset_h5ad} \\
    --gene-list {gene_list} \\
    --output {template_out} \\
    --cell-type-col {ct_col} \\
    --n-cells {n_cells} \\
    --seed {seed}

if [ $? -ne 0 ]; then echo "Preprocessing failed!"; exit 1; fi

echo "Step 2: Inference with {model['name']}..."
state tx infer \\
    --adata {template_out} \\
    --output {pred_out} \\
    --model-dir {model['model_dir']} \\
    --checkpoint {model['checkpoint']} \\
    --pert-col target_gene \\
    --celltype-col cell_type \\
    --batch-col batch_var \\
    --seed {seed}

if [ $? -ne 0 ]; then echo "Inference failed!"; exit 1; fi

echo "Done! Output: {pred_out}"
echo "End: $(date)"
"""
                job_path = os.path.join(slurm_dir, f"{job_name}.sh")
                if not dry_run:
                    with open(job_path, "w") as f:
                        f.write(script_content)
                job_count += 1

    return job_count


# ── Submit ───────────────────────────────────────────────────────────────────


def submit(cfg, dry_run=False):
    """Submit SLURM jobs with rate limiting."""
    output_dir = cfg["output_dir"]
    slurm_dir = os.path.join(output_dir, "slurm_jobs")
    max_concurrent = cfg["slurm"].get("max_concurrent", 4)

    jobs = sorted(Path(slurm_dir).glob("*.sh"))
    if not jobs:
        print(f"No SLURM scripts found in {slurm_dir}/")
        return

    # Skip jobs that already have predictions
    pending = []
    for job in jobs:
        # Infer prediction path from job name
        parts = job.stem.rsplit("_", 1)  # e.g., adult_gut_eecs_batch000
        # Check if predictions exist by looking for the prediction file
        # Parse the script to find the prediction output path
        content = job.read_text()
        pred_line = [l for l in content.splitlines() if "--output" in l and "predictions_" in l]
        if pred_line:
            pred_path = pred_line[0].strip().rstrip(" \\")
            if os.path.exists(pred_path):
                print(f"  SKIP {job.stem} (predictions exist)")
                continue
        pending.append(job)

    total = len(pending)
    skipped = len(jobs) - total
    if skipped:
        print(f"Skipping {skipped} jobs with existing predictions")
    print(f"Submitting {total} jobs, {max_concurrent} at a time")
    if dry_run:
        for j in pending:
            print(f"  [dry run] sbatch {j}")
        return

    i = 0
    while i < total:
        batch_end = min(i + max_concurrent, total)
        job_ids = []

        print(f"\n--- Batch {i + 1}-{batch_end} of {total} ---")
        for j in range(i, batch_end):
            result = subprocess.run(
                ["sbatch", "--parsable", str(pending[j])],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                jid = result.stdout.strip()
                job_ids.append(jid)
                print(f"  Submitted: {pending[j].stem} (job {jid})")
            else:
                print(f"  FAILED: {pending[j].stem}: {result.stderr.strip()}")

        # Wait for all jobs in this batch
        print(f"  Waiting for {len(job_ids)} jobs...")
        for jid in job_ids:
            while True:
                check = subprocess.run(
                    ["squeue", "-j", jid], capture_output=True, text=True,
                )
                if jid not in check.stdout:
                    break
                time.sleep(60)
            print(f"  Job {jid} finished")

        i = batch_end

    print(f"\nAll {total} jobs complete!")


# ── Status ───────────────────────────────────────────────────────────────────


def status(cfg):
    """Show pipeline status: which batches have predictions."""
    output_dir = cfg["output_dir"]
    total = 0
    done = 0
    missing = []

    for ds_name in cfg["datasets"]:
        ds_dir = os.path.join(output_dir, ds_name)
        if not os.path.isdir(ds_dir):
            continue
        for ct_name in sorted(os.listdir(ds_dir)):
            batches_dir = os.path.join(ds_dir, ct_name, "batches")
            if not os.path.isdir(batches_dir):
                continue
            for f in sorted(os.listdir(batches_dir)):
                if f.startswith("gene_list_batch") and f.endswith(".txt"):
                    batch = f.replace("gene_list_", "").replace(".txt", "")
                    pred = os.path.join(batches_dir, f"predictions_{batch}.h5ad")
                    total += 1
                    if os.path.exists(pred):
                        done += 1
                    else:
                        missing.append(f"{ds_name}/{ct_name}/{batch}")

    print(f"Pipeline status: {done}/{total} batches complete")
    if missing:
        print(f"\nMissing predictions ({len(missing)}):")
        for m in missing[:20]:
            print(f"  {m}")
        if len(missing) > 20:
            print(f"  ... and {len(missing) - 20} more")


# ── Main ─────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Config-driven perturbation pipeline orchestrator"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # prepare
    p_prep = sub.add_parser("prepare", help="Prepare all cell types + generate SLURM jobs")
    p_prep.add_argument("--config", "-c", required=True, help="YAML config file")
    p_prep.add_argument("--dataset", "-d", nargs="*", help="Only process these datasets")
    p_prep.add_argument("--dry-run", action="store_true", help="Show what would be done")

    # submit
    p_sub = sub.add_parser("submit", help="Submit SLURM jobs with rate limiting")
    p_sub.add_argument("--config", "-c", required=True, help="YAML config file")
    p_sub.add_argument("--dry-run", action="store_true", help="Show what would be submitted")

    # status
    p_stat = sub.add_parser("status", help="Show pipeline completion status")
    p_stat.add_argument("--config", "-c", required=True, help="YAML config file")

    args = parser.parse_args()
    cfg = load_config(args.config)

    if args.command == "prepare":
        prepare(cfg, dataset_filter=args.dataset, dry_run=args.dry_run)
    elif args.command == "submit":
        submit(cfg, dry_run=args.dry_run)
    elif args.command == "status":
        status(cfg)


if __name__ == "__main__":
    main()
