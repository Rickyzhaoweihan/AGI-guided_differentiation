#!/usr/bin/env python3
"""
Run in-silico perturbation inference using trained STATE model.

This script:
1. Loads the preprocessed inference template
2. Runs inference using the trained STATE model
3. Saves the predicted perturbation effects

The output contains predicted gene expression for each cell under each
perturbation condition.

Usage:
    python run_insilico_perturbation.py \
        --input data/inference_template.h5ad \
        --output results/perturbation_predictions.h5ad \
        --model-dir results/state_sm_output/<run_name> \
        --checkpoint results/state_sm_output/<run_name>/checkpoints/best.ckpt
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_inference(args):
    """Run state tx infer command."""

    print("=" * 80)
    print("IN-SILICO PERTURBATION INFERENCE")
    print("=" * 80)
    print(f"\nInput template: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Model directory: {args.model_dir}")
    print(f"Checkpoint: {args.checkpoint}")

    # Build command
    cmd = [
        "python", "-m", "state", "tx", "infer",
        "--adata", args.input,
        "--output", args.output,
        "--model-dir", args.model_dir,
        "--checkpoint", args.checkpoint,
        "--pert-col", "target_gene",
        "--celltype-col", "cell_type",
        "--batch-col", "batch_var",
        "--seed", str(args.seed),
    ]

    if args.max_set_len:
        cmd.extend(["--max-set-len", str(args.max_set_len)])

    if args.quiet:
        cmd.append("--quiet")

    print(f"\nCommand:")
    print(" ".join(cmd))
    print("\n" + "=" * 80)
    print("Starting inference...")
    print("=" * 80 + "\n")

    try:
        result = subprocess.run(cmd, check=True)

        print("\n" + "=" * 80)
        print("INFERENCE COMPLETE")
        print("=" * 80)
        print(f"\nOutput saved to: {args.output}")
        print(f"\nNext steps:")
        print(f"  1. Load predictions: ad.read_h5ad('{args.output}')")
        print(f"  2. Compare predicted vs control expression")
        print(f"  3. Analyze perturbation effects by cell type")

        return 0

    except subprocess.CalledProcessError as e:
        print(f"\nInference failed with error code: {e.returncode}")
        return e.returncode
    except FileNotFoundError:
        print("\nCould not find 'state' command. Make sure the environment is activated.")
        return 1


def main():
    parser = argparse.ArgumentParser(
        description="Run in-silico perturbation inference using trained STATE model"
    )
    parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Input template h5ad file"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output predictions h5ad file"
    )
    parser.add_argument(
        "--model-dir",
        type=str,
        required=True,
        help="Model training directory"
    )
    parser.add_argument(
        "--checkpoint",
        type=str,
        required=True,
        help="Model checkpoint path"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--max-set-len",
        type=int,
        default=None,
        help="Maximum set length per forward pass (optional)"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce logging verbosity"
    )

    args = parser.parse_args()

    # Validate paths
    if not Path(args.input).exists():
        print(f"Input file not found: {args.input}")
        print(f"\nPlease run preprocessing first:")
        print(f"  python preprocess_for_inference.py --input your_data.h5ad --output {args.input}")
        sys.exit(1)

    if not Path(args.model_dir).exists():
        print(f"Model directory not found: {args.model_dir}")
        sys.exit(1)

    if not Path(args.checkpoint).exists():
        print(f"Checkpoint not found: {args.checkpoint}")
        sys.exit(1)

    # Create output directory
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Run inference
    exit_code = run_inference(args)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
