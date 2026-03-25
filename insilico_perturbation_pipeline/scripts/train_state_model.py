#!/usr/bin/env python3
"""
Train STATE model with ESM2 perturbation features.

This script trains the STATE (State Transition) model for perturbation prediction
using ESM2 protein embeddings to encode genetic perturbations.

Model variants:
- state_sm: Small model (~162M parameters), faster training, suitable for testing
- state:    Base model (~500M parameters), balanced performance
- state_lg: Large model (~1.1B parameters), best quality, requires more compute

Usage:
    # Quick test with competition subset
    python train_state_model.py --use-starter --model state_sm --max-steps 1000

    # Full training with small model
    python train_state_model.py --use-starter --model state_sm --max-steps 40000

    # Production training with large model
    python train_state_model.py --model state_lg --max-steps 200000
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path


# Default paths - adjust these for your setup
SCRIPT_DIR = Path(__file__).parent.resolve()
DATA_DIR = SCRIPT_DIR.parent / "data"
CONFIGS_DIR = SCRIPT_DIR.parent / "configs"
RESULTS_DIR = SCRIPT_DIR.parent / "results"


def check_gpu():
    """Check if GPU is available."""
    try:
        import torch
        if torch.cuda.is_available():
            return True, torch.cuda.get_device_name(0)
        return False, None
    except ImportError:
        return None, "torch not installed"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train STATE model with ESM perturbation features",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Data configuration
    data_group = parser.add_argument_group('Data Configuration')
    data_group.add_argument(
        "--use-starter",
        action="store_true",
        help="Use starter.toml (competition subset, ~40K cells) for faster testing"
    )
    data_group.add_argument(
        "--config",
        type=str,
        default=None,
        help="Custom path to TOML config file"
    )
    data_group.add_argument(
        "--pert-features",
        type=str,
        default=None,
        help="Path to ESM2 perturbation features file (.pt)"
    )
    data_group.add_argument(
        "--pert-col",
        type=str,
        default="target_gene",
        help="Column name for perturbations (default: target_gene)"
    )
    data_group.add_argument(
        "--control-pert",
        type=str,
        default="non-targeting",
        help="Control perturbation name (default: non-targeting)"
    )
    data_group.add_argument(
        "--cell-type-key",
        type=str,
        default="cell_type",
        help="Column name for cell types (default: cell_type)"
    )
    data_group.add_argument(
        "--batch-col",
        type=str,
        default="batch_var",
        help="Column name for batch (default: batch_var)"
    )
    data_group.add_argument(
        "--num-workers",
        type=int,
        default=8,
        help="Number of data loading workers (default: 8)"
    )

    # Model configuration
    model_group = parser.add_argument_group('Model Configuration')
    model_group.add_argument(
        "--model",
        type=str,
        default="state_sm",
        choices=["state", "state_sm", "state_lg"],
        help="Model size: state_sm (~162M), state (~500M), state_lg (~1.1B)"
    )

    # Training configuration
    train_group = parser.add_argument_group('Training Configuration')
    train_group.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for checkpoints and logs"
    )
    train_group.add_argument(
        "--name",
        type=str,
        default=None,
        help="Run name for logging (auto-generated if not specified)"
    )
    train_group.add_argument(
        "--max-steps",
        type=int,
        default=40000,
        help="Maximum training steps (default: 40000)"
    )
    train_group.add_argument(
        "--ckpt-every-n-steps",
        type=int,
        default=10000,
        help="Save checkpoint every N steps (default: 10000)"
    )
    train_group.add_argument(
        "--val-freq",
        type=int,
        default=1000,
        help="Validate every N steps (default: 1000)"
    )
    train_group.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="Training batch size (default: from model config)"
    )
    train_group.add_argument(
        "--learning-rate",
        type=float,
        default=None,
        help="Learning rate (default: from model config)"
    )

    # W&B configuration
    wandb_group = parser.add_argument_group('Weights & Biases Configuration')
    wandb_group.add_argument(
        "--no-wandb",
        action="store_true",
        help="Disable Weights & Biases logging"
    )
    wandb_group.add_argument(
        "--wandb-project",
        type=str,
        default="state-perturbation",
        help="W&B project name"
    )
    wandb_group.add_argument(
        "--wandb-tags",
        type=str,
        nargs="+",
        default=["esm-features"],
        help="W&B tags"
    )

    # Additional options
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the command without executing it"
    )
    parser.add_argument(
        "hydra_overrides",
        nargs="*",
        help="Additional Hydra overrides (e.g., training.lr=1e-4)"
    )

    return parser.parse_args()


def get_config_path(args):
    """Determine which config file to use."""
    if args.config:
        return Path(args.config)
    elif args.use_starter:
        return CONFIGS_DIR / "starter.toml"
    else:
        return CONFIGS_DIR / "full_dataset.toml"


def get_pert_features_path(args):
    """Determine perturbation features file path."""
    if args.pert_features:
        return Path(args.pert_features)
    return DATA_DIR / "ESM2_pert_features.pt"


def get_output_dir(args):
    """Determine output directory."""
    if args.output_dir:
        return Path(args.output_dir)
    return RESULTS_DIR / f"{args.model}_output"


def verify_files(args, config_path, pert_features_path):
    """Verify that required files exist."""
    errors = []

    if not config_path.exists():
        errors.append(f"Config file not found: {config_path}")

    if not pert_features_path.exists():
        errors.append(f"Perturbation features file not found: {pert_features_path}")

    return errors


def build_command(args, config_path, pert_features_path, output_dir):
    """Build the state tx train command."""
    cmd_parts = ["state", "tx", "train"]

    # Data configuration
    cmd_parts.extend([
        f"data.kwargs.toml_config_path={config_path}",
        f"data.kwargs.perturbation_features_file={pert_features_path}",
        f"data.kwargs.num_workers={args.num_workers}",
        f"data.kwargs.batch_col={args.batch_col}",
        f"data.kwargs.pert_col={args.pert_col}",
        f"data.kwargs.cell_type_key={args.cell_type_key}",
        f"data.kwargs.control_pert={args.control_pert}",
    ])

    # Training configuration
    cmd_parts.extend([
        f"training.max_steps={args.max_steps}",
        f"training.ckpt_every_n_steps={args.ckpt_every_n_steps}",
        f"training.val_freq={args.val_freq}",
    ])

    if args.batch_size is not None:
        cmd_parts.append(f"training.batch_size={args.batch_size}")

    if args.learning_rate is not None:
        cmd_parts.append(f"training.lr={args.learning_rate}")

    # Model configuration
    cmd_parts.append(f"model={args.model}")

    # Generate run name if not provided
    run_name = args.name
    if run_name is None:
        dataset_type = "starter" if args.use_starter else "full"
        run_name = f"{args.model}_{dataset_type}_{int(time.time())}"

    # Output configuration
    cmd_parts.extend([
        f"output_dir={output_dir}",
        f"name={run_name}",
    ])

    # W&B configuration
    if args.no_wandb:
        cmd_parts.append("use_wandb=false")
    else:
        cmd_parts.append("use_wandb=true")
        cmd_parts.append(f"wandb.project={args.wandb_project}")
        if args.wandb_tags:
            tags_str = "[" + ",".join(args.wandb_tags) + "]"
            cmd_parts.append(f"wandb.tags={tags_str}")

    # Add any additional overrides
    if args.hydra_overrides:
        cmd_parts.extend(args.hydra_overrides)

    return cmd_parts, run_name


def print_config(args, config_path, pert_features_path, output_dir, run_name):
    """Print configuration summary."""
    print("=" * 70)
    print("STATE TRAINING CONFIGURATION")
    print("=" * 70)

    gpu_available, gpu_name = check_gpu()
    if gpu_available:
        print(f"\nGPU: {gpu_name}")
    elif gpu_available is False:
        print(f"\nWARNING: No GPU detected! Training will be very slow.")

    print(f"\nDataset:")
    print(f"  Config: {config_path.name}")
    print(f"  Pert features: {pert_features_path.name}")

    print(f"\nModel: {args.model}")

    print(f"\nTraining:")
    print(f"  Run name:        {run_name}")
    print(f"  Max steps:       {args.max_steps:,}")
    print(f"  Val frequency:   {args.val_freq:,}")
    print(f"  Ckpt frequency:  {args.ckpt_every_n_steps:,}")
    print(f"  Output dir:      {output_dir}")
    if args.batch_size:
        print(f"  Batch size:      {args.batch_size}")

    print(f"\nW&B: {'Disabled' if args.no_wandb else args.wandb_project}")
    print("=" * 70)


def main():
    args = parse_args()

    config_path = get_config_path(args)
    pert_features_path = get_pert_features_path(args)
    output_dir = get_output_dir(args)

    # Verify files exist
    errors = verify_files(args, config_path, pert_features_path)
    if errors:
        print("\nError: Missing required files:\n")
        for error in errors:
            print(f"  {error}")
        sys.exit(1)

    # Build command
    cmd_parts, run_name = build_command(args, config_path, pert_features_path, output_dir)

    # Print configuration
    print_config(args, config_path, pert_features_path, output_dir, run_name)

    # Print command
    print("\nCommand:")
    print("-" * 70)
    print(" \\\n  ".join(cmd_parts))
    print("-" * 70)

    if args.dry_run:
        print("\nDry run complete. Command not executed.")
        return 0

    # Execute command
    print("\nStarting training...\n")
    start_time = time.time()

    result = subprocess.run(cmd_parts)

    elapsed = time.time() - start_time
    elapsed_str = f"{elapsed/3600:.1f}h" if elapsed > 3600 else f"{elapsed/60:.1f}m"

    print("\n" + "=" * 70)
    if result.returncode == 0:
        print(f"Training completed successfully! (Time: {elapsed_str})")
        print(f"Checkpoints saved to: {output_dir}/{run_name}")
    else:
        print(f"Training failed with exit code: {result.returncode}")
    print("=" * 70)

    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
