#!/bin/bash
#############################################################################
# SLURM Template: In-Silico Perturbation Inference
#
# Fill in the CONFIGURATION section below with your specific paths.
# Preprocessing step is optional (comment out if template already exists).
#############################################################################

#SBATCH --mail-type=BEGIN,END
#SBATCH --job-name=insilico_pert
#SBATCH --partition=drjieliu-h100,drjieliu-l40s,drjieliu-v100
#SBATCH --time=4:00:00
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem-per-gpu=128GB
#SBATCH --output=slurm_logs/insilico_pert_%j.out
#SBATCH --error=slurm_logs/insilico_pert_%j.err

#############################################################################
# CONFIGURATION - FILL THESE IN
#############################################################################

# Data paths
INPUT_H5AD=""           # e.g., /path/to/fetalgut.h5ad
GENE_LIST=""            # e.g., /path/to/Genelistgut.xlsx
TEMPLATE_OUTPUT=""      # e.g., /path/to/inference_template.h5ad
PREDICTION_OUTPUT=""    # e.g., /path/to/perturbation_predictions.h5ad

# Model paths (choose one)
MODEL_DIR=""            # e.g., output_state_sm_baseline/state_sm_40k_20251202_001509
CHECKPOINT=""           # e.g., ${MODEL_DIR}/checkpoints/step=step=18000-val_loss=val_loss=1.7692.ckpt

# Preprocessing params
CELL_TYPE_COL="celltype"  # Column name for cell type in .obs
N_CELLS=10000
SEED=42

# Path to preprocessing script (relative to repo root)
PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
PREPROCESS_SCRIPT="${PIPELINE_DIR}/scripts/preprocess_for_inference.py"

#############################################################################
# ENVIRONMENT SETUP
#############################################################################

echo "=============================================="
echo "In-Silico Perturbation Inference"
echo "Job ID: $SLURM_JOB_ID | Node: $SLURMD_NODENAME"
echo "Start: $(date)"
echo "=============================================="

cd /nfs/turbo/umms-drjieliu/usr/zheyuz/state_with_esm
mkdir -p slurm_logs

source activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env

echo "Python: $(which python)"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv
echo ""

#############################################################################
# STEP 1: PREPROCESSING (uncomment if needed)
#############################################################################

# python "$PREPROCESS_SCRIPT" \
#     --input "$INPUT_H5AD" \
#     --gene-list "$GENE_LIST" \
#     --output "$TEMPLATE_OUTPUT" \
#     --cell-type-col "$CELL_TYPE_COL" \
#     --n-cells $N_CELLS \
#     --seed $SEED
#
# if [ $? -ne 0 ]; then
#     echo "Preprocessing failed!"
#     exit 1
# fi

#############################################################################
# STEP 2: INFERENCE
#############################################################################

echo "Running inference..."
echo "  Template: $TEMPLATE_OUTPUT"
echo "  Model: $MODEL_DIR"

state tx infer \
    --adata "$TEMPLATE_OUTPUT" \
    --output "$PREDICTION_OUTPUT" \
    --model-dir "$MODEL_DIR" \
    --checkpoint "$CHECKPOINT" \
    --pert-col "target_gene" \
    --celltype-col "cell_type" \
    --batch-col "batch_var" \
    --seed $SEED

if [ $? -ne 0 ]; then
    echo "Inference failed!"
    exit 1
fi

echo ""
echo "Done! Output: $PREDICTION_OUTPUT"
echo "End: $(date)"
