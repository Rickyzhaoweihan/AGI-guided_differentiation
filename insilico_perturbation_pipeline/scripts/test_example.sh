#!/bin/bash
#############################################################################
# Test script: runs the full pipeline on the bundled example data
# (fetal heart + 20 cardiac TFs, state_sm checkpoint)
#
# Submit with:
#   sbatch insilico_perturbation_pipeline/scripts/test_example.sh
#
# Run from repo root:
#   /nfs/turbo/umms-drjieliu/usr/rickyhan/AGI-guided_differentiation/
#############################################################################

#SBATCH --job-name=pert_test
#SBATCH --partition=drjieliu-h100,drjieliu-l40s,drjieliu-v100
#SBATCH --time=12:00:00
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem-per-gpu=128GB
#SBATCH --output=insilico_perturbation_pipeline/logs/test_%j.out
#SBATCH --error=insilico_perturbation_pipeline/logs/test_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL

set -e

REPO_ROOT="/nfs/turbo/umms-drjieliu/usr/rickyhan/AGI-guided_differentiation"
PIPELINE="${REPO_ROOT}/insilico_perturbation_pipeline"

INPUT_H5AD="${PIPELINE}/data/example_fetal_heart.h5ad"
GENE_LIST="${PIPELINE}/data/example_cardiac_tfs.csv"
TEMPLATE="${PIPELINE}/data/test_inference_template.h5ad"
PREDICTIONS="${PIPELINE}/data/test_perturbation_predictions.h5ad"
CHECKPOINT="${PIPELINE}/models/checkpoints/step=step=18000-val_loss=val_loss=1.7692.ckpt"

mkdir -p "${PIPELINE}/logs"

echo "=============================================="
echo "In-Silico Perturbation — Example Test Run"
echo "Job ID: $SLURM_JOB_ID | Node: $SLURMD_NODENAME"
echo "Start: $(date)"
echo "=============================================="

source activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
echo "Python: $(which python)"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv
echo ""

# Step 1: Preprocess
echo "[1/3] Preprocessing..."
python "${PIPELINE}/scripts/preprocess_for_inference.py" \
    --input "$INPUT_H5AD" \
    --gene-list "$GENE_LIST" \
    --output "$TEMPLATE" \
    --cell-type-col celltype \
    --n-cells 10000 \
    --seed 42
echo "Preprocessing done."
echo ""

# Step 2: Inference
echo "[2/3] Running inference..."
python "${PIPELINE}/scripts/run_insilico_perturbation.py" \
    --input "$TEMPLATE" \
    --output "$PREDICTIONS" \
    --model-dir "${PIPELINE}/models" \
    --checkpoint "$CHECKPOINT" \
    --seed 42
echo "Inference done."
echo ""

# Step 3: Export to CSV
echo "[3/3] Exporting predictions to CSV..."
python "${PIPELINE}/scripts/extract_to_csv.py" \
    --predictions "$PREDICTIONS" \
    --output-dir "${PIPELINE}/data/test_perturbation_csvs"
echo "Export done."

echo ""
echo "=============================================="
echo "All steps complete. End: $(date)"
echo "Outputs:"
echo "  Template:    $TEMPLATE"
echo "  Predictions: $PREDICTIONS"
echo "  CSVs:        ${PIPELINE}/data/test_perturbation_csvs/"
echo "=============================================="
