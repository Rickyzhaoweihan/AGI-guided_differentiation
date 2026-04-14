#!/bin/bash
#SBATCH --job-name=hes1_synthesis_demo
#SBATCH --partition=drjieliu-l40s,drjieliu-v100
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --output=slurm_logs/hes1_demo_%j.out
#SBATCH --error=slurm_logs/hes1_demo_%j.err
#SBATCH --mail-type=END,FAIL

# Full synthesis demo: HES1 knockdown in adult gut stem cells
# Runs all 7 phases including literature search

echo "=============================================="
echo "Synthesis Demo: HES1 in Adult Gut Stem Cells"
echo "Job: ${SLURM_JOB_ID} | Node: ${SLURMD_NODENAME}"
echo "Start: $(date)"
echo "=============================================="

source activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env

cd /nfs/turbo/umms-drjieliu/proj/in_silico_pert/agi_guide_differentation_workspace/github_overall_workflow

PRED="/nfs/turbo/umms-drjieliu/proj/in_silico_pert/agi_guide_differentation_workspace/perturbation_workspace/runs/collab_filtered_v1/adult_gut/stem_cells/batches/predictions_batch000.h5ad"
OUTPUT="demo_results/hes1_stem_cells"
CONFIG="insilico_perturbation_pipeline/configs/expert_knowledge_example.yaml"

mkdir -p "$OUTPUT"
mkdir -p slurm_logs

echo ""
echo "Running full synthesis pipeline..."
echo "  Predictions: $PRED"
echo "  Config: $CONFIG"
echo "  Output: $OUTPUT"
echo ""

python insilico_perturbation_pipeline/scripts/analyze_perturbation.py \
    --predictions "$PRED" \
    --config "$CONFIG" \
    --target-gene HES1 \
    --output-dir "$OUTPUT" \
    --tissue "adult gut" \
    --cell-type "Stem cells"

echo ""
echo "=============================================="
echo "Demo complete!"
echo "Report: ${OUTPUT}/HES1/report.md"
echo "End: $(date)"
echo "=============================================="
