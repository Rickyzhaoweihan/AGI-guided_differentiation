# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **self-contained in-silico gene perturbation prediction pipeline**. Given scRNA-seq data and a list of target genes, it predicts how each gene perturbation would affect cell state (gene expression) using a bundled pre-trained State model.

## Conda Environment

```bash
conda activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
```

## Directory Structure

```
insilico_perturbation_pipeline/
├── models/
│   └── checkpoints/
│       └── step=step=18000-val_loss=val_loss=1.7692.ckpt  ← bundled state_sm model
├── configs/
│   ├── starter.toml       ← training config (subset data)
│   └── full_dataset.toml  ← training config (full dataset)
├── data/                  ← place your .h5ad input here; outputs go here too
├── scripts/
│   ├── preprocess_for_inference.py   ← Step 1: prep data
│   ├── run_insilico_perturbation.py  ← Step 2: run inference (Python wrapper)
│   ├── run_inference.sh              ← Step 2 alt: SLURM job submission
│   ├── extract_to_csv.py             ← Step 3: split predictions to CSV
│   ├── train_state_model.py          ← optional: train your own model
│   └── create_esm_pert_features.py   ← optional: generate ESM2 embeddings
└── references/
    └── models.md          ← registry of available models and paths
```

## Running the Pipeline

### Step 1: Preprocess
```bash
python insilico_perturbation_pipeline/scripts/preprocess_for_inference.py \
    --input <YOUR_DATA.h5ad> \
    --gene-list <YOUR_GENE_LIST.xlsx> \
    --output insilico_perturbation_pipeline/data/inference_template.h5ad \
    --cell-type-col <CELL_TYPE_COLUMN> \
    --n-cells 10000 \
    --seed 42
```

### Step 2: Run Inference (Python)
```bash
python insilico_perturbation_pipeline/scripts/run_insilico_perturbation.py \
    --input insilico_perturbation_pipeline/data/inference_template.h5ad \
    --output insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --model-dir insilico_perturbation_pipeline/models \
    --checkpoint "insilico_perturbation_pipeline/models/checkpoints/step=step=18000-val_loss=val_loss=1.7692.ckpt"
```

Or via SLURM (edit paths in the script first):
```bash
sbatch insilico_perturbation_pipeline/scripts/run_inference.sh
```

### Step 3: Export to CSV
```bash
python insilico_perturbation_pipeline/scripts/extract_to_csv.py \
    --predictions insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --output-dir insilico_perturbation_pipeline/data/perturbation_csvs
```

## Architecture

The pipeline wraps the external `state` CLI tool (`python -m state tx infer`). Key design points:

- **`preprocess_for_inference.py`**: Auto-detects normalization, stratified-samples cells to preserve cell type proportions, zero-pads/reorders to the model's expected 18,080-gene format, creates one copy of cells per perturbation condition + non-targeting controls
- **`run_insilico_perturbation.py`**: Python wrapper that calls `state tx infer` with correct column names (`target_gene`, `cell_type`, `batch_var`)
- **`extract_to_csv.py`**: Splits the large output h5ad into one CSV per perturbation gene

**Bundled model**: `state_sm` (small, 40k steps, trained on Replogle-Nadig data) — works on any GPU. For `state_base` (H100 only), see `references/models.md` for the external path.

## Key Constraints

- Model expects exactly **18,080 genes** — preprocessing handles alignment automatically
- Matrices must be **dense** (not sparse); memory ≈ `n_cells × (n_perturbations + 1) × 18,080 × 4 bytes`
- Input h5ad must include cell type labels (specify column with `--cell-type-col`)
- Inference template must include non-targeting control cells (created automatically)

## Example Data

Two example files are bundled for a complete end-to-end test:
```
insilico_perturbation_pipeline/data/example_fetal_heart.h5ad   ← fetal heart scRNA-seq (Jiajun/Chen Lab, 28,892 cells)
insilico_perturbation_pipeline/data/example_cardiac_tfs.csv    ← 20 cardiac transcription factors (TBX5, GATA4, MEF2C, ...)
```
Cell type column for this dataset: `celltype` (no underscore)

## Resource Requirements

| Step | RAM | GPU | Time |
|------|-----|-----|------|
| Preprocess | 256 GB | None | Minutes |
| state_sm inference | 256 GB | Any GPU | 8–20 hrs |
| Export CSV | Low | None | Minutes |

## Skill

Use `/run-perturbation` to get step-by-step guidance running the pipeline end-to-end.
