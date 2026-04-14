# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **multi-agent system for AGI-guided cell differentiation**, centered around an **in-silico gene perturbation prediction pipeline**. Given scRNA-seq data and target genes (from DEG lists or curated TF sets), it predicts how each gene perturbation would affect cell state using trained State models.

Two modes of operation:
1. **Config-driven batch mode** — run all datasets × all cell types from a YAML config (recommended)
2. **Manual single-run mode** — run one gene list against one dataset

## Conda Environment

```bash
source activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
```

## Directory Structure

```
AGI-guided_differentiation/
├── insilico_perturbation_pipeline/
│   ├── scripts/
│   │   ├── run_full_pipeline.py         ← Config-driven orchestrator (prepare/submit/status)
│   │   ├── prepare_celltype_run.py      ← Subset h5ad + filter DEGs per cell type
│   │   ├── preprocess_for_inference.py  ← Normalize, align genes, create template
│   │   ├── run_insilico_perturbation.py ← Python wrapper for state tx infer
│   │   ├── run_inference.sh             ← SLURM job template (single run)
│   │   ├── extract_to_csv.py            ← Split predictions to CSV
│   │   ├── train_state_model.py         ← Optional: train your own model
│   │   └── create_esm_pert_features.py  ← Optional: generate ESM2 embeddings
│   ├── configs/
│   │   ├── collab_filtered_v1.yaml      ← Demo: 4 datasets, 46 cell types, 81 batches
│   │   ├── starter.toml                 ← Training config (subset data)
│   │   └── full_dataset.toml            ← Training config (full dataset)
│   ├── models/                          ← Bundled state_sm model + gene_names.csv
│   ├── data/                            ← Input/output data
│   └── references/models.md             ← Model registry
├── GSEA_skill/                          ← Pathway enrichment analysis
├── hirn_publication_retrieval/          ← HIRN literature retrieval
└── PLAN.md                             ← Multi-agent architecture plan
```

## Running the Pipeline

### Config-Driven Batch Mode (Recommended)

For running perturbation predictions across all cell types in multiple datasets:

```bash
# 1. Prepare all datasets × cell types + generate SLURM jobs
python insilico_perturbation_pipeline/scripts/run_full_pipeline.py prepare \
    --config insilico_perturbation_pipeline/configs/collab_filtered_v1.yaml

# 2. Check status
python insilico_perturbation_pipeline/scripts/run_full_pipeline.py status \
    --config insilico_perturbation_pipeline/configs/collab_filtered_v1.yaml

# 3. Submit SLURM jobs (rate-limited, skips completed batches)
python insilico_perturbation_pipeline/scripts/run_full_pipeline.py submit \
    --config insilico_perturbation_pipeline/configs/collab_filtered_v1.yaml

# Prepare a single dataset only
python insilico_perturbation_pipeline/scripts/run_full_pipeline.py prepare \
    --config insilico_perturbation_pipeline/configs/collab_filtered_v1.yaml \
    --dataset adult_gut

# Dry run (show what would be done)
python insilico_perturbation_pipeline/scripts/run_full_pipeline.py prepare \
    --config insilico_perturbation_pipeline/configs/collab_filtered_v1.yaml --dry-run
```

See `configs/collab_filtered_v1.yaml` for config format. Key fields: `data_dir`, `output_dir`, `defaults` (n_cells, batch_size, filters), `model`, `slurm`, and `datasets` (h5ad, deg_csv, cell_types per dataset).

### Manual Single-Run Mode

For running one gene list against one dataset (collaborator-facing):

```bash
# Step 1: Preprocess
python insilico_perturbation_pipeline/scripts/preprocess_for_inference.py \
    --input <YOUR_DATA.h5ad> \
    --gene-list <YOUR_GENE_LIST.xlsx> \
    --output data/inference_template.h5ad \
    --cell-type-col <CELL_TYPE_COLUMN> \
    --n-cells 10000

# Step 2: Run inference via SLURM
sbatch insilico_perturbation_pipeline/scripts/run_inference.sh

# Step 3: Export to CSV
python insilico_perturbation_pipeline/scripts/extract_to_csv.py \
    --predictions data/perturbation_predictions.h5ad \
    --output-dir data/perturbation_csvs
```

## Architecture

The pipeline wraps the external `state` CLI tool (`python -m state tx infer`). Key design:

- **`run_full_pipeline.py`**: Reads YAML config, loops all datasets × cell types, calls prepare + preprocess for each, generates SLURM jobs, submits with rate limiting
- **`prepare_celltype_run.py`**: Subsets h5ad to one cell type, filters DEGs, batches genes into groups of 500
- **`preprocess_for_inference.py`**: Auto-detects normalization, stratified-samples cells, aligns to model's 18,080-gene format, creates inference template
- **`extract_to_csv.py`**: Splits large output h5ad into per-perturbation CSV files

**Bundled model**: `state_sm` (small, 40k steps, trained on Replogle-Nadig data) — works on any GPU. For `state_base` (H100 only), see `references/models.md`.

## Key Constraints

- Model expects exactly **18,080 genes** — preprocessing handles alignment automatically
- Matrices must be **dense** (not sparse); memory ≈ `n_cells × (n_perturbations + 1) × 18,080 × 4 bytes`
- Input h5ad must include cell type labels (specify column with `--cell-type-col`)
- Use `source activate` (not `conda activate`) for SLURM jobs
- `state` CLI must be installed via `pip install -e .` from the State repo

## GSEA Analysis

After predictions complete, run pathway enrichment:

```bash
python GSEA_skill/run_gsea.py \
    --predictions perturbation_predictions.h5ad \
    --gene-sets MSigDB_Hallmark_2020 \
    --output-dir results/gsea/
```

## Resource Requirements

| Step | RAM | GPU | Time |
|------|-----|-----|------|
| Prepare (per cell type) | 16 GB | None | Minutes |
| Preprocess (per batch) | 256 GB | None | Minutes |
| state_sm inference | 256 GB | Any GPU | 8–20 hrs |
| state_base inference | 256 GB | H100 only | 8–20 hrs |
| Export CSV | Low | None | Minutes |
