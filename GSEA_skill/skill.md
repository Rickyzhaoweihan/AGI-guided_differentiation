# Skill: gsea-analysis

Run GSEA and differential gene expression analysis on in-silico perturbation predictions.

## Input

The output of the `run-perturbation` skill:
- `perturbation_predictions.h5ad` — predicted gene expression for each cell under each perturbation + non-targeting control
- OR the per-gene CSVs from `extract_to_csv.py`

## What this skill does

For each perturbed gene, computes **predicted log2 fold change** (perturbed vs. non-targeting control) across all cells or per cell type, then runs **pre-ranked GSEA** using `gseapy` against MSigDB gene sets (Hallmark, GO, KEGG).

## Environment

```bash
conda activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
pip install gseapy  # if not already installed
```

## Run GSEA on perturbation predictions

```bash
python GSEA_skill/run_gsea.py \
    --predictions insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --output-dir GSEA_skill/results \
    --gene-sets MSigDB_Hallmark_2020   # or GO_Biological_Process_2023, KEGG_2021_Human
```

### Run on a specific perturbation gene
```bash
python GSEA_skill/run_gsea.py \
    --predictions insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --output-dir GSEA_skill/results \
    --gene-sets MSigDB_Hallmark_2020 \
    --target-gene TBX5
```

### Run per cell type
```bash
python GSEA_skill/run_gsea.py \
    --predictions insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --output-dir GSEA_skill/results \
    --gene-sets MSigDB_Hallmark_2020 \
    --per-cell-type
```

## Output structure

```
GSEA_skill/results/
├── TBX5/
│   ├── log2fc.csv                  ← ranked gene list (gene, log2fc, mean_pert, mean_ctrl)
│   ├── gsea_MSigDB_Hallmark_2020/  ← gseapy output (enrichment scores, plots)
│   └── volcano.png                 ← volcano plot of predicted DE
├── GATA4/
│   └── ...
└── summary_all_perturbations.csv   ← top pathway per perturbation, all genes
```

## How log2FC is computed

The predictions are **log-normalized expression** (log1p scale). Log2FC is computed as:

```
log2FC(gene) = mean_perturbed(gene) - mean_control(gene)   # difference in log space ≈ log2FC
```

This gives a ranked list of all 18,080 genes for pre-ranked GSEA.

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gene-sets` | `MSigDB_Hallmark_2020` | Gene set library (any Enrichr library name) |
| `--target-gene` | all | Run for one specific perturbation |
| `--per-cell-type` | False | Separate GSEA per cell type |
| `--min-abs-lfc` | 0.1 | Min log2FC to label on volcano plot |
| `--n-top-pathways` | 20 | Pathways to show in summary |

## Available gene set libraries (Enrichr)

- `MSigDB_Hallmark_2020` — 50 well-defined hallmark gene sets
- `GO_Biological_Process_2023` — GO biological process terms
- `KEGG_2021_Human` — KEGG pathways
- `Reactome_2022` — Reactome pathways
- `TF_Perturbations_Followed_by_Expression` — TF perturbation signatures (directly relevant)

## Connecting the two skills

```bash
# Step 1: Run perturbation pipeline
python insilico_perturbation_pipeline/scripts/run_insilico_perturbation.py \
    --input  insilico_perturbation_pipeline/data/inference_template.h5ad \
    --output insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --model-dir  insilico_perturbation_pipeline/models \
    --checkpoint "insilico_perturbation_pipeline/models/checkpoints/step=step=18000-val_loss=val_loss=1.7692.ckpt"

# Step 2: Run GSEA on the predictions
python GSEA_skill/run_gsea.py \
    --predictions insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --output-dir  GSEA_skill/results
```
