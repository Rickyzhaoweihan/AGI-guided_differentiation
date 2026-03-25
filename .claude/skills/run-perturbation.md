# Skill: run-perturbation

Run the in-silico gene perturbation pipeline end-to-end.

## What this skill does

When invoked, guide the user through the full pipeline:
1. **Preprocess** their scRNA-seq `.h5ad` + gene list → inference template
2. **Run inference** using the bundled `state_sm` checkpoint
3. **Export** predictions to per-gene CSVs

## Bundled checkpoint (self-contained)

```
insilico_perturbation_pipeline/models/checkpoints/step=step=18000-val_loss=val_loss=1.7692.ckpt
```
This is the `state_sm` model (40k training steps). No external model path needed.

## Step-by-step

### Activate environment
```bash
conda activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
```

### Step 1 – Preprocess
```bash
python insilico_perturbation_pipeline/scripts/preprocess_for_inference.py \
    --input <YOUR_DATA.h5ad> \
    --gene-list <YOUR_GENE_LIST.xlsx> \
    --output insilico_perturbation_pipeline/data/inference_template.h5ad \
    --cell-type-col <CELL_TYPE_COLUMN> \
    --n-cells 10000 \
    --seed 42
```
Accepts `.xlsx`, `.csv`, `.tsv`, or `.txt` for gene lists. Auto-detects normalization state.

### Step 2 – Run inference
```bash
python insilico_perturbation_pipeline/scripts/run_insilico_perturbation.py \
    --input insilico_perturbation_pipeline/data/inference_template.h5ad \
    --output insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --model-dir insilico_perturbation_pipeline/models \
    --checkpoint "insilico_perturbation_pipeline/models/checkpoints/step=step=18000-val_loss=val_loss=1.7692.ckpt"
```

Or via SLURM (edit paths first):
```bash
sbatch insilico_perturbation_pipeline/scripts/run_inference.sh
```

### Step 3 – Export to CSV
```bash
python insilico_perturbation_pipeline/scripts/extract_to_csv.py \
    --predictions insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --output-dir insilico_perturbation_pipeline/data/perturbation_csvs
```

## Resource requirements

| Step | RAM | GPU | Time |
|------|-----|-----|------|
| Preprocess | 256 GB | None | Minutes |
| Inference (state_sm) | 256 GB | Any GPU | 8–20 hrs |
| Export CSV | Low | None | Minutes |

## Example (fully bundled — run this to test end-to-end)

```bash
# Step 1 – Preprocess (fetal heart + 20 cardiac TFs)
python insilico_perturbation_pipeline/scripts/preprocess_for_inference.py \
    --input insilico_perturbation_pipeline/data/example_fetal_heart.h5ad \
    --gene-list insilico_perturbation_pipeline/data/example_cardiac_tfs.csv \
    --output insilico_perturbation_pipeline/data/inference_template.h5ad \
    --cell-type-col celltype \
    --n-cells 10000 \
    --seed 42

# Step 2 – Inference
python insilico_perturbation_pipeline/scripts/run_insilico_perturbation.py \
    --input insilico_perturbation_pipeline/data/inference_template.h5ad \
    --output insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --model-dir insilico_perturbation_pipeline/models \
    --checkpoint "insilico_perturbation_pipeline/models/checkpoints/step=step=18000-val_loss=val_loss=1.7692.ckpt"

# Step 3 – Export to CSV
python insilico_perturbation_pipeline/scripts/extract_to_csv.py \
    --predictions insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --output-dir insilico_perturbation_pipeline/data/perturbation_csvs
```

## Notes

- Model expects exactly 18,080 genes — preprocessing handles alignment automatically
- Input data must have cell type labels; specify the column name with `--cell-type-col`
- The bundled example (`example_fetal_heart.h5ad`) uses `cell_type` as the cell type column
- For H100-only `state_base` model, see `references/models.md`
