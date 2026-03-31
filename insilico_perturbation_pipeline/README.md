# In-Silico Perturbation Pipeline

Predict gene perturbation effects on cell state using trained State models. Supports both config-driven batch runs (all datasets x all cell types) and manual single-run mode.

## Quick Start: Config-Driven Batch Mode

Define datasets, cell types, and parameters in a YAML config, then run everything:

```bash
source activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env

# Prepare all cell types + generate SLURM jobs
python scripts/run_full_pipeline.py prepare --config configs/collab_filtered_v1.yaml

# Check what's done
python scripts/run_full_pipeline.py status --config configs/collab_filtered_v1.yaml

# Submit SLURM jobs (4 at a time, auto-skips completed)
python scripts/run_full_pipeline.py submit --config configs/collab_filtered_v1.yaml
```

### Config Format

```yaml
run_name: my_run
data_dir: /path/to/data
output_dir: /path/to/runs/my_run

defaults:
  n_cells: 500          # cells per cell type
  batch_size: 500        # genes per SLURM job
  min_log2fc: 0.0        # DEG filter
  max_padj: 1.0          # DEG filter
  cell_type_col: ident
  seed: 42

model:
  name: state_sm
  model_dir: /path/to/model
  checkpoint: checkpoints/best.ckpt

slurm:
  partition: drjieliu-h100,drjieliu-l40s,drjieliu-v100
  time: "20:00:00"
  mem: 256GB
  max_concurrent: 4
  conda_env: /path/to/envs/state_env
  state_repo: /path/to/state_with_esm

datasets:
  adult_gut:
    h5ad: adult_gut/adultgut.h5ad        # relative to data_dir
    deg_csv: adult_gut/adultgutDEG.csv   # relative to data_dir
    cell_types: ["Stem cells", "TA cells", "Enterocytes", ...]
```

See `configs/collab_filtered_v1.yaml` for a complete working example (4 datasets, 46 cell types).

### Output Structure

```
runs/{run_name}/
├── {dataset}/
│   └── {cell_type}/
│       ├── {cell_type}_subset.h5ad       ← subsetted cells
│       ├── filtered_degs.csv             ← DEG table after filtering
│       ├── gene_list_all.txt             ← all filtered genes
│       └── batches/
│           ├── gene_list_batch000.txt    ← genes for this batch
│           ├── template_batch000.h5ad    ← inference template (SLURM creates)
│           └── predictions_batch000.h5ad ← predictions (SLURM creates)
├── slurm_jobs/                           ← one .sh per batch
└── slurm_logs/
```

## Quick Start: Manual Single-Run Mode

For a single h5ad + gene list (e.g., collaborator data):

```bash
# Step 1: Preprocess
python scripts/preprocess_for_inference.py \
    --input data.h5ad \
    --gene-list genes.xlsx \
    --output inference_template.h5ad \
    --cell-type-col celltype \
    --n-cells 10000

# Step 2: Inference (SLURM or direct)
sbatch scripts/run_inference.sh
# or
python scripts/run_insilico_perturbation.py \
    --input inference_template.h5ad \
    --output perturbation_predictions.h5ad \
    --model-dir models \
    --checkpoint "models/checkpoints/step=step=18000-val_loss=val_loss=1.7692.ckpt"

# Step 3: Export
python scripts/extract_to_csv.py \
    --predictions perturbation_predictions.h5ad \
    --output-dir perturbation_csvs
```

## Scripts

| Script | Purpose |
|--------|---------|
| `run_full_pipeline.py` | Config-driven orchestrator: prepare, submit, status |
| `prepare_celltype_run.py` | Subset h5ad + filter DEGs for one cell type, batch genes |
| `preprocess_for_inference.py` | Normalize, align to 18,080 genes, create inference template |
| `run_insilico_perturbation.py` | Python wrapper for `state tx infer` |
| `run_inference.sh` | SLURM job template (single run) |
| `extract_to_csv.py` | Split prediction h5ad into per-perturbation CSVs |
| `train_state_model.py` | Train a new State model (optional) |
| `create_esm_pert_features.py` | Generate ESM2 embeddings (optional) |

## Models

| Model | GPU | Notes |
|-------|-----|-------|
| `state_sm` (bundled) | Any GPU | Faster, good for initial runs |
| `state_base` (external) | H100 only | Better predictions, slower |

See [references/models.md](references/models.md) for paths and checkpoints.

## Key Constraints

- Model expects exactly **18,080 genes** — preprocessing handles alignment
- Dense float32 matrices: memory = `n_cells x (n_perts + 1) x 18,080 x 4 bytes`
- Input must be raw counts or log-normalized (auto-detected)
- Use `source activate` (not `conda activate`) in SLURM
- `state` CLI: `pip install -e .` from State repo (not `uv tool install`)
