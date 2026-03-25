# In-Silico Perturbation Pipeline (Chen Lab Collaboration)

Predict gene perturbation effects on cell state using trained State models.
Data comes from Shuibing Chen's lab at Cornell — fetal tissue scRNA-seq from various organs/donors.

## Workflow

### Step 1: Gather Inputs

Required inputs:
1. **h5ad file** — scRNA-seq data (any tissue/donor)
2. **Gene list** — perturbation targets (.xlsx, .csv, .tsv, or .txt)
3. **Cell type column** — name in `.obs` (inspect data to find it)
4. **Model** — which State model to use. See [references/models.md](references/models.md)
5. **Output directory** — where to save all outputs

### Step 2: Explore Input Data

Before preprocessing, inspect the h5ad:

```python
import anndata as ad
adata = ad.read_h5ad("path/to/data.h5ad")
print(f"Shape: {adata.shape}")
print(f"Obs columns: {adata.obs.columns.tolist()}")
print(f"Var columns: {adata.var.columns.tolist()}")
print(f"Cell types: {adata.obs['<cell_type_col>'].value_counts()}")
```

Check:
- What column has cell type labels? (common: `celltype`, `cell_type`, `CellType`)
- What column has gene names in `.var`? (common: `var_names`, `features`, `gene_name`)
- **Data format**: Must be raw counts or log-normalized. Check for: integer values + high max (raw counts), negative values (scaled — BAD), float values with max ~10 (log-normalized). The model was trained on log-normalized data (`normalize_total` + `log1p`).
- How many genes? Need full transcriptome (not just HVGs) for best results — model uses 18,080 genes
- How many cells? If >15k, subsampling is recommended

Also inspect the gene list file to confirm format and gene count.

### Step 3: Preprocess

Run the preprocessing script:

```bash
python scripts/preprocess_for_inference.py \
    --input /path/to/data.h5ad \
    --gene-list /path/to/gene_list.xlsx \
    --output /path/to/output_dir/inference_template.h5ad \
    --cell-type-col celltype \
    --n-cells 10000 \
    --seed 42
```

The script handles:
- **Normalization**: Auto-detects raw counts and applies `normalize_total` + `log1p` (matching model training). Warns on scaled/SCT data. Use `--skip-normalize` if data is already log-normalized.
- Stratified sampling preserving cell type proportions
- Gene alignment to the model's 18,080-gene format (missing genes zero-padded)
- Creating control cells (non-targeting) + perturbation entries
- Validation that perturbation targets exist in model gene set

**Output**: inference template h5ad with shape `(n_cells * (n_perturbations + 1), 18080)`

### Step 4: Validate Gene List

Before running inference, validate the gene list:
- Check overlap with ESM2 embeddings (19,790 genes) and model perturbation map
- Check overlap with 18,080 expression gene set
- Flag any mouse gene symbols from JASPAR — map to human orthologs (e.g., RHOX11 → RHOXF1)
- Genes not in ESM2 will fall back to default encoding (less meaningful predictions)

### Step 5: Run Inference via SLURM

Edit `scripts/run_inference.sh` with your paths and submit:

```bash
sbatch scripts/run_inference.sh
```

Key `state tx infer` parameters:
- `--adata` — the inference template from Step 3
- `--output` — where to save predictions
- `--model-dir` — model training directory
- `--checkpoint` — specific checkpoint file
- `--pert-col target_gene` — always this (set by preprocessing)
- `--celltype-col cell_type` — always this (set by preprocessing)
- `--batch-col batch_var` — always this (set by preprocessing)

SLURM settings (for ~100 perturbations x 10k cells):
- **Preprocessing**: `--mem=256GB` (dense template is n_cells x n_conditions x 18080 x 4 bytes)
- **state_sm inference**: any GPU partition, 8-20h, `--mem=256GB`
- **state_base inference**: H100 only (`--partition=drjieliu-h100`), 20h, `--mem=256GB`
- Always use `source activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env`
- `state` CLI must be installed via `pip install -e .` from repo root (NOT `uv tool install`)

### Step 6: Split Predictions

The full prediction h5ad (~69 GB for 10k cells x 100 perts) is too large for collaborators. Split into per-perturbation files:

```python
import anndata as ad, os
pred = ad.read_h5ad('perturbation_predictions.h5ad')
os.makedirs('predictions_<model>/', exist_ok=True)

# Control baseline
ctrl = pred[pred.obs['target_gene'] == 'non-targeting'].copy()
ctrl.write_h5ad('predictions_<model>/control_baseline.h5ad')

# Each perturbation
for p in pred.obs['target_gene'].unique():
    if p != 'non-targeting':
        pred[pred.obs['target_gene'] == p].copy().write_h5ad(f'predictions_<model>/{p}.h5ad')
```

Result: ~726 MB per file, loadable on any laptop.

Alternatively, use `scripts/extract_to_csv.py` to export as CSV files:

```bash
python scripts/extract_to_csv.py \
    --predictions perturbation_predictions.h5ad \
    --output-dir perturbation_csv_files
```

### Step 7: Create Analysis Notebooks & Tutorials

Provide collaborators with:
1. **Jupyter notebooks** (`analyze_predictions_<model>.ipynb`) — load from split h5ad files, compute DE on the fly
2. **R Markdown tutorials** (`analyze_predictions_tutorial_<model>.Rmd`) — same analyses in R using `anndata` R package

Both should work from split files (~2-4 GB RAM peak). See Jiajun and Jeya examples for templates.

## Critical Details

- **Normalization**: The model was trained on log-normalized data (`sc.pp.normalize_total` + `sc.pp.log1p`). The preprocessing script auto-detects raw counts and normalizes. The `state tx infer` pipeline does NOT normalize internally — data must be pre-normalized. Ask collaborators for **raw counts**; we handle normalization.
- **Perturbation encoding**: Both models use **ESM2 protein embeddings** (not one-hot) for perturbation encoding. Config: `pert_rep: onehot` with `perturbation_features_file: ESM2_pert_features.pt`. Coverage: ~19,790 human genes.
- **Control cells are required**: Template MUST include `non-targeting` control cells. The preprocessing script handles this automatically.
- **Gene alignment**: Model expects exactly 18,080 genes in a specific order. Preprocessing aligns to `competition_support_set/gene_names.csv`.
- **Full transcriptome needed**: HVG-only data (e.g., 2,000 genes) -> ~16,000 zero-padded genes. Ask collaborators for all genes before HVG filtering.
- **Obs column names after preprocessing**: Always `target_gene`, `cell_type`, `batch_var` regardless of original column names.
- **Dense matrix**: The inference template uses dense float32 matrices (not sparse). Memory: `n_cells x (n_perts + 1) x 18080 x 4 bytes`.
- **JASPAR gene lists**: May contain mouse gene symbols. Check and map to human orthologs before running.
- **`state` CLI**: Must be installed in the conda env via `pip install -e .` from the state repo root. The old `uv tool install` symlink is broken.

## Known Issues

- **Jiajun's prior runs used raw counts without normalization** — predictions may be suboptimal. Consider re-running with the updated preprocessing script.
- **`sc.pp.normalize_total()` with default `target_sum=None`** normalizes to median library size, not 10,000. This matches training if training also used the default.

## Prior Runs

| Dataset | Tissue | Perturbations | Models Used | Output Dir |
|---------|--------|---------------|-------------|------------|
| Jiajun | Fetal heart | 20 cardiac TFs | state_sm, state_base | `jiajun_data/` |
| Jeya | Fetal gut | 101 gut TFs (JASPAR) | state_sm, state_base | `jeya_data_and_pred/` |

### Jeya Run Details (2026-03)
- Input: `fetalgut_latest.h5ad` (71,674 cells x 36,601 genes, raw counts)
- Gene list: `Genelistgut_updated.xlsx` (101 TFs, RHOX11->RHOXF1 mapped)
- Cell types: Enterocytes (56k), TA (7.7k), Stem Cells (4.4k), EECs (1.7k), Goblet Cells (1.7k)
- Subsampled to ~10k cells, normalized with `normalize_total` + `log1p`
- 3 genes not in expression output: BHLHA15, ESR2, SOX11
- Predictions split into `predictions_sm/` and `predictions_base/` (102 files each, ~726 MB per file)
- Analysis notebooks: `analyze_predictions_sm.ipynb`, `analyze_predictions_base.ipynb`
- R tutorials: `analyze_predictions_tutorial_sm.Rmd`, `analyze_predictions_tutorial_base.Rmd`
