#!/usr/bin/env python3
"""
Generalized preprocessing for in-silico perturbation inference with State model.

Accepts gene lists from various file formats (xlsx, csv, txt) and h5ad data
from any tissue/donor. Performs normalization (normalize_total + log1p),
stratified sampling, gene alignment to the model's 18,080-gene format,
and creates the inference template.

Usage:
    python preprocess_for_inference.py \
        --input data.h5ad \
        --gene-list genes.xlsx \
        --output inference_template.h5ad \
        --cell-type-col celltype \
        --n-cells 10000
"""

import argparse
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
from tqdm import tqdm


# Default path to the 18,080 target genes from the competition gene list
DEFAULT_GENE_NAMES = "/nfs/turbo/umms-drjieliu/usr/rickyhan/AGI-guided_differentiation/insilico_perturbation_pipeline/models/gene_names.csv"


def load_gene_list(gene_list_path):
    """
    Load perturbation target gene list from various formats.

    Supported formats:
    - .xlsx: reads first column, skips header if present
    - .csv: reads first column, skips header if present
    - .txt: one gene per line
    - .tsv: reads first column, skips header if present

    Returns list of gene symbols (strings).
    """
    path = Path(gene_list_path)
    suffix = path.suffix.lower()

    if suffix == '.xlsx':
        df = pd.read_excel(path)
        # Use first column, drop NaN
        genes = df.iloc[:, 0].dropna().astype(str).tolist()
    elif suffix in ('.csv', '.tsv'):
        sep = '\t' if suffix == '.tsv' else ','
        df = pd.read_csv(path, sep=sep)
        genes = df.iloc[:, 0].dropna().astype(str).tolist()
    elif suffix == '.txt':
        with open(path) as f:
            genes = [line.strip() for line in f if line.strip()]
    else:
        raise ValueError(f"Unsupported gene list format: {suffix}. Use .xlsx, .csv, .tsv, or .txt")

    # Filter out potential header values
    # If first value looks like a header (contains 'gene', 'name', 'symbol'), skip it
    if genes and any(kw in genes[0].lower() for kw in ['gene', 'name', 'symbol', 'target']):
        print(f"  Skipping likely header row: '{genes[0]}'")
        genes = genes[1:]

    # Remove duplicates while preserving order
    seen = set()
    unique_genes = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)

    print(f"  Loaded {len(unique_genes)} unique perturbation targets from {path.name}")
    return unique_genes


def load_target_genes(gene_names_path=None):
    """Load the 18,080 target genes that the model expects."""
    path = gene_names_path or DEFAULT_GENE_NAMES
    gene_names_df = pd.read_csv(path, header=None)
    target_genes = gene_names_df[0].tolist()
    print(f"Loaded {len(target_genes)} model target genes")
    return target_genes


def stratified_sample(adata, n_cells, cell_type_col='celltype', random_state=42):
    """Stratified sampling preserving cell type proportions."""
    print(f"\nStratified sampling: {adata.n_obs} -> {n_cells} cells")
    print(f"  Cell type column: {cell_type_col}")

    cell_types = adata.obs[cell_type_col].values
    unique_types = np.unique(cell_types)
    print(f"  Cell types: {len(unique_types)}")

    # Calculate target counts
    target_counts = {}
    for ct in unique_types:
        count = np.sum(cell_types == ct)
        prop = count / len(cell_types)
        target_counts[ct] = max(1, int(n_cells * prop))

    # Adjust if total exceeds n_cells
    total_target = sum(target_counts.values())
    if total_target > n_cells:
        sorted_types = sorted(target_counts.keys(), key=lambda x: target_counts[x], reverse=True)
        excess = total_target - n_cells
        for ct in sorted_types:
            if excess <= 0:
                break
            reduce_by = min(excess, target_counts[ct] - 1)
            target_counts[ct] -= reduce_by
            excess -= reduce_by

    # Cap at available count
    for ct in unique_types:
        available = np.sum(cell_types == ct)
        target_counts[ct] = min(target_counts[ct], available)

    # Sample
    np.random.seed(random_state)
    selected_indices = []
    for ct in unique_types:
        ct_indices = np.where(cell_types == ct)[0]
        n_sample = target_counts[ct]
        if n_sample >= len(ct_indices):
            selected_indices.extend(ct_indices)
        else:
            sampled = np.random.choice(ct_indices, size=n_sample, replace=False)
            selected_indices.extend(sampled)

    selected_indices = np.array(sorted(selected_indices))
    sampled = adata[selected_indices].copy()

    print(f"  Sampled {sampled.n_obs} cells")
    for ct in np.unique(sampled.obs[cell_type_col].values):
        n = np.sum(sampled.obs[cell_type_col].values == ct)
        print(f"    {ct}: {n}")

    return sampled


def align_genes_to_target(adata, target_genes):
    """Align AnnData genes to model's target gene list. Missing genes filled with zeros."""
    print(f"\nAligning genes: {adata.n_vars} -> {len(target_genes)}")

    # Determine gene name column
    if 'features' in adata.var.columns:
        current_genes = adata.var['features'].tolist()
    elif 'gene_name' in adata.var.columns:
        current_genes = adata.var['gene_name'].tolist()
    else:
        current_genes = adata.var_names.tolist()

    common = set(current_genes) & set(target_genes)
    missing = set(target_genes) - set(current_genes)
    print(f"  Common: {len(common)}, Missing (zero-padded): {len(missing)}")

    gene_to_idx = {gene: idx for idx, gene in enumerate(current_genes)}

    X = adata.X
    if hasattr(X, 'toarray'):
        X_dense = X.toarray()
    else:
        X_dense = np.array(X)

    X_aligned = np.zeros((adata.n_obs, len(target_genes)), dtype=np.float32)
    for target_idx, gene in enumerate(tqdm(target_genes, desc="  Aligning")):
        if gene in gene_to_idx:
            X_aligned[:, target_idx] = X_dense[:, gene_to_idx[gene]]

    return X_aligned


def create_inference_template(adata, X_aligned, target_genes, perturbation_targets, cell_type_col='celltype', normalized=False):
    """Create inference template with control cells + perturbation entries.

    Uses pre-allocated matrix to avoid OOM from appending copies + vstack.
    """
    n_cells = adata.n_obs
    n_conditions = len(perturbation_targets) + 1  # +1 for control
    n_total = n_cells * n_conditions
    n_genes = len(target_genes)
    mem_gb = n_total * n_genes * 4 / 1e9
    print(f"\nCreating template: {n_cells} cells x {n_conditions} conditions = {n_total} entries")
    print(f"  Estimated matrix size: {mem_gb:.1f} GB")

    original_cell_types = adata.obs[cell_type_col].values
    original_barcodes = adata.obs_names.tolist()
    ct_strs = [str(ct) for ct in original_cell_types]

    # Pre-allocate the full matrix (avoids 2x memory from append + vstack)
    print("  Pre-allocating matrix...")
    X_final = np.empty((n_total, n_genes), dtype=np.float32)

    all_target_genes = []
    all_batch_vars = []
    all_cell_types = []
    all_original_barcodes = []
    all_original_cell_types = []

    # All conditions: control first, then perturbations
    conditions = ["non-targeting"] + list(perturbation_targets)
    for i, cond in enumerate(tqdm(conditions, desc="  Building template")):
        start = i * n_cells
        end = start + n_cells
        X_final[start:end, :] = X_aligned

        all_target_genes.extend([cond] * n_cells)
        batch_name = "batch_control" if cond == "non-targeting" else f"batch_{cond}"
        all_batch_vars.extend([batch_name] * n_cells)
        all_cell_types.extend(ct_strs)
        all_original_barcodes.extend(original_barcodes)
        all_original_cell_types.extend(ct_strs)

    obs_df = pd.DataFrame({
        'target_gene': pd.Categorical(all_target_genes),
        'batch_var': pd.Categorical(all_batch_vars),
        'cell_type': pd.Categorical(all_cell_types),
        'original_barcode': all_original_barcodes,
        'original_cell_type': all_original_cell_types
    })

    var_df = pd.DataFrame(index=target_genes)
    var_df.index.name = 'gene_names'

    obs_names = [f"{bc}_{tg}_{i}" for i, (bc, tg) in enumerate(zip(all_original_barcodes, all_target_genes))]
    obs_df.index = obs_names

    uns = {'log1p': {'base': None}} if normalized else {}
    template = ad.AnnData(X=X_final, obs=obs_df, var=var_df, uns=uns)

    n_controls = (template.obs['target_gene'] == 'non-targeting').sum()
    print(f"  Template shape: {template.shape}")
    print(f"  Control cells: {n_controls}, Perturbation cells: {template.n_obs - n_controls}")

    if n_controls == 0:
        raise ValueError("No control cells created! Inference will fail.")

    return template


def main():
    parser = argparse.ArgumentParser(description="Preprocess data for in-silico perturbation inference")
    parser.add_argument("--input", "-i", required=True, help="Input h5ad file")
    parser.add_argument("--gene-list", "-g", required=True, help="Perturbation target gene list (.xlsx/.csv/.tsv/.txt)")
    parser.add_argument("--output", "-o", required=True, help="Output inference template h5ad")
    parser.add_argument("--gene-names", default=None, help=f"Model target genes CSV (default: {DEFAULT_GENE_NAMES})")
    parser.add_argument("--n-cells", "-n", type=int, default=10000, help="Number of cells to sample (default: 10000)")
    parser.add_argument("--cell-type-col", default="celltype", help="Cell type column in obs (default: celltype)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    parser.add_argument("--skip-normalize", action="store_true",
                        help="Skip normalization (use if data is already log-normalized)")

    args = parser.parse_args()

    print("=" * 70)
    print("IN-SILICO PERTURBATION PREPROCESSING")
    print("=" * 70)

    # Load perturbation targets
    perturbation_targets = load_gene_list(args.gene_list)

    # Load model target genes
    target_genes = load_target_genes(args.gene_names)

    # Validate perturbation targets exist in model genes
    missing_perts = [g for g in perturbation_targets if g not in set(target_genes)]
    if missing_perts:
        print(f"\n  WARNING: {len(missing_perts)} perturbation targets not in model gene set:")
        for g in missing_perts:
            print(f"    - {g}")
        print("  These genes will still be included but may have limited prediction quality.")

    # Load data
    print(f"\nLoading {args.input}...")
    adata = ad.read_h5ad(args.input)
    print(f"  Shape: {adata.shape}")
    print(f"  Obs columns: {adata.obs.columns.tolist()}")

    # Normalize: normalize_total + log1p (matching State model training)
    normalized = False
    if args.skip_normalize:
        print("\n  Skipping normalization (--skip-normalize flag set)")
    else:
        # Coerce to CSR if sparse (normalize_total expects CSR, not CSC)
        import scipy.sparse as sp
        if sp.issparse(adata.X) and not sp.isspmatrix_csr(adata.X):
            print(f"\n  Converting sparse matrix from {type(adata.X).__name__} to CSR")
            adata.X = adata.X.tocsr()

        X = adata.X
        if hasattr(X, 'toarray'):
            sample = X[:min(500, adata.n_obs)].toarray()
        else:
            sample = np.array(X[:min(500, adata.n_obs)])

        has_negative = sample.min() < 0
        nz = sample[sample != 0]
        is_integer = np.allclose(nz[:1000], np.round(nz[:1000])) if len(nz) > 0 else True
        max_val = sample.max()
        frac_zero = (sample == 0).mean()

        # Detection logic:
        # - Negative values -> scaled/SCT (cannot normalize, warn)
        # - Integer values + high sparsity -> raw counts (normalize)
        # - Float values + low max + no negatives -> likely log-normalized (skip)
        # - Ambiguous -> warn and ask user to use --skip-normalize if needed
        if has_negative:
            print(f"\n  WARNING: Data has negative values (min={sample.min():.3f}).")
            print("  This looks like scaled/SCT data, not raw counts.")
            print("  Normalization will be skipped, but predictions may be affected.")
            print("  Consider providing raw counts instead.")
        elif is_integer and frac_zero > 0.5:
            print(f"\n  Detected raw counts (integer values, max={max_val:.0f}, {frac_zero*100:.0f}% zeros)")
            print("  Applying: sc.pp.normalize_total() + sc.pp.log1p()")
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            normalized = True
            print("  Normalization complete.")
        elif not is_integer and max_val < 20:
            print(f"\n  Data appears already log-normalized (float values, max={max_val:.2f})")
            print("  Skipping normalization. Use --skip-normalize to suppress this check.")
        else:
            print(f"\n  WARNING: Uncertain data format (max={max_val:.2f}, integer={is_integer}, zeros={frac_zero*100:.0f}%)")
            print("  Cannot safely auto-detect. Please specify:")
            print("    --skip-normalize  if data is already log-normalized")
            print("  Without flag, assuming raw counts and normalizing.")
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            normalized = True
            print("  Normalization applied.")

    # Stratified sampling
    sampled = stratified_sample(adata, args.n_cells, args.cell_type_col, args.seed)

    # Align genes
    X_aligned = align_genes_to_target(sampled, target_genes)

    # Create template
    template = create_inference_template(sampled, X_aligned, target_genes, perturbation_targets, args.cell_type_col, normalized=normalized)

    # Save
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nSaving to {args.output}...")
    template.write_h5ad(args.output)

    # Save metadata
    metadata_path = output_path.parent / f"{output_path.stem}_metadata.csv"
    cell_type_info = sampled.obs[args.cell_type_col].value_counts().reset_index()
    cell_type_info.columns = ['cell_type_id', 'n_cells']
    cell_type_info.to_csv(metadata_path, index=False)

    print(f"\nDone! Template: {template.shape}")
    print(f"  Metadata: {metadata_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
