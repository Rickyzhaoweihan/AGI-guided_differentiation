#!/usr/bin/env python3
"""Prepare a per-cell-type perturbation run from h5ad + DEG CSV.

Given an h5ad file and a Seurat FindAllMarkers DEG CSV, this script:
1. Subsets the h5ad to cells of a specified cell type
2. Filters DEGs for that cell type
3. Applies optional filters (p_val_adj, log2FC, top-N)
4. Exports a flat gene list compatible with the preprocessing script
5. Saves the subset h5ad

Usage:
    python prepare_celltype_run.py \
        --h5ad /path/to/data.h5ad \
        --deg-csv /path/to/DEGs.csv \
        --cell-type "Stem cells" \
        --cell-type-col ident \
        --output-dir /path/to/output/stem_cells/ \
        --top-n 50 \
        --min-log2fc 0.5 \
        --max-padj 0.05
"""

import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare per-cell-type perturbation run")
    parser.add_argument("--h5ad", required=True, help="Input h5ad file")
    parser.add_argument("--deg-csv", required=True, help="Seurat FindAllMarkers DEG CSV")
    parser.add_argument("--cell-type", required=True, help="Cell type to extract (must match DEG cluster name)")
    parser.add_argument("--cell-type-col", default="ident", help="Column in h5ad .obs for cell type (default: ident)")
    parser.add_argument("--output-dir", required=True, help="Output directory for subset h5ad and gene list")
    parser.add_argument("--top-n", type=int, default=None, help="Keep only top N genes by abs(log2FC)")
    parser.add_argument("--min-log2fc", type=float, default=1.0, help="Minimum abs(avg_log2FC) (default: 1.0)")
    parser.add_argument("--max-padj", type=float, default=0.05, help="Maximum adjusted p-value (default: 0.05)")
    parser.add_argument("--upregulated-only", action="store_true", help="Only keep upregulated genes (positive log2FC)")
    parser.add_argument("--batch-size", type=int, default=500, help="Max genes per batch (default: 500)")
    parser.add_argument("--gene-names-csv", default=None,
                        help="Model gene names CSV for validation. Default: competition_support_set path")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load DEG CSV ---
    print(f"Loading DEGs from {args.deg_csv}")
    deg = pd.read_csv(args.deg_csv)
    assert "cluster" in deg.columns and "gene" in deg.columns, \
        f"DEG CSV must have 'cluster' and 'gene' columns. Found: {deg.columns.tolist()}"

    available_clusters = sorted(deg["cluster"].unique())
    if args.cell_type not in available_clusters:
        print(f"ERROR: Cell type '{args.cell_type}' not found in DEG clusters.")
        print(f"Available: {available_clusters}")
        sys.exit(1)

    # Filter to this cell type
    ct_deg = deg[deg["cluster"] == args.cell_type].copy()
    print(f"  {len(ct_deg)} DEGs for '{args.cell_type}'")

    # Apply filters
    if "p_val_adj" in ct_deg.columns:
        ct_deg = ct_deg[ct_deg["p_val_adj"] <= args.max_padj]
        print(f"  {len(ct_deg)} after p_val_adj <= {args.max_padj}")

    if "avg_log2FC" in ct_deg.columns:
        if args.upregulated_only:
            ct_deg = ct_deg[ct_deg["avg_log2FC"] >= args.min_log2fc]
            print(f"  {len(ct_deg)} after avg_log2FC >= {args.min_log2fc} (upregulated only)")
        else:
            ct_deg = ct_deg[ct_deg["avg_log2FC"].abs() >= args.min_log2fc]
            print(f"  {len(ct_deg)} after |avg_log2FC| >= {args.min_log2fc}")

    # Sort by absolute fold change
    if "avg_log2FC" in ct_deg.columns:
        ct_deg = ct_deg.sort_values("avg_log2FC", key=abs, ascending=False)

    if args.top_n and len(ct_deg) > args.top_n:
        ct_deg = ct_deg.head(args.top_n)
        print(f"  Kept top {args.top_n} by |avg_log2FC|")

    genes = ct_deg["gene"].unique().tolist()
    print(f"  Final: {len(genes)} unique genes")

    if len(genes) == 0:
        print("ERROR: No genes passed filters. Adjust --min-log2fc, --max-padj, or --top-n.")
        sys.exit(1)

    # --- Validate against model gene set ---
    gene_names_path = args.gene_names_csv
    if gene_names_path is None:
        # Try bundled gene_names.csv first, then NFS fallback
        bundled = os.path.join(os.path.dirname(__file__), "..", "models", "gene_names.csv")
        nfs = "/nfs/turbo/umms-drjieliu/usr/zheyuz/state_with_esm/competition_support_set/competition_support_set/gene_names.csv"
        if os.path.exists(bundled):
            gene_names_path = bundled
        elif os.path.exists(nfs):
            gene_names_path = nfs

    if gene_names_path and os.path.exists(gene_names_path):
        model_genes = set(pd.read_csv(gene_names_path, header=None)[0].tolist())
        in_model = [g for g in genes if g in model_genes]
        not_in_model = [g for g in genes if g not in model_genes]
        print(f"  {len(in_model)}/{len(genes)} genes in model's 18,080 gene set")
        if not_in_model:
            print(f"  Not in model (will be zero-padded in expression): {not_in_model[:10]}{'...' if len(not_in_model) > 10 else ''}")

    # --- Save gene lists (batched) ---
    # Save full filtered DEG table for reference
    deg_path = os.path.join(args.output_dir, "filtered_degs.csv")
    ct_deg.to_csv(deg_path, index=False)
    print(f"  Filtered DEG table saved: {deg_path}")

    # Split genes into batches
    batch_size = args.batch_size
    n_batches = (len(genes) + batch_size - 1) // batch_size
    batches_dir = os.path.join(args.output_dir, "batches")
    os.makedirs(batches_dir, exist_ok=True)

    gene_list_paths = []
    for i in range(n_batches):
        batch_genes = genes[i * batch_size : (i + 1) * batch_size]
        batch_path = os.path.join(batches_dir, f"gene_list_batch{i:03d}.txt")
        with open(batch_path, "w") as f:
            for g in batch_genes:
                f.write(g + "\n")
        gene_list_paths.append(batch_path)

    # Also save full gene list for reference
    full_list_path = os.path.join(args.output_dir, "gene_list_all.txt")
    with open(full_list_path, "w") as f:
        for g in genes:
            f.write(g + "\n")

    print(f"  {len(genes)} genes split into {n_batches} batches of <={batch_size}")
    print(f"  Batch gene lists saved in: {batches_dir}/")

    # --- Subset h5ad ---
    print(f"\nLoading h5ad from {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad)
    print(f"  Full dataset: {adata.shape}")

    if args.cell_type_col not in adata.obs.columns:
        print(f"ERROR: Column '{args.cell_type_col}' not in .obs. Available: {adata.obs.columns.tolist()}")
        sys.exit(1)

    # Find matching cell type (handle case differences)
    obs_types = adata.obs[args.cell_type_col].unique()
    exact_match = args.cell_type in obs_types
    if not exact_match:
        # Try case-insensitive match
        lower_map = {t.lower(): t for t in obs_types}
        if args.cell_type.lower() in lower_map:
            matched = lower_map[args.cell_type.lower()]
            print(f"  Case-insensitive match: '{args.cell_type}' -> '{matched}'")
            mask = adata.obs[args.cell_type_col] == matched
        else:
            print(f"ERROR: Cell type '{args.cell_type}' not found in .obs['{args.cell_type_col}'].")
            print(f"Available: {sorted(obs_types)}")
            sys.exit(1)
    else:
        mask = adata.obs[args.cell_type_col] == args.cell_type

    subset = adata[mask].copy()
    print(f"  Subset '{args.cell_type}': {subset.shape[0]} cells")

    # Save subset
    safe_name = args.cell_type.replace(" ", "_").replace("/", "_").lower()
    subset_path = os.path.join(args.output_dir, f"{safe_name}_subset.h5ad")
    subset.write_h5ad(subset_path)
    print(f"  Subset h5ad saved: {subset_path}")

    # --- Summary ---
    n_cells = min(subset.shape[0], 1000)
    preprocess_script = os.path.join(os.path.dirname(__file__), "preprocess_for_inference.py")
    print(f"\n{'='*60}")
    print(f"Ready for preprocessing. {n_batches} batch(es) to run:")
    print(f"  Cell type subset: {subset_path}")
    print(f"  Cells: {n_cells} (from {subset.shape[0]} available)")
    print(f"  Genes: {len(genes)} total in {n_batches} batches of <={batch_size}")
    for i, bp in enumerate(gene_list_paths):
        batch_n = len(open(bp).readlines())
        template_out = os.path.join(args.output_dir, f"batches/template_batch{i:03d}.h5ad")
        pred_out = os.path.join(args.output_dir, f"batches/predictions_batch{i:03d}.h5ad")
        mem_gb = n_cells * (batch_n + 1) * 18080 * 4 / 1e9
        print(f"\n  --- Batch {i} ({batch_n} genes, ~{mem_gb:.1f} GB template) ---")
        print(f"  python {preprocess_script} \\")
        print(f"      --input {subset_path} \\")
        print(f"      --gene-list {bp} \\")
        print(f"      --output {template_out} \\")
        print(f"      --cell-type-col {args.cell_type_col} \\")
        print(f"      --n-cells {n_cells} \\")
        print(f"      --seed {args.seed}")
    print(f"\n{'='*60}")


if __name__ == "__main__":
    main()
