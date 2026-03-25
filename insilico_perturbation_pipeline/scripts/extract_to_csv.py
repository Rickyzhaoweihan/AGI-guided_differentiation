#!/usr/bin/env python3
"""
Extract perturbation predictions to individual CSV files.

Creates one CSV per perturbation condition + one control baseline CSV.
Each CSV has cells as rows, genes as columns, with optional metadata.

Usage:
    python extract_to_csv.py \
        --predictions perturbation_predictions.h5ad \
        --output-dir perturbation_csv_files
"""

import argparse
import anndata as ad
import pandas as pd
from pathlib import Path
from tqdm import tqdm


def extract_perturbations_to_csv(predictions_file, output_dir, include_metadata=True):
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    print(f"Loading {predictions_file}...")
    pred = ad.read_h5ad(predictions_file)
    print(f"  {pred.n_obs:,} predictions, {pred.n_vars:,} genes")

    gene_names = pred.var_names.tolist()
    perturbations = sorted(pred.obs['target_gene'].unique())
    print(f"  {len(perturbations)} conditions")

    for pert in tqdm(perturbations, desc="Extracting"):
        pert_data = pred[pred.obs['target_gene'] == pert].copy()

        if hasattr(pert_data.X, 'toarray'):
            expr = pert_data.X.toarray()
        else:
            expr = pert_data.X

        df = pd.DataFrame(expr, columns=gene_names, index=pert_data.obs_names)

        if include_metadata:
            df.insert(0, 'cell_type', pert_data.obs['cell_type'].values)
            df.insert(1, 'original_barcode', pert_data.obs['original_barcode'].values)
            if 'original_cell_type' in pert_data.obs.columns:
                df.insert(2, 'original_cell_type', pert_data.obs['original_cell_type'].values)

        filename = 'control_baseline.csv' if pert == 'non-targeting' else f'{pert}_perturbation.csv'
        df.to_csv(output_path / filename, index_label='cell_id')

    print(f"Done! {len(perturbations)} files in {output_dir}/")


def main():
    parser = argparse.ArgumentParser(description="Extract perturbation predictions to CSV files")
    parser.add_argument('--predictions', '-p', required=True, help='Predictions h5ad file')
    parser.add_argument('--output-dir', '-o', default='perturbation_csv_files', help='Output directory')
    parser.add_argument('--no-metadata', action='store_true', help='Exclude metadata columns')
    args = parser.parse_args()

    extract_perturbations_to_csv(args.predictions, args.output_dir, not args.no_metadata)


if __name__ == '__main__':
    main()
