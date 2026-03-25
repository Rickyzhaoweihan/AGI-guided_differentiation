#!/usr/bin/env python
"""
Generate ESM2 perturbation features file for genetic perturbations.

This script creates a dictionary mapping gene names to their ESM2 embeddings,
which can be used with the STATE data module for training.

The ESM2 embeddings are 5,120-dimensional vectors from the ESM2-3B protein
language model, providing rich representations of protein function.

Usage:
    # From gene list file
    python create_esm_pert_features.py \
        --gene-list genes.txt \
        --output ESM2_pert_features.pt

    # From TOML config (extracts genes from datasets)
    python create_esm_pert_features.py \
        --from-toml competition.toml \
        --output ESM2_pert_features.pt

Output format:
    Dictionary saved as PyTorch file (.pt):
    {
        "GENE1": torch.Tensor([5120]),
        "GENE2": torch.Tensor([5120]),
        "non-targeting": torch.zeros([5120]),  # Control
        ...
    }
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, Set

import torch
import toml
import h5py
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_esm_embeddings(esm_path: str) -> Dict[str, torch.Tensor]:
    """
    Load ESM2 embeddings from file.

    Args:
        esm_path: Path to ESM embeddings file (dict format)

    Returns:
        Dictionary mapping gene names to ESM embeddings
    """
    logger.info(f"Loading ESM embeddings from {esm_path}")
    esm_data = torch.load(esm_path, weights_only=False)

    if isinstance(esm_data, dict):
        logger.info(f"Loaded {len(esm_data)} gene embeddings")
        first_emb = next(iter(esm_data.values()))
        logger.info(f"Embedding dimension: {first_emb.shape}")
        return esm_data
    else:
        raise ValueError(f"Expected dict format, got {type(esm_data)}")


def extract_genes_from_toml(toml_path: str) -> Set[str]:
    """
    Extract all unique gene names from datasets specified in TOML config.

    Args:
        toml_path: Path to TOML configuration file

    Returns:
        Set of unique gene names
    """
    logger.info(f"Reading dataset configuration from {toml_path}")
    config = toml.load(toml_path)

    genes = set()

    if 'datasets' not in config:
        raise ValueError("TOML file must contain [datasets] section")

    dataset_paths = config['datasets']
    logger.info(f"Found {len(dataset_paths)} datasets")

    for dataset_name, dataset_path in dataset_paths.items():
        logger.info(f"Processing dataset: {dataset_name} at {dataset_path}")

        dataset_dir = Path(dataset_path)
        if not dataset_dir.exists():
            logger.warning(f"Dataset directory does not exist: {dataset_path}")
            continue

        h5_files = list(dataset_dir.glob("*.h5")) + list(dataset_dir.glob("*.h5ad"))
        logger.info(f"  Found {len(h5_files)} h5/h5ad files")

        for h5_file in h5_files:
            logger.info(f"  Reading: {h5_file.name}")
            try:
                with h5py.File(h5_file, 'r') as f:
                    pert_col_options = [
                        'obs/target_gene', 'obs/gene',
                        'obs/pert', 'obs/perturbation'
                    ]

                    for pert_col in pert_col_options:
                        if pert_col in f:
                            pert_data = f[pert_col]

                            if f'{pert_col}/categories' in f:
                                categories = f[f'{pert_col}/categories'][:].astype(str)
                                file_genes = set(categories)
                            else:
                                file_genes = set(pert_data[:].astype(str))

                            genes.update(file_genes)
                            logger.info(f"    Found {len(file_genes)} genes from {pert_col}")
                            break
                    else:
                        logger.warning(f"    No perturbation column found in {h5_file.name}")

            except Exception as e:
                logger.error(f"    Error reading {h5_file.name}: {e}")

    # Remove control perturbations
    control_perts = {'non-targeting', 'DMSO', 'DMSO_TF', 'ctrl', 'control', 'NT', 'neg'}
    genes = genes - control_perts

    logger.info(f"Total unique genes found: {len(genes)}")
    return genes


def extract_genes_from_file(gene_file: str) -> Set[str]:
    """
    Extract gene names from a text file (one gene per line).

    Args:
        gene_file: Path to text file with gene names

    Returns:
        Set of gene names
    """
    logger.info(f"Reading genes from {gene_file}")
    genes = set()

    with open(gene_file, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.add(gene)

    logger.info(f"Loaded {len(genes)} genes from file")
    return genes


def create_pert_features(
    genes: Set[str],
    esm_embeddings: Dict[str, torch.Tensor],
    include_control: bool = True
) -> Dict[str, torch.Tensor]:
    """
    Create perturbation features dictionary from gene list and ESM embeddings.

    Args:
        genes: Set of gene names to include
        esm_embeddings: Dictionary of ESM embeddings
        include_control: If True, add zero embedding for control perturbations

    Returns:
        Dictionary mapping gene names to ESM embeddings
    """
    pert_features = {}
    missing_genes = []

    first_emb = next(iter(esm_embeddings.values()))
    emb_dim = first_emb.shape[0]

    # Add control perturbations with zero embedding
    if include_control:
        control_names = ['non-targeting', 'DMSO_TF', 'DMSO', 'ctrl']
        for control in control_names:
            pert_features[control] = torch.zeros(emb_dim)
        logger.info(f"Added {len(control_names)} control perturbations with zero embeddings")

    # Map genes to embeddings
    for gene in tqdm(genes, desc="Creating perturbation features"):
        if gene in esm_embeddings:
            pert_features[gene] = esm_embeddings[gene]
        else:
            missing_genes.append(gene)
            pert_features[gene] = torch.zeros(emb_dim)

    logger.info(f"Created perturbation features for {len(pert_features)} genes/perturbations")

    if missing_genes:
        logger.warning(f"Missing ESM embeddings for {len(missing_genes)} genes:")
        logger.warning(f"  {missing_genes[:10]}{'...' if len(missing_genes) > 10 else ''}")

    return pert_features


def main():
    parser = argparse.ArgumentParser(
        description="Create ESM2 perturbation features file",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--gene-list',
        type=str,
        help='Text file with gene names (one per line)'
    )
    input_group.add_argument(
        '--from-toml',
        type=str,
        help='TOML config file - extract genes from datasets'
    )

    # ESM embeddings path
    parser.add_argument(
        '--esm-embeddings',
        type=str,
        default='/large_storage/ctc/ML/data/cell/misc/Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM2.pt',
        help='Path to ESM2 embeddings file'
    )

    # Output
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output path for perturbation features file (.pt)'
    )

    parser.add_argument(
        '--no-control',
        action='store_true',
        help='Do not add control perturbation with zero embedding'
    )

    args = parser.parse_args()

    # Load ESM embeddings
    try:
        esm_embeddings = load_esm_embeddings(args.esm_embeddings)
    except Exception as e:
        logger.error(f"Failed to load ESM embeddings: {e}")
        raise SystemExit(1)

    # Get gene list
    if args.gene_list:
        genes = extract_genes_from_file(args.gene_list)
    else:
        genes = extract_genes_from_toml(args.from_toml)

    if not genes:
        logger.error("No genes found!")
        raise SystemExit(1)

    # Create perturbation features
    pert_features = create_pert_features(
        genes,
        esm_embeddings,
        include_control=not args.no_control
    )

    # Save
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Saving perturbation features to {output_path}")
    torch.save(pert_features, output_path)

    logger.info("Successfully created perturbation features file!")
    logger.info(f"  - Genes: {len(pert_features)}")
    logger.info(f"  - Embedding dim: {next(iter(pert_features.values())).shape[0]}")
    logger.info(f"\nUse with STATE training:")
    logger.info(f'  data.kwargs.perturbation_features_file="{output_path}"')


if __name__ == '__main__':
    main()
