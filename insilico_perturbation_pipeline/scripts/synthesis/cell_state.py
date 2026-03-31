"""Phase 5: Phenotypic interpretation via cell state scoring.

Scores curated marker gene sets in perturbed vs. control cells.
Pro-apoptotic and anti-apoptotic are scored separately (unsigned score_genes
cannot handle opposing genes in one set).
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


# Built-in marker sets — well-established, tissue-agnostic
BUILTIN_MARKERS = {
    "proliferation": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM6", "CDK1", "CCNB1", "CCNA2"],
    "pro_apoptotic": ["BAX", "BAK1", "BID", "CASP3", "CASP9", "CYCS", "TP53", "PMAIP1"],
    "anti_apoptotic": ["BCL2", "BCL2L1", "MCL1", "BIRC5", "XIAP"],
    "stress_upr": ["HSPA1A", "HSPA1B", "HSP90AA1", "DDIT3", "ATF4", "XBP1"],
}


def run_cell_state_scoring(X_ctrl, X_pert, gene_names, custom_gene_sets=None):
    """Score cell states in control vs. perturbed cells.

    Args:
        X_ctrl, X_pert: np.ndarray (n_cells, n_genes) matched by row
        gene_names: list of gene names
        custom_gene_sets: dict of {name: [gene_list]} from expert YAML

    Returns:
        pd.DataFrame with state, ctrl_score, pert_score, delta, interpretation
    """
    # Combine all gene sets
    all_sets = dict(BUILTIN_MARKERS)
    if custom_gene_sets:
        all_sets.update(custom_gene_sets)

    # Filter to genes present in the data
    gene_set = set(gene_names)
    filtered_sets = {}
    for name, genes in all_sets.items():
        present = [g for g in genes if g in gene_set]
        if len(present) >= 3:  # Need at least 3 genes for meaningful scoring
            filtered_sets[name] = present

    if not filtered_sets:
        return pd.DataFrame(columns=["state", "ctrl_score", "pert_score", "delta", "interpretation"])

    # Create temporary AnnData for scoring
    X_combined = np.vstack([X_ctrl, X_pert])
    n_ctrl = X_ctrl.shape[0]
    labels = ["ctrl"] * n_ctrl + ["pert"] * X_pert.shape[0]

    adata = ad.AnnData(
        X=X_combined,
        var=pd.DataFrame(index=gene_names),
        obs=pd.DataFrame({"group": labels}),
    )

    # Score each gene set
    results = []
    for name, genes in filtered_sets.items():
        score_key = f"score_{name}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=score_key, use_raw=False)

        ctrl_scores = adata.obs.loc[adata.obs["group"] == "ctrl", score_key].values
        pert_scores = adata.obs.loc[adata.obs["group"] == "pert", score_key].values

        ctrl_mean = ctrl_scores.mean()
        pert_mean = pert_scores.mean()
        delta = pert_mean - ctrl_mean

        interpretation = _interpret_delta(name, delta)

        results.append({
            "state": name,
            "ctrl_score": round(ctrl_mean, 4),
            "pert_score": round(pert_mean, 4),
            "delta": round(delta, 4),
            "n_genes_scored": len(genes),
            "interpretation": interpretation,
        })

    return pd.DataFrame(results)


def _interpret_delta(state_name, delta, threshold=0.05):
    """Generate human-readable interpretation of score delta."""
    if abs(delta) < threshold:
        return "no change"

    direction = "increased" if delta > 0 else "decreased"
    magnitude = "slightly" if abs(delta) < 0.1 else ("moderately" if abs(delta) < 0.3 else "strongly")

    interpretations = {
        "proliferation": f"Proliferation {magnitude} {direction}",
        "pro_apoptotic": f"Pro-apoptotic program {magnitude} {direction}",
        "anti_apoptotic": f"Anti-apoptotic (survival) program {magnitude} {direction}",
        "stress_upr": f"Stress/UPR response {magnitude} {direction}",
    }

    return interpretations.get(state_name, f"{state_name} {magnitude} {direction}")
