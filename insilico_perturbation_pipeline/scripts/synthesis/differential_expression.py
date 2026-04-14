"""Phase 2: Differential expression with paired per-cell deltas.

Exploits the paired structure of the prediction template: each original cell
appears once under non-targeting control and once per perturbation, linked by
original_barcode. Computes per-cell deltas, true log2FC from natural-scale
means, and FDR-corrected significance via paired t-test.
"""

import numpy as np
import pandas as pd
from scipy import stats


def compute_de(X_ctrl, X_pert, gene_names, n_bootstrap=0):
    """Compute differential expression with paired structure.

    Args:
        X_ctrl: np.ndarray (n_cells, n_genes) — control expression (log1p-normalized)
        X_pert: np.ndarray (n_cells, n_genes) — perturbed expression (matched by row)
        gene_names: list of gene names
        n_bootstrap: number of bootstrap iterations for CI (0 = skip, use t-test CI)

    Returns:
        pd.DataFrame with columns: gene, log2fc, log2fc_ci_low, log2fc_ci_high,
                                   pval, fdr, mean_pert, mean_ctrl, is_significant
    """
    n_cells, n_genes = X_ctrl.shape

    # Natural-scale means for true log2FC
    nat_ctrl = np.expm1(X_ctrl)  # reverse log1p
    nat_pert = np.expm1(X_pert)

    mean_nat_ctrl = nat_ctrl.mean(axis=0)
    mean_nat_pert = nat_pert.mean(axis=0)

    # True log2FC: log2((mean_nat_pert + 1) / (mean_nat_ctrl + 1))
    log2fc = np.log2(mean_nat_pert + 1) - np.log2(mean_nat_ctrl + 1)

    # Log-space means (for reference)
    mean_log_ctrl = X_ctrl.mean(axis=0)
    mean_log_pert = X_pert.mean(axis=0)

    # Paired per-cell deltas in natural scale
    deltas = nat_pert - nat_ctrl  # (n_cells, n_genes)

    # Paired t-test: is the mean delta significantly different from 0?
    # Using scipy.stats.ttest_1samp on the deltas
    t_stats, pvals = stats.ttest_1samp(deltas, popmean=0, axis=0)

    # Handle constant columns (zero variance → nan p-value)
    pvals = np.nan_to_num(pvals, nan=1.0)

    # BH FDR correction
    fdr = _fdr_bh(pvals)

    # Confidence intervals
    if n_bootstrap > 0:
        ci_low, ci_high = _bootstrap_ci(nat_ctrl, nat_pert, n_bootstrap)
    else:
        # Use t-distribution CI from the paired deltas
        ci_low, ci_high = _ttest_ci(deltas, log2fc, mean_nat_ctrl, mean_nat_pert)

    # Significance: both effect size and statistical test
    is_significant = (np.abs(log2fc) > 0.5) & (fdr < 0.05)

    result = pd.DataFrame({
        "gene": gene_names,
        "log2fc": log2fc,
        "log2fc_ci_low": ci_low,
        "log2fc_ci_high": ci_high,
        "pval": pvals,
        "fdr": fdr,
        "mean_pert": mean_log_pert,
        "mean_ctrl": mean_log_ctrl,
        "is_significant": is_significant,
    }).set_index("gene").sort_values("log2fc", key=abs, ascending=False)

    return result


def extract_paired_matrices(pred, target_gene):
    """Extract matched control and perturbed expression matrices.

    Args:
        pred: AnnData with target_gene and original_barcode in obs
        target_gene: name of the perturbation to extract

    Returns:
        X_ctrl, X_pert: np.ndarrays with rows matched by original_barcode
        gene_names: list of gene names
    """
    ctrl_mask = pred.obs["target_gene"] == "non-targeting"
    pert_mask = pred.obs["target_gene"] == target_gene

    ctrl = pred[ctrl_mask]
    pert = pred[pert_mask]

    # Match by original_barcode
    ctrl_barcodes = ctrl.obs["original_barcode"].values
    pert_barcodes = pert.obs["original_barcode"].values

    # Find common barcodes in same order
    common = set(ctrl_barcodes) & set(pert_barcodes)
    if not common:
        raise ValueError(f"No matching barcodes between control and {target_gene}")

    # Build index maps
    ctrl_idx = {bc: i for i, bc in enumerate(ctrl_barcodes)}
    pert_idx = {bc: i for i, bc in enumerate(pert_barcodes)}

    # Order by sorted common barcodes for reproducibility
    ordered = sorted(common)
    ctrl_order = [ctrl_idx[bc] for bc in ordered]
    pert_order = [pert_idx[bc] for bc in ordered]

    X_ctrl = ctrl.X[ctrl_order]
    X_pert = pert.X[pert_order]

    if hasattr(X_ctrl, "toarray"):
        X_ctrl = X_ctrl.toarray()
        X_pert = X_pert.toarray()

    X_ctrl = np.asarray(X_ctrl, dtype=np.float32)
    X_pert = np.asarray(X_pert, dtype=np.float32)

    return X_ctrl, X_pert, pred.var_names.tolist()


def _fdr_bh(pvals):
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    ranked = np.argsort(pvals)
    fdr = np.ones(n)
    for i, idx in enumerate(ranked):
        fdr[idx] = pvals[idx] * n / (i + 1)
    # Enforce monotonicity (from largest rank down)
    for i in range(n - 2, -1, -1):
        idx = ranked[i]
        idx_next = ranked[i + 1]
        fdr[idx] = min(fdr[idx], fdr[idx_next])
    fdr = np.clip(fdr, 0, 1)
    return fdr


def _ttest_ci(deltas, log2fc, mean_nat_ctrl, mean_nat_pert, alpha=0.05):
    """Approximate CI for log2FC using t-distribution on paired deltas."""
    n = deltas.shape[0]
    se = deltas.std(axis=0, ddof=1) / np.sqrt(n)
    t_crit = stats.t.ppf(1 - alpha / 2, df=n - 1)

    # CI on mean delta in natural scale
    delta_mean = deltas.mean(axis=0)
    ci_low_nat = delta_mean - t_crit * se
    ci_high_nat = delta_mean + t_crit * se

    # Convert to log2FC scale (approximate via delta method)
    # log2fc ≈ log2(mean_nat_pert + 1) - log2(mean_nat_ctrl + 1)
    # Perturb the natural-scale mean by CI bounds
    # Clamp to avoid log2 of negative values
    ci_low = np.log2(np.maximum(mean_nat_ctrl + ci_low_nat + 1, 1e-10)) - np.log2(mean_nat_ctrl + 1)
    ci_high = np.log2(np.maximum(mean_nat_ctrl + ci_high_nat + 1, 1e-10)) - np.log2(mean_nat_ctrl + 1)

    return ci_low, ci_high


def _bootstrap_ci(nat_ctrl, nat_pert, n_bootstrap, alpha=0.05):
    """Bootstrap CI for log2FC."""
    n_cells = nat_ctrl.shape[0]
    rng = np.random.default_rng(42)

    boot_log2fc = np.zeros((n_bootstrap, nat_ctrl.shape[1]), dtype=np.float32)
    for i in range(n_bootstrap):
        idx = rng.choice(n_cells, size=n_cells, replace=True)
        mc = nat_ctrl[idx].mean(axis=0)
        mp = nat_pert[idx].mean(axis=0)
        boot_log2fc[i] = np.log2(mp + 1) - np.log2(mc + 1)

    ci_low = np.percentile(boot_log2fc, 100 * alpha / 2, axis=0)
    ci_high = np.percentile(boot_log2fc, 100 * (1 - alpha / 2), axis=0)
    return ci_low, ci_high
