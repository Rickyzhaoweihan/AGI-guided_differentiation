#!/usr/bin/env python3
"""
GSEA analysis on in-silico perturbation predictions.

Takes perturbation_predictions.h5ad (output of run_insilico_perturbation.py),
computes predicted log2FC per gene (perturbed vs. non-targeting control),
and runs pre-ranked GSEA using gseapy.

Usage:
    python GSEA_skill/run_gsea.py \
        --predictions insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
        --output-dir  GSEA_skill/results \
        --gene-sets   MSigDB_Hallmark_2020
"""

import argparse
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


# ── Core analysis ─────────────────────────────────────────────────────────────

def compute_log2fc(pred, target_gene, cell_type=None):
    """
    Compute per-gene expression changes: perturbed vs. non-targeting control.

    Predictions are log1p-normalized. We compute:
    - log2fc: true log2 fold change from natural-scale means:
              log2(mean(expm1(pert)) + 1) - log2(mean(expm1(ctrl)) + 1)
    - log_diff: mean log-expression difference (used as GSEA ranking statistic)

    Returns a DataFrame indexed by gene name, sorted by log_diff descending.
    """
    mask_ctrl = pred.obs["target_gene"] == "non-targeting"
    mask_pert = pred.obs["target_gene"] == target_gene

    if cell_type is not None:
        ct_col = "cell_type" if "cell_type" in pred.obs.columns else "original_cell_type"
        mask_ctrl &= pred.obs[ct_col] == cell_type
        mask_pert &= pred.obs[ct_col] == cell_type

    if mask_ctrl.sum() == 0 or mask_pert.sum() == 0:
        return None

    X_ctrl = pred[mask_ctrl].X
    X_pert = pred[mask_pert].X

    if hasattr(X_ctrl, "toarray"):
        X_ctrl = X_ctrl.toarray()
        X_pert = X_pert.toarray()

    # Mean in log-space (for GSEA ranking — any monotonic statistic works)
    mean_log_ctrl = np.array(X_ctrl.mean(axis=0)).flatten()
    mean_log_pert = np.array(X_pert.mean(axis=0)).flatten()
    log_diff = mean_log_pert - mean_log_ctrl

    # True log2FC from natural-scale means: log2((mean(expm1(pert))+1) / (mean(expm1(ctrl))+1))
    mean_nat_ctrl = np.array(np.expm1(X_ctrl).mean(axis=0)).flatten()
    mean_nat_pert = np.array(np.expm1(X_pert).mean(axis=0)).flatten()
    log2fc = np.log2(mean_nat_pert + 1) - np.log2(mean_nat_ctrl + 1)

    gene_names = pred.var_names.tolist()
    result = pd.DataFrame({
        "gene":      gene_names,
        "log2fc":    log2fc,
        "log_diff":  log_diff,
        "mean_pert": mean_log_pert,
        "mean_ctrl": mean_log_ctrl,
    }).set_index("gene").sort_values("log_diff", ascending=False)

    return result


def run_gsea_preranked(ranked_genes, gene_sets, output_dir, label):
    """Run gseapy prerank GSEA and save results."""
    import gseapy as gp

    rnk = ranked_genes["log2fc"]

    out = Path(output_dir) / f"gsea_{gene_sets}"
    out.mkdir(parents=True, exist_ok=True)

    res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        outdir=str(out),
        min_size=10,
        max_size=500,
        permutation_num=1000,
        seed=42,
        verbose=False,
        no_plot=False,
    )

    return res.res2d


def plot_volcano(ranked_genes, target_gene, output_path, min_abs_lfc=0.1, n_label=20):
    """MA-style plot: true log2FC vs. mean expression."""
    df = ranked_genes.copy()
    df["significant"] = df["log2fc"].abs() > min_abs_lfc

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df["mean_ctrl"], df["log2fc"], c="lightgrey", s=5, alpha=0.5, rasterized=True)

    sig = df[df["significant"]]
    up = sig[sig["log2fc"] > 0]
    dn = sig[sig["log2fc"] < 0]
    ax.scatter(up["mean_ctrl"], up["log2fc"], c="tomato",    s=8, alpha=0.7, label=f"Up ({len(up)})")
    ax.scatter(dn["mean_ctrl"], dn["log2fc"], c="steelblue", s=8, alpha=0.7, label=f"Down ({len(dn)})")

    # Label top genes
    top = pd.concat([
        df.nlargest(n_label // 2, "log2fc"),
        df.nsmallest(n_label // 2, "log2fc"),
    ])
    for gene, row in top.iterrows():
        ax.annotate(gene, (row["mean_ctrl"], row["log2fc"]),
                    fontsize=6, alpha=0.8,
                    xytext=(3, 0), textcoords="offset points")

    ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Mean control expression (log-norm)")
    ax.set_ylabel("log2 fold change (natural-scale means)")
    ax.set_title(f"Predicted DE: {target_gene} perturbation")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="GSEA on in-silico perturbation predictions"
    )
    parser.add_argument("--predictions", "-p", required=True,
                        help="perturbation_predictions.h5ad")
    parser.add_argument("--output-dir", "-o", default="GSEA_skill/results",
                        help="Output directory (default: GSEA_skill/results)")
    parser.add_argument("--gene-sets", default="MSigDB_Hallmark_2020",
                        help="Enrichr gene set library name (default: MSigDB_Hallmark_2020)")
    parser.add_argument("--target-gene", default=None,
                        help="Run for one perturbation gene only (default: all)")
    parser.add_argument("--per-cell-type", action="store_true",
                        help="Run GSEA separately per cell type")
    parser.add_argument("--min-abs-lfc", type=float, default=0.1,
                        help="Min |log2FC| to highlight on volcano (default: 0.1)")
    parser.add_argument("--n-top-pathways", type=int, default=20,
                        help="Top N pathways to include in summary (default: 20)")
    args = parser.parse_args()

    out_root = Path(args.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.predictions}...")
    pred = ad.read_h5ad(args.predictions)
    print(f"  {pred.n_obs:,} cells × {pred.n_vars:,} genes")

    perturbations = [g for g in pred.obs["target_gene"].unique() if g != "non-targeting"]
    if args.target_gene:
        if args.target_gene not in perturbations:
            print(f"ERROR: '{args.target_gene}' not found. Available: {perturbations}")
            return
        perturbations = [args.target_gene]

    ct_col = "cell_type" if "cell_type" in pred.obs.columns else "original_cell_type"
    cell_types = sorted(pred.obs[ct_col].unique()) if args.per_cell_type else [None]

    all_summaries = []

    for gene in perturbations:
        for ct in cell_types:
            label = f"{gene}" + (f"_{ct}" if ct else "")
            print(f"\n── {label} ──")

            out_dir = out_root / label
            out_dir.mkdir(parents=True, exist_ok=True)

            # Compute log2FC
            ranked = compute_log2fc(pred, gene, cell_type=ct)
            if ranked is None:
                print(f"  Skipping: no cells for this condition")
                continue

            ranked.to_csv(out_dir / "de_results.csv")
            print(f"  DE computed: {(ranked['log2fc'].abs() > args.min_abs_lfc).sum()} genes |log2FC| > {args.min_abs_lfc}")

            # Volcano plot
            plot_volcano(ranked, label, out_dir / "volcano.png", args.min_abs_lfc)
            print(f"  Volcano plot saved")

            # GSEA
            try:
                gsea_res = run_gsea_preranked(ranked, args.gene_sets, out_dir, label)
                gsea_res.to_csv(out_dir / f"gsea_{args.gene_sets}_results.csv")

                sig = gsea_res[gsea_res["FDR q-val"] < 0.25].sort_values("NES", ascending=False)
                print(f"  GSEA: {len(sig)} significant pathways (FDR < 0.25)")
                if len(sig) > 0:
                    print(f"  Top pathway: {sig.index[0]}  NES={sig['NES'].iloc[0]:.2f}")

                # Only include FDR-significant pathways in the summary
                sig_for_summary = gsea_res[gsea_res["FDR q-val"] < 0.25]
                for _, row in sig_for_summary.nlargest(args.n_top_pathways, "NES").iterrows():
                    all_summaries.append({
                        "perturbation": gene,
                        "cell_type": ct or "all",
                        "pathway": row.name if hasattr(row, "name") else _,
                        "NES": row["NES"],
                        "pval": row["NOM p-val"],
                        "fdr": row["FDR q-val"],
                    })

            except ImportError:
                print("  gseapy not installed — skipping GSEA. Run: pip install gseapy")
            except Exception as e:
                print(f"  GSEA failed: {e}")

    # Write summary
    if all_summaries:
        summary_df = pd.DataFrame(all_summaries)
        summary_path = out_root / "summary_all_perturbations.csv"
        summary_df.to_csv(summary_path, index=False)
        print(f"\nSummary saved: {summary_path}")

    print(f"\nDone. Results in {out_root}/")


if __name__ == "__main__":
    main()
