"""Phase 3a-b: GSEA (pre-ranked) and Over-Representation Analysis (ORA).

Uses gseapy for GSEA and decoupleR for ORA.
"""

import numpy as np
import pandas as pd


def run_gsea(de_results, gene_sets="MSigDB_Hallmark_2020", output_dir=None,
             min_size=10, max_size=500, permutation_num=1000):
    """Run pre-ranked GSEA on log2FC values.

    Args:
        de_results: DataFrame from Phase 2 with 'log2fc' column
        gene_sets: Enrichr gene set library name
        output_dir: directory to save gseapy output (None = temp dir)
        min_size, max_size, permutation_num: gseapy parameters

    Returns:
        pd.DataFrame with GSEA results (Term, NES, pval, FDR, etc.)
    """
    import gseapy as gp

    rnk = de_results["log2fc"].dropna()
    if rnk.empty:
        return pd.DataFrame()

    outdir = str(output_dir) if output_dir else None

    res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        outdir=outdir,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        seed=42,
        verbose=False,
        no_plot=True,
    )

    result = res.res2d.copy()
    result["significant"] = result["FDR q-val"] < 0.25
    return result


def run_ora(de_results, log2fc_threshold=0.5, fdr_threshold=0.05):
    """Run Over-Representation Analysis on significant DE genes using decoupleR.

    Args:
        de_results: DataFrame from Phase 2 with 'log2fc', 'fdr', 'is_significant'

    Returns:
        dict with 'up' and 'down' ORA results DataFrames
    """
    import decoupler as dc

    sig = de_results[de_results["is_significant"]]
    if sig.empty:
        # Fall back to top genes by effect size if nothing passes FDR
        sig = de_results.nlargest(100, "log2fc")

    # Get MSigDB Hallmark gene sets
    try:
        msigdb = dc.op.hallmark()
    except Exception:
        return {"up": pd.DataFrame(), "down": pd.DataFrame(), "error": "Failed to fetch hallmark gene sets"}

    # Upregulated genes
    up_genes = sig[sig["log2fc"] > 0].index.tolist()
    # Downregulated genes
    down_genes = sig[sig["log2fc"] < 0].index.tolist()

    results = {}
    for direction, genes in [("up", up_genes), ("down", down_genes)]:
        if not genes:
            results[direction] = pd.DataFrame()
            continue

        try:
            ora_res = dc.mt.ora(
                mat=pd.DataFrame(
                    np.ones((1, len(genes))),
                    columns=genes,
                ),
                net=msigdb,
            )
            # Extract results
            if hasattr(ora_res, 'shape') and ora_res.shape[0] > 0:
                results[direction] = ora_res
            else:
                results[direction] = pd.DataFrame()
        except Exception as e:
            results[direction] = pd.DataFrame()
            results[f"{direction}_error"] = str(e)

    return results


def run_multi_gsea(de_results, output_dir=None):
    """Run GSEA against multiple gene set libraries.

    Returns:
        dict mapping library name → GSEA results DataFrame
    """
    libraries = [
        "MSigDB_Hallmark_2020",
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
    ]

    results = {}
    for lib in libraries:
        try:
            sub_dir = output_dir / f"gsea_{lib}" if output_dir else None
            res = run_gsea(de_results, gene_sets=lib, output_dir=sub_dir)
            # Only keep significant results for summary
            sig = res[res["significant"]] if "significant" in res.columns else res
            results[lib] = sig
        except Exception as e:
            results[lib] = pd.DataFrame()
            print(f"  GSEA {lib} failed: {e}")

    return results
