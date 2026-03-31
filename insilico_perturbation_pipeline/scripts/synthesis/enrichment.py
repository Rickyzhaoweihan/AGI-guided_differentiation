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

    Uses the full DE result vector with dc.mt.ora, which internally selects
    top/bottom features for enrichment testing against the full gene background.

    Args:
        de_results: DataFrame from Phase 2 with 'log2fc', 'fdr', 'is_significant'

    Returns:
        dict with 'result' DataFrame from ORA
    """
    import decoupler as dc
    import anndata as adat

    # Get MSigDB Hallmark gene sets
    try:
        msigdb = dc.op.hallmark()
    except Exception:
        return {"error": "Failed to fetch hallmark gene sets"}

    # Pass the full log2fc vector as an AnnData — decoupleR ORA selects
    # top positive and bottom negative features internally
    mat = adat.AnnData(
        X=de_results["log2fc"].values.reshape(1, -1).astype(np.float32),
        var=pd.DataFrame(index=de_results.index),
        obs=pd.DataFrame(index=["perturbation"]),
    )

    try:
        dc.mt.ora(mat, net=msigdb)
        # Extract results from obsm
        est_key = next((k for k in mat.obsm if "score" in k and "ora" in k), None)
        pval_key = next((k for k in mat.obsm if "padj" in k and "ora" in k), None)

        if est_key:
            est_data = mat.obsm[est_key]
            if hasattr(est_data, 'iloc'):
                scores = est_data.iloc[0].sort_values(ascending=False)
            else:
                scores = pd.Series(est_data[0], dtype=float).sort_values(ascending=False)

            pvals = None
            if pval_key:
                pv_data = mat.obsm[pval_key]
                if hasattr(pv_data, 'iloc'):
                    pvals = pv_data.iloc[0]
                else:
                    pvals = pd.Series(pv_data[0], index=scores.index, dtype=float)

            result_df = pd.DataFrame({
                "term": scores.index,
                "score": scores.values,
                "padj": pvals.reindex(scores.index).values if pvals is not None else np.nan,
            }).set_index("term")
            return {"result": result_df}
        else:
            return {"error": f"No ORA results in obsm: {list(mat.obsm.keys())}"}
    except Exception as e:
        return {"error": str(e)}


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
