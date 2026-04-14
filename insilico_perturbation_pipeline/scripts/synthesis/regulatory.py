"""Phase 3c-d: TF activity inference (CollecTRI + ULM) and pathway activity (PROGENy).

Uses decoupleR for both analyses.
"""

import numpy as np
import pandas as pd


def run_tf_activity(de_results, organism="human"):
    """Infer TF activity from expression changes using decoupleR ULM + CollecTRI.

    Args:
        de_results: DataFrame from Phase 2 with 'log2fc' column
        organism: 'human' or 'mouse'

    Returns:
        pd.DataFrame with columns: tf, activity_score, pval, significant
    """
    import decoupler as dc

    # Fetch CollecTRI TF-target network
    try:
        collectri = dc.op.collectri(organism=organism)
    except Exception as e:
        print(f"  Failed to fetch CollecTRI: {e}")
        return pd.DataFrame()

    # decoupleR v2 expects AnnData input
    import anndata as adat
    mat = adat.AnnData(
        X=de_results["log2fc"].values.reshape(1, -1).astype(np.float32),
        var=pd.DataFrame(index=de_results.index),
        obs=pd.DataFrame(index=["perturbation"]),
    )

    try:
        dc.mt.ulm(mat, net=collectri)
        # decoupleR v2.1 stores results as score_ulm / padj_ulm in obsm
        est_key = next((k for k in mat.obsm if "score" in k and "ulm" in k), None)
        pval_key = next((k for k in mat.obsm if "padj" in k and "ulm" in k), None)

        if est_key is None:
            print(f"  ULM: no score key in obsm: {list(mat.obsm.keys())}")
            return pd.DataFrame()

        est_data = mat.obsm[est_key]
        # Column names stored in the DataFrame if it's a DataFrame, or in uns
        if hasattr(est_data, 'columns'):
            col_names = est_data.columns.tolist()
            est_vals = est_data.iloc[0]
        else:
            col_names = [f"TF_{i}" for i in range(est_data.shape[1])]
            est_vals = pd.Series(est_data[0], index=col_names)

        tf_scores = est_vals.sort_values(key=abs, ascending=False)

        if pval_key is not None:
            pv_data = mat.obsm[pval_key]
            if hasattr(pv_data, 'iloc'):
                tf_pvals = pv_data.iloc[0].reindex(tf_scores.index)
            else:
                tf_pvals = pd.Series(pv_data[0], index=col_names).reindex(tf_scores.index)
        else:
            tf_pvals = pd.Series(np.nan, index=tf_scores.index)

        tf_df = pd.DataFrame({
            "tf": tf_scores.index,
            "activity_score": tf_scores.values.astype(float),
            "pval": tf_pvals.reindex(tf_scores.index).values.astype(float),
        }).set_index("tf")
        tf_df["significant"] = tf_df["pval"] < 0.05
        return tf_df
    except Exception as e:
        print(f"  ULM failed: {e}")
        import traceback; traceback.print_exc()

    return pd.DataFrame()


def run_pathway_activity(de_results, organism="human"):
    """Score pathway activity using decoupleR + PROGENy model.

    PROGENy scores 14 signaling pathways based on downstream transcriptional
    footprints (not pathway membership), capturing pathway *activity*.

    Args:
        de_results: DataFrame from Phase 2 with 'log2fc' column
        organism: 'human' or 'mouse'

    Returns:
        pd.DataFrame with columns: pathway, activity_score, pval
    """
    import decoupler as dc

    # Fetch PROGENy model
    try:
        progeny = dc.op.progeny(organism=organism, top=500)
    except Exception as e:
        print(f"  Failed to fetch PROGENy model: {e}")
        return pd.DataFrame()

    import anndata as adat
    mat = adat.AnnData(
        X=de_results["log2fc"].values.reshape(1, -1).astype(np.float32),
        var=pd.DataFrame(index=de_results.index),
        obs=pd.DataFrame(index=["perturbation"]),
    )

    try:
        dc.mt.mlm(mat, net=progeny)
        est_key = next((k for k in mat.obsm if "score" in k and "mlm" in k), None)
        pval_key = next((k for k in mat.obsm if "padj" in k and "mlm" in k), None)

        if est_key is None:
            print(f"  MLM: no score key in obsm: {list(mat.obsm.keys())}")
            return pd.DataFrame()

        est_data = mat.obsm[est_key]
        if hasattr(est_data, 'columns'):
            pw_scores = est_data.iloc[0].sort_values(key=abs, ascending=False)
        else:
            pw_scores = pd.Series(est_data[0], dtype=float).sort_values(key=abs, ascending=False)

        if pval_key is not None:
            pv_data = mat.obsm[pval_key]
            if hasattr(pv_data, 'iloc'):
                pw_pvals = pv_data.iloc[0].reindex(pw_scores.index)
            else:
                pw_pvals = pd.Series(pv_data[0], index=pw_scores.index, dtype=float)
        else:
            pw_pvals = pd.Series(np.nan, index=pw_scores.index)

        pw_df = pd.DataFrame({
            "pathway": pw_scores.index,
            "activity_score": pw_scores.values.astype(float),
            "pval": pw_pvals.values.astype(float),
        }).set_index("pathway")
        return pw_df
    except Exception as e:
        print(f"  MLM failed: {e}")
        import traceback; traceback.print_exc()

    return pd.DataFrame()


def get_collectri_targets(gene, organism="human"):
    """Get known TF targets from CollecTRI for sanity checking.

    Returns:
        dict with 'activates' and 'represses' lists, or None if gene not found
    """
    import decoupler as dc

    try:
        collectri = dc.op.collectri(organism=organism)
    except Exception:
        return None

    # CollecTRI has columns: source, target, weight (1 = activation, -1 = repression)
    tf_rows = collectri[collectri["source"] == gene]
    if tf_rows.empty:
        return None

    activates = tf_rows[tf_rows["weight"] > 0]["target"].tolist()
    represses = tf_rows[tf_rows["weight"] < 0]["target"].tolist()

    return {"activates": activates, "represses": represses}
