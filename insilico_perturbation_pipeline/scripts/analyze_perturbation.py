#!/usr/bin/env python3
"""Perturbation interpretation and synthesis pipeline.

Analyzes in-silico perturbation predictions through 7 phases:
  1. Sanity checks (tiered QC)
  2. Differential expression (paired per-cell deltas)
  3. Functional enrichment (GSEA + ORA + TF activity + PROGENy)
  5. Cell state scoring (proliferation, apoptosis, stress)
  6. Literature integration (PubMed)
  7. Multi-evidence synthesis (structured report)

Usage:
    python analyze_perturbation.py \
        --predictions predictions_batch000.h5ad \
        --output-dir results/adult_gut/stem_cells/ \
        --tissue "adult gut" \
        --cell-type "Stem cells"

    # With expert knowledge
    python analyze_perturbation.py \
        --predictions predictions_batch000.h5ad \
        --config expert_knowledge.yaml \
        --output-dir results/ \
        --tissue "adult gut" \
        --cell-type "Stem cells"

    # Single gene, no literature (faster)
    python analyze_perturbation.py \
        --predictions predictions_batch000.h5ad \
        --target-gene CDX2 \
        --no-literature \
        --output-dir results/
"""

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad
import numpy as np
import pandas as pd

# Add synthesis module to path
sys.path.insert(0, str(Path(__file__).parent))

from synthesis.sanity_checks import run_sanity_checks
from synthesis.differential_expression import compute_de, extract_paired_matrices
from synthesis.enrichment import run_gsea, run_multi_gsea, run_ora
from synthesis.regulatory import run_tf_activity, run_pathway_activity, get_collectri_targets
from synthesis.cell_state import run_cell_state_scoring
from synthesis.literature import search_for_perturbation, get_gene_summary
from synthesis.report import generate_report, save_results


def load_expert_config(config_path):
    """Load expert knowledge YAML config."""
    if config_path is None:
        return {}
    try:
        import yaml
        with open(config_path) as f:
            return yaml.safe_load(f) or {}
    except Exception as e:
        print(f"  Warning: Failed to load config {config_path}: {e}")
        return {}


def analyze_one_perturbation(pred, gene, args, expert_config):
    """Run full analysis pipeline for one perturbation gene.

    Returns:
        dict with all results, or None if extraction fails
    """
    tissue = args.tissue or "unknown"
    cell_type = args.cell_type or "unknown"
    force = args.force_full_report

    print(f"\n{'='*60}")
    print(f"  {gene} knockdown in {cell_type}")
    print(f"{'='*60}")

    # ── Extract paired matrices ──────────────────────────────────
    print("  Extracting paired cells...", end=" ", flush=True)
    try:
        X_ctrl, X_pert, gene_names = extract_paired_matrices(pred, gene)
        print(f"{X_ctrl.shape[0]} matched cells")
    except ValueError as e:
        print(f"FAILED: {e}")
        return None

    # ── Phase 2: Differential Expression ─────────────────────────
    print("  Phase 2: Differential expression...", end=" ", flush=True)
    de_results = compute_de(X_ctrl, X_pert, gene_names)
    n_sig = de_results["is_significant"].sum()
    print(f"{n_sig} significant genes")

    # ── Phase 1: Sanity Checks ───────────────────────────────────
    print("  Phase 1: Sanity checks...", end=" ", flush=True)
    # Get known targets from expert config and/or CollecTRI
    expert_targets = (expert_config.get("known_targets", {}) or {}).get(gene)
    collectri_targets = None
    if not args.no_literature:
        try:
            collectri_targets = get_collectri_targets(gene)
        except Exception:
            pass

    qc_result = run_sanity_checks(
        de_results["log2fc"], gene, gene_names,
        known_targets=expert_targets,
        collectri_targets=collectri_targets,
    )
    print(f"{qc_result['tier']}")
    for check in qc_result["checks"]:
        print(f"    {check['check']}: {check['result']}")

    # Check if we should stop
    if qc_result["tier"] == "FAIL" and not force:
        print(f"  QC FAILED — generating abbreviated report only")
        report = generate_report(
            gene, cell_type, tissue, qc_result, de_results,
            None, None, None, None, None, None, None,
        )
        out = save_results(
            args.output_dir, gene, de_results, qc_result,
            None, None, None, None, None, None, report,
        )
        print(f"  Saved to {out}/")
        return {"tier": "FAIL", "output_dir": str(out)}

    # ── Phase 3a: GSEA ───────────────────────────────────────────
    gsea_results = None
    if not args.no_gsea:
        print("  Phase 3a: GSEA...", end=" ", flush=True)
        try:
            gsea_dir = Path(args.output_dir) / gene / "gsea"
            gsea_results = run_multi_gsea(de_results, output_dir=gsea_dir)
            total = sum(len(df) for df in gsea_results.values() if df is not None and not df.empty)
            print(f"{total} significant pathways")
        except Exception as e:
            print(f"failed: {e}")

    # ── Phase 3b: ORA ────────────────────────────────────────────
    ora_results = None
    print("  Phase 3b: ORA...", end=" ", flush=True)
    try:
        ora_results = run_ora(de_results)
        print("done")
    except Exception as e:
        print(f"failed: {e}")

    # ── Phase 3c: TF Activity ────────────────────────────────────
    tf_activity = None
    print("  Phase 3c: TF activity...", end=" ", flush=True)
    try:
        tf_activity = run_tf_activity(de_results)
        n_tf = len(tf_activity) if tf_activity is not None else 0
        print(f"{n_tf} TFs scored")
    except Exception as e:
        print(f"failed: {e}")

    # ── Phase 3d: Pathway Activity (PROGENy) ─────────────────────
    pathway_activity = None
    print("  Phase 3d: PROGENy pathway activity...", end=" ", flush=True)
    try:
        pathway_activity = run_pathway_activity(de_results)
        n_pw = len(pathway_activity) if pathway_activity is not None else 0
        print(f"{n_pw} pathways scored")
    except Exception as e:
        print(f"failed: {e}")

    # ── Phase 5: Cell State Scoring ──────────────────────────────
    cell_state_scores = None
    print("  Phase 5: Cell state scoring...", end=" ", flush=True)
    try:
        custom_sets = expert_config.get("custom_gene_sets")
        cell_state_scores = run_cell_state_scoring(X_ctrl, X_pert, gene_names, custom_sets)
        n_changed = len(cell_state_scores[cell_state_scores["delta"].abs() > 0.05]) \
            if cell_state_scores is not None else 0
        print(f"{n_changed} states changed")
    except Exception as e:
        print(f"failed: {e}")

    # ── Phase 6: Literature ──────────────────────────────────────
    literature_results = None
    gene_summary = None
    if not args.no_literature:
        print("  Phase 6: Literature search...", end=" ", flush=True)
        try:
            # Get top pathways for literature queries
            top_pathways = []
            if gsea_results:
                for df in gsea_results.values():
                    if df is not None and not df.empty and "NES" in df.columns:
                        df_copy = df.copy()
                        df_copy["NES"] = pd.to_numeric(df_copy["NES"], errors="coerce")
                        top = df_copy.dropna(subset=["NES"]).nlargest(2, "NES")
                        if "Term" in top.columns:
                            top_pathways.extend(top["Term"].tolist())
                        else:
                            top_pathways.extend(top.index.tolist())

            literature_results = search_for_perturbation(
                gene, tissue, cell_type,
                top_pathways=top_pathways[:4],
                api_key=os.environ.get("NCBI_API_KEY"),
            )
            gene_summary = get_gene_summary(gene, api_key=os.environ.get("NCBI_API_KEY"))
            total_articles = sum(
                len([a for a in v.get("articles", []) if "error" not in a])
                for v in literature_results.values()
            )
            print(f"{total_articles} articles found")
        except Exception as e:
            print(f"failed: {e}")
    else:
        print("  Phase 6: Skipped (--no-literature)")

    # ── Phase 7: Report Generation ───────────────────────────────
    print("  Phase 7: Generating report...", end=" ", flush=True)
    report = generate_report(
        gene, cell_type, tissue, qc_result, de_results,
        gsea_results, ora_results, tf_activity, pathway_activity,
        cell_state_scores, literature_results, gene_summary,
    )
    out = save_results(
        args.output_dir, gene, de_results, qc_result, gsea_results,
        tf_activity, pathway_activity, cell_state_scores,
        literature_results, gene_summary, report,
    )
    print(f"done → {out}/")

    return {
        "tier": qc_result["tier"],
        "n_sig_genes": int(n_sig),
        "output_dir": str(out),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Perturbation interpretation and synthesis pipeline"
    )
    parser.add_argument("--predictions", "-p", required=True,
                        help="Predictions h5ad file (output of state tx infer)")
    parser.add_argument("--config", "-c", default=None,
                        help="Expert knowledge YAML config (optional)")
    parser.add_argument("--output-dir", "-o", required=True,
                        help="Output directory for reports and results")
    parser.add_argument("--tissue", "-t", default=None,
                        help="Tissue context (e.g., 'adult gut', 'human islets')")
    parser.add_argument("--cell-type", default=None,
                        help="Cell type (e.g., 'Stem cells')")
    parser.add_argument("--target-gene", "-g", default=None, nargs="*",
                        help="Analyze specific gene(s) only (default: all)")
    parser.add_argument("--no-literature", action="store_true",
                        help="Skip literature search (faster, offline)")
    parser.add_argument("--no-gsea", action="store_true",
                        help="Skip GSEA (faster)")
    parser.add_argument("--force-full-report", action="store_true",
                        help="Generate full report even for QC-failed perturbations")
    parser.add_argument("--n-bootstrap", type=int, default=0,
                        help="Bootstrap iterations for CI (0 = use t-test CI, faster)")

    args = parser.parse_args()

    # Load data
    print(f"Loading {args.predictions}...")
    pred = ad.read_h5ad(args.predictions)
    print(f"  {pred.n_obs:,} cells x {pred.n_vars:,} genes")

    # Load expert config
    expert_config = load_expert_config(args.config)
    if expert_config:
        print(f"  Expert config: {args.config}")
        if "custom_gene_sets" in expert_config:
            print(f"    Custom gene sets: {list(expert_config['custom_gene_sets'].keys())}")
        if "known_targets" in expert_config:
            print(f"    Known targets for: {list(expert_config['known_targets'].keys())}")

    # Determine which perturbations to analyze
    all_perts = [g for g in pred.obs["target_gene"].unique() if g != "non-targeting"]
    if args.target_gene:
        perts = [g for g in args.target_gene if g in all_perts]
        missing = [g for g in args.target_gene if g not in all_perts]
        if missing:
            print(f"  Warning: genes not found in predictions: {missing}")
    else:
        # Focus genes from expert config, or all
        focus = expert_config.get("focus_genes")
        if focus:
            perts = [g for g in focus if g in all_perts]
            print(f"  Using focus genes from config: {perts}")
        else:
            perts = all_perts

    print(f"\nAnalyzing {len(perts)} perturbations...")
    os.makedirs(args.output_dir, exist_ok=True)

    # Run analysis for each perturbation
    summary = []
    for gene in perts:
        result = analyze_one_perturbation(pred, gene, args, expert_config)
        if result:
            result["gene"] = gene
            summary.append(result)

    # Write cross-perturbation summary
    if summary:
        summary_df = pd.DataFrame(summary)
        summary_path = os.path.join(args.output_dir, "summary.csv")
        summary_df.to_csv(summary_path, index=False)
        print(f"\n{'='*60}")
        print(f"Summary: {len(summary)} perturbations analyzed")
        print(f"  PASS: {sum(1 for s in summary if s['tier'] == 'PASS')}")
        print(f"  WARN: {sum(1 for s in summary if s['tier'] == 'WARN')}")
        print(f"  FAIL: {sum(1 for s in summary if s['tier'] == 'FAIL')}")
        print(f"Results: {args.output_dir}")
        print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
