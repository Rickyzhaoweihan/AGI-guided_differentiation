"""Phase 7: Multi-evidence synthesis and structured report generation.

Generates markdown reports per (perturbation, cell_type) with convergence
assessment across all analysis phases. No composite score — evidence
dimensions presented separately.
"""

import json
from pathlib import Path

import pandas as pd


def generate_report(gene, cell_type, tissue, qc_result, de_results,
                    gsea_results, ora_results, tf_activity, pathway_activity,
                    cell_state_scores, literature_results, gene_summary):
    """Generate a structured markdown report for one perturbation.

    Returns:
        str: markdown-formatted report
    """
    tier = qc_result["tier"]
    sections = []

    # Header
    sections.append(f"# Perturbation Report: {gene} knockdown in {cell_type}\n")
    sections.append(f"**Tissue:** {tissue}  ")
    sections.append(f"**QC Status:** {tier}\n")

    # QC banner
    if tier == "FAIL":
        sections.append("> **WARNING: QC FAILED** — Predictions flagged as unreliable. "
                       "Only raw DE results shown below. No synthesis narrative generated.\n")
    elif tier == "WARN":
        warns = [c for c in qc_result["checks"] if c["result"] == "WARN"]
        for w in warns:
            sections.append(f"> **Warning ({w['check']}):** {w['details']}\n")

    # Sanity checks table
    sections.append("## Sanity Checks\n")
    sections.append("| Check | Result | Details |")
    sections.append("|-------|--------|---------|")
    for check in qc_result["checks"]:
        sections.append(f"| {check['check']} | {check['result']} | {check['details'][:100]} |")
    sections.append("")

    # Differential expression
    sections.append("## Differential Expression\n")
    if de_results is not None and not de_results.empty:
        n_sig = de_results["is_significant"].sum() if "is_significant" in de_results.columns else 0
        sections.append(f"- **{n_sig}** genes significantly changed (|log2FC| > 0.5, FDR < 0.05)\n")

        # Strongest up (largest positive log2fc) and down (most negative log2fc)
        top_up = de_results[de_results["log2fc"] > 0].nlargest(10, "log2fc")
        top_down = de_results[de_results["log2fc"] < 0].nsmallest(10, "log2fc")

        if not top_up.empty:
            sections.append("**Top upregulated:**\n")
            sections.append("| Gene | log2FC | 95% CI | FDR |")
            sections.append("|------|--------|--------|-----|")
            for g, row in top_up.iterrows():
                ci = f"[{row.get('log2fc_ci_low', 0):.3f}, {row.get('log2fc_ci_high', 0):.3f}]"
                fdr = f"{row.get('fdr', 1):.2e}" if row.get('fdr', 1) < 1 else "ns"
                sig = " *" if row.get("is_significant", False) else ""
                sections.append(f"| {g}{sig} | +{row['log2fc']:.3f} | {ci} | {fdr} |")
            sections.append("")

        if not top_down.empty:
            sections.append("**Top downregulated:**\n")
            sections.append("| Gene | log2FC | 95% CI | FDR |")
            sections.append("|------|--------|--------|-----|")
            for g, row in top_down.iterrows():
                ci = f"[{row.get('log2fc_ci_low', 0):.3f}, {row.get('log2fc_ci_high', 0):.3f}]"
                fdr = f"{row.get('fdr', 1):.2e}" if row.get('fdr', 1) < 1 else "ns"
                sig = " *" if row.get("is_significant", False) else ""
                sections.append(f"| {g}{sig} | {row['log2fc']:.3f} | {ci} | {fdr} |")
            sections.append("")

    # Stop here for QC-failed perturbations
    if tier == "FAIL":
        sections.append("---\n*Report truncated due to QC failure. "
                       "Use `--force-full-report` for exploratory analysis.*\n")
        return "\n".join(sections)

    # GSEA results
    sections.append("## Pathway Enrichment\n")
    sections.append("### GSEA (pre-ranked, FDR < 0.25)\n")
    if gsea_results:
        for lib_name, gsea_df in gsea_results.items():
            if gsea_df is not None and not gsea_df.empty:
                gsea_df["NES"] = pd.to_numeric(gsea_df["NES"], errors="coerce")
                gsea_df["FDR q-val"] = pd.to_numeric(gsea_df["FDR q-val"], errors="coerce")
                sig_df = gsea_df[gsea_df["FDR q-val"] < 0.25].dropna(subset=["NES"])
                if not sig_df.empty:
                    sections.append(f"**{lib_name}** ({len(sig_df)} significant):\n")
                    # Show top positive AND top negative NES separately
                    pos_sig = sig_df[sig_df["NES"] > 0].nlargest(5, "NES")
                    neg_sig = sig_df[sig_df["NES"] < 0].nsmallest(5, "NES")
                    top_both = pd.concat([pos_sig, neg_sig]).sort_values("NES", key=abs, ascending=False)
                    sections.append("| Pathway | NES | FDR | Direction |")
                    sections.append("|---------|-----|-----|-----------|")
                    for _, row in top_both.iterrows():
                        name = row.name if hasattr(row, "name") and isinstance(row.name, str) else row.get("Term", "")
                        direction = "UP" if row["NES"] > 0 else "DOWN"
                        sections.append(f"| {name} | {row['NES']:.2f} | {row['FDR q-val']:.3f} | {direction} |")
                    sections.append("")
                else:
                    sections.append(f"**{lib_name}**: No significant pathways (FDR < 0.25)\n")
    else:
        sections.append("GSEA not run.\n")

    # ORA results
    if ora_results:
        sections.append("### Over-Representation Analysis\n")
        for direction in ["up", "down"]:
            ora_df = ora_results.get(direction)
            if ora_df is not None and not ora_df.empty:
                sections.append(f"**{direction.capitalize()}regulated gene ORA:**\n")
                sections.append("(ORA results available in output files)\n")
        sections.append("")

    # TF activity
    sections.append("## Regulatory Analysis\n")
    sections.append("### TF Activity (decoupleR ULM + CollecTRI)\n")
    if tf_activity is not None and not tf_activity.empty:
        top_tfs = tf_activity.head(15)
        sections.append("| TF | Activity Score | p-value |")
        sections.append("|----|---------------|---------|")
        for tf, row in top_tfs.iterrows():
            pval = f"{row['pval']:.2e}" if pd.notna(row['pval']) else "n/a"
            direction = "activated" if row["activity_score"] > 0 else "repressed"
            sections.append(f"| {tf} | {row['activity_score']:+.3f} ({direction}) | {pval} |")
        sections.append("")
    else:
        sections.append("TF activity analysis not available.\n")

    # Pathway activity (PROGENy)
    sections.append("### Pathway Activity (PROGENy)\n")
    if pathway_activity is not None and not pathway_activity.empty:
        sections.append("| Pathway | Activity Score | p-value |")
        sections.append("|---------|---------------|---------|")
        for pw, row in pathway_activity.iterrows():
            pval = f"{row['pval']:.2e}" if pd.notna(row['pval']) else "n/a"
            sections.append(f"| {pw} | {row['activity_score']:+.3f} | {pval} |")
        sections.append("")
    else:
        sections.append("PROGENy pathway analysis not available.\n")

    # Cell state impact
    sections.append("## Cell State Impact\n")
    sections.append("*Note: pro-apoptotic and anti-apoptotic scored separately*\n")
    if cell_state_scores is not None and not cell_state_scores.empty:
        sections.append("| State | Control | Perturbed | Delta | Interpretation |")
        sections.append("|-------|---------|-----------|-------|----------------|")
        for _, row in cell_state_scores.iterrows():
            sections.append(
                f"| {row['state']} | {row['ctrl_score']:.3f} | {row['pert_score']:.3f} "
                f"| {row['delta']:+.3f} | {row['interpretation']} |"
            )
        sections.append("")
    else:
        sections.append("Cell state scoring not available.\n")

    # Literature
    sections.append("## Literature Context\n")
    if gene_summary and "description" in gene_summary:
        sections.append(f"**{gene}**: {gene_summary.get('description', 'N/A')}\n")
        if gene_summary.get("summary"):
            sections.append(f"> {gene_summary['summary'][:300]}\n")

    if literature_results:
        for query_name, query_data in literature_results.items():
            articles = query_data.get("articles", [])
            valid = [a for a in articles if "error" not in a]
            if valid:
                sections.append(f"**{query_name}** ({len(valid)} results):\n")
                for a in valid[:3]:
                    sections.append(
                        f"- [{a['title']}](https://pubmed.ncbi.nlm.nih.gov/{a['pmid']}/) "
                        f"({a.get('authors', '')}, {a.get('year', '')})"
                    )
                    if a.get("snippet"):
                        sections.append(f"  > {a['snippet'][:200]}...")
                sections.append("")
    else:
        sections.append("Literature search not performed.\n")

    # Evidence convergence
    sections.append("## Evidence Convergence\n")
    convergence = _assess_convergence(
        qc_result, de_results, gsea_results, tf_activity,
        pathway_activity, cell_state_scores, literature_results
    )
    sections.append("| Dimension | Signal | Details |")
    sections.append("|-----------|--------|---------|")
    for dim in convergence:
        sections.append(f"| {dim['dimension']} | {dim['signal']} | {dim['details']} |")
    sections.append("")

    return "\n".join(sections)


def _assess_convergence(qc_result, de_results, gsea_results, tf_activity,
                        pathway_activity, cell_state_scores, literature_results):
    """Assess evidence convergence across all phases."""
    convergence = []

    # Mechanistic consistency (Phase 1+2)
    target_check = next((c for c in qc_result["checks"] if c["check"] == "known_targets"), None)
    if target_check and target_check["result"] != "SKIP":
        convergence.append({
            "dimension": "Mechanistic",
            "signal": "+" if target_check["result"] == "PASS" else "-",
            "details": target_check["details"][:80],
        })
    else:
        convergence.append({"dimension": "Mechanistic", "signal": "?", "details": "No known targets to validate"})

    # Functional consistency (Phase 3a-b)
    if gsea_results:
        total_sig = 0
        for df in gsea_results.values():
            if df is not None and not df.empty and "FDR q-val" in df.columns:
                fdr = pd.to_numeric(df["FDR q-val"], errors="coerce")
                total_sig += int((fdr < 0.25).sum())
        convergence.append({
            "dimension": "Functional",
            "signal": "+" if total_sig > 0 else "-",
            "details": f"{total_sig} significant pathways across libraries",
        })
    else:
        convergence.append({"dimension": "Functional", "signal": "?", "details": "GSEA not run"})

    # Regulatory consistency (Phase 3c) — require statistical significance
    if tf_activity is not None and not tf_activity.empty:
        sig_tfs = tf_activity[tf_activity["pval"] < 0.05] if "pval" in tf_activity.columns else pd.DataFrame()
        if not sig_tfs.empty:
            top_tf = sig_tfs.index[0]
            top_score = sig_tfs.iloc[0]["activity_score"]
            convergence.append({
                "dimension": "Regulatory",
                "signal": "+",
                "details": f"{len(sig_tfs)} significant TFs (p<0.05). Top: {top_tf} ({top_score:+.2f})",
            })
        else:
            convergence.append({
                "dimension": "Regulatory",
                "signal": "~",
                "details": f"TF scores computed but none significant (p<0.05)",
            })
    else:
        convergence.append({"dimension": "Regulatory", "signal": "?", "details": "TF activity not available"})

    # Pathway consistency (Phase 3d) — require statistical significance
    if pathway_activity is not None and not pathway_activity.empty:
        sig_pw = pathway_activity[pathway_activity["pval"] < 0.05] if "pval" in pathway_activity.columns else pd.DataFrame()
        if not sig_pw.empty:
            convergence.append({
                "dimension": "Pathway",
                "signal": "+",
                "details": f"{len(sig_pw)} significant pathways (p<0.05)",
            })
        else:
            convergence.append({
                "dimension": "Pathway",
                "signal": "~",
                "details": "Pathway scores computed but none significant (p<0.05)",
            })
    else:
        convergence.append({"dimension": "Pathway", "signal": "?", "details": "PROGENy not available"})

    # Phenotypic consistency (Phase 5)
    if cell_state_scores is not None and not cell_state_scores.empty:
        changed = cell_state_scores[cell_state_scores["delta"].abs() > 0.05]
        convergence.append({
            "dimension": "Phenotypic",
            "signal": "+" if len(changed) > 0 else "~",
            "details": f"{len(changed)} cell states changed",
        })
    else:
        convergence.append({"dimension": "Phenotypic", "signal": "?", "details": "Cell state scoring not available"})

    # Literature consistency (Phase 6)
    if literature_results:
        total_articles = sum(
            len([a for a in v.get("articles", []) if "error" not in a])
            for v in literature_results.values()
        )
        convergence.append({
            "dimension": "Literature",
            "signal": "+" if total_articles > 0 else "-",
            "details": f"{total_articles} relevant publications found",
        })
    else:
        convergence.append({"dimension": "Literature", "signal": "?", "details": "Literature search not performed"})

    return convergence


def save_results(output_dir, gene, de_results, qc_result, gsea_results,
                 tf_activity, pathway_activity, cell_state_scores,
                 literature_results, gene_summary, report_text):
    """Save all analysis results to files."""
    out = Path(output_dir) / gene
    out.mkdir(parents=True, exist_ok=True)

    # DE results
    if de_results is not None:
        de_results.to_csv(out / "de_results.csv")

    # QC
    with open(out / "sanity_checks.json", "w") as f:
        json.dump(qc_result, f, indent=2, default=str)

    # GSEA
    if gsea_results:
        for lib_name, df in gsea_results.items():
            if df is not None and not df.empty:
                df.to_csv(out / f"gsea_{lib_name}.csv")

    # TF activity
    if tf_activity is not None and not tf_activity.empty:
        tf_activity.to_csv(out / "tf_activity.csv")

    # Pathway activity
    if pathway_activity is not None and not pathway_activity.empty:
        pathway_activity.to_csv(out / "pathway_activity.csv")

    # Cell state
    if cell_state_scores is not None and not cell_state_scores.empty:
        cell_state_scores.to_csv(out / "cell_state_scores.csv", index=False)

    # Literature
    if literature_results:
        with open(out / "literature_results.json", "w") as f:
            json.dump(literature_results, f, indent=2, default=str)

    # Gene summary
    if gene_summary:
        with open(out / "gene_summary.json", "w") as f:
            json.dump(gene_summary, f, indent=2, default=str)

    # Report
    with open(out / "report.md", "w") as f:
        f.write(report_text)

    return out
