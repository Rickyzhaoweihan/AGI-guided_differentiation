"""Phase 1: Sanity checks on perturbation predictions.

Tiered QC:
  FAIL — housekeeping instability or zero effect → block synthesis
  WARN — self-knockdown or known-target issues → full report with warnings
  PASS — all checks pass → full report
"""

import numpy as np

HOUSEKEEPING_GENES = ["ACTB", "GAPDH", "RPL13A", "UBC", "B2M", "HPRT1", "TBP"]
HK_THRESHOLD = 0.5
EFFECT_THRESHOLD = 0.001


def run_sanity_checks(log2fc_series, perturbed_gene, gene_names,
                      known_targets=None, collectri_targets=None):
    """Run all sanity checks and return tiered result.

    Args:
        log2fc_series: pd.Series indexed by gene name with log2FC values
        perturbed_gene: name of the perturbed gene
        gene_names: list of all gene names
        known_targets: dict with 'activates' and 'represses' lists (from expert YAML)
        collectri_targets: dict with 'activates' and 'represses' from CollecTRI

    Returns:
        dict with 'tier' (FAIL/WARN/PASS), 'checks' list, 'details' dict
    """
    checks = []
    tier = "PASS"

    # Check 1: Housekeeping stability
    hk_present = [g for g in HOUSEKEEPING_GENES if g in log2fc_series.index]
    if hk_present:
        hk_fc = log2fc_series[hk_present]
        hk_unstable = hk_fc[hk_fc.abs() > HK_THRESHOLD]
        if len(hk_unstable) >= 2:
            tier = "FAIL"
            checks.append({
                "check": "housekeeping_stability",
                "result": "FAIL",
                "details": f"{len(hk_unstable)}/{len(hk_present)} HK genes unstable: "
                           f"{dict(hk_unstable.round(3))}"
            })
        else:
            checks.append({
                "check": "housekeeping_stability",
                "result": "PASS",
                "details": f"All {len(hk_present)} HK genes stable (max |log2FC| = {hk_fc.abs().max():.4f})"
            })
    else:
        checks.append({
            "check": "housekeeping_stability",
            "result": "SKIP",
            "details": "No housekeeping genes found in gene set"
        })

    # Check 2: Effect magnitude
    median_abs_fc = log2fc_series.abs().median()
    if median_abs_fc < EFFECT_THRESHOLD:
        tier = "FAIL"
        checks.append({
            "check": "effect_magnitude",
            "result": "FAIL",
            "details": f"Median |log2FC| = {median_abs_fc:.6f} < {EFFECT_THRESHOLD} — no predicted effect"
        })
    else:
        checks.append({
            "check": "effect_magnitude",
            "result": "PASS",
            "details": f"Median |log2FC| = {median_abs_fc:.4f}"
        })

    # Check 3: Self-knockdown
    if perturbed_gene in log2fc_series.index:
        self_fc = log2fc_series[perturbed_gene]
        if self_fc >= 0:
            if tier != "FAIL":
                tier = "WARN"
            checks.append({
                "check": "self_knockdown",
                "result": "WARN",
                "details": f"{perturbed_gene} log2FC = {self_fc:.4f} (expected < 0). "
                           "May indicate model didn't encode self-suppression, or gene has autoregulatory feedback."
            })
        else:
            checks.append({
                "check": "self_knockdown",
                "result": "PASS",
                "details": f"{perturbed_gene} log2FC = {self_fc:.4f} (correctly downregulated)"
            })
    else:
        checks.append({
            "check": "self_knockdown",
            "result": "SKIP",
            "details": f"{perturbed_gene} not in model gene set"
        })

    # Check 4: Known targets
    all_targets = _merge_targets(known_targets, collectri_targets)
    if all_targets:
        activated = all_targets.get("activates", [])
        repressed = all_targets.get("represses", [])
        n_correct = 0
        n_total = 0
        target_details = []

        for gene in activated:
            if gene in log2fc_series.index:
                fc = log2fc_series[gene]
                n_total += 1
                correct = fc < 0  # knockdown of activator → target should go down
                if correct:
                    n_correct += 1
                target_details.append(f"{gene}: {fc:+.3f} ({'correct' if correct else 'unexpected'})")

        for gene in repressed:
            if gene in log2fc_series.index:
                fc = log2fc_series[gene]
                n_total += 1
                correct = fc > 0  # knockdown of repressor → target should go up
                if correct:
                    n_correct += 1
                target_details.append(f"{gene}: {fc:+.3f} ({'correct' if correct else 'unexpected'})")

        if n_total > 0:
            frac_correct = n_correct / n_total
            if frac_correct < 0.25:
                if tier != "FAIL":
                    tier = "WARN"
                result = "WARN"
            else:
                result = "PASS"
            checks.append({
                "check": "known_targets",
                "result": result,
                "details": f"{n_correct}/{n_total} ({frac_correct:.0%}) targets responded as expected. "
                           + "; ".join(target_details[:10])
            })
        else:
            checks.append({
                "check": "known_targets",
                "result": "SKIP",
                "details": "No known targets found in gene set"
            })
    else:
        checks.append({
            "check": "known_targets",
            "result": "SKIP",
            "details": "No known targets provided"
        })

    return {"tier": tier, "checks": checks}


def _merge_targets(expert_targets, collectri_targets):
    """Merge targets from expert YAML and CollecTRI."""
    merged = {"activates": [], "represses": []}
    if expert_targets:
        merged["activates"].extend(expert_targets.get("activates", []))
        merged["represses"].extend(expert_targets.get("represses", []))
    if collectri_targets:
        merged["activates"].extend(collectri_targets.get("activates", []))
        merged["represses"].extend(collectri_targets.get("represses", []))
    # Deduplicate
    merged["activates"] = list(set(merged["activates"]))
    merged["represses"] = list(set(merged["represses"]))
    if not merged["activates"] and not merged["represses"]:
        return None
    return merged
