# Synthesis Pipeline Plan: Perturbation Interpretation & Biological Analysis

## Goal

Given in-silico perturbation predictions (predicted gene expression after gene knockdown), automatically generate a structured biological interpretation report that integrates multiple lines of evidence — the way an experienced biologist would.

## Inputs

1. **Prediction h5ad** — output of `state tx infer` (cells × 18,080 genes, with `target_gene` and `cell_type` in `.obs`)
2. **Expert knowledge YAML** (optional) — custom gene sets, known targets, hypotheses, focus areas
3. **Dataset context** — which tissue/organ, which cell type

## Output

A **structured markdown report** per perturbation per cell type, containing:
- Sanity check results (pass/fail)
- Top DE genes (true log2FC) with known-target annotations
- Pathway enrichment (GSEA + ORA) results
- TF activity changes (decoupleR)
- Pathway activity scores (PROGENy)
- Cell state scores (proliferation, apoptosis, differentiation, stress)
- Literature references for top hits
- Evidence convergence summary

Plus a **cross-perturbation summary table** (CSV) with one row per perturbation.

---

## Architecture

```
analyze_perturbation.py --predictions pred.h5ad --config expert.yaml --output-dir results/

Phase 1: Sanity Checks
  └─► sanity_report (pass/fail per check)

Phase 2: Differential Expression
  └─► de_results.csv (ranked gene list with log2fc, log_diff, means)

Phase 3: Functional Enrichment
  ├─► 3a. GSEA (gseapy preranked with log2fc)
  ├─► 3b. ORA (decoupleR dc.mt.ora — top DE genes vs GO/Reactome)
  ├─► 3c. TF activity (decoupleR dc.mt.ulm + collectri network)
  ├─► 3d. Pathway activity (decoupleR PROGENy — 14 signaling pathways)
  └─► 3e. Custom gene set scoring (scanpy score_genes from expert YAML)

Phase 5: Phenotypic Interpretation
  ├─► Cell state scores (proliferation, apoptosis, differentiation, stress)
  └─► Cross-cell-type comparison (same perturbation, different cell types)

Phase 6: Literature Integration
  ├─► PubMed search per top gene + pathway
  ├─► HIRN search (for islet-relevant datasets)
  └─► Gene annotation (Ensembl/UniProt via APIs)

Phase 7: Multi-Evidence Synthesis
  ├─► Convergence assessment across all phases
  ├─► Structured report generation (markdown)
  └─► Cross-perturbation summary table (CSV)
```

Phase 4 (PPI network analysis) is deferred to v2.

---

## Phase 1: Sanity Checks

**Purpose:** Catch prediction failures before wasting analysis time.

**Checks:**
1. **Housekeeping stability**: Are canonical housekeeping genes (ACTB, GAPDH, RPL13A, UBC, B2M) stable (|log2FC| < 0.5)?
2. **Effect magnitude**: Is the median absolute log2FC > some minimum threshold (e.g., 0.001)?
   - If all genes have near-zero changes, the perturbation had no predicted effect
3. **Self-knockdown**: Is the perturbed gene itself downregulated? (log2FC < 0)
   - Note: not always expected — depends on autoregulatory feedback, whether the gene is expressed in this cell type, and whether the model learned self-suppression
4. **Known targets** (if provided in expert YAML): Are known direct targets of the perturbed gene changing in the expected direction?
   - Uses `expert_knowledge.yaml → known_targets` mapping
   - Also checks CollecTRI database (via decoupleR) for known TF→target relationships

**Tiered QC outcome:**

| Tier | Condition | Action |
|------|-----------|--------|
| **FAIL** | Multiple housekeeping genes unstable (\|log2FC\| > 0.5) | Block synthesis — abbreviated report only (sanity results + raw DE list), header: "QC FAILED — predictions unreliable" |
| **FAIL** | No effect (median \|log2FC\| < 0.001 across all genes) | Block synthesis — same abbreviated report |
| **WARN** | Self-knockdown not detected (gene expressed but doesn't decrease) | Full report generated with prominent warning banner |
| **WARN** | Known targets don't respond (< 25% move in expected direction) | Full report with warning — could be novel biology |
| **PASS** | All checks pass | Full report, no caveats |

QC-failed perturbations get no synthesis narrative, no literature integration, no convergence assessment — only raw data. This prevents polished biological stories from being generated for mechanically broken predictions.

An `--force-full-report` flag allows override for exploratory analysis, but all downstream artifacts are watermarked as QC-failed.

**Output:** `sanity_checks.json` with tier (FAIL/WARN/PASS), per-check results, and details.

---

## Phase 2: Differential Expression

**Purpose:** Extract the ranked gene list that feeds all downstream analyses, with proper uncertainty estimates.

**Key design:** The prediction template creates matched copies of each original cell across non-targeting and every perturbation (linked by `original_barcode`). This paired structure should be exploited, not collapsed into a single mean (Squair et al., Nat Commun 2021).

**Method:**
1. **Paired per-cell deltas:** For each cell (matched by `original_barcode`), compute `delta = expm1(pert) - expm1(ctrl)` in natural scale
2. **Aggregate:** Mean delta per gene across all cells
3. **True log2FC:** `log2(mean_nat_pert + 1) - log2(mean_nat_ctrl + 1)` from natural-scale means
4. **Uncertainty (bootstrap CI):** Resample cells with replacement (1000 iterations), compute log2FC per resample, report 95% CI per gene
5. **Statistical test:** For each gene, test whether the paired deltas are significantly different from zero (one-sample t-test or Wilcoxon on deltas). Apply BH FDR correction.

**"Significant" means tested:** A gene is only labeled "significant" if it passes both:
- |log2FC| > threshold (e.g., 0.5) — biological effect size
- FDR-adjusted p-value < 0.05 — statistical significance from paired test

**Output:** `de_results.csv` with columns: `gene, log2fc, log2fc_ci_low, log2fc_ci_high, pval, fdr, mean_pert, mean_ctrl, is_significant`

**Additional annotations:**
- Flag genes that are known targets of the perturbed gene (from CollecTRI or expert YAML)
- Flag marker genes for the cell type being analyzed
- Flag housekeeping genes

**Performance note:** With 500 cells × 1000 bootstrap iterations × 18,080 genes, this is ~36 GB of computation. Can be optimized by vectorizing the bootstrap (resample indices, matrix multiply) rather than looping. For the demo, a simpler approach (t-test on paired deltas, no bootstrap) is acceptable.

---

## Phase 3: Functional Enrichment

### 3a. GSEA (Pre-ranked)

**Tool:** gseapy (already in repo)
**Ranking:** True log2FC (natural-scale)
**Gene set libraries:**
- MSigDB Hallmark (50 pathway signatures) — always run
- GO Biological Process 2023 — always run
- KEGG 2021 Human — always run
- TF Perturbations Followed by Expression — if available (directly relevant)

**Output:** GSEA results with NES, p-value, FDR. Only FDR < 0.25 included in summary.

### 3b. Over-Representation Analysis (ORA)

**Tool:** decoupleR `dc.mt.ora`
**Input:** Top DE genes (|log2FC| > threshold, e.g., 0.5)
**Gene sets:** MSigDB Hallmark, GO BP, Reactome
**Purpose:** Complementary to GSEA — tests whether top DE genes are enriched in specific pathways. Biologists often run both.

**Output:** ORA results table with odds ratio, p-value, FDR.

### 3c. TF Activity Inference

**Tool:** decoupleR `dc.mt.ulm` (Univariate Linear Model)
**Network:** CollecTRI (curated TF→target database via `dc.op.collectri()`)
**Purpose:** Given the expression changes, which TFs are functionally activated or repressed? This is more informative than raw DE because it infers regulatory state from target gene behavior.

**Method:**
1. Fetch CollecTRI network (TF→target with sign: activation/repression)
2. Run ULM on the DE results (log2FC as input, TF activity as output)
3. Rank TFs by activity score

**Output:** `tf_activity.csv` — TF name, activity score, p-value. Shows which TFs are functionally affected by the perturbation.

**Biological value:** If we knock down TF_X and TF_Y's activity drops, that suggests TF_X → TF_Y regulatory relationship. If TF_Y's known function matches an enriched pathway, that's convergent evidence.

### 3d. Pathway Activity Scoring (PROGENy)

**Tool:** decoupleR with PROGENy model (`dc.op.progeny()`)
**Purpose:** Score 14 cancer-relevant signaling pathways (EGFR, MAPK, PI3K, TNFa, NFkB, Hypoxia, JAK-STAT, TGFb, Trail, p53, Androgen, Estrogen, VEGF, WNT) based on expression changes of their footprint genes.
**Advantage over GSEA:** PROGENy uses downstream transcriptional footprints (not pathway membership), so it captures pathway *activity* rather than just gene overlap.

**Method:**
1. Fetch PROGENy model via `dc.op.progeny()`
2. Run MLM (Multivariate Linear Model) on DE results
3. Report pathway activity scores

**Output:** `pathway_activity.csv` — pathway name, activity score, p-value.

### 3e. Custom Gene Set Scoring

**Tool:** scanpy `sc.tl.score_genes()`
**Input:** Custom gene sets from `expert_knowledge.yaml → custom_gene_sets`
**Purpose:** Let biologists test their own hypotheses — e.g., "does this perturbation affect my stem cell markers?"

**Method:**
1. For each custom gene set, compute score in perturbed vs. control cells
2. Report delta score (perturbed - control)

**Output:** `custom_scores.csv` — gene set name, control score, perturbed score, delta.

---

## Phase 5: Phenotypic Interpretation

**Purpose:** Translate molecular changes into cellular phenotype predictions.

### Cell State Scoring

**Tool:** scanpy `sc.tl.score_genes()` with curated marker gene sets

**Built-in marker sets** (hardcoded, well-established):
- **Proliferation:** MKI67, TOP2A, PCNA, MCM2, MCM6, CDK1, CCNB1, CCNA2
- **Pro-apoptotic:** BAX, BAK1, BID, CASP3, CASP9, CYCS, TP53, PMAIP1
- **Anti-apoptotic:** BCL2, BCL2L1, MCL1, BIRC5, XIAP
- **Stress response (UPR):** HSPA1A, HSPA1B, HSP90AA1, DDIT3, ATF4, XBP1
- **Stemness** (context-dependent, may override from expert YAML)

**Important:** Pro-apoptotic and anti-apoptotic programs are scored separately, not combined. `score_genes()` is unsigned — mixing opposing genes (e.g., BAX + BCL2) would produce a meaningless composite. The report presents both scores and interprets the balance: e.g., "pro-apoptotic up, anti-apoptotic down → net pro-apoptotic shift."

**Method:**
1. Score each marker set in perturbed vs. control cells (using paired per-cell structure from Phase 2)
2. Report delta for each state with bootstrap CI

**Output:** `cell_state_scores.csv` — state, control_score, perturbed_score, delta, ci_low, ci_high, interpretation.

### Cross-Cell-Type Comparison

If predictions exist for the same perturbation across multiple cell types (which they do — we run per cell type):
- Compare log2FC profiles across cell types (correlation)
- Identify cell-type-specific vs. universal effects
- Flag perturbations that have opposite effects in different cell types

**Output:** `cross_celltype_comparison.csv` — perturbation, cell_type_A, cell_type_B, correlation, n_shared_de_genes.

---

## Phase 6: Literature Integration

**Purpose:** Cross-reference computational predictions with published biology.

### 6a. PubMed Search

**Tool:** NCBI E-utilities API (via requests)
**Queries per perturbation** (properly parenthesized, field-tagged):
1. `({gene}[tiab]) AND ({cell_type}[tiab] OR {tissue}[tiab])` — direct evidence in this context
2. `({gene}[tiab]) AND ({top_pathway}[tiab])` — pathway context
3. `({gene}[tiab]) AND (knockdown[tiab] OR knockout[tiab] OR perturbation[tiab]) AND ({tissue}[tiab])` — perturbation studies in relevant tissue

All queries constrained to `"Homo sapiens"[Mesh]` when working with human datasets.

**Gene symbol disambiguation:** Use official HGNC symbol. For ambiguous symbols (e.g., ACE, REST), add the full gene name as OR term.

**Return:** Top 5 abstracts per query (title, PMID, snippet).

**Rate limiting:** 3 requests/second without API key, 10 with `NCBI_API_KEY`.

### 6b. HIRN Search (islet datasets only)

**Tool:** hirn_publication_retrieval scripts
**When:** Only for human islets dataset
**Queries:** Same as PubMed but searched against HIRN publication corpus
**Return:** Top 5 relevant text chunks with PMIDs.

### 6c. Gene Annotation

**Tool:** Ensembl REST API or mygene.info
**For:** Top 20 DE genes + the perturbed gene itself
**Return:** Gene description, known function, associated diseases, GO terms.

**Output:** `literature_results.json` — structured per-gene and per-pathway literature hits.

---

## Phase 7: Multi-Evidence Synthesis

**Purpose:** Integrate all evidence into a coherent biological narrative.

### Convergence Assessment

For each perturbation, evaluate evidence across 6 dimensions:

| Dimension | Source | Positive signal |
|-----------|--------|-----------------|
| Mechanistic consistency | Phase 1 + 2 | Known targets respond correctly |
| Functional consistency | Phase 3a-b | Enriched pathways match gene's known function |
| Regulatory consistency | Phase 3c | TF activity changes are mechanistically sensible |
| Pathway consistency | Phase 3d | PROGENy pathway scores align with GSEA |
| Phenotypic consistency | Phase 5 | Cell state shifts match expected biology |
| Literature consistency | Phase 6 | Published evidence supports predicted effects |

**No composite score.** Each dimension reported separately. The synthesis narrative (below) contextualizes them.

### Structured Report Generation

**Format:** Markdown, one report per (perturbation, cell_type).

**Template:**
```markdown
# Perturbation Report: {GENE} knockdown in {CELL_TYPE}

## QC Status: {PASS / WARN / FAIL}
{If WARN: banner explaining which checks triggered warnings}
{If FAIL: "Predictions flagged as unreliable. Only raw DE shown below."}

## Summary
{2-3 sentence LLM-generated summary integrating key findings}
{Only generated for PASS/WARN perturbations}

## Sanity Checks
{table: check, result, details}

## Differential Expression
- {N} genes significantly changed (|log2FC| > 0.5, FDR < 0.05 from paired test)
- Top upregulated: {gene1} (+{fc}, 95% CI [{lo}, {hi}]), ...
- Top downregulated: {gene1} ({fc}, 95% CI [{lo}, {hi}]), ...
- Known target validation: {N}/{M} known targets responded as expected

## Pathway Enrichment
### GSEA (pre-ranked, FDR < 0.25)
{table: pathway, NES, FDR}

### Over-Representation Analysis
{table: pathway, odds_ratio, FDR}

## Regulatory Analysis
### TF Activity Changes (decoupleR ULM + CollecTRI)
{table: TF, activity_score, p-value, interpretation}

### Pathway Activity (PROGENy)
{table: pathway, activity_score, p-value}

## Cell State Impact
{table: state, delta_score, CI, interpretation}
Note: pro-apoptotic and anti-apoptotic scored separately

## Literature Context
{Per top gene/pathway: title, PMID, key finding}

## Evidence Convergence
{Table showing 6 dimensions, each with signal strength}

## Interpretation
{LLM-generated paragraph synthesizing all evidence into a biological narrative.
What does this perturbation do to the cell? Is the effect consistent with
known biology? What is surprising or novel? What follow-up experiments
would test the predictions?}
```

### Cross-Perturbation Summary Table

One CSV summarizing all perturbations analyzed:

| perturbation | cell_type | n_de_genes | top_pathway | top_pathway_NES | top_tf_affected | proliferation_delta | apoptosis_delta | sanity_pass | n_literature_hits |
|---|---|---|---|---|---|---|---|---|---|

---

## Expert Knowledge YAML Format

```yaml
# expert_knowledge.yaml — optional, all fields optional

custom_gene_sets:
  gut_stem_markers: [LGR5, OLFM4, ASCL2, SOX9, SMOC2]
  paneth_markers: [LYZ, DEFA5, DEFA6, MMP7]
  proliferation: [MKI67, TOP2A, PCNA, MCM6]

known_targets:
  CDX2: {activates: [MUC2, LYZ, FABP1, VIL1], represses: []}
  HNF4A: {activates: [APOB, APOA4, CUBN], represses: []}

hypotheses:
  - "CDX2 knockdown should reduce enterocyte markers and increase stem markers"
  - "GATA4 should affect epithelial differentiation in gut but not islets"

context: |
  Adult gut dataset from healthy donors.
  Studying TFs that maintain intestinal epithelial cell identity.

focus_genes: [CDX2, HNF4A, GATA4, KLF4]
focus_pathways: ["Wnt signaling", "Notch signaling"]
```

---

## Implementation Plan

### File Structure

```
insilico_perturbation_pipeline/
├── scripts/
│   ├── analyze_perturbation.py      ← NEW: main entry point
│   ├── synthesis/
│   │   ├── __init__.py
│   │   ├── sanity_checks.py         ← Phase 1
│   │   ├── differential_expression.py ← Phase 2 (refactored from run_gsea.py)
│   │   ├── enrichment.py            ← Phase 3a-b (GSEA + ORA)
│   │   ├── regulatory.py            ← Phase 3c-d (TF activity + PROGENy)
│   │   ├── cell_state.py            ← Phase 5
│   │   ├── literature.py            ← Phase 6
│   │   └── report.py                ← Phase 7 (report generation)
│   └── ... (existing scripts)
├── configs/
│   ├── expert_knowledge_example.yaml ← NEW: example expert input
│   └── ... (existing configs)
```

### Dependencies (all in state_env)

| Package | Version | Use |
|---------|---------|-----|
| decoupler | 2.1.4 | TF activity, PROGENy, ORA |
| gseapy | 1.1.12 | Pre-ranked GSEA |
| scanpy | 1.11.5 | Gene scoring, DE |
| networkx | 3.5 | Future: PPI analysis |
| requests | 2.32.5 | PubMed API, Ensembl API |
| pyyaml | 6.0.3 | Config parsing |

### CLI Usage

```bash
# Analyze one prediction file (one cell type, one batch)
python scripts/analyze_perturbation.py \
    --predictions runs/adult_gut/stem_cells/batches/predictions_batch000.h5ad \
    --config configs/expert_knowledge_example.yaml \
    --output-dir results/adult_gut/stem_cells/ \
    --tissue "adult gut" \
    --cell-type "Stem cells"

# Analyze specific perturbation only
python scripts/analyze_perturbation.py \
    --predictions predictions.h5ad \
    --target-gene CDX2 \
    --output-dir results/

# Skip literature search (faster, offline)
python scripts/analyze_perturbation.py \
    --predictions predictions.h5ad \
    --no-literature \
    --output-dir results/
```

### Demo Plan

1. Pick: adult gut stem cells, batch000 predictions
2. Create example `expert_knowledge.yaml` with gut-relevant gene sets and known targets
3. Run `analyze_perturbation.py` for one perturbation (e.g., a well-known gut TF)
4. Inspect the generated report
5. Iterate on report format based on what's informative

---

## Design Decisions

### Why no composite score
Multiplying NES × literature_count × TF_confidence creates a single number that:
- Mixes incompatible quantities
- Biases toward well-studied genes
- Gives false quantitative confidence
Instead: present evidence dimensions separately, let the narrative synthesis contextualize.

### Why log2FC from natural-scale means for GSEA
Mean-of-logs ≠ log-of-means (Jensen's inequality). For sparse scRNA-seq data with high dropout, the difference is material (Booeshaghi & Pachter 2021). Natural-scale log2FC is the statistic biologists expect and interpret.

### Why decoupleR over manual implementation
decoupleR provides validated, published methods (ULM, VIPER, ORA, GSEA) with curated regulatory networks (CollecTRI, PROGENy, DoRothEA). Reimplementing these would introduce bugs and lack the validation.

### Why separate report files per phase
Each phase produces its own output file (CSV/JSON). The final report reads these to generate the synthesis. This means:
- Individual phases can be rerun independently
- Results are inspectable before synthesis
- The LLM synthesis step can be rerun without recomputing everything

### Why HITL checkpoint after Phase 3
After enrichment results are available, the biologist may want to:
- Add custom gene sets they hadn't thought of before
- Adjust focus areas based on what pathways came up
- Flag known artifacts or expected results to skip in the narrative
The system should support re-running Phase 7 with updated expert input without recomputing Phases 1-5.

---

## Future Work (v2)

- **Phase 4: PPI network analysis** — STRING API integration, hub gene identification, perturbation propagation visualization
- **Trajectory impact** — CellRank integration for differentiation trajectory analysis
- **Connectivity Map** — match predicted signatures to L1000 drug perturbation database
- **Cross-dataset synthesis** — compare same perturbation across adult gut, fetal gut, islets, hypothalamus
- **Automated batch analysis** — run `analyze_perturbation.py` across all 81 batches via SLURM
- **Interactive dashboard** — Plotly-based visualization of results across perturbations and cell types
