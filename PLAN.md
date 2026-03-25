# AGI-Guided Differentiation Agent — Concrete Plan

## Vision

A **multi-agent conversational system** that guides researchers through cell differentiation
strategy design. Given a natural language query (e.g., *"help me differentiate iPSCs into
pancreatic beta cells"*), the system mines literature, discovers and downloads relevant
scRNA-seq data from public repositories, runs in-silico TF perturbation, performs pathway
analysis, and delivers a concrete differentiation protocol — keeping the human in the loop
at every major decision point.

---

## High-Level Multi-Agent Architecture

```
╔══════════════════════════════════════════════════════════════════════════════════╗
║                          USER (researcher)                                       ║
║   "Differentiate iPSCs into intestinal enterocytes using transcription factors"  ║
╚══════════════════════╦═══════════════════════════════════════════════════════════╝
                       │  natural language query
                       ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│                        ORCHESTRATOR AGENT                                        │
│   • Parses query → organ / target cell type / TFs of interest / constraints      │
│   • Maintains conversation state and human-in-the-loop checkpoints               │
│   • Routes to sub-agents, collects results, synthesizes final output             │
│   • Asks clarifying questions when query is ambiguous                            │
└───────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────────────┘
        │          │          │          │          │          │
        ▼          ▼          ▼          ▼          ▼          ▼
   [Agent 1]  [Agent 2]  [Agent 3]  [Agent 4]  [Agent 5]  [Agent 6]
  Literature   TF Intel  Dataset   Data Acq.  Preprocess  Perturbation
    Mining    Agent      Discovery  Agent      Agent       + Analysis
   Agent                Agent                             Agent
```

---

## Agent Definitions

### Agent 0 — Orchestrator

**Role**: The brain. Interprets user intent, sequences agent calls, presents results,
asks for human confirmation before each major step.

**Inputs**: Natural language query
**Outputs**: Final protocol report (narrative + table)
**Tools**: All sub-agents, conversation memory

**Human-in-the-loop checkpoints** (HITL — marked 🛑 throughout):

```
🛑 HITL-1  After query parsing     → confirm: target cell type, organ, TFs to test
🛑 HITL-2  After literature mining → confirm: TF shortlist before lookup
🛑 HITL-3  After dataset discovery → user selects which dataset(s) to download
🛑 HITL-4  After preprocessing     → user confirms cell type annotations
🛑 HITL-5  After SLURM submission  → user triggers analysis after GPU job completes
🛑 HITL-6  Before final report     → user can re-run with different TFs or datasets
```

---

### Agent 1 — Literature Mining Agent

**Role**: Find relevant papers and extract candidate transcription factors.

**Tools**:
- HIRN Literature Retrieve skill (existing) — for islet/pancreas context
- PubMed eutils API (`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`) — general search
- PubMed abstract fetcher — extract TF mentions from abstracts

**Flow**:
```
User query (organ + cell type)
        │
        ├──► HIRN skill search (if pancreas/islet relevant)
        │
        ├──► PubMed search: "{target cell type} transcription factor differentiation"
        │    PubMed search: "{organ} scRNA-seq cell identity"
        │
        ├──► Fetch top-N abstracts
        │
        └──► Extract TF candidates (NER on gene symbols + context)
             → ranked list: [CDX2, GATA4, HNF4A, VIL1, ...] + supporting PMIDs

                    🛑 HITL-2: present TF candidates to user → user confirms/edits list
```

**Output**: `{"tfs": ["CDX2", "GATA4", "HNF4A"], "evidence": [{pmid, title, snippet}]}`

---

### Agent 2 — TF Intelligence Agent

**Role**: For each confirmed TF, retrieve known target genes and biological context.

**Tools**:
- ChEA3 REST API (`https://maayanlab.cloud/chea3/api/enrich/`) — no auth, free
- DoRothEA via decoupleR Python API (confidence levels A/B/C)
- Gene Ontology API — validate TF biological role

**Flow**:
```
Confirmed TF list [CDX2, GATA4, HNF4A]
        │
        ├──► ChEA3: for each TF → ranked target gene list (7 evidence libraries)
        │
        ├──► Filter: keep targets with evidence in ≥2 libraries
        │
        └──► Cross-reference: TF targets ∩ genes differentially expressed
             in target cell type (from literature)
             → final gene list for perturbation input
```

**Output**: `{"CDX2": ["LYZ", "MUC2", "FABP1", ...], "GATA4": [...], ...}`

---

### Agent 3 — Dataset Discovery Agent

**Role**: Dynamically search public repositories for scRNA-seq datasets matching the
user's organ/tissue/developmental stage. No pre-existing catalog — always discovered fresh.

**Tools**:
- **GEO API** (`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi`)
  - Query: `"{organ}" AND "single cell RNA" AND "Homo sapiens"[Organism]`
  - Fetch GSE metadata: title, organism, n_samples, platform, supplementary files
- **ArrayExpress API** (`https://www.ebi.ac.uk/biostudies/api/v1/search`)
  - Query by tissue keyword, filter for RNA-seq, check for processed data availability
- **HuBMAP API** (`https://search.api.hubmapconsortium.org/v3/search`)
  - For adult human tissue atlases

**Flow**:
```
Organ / target cell type
        │
        ├──► GEO esearch: top-20 GSE accessions ranked by:
        │    - recency (newer = better protocols)
        │    - n_cells (larger = more robust clustering)
        │    - has processed h5ad/matrix in supplementary files
        │
        ├──► ArrayExpress search: top-10 E-MTAB accessions
        │
        ├──► Format as comparison table:
        │    Accession | Title | N_cells | Stage | Platform | Has_h5ad
        │    GSE158702 | Intestinal develop... | 50k | Fetal 8-20 PCW | 10x | Yes
        │    E-MTAB-9543 | Human intestine atlas | 400k | Adult | 10x | Yes
        │    ...
        │
        └──► 🛑 HITL-3: present table → user picks 1-2 datasets

```

**Output**: `{"selected": ["GSE158702"], "metadata": {...}}`

**Key design**: Results are never cached as a static file. Every run re-queries the
APIs so the agent always finds the most current public data.

---

### Agent 4 — Data Acquisition Agent

**Role**: Download selected dataset(s) and build a local h5ad file.

**Tools**: `requests`, `urllib`, `gzip`, `scipy.io`, `anndata`, `scanpy`

**Three-tier download strategy** (tried in order):

```
Tier 1 — Direct h5ad (preferred)
  NCBI eutils → list supplementary files for GSE accession
  → find *.h5ad or *_adata.h5 → wget → done ✓

Tier 2 — Count matrix reconstruction
  GEO SOFT file → parse sample metadata
  → download matrix.mtx.gz + barcodes.tsv.gz + features.tsv.gz
  → scipy.io.mmread → anndata.AnnData → save as h5ad ✓

Tier 3 — ArrayExpress processed files
  EBI BioStudies API → find processed_data/ directory
  → download tab-separated count matrix → build h5ad ✓

Tier 4 — SKIP (no raw FASTQ processing)
  If only SRA FASTQs available → inform user, skip this dataset
```

**Output**: `data/downloaded/{accession}.h5ad`

---

### Agent 5 — Preprocessing Agent

**Role**: Standard scanpy preprocessing + automated cell type annotation.

**Environment**: `state_env` (scanpy 1.11.5 already installed)

**Pipeline**:
```
raw h5ad
  │
  ├── QC filtering
  │   - min_genes=200, max_genes=6000
  │   - mito_fraction < 20%
  │   - doublet removal (scrublet)
  │
  ├── Normalization
  │   - normalize_total(target_sum=1e4)
  │   - log1p
  │
  ├── Feature selection
  │   - highly_variable_genes(n_top_genes=3000)
  │
  ├── Dimensionality reduction
  │   - PCA(n_comps=50) → neighbors → UMAP
  │
  ├── Clustering
  │   - Leiden (resolution=0.5 default, tunable)
  │
  └── Cell type annotation (marker scoring, no celltypist needed)
      - scanpy.tl.score_genes() per known cell type marker set
      - Marker sets: fetched dynamically from CellMarker 2.0 API
        for the user's organ of interest
      - Each cluster → top-scoring cell type label
      - Outputs: UMAP plot + cluster→celltype table

              🛑 HITL-4: show UMAP + annotation table → user corrects any labels
```

**Output**: `data/preprocessed/{accession}_annotated.h5ad` + `umap.png`

---

### Agent 6 — Perturbation + Analysis Agent

**Role**: Run in-silico perturbation and GSEA; synthesize the final protocol.

**Sub-steps**:

```
Step A — Preprocess for inference (existing script)
  preprocess_for_inference.py
  → inference_template.h5ad (18,080 genes, stratified cells, all TF conditions)

Step B — SLURM job submission
  → generate run_inference.sh with correct paths
  → print sbatch command for user
  🛑 HITL-5: user runs sbatch, waits 8-20hrs, then resumes conversation

Step C — GSEA (existing script, runs after GPU job)
  run_gsea.py → log2FC + volcano plots + pathway enrichment per TF

Step D — Synthesis
  For each TF:
    - Top upregulated pathways (NES > 0, FDR < 0.25)
    - Top target genes activated
    - Literature support (from Agent 1)
  → rank TFs by "differentiation potential score":
      score = NES(target pathway) × literature_evidence_count × ChEA3_confidence

Step E — Report generation (both formats)
  Narrative: paragraph summary per TF with citations
  Table:     | TF | Top Pathway | NES | Key Targets | Recommended Timing |
```

**Output**: `results/protocol_report.md` + `results/protocol_table.csv`

---

## Complete Data Flow

```
                         User Query
                             │
                    ┌────────▼────────┐
                    │  Orchestrator   │
                    │    Agent 0      │
                    └───┬────┬────┬───┘
                        │    │    │
              ┌─────────┘    │    └──────────┐
              ▼              ▼               ▼
         ┌─────────┐   ┌─────────┐    ┌──────────┐
         │Agent 1  │   │Agent 2  │    │ Agent 3  │
         │Lit Mine │   │TF Intel │    │ Dataset  │
         │         │   │         │    │Discovery │
         └────┬────┘   └────┬────┘    └────┬─────┘
              │             │              │
              │    TF list  │  targets     │  accession
              └──────┬──────┘      ┌───────┘
                     │             │
                 🛑 HITL-2    🛑 HITL-3
                     │             │
                     │        ┌────▼─────┐
                     │        │ Agent 4  │
                     │        │Data Acq. │
                     │        └────┬─────┘
                     │             │ raw h5ad
                     │        ┌────▼─────┐
                     │        │ Agent 5  │
                     │        │Preprocess│
                     │        └────┬─────┘
                     │             │ annotated h5ad
                     │         🛑 HITL-4
                     │             │
                     └──────┬──────┘
                            │ TFs + annotated h5ad
                       ┌────▼──────┐
                       │  Agent 6  │
                       │Perturbation│
                       │+ Analysis │
                       └────┬──────┘
                            │
                        🛑 HITL-5 (SLURM wait)
                            │
                       ┌────▼──────┐
                       │  REPORT   │
                       │ narrative │
                       │  + table  │
                       └───────────┘
```

---

## Directory Structure (to build)

```
AGI-guided_differentiation/
├── PLAN.md                          ← this file
├── CLAUDE.md                        ← project instructions
│
├── agents/                          ← agent definitions + scripts
│   ├── orchestrator/
│   │   └── AGENT.md                 ← orchestrator prompt + routing logic
│   ├── literature_mining/
│   │   ├── AGENT.md
│   │   └── scripts/
│   │       ├── pubmed_search.py     ← NCBI eutils search + abstract fetch
│   │       └── extract_tfs.py       ← TF mention extraction from abstracts
│   ├── tf_intelligence/
│   │   ├── AGENT.md
│   │   └── scripts/
│   │       └── chea3_lookup.py      ← ChEA3 REST API query
│   ├── dataset_discovery/
│   │   ├── AGENT.md
│   │   └── scripts/
│   │       ├── search_geo.py        ← NCBI eutils GEO search
│   │       ├── search_ebi.py        ← ArrayExpress/BioStudies API
│   │       └── rank_datasets.py     ← score by recency, size, data availability
│   ├── data_acquisition/
│   │   ├── AGENT.md
│   │   └── scripts/
│   │       ├── download_geo.py      ← Tier 1/2: h5ad or count matrix from GEO
│   │       └── download_ebi.py      ← Tier 3: ArrayExpress processed files
│   ├── preprocessing/
│   │   ├── AGENT.md
│   │   └── scripts/
│   │       ├── preprocess.py        ← scanpy QC/norm/cluster pipeline
│   │       ├── annotate.py          ← CellMarker API + score_genes annotation
│   │       └── cellmarker_api.py    ← fetch marker sets from CellMarker 2.0
│   └── perturbation_analysis/
│       ├── AGENT.md
│       └── scripts/
│           └── synthesize.py        ← combine GSEA results → protocol report
│
├── insilico_perturbation_pipeline/  ← existing (unchanged)
├── GSEA_skill/                      ← existing (unchanged)
├── hirn_publication_retrieval/      ← existing (unchanged)
│
├── data/
│   ├── downloaded/                  ← raw h5ads from GEO/EBI
│   ├── preprocessed/                ← annotated h5ads
│   └── perturbation/                ← STATE model outputs
│
└── results/
    ├── protocol_report.md           ← narrative output
    └── protocol_table.csv           ← structured output
```

---

## Build Order

### Phase 1 — Fast path (no GPU, 1-2 weeks)

| # | Task | Agent | Notes |
|---|------|-------|-------|
| 1 | `pubmed_search.py` | Lit Mining | NCBI eutils, free |
| 2 | `extract_tfs.py` | Lit Mining | regex + gene symbol dict |
| 3 | `chea3_lookup.py` | TF Intel | POST to maayanlab.cloud |
| 4 | `search_geo.py` | Dataset Discovery | eutils esearch + efetch |
| 5 | `search_ebi.py` | Dataset Discovery | BioStudies REST API |
| 6 | `rank_datasets.py` | Dataset Discovery | score + format table for HITL |
| 7 | `AGENT.md` files | All | prompts + routing instructions |
| 8 | Orchestrator AGENT.md | Orchestrator | HITL flow, checkpoint logic |

### Phase 2 — Data + compute (1-2 weeks)

| # | Task | Agent | Notes |
|---|------|-------|-------|
| 9 | `download_geo.py` | Data Acq. | Tier 1+2 strategy |
| 10 | `download_ebi.py` | Data Acq. | Tier 3 strategy |
| 11 | `preprocess.py` | Preprocess | scanpy pipeline |
| 12 | `cellmarker_api.py` | Preprocess | fetch marker sets dynamically |
| 13 | `annotate.py` | Preprocess | score_genes + HITL annotation |

### Phase 3 — Synthesis (1 week)

| # | Task | Agent | Notes |
|---|------|-------|-------|
| 14 | SLURM script generator | Perturbation | parameterized from HITL selections |
| 15 | `synthesize.py` | Perturbation | GSEA → score TFs → report |
| 16 | Report templates | Perturbation | markdown + CSV formatters |
| 17 | End-to-end test | All | example: iPSC → beta cell |

---

## Key Design Principles

1. **No static catalog**: Dataset discovery always runs live against GEO/EBI APIs.
   The agent's knowledge of available data is always current.

2. **Human always in the loop**: The system never makes consequential decisions
   autonomously (which TFs to test, which dataset to download, how to annotate
   cell types). It presents options with evidence and waits for confirmation.

3. **Graceful degradation**: If a dataset has no processed h5ad (Tier 4), the agent
   skips it and presents the next best option — never blocks on unavailable data.

4. **Async compute**: GPU inference is handled via SLURM. The agent generates the
   submission script, the user runs it, and resumes the conversation with results.

5. **Reproducibility**: Every run logs the exact accessions, TF lists, parameters,
   and script versions used, so the protocol is fully reproducible.

6. **Modular agents**: Each agent can be invoked independently. A user who already
   has preprocessed data can start at Agent 6. A user who already knows their TFs
   can skip Agents 1-2.

---

## Example Session (end-to-end)

```
User:   "I want to differentiate human iPSCs into intestinal enterocytes.
         I'm interested in CDX2 and maybe a few other TFs."

Agent:  [Agent 1] Searching PubMed for intestinal enterocyte differentiation TFs...
        Found 847 papers. Top TF candidates: CDX2, GATA4, HNF4A, VIL1, KLF4
        Supporting evidence: [PMID:31348891, PMID:33406409, ...]

        🛑 HITL-2: "Shall I proceed with CDX2, GATA4, HNF4A? Add/remove any?"

User:   "Yes, also add ISX."

Agent:  [Agent 2] Fetching ChEA3 targets for CDX2, GATA4, HNF4A, ISX...
        CDX2 targets (n=312): MUC2, LYZ, FABP1, CLCA4...
        GATA4 targets (n=287): HAND1, NKX2-5, MYH6...

        [Agent 3] Searching GEO + ArrayExpress for intestinal scRNA-seq data...

        Found 8 relevant datasets:
        ┌────────────┬──────────────────────────────┬────────┬──────────┬────────┐
        │ Accession  │ Title                        │ Cells  │ Stage    │ h5ad?  │
        ├────────────┼──────────────────────────────┼────────┼──────────┼────────┤
        │ GSE158702  │ Spatiotemporal intestinal dev│ 50,517 │ Fetal    │ Yes    │
        │ E-MTAB-9543│ Cells of human intestinal... │ 428k   │ Adult    │ Yes    │
        │ GSE125970  │ Single-cell transcriptome... │ 9,261  │ Fetal    │ No     │
        └────────────┴──────────────────────────────┴────────┴──────────┴────────┘

        🛑 HITL-3: "Which dataset would you like to use?"

User:   "GSE158702 — fetal makes more sense for iPSC comparison."

Agent:  [Agent 4] Downloading GSE158702... found h5ad in supplementary files.
        Saved to data/downloaded/GSE158702.h5ad (2.1 GB)

        [Agent 5] Preprocessing...
        After QC: 47,382 cells retained. Running UMAP + Leiden clustering...
        Annotating 18 clusters using CellMarker intestine marker sets:
        Cluster 0 → Enterocyte (score 0.82)
        Cluster 1 → Goblet cell (score 0.79)
        Cluster 3 → Enteroendocrine (score 0.71)
        ...

        🛑 HITL-4: [shows UMAP] "Do these annotations look correct?"

User:   "Cluster 5 should be Paneth cells, not stem cells."

Agent:  Updated. Saving annotated h5ad.

        [Agent 6] Generating inference template with 4 TFs + non-targeting control...
        Cells per condition: ~9,476 (stratified by cell type)
        Memory estimate: 9,476 × 5 × 18,080 × 4 bytes ≈ 3.4 GB

        Ready to run inference. Here is your SLURM script:
        [run_inference_GSE158702_CDX2_GATA4_HNF4A_ISX.sh]

        🛑 HITL-5: "Run: sbatch run_inference_GSE158702_CDX2_GATA4_HNF4A_ISX.sh
                    Then come back when the job completes (~12 hrs)."

User:   [12 hours later] "Job is done."

Agent:  [Agent 6] Running GSEA on perturbation predictions...

        ╔══════════════════════════════════════════════════════════════════╗
        ║              DIFFERENTIATION PROTOCOL RECOMMENDATION            ║
        ╠════════╦════════════════════════════╦══════╦════════════════════╣
        ║ TF     ║ Top Activated Pathway      ║  NES ║ Key Target Genes   ║
        ╠════════╬════════════════════════════╬══════╬════════════════════╣
        ║ CDX2   ║ Intestinal absorptive cells║ 2.41 ║ MUC2, LYZ, FABP1  ║
        ║ GATA4  ║ Epithelial differentiation ║ 1.98 ║ HNF1A, ONECUT1    ║
        ║ HNF4A  ║ Lipid metabolism           ║ 1.87 ║ APOB, APOA4, CUBN ║
        ║ ISX    ║ Carotenoid absorption      ║ 1.43 ║ BCO1, SCARB1      ║
        ╚════════╩════════════════════════════╩══════╩════════════════════╝

        Recommended protocol: Overexpress CDX2 (primary driver) + HNF4A
        (metabolic maturation) in combination. GATA4 at earlier stage.
        Full narrative report saved to results/protocol_report.md
```

---

## Dependencies Summary

| Package | Used by | Status in state_env |
|---------|---------|---------------------|
| `scanpy` | Agent 5, 6 | ✓ installed (1.11.5) |
| `anndata` | Agent 4, 5 | ✓ installed |
| `requests` | Agents 1,2,3,4 | ✓ stdlib-adjacent |
| `scipy` | Agent 4 (mtx read) | ✓ installed |
| `gseapy` | Agent 6 (GSEA) | ✓ (existing skill) |
| `matplotlib` | Agent 5, 6 | ✓ installed |
| `pandas` | All | ✓ installed |
| `scrublet` | Agent 5 (doublets) | ⚠ may need pip install |
| `celltypist` | — | ✗ not used (replaced by score_genes) |
| `GEOparse` | — | ✗ not used (direct API instead) |
| `ffq` | — | ✗ not used (direct API instead) |

No new heavy dependencies required. All data access via standard REST APIs + `requests`.
