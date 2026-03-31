# AGI-Guided Differentiation

A multi-agent AI system that guides researchers through cell differentiation strategy design. Given a natural language query (e.g., *"help me differentiate iPSCs into pancreatic beta cells"*), the system mines literature, discovers public scRNA-seq datasets, runs in-silico transcription factor (TF) perturbation, performs pathway analysis, and delivers a concrete differentiation protocol — keeping the human in the loop at every major decision point.

## Overview

The system is built around a pipeline of specialized agents:

```
User Query (natural language)
        │
   Orchestrator Agent
        │
   ┌────┼────┬────────┬────────┬─────────────┐
   ▼    ▼    ▼        ▼        ▼             ▼
 Lit   TF  Dataset  Data    Preprocess  Perturbation
Mining Intel Discovery Acquisition           + Analysis
```

| Agent | Role |
|-------|------|
| **Orchestrator** | Parses query, sequences agents, manages human-in-the-loop checkpoints |
| **Literature Mining** | Searches PubMed for TF candidates relevant to the target cell type |
| **TF Intelligence** | Queries ChEA3 / DoRothEA for known TF target genes |
| **Dataset Discovery** | Searches GEO and ArrayExpress for relevant scRNA-seq datasets |
| **Data Acquisition** | Downloads selected datasets (h5ad, count matrix, or ArrayExpress) |
| **Preprocessing** | QC, normalization, clustering, automated cell type annotation |
| **Perturbation + Analysis** | Runs STATE model inference + GSEA, synthesizes protocol report |

### Human-in-the-Loop Checkpoints

The system pauses for user confirmation at six key decision points:

| # | When | What the user decides |
|---|------|-----------------------|
| HITL-1 | After query parsing | Confirm target cell type, organ, TFs to test |
| HITL-2 | After literature mining | Edit/approve TF candidate list |
| HITL-3 | After dataset discovery | Choose which dataset(s) to download |
| HITL-4 | After preprocessing | Correct cell type annotations on UMAP |
| HITL-5 | After SLURM submission | Trigger analysis after GPU job completes |
| HITL-6 | Before final report | Re-run with different TFs or datasets |

---

## Repository Structure

```
AGI-guided_differentiation/
├── PLAN.md                              # Full multi-agent architecture plan
├── CLAUDE.md                            # Instructions for Claude Code
│
├── insilico_perturbation_pipeline/      # Perturbation prediction pipeline
│   ├── scripts/
│   │   ├── run_full_pipeline.py         # Config-driven orchestrator (prepare/submit/status)
│   │   ├── prepare_celltype_run.py      # Subset h5ad + filter DEGs per cell type
│   │   ├── preprocess_for_inference.py  # Normalize, align genes, create template
│   │   ├── run_insilico_perturbation.py # Python wrapper for state tx infer
│   │   ├── run_inference.sh             # SLURM job template (single run)
│   │   ├── extract_to_csv.py            # Split predictions to CSV
│   │   ├── train_state_model.py         # Optional: train your own model
│   │   └── create_esm_pert_features.py  # Optional: generate ESM2 embeddings
│   ├── configs/
│   │   ├── collab_filtered_v1.yaml      # Demo: 4 datasets, 46 cell types
│   │   ├── starter.toml                 # Training config (subset data)
│   │   └── full_dataset.toml            # Training config (full dataset)
│   ├── models/                          # Bundled state_sm model + gene_names.csv
│   ├── data/                            # Place input h5ad here; outputs go here
│   └── references/models.md             # Registry of available STATE models
│
├── GSEA_skill/                          # Gene Set Enrichment Analysis skill
│   ├── run_gsea.py                      # Compute log2FC + run gseapy prerank
│   └── test_results/TBX5/              # Example GSEA output (TBX5 perturbation)
│
└── hirn_publication_retrieval/          # HIRN literature retrieval
```

---

## In-Silico Perturbation Pipeline

Predicts how gene perturbations would shift cell state (gene expression) using pre-trained [STATE](https://github.com/genentech/state) models.

### Requirements

```bash
source activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
```

### Quick Start: Config-Driven Batch Mode

Run perturbation predictions across all cell types in multiple datasets from a YAML config:

```bash
# Prepare all datasets x cell types + generate SLURM jobs
python insilico_perturbation_pipeline/scripts/run_full_pipeline.py prepare \
    --config insilico_perturbation_pipeline/configs/collab_filtered_v1.yaml

# Check completion status
python insilico_perturbation_pipeline/scripts/run_full_pipeline.py status \
    --config insilico_perturbation_pipeline/configs/collab_filtered_v1.yaml

# Submit SLURM jobs (rate-limited, skips completed batches)
python insilico_perturbation_pipeline/scripts/run_full_pipeline.py submit \
    --config insilico_perturbation_pipeline/configs/collab_filtered_v1.yaml
```

See `configs/collab_filtered_v1.yaml` for a working example (4 datasets, 46 cell types, 81 batches).

### Quick Start: Manual Single Run

For a single h5ad + gene list (e.g., collaborator data):

```bash
# Step 1: Preprocess
python insilico_perturbation_pipeline/scripts/preprocess_for_inference.py \
    --input data/example_fetal_heart.h5ad \
    --gene-list data/example_cardiac_tfs.csv \
    --output data/inference_template.h5ad \
    --cell-type-col celltype --n-cells 10000

# Step 2: Inference (SLURM or direct)
sbatch insilico_perturbation_pipeline/scripts/run_inference.sh

# Step 3: Export
python insilico_perturbation_pipeline/scripts/extract_to_csv.py \
    --predictions data/perturbation_predictions.h5ad \
    --output-dir data/perturbation_csvs
```

### Resource Requirements

| Step | RAM | GPU | Time |
|------|-----|-----|------|
| Preprocess | 256 GB | None | Minutes |
| `state_sm` inference | 256 GB | Any GPU | 8–20 hrs |
| Export CSV | Low | None | Minutes |

### Key Constraints

- Model expects exactly **18,080 genes** — preprocessing handles alignment automatically
- Input must be **log-normalized** (`normalize_total` + `log1p`); raw counts are auto-detected and normalized
- Matrices must be **dense** (not sparse); memory ≈ `n_cells × (n_perturbations + 1) × 18,080 × 4 bytes`
- Input h5ad must include cell type labels (specify column with `--cell-type-col`)

---

## GSEA Skill

Computes log2 fold change between a perturbed condition and non-targeting control, then runs gseapy prerank GSEA against any MSigDB gene set collection.

```bash
python GSEA_skill/run_gsea.py \
    --predictions perturbation_predictions.h5ad \
    --gene <TF_NAME> \
    --gene-set "MSigDB_Hallmark_2020" \
    --output-dir results/gsea/
```

Example output (TBX5 perturbation): `GSEA_skill/test_results/TBX5/`

---

## Example Session

```
User:   "I want to differentiate human iPSCs into intestinal enterocytes.
         I'm interested in CDX2 and maybe a few other TFs."

Agent:  [Literature Mining] Found TF candidates: CDX2, GATA4, HNF4A, VIL1, KLF4
        HITL-2: "Shall I proceed with these? Add/remove any?"

User:   "Yes, also add ISX."

Agent:  [TF Intelligence] CDX2 targets: MUC2, LYZ, FABP1 (n=312); GATA4: HAND1, NKX2-5...
        [Dataset Discovery] Found 8 relevant datasets on GEO + ArrayExpress:
          GSE158702 | Fetal intestinal development | 50k cells | h5ad available
          E-MTAB-9543 | Human intestine atlas      | 428k cells | h5ad available
        HITL-3: "Which dataset would you like to use?"

User:   "GSE158702."

Agent:  [Data Acquisition] Downloading GSE158702...
        [Preprocessing] 47,382 cells after QC. UMAP + cell type annotation done.
        HITL-4: [shows UMAP] "Do these annotations look correct?"

User:   "Cluster 5 should be Paneth cells."

Agent:  [Perturbation] Generating inference template (4 TFs, ~10k cells)...
        HITL-5: "Run: sbatch run_inference_GSE158702_CDX2_GATA4_HNF4A_ISX.sh (~12 hrs)"

User:   [12 hours later] "Job is done."

Agent:  ╔═══════╦══════════════════════════╦══════╦══════════════════╗
        ║ TF    ║ Top Activated Pathway    ║  NES ║ Key Targets      ║
        ╠═══════╬══════════════════════════╬══════╬══════════════════╣
        ║ CDX2  ║ Intestinal absorptive    ║ 2.41 ║ MUC2, LYZ, FABP1 ║
        ║ GATA4 ║ Epithelial diff.         ║ 1.98 ║ HNF1A, ONECUT1   ║
        ║ HNF4A ║ Lipid metabolism         ║ 1.87 ║ APOB, APOA4      ║
        ║ ISX   ║ Carotenoid absorption    ║ 1.43 ║ BCO1, SCARB1     ║
        ╚═══════╩══════════════════════════╩══════╩══════════════════╝
        Recommendation: CDX2 (primary) + HNF4A (metabolic maturation).
```

---

## Citation / Acknowledgements

- STATE model: [Genentech/state](https://github.com/genentech/state)
- Example fetal heart data: Jiajun Chen (Chen Lab, Cornell)
- GSEA: [gseapy](https://gseapy.readthedocs.io/)
