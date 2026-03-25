# Skill: gsea-analysis

Run GSEA on in-silico perturbation predictions. See full documentation and usage in `GSEA_skill/skill.md`.

## Quick start

```bash
conda activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
pip install gseapy  # first time only

python GSEA_skill/run_gsea.py \
    --predictions insilico_perturbation_pipeline/data/perturbation_predictions.h5ad \
    --output-dir  GSEA_skill/results \
    --gene-sets   MSigDB_Hallmark_2020
```

## Input
- `perturbation_predictions.h5ad` — output of `run-perturbation` skill

## Output
- `GSEA_skill/results/<GENE>/log2fc.csv` — ranked gene list
- `GSEA_skill/results/<GENE>/volcano.png` — predicted DE plot
- `GSEA_skill/results/<GENE>/gsea_*/` — gseapy GSEA output
- `GSEA_skill/results/summary_all_perturbations.csv` — top pathways across all genes

## Key options
- `--target-gene TBX5` — run for one gene only
- `--per-cell-type` — separate GSEA per cell type
- `--gene-sets KEGG_2021_Human` — use KEGG instead of Hallmark

For full details: `GSEA_skill/skill.md`
