# Calibration / Ground-Truth Notes

## Source

`docs/originals/AGImodel.xlsx` (collaborator-curated, provided 2026-04-13).

The spreadsheet defines, for 17 cell types across 2 tissues (gut and human islets),
the **expected** answers for the three deliverable streams from the
[analysis pipeline plan](../../docs/analysis_pipeline_plan.md):

1. **Stream 1**: Transcription factors that should be tested in differentiation
2. **Stream 2**: Pathway → chemical/growth factor mapping for differentiation media
3. **Stream 3**: ChemPerturb-seq hits (5 per cell type)

These act as ground truth for evaluating Phase 7 (perturbation scoring) and
Phase 8 (chemical aggregation), and for benchmarking different perturbation
models (state_sm, state_base, future variants).

## Extracted YAMLs

| File | Stream | Use |
|------|--------|-----|
| `expected_tfs_by_celltype.yaml` | 1 | Phase 7 positive-control TF list per cell type. Compute top-K recall. |
| `differentiation_compounds.yaml` | 2 | Phase 8 curated compound table. Replaces the manual ~40-compound table the plan said we would build. |

Stream 3 is intentionally **not** extracted into a YAML — see "Caveats" below.

## Cell Type Coverage

**Gut (11 cell types)**, training on Adult or Fetal+Adult gut:
D, Enterochromaffin, Enterocytes, Enteroendocrine, Goblet, I, K, L, Stem, TA, X

**Human islets (6 cell types)**, training on HPAP:
Acinar, Alpha, Beta, Delta, Ductal, Enterochromaffin

## Caveats

### Spreadsheet structure

The xlsx uses **merged-header columns**. The column labeled
`2. Chemicals/growth factors from pathway analysis...` actually contains the
**pathway names**; the unlabeled column to its right (`Unnamed: 5` in pandas)
contains the **chemicals**. Verify against row 2 (sub-header) before parsing —
the first parser pass got it backwards and produced nonsense entries like
`pathway: "DAPT, Dibenzazepine"`.

### Stream 3 (ChemPerturb-seq hits) is suspect

11 of 17 cell types have **identical** ChemPerturb-seq hits:
`ozanimod, ketotifen_fumarate, duloxetine_hydrochloride, deflazacort, dimethyl_fumarate`.

Only Alpha, Beta, Delta, and (separately) human-islet Enterochromaffin have
unique values. This pattern suggests the column is either:

- A placeholder waiting for the collaborator to fill in
- Default top-5 hits from a baseline analysis the collaborator already ran
  (possibly using prior code we haven't seen)
- A join error in the spreadsheet population

**Action:** Do not calibrate Phase 8 / Stream 3 against this column until
the collaborator clarifies. We extracted it into the parsed JSON for
reference but did not promote it to a YAML.

## How To Use For Evaluation

```python
import yaml
from pathlib import Path

config_dir = Path("insilico_perturbation_pipeline/configs")
expected = yaml.safe_load((config_dir / "expected_tfs_by_celltype.yaml").read_text())

beta = expected["cell_types"]["human_islets/beta_cells"]
print(beta["expected_tfs"])
# ['PDX1', 'NKX6-1', 'MAFA', 'NEUROD1', 'ISL1', 'RFX6', 'PAX6', 'NKX2-2', 'UCN3']
```

Then compute the top-K recall metric in Phase 7 evaluation:
```
top_k_recall = |expected_tfs ∩ top_K_predicted| / |expected_tfs|
```

Beta cells are recommended as the **first cell type to evaluate** because the
expected TFs are completely unambiguous (every beta-cell biology textbook
contains PDX1, NKX6-1, MAFA), and we already have HPAP predictions in
`runs/collab_filtered_v1/human_islets/beta/`.
