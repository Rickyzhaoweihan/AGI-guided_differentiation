# Perturbation Report: HES1 knockdown in Stem cells

**Tissue:** adult gut  
**QC Status:** WARN

> **Warning (self_knockdown):** HES1 log2FC = 0.0121 (expected < 0). May indicate model didn't encode self-suppression, or gene has autoregulatory feedback.

## Sanity Checks

| Check | Result | Details |
|-------|--------|---------|
| housekeeping_stability | PASS | All 5 HK genes stable (max |log2FC| = 0.4449) |
| effect_magnitude | PASS | Median |log2FC| = 0.0014 |
| self_knockdown | WARN | HES1 log2FC = 0.0121 (expected < 0). May indicate model didn't encode self-suppression, or gene has  |
| known_targets | PASS | 29/91 (32%) targets responded as expected. BCL2: +0.000 (unexpected); COL1A2: +0.001 (unexpected); S |

## Differential Expression

- **558** genes significantly changed (|log2FC| > 0.5, FDR < 0.05)

**Top upregulated:**

| Gene | log2FC | 95% CI | FDR |
|------|--------|--------|-----|
| MT-ATP6 | +6.222 | [-57.880, 7.483] | 6.07e-01 |
| MT-CO3 | +5.517 | [-56.811, 6.793] | 6.24e-01 |
| MT-ND4 | +5.148 | [-54.220, 6.409] | 6.17e-01 |
| MT-CO2 | +5.061 | [-55.145, 6.315] | 6.13e-01 |
| MT-ND3 | +4.976 | [-52.368, 6.221] | 6.07e-01 |
| CCDC115 * | +3.074 | [2.056, 3.665] | 9.16e-03 |
| ERCC6L2 * | +2.638 | [1.640, 3.222] | 1.28e-02 |
| WBP2 | +2.179 | [0.509, 2.932] | 1.65e-01 |
| SOX18 | +2.097 | [0.866, 2.752] | 7.15e-02 |
| SOCS4 * | +2.085 | [1.382, 2.556] | 2.11e-03 |

**Top downregulated:**

| Gene | log2FC | 95% CI | FDR |
|------|--------|--------|-----|
| TIMM29 | -2.292 | [-35.534, 0.282] | 5.06e-01 |
| ZNF800 | -1.768 | [-35.633, 0.203] | 4.57e-01 |
| HMGB2 * | -1.586 | [-6.038, -0.619] | 9.61e-04 |
| SET * | -1.365 | [-3.107, -0.599] | 3.25e-04 |
| S100A6 * | -1.343 | [-1.957, -0.913] | 2.26e-14 |
| MDH1 * | -1.282 | [-1.885, -0.858] | 6.00e-13 |
| SNRPD1 * | -1.280 | [-1.871, -0.862] | 2.78e-13 |
| ALB * | -1.265 | [-1.788, -0.882] | 1.95e-15 |
| PHLDA1 * | -1.215 | [-2.021, -0.700] | 2.41e-07 |
| APOC1 * | -1.214 | [-1.635, -0.888] | 3.42e-19 |

## Pathway Enrichment

### GSEA (pre-ranked, FDR < 0.25)

**MSigDB_Hallmark_2020** (5 significant):

| Pathway | NES | FDR | Direction |
|---------|-----|-----|-----------|
| Pancreas Beta Cells | -1.73 | 0.006 | DOWN |
| G2-M Checkpoint | -1.50 | 0.085 | DOWN |
| Myc Targets V1 | -1.46 | 0.087 | DOWN |
| E2F Targets | -1.29 | 0.223 | DOWN |
| Coagulation | -1.29 | 0.186 | DOWN |

**GO_Biological_Process_2023** (46 significant):

| Pathway | NES | FDR | Direction |
|---------|-----|-----|-----------|
| Potassium Ion Transport (GO:0006813) | -1.97 | 0.000 | DOWN |
| B Cell Proliferation (GO:0042100) | -1.89 | 0.001 | DOWN |
| Monocyte Chemotaxis (GO:0002548) | -1.82 | 0.008 | DOWN |
| Calcium-Dependent Cell-Cell Adhesion Via Plasma Membrane Cell Adhesion Molecules (GO:0016339) | -1.74 | 0.053 | DOWN |
| Defense Response To Gram-positive Bacterium (GO:0050830) | -1.72 | 0.073 | DOWN |
| Neuropeptide Signaling Pathway (GO:0007218) | 1.48 | 0.091 | UP |
| Regulation Of AMPA Receptor Activity (GO:2000311) | 1.45 | 0.221 | UP |
| Negative Regulation Of Notch Signaling Pathway (GO:0045746) | 1.43 | 0.243 | UP |

**KEGG_2021_Human** (7 significant):

| Pathway | NES | FDR | Direction |
|---------|-----|-----|-----------|
| Systemic lupus erythematosus | -1.85 | 0.000 | DOWN |
| Type I diabetes mellitus | -1.72 | 0.013 | DOWN |
| Malaria | -1.71 | 0.010 | DOWN |
| Hematopoietic cell lineage | -1.62 | 0.051 | DOWN |
| Primary immunodeficiency | -1.57 | 0.083 | DOWN |

### Over-Representation Analysis


## Regulatory Analysis

### TF Activity (decoupleR ULM + CollecTRI)

| TF | Activity Score | p-value |
|----|---------------|---------|
| ZBED1 | +5.504 (activated) | 2.43e-05 |
| HSF4 | -5.109 (repressed) | 1.05e-04 |
| TBPL1 | +4.832 (activated) | 2.92e-04 |
| ARID4B | -4.364 (repressed) | 1.87e-03 |
| VEZF1 | +4.335 (activated) | 1.87e-03 |
| RELA | -3.651 (repressed) | 2.64e-02 |
| NFKB | -3.606 (repressed) | 2.64e-02 |
| AHRR | -3.591 (repressed) | 2.64e-02 |
| ONECUT1 | -3.476 (repressed) | 3.50e-02 |
| NR1H3 | -3.457 (repressed) | 3.50e-02 |
| STAT3 | -3.289 (repressed) | 5.85e-02 |
| TTF1 | +3.085 (activated) | 9.89e-02 |
| HNF1A | -3.071 (repressed) | 9.89e-02 |
| PIN1 | -3.046 (repressed) | 9.89e-02 |
| ZNF148 | -3.046 (repressed) | 9.89e-02 |

### Pathway Activity (PROGENy)

| Pathway | Activity Score | p-value |
|---------|---------------|---------|
| p53 | -11.987 | 0.00e+00 |
| Hypoxia | -3.179 | 1.48e-03 |
| EGFR | -2.714 | 6.66e-03 |
| MAPK | -2.473 | 1.34e-02 |
| Estrogen | -1.531 | 1.26e-01 |
| WNT | -1.341 | 1.80e-01 |
| TNFa | -1.267 | 2.05e-01 |
| PI3K | -1.188 | 2.35e-01 |
| Trail | -1.154 | 2.49e-01 |
| Androgen | +0.857 | 3.91e-01 |
| NFkB | +0.633 | 5.27e-01 |
| JAK-STAT | -0.529 | 5.97e-01 |
| TGFb | -0.269 | 7.88e-01 |
| VEGF | +0.173 | 8.63e-01 |

## Cell State Impact

*Note: pro-apoptotic and anti-apoptotic scored separately*

| State | Control | Perturbed | Delta | Interpretation |
|-------|---------|-----------|-------|----------------|
| proliferation | -0.003 | -0.072 | -0.069 | Proliferation slightly decreased |
| pro_apoptotic | -0.073 | -0.111 | -0.038 | no change |
| anti_apoptotic | -0.009 | -0.021 | -0.011 | no change |
| stress_upr | 0.543 | 0.457 | -0.086 | Stress/UPR response slightly decreased |
| gut_stem_markers | 0.001 | -0.002 | -0.002 | no change |
| enterocyte_markers | 0.001 | -0.008 | -0.009 | no change |
| goblet_markers | -0.000 | 0.000 | +0.000 | no change |
| paneth_markers | -0.010 | -0.012 | -0.002 | no change |
| enteroendocrine_markers | 0.000 | -0.000 | -0.001 | no change |
| wnt_targets | 0.030 | -0.116 | -0.146 | wnt_targets moderately decreased |
| notch_targets | -0.001 | 0.001 | +0.002 | no change |

## Literature Context

**HES1**: hes family bHLH transcription factor 1

> This protein belongs to the basic helix-loop-helix family of transcription factors. It is a transcriptional repressor of genes that require a bHLH protein for their transcription. The protein has a particular type of basic domain that contains a helix interrupting protein that binds to the N-box rat

**gene_in_context** (5 results):

- [Modelling human hepato-biliary-pancreatic organogenesis from the foregut-midgut boundary.](https://pubmed.ncbi.nlm.nih.gov/31554966/) (Koike, Iwasawa, Ouchi et al., 2019)
  > Organogenesis is a complex and interconnected process that is orchestrated by multiple boundary tissue interactions...
- [Cell cycle arrest determines adult neural stem cell ontogeny by an embryonic Notch-nonoscillatory Hey1 module.](https://pubmed.ncbi.nlm.nih.gov/34772946/) (Harada, Yamada, Imayoshi et al., 2021)
  > Quiescent neural stem cells (NSCs) in the adult mouse brain are the source of neurogenesis that regulates innate and adaptive behaviors. Adult NSCs in the subventricular zone are derived from a subpop...
- [HES1 revitalizes the functionality of aged adipose-derived stem cells by inhibiting the transcription of STAT1.](https://pubmed.ncbi.nlm.nih.gov/39501364/) (Li, Ren, Yan et al., 2024)
  > The effectiveness of adipose-derived stem cells (ADSCs) in therapy diminishes with age. It has been reported that transcription factors (TFs) play a crucial role in the aging and functionality of stem...

**gene_perturbation** (5 results):

- [Endothelial POFUT1 controls injury-induced liver fibrosis by repressing fibrinogen synthesis.](https://pubmed.ncbi.nlm.nih.gov/38460791/) (He, Luo, Ma et al., 2024)
  > NOTCH signaling in liver sinusoidal endothelial cells (LSECs) regulates liver fibrosis, a pathological feature of chronic liver diseases. POFUT1 is an essential regulator of NOTCH signaling. Here, we ...
- [Impaired Glycosylation of Gastric Mucins Drives Gastric Tumorigenesis and Serves as a Novel Therapeutic Target.](https://pubmed.ncbi.nlm.nih.gov/38583723/) (Arai, Hayakawa, Tateno et al., 2024)
  > Gastric cancer is often accompanied by a loss of mucin 6 (MUC6), but its pathogenic role in gastric carcinogenesis remains unclear....
- [Myt1l safeguards neuronal identity by actively repressing many non-neuronal fates.](https://pubmed.ncbi.nlm.nih.gov/28379941/) (Mall, Kareta, Chanda et al., 2017)
  > Normal differentiation and induced reprogramming require the activation of target cell programs and silencing of donor cell programs. In reprogramming, the same factors are often used to reprogram man...

**gene_function** (5 results):

- [A Single-Cell Transcriptomic Atlas of Human Skin Aging.](https://pubmed.ncbi.nlm.nih.gov/33238152/) (Zou, Long, Zhao et al., 2021)
  > Skin undergoes constant self-renewal, and its functional decline is a visible consequence of aging. Understanding human skin aging requires in-depth knowledge of the molecular and functional propertie...
- [The cloning and activity of human Hes1 gene promoter.](https://pubmed.ncbi.nlm.nih.gov/29257279/) (Lu, Jiang, Gao, 2018)
  > The aim of the current study was to obtain and analyze the activity of the human Hes1 gene promoter. The genomic DNA of human HeLa cell was used as template, polymerase chain reaction (PCR) was used t...
- [The pivotal role of the Hes1/Piezo1 pathway in the pathophysiology of glucocorticoid-induced osteoporosis.](https://pubmed.ncbi.nlm.nih.gov/39641269/) (Ochiai, Etani, Noguchi et al., 2024)
  > Glucocorticoid-induced osteoporosis (GIOP) lacks fully effective treatments. This study investigated the role of Piezo1, a mechanosensitive ion channel component 1, in GIOP. We found reduced Piezo1 ex...

**gene_pathway_0** (2 results):

- [Endosulfan inhibits proliferation through the Notch signaling pathway in human umbilical vein endothelial cells.](https://pubmed.ncbi.nlm.nih.gov/27939630/) (Wei, Zhang, Ren et al., 2017)
  > Our previous research showed that endosulfan triggers the extrinsic coagulation pathway by damaging endothelial cells and causes hypercoagulation of blood. To identify the mechanism of endosulfan-impa...
- [Comprehensive transcriptomics and proteomics analysis of neointima formation in human saphenous vein: implications for bypass graft disease.](https://pubmed.ncbi.nlm.nih.gov/41569144/) (Kim, Goo, Shi et al., 2026)
  > Human saphenous veins (SVs) are widely used as grafts in coronary artery bypass (CABG) surgery but often fail due to neointima formation. Little is known, however, regarding the cellular, transcriptom...

## Evidence Convergence

| Dimension | Signal | Details |
|-----------|--------|---------|
| Mechanistic | + | 29/91 (32%) targets responded as expected. BCL2: +0.000 (unexpected); COL1A2: +0 |
| Functional | + | 58 significant pathways across libraries |
| Regulatory | + | 10 significant TFs (p<0.05). Top: ZBED1 (+5.50) |
| Pathway | + | 4 significant pathways (p<0.05) |
| Phenotypic | + | 3 cell states changed |
| Literature | + | 17 relevant publications found |
