# acquire-datasets skill

Use this skill when the user wants to download scRNA-seq datasets from a catalog file (CSV or xlsx).

## What this skill does

Reads a catalog CSV/xlsx listing manuscript names and database accessions, then downloads each publicly accessible dataset as an h5ad file.

Supported databases:
- **GEO** (GSE*) — NCBI, fully public
- **EBI ArrayExpress** (E-MTAB-*) — EBI BioStudies, fully public
- **HuBMAP** (HBM*) — public datasets (protected ones need a token)
- **EGA / CNCB / dbGaP / PANC DB** — controlled access, reported but not downloaded

## How to invoke

When the user says something like:
- "download the datasets from this spreadsheet"
- "fetch the data from my catalog"
- "acquire the scRNA-seq datasets listed in this CSV"
- `/acquire-datasets`

Follow these steps:

### Step 1 — Identify the catalog file

Ask if not provided. Accepts `.csv` (preferred) or `.xlsx`.

Required columns (detected automatically by name):
- A title/manuscript column
- An accession/data resource column (GSE*, E-MTAB-*, HBM*, etc.)
- A URL/link column (used as fallback)

### Step 2 — Run a dry run first

```bash
conda activate /nfs/turbo/umms-drjieliu/usr/zheyuz/miniforge-pypy3/envs/state_env
cd /nfs/turbo/umms-drjieliu/usr/rickyhan/AGI-guided_differentiation

python agents/data_acquisition/acquire_datasets.py \
    --catalog <PATH_TO_CATALOG> \
    --output-dir data/downloaded \
    --dry-run
```

Show the user the summary table of what will be downloaded vs skipped.

### Step 3 — Confirm with user

Tell the user:
- How many datasets are publicly accessible and will be downloaded
- Which accessions are controlled-access (EGA, dbGaP, etc.) and will be skipped
- Estimated disk space (rough: 100–500 MB per dataset)
- Ask: "Shall I proceed with downloading all public datasets?"

### Step 4 — Run the download

```bash
python agents/data_acquisition/acquire_datasets.py \
    --catalog <PATH_TO_CATALOG> \
    --output-dir data/downloaded
```

To download a single specific accession:
```bash
python agents/data_acquisition/acquire_datasets.py \
    --catalog <PATH_TO_CATALOG> \
    --output-dir data/downloaded \
    --accession GSE158702
```

### Step 5 — Report results

After completion, show the user:
- Which datasets were successfully downloaded and their file paths
- Which failed and why (e.g. FASTQ-only on GEO)
- Which were skipped (controlled access) with links to apply for access

## Output

All h5ad files saved to `data/downloaded/{accession}.h5ad`

## Troubleshooting

- **GEO dataset has only FASTQs**: The script will report this. The user must download from SRA manually.
- **EBI dataset not found**: Check if the accession is correct; some older E-MTAB studies use different URL formats.
- **HuBMAP requires token**: Run with `--hubmap-token <TOKEN>`. Token obtained from https://portal.hubmapconsortium.org (sign in → profile → token).
- **Large dataset (>10 GB)**: Downloads stream in chunks. Ensure sufficient disk space before starting.
