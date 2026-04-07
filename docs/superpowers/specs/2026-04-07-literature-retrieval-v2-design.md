# Literature Retrieval v2 — Design Spec

## Problem

The current `mine_literature.py` pipeline has poor precision (~3-5%) and moderate recall (~35-70%) when evaluated against two curated ground truth tables of scRNA-seq dataset papers (islet: 27 papers, gut: 26 papers).

**Root causes:**
- Query expansion generates too many sub-queries (disease terms, synonyms) that flood results with irrelevant papers
- No filter for whether a paper actually produced scRNA-seq data vs merely mentioning it
- Misses multi-organ atlases, embryo/stem cell papers, and organoid papers where scRNA-seq isn't in the abstract
- PubMed-only retrieval misses the data-repository angle entirely

## Design

A two-channel retrieval pipeline with LLM-based relevance filtering.

### Architecture

```
User query (e.g. "human gut scRNA-seq")
        │
        ├──► Channel 1: PubMed broad retrieval
        │       - Improved query expansion (add atlas, organoid, embryo, iPSC categories)
        │       - Remove noisy disease sub-queries (or cap them)
        │       - Fetch title + abstract + metadata via eutils
        │
        ├──► Channel 2: GEO + ArrayExpress direct search
        │       - Search NCBI GDS (gds database) for scRNA-seq datasets matching tissue
        │       - Search EBI ArrayExpress for matching experiments
        │       - Resolve dataset → PMID links
        │       - Fetch paper metadata via eutils
        │
        ▼
   Merge & deduplicate by PMID
        │
        ▼
   LLM relevance filter (Claude API)
        - Input: title + abstract for each paper
        - Prompt: "Does this paper report/re-analyze scRNA-seq data
          that includes [tissue]? Consider multi-organ atlases and
          developmental lineage."
        - Output: yes/no + confidence (high/medium/low) + reasoning
        - Batch abstracts to minimize API calls
        │
        ▼
   Output: filtered CSV
        - Columns: PMID, Title, Abstract, Year, Journal, Authors, DOI,
          source (PubMed/GEO/ArrayExpress), GEO_accession, ArrayExpress_accession,
          LLM_relevance (yes/no), LLM_confidence (high/medium/low), LLM_reasoning
```

### Channel 1: Improved PubMed Retrieval

Refactor `expand_query()` to generate queries along these axes:

1. **Core query**: user query as-is
2. **Technology variants**: scRNA-seq, snRNA-seq, single-cell transcriptome, single-nucleus RNA
3. **Tissue synonyms**: same as current, but capped at top 3 synonyms
4. **Organoid/iPSC**: `"[tissue] organoid single-cell"`, `"pluripotent stem cell [tissue] differentiation"`
5. **Atlas/multi-organ**: `"human cell atlas single-cell"`, `"Tabula Sapiens"`, `"fetal gene expression atlas"`
6. **Developmental**: `"[tissue] fetal development single-cell"`, `"embryo organogenesis single-cell"`

**Removed**: Disease-specific sub-queries (IBD, Crohn, T2D) — these added 100+ irrelevant papers. If the user wants disease context, they can include it in the base query.

Max results per sub-query: 50 (unchanged).

### Channel 2: GEO + ArrayExpress Search

**GEO (NCBI GDS database)**:
- Use eutils `esearch` on `db=gds` with query like `"human islet" AND ("scRNA-seq" OR "single cell RNA")`
- Parse GDS records to extract: GEO accession (GSE), linked PMID, title, summary
- Also search `db=gds` for `"10x Genomics" AND "[tissue]"` to catch platform-tagged datasets

**ArrayExpress (EBI)**:
- Use the ArrayExpress REST API: `https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?keywords=[query]`
- Filter by experiment type containing "RNA-seq of coding RNA from single cells"
- Extract: accession (E-MTAB-*), title, linked publications

Both channels resolve to PMIDs via NCBI ID converter or direct field extraction.

### LLM Filter

**Model**: Claude (via Anthropic API, model configurable, default `claude-sonnet-4-20250514`)

**Batching**: Send 10 abstracts per API call to reduce overhead. Each batch is a single prompt with structured output.

**Prompt template**:
```
You are a bioinformatics expert. For each paper below, determine whether it
reports or re-analyzes single-cell or single-nucleus RNA-seq data that includes
{tissue_description}.

Consider:
- Multi-organ atlases that include this tissue
- Developmental biology: embryo/fetal datasets where this tissue lineage originates
- Organoid or iPSC-derived models of this tissue
- Papers that generated original data vs reviews/methods papers that don't

For each paper, respond with:
- relevant: yes or no
- confidence: high, medium, or low
- reasoning: one sentence explanation

Papers:
{papers_json}
```

**Output parsing**: Expect JSON array response. Papers marked `yes` with `high` or `medium` confidence are included in final output.

**Cost estimate**: ~500 abstracts / 10 per batch = 50 API calls. At ~1000 tokens per call, roughly $0.02-0.05 total per query.

**API key**: Passed via `--api-key` CLI argument.

### CLI Interface

```bash
python mine_literature_v2.py \
    --query "human pancreatic islet scRNA-seq" \
    --api-key "sk-ant-..." \
    --output results.csv \
    --max 50 \
    --model claude-sonnet-4-20250514 \
    --no-geo           # skip GEO/ArrayExpress channel
    --no-filter        # skip LLM filter (return raw merged results)
    --use-hirn         # force HIRN search
    --verbose          # print detailed progress
```

### Output Format

CSV columns:
- `PMID` — PubMed ID
- `Title` — paper title
- `Abstract` — abstract text
- `Year` — publication year
- `Journal` — journal name
- `Authors` — first 3 authors + "et al."
- `DOI` — digital object identifier
- `source` — comma-separated: PubMed, GEO, ArrayExpress, HIRN
- `GEO_accession` — GSE ID if found via GEO channel
- `ArrayExpress_accession` — E-MTAB ID if found via ArrayExpress channel
- `LLM_relevant` — yes/no
- `LLM_confidence` — high/medium/low
- `LLM_reasoning` — one-sentence explanation

### File Structure

```
agents/literature_mining/
├── mine_literature_v2.py      ← new main orchestrator
├── pubmed_search.py           ← unchanged (reused)
├── hirn_search.py             ← unchanged (reused)
├── geo_search.py              ← new: GEO/ArrayExpress retrieval
├── llm_filter.py              ← new: Claude API relevance filter
└── query_expansion.py         ← new: refactored query expansion (extracted from mine_literature.py)
```

### Evaluation

After implementation, re-run against both ground truth tables:
- Target precision: >40% (up from 3-5%)
- Target recall: >80% (up from 35-70%)
- Target F1: >50% (up from 6-9%)

### Dependencies

- `anthropic` Python SDK (for LLM filter)
- Existing: `urllib`, `xml.etree`, `csv`, `json` (all stdlib)

### Not in scope

- Structured metadata extraction from abstracts (separate tool)
- Full-text retrieval/parsing (future step after retrieval is validated)
- XLSX output format (CSV only; user converts as needed)
