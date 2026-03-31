# HIRN Literature Retrieve

A Claude Code plugin for searching and retrieving relevant passages from articles published by the [Human Islet Research Network (HIRN)](https://hirnetwork.org). Covers ~1,160 publications across 9 HIRN consortia.

## Features

- **Search HIRN publications** — keyword search across all HIRN published works with optional consortium filtering
- **Resolve article identifiers** — batch PMID-to-PMCID resolution via the NCBI ID Converter API
- **Fetch full text** — retrieve structured article text from PubMed Central Open Access (BioC JSON)
- **Chunk and rank passages** — split articles into retrievable segments and rank by BM25 relevance
- **Automatic caching** — file-based cache with TTL to minimize redundant API calls
- **Rate limiting** — token-bucket rate limiter respecting NCBI API limits

## Installation

Requires Claude Code v1.0.33 or later.

### From GitHub

```bash
/plugin marketplace add yuanhao96/hirn_publication_retrieval
/plugin install hirn-literature-retrieve@hirn-literature-retrieve
```

### From a Local Clone

```bash
git clone git@github.com:yuanhao96/hirn_publication_retrieval.git
/plugin marketplace add ./hirn_publication_retrieval
/plugin install hirn-literature-retrieve@hirn-literature-retrieve
```

After installation, the skill is automatically available when Claude Code detects HIRN-related queries. No external Python dependencies are required — the skill uses only the Python standard library.

## Configuration

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `NCBI_API_KEY` | No | None | [NCBI API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) for 10 requests/sec (vs 3 without) |

## Usage

The skill activates automatically when Claude Code encounters HIRN-related queries. You can also invoke it explicitly:

### Search HIRN publications

```bash
python -c "
import json
from scripts.scrape_hirn import fetch_hirn_publications, search_publications

pubs = fetch_hirn_publications()
results = search_publications(pubs, query='islet autoimmunity', max_results=10)
print(json.dumps(results, indent=2))
"
```

### Filter by consortium

```bash
python -c "
from scripts.scrape_hirn import fetch_hirn_publications, search_publications

pubs = fetch_hirn_publications()
results = search_publications(pubs, query='beta cell', consortium='HPAC')
for r in results:
    print(f'PMID {r[\"pmid\"]}: {r[\"title\"]}')
"
```

Available consortia: `CBDS`, `CHIB`, `CMAD`, `CMAI`, `CTAR`, `HIREC`, `HPAC`, `Opportunity Pool`, `PanKbase`

### End-to-end: query to relevant passages

```bash
python -c "
import json
from scripts.scrape_hirn import fetch_hirn_publications, search_publications
from scripts.resolve_ids import resolve_pmcids
from scripts.fetch_fulltext import fetch_fulltext
from scripts.chunk_text import chunk_passages
from scripts.search_chunks import search_chunks

query = 'islet autoimmunity HLA markers'

# Fetch index, search titles, resolve PMCIDs
pubs = fetch_hirn_publications()
matches = search_publications(pubs, query=query, max_results=5)
pmcid_map = resolve_pmcids([p['pmid'] for p in matches])

# Fetch full text and search chunks
all_hits = []
for pub in matches:
    pmcid = pmcid_map.get(pub['pmid'])
    if not pmcid:
        continue
    ft = fetch_fulltext(pmcid)
    if not ft['success']:
        continue
    chunks = chunk_passages(ft['passages'])
    hits = search_chunks(chunks, query=query, top_k=3)
    for h in hits:
        h['pmid'] = pub['pmid']
        h['pmcid'] = pmcid
        h['article_title'] = pub['title']
    all_hits.extend(hits)

all_hits.sort(key=lambda x: x['score'], reverse=True)
print(json.dumps(all_hits[:10], indent=2))
"
```

### Output format

Each returned chunk contains:

```json
{
  "text": "The passage text...",
  "section": "RESULTS",
  "score": 3.14,
  "pmid": "39630627",
  "pmcid": "PMC11615173",
  "article_title": "Islet autoimmunity and HLA markers...",
  "chunk_index": 4
}
```

## Pipeline

```
User Query
    |
    v
 Fetch HIRN Index -----> AJAX scrape of hirnetwork.org (~1,160 pubs)
    |                     Cached 24 hours
    v
 Search Titles ---------> Keyword matching, optional consortium filter
    |
    v
 Resolve PMCIDs -------> NCBI ID Converter API (batch, cached 30 days)
    |
    v
 Fetch Full Text ------> PMC BioC JSON API (cached 30 days)
    |
    v
 Chunk & Rank ----------> Section/paragraph chunking + BM25 scoring
    |
    v
 Relevant passages with source citations
```

## Project Structure

```
.claude-plugin/
  plugin.json                       # Plugin manifest
skills/hirn-literature-retrieve/
  SKILL.md                          # Skill definition (auto-discovered)
  scripts/
    scrape_hirn.py                  # Fetch HIRN publications via AJAX
    resolve_ids.py                  # PMID -> PMCID resolution
    fetch_fulltext.py               # PMC full text via BioC JSON
    chunk_text.py                   # Section/paragraph chunking
    search_chunks.py                # BM25 ranking
    utils/
      rate_limiter.py               # Token-bucket rate limiter
      cache_manager.py              # File-based cache with TTL
      html_parser.py                # Parse HIRN publication HTML
  data/cache/                       # Runtime cache (gitignored)
tests/
  test_*.py                         # 43 tests covering all modules
  conftest.py                       # Shared fixtures and sample data
conftest.py                         # Root conftest (adds skill dir to sys.path)
```

## Testing

```bash
python -m pytest tests/ -v
```

All 43 tests run with mocked HTTP — no network access required.

## Limitations

- **PMC Open Access only** — full text is available for articles in the PMC Open Access subset (~50-60% of HIRN publications). Others return title/metadata only.
- **Title-based initial search** — the HIRN listing does not include abstracts, so the first search step matches on titles. Full-text chunk search provides deeper retrieval.
- **Rate limits** — NCBI APIs allow 3 requests/second without an API key (10 with one). Large batch operations are handled automatically but may take time.

## Data Sources

- [HIRN Published Works](https://hirnetwork.org/published-works) — WordPress AJAX endpoint
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) — article metadata
- [NCBI ID Converter](https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/) — PMID/DOI to PMCID mapping
- [PMC BioC API](https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi) — structured full text

## License

MIT
