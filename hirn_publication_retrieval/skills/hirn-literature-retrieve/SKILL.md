---
name: hirn-literature-retrieve
description: >-
  This skill should be used when the user asks to "search HIRN publications",
  "find HIRN research", "retrieve HIRN articles", "search Human Islet Research Network",
  "find papers from HIRN consortia", or mentions HIRN-specific terms like CBDS, CHIB,
  CMAD, CMAI, CTAR, HIREC, HPAC, PanKbase, or Opportunity Pool in the context of
  literature search. Also triggers when the user asks to "search islet research network
  publications" or "find HIRN published works".
---

# HIRN Literature Retrieve

Search and retrieve relevant passages from ~1,160 articles published by the Human Islet Research Network (HIRN) across 9 consortia.

## Pipeline Overview

```
User Query
    |
    v
 Fetch HIRN Index -----> AJAX scrape of hirnetwork.org (~1,160 pubs, cached 24h)
    |
    v
 Search Titles ---------> Keyword matching, optional consortium filter
    |
    v
 Resolve PMCIDs -------> NCBI ID Converter API (batch, cached 30 days)
    |
    v
 Fetch Full Text ------> PMC BioC JSON API (Open Access only, cached 30 days)
    |
    v
 Chunk & Rank ----------> Section/paragraph chunking + BM25 scoring
    |
    v
 Relevant passages with source citations
```

## Rate Limiting and Caching

All NCBI API calls are automatically rate-limited:
- Without `NCBI_API_KEY`: 3 requests/second
- With `NCBI_API_KEY` env var: 10 requests/second
- HIRN WordPress AJAX: 3 requests/second

File-based cache in `data/cache/` with TTLs: index (24h), PMCIDs (30d), full text (30d).

## Workflow

### Step 1: Fetch the HIRN publication index

```bash
cd {SKILL_DIR}
python -c "
import json
from scripts.scrape_hirn import fetch_hirn_publications
pubs = fetch_hirn_publications()
print(json.dumps({'count': len(pubs), 'sample': pubs[:3]}, indent=2))
"
```

Fetches all ~1,160 HIRN publications. Cached for 24 hours after first call.

### Step 2: Search publications by keyword

```bash
cd {SKILL_DIR}
python -c "
import json
from scripts.scrape_hirn import fetch_hirn_publications, search_publications
pubs = fetch_hirn_publications()
results = search_publications(pubs, query='USER_QUERY_HERE', max_results=10)
print(json.dumps(results, indent=2))
"
```

To filter by consortium, add `consortium='CBDS'` (or any consortium name: CBDS, CHIB, CMAD, CMAI, CTAR, HIREC, HPAC, Opportunity Pool, PanKbase).

### Step 3: Resolve PMCIDs for matched articles

```bash
cd {SKILL_DIR}
python -c "
import json
from scripts.resolve_ids import resolve_pmcids
pmids = ['PMID1', 'PMID2']  # from step 2 results
mapping = resolve_pmcids(pmids)
print(json.dumps(mapping, indent=2))
"
```

### Step 4: Fetch full text from PMC

```bash
cd {SKILL_DIR}
python -c "
import json
from scripts.fetch_fulltext import fetch_fulltext
result = fetch_fulltext('PMCID_HERE')
if result['success']:
    print(f'Passages: {len(result[\"passages\"])}')
    for p in result['passages'][:5]:
        print(f'  [{p[\"section\"]}] {p[\"text\"][:100]}...')
else:
    print(f'Error: {result[\"error\"][\"message\"]}')
"
```

### Step 5: Chunk and search for relevant passages

```bash
cd {SKILL_DIR}
python -c "
import json
from scripts.fetch_fulltext import fetch_fulltext
from scripts.chunk_text import chunk_passages
from scripts.search_chunks import search_chunks

result = fetch_fulltext('PMCID_HERE')
if result['success']:
    chunks = chunk_passages(result['passages'])
    hits = search_chunks(chunks, query='USER_QUERY_HERE', top_k=5)
    for h in hits:
        print(f'[{h[\"section\"]}] (score={h[\"score\"]})')
        print(f'  {h[\"text\"][:200]}...')
        print()
"
```

## Complete End-to-End Workflow

For a full search-to-chunks pipeline, run all steps together:

```bash
cd {SKILL_DIR}
python -c "
import json
from scripts.scrape_hirn import fetch_hirn_publications, search_publications
from scripts.resolve_ids import resolve_pmcids
from scripts.fetch_fulltext import fetch_fulltext
from scripts.chunk_text import chunk_passages
from scripts.search_chunks import search_chunks

query = 'USER_QUERY_HERE'
consortium = None  # or 'CBDS', 'HPAC', etc.

# 1. Fetch and search HIRN index
pubs = fetch_hirn_publications()
matches = search_publications(pubs, query=query, consortium=consortium, max_results=5)
print(f'Found {len(matches)} matching publications')

# 2. Resolve PMCIDs
pmids = [p['pmid'] for p in matches if p.get('pmid')]
pmcid_map = resolve_pmcids(pmids)

# 3. Fetch full text and search chunks for each article
all_results = []
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
        h['doi'] = pub.get('doi', '')
        h['authors'] = pub.get('authors', '')
    all_results.extend(hits)

# Sort all chunks by score
all_results.sort(key=lambda x: x['score'], reverse=True)

# Output top results
print(json.dumps(all_results[:10], indent=2))
"
```

## Output Format

Each returned chunk contains:

```json
{
    "text": "The passage text...",
    "section": "RESULTS",
    "type": "paragraph",
    "offset": 800,
    "chunk_index": 4,
    "score": 3.1415,
    "pmid": "39630627",
    "pmcid": "PMC11615173",
    "article_title": "Islet autoimmunity and HLA markers...",
    "doi": "10.1007/s00125-024-06244-0",
    "authors": "Jacobsen LM, Bocchino LE, et al."
}
```

## Reference Formatting

When citing retrieved passages, use NLM/Vancouver format:

```
AuthorLastName Initials, et al. Article title. Journal. Year;Vol(Issue):Pages. doi:DOI. PMID: XXXXX.
```

## Limitations

- **PMC Open Access only**: Full text available only for PMC Open Access articles (~50-60% of HIRN publications). Others return title/metadata only.
- **Title-based initial search**: The HIRN listing does not include abstracts; initial matching is on titles. Full-text chunk search provides deeper retrieval.
- **Rate limits**: NCBI APIs are rate-limited. Handled automatically but large batches may take time.

## Environment Variables

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `NCBI_API_KEY` | No | None | NCBI API key for 10 req/sec (vs 3) |
