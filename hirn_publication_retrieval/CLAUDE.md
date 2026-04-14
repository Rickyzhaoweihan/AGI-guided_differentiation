# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

A Claude Code skill for retrieving published articles from the **Human Islet Research Network (HIRN)** and locating relevant text chunks from those articles via **PubMed Central (PMC)**.

Workflow: HIRN published works → extract PMIDs/DOIs → fetch full text from PMC → chunk and search for relevant passages.

## Data Sources

### HIRN Published Works
- URL: https://hirnetwork.org/published-works
- ~1,160 publications across consortia (CBDS, CHIB, CMAD, CMAI, CTAR, HIREC, HPAC, Opportunity Pool, PanKbase)
- WordPress site using Ajax Load More plugin (v7.8.2)
- AJAX endpoint: `https://hirnetwork.org/2021/wp-admin/admin-ajax.php`
- Post type: `publication`, taxonomy: `publication_category`, term: `published_works`
- Posts load 250 at a time; publication metadata (PMIDs, DOIs, titles, authors) is rendered dynamically

### PubMed Central (PMC)
- NCBI E-utilities for article metadata: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`
- PMC Open Access full-text API: `https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{PMCID}/unicode`
- BioC format returns structured sections (title, abstract, body paragraphs) with offset metadata
- Rate limit: 3 requests/second without API key, 10/second with NCBI API key (`NCBI_API_KEY` env var)

## Architecture Notes

This project is a **Claude Code plugin** containing a skill for HIRN literature retrieval. It follows the standard plugin layout:

```
.claude-plugin/plugin.json          # Plugin manifest
skills/hirn-literature-retrieve/    # Skill directory
  SKILL.md                          # Skill definition (auto-discovered)
  scripts/                          # Python scripts for the pipeline
  data/cache/                       # Runtime file cache (gitignored)
tests/                              # Test suite (43 tests, mocked HTTP)
conftest.py                         # Root conftest adding skill dir to sys.path
```

The core pipeline:
1. **Scrape/fetch** HIRN publication list (handle AJAX pagination and consortium filtering)
2. **Extract identifiers** (PMID, DOI, PMCID) from each publication entry
3. **Resolve PMCIDs** — use NCBI ID Converter API (`https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/`) to map PMIDs/DOIs → PMCIDs
4. **Fetch full text** from PMC Open Access subset via BioC JSON API
5. **Chunk text** into retrievable passages (by section, paragraph, or sliding window)
6. **Search/rank chunks** against a user query for relevant passage retrieval
