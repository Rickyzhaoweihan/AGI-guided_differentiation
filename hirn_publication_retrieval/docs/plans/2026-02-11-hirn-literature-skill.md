# HIRN Literature Retrieve Skill — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a Claude Code skill that searches HIRN-published articles, fetches their full text from PMC, and returns relevant text chunks for a user query.

**Architecture:** A set of Python scripts (stdlib only, no external deps) invoked by Claude Code via Bash. The pipeline: scrape HIRN publication index via AJAX → keyword-search titles → resolve PMCIDs via NCBI ID Converter → fetch full text via BioC JSON API → chunk by section/paragraph → rank chunks by BM25 score → return top results. All API calls use token-bucket rate limiting and file-based TTL caching.

**Tech Stack:** Python 3.10+ (stdlib only: `urllib`, `json`, `html.parser`, `re`, `math`, `hashlib`, `pathlib`, `logging`, `threading`), pytest for tests.

---

## File Structure

```
hirn_literature_retrieve/
├── CLAUDE.md                        # Existing
├── SKILL.md                         # Skill definition for Claude Code
├── scripts/
│   ├── __init__.py
│   ├── scrape_hirn.py               # Fetch HIRN publications via AJAX
│   ├── resolve_ids.py               # PMID/DOI → PMCID resolution
│   ├── fetch_fulltext.py            # PMC full text via BioC JSON
│   ├── chunk_text.py                # Section/paragraph chunking
│   ├── search_chunks.py             # BM25 chunk ranking
│   └── utils/
│       ├── __init__.py
│       ├── rate_limiter.py          # Token bucket rate limiter
│       ├── cache_manager.py         # File-based cache with TTL
│       └── html_parser.py           # Parse HIRN publication HTML
├── tests/
│   ├── __init__.py
│   ├── conftest.py                  # Shared fixtures + sample data
│   ├── test_rate_limiter.py
│   ├── test_cache_manager.py
│   ├── test_html_parser.py
│   ├── test_scrape_hirn.py
│   ├── test_resolve_ids.py
│   ├── test_fetch_fulltext.py
│   ├── test_chunk_text.py
│   └── test_search_chunks.py
└── data/
    └── cache/                       # Runtime cache (gitignored)
```

---

## Task 1: Project Scaffolding

**Files:**
- Create: `scripts/__init__.py`
- Create: `scripts/utils/__init__.py`
- Create: `tests/__init__.py`
- Create: `tests/conftest.py`
- Create: `.gitignore`
- Create: `data/cache/.gitkeep`

**Step 1: Create directory structure**

```bash
mkdir -p scripts/utils tests data/cache
```

**Step 2: Create package init files**

`scripts/__init__.py`:
```python
"""HIRN Literature Retrieve — scripts package."""
```

`scripts/utils/__init__.py`:
```python
"""Utility modules for HIRN Literature Retrieve."""
```

`tests/__init__.py`:
```python
"""Tests for HIRN Literature Retrieve."""
```

**Step 3: Create .gitignore**

`.gitignore`:
```
__pycache__/
*.pyc
data/cache/*
!data/cache/.gitkeep
.pytest_cache/
*.egg-info/
```

**Step 4: Create conftest.py with sample HIRN HTML fixture**

`tests/conftest.py`:
```python
"""Shared test fixtures for HIRN Literature Retrieve."""
import pytest


SAMPLE_HIRN_HTML = """
<div class="row">
    <div class="col-sm-2">12/03/2025</div>
    <div class="col-sm-2">CBDS</div>
    <div class="col-sm-4">
        <p class="pub-title">
            <a href="https://pubmed.ncbi.nlm.nih.gov/39630627" target="_blank">
                Islet autoimmunity and HLA markers of presymptomatic and clinical type 1 diabetes.
            </a>
        </p>
        <p>
            Jacobsen LM, Bocchino LE, Evans-Molina C, et al.<br />
            Diabetologia. 2024 Dec;67(12):2611-2637.<br />
            doi: 10.1007/s00125-024-06244-0.<br />
            PMID:39630627
        </p>
    </div>
</div>
<div class="row">
    <div class="col-sm-2">11/15/2025</div>
    <div class="col-sm-2">HPAC, PanKbase</div>
    <div class="col-sm-4">
        <p class="pub-title">
            <a href="https://pubmed.ncbi.nlm.nih.gov/38012345" target="_blank">
                Single-cell transcriptomics of human pancreatic islets.
            </a>
        </p>
        <p>
            Smith J, Doe A, Johnson B.<br />
            Nature. 2024 Nov;615(7950):123-130.<br />
            doi: 10.1038/s41586-024-00001-1.<br />
            PMID:38012345
        </p>
    </div>
</div>
""".strip()


SAMPLE_AJAX_RESPONSE = {
    "html": SAMPLE_HIRN_HTML,
    "meta": {
        "postcount": 2,
        "totalposts": 2,
        "debug": False,
    },
}


SAMPLE_BIOC_JSON = {
    "source": "Auto-CuratedFull",
    "date": "20240101",
    "key": "autocuratedfull.key",
    "documents": [
        {
            "id": "39630627",
            "passages": [
                {
                    "infons": {"section_type": "TITLE", "type": "front"},
                    "offset": 0,
                    "text": "Islet autoimmunity and HLA markers of presymptomatic and clinical type 1 diabetes.",
                },
                {
                    "infons": {"section_type": "ABSTRACT", "type": "abstract"},
                    "offset": 83,
                    "text": "Type 1 diabetes is an autoimmune disease that destroys pancreatic beta cells. "
                    "Islet autoantibodies are the primary biomarkers for disease prediction. "
                    "HLA genotyping provides additional risk stratification.",
                },
                {
                    "infons": {"section_type": "INTRO", "type": "paragraph"},
                    "offset": 300,
                    "text": "The natural history of type 1 diabetes begins with genetic susceptibility, "
                    "progresses through islet autoimmunity, and culminates in clinical disease. "
                    "Understanding this progression is critical for intervention strategies.",
                },
                {
                    "infons": {"section_type": "METHODS", "type": "paragraph"},
                    "offset": 550,
                    "text": "We analyzed data from the TrialNet Pathway to Prevention study cohort. "
                    "Participants were screened for islet autoantibodies and HLA genotyped. "
                    "Follow-up extended to 15 years from initial autoantibody detection.",
                },
                {
                    "infons": {"section_type": "RESULTS", "type": "paragraph"},
                    "offset": 800,
                    "text": "Among 1,500 autoantibody-positive relatives, 45% progressed to clinical diabetes. "
                    "HLA-DR3/DR4-DQ8 carriers had the highest risk. "
                    "Multiple autoantibodies predicted faster progression.",
                },
                {
                    "infons": {"section_type": "DISCUSS", "type": "paragraph"},
                    "offset": 1050,
                    "text": "Our findings confirm that combined autoantibody and HLA profiling improves prediction "
                    "of type 1 diabetes progression. This has implications for clinical trial enrollment "
                    "and early intervention strategies.",
                },
            ],
        }
    ],
}


SAMPLE_ID_CONVERTER_RESPONSE = {
    "status": "ok",
    "responseDate": "2024-01-01 00:00:00",
    "request": "ids=39630627&format=json",
    "records": [
        {
            "pmid": "39630627",
            "pmcid": "PMC11615173",
            "doi": "10.1007/s00125-024-06244-0",
            "status": "ok",
        }
    ],
}


@pytest.fixture
def sample_hirn_html():
    return SAMPLE_HIRN_HTML


@pytest.fixture
def sample_ajax_response():
    return SAMPLE_AJAX_RESPONSE


@pytest.fixture
def sample_bioc_json():
    return SAMPLE_BIOC_JSON


@pytest.fixture
def sample_id_converter_response():
    return SAMPLE_ID_CONVERTER_RESPONSE


@pytest.fixture
def tmp_cache_dir(tmp_path):
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    return cache_dir
```

**Step 5: Create data/cache/.gitkeep**

```bash
touch data/cache/.gitkeep
```

**Step 6: Verify structure**

Run: `find . -type f | grep -v __pycache__ | sort`

Expected: all files listed above present.

**Step 7: Commit**

```bash
git init
git add .gitignore CLAUDE.md scripts/__init__.py scripts/utils/__init__.py tests/__init__.py tests/conftest.py data/cache/.gitkeep
git commit -m "feat: project scaffolding with test fixtures"
```

---

## Task 2: Rate Limiter

**Files:**
- Create: `scripts/utils/rate_limiter.py`
- Create: `tests/test_rate_limiter.py`

**Step 1: Write the failing tests**

`tests/test_rate_limiter.py`:
```python
"""Tests for token-bucket rate limiter."""
import time
from scripts.utils.rate_limiter import RateLimiter


def test_limiter_allows_burst_up_to_capacity():
    limiter = RateLimiter(rate=10.0, capacity=3)
    # Should allow 3 requests immediately (burst capacity)
    for _ in range(3):
        wait = limiter.wait()
        assert wait < 0.05  # near-instant


def test_limiter_throttles_after_burst():
    limiter = RateLimiter(rate=10.0, capacity=1)
    limiter.wait()  # consume the 1 token
    start = time.monotonic()
    limiter.wait()  # must wait ~0.1s for refill
    elapsed = time.monotonic() - start
    assert elapsed >= 0.08  # at least ~80ms


def test_limiter_default_from_env_no_key(monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    limiter = RateLimiter.default()
    assert limiter.rate == 3.0


def test_limiter_default_from_env_with_key(monkeypatch):
    monkeypatch.setenv("NCBI_API_KEY", "fake_key")
    limiter = RateLimiter.default()
    assert limiter.rate == 10.0
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_rate_limiter.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'scripts.utils.rate_limiter'`

**Step 3: Write minimal implementation**

`scripts/utils/rate_limiter.py`:
```python
"""Token-bucket rate limiter for NCBI API calls."""
from __future__ import annotations

import os
import threading
import time


class RateLimiter:
    """Thread-safe token-bucket rate limiter.

    Args:
        rate: Maximum sustained requests per second.
        capacity: Burst capacity (tokens available immediately).
    """

    def __init__(self, rate: float = 3.0, capacity: int = 1) -> None:
        self.rate = rate
        self.capacity = capacity
        self._tokens = float(capacity)
        self._last = time.monotonic()
        self._lock = threading.Lock()

    @classmethod
    def default(cls) -> RateLimiter:
        """Create a limiter configured from NCBI_API_KEY env var."""
        has_key = bool(os.environ.get("NCBI_API_KEY"))
        rate = 10.0 if has_key else 3.0
        return cls(rate=rate, capacity=1)

    def wait(self) -> float:
        """Block until a token is available. Returns seconds waited."""
        with self._lock:
            now = time.monotonic()
            elapsed = now - self._last
            self._tokens = min(self.capacity, self._tokens + elapsed * self.rate)
            self._last = now

            if self._tokens >= 1.0:
                self._tokens -= 1.0
                return 0.0

            deficit = 1.0 - self._tokens
            sleep_time = deficit / self.rate
            self._tokens = 0.0

        time.sleep(sleep_time)

        with self._lock:
            self._last = time.monotonic()

        return sleep_time
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_rate_limiter.py -v`
Expected: 4 passed

**Step 5: Commit**

```bash
git add scripts/utils/rate_limiter.py tests/test_rate_limiter.py
git commit -m "feat: token-bucket rate limiter for NCBI API"
```

---

## Task 3: Cache Manager

**Files:**
- Create: `scripts/utils/cache_manager.py`
- Create: `tests/test_cache_manager.py`

**Step 1: Write the failing tests**

`tests/test_cache_manager.py`:
```python
"""Tests for file-based cache manager."""
import json
import time
from scripts.utils.cache_manager import CacheManager


def test_cache_set_and_get(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    cache.set("test_key", {"data": "hello"}, ttl=60)
    result = cache.get("test_key")
    assert result == {"data": "hello"}


def test_cache_miss_returns_none(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    assert cache.get("nonexistent") is None


def test_cache_expired_returns_none(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    cache.set("expire_me", {"data": 1}, ttl=0)
    time.sleep(0.05)
    assert cache.get("expire_me") is None


def test_cache_cleanup_removes_expired(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    cache.set("old", {"data": 1}, ttl=0)
    cache.set("fresh", {"data": 2}, ttl=3600)
    time.sleep(0.05)
    removed = cache.cleanup()
    assert removed >= 1
    assert cache.get("old") is None
    assert cache.get("fresh") == {"data": 2}


def test_cache_key_hashing_is_deterministic(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    cache.set("same_key", {"v": 1}, ttl=60)
    cache.set("same_key", {"v": 2}, ttl=60)
    assert cache.get("same_key") == {"v": 2}
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_cache_manager.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write minimal implementation**

`scripts/utils/cache_manager.py`:
```python
"""File-based cache with TTL expiration."""
from __future__ import annotations

import hashlib
import json
import logging
import threading
import time
from pathlib import Path

logger = logging.getLogger(__name__)

# Default TTLs in seconds
TTL_INDEX = 86400       # 1 day  — HIRN publication index
TTL_METADATA = 2592000  # 30 days — article metadata
TTL_FULLTEXT = 2592000  # 30 days — article full text
TTL_IDS = 2592000       # 30 days — ID mappings


class CacheManager:
    """Thread-safe file-based JSON cache with TTL expiration.

    Args:
        cache_dir: Directory to store cache files.
    """

    def __init__(self, cache_dir: str | Path) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def _key_path(self, key: str) -> Path:
        hashed = hashlib.md5(key.encode()).hexdigest()
        return self.cache_dir / f"{hashed}.json"

    def get(self, key: str) -> dict | list | None:
        """Retrieve a cached value. Returns None on miss or expiry."""
        path = self._key_path(key)
        if not path.exists():
            return None
        try:
            with self._lock:
                data = json.loads(path.read_text())
            if time.time() > data["expires_at"]:
                path.unlink(missing_ok=True)
                return None
            return data["value"]
        except (json.JSONDecodeError, KeyError):
            path.unlink(missing_ok=True)
            return None

    def set(self, key: str, value: dict | list, ttl: int = 3600) -> None:
        """Store a value with a TTL in seconds."""
        path = self._key_path(key)
        envelope = {
            "key": key,
            "value": value,
            "created_at": time.time(),
            "expires_at": time.time() + ttl,
        }
        with self._lock:
            path.write_text(json.dumps(envelope))

    def cleanup(self) -> int:
        """Remove all expired entries. Returns count removed."""
        removed = 0
        for path in self.cache_dir.glob("*.json"):
            try:
                data = json.loads(path.read_text())
                if time.time() > data["expires_at"]:
                    path.unlink()
                    removed += 1
            except (json.JSONDecodeError, KeyError):
                path.unlink()
                removed += 1
        return removed
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_cache_manager.py -v`
Expected: 5 passed

**Step 5: Commit**

```bash
git add scripts/utils/cache_manager.py tests/test_cache_manager.py
git commit -m "feat: file-based cache manager with TTL expiration"
```

---

## Task 4: HTML Parser for HIRN Publications

**Files:**
- Create: `scripts/utils/html_parser.py`
- Create: `tests/test_html_parser.py`

**Step 1: Write the failing tests**

`tests/test_html_parser.py`:
```python
"""Tests for HIRN publication HTML parser."""
from scripts.utils.html_parser import parse_hirn_publications
from tests.conftest import SAMPLE_HIRN_HTML


def test_parse_returns_list():
    pubs = parse_hirn_publications(SAMPLE_HIRN_HTML)
    assert isinstance(pubs, list)
    assert len(pubs) == 2


def test_parse_extracts_pmid():
    pubs = parse_hirn_publications(SAMPLE_HIRN_HTML)
    assert pubs[0]["pmid"] == "39630627"
    assert pubs[1]["pmid"] == "38012345"


def test_parse_extracts_doi():
    pubs = parse_hirn_publications(SAMPLE_HIRN_HTML)
    assert pubs[0]["doi"] == "10.1007/s00125-024-06244-0"


def test_parse_extracts_title():
    pubs = parse_hirn_publications(SAMPLE_HIRN_HTML)
    assert "Islet autoimmunity" in pubs[0]["title"]


def test_parse_extracts_authors():
    pubs = parse_hirn_publications(SAMPLE_HIRN_HTML)
    assert "Jacobsen LM" in pubs[0]["authors"]


def test_parse_extracts_consortia():
    pubs = parse_hirn_publications(SAMPLE_HIRN_HTML)
    assert pubs[0]["consortia"] == ["CBDS"]
    assert pubs[1]["consortia"] == ["HPAC", "PanKbase"]


def test_parse_extracts_date():
    pubs = parse_hirn_publications(SAMPLE_HIRN_HTML)
    assert pubs[0]["date"] == "12/03/2025"


def test_parse_extracts_journal():
    pubs = parse_hirn_publications(SAMPLE_HIRN_HTML)
    assert "Diabetologia" in pubs[0]["journal"]


def test_parse_empty_html():
    pubs = parse_hirn_publications("")
    assert pubs == []
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_html_parser.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write minimal implementation**

`scripts/utils/html_parser.py`:
```python
"""Parse HIRN publication entries from AJAX HTML responses."""
from __future__ import annotations

import re
from html.parser import HTMLParser


class _HIRNParser(HTMLParser):
    """State-machine HTML parser for HIRN publication rows."""

    def __init__(self) -> None:
        super().__init__()
        self.publications: list[dict] = []
        self._in_row = False
        self._col_index = 0
        self._in_pub_title = False
        self._in_link = False
        self._in_citation_p = False
        self._current: dict = {}
        self._text_buf: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr_dict = dict(attrs)
        classes = (attr_dict.get("class") or "").split()

        if tag == "div" and "row" in classes:
            self._in_row = True
            self._col_index = 0
            self._current = {
                "pmid": "",
                "doi": "",
                "title": "",
                "authors": "",
                "journal": "",
                "date": "",
                "consortia": [],
            }

        elif tag == "div" and self._in_row and "col-sm-2" in classes:
            self._col_index += 1
            self._text_buf = []

        elif tag == "div" and self._in_row and "col-sm-4" in classes:
            self._col_index = 3
            self._text_buf = []

        elif tag == "p" and self._in_row and "pub-title" in classes:
            self._in_pub_title = True
            self._text_buf = []

        elif tag == "a" and self._in_pub_title:
            self._in_link = True
            href = attr_dict.get("href", "")
            pmid_match = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", href)
            if pmid_match:
                self._current["pmid"] = pmid_match.group(1)

        elif tag == "p" and self._in_row and self._col_index == 3 and not self._in_pub_title:
            self._in_citation_p = True
            self._text_buf = []

    def handle_endtag(self, tag: str) -> None:
        if tag == "a" and self._in_link:
            self._in_link = False
            self._current["title"] = "".join(self._text_buf).strip()
            self._text_buf = []

        elif tag == "p" and self._in_pub_title:
            self._in_pub_title = False

        elif tag == "p" and self._in_citation_p:
            self._in_citation_p = False
            citation_text = "".join(self._text_buf).strip()
            self._parse_citation(citation_text)

        elif tag == "div" and self._in_row:
            if self._col_index == 1:
                self._current["date"] = "".join(self._text_buf).strip()
            elif self._col_index == 2:
                raw = "".join(self._text_buf).strip()
                self._current["consortia"] = [
                    c.strip() for c in raw.split(",") if c.strip()
                ]

        # Detect end of a row: the second top-level </div> after a row start
        # We rely on the fact that the row div closes after all its children
        if tag == "div" and self._in_row and self._current.get("pmid"):
            # Check if we've accumulated enough data
            if self._current.get("title"):
                self.publications.append(self._current.copy())
                self._in_row = False

    def handle_data(self, data: str) -> None:
        if self._in_link or self._in_citation_p:
            self._text_buf.append(data)
        elif self._in_row and self._col_index in (1, 2):
            self._text_buf.append(data)

    def _parse_citation(self, text: str) -> None:
        # Extract DOI
        doi_match = re.search(r"doi:\s*(10\.\S+?)[\s.]", text)
        if doi_match:
            self._current["doi"] = doi_match.group(1).rstrip(".")

        # Extract PMID (backup if not from href)
        pmid_match = re.search(r"PMID:\s*(\d+)", text)
        if pmid_match and not self._current["pmid"]:
            self._current["pmid"] = pmid_match.group(1)

        # Extract authors (text before first journal-like pattern)
        # Authors end at first period followed by journal name
        parts = text.split(".\n")
        if not parts:
            parts = text.split(". ")
        if parts:
            self._current["authors"] = parts[0].strip().rstrip(".")

        # Extract journal (after authors, before year)
        journal_match = re.search(
            r"\.\s*([A-Z][^.]+?)\.\s*\d{4}", text
        )
        if journal_match:
            self._current["journal"] = journal_match.group(1).strip()


def parse_hirn_publications(html: str) -> list[dict]:
    """Parse HIRN AJAX HTML response into publication dicts.

    Args:
        html: Raw HTML string from HIRN AJAX endpoint.

    Returns:
        List of dicts with keys: pmid, doi, title, authors,
        journal, date, consortia.
    """
    if not html or not html.strip():
        return []
    parser = _HIRNParser()
    parser.feed(html)
    return parser.publications
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_html_parser.py -v`
Expected: 9 passed

**Step 5: Commit**

```bash
git add scripts/utils/html_parser.py tests/test_html_parser.py
git commit -m "feat: HTML parser for HIRN publication entries"
```

---

## Task 5: HIRN Scraper

**Files:**
- Create: `scripts/scrape_hirn.py`
- Create: `tests/test_scrape_hirn.py`

**Step 1: Write the failing tests**

`tests/test_scrape_hirn.py`:
```python
"""Tests for HIRN publications scraper."""
import json
from unittest.mock import patch, MagicMock
from scripts.scrape_hirn import fetch_hirn_publications, search_publications
from tests.conftest import SAMPLE_AJAX_RESPONSE


def _mock_urlopen(response_data):
    """Create a mock for urllib.request.urlopen."""
    mock_resp = MagicMock()
    mock_resp.read.return_value = json.dumps(response_data).encode()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    return mock_resp


@patch("scripts.scrape_hirn.urlopen")
def test_fetch_returns_publications(mock_urlopen, tmp_cache_dir):
    # Simulate single page (totalposts == postcount)
    mock_urlopen.return_value = _mock_urlopen(SAMPLE_AJAX_RESPONSE)
    pubs = fetch_hirn_publications(cache_dir=str(tmp_cache_dir))
    assert len(pubs) == 2
    assert pubs[0]["pmid"] == "39630627"


@patch("scripts.scrape_hirn.urlopen")
def test_fetch_caches_result(mock_urlopen, tmp_cache_dir):
    mock_urlopen.return_value = _mock_urlopen(SAMPLE_AJAX_RESPONSE)
    # First call fetches from network
    pubs1 = fetch_hirn_publications(cache_dir=str(tmp_cache_dir))
    # Second call should use cache (no new urlopen call)
    pubs2 = fetch_hirn_publications(cache_dir=str(tmp_cache_dir))
    assert pubs1 == pubs2
    assert mock_urlopen.call_count == 1


def test_search_publications_by_keyword():
    pubs = [
        {"pmid": "1", "title": "Islet autoimmunity in type 1 diabetes", "consortia": ["CBDS"]},
        {"pmid": "2", "title": "Single-cell transcriptomics of pancreas", "consortia": ["HPAC"]},
        {"pmid": "3", "title": "Beta cell regeneration methods", "consortia": ["CTAR"]},
    ]
    results = search_publications(pubs, query="islet autoimmunity")
    assert len(results) > 0
    assert results[0]["pmid"] == "1"  # Best match first


def test_search_publications_filter_consortium():
    pubs = [
        {"pmid": "1", "title": "Islet autoimmunity", "consortia": ["CBDS"]},
        {"pmid": "2", "title": "Islet function", "consortia": ["HPAC"]},
    ]
    results = search_publications(pubs, query="islet", consortium="HPAC")
    assert len(results) == 1
    assert results[0]["pmid"] == "2"


def test_search_publications_no_match():
    pubs = [{"pmid": "1", "title": "Beta cell study", "consortia": ["CBDS"]}]
    results = search_publications(pubs, query="quantum computing")
    assert results == []
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_scrape_hirn.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write minimal implementation**

`scripts/scrape_hirn.py`:
```python
"""Fetch and search HIRN published works via AJAX endpoint."""
from __future__ import annotations

import json
import logging
import re
import urllib.parse
import urllib.request
from urllib.request import urlopen

from scripts.utils.cache_manager import CacheManager, TTL_INDEX
from scripts.utils.html_parser import parse_hirn_publications
from scripts.utils.rate_limiter import RateLimiter

logger = logging.getLogger(__name__)

AJAX_URL = "https://hirnetwork.org/2021/wp-admin/admin-ajax.php"
POSTS_PER_PAGE = 250
DEFAULT_CACHE_DIR = "data/cache"


def fetch_hirn_publications(
    cache_dir: str = DEFAULT_CACHE_DIR,
    consortium: str | None = None,
) -> list[dict]:
    """Fetch all HIRN publications, paginating through AJAX endpoint.

    Args:
        cache_dir: Directory for file cache.
        consortium: Optional consortium slug to filter (e.g. 'cbds').

    Returns:
        List of publication dicts with keys: pmid, doi, title,
        authors, journal, date, consortia.
    """
    cache = CacheManager(cache_dir)
    cache_key = f"hirn_index_{consortium or 'all'}"
    cached = cache.get(cache_key)
    if cached is not None:
        logger.info("Using cached HIRN publication index")
        return cached

    limiter = RateLimiter(rate=3.0, capacity=1)  # gentle on WordPress
    all_pubs: list[dict] = []
    page = 0

    while True:
        params = {
            "action": "alm_get_posts",
            "post_type": "publication",
            "posts_per_page": str(POSTS_PER_PAGE),
            "page": str(page),
        }
        if consortium:
            params["taxonomy"] = "publication_consortium"
            params["taxonomy_terms"] = consortium.lower()
        else:
            params["taxonomy"] = "publication_category"
            params["taxonomy_terms"] = "published_works"

        url = f"{AJAX_URL}?{urllib.parse.urlencode(params)}"
        limiter.wait()

        logger.info("Fetching HIRN publications page %d", page)
        with urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))

        html = data.get("html", "")
        meta = data.get("meta", {})
        pubs = parse_hirn_publications(html)
        all_pubs.extend(pubs)

        total = meta.get("totalposts", 0)
        if not pubs or len(all_pubs) >= total:
            break
        page += 1

    logger.info("Fetched %d HIRN publications", len(all_pubs))
    cache.set(cache_key, all_pubs, ttl=TTL_INDEX)
    return all_pubs


def search_publications(
    publications: list[dict],
    query: str,
    consortium: str | None = None,
    max_results: int = 20,
) -> list[dict]:
    """Search publications by keyword matching on titles.

    Args:
        publications: List of publication dicts.
        query: Search query string.
        consortium: Optional consortium filter.
        max_results: Maximum results to return.

    Returns:
        List of matching publications, sorted by relevance (best first).
    """
    if consortium:
        publications = [
            p for p in publications
            if consortium.upper() in [c.upper() for c in p.get("consortia", [])]
        ]

    query_terms = re.findall(r"\w+", query.lower())
    if not query_terms:
        return []

    scored: list[tuple[float, dict]] = []
    for pub in publications:
        title_lower = pub.get("title", "").lower()
        # Score: fraction of query terms found in title
        matches = sum(1 for t in query_terms if t in title_lower)
        if matches > 0:
            score = matches / len(query_terms)
            scored.append((score, pub))

    scored.sort(key=lambda x: x[0], reverse=True)
    return [pub for _, pub in scored[:max_results]]
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_scrape_hirn.py -v`
Expected: 5 passed

**Step 5: Commit**

```bash
git add scripts/scrape_hirn.py tests/test_scrape_hirn.py
git commit -m "feat: HIRN AJAX scraper with keyword search"
```

---

## Task 6: ID Resolver (PMID/DOI → PMCID)

**Files:**
- Create: `scripts/resolve_ids.py`
- Create: `tests/test_resolve_ids.py`

**Step 1: Write the failing tests**

`tests/test_resolve_ids.py`:
```python
"""Tests for NCBI ID resolver."""
import json
from unittest.mock import patch, MagicMock
from scripts.resolve_ids import resolve_pmcids
from tests.conftest import SAMPLE_ID_CONVERTER_RESPONSE


def _mock_urlopen(response_data):
    mock_resp = MagicMock()
    mock_resp.read.return_value = json.dumps(response_data).encode()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    return mock_resp


@patch("scripts.resolve_ids.urlopen")
def test_resolve_single_pmid(mock_urlopen, tmp_cache_dir):
    mock_urlopen.return_value = _mock_urlopen(SAMPLE_ID_CONVERTER_RESPONSE)
    result = resolve_pmcids(["39630627"], cache_dir=str(tmp_cache_dir))
    assert result["39630627"] == "PMC11615173"


@patch("scripts.resolve_ids.urlopen")
def test_resolve_missing_pmcid(mock_urlopen, tmp_cache_dir):
    response = {
        "status": "ok",
        "records": [
            {"pmid": "99999999", "status": "error", "errmsg": "not found"}
        ],
    }
    mock_urlopen.return_value = _mock_urlopen(response)
    result = resolve_pmcids(["99999999"], cache_dir=str(tmp_cache_dir))
    assert result.get("99999999") is None


@patch("scripts.resolve_ids.urlopen")
def test_resolve_batches_large_lists(mock_urlopen, tmp_cache_dir):
    # 250 PMIDs should be split into batches of 200
    response = {"status": "ok", "records": []}
    mock_urlopen.return_value = _mock_urlopen(response)
    pmids = [str(i) for i in range(250)]
    resolve_pmcids(pmids, cache_dir=str(tmp_cache_dir))
    assert mock_urlopen.call_count == 2  # ceil(250/200) = 2


@patch("scripts.resolve_ids.urlopen")
def test_resolve_uses_cache(mock_urlopen, tmp_cache_dir):
    mock_urlopen.return_value = _mock_urlopen(SAMPLE_ID_CONVERTER_RESPONSE)
    resolve_pmcids(["39630627"], cache_dir=str(tmp_cache_dir))
    resolve_pmcids(["39630627"], cache_dir=str(tmp_cache_dir))
    assert mock_urlopen.call_count == 1  # only 1 network call
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_resolve_ids.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write minimal implementation**

`scripts/resolve_ids.py`:
```python
"""Resolve PMIDs/DOIs to PMCIDs via NCBI ID Converter API."""
from __future__ import annotations

import json
import logging
import urllib.parse
from urllib.request import urlopen

from scripts.utils.cache_manager import CacheManager, TTL_IDS
from scripts.utils.rate_limiter import RateLimiter

logger = logging.getLogger(__name__)

ID_CONVERTER_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
BATCH_SIZE = 200
DEFAULT_CACHE_DIR = "data/cache"


def resolve_pmcids(
    pmids: list[str],
    cache_dir: str = DEFAULT_CACHE_DIR,
) -> dict[str, str | None]:
    """Resolve a list of PMIDs to PMCIDs.

    Args:
        pmids: List of PMID strings.
        cache_dir: Directory for file cache.

    Returns:
        Dict mapping PMID → PMCID (or None if not in PMC).
    """
    cache = CacheManager(cache_dir)
    limiter = RateLimiter.default()
    result: dict[str, str | None] = {}
    uncached: list[str] = []

    # Check cache first
    for pmid in pmids:
        cached = cache.get(f"pmcid_{pmid}")
        if cached is not None:
            result[pmid] = cached.get("pmcid")
        else:
            uncached.append(pmid)

    if not uncached:
        return result

    # Batch resolve uncached PMIDs
    for i in range(0, len(uncached), BATCH_SIZE):
        batch = uncached[i : i + BATCH_SIZE]
        params = {
            "ids": ",".join(batch),
            "format": "json",
            "tool": "hirn-literature-retrieve",
            "email": "tool@example.com",
        }
        url = f"{ID_CONVERTER_URL}?{urllib.parse.urlencode(params)}"
        limiter.wait()

        logger.info("Resolving PMCIDs for %d PMIDs", len(batch))
        with urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))

        for record in data.get("records", []):
            pmid = record.get("pmid", "")
            pmcid = record.get("pmcid") if record.get("status") == "ok" else None
            result[pmid] = pmcid
            cache.set(f"pmcid_{pmid}", {"pmcid": pmcid}, ttl=TTL_IDS)

    return result
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_resolve_ids.py -v`
Expected: 4 passed

**Step 5: Commit**

```bash
git add scripts/resolve_ids.py tests/test_resolve_ids.py
git commit -m "feat: NCBI ID converter for PMID to PMCID resolution"
```

---

## Task 7: Full Text Fetcher

**Files:**
- Create: `scripts/fetch_fulltext.py`
- Create: `tests/test_fetch_fulltext.py`

**Step 1: Write the failing tests**

`tests/test_fetch_fulltext.py`:
```python
"""Tests for PMC full text fetcher."""
import json
from unittest.mock import patch, MagicMock
from scripts.fetch_fulltext import fetch_fulltext, parse_bioc_passages
from tests.conftest import SAMPLE_BIOC_JSON


def _mock_urlopen(response_data):
    mock_resp = MagicMock()
    mock_resp.read.return_value = json.dumps(response_data).encode()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    return mock_resp


def test_parse_bioc_passages():
    passages = parse_bioc_passages(SAMPLE_BIOC_JSON)
    assert len(passages) == 6
    assert passages[0]["section"] == "TITLE"
    assert "Islet autoimmunity" in passages[0]["text"]
    assert passages[1]["section"] == "ABSTRACT"


@patch("scripts.fetch_fulltext.urlopen")
def test_fetch_fulltext_success(mock_urlopen, tmp_cache_dir):
    mock_urlopen.return_value = _mock_urlopen(SAMPLE_BIOC_JSON)
    result = fetch_fulltext("PMC11615173", cache_dir=str(tmp_cache_dir))
    assert result["success"] is True
    assert len(result["passages"]) == 6
    assert result["pmcid"] == "PMC11615173"


@patch("scripts.fetch_fulltext.urlopen")
def test_fetch_fulltext_caches(mock_urlopen, tmp_cache_dir):
    mock_urlopen.return_value = _mock_urlopen(SAMPLE_BIOC_JSON)
    fetch_fulltext("PMC11615173", cache_dir=str(tmp_cache_dir))
    fetch_fulltext("PMC11615173", cache_dir=str(tmp_cache_dir))
    assert mock_urlopen.call_count == 1


@patch("scripts.fetch_fulltext.urlopen")
def test_fetch_fulltext_not_found(mock_urlopen, tmp_cache_dir):
    import urllib.error
    mock_urlopen.side_effect = urllib.error.HTTPError(
        url="", code=404, msg="Not Found", hdrs=None, fp=None
    )
    result = fetch_fulltext("PMC00000000", cache_dir=str(tmp_cache_dir))
    assert result["success"] is False
    assert result["error"]["code"] == "NOT_FOUND"
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_fetch_fulltext.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write minimal implementation**

`scripts/fetch_fulltext.py`:
```python
"""Fetch full text from PMC Open Access via BioC JSON API."""
from __future__ import annotations

import json
import logging
import urllib.error
from urllib.request import urlopen

from scripts.utils.cache_manager import CacheManager, TTL_FULLTEXT
from scripts.utils.rate_limiter import RateLimiter

logger = logging.getLogger(__name__)

BIOC_URL = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{pmcid}/unicode"
DEFAULT_CACHE_DIR = "data/cache"


def parse_bioc_passages(bioc_data: dict) -> list[dict]:
    """Extract passages from BioC JSON response.

    Args:
        bioc_data: Parsed BioC JSON dict.

    Returns:
        List of passage dicts with keys: section, type, offset, text.
    """
    passages = []
    for doc in bioc_data.get("documents", []):
        for passage in doc.get("passages", []):
            infons = passage.get("infons", {})
            passages.append({
                "section": infons.get("section_type", "UNKNOWN"),
                "type": infons.get("type", "unknown"),
                "offset": passage.get("offset", 0),
                "text": passage.get("text", ""),
            })
    return passages


def fetch_fulltext(
    pmcid: str,
    cache_dir: str = DEFAULT_CACHE_DIR,
) -> dict:
    """Fetch full text for a PMC article.

    Args:
        pmcid: PMC ID (e.g. 'PMC11615173').
        cache_dir: Directory for file cache.

    Returns:
        Dict with keys: success, pmcid, passages, error.
    """
    cache = CacheManager(cache_dir)
    cache_key = f"fulltext_{pmcid}"
    cached = cache.get(cache_key)
    if cached is not None:
        logger.info("Cache hit for full text: %s", pmcid)
        return cached

    limiter = RateLimiter.default()
    url = BIOC_URL.format(pmcid=pmcid)
    limiter.wait()

    try:
        logger.info("Fetching full text for %s", pmcid)
        with urlopen(url, timeout=60) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return {
                "success": False,
                "pmcid": pmcid,
                "passages": [],
                "error": {
                    "code": "NOT_FOUND",
                    "message": f"Full text not available in PMC Open Access for {pmcid}",
                },
            }
        return {
            "success": False,
            "pmcid": pmcid,
            "passages": [],
            "error": {"code": "API_ERROR", "message": str(e)},
        }
    except urllib.error.URLError as e:
        return {
            "success": False,
            "pmcid": pmcid,
            "passages": [],
            "error": {"code": "NETWORK_ERROR", "message": str(e)},
        }

    passages = parse_bioc_passages(data)
    result = {"success": True, "pmcid": pmcid, "passages": passages}
    cache.set(cache_key, result, ttl=TTL_FULLTEXT)
    return result
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_fetch_fulltext.py -v`
Expected: 4 passed

**Step 5: Commit**

```bash
git add scripts/fetch_fulltext.py tests/test_fetch_fulltext.py
git commit -m "feat: PMC full text fetcher via BioC JSON API"
```

---

## Task 8: Text Chunker

**Files:**
- Create: `scripts/chunk_text.py`
- Create: `tests/test_chunk_text.py`

**Step 1: Write the failing tests**

`tests/test_chunk_text.py`:
```python
"""Tests for text chunking."""
from scripts.chunk_text import chunk_passages


def test_chunk_short_passages_stay_intact():
    passages = [
        {"section": "INTRO", "text": "Short paragraph.", "offset": 0, "type": "paragraph"},
    ]
    chunks = chunk_passages(passages, max_chunk_size=500)
    assert len(chunks) == 1
    assert chunks[0]["text"] == "Short paragraph."
    assert chunks[0]["section"] == "INTRO"


def test_chunk_long_passage_splits():
    long_text = "word " * 200  # ~1000 chars
    passages = [
        {"section": "RESULTS", "text": long_text.strip(), "offset": 0, "type": "paragraph"},
    ]
    chunks = chunk_passages(passages, max_chunk_size=300, overlap=50)
    assert len(chunks) > 1
    # Each chunk should be <= max_chunk_size
    for c in chunks:
        assert len(c["text"]) <= 350  # allow small overflow at word boundaries


def test_chunk_preserves_section_metadata():
    passages = [
        {"section": "METHODS", "text": "We used method A.", "offset": 100, "type": "paragraph"},
        {"section": "RESULTS", "text": "We found result B.", "offset": 200, "type": "paragraph"},
    ]
    chunks = chunk_passages(passages)
    assert chunks[0]["section"] == "METHODS"
    assert chunks[1]["section"] == "RESULTS"


def test_chunk_assigns_index():
    passages = [
        {"section": "INTRO", "text": "Para one.", "offset": 0, "type": "paragraph"},
        {"section": "INTRO", "text": "Para two.", "offset": 50, "type": "paragraph"},
    ]
    chunks = chunk_passages(passages)
    assert chunks[0]["chunk_index"] == 0
    assert chunks[1]["chunk_index"] == 1


def test_chunk_empty_passages():
    assert chunk_passages([]) == []


def test_chunk_skips_empty_text():
    passages = [
        {"section": "TITLE", "text": "", "offset": 0, "type": "front"},
        {"section": "INTRO", "text": "Content here.", "offset": 10, "type": "paragraph"},
    ]
    chunks = chunk_passages(passages)
    assert len(chunks) == 1
    assert chunks[0]["section"] == "INTRO"
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_chunk_text.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write minimal implementation**

`scripts/chunk_text.py`:
```python
"""Chunk article passages into retrievable text segments."""
from __future__ import annotations


def chunk_passages(
    passages: list[dict],
    max_chunk_size: int = 500,
    overlap: int = 50,
) -> list[dict]:
    """Split passages into chunks suitable for retrieval.

    Short passages are kept intact. Long passages are split at word
    boundaries with configurable overlap.

    Args:
        passages: List of passage dicts from parse_bioc_passages().
        max_chunk_size: Maximum characters per chunk.
        overlap: Character overlap between consecutive chunks of a split passage.

    Returns:
        List of chunk dicts with keys: text, section, type, offset,
        chunk_index.
    """
    chunks: list[dict] = []
    idx = 0

    for passage in passages:
        text = passage.get("text", "").strip()
        if not text:
            continue

        section = passage.get("section", "UNKNOWN")
        ptype = passage.get("type", "unknown")
        offset = passage.get("offset", 0)

        if len(text) <= max_chunk_size:
            chunks.append({
                "text": text,
                "section": section,
                "type": ptype,
                "offset": offset,
                "chunk_index": idx,
            })
            idx += 1
        else:
            # Split long passage at word boundaries
            start = 0
            while start < len(text):
                end = start + max_chunk_size
                if end < len(text):
                    # Find last space before end
                    space_pos = text.rfind(" ", start, end)
                    if space_pos > start:
                        end = space_pos
                chunk_text = text[start:end].strip()
                if chunk_text:
                    chunks.append({
                        "text": chunk_text,
                        "section": section,
                        "type": ptype,
                        "offset": offset + start,
                        "chunk_index": idx,
                    })
                    idx += 1
                start = end - overlap if end < len(text) else len(text)

    return chunks
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_chunk_text.py -v`
Expected: 6 passed

**Step 5: Commit**

```bash
git add scripts/chunk_text.py tests/test_chunk_text.py
git commit -m "feat: text chunker with overlap for passage retrieval"
```

---

## Task 9: Chunk Searcher (BM25 Ranking)

**Files:**
- Create: `scripts/search_chunks.py`
- Create: `tests/test_search_chunks.py`

**Step 1: Write the failing tests**

`tests/test_search_chunks.py`:
```python
"""Tests for BM25 chunk search/ranking."""
from scripts.search_chunks import search_chunks


def test_search_ranks_relevant_first():
    chunks = [
        {"text": "The weather is sunny today.", "section": "INTRO", "chunk_index": 0},
        {"text": "Islet autoantibodies predict type 1 diabetes progression.", "section": "RESULTS", "chunk_index": 1},
        {"text": "Beta cells produce insulin in the pancreas.", "section": "INTRO", "chunk_index": 2},
    ]
    results = search_chunks(chunks, query="islet autoantibodies diabetes")
    assert len(results) > 0
    assert results[0]["chunk_index"] == 1  # most relevant


def test_search_returns_top_k():
    chunks = [
        {"text": f"Diabetes topic sentence number {i}.", "section": "BODY", "chunk_index": i}
        for i in range(20)
    ]
    results = search_chunks(chunks, query="diabetes", top_k=5)
    assert len(results) <= 5


def test_search_no_match():
    chunks = [
        {"text": "The sun rises in the east.", "section": "INTRO", "chunk_index": 0},
    ]
    results = search_chunks(chunks, query="quantum entanglement")
    assert results == []


def test_search_includes_score():
    chunks = [
        {"text": "Islet autoimmunity is a marker.", "section": "RESULTS", "chunk_index": 0},
    ]
    results = search_chunks(chunks, query="islet autoimmunity")
    assert "score" in results[0]
    assert results[0]["score"] > 0


def test_search_empty_input():
    assert search_chunks([], query="anything") == []
    chunks = [{"text": "Some text.", "section": "INTRO", "chunk_index": 0}]
    assert search_chunks(chunks, query="") == []
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_search_chunks.py -v`
Expected: FAIL — `ModuleNotFoundError`

**Step 3: Write minimal implementation**

`scripts/search_chunks.py`:
```python
"""BM25 search and ranking for text chunks."""
from __future__ import annotations

import math
import re


def _tokenize(text: str) -> list[str]:
    """Lowercase word tokenization."""
    return re.findall(r"\w+", text.lower())


def search_chunks(
    chunks: list[dict],
    query: str,
    top_k: int = 10,
    k1: float = 1.5,
    b: float = 0.75,
) -> list[dict]:
    """Rank chunks by BM25 relevance to a query.

    Args:
        chunks: List of chunk dicts (must have 'text' key).
        query: Search query string.
        top_k: Number of top results to return.
        k1: BM25 term frequency saturation parameter.
        b: BM25 length normalization parameter.

    Returns:
        Top-k chunks sorted by score (descending), each with added 'score' key.
    """
    if not chunks or not query or not query.strip():
        return []

    query_terms = _tokenize(query)
    if not query_terms:
        return []

    # Tokenize all chunks
    doc_tokens = [_tokenize(c.get("text", "")) for c in chunks]
    n = len(doc_tokens)
    avgdl = sum(len(d) for d in doc_tokens) / n if n else 1

    # Document frequency for each query term
    df: dict[str, int] = {}
    for term in set(query_terms):
        df[term] = sum(1 for d in doc_tokens if term in d)

    # Score each chunk
    scored: list[tuple[float, int]] = []
    for i, tokens in enumerate(doc_tokens):
        score = 0.0
        dl = len(tokens)
        tf_map: dict[str, int] = {}
        for t in tokens:
            tf_map[t] = tf_map.get(t, 0) + 1

        for term in query_terms:
            if term not in tf_map:
                continue
            tf = tf_map[term]
            idf = math.log((n - df[term] + 0.5) / (df[term] + 0.5) + 1.0)
            tf_norm = (tf * (k1 + 1)) / (tf + k1 * (1 - b + b * dl / avgdl))
            score += idf * tf_norm

        if score > 0:
            scored.append((score, i))

    scored.sort(key=lambda x: x[0], reverse=True)
    results = []
    for score, i in scored[:top_k]:
        chunk = dict(chunks[i])
        chunk["score"] = round(score, 4)
        results.append(chunk)

    return results
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_search_chunks.py -v`
Expected: 5 passed

**Step 5: Commit**

```bash
git add scripts/search_chunks.py tests/test_search_chunks.py
git commit -m "feat: BM25 chunk search and ranking"
```

---

## Task 10: SKILL.md — Skill Definition

**Files:**
- Create: `SKILL.md`

**Step 1: Write SKILL.md**

`SKILL.md`:
```markdown
---
name: hirn-literature-retrieve
description: Search HIRN (Human Islet Research Network) published articles, retrieve metadata and full text from PubMed Central, and locate relevant text chunks. Covers ~1,160 publications across 9 HIRN consortia (CBDS, CHIB, CMAD, CMAI, CTAR, HIREC, HPAC, Opportunity Pool, PanKbase).
---

# HIRN Literature Retrieve

Search and retrieve relevant passages from articles published by the Human Islet Research Network (HIRN).

## When to Use This Skill

Use this skill when:
- Searching for HIRN-published research on islet biology, type 1 diabetes, beta cells, or related topics
- Retrieving full text and relevant passages from HIRN consortium publications
- Looking for HIRN research by consortium (CBDS, CHIB, CMAD, CMAI, CTAR, HIREC, HPAC, Opportunity Pool, PanKbase)
- Needing specific passages or evidence from HIRN literature to support analysis

### Keywords That Trigger Activation

**HIRN Keywords**: HIRN, Human Islet Research Network, HIRN publications, HIRN published works, islet research network

**Consortium Keywords**: CBDS, CHIB, CMAD, CMAI, CTAR, HIREC, HPAC, PanKbase, Opportunity Pool

**Topic Keywords**: islet autoimmunity, beta cell, pancreatic islet, type 1 diabetes HIRN, islet transplantation

## How It Works

### Pipeline Overview

```
User Query
    │
    ▼
┌────────────────────┐
│ 1. Fetch HIRN Index │ ──▶ AJAX scrape of hirnetwork.org/published-works
└────────────────────┘     (cached for 24 hours)
    │
    ▼
┌────────────────────┐
│ 2. Search Titles    │ ──▶ Keyword matching on publication titles
└────────────────────┘     Optional consortium filter
    │
    ▼
┌────────────────────┐
│ 3. Resolve PMCIDs   │ ──▶ NCBI ID Converter API (PMID → PMCID)
└────────────────────┘     (cached for 30 days)
    │
    ▼
┌────────────────────┐
│ 4. Fetch Full Text  │ ──▶ PMC BioC JSON API (Open Access only)
└────────────────────┘     (cached for 30 days)
    │
    ▼
┌────────────────────┐
│ 5. Chunk & Search   │ ──▶ Section/paragraph chunking + BM25 ranking
└────────────────────┘
    │
    ▼
  Relevant passages with citations
```

### Rate Limiting

All NCBI API calls are rate-limited:
- Without `NCBI_API_KEY`: 3 requests/second
- With `NCBI_API_KEY` env var: 10 requests/second
- HIRN WordPress AJAX: 3 requests/second (gentle)

### Caching

File-based cache in `data/cache/` with TTLs:
- HIRN publication index: 24 hours
- PMCID mappings: 30 days
- Full text: 30 days

## Workflow: Search and Retrieve HIRN Literature

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

This fetches all ~1,160 HIRN publications. Cached for 24 hours after first call.

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

To filter by consortium, add `consortium='CBDS'` (or any consortium name).

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

- **PMC Open Access only**: Full text is available only for articles in the PMC Open Access subset (~50-60% of HIRN publications). For others, only title/metadata from the HIRN listing is available.
- **Title-based search**: Initial publication matching uses title keywords only (abstracts are not in the HIRN listing). For more precise searches, use the full-text chunk search on fetched articles.
- **Rate limits**: NCBI APIs are rate-limited. The skill handles this automatically but large batch operations may take time.

## Environment Variables

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `NCBI_API_KEY` | No | None | NCBI API key for 10 req/sec (vs 3) |
```

**Step 2: Verify SKILL.md is valid YAML frontmatter**

Run: `python -c "import yaml; yaml.safe_load(open('SKILL.md').read().split('---')[1]); print('Valid')"`

If `yaml` not available, visually verify the frontmatter is correct.

**Step 3: Commit**

```bash
git add SKILL.md
git commit -m "feat: SKILL.md skill definition for HIRN literature retrieve"
```

---

## Task 11: Integration Test

**Files:**
- Create: `tests/test_integration.py`

**Step 1: Write the integration test**

`tests/test_integration.py`:
```python
"""Integration test: end-to-end pipeline with mocked HTTP."""
import json
from unittest.mock import patch, MagicMock, call
from tests.conftest import (
    SAMPLE_AJAX_RESPONSE,
    SAMPLE_ID_CONVERTER_RESPONSE,
    SAMPLE_BIOC_JSON,
)
from scripts.scrape_hirn import fetch_hirn_publications, search_publications
from scripts.resolve_ids import resolve_pmcids
from scripts.fetch_fulltext import fetch_fulltext
from scripts.chunk_text import chunk_passages
from scripts.search_chunks import search_chunks


def _mock_urlopen_factory(responses):
    """Return a side_effect function that returns different responses per call."""
    call_count = [0]

    def side_effect(url, **kwargs):
        mock_resp = MagicMock()
        idx = min(call_count[0], len(responses) - 1)
        mock_resp.read.return_value = json.dumps(responses[idx]).encode()
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)
        call_count[0] += 1
        return mock_resp

    return side_effect


def test_end_to_end_pipeline(tmp_cache_dir):
    """Full pipeline: scrape → search → resolve → fetch → chunk → search."""
    cache_dir = str(tmp_cache_dir)

    # Mock all HTTP calls
    with patch("scripts.scrape_hirn.urlopen") as mock_scrape, \
         patch("scripts.resolve_ids.urlopen") as mock_resolve, \
         patch("scripts.fetch_fulltext.urlopen") as mock_fetch:

        # Setup mocks
        scrape_resp = MagicMock()
        scrape_resp.read.return_value = json.dumps(SAMPLE_AJAX_RESPONSE).encode()
        scrape_resp.__enter__ = lambda s: s
        scrape_resp.__exit__ = MagicMock(return_value=False)
        mock_scrape.return_value = scrape_resp

        resolve_resp = MagicMock()
        resolve_resp.read.return_value = json.dumps(SAMPLE_ID_CONVERTER_RESPONSE).encode()
        resolve_resp.__enter__ = lambda s: s
        resolve_resp.__exit__ = MagicMock(return_value=False)
        mock_resolve.return_value = resolve_resp

        fetch_resp = MagicMock()
        fetch_resp.read.return_value = json.dumps(SAMPLE_BIOC_JSON).encode()
        fetch_resp.__enter__ = lambda s: s
        fetch_resp.__exit__ = MagicMock(return_value=False)
        mock_fetch.return_value = fetch_resp

        # 1. Fetch HIRN index
        pubs = fetch_hirn_publications(cache_dir=cache_dir)
        assert len(pubs) == 2

        # 2. Search by keyword
        matches = search_publications(pubs, query="islet autoimmunity")
        assert len(matches) >= 1
        assert matches[0]["pmid"] == "39630627"

        # 3. Resolve PMCIDs
        pmids = [m["pmid"] for m in matches]
        pmcid_map = resolve_pmcids(pmids, cache_dir=cache_dir)
        assert pmcid_map["39630627"] == "PMC11615173"

        # 4. Fetch full text
        ft = fetch_fulltext("PMC11615173", cache_dir=cache_dir)
        assert ft["success"] is True
        assert len(ft["passages"]) == 6

        # 5. Chunk passages
        chunks = chunk_passages(ft["passages"])
        assert len(chunks) >= 5  # at least one chunk per non-empty passage

        # 6. Search chunks
        hits = search_chunks(chunks, query="islet autoantibodies HLA prediction")
        assert len(hits) > 0
        # Top hit should be from RESULTS or ABSTRACT (where these terms appear)
        assert hits[0]["section"] in ("RESULTS", "ABSTRACT", "DISCUSS")
        assert hits[0]["score"] > 0
```

**Step 2: Run the integration test**

Run: `python -m pytest tests/test_integration.py -v`
Expected: 1 passed

**Step 3: Run full test suite**

Run: `python -m pytest tests/ -v`
Expected: All tests pass (25+ tests)

**Step 4: Commit**

```bash
git add tests/test_integration.py
git commit -m "test: end-to-end integration test for full pipeline"
```

---

## Task 12: Final Cleanup and README

**Files:**
- Modify: `.gitignore` (if needed)
- Verify all tests pass

**Step 1: Run full test suite one last time**

Run: `python -m pytest tests/ -v --tb=short`
Expected: All tests pass

**Step 2: Verify file structure matches plan**

Run: `find . -type f | grep -v __pycache__ | grep -v '.git/' | sort`

Expected output should match the file structure from the top of this plan.

**Step 3: Final commit (if any cleanup was needed)**

```bash
git add -A
git status
# Only commit if there are changes
git commit -m "chore: final cleanup"
```

---

## Summary

| Task | Component | Tests | Lines (approx) |
|------|-----------|-------|-----------------|
| 1 | Scaffolding + fixtures | conftest.py | ~120 |
| 2 | Rate limiter | 4 tests | ~50 |
| 3 | Cache manager | 5 tests | ~70 |
| 4 | HTML parser | 9 tests | ~120 |
| 5 | HIRN scraper | 5 tests | ~90 |
| 6 | ID resolver | 4 tests | ~60 |
| 7 | Full text fetcher | 4 tests | ~80 |
| 8 | Text chunker | 6 tests | ~70 |
| 9 | Chunk searcher (BM25) | 5 tests | ~70 |
| 10 | SKILL.md | — | ~200 |
| 11 | Integration test | 1 test | ~70 |
| 12 | Final cleanup | — | — |
| **Total** | **12 tasks** | **43 tests** | **~1000** |
