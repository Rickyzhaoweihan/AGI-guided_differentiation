"""Integration test: end-to-end pipeline with mocked HTTP."""
import json
from unittest.mock import patch, MagicMock
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


def test_end_to_end_pipeline(tmp_cache_dir):
    """Full pipeline: scrape -> search -> resolve -> fetch -> chunk -> search."""
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
