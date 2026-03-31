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
