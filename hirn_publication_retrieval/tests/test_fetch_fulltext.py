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
