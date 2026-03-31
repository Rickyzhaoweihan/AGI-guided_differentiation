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
