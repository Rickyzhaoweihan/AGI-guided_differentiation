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
