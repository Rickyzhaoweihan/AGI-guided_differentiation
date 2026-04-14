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
