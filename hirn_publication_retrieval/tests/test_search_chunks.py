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
