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
