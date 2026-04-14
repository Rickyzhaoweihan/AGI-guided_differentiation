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
