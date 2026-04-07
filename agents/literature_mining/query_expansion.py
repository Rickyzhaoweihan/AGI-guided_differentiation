#!/usr/bin/env python3
"""
query_expansion.py — Generate expanded PubMed queries for scRNA-seq dataset retrieval.

Given a user query (e.g. "human gut scRNA-seq"), produces multiple sub-queries
covering technology variants, tissue synonyms, organoid/iPSC, atlas, and
developmental categories. Disease-specific queries are NOT generated (user can
include disease terms in their base query if desired).
"""

from __future__ import annotations

import re


# Map common tissue keywords to synonyms for query expansion
TISSUE_SYNONYMS: dict[str, list[str]] = {
    "gut":       ["intestine", "intestinal", "colon", "duodenum", "ileum", "enterocyte", "gastrointestinal"],
    "intestin":  ["gut", "colon", "duodenum", "ileum", "enterocyte", "gastrointestinal", "bowel"],
    "colon":     ["intestine", "intestinal", "gut", "bowel"],
    "pancrea":   ["islet", "islets", "beta cell", "endocrine pancreas"],
    "islet":     ["pancreas", "pancreatic", "beta cell", "endocrine pancreas"],
    "heart":     ["cardiac", "cardiomyocyte", "myocardium"],
    "liver":     ["hepatocyte", "hepatic"],
    "lung":      ["pulmonary", "airway", "alveolar"],
    "kidney":    ["renal", "nephron"],
    "brain":     ["cerebral", "cortex", "neuron", "neural"],
}

# Regex to strip technology terms from the query to get the "base" tissue terms
_TECH_RE = re.compile(
    r"(?i)\b(scrna-?seq|snrna-?seq|single[- ]cell|single[- ]nucleus|"
    r"RNA[- ]seq(?:uencing)?|transcriptom\w*|scRNA|10x|chromium)\b"
)


def expand_query(query: str) -> list[str]:
    """
    Expand a user query into multiple PubMed sub-queries to maximize recall.

    Categories:
      1. Core query (as-is)
      2. Technology variants (scRNA-seq, snRNA-seq, single-cell transcriptome)
      3. Tissue synonyms (capped at top 3)
      4. Organoid / iPSC / stem cell differentiation
      5. Atlas / multi-organ
      6. Developmental (fetal, embryo, organogenesis)

    Returns deduplicated list of query strings.
    """
    queries: list[str] = [query]
    q_lower = query.lower()

    # Extract base tissue terms (remove technology keywords)
    base = _TECH_RE.sub("", query).strip()
    base = re.sub(r"\s+", " ", base).strip()

    # --- 1. Technology variants ---
    if base and len(base) > 3:
        queries.append(f"{base} single-cell RNA-seq")
        queries.append(f"{base} single-nucleus RNA-seq")
        queries.append(f"{base} single-cell transcriptome")

    # --- 2. Tissue synonyms (top 3) ---
    for key, syns in TISSUE_SYNONYMS.items():
        if key in q_lower:
            for syn in syns[:3]:
                if syn.lower() not in q_lower:
                    queries.append(f"human {syn} single-cell RNA-seq")
            break  # only match first tissue

    # --- 3. Organoid / iPSC ---
    if base and len(base) > 3:
        queries.append(f"{base} organoid single-cell RNA-seq")
        queries.append(f"human pluripotent stem cell {base} differentiation single-cell")
        queries.append(f"iPSC {base} single-cell RNA-seq")

    # --- 4. Atlas / multi-organ ---
    queries.append("human cell atlas single-cell transcriptome")
    queries.append("Tabula Sapiens single-cell")
    queries.append("human fetal gene expression atlas single-cell")
    if base and len(base) > 3:
        queries.append(f"{base} cell atlas single-cell")

    # --- 5. Developmental ---
    if base and len(base) > 3:
        queries.append(f"{base} fetal development single-cell RNA-seq")
        queries.append(f"human embryo organogenesis single-cell {base}")
        queries.append(f"{base} embryo single-cell transcriptome")

    # --- Deduplicate preserving order ---
    seen: set[str] = set()
    unique: list[str] = []
    for q in queries:
        q_norm = re.sub(r"\s+", " ", q.strip().lower())
        if q_norm not in seen:
            seen.add(q_norm)
            unique.append(q.strip())

    return unique
