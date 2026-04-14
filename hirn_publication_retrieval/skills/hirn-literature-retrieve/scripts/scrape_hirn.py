"""Fetch and search HIRN published works via AJAX endpoint."""
from __future__ import annotations

import json
import logging
import re
import urllib.parse
import urllib.request
from urllib.request import urlopen

from scripts.utils.cache_manager import CacheManager, TTL_INDEX
from scripts.utils.html_parser import parse_hirn_publications
from scripts.utils.rate_limiter import RateLimiter

logger = logging.getLogger(__name__)

AJAX_URL = "https://hirnetwork.org/2021/wp-admin/admin-ajax.php"
POSTS_PER_PAGE = 250
DEFAULT_CACHE_DIR = "data/cache"


def fetch_hirn_publications(
    cache_dir: str = DEFAULT_CACHE_DIR,
    consortium: str | None = None,
) -> list[dict]:
    """Fetch all HIRN publications, paginating through AJAX endpoint.

    Args:
        cache_dir: Directory for file cache.
        consortium: Optional consortium slug to filter (e.g. 'cbds').

    Returns:
        List of publication dicts with keys: pmid, doi, title,
        authors, journal, date, consortia.
    """
    cache = CacheManager(cache_dir)
    cache_key = f"hirn_index_{consortium or 'all'}"
    cached = cache.get(cache_key)
    if cached is not None:
        logger.info("Using cached HIRN publication index")
        return cached

    limiter = RateLimiter(rate=3.0, capacity=1)  # gentle on WordPress
    all_pubs: list[dict] = []
    page = 0

    while True:
        params = {
            "action": "alm_get_posts",
            "post_type": "publication",
            "posts_per_page": str(POSTS_PER_PAGE),
            "page": str(page),
        }
        if consortium:
            params["taxonomy"] = "publication_consortium"
            params["taxonomy_terms"] = consortium.lower()
        else:
            params["taxonomy"] = "publication_category"
            params["taxonomy_terms"] = "published_works"

        url = f"{AJAX_URL}?{urllib.parse.urlencode(params)}"
        limiter.wait()

        logger.info("Fetching HIRN publications page %d", page)
        with urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))

        html = data.get("html", "")
        meta = data.get("meta", {})
        pubs = parse_hirn_publications(html)
        all_pubs.extend(pubs)

        total = meta.get("totalposts", 0)
        if not pubs or len(all_pubs) >= total:
            break
        page += 1

    logger.info("Fetched %d HIRN publications", len(all_pubs))
    cache.set(cache_key, all_pubs, ttl=TTL_INDEX)
    return all_pubs


def search_publications(
    publications: list[dict],
    query: str,
    consortium: str | None = None,
    max_results: int = 20,
) -> list[dict]:
    """Search publications by keyword matching on titles.

    Args:
        publications: List of publication dicts.
        query: Search query string.
        consortium: Optional consortium filter.
        max_results: Maximum results to return.

    Returns:
        List of matching publications, sorted by relevance (best first).
    """
    if consortium:
        publications = [
            p for p in publications
            if consortium.upper() in [c.upper() for c in p.get("consortia", [])]
        ]

    query_terms = re.findall(r"\w+", query.lower())
    if not query_terms:
        return []

    scored: list[tuple[float, dict]] = []
    for pub in publications:
        title_lower = pub.get("title", "").lower()
        # Score: fraction of query terms found in title
        matches = sum(1 for t in query_terms if t in title_lower)
        if matches > 0:
            score = matches / len(query_terms)
            scored.append((score, pub))

    scored.sort(key=lambda x: x[0], reverse=True)
    return [pub for _, pub in scored[:max_results]]
