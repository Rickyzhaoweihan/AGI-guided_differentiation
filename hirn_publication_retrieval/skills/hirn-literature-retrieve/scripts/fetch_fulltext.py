"""Fetch full text from PMC Open Access via BioC JSON API."""
from __future__ import annotations

import json
import logging
import urllib.error
from urllib.request import urlopen

from scripts.utils.cache_manager import CacheManager, TTL_FULLTEXT
from scripts.utils.rate_limiter import RateLimiter

logger = logging.getLogger(__name__)

BIOC_URL = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{pmcid}/unicode"
DEFAULT_CACHE_DIR = "data/cache"


def parse_bioc_passages(bioc_data: dict | list) -> list[dict]:
    """Extract passages from BioC JSON response.

    Args:
        bioc_data: Parsed BioC JSON — either a dict with "documents" key
            or a list containing such a dict (both formats are used by the API).

    Returns:
        List of passage dicts with keys: section, type, offset, text.
    """
    # API sometimes returns a list wrapping the document object
    if isinstance(bioc_data, list):
        bioc_data = bioc_data[0] if bioc_data else {}
    passages = []
    for doc in bioc_data.get("documents", []):
        for passage in doc.get("passages", []):
            infons = passage.get("infons", {})
            passages.append({
                "section": infons.get("section_type", "UNKNOWN"),
                "type": infons.get("type", "unknown"),
                "offset": passage.get("offset", 0),
                "text": passage.get("text", ""),
            })
    return passages


def fetch_fulltext(
    pmcid: str,
    cache_dir: str = DEFAULT_CACHE_DIR,
) -> dict:
    """Fetch full text for a PMC article.

    Args:
        pmcid: PMC ID (e.g. 'PMC11615173').
        cache_dir: Directory for file cache.

    Returns:
        Dict with keys: success, pmcid, passages, error.
    """
    cache = CacheManager(cache_dir)
    cache_key = f"fulltext_{pmcid}"
    cached = cache.get(cache_key)
    if cached is not None:
        logger.info("Cache hit for full text: %s", pmcid)
        return cached

    limiter = RateLimiter.default()
    url = BIOC_URL.format(pmcid=pmcid)
    limiter.wait()

    try:
        logger.info("Fetching full text for %s", pmcid)
        with urlopen(url, timeout=60) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return {
                "success": False,
                "pmcid": pmcid,
                "passages": [],
                "error": {
                    "code": "NOT_FOUND",
                    "message": f"Full text not available in PMC Open Access for {pmcid}",
                },
            }
        return {
            "success": False,
            "pmcid": pmcid,
            "passages": [],
            "error": {"code": "API_ERROR", "message": str(e)},
        }
    except urllib.error.URLError as e:
        return {
            "success": False,
            "pmcid": pmcid,
            "passages": [],
            "error": {"code": "NETWORK_ERROR", "message": str(e)},
        }
    except json.JSONDecodeError as e:
        return {
            "success": False,
            "pmcid": pmcid,
            "passages": [],
            "error": {"code": "PARSE_ERROR", "message": f"Invalid JSON from BioC API for {pmcid}: {e}"},
        }

    passages = parse_bioc_passages(data)
    result = {"success": True, "pmcid": pmcid, "passages": passages}
    cache.set(cache_key, result, ttl=TTL_FULLTEXT)
    return result
