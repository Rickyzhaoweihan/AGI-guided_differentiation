"""Resolve PMIDs/DOIs to PMCIDs via NCBI ID Converter API."""
from __future__ import annotations

import json
import logging
import urllib.parse
from urllib.request import urlopen

from scripts.utils.cache_manager import CacheManager, TTL_IDS
from scripts.utils.rate_limiter import RateLimiter

logger = logging.getLogger(__name__)

ID_CONVERTER_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
BATCH_SIZE = 200
DEFAULT_CACHE_DIR = "data/cache"


def resolve_pmcids(
    pmids: list[str],
    cache_dir: str = DEFAULT_CACHE_DIR,
) -> dict[str, str | None]:
    """Resolve a list of PMIDs to PMCIDs.

    Args:
        pmids: List of PMID strings.
        cache_dir: Directory for file cache.

    Returns:
        Dict mapping PMID -> PMCID (or None if not in PMC).
    """
    cache = CacheManager(cache_dir)
    limiter = RateLimiter.default()
    result: dict[str, str | None] = {}
    uncached: list[str] = []

    # Check cache first
    for pmid in pmids:
        cached = cache.get(f"pmcid_{pmid}")
        if cached is not None:
            result[pmid] = cached.get("pmcid")
        else:
            uncached.append(pmid)

    if not uncached:
        return result

    # Batch resolve uncached PMIDs
    for i in range(0, len(uncached), BATCH_SIZE):
        batch = uncached[i : i + BATCH_SIZE]
        params = {
            "ids": ",".join(batch),
            "format": "json",
            "tool": "hirn-literature-retrieve",
            "email": "tool@example.com",
        }
        url = f"{ID_CONVERTER_URL}?{urllib.parse.urlencode(params)}"
        limiter.wait()

        logger.info("Resolving PMCIDs for %d PMIDs", len(batch))
        with urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))

        for record in data.get("records", []):
            pmid = str(record.get("pmid", record.get("requested-id", "")))
            is_error = record.get("status") == "error"
            pmcid = record.get("pmcid") if not is_error else None
            result[pmid] = pmcid
            cache.set(f"pmcid_{pmid}", {"pmcid": pmcid}, ttl=TTL_IDS)

    return result
