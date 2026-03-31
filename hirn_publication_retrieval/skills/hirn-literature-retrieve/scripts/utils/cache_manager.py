"""File-based cache with TTL expiration."""
from __future__ import annotations

import hashlib
import json
import logging
import threading
import time
from pathlib import Path

logger = logging.getLogger(__name__)

# Default TTLs in seconds
TTL_INDEX = 86400       # 1 day  — HIRN publication index
TTL_METADATA = 2592000  # 30 days — article metadata
TTL_FULLTEXT = 2592000  # 30 days — article full text
TTL_IDS = 2592000       # 30 days — ID mappings


class CacheManager:
    """Thread-safe file-based JSON cache with TTL expiration.

    Args:
        cache_dir: Directory to store cache files.
    """

    def __init__(self, cache_dir: str | Path) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def _key_path(self, key: str) -> Path:
        hashed = hashlib.md5(key.encode()).hexdigest()
        return self.cache_dir / f"{hashed}.json"

    def get(self, key: str) -> dict | list | None:
        """Retrieve a cached value. Returns None on miss or expiry."""
        path = self._key_path(key)
        if not path.exists():
            return None
        try:
            with self._lock:
                data = json.loads(path.read_text())
            if time.time() > data["expires_at"]:
                path.unlink(missing_ok=True)
                return None
            return data["value"]
        except (json.JSONDecodeError, KeyError):
            path.unlink(missing_ok=True)
            return None

    def set(self, key: str, value: dict | list, ttl: int = 3600) -> None:
        """Store a value with a TTL in seconds."""
        path = self._key_path(key)
        envelope = {
            "key": key,
            "value": value,
            "created_at": time.time(),
            "expires_at": time.time() + ttl,
        }
        with self._lock:
            path.write_text(json.dumps(envelope))

    def cleanup(self) -> int:
        """Remove all expired entries. Returns count removed."""
        removed = 0
        for path in self.cache_dir.glob("*.json"):
            try:
                data = json.loads(path.read_text())
                if time.time() > data["expires_at"]:
                    path.unlink()
                    removed += 1
            except (json.JSONDecodeError, KeyError):
                path.unlink()
                removed += 1
        return removed
