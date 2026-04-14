"""Tests for file-based cache manager."""
import json
import time
from scripts.utils.cache_manager import CacheManager


def test_cache_set_and_get(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    cache.set("test_key", {"data": "hello"}, ttl=60)
    result = cache.get("test_key")
    assert result == {"data": "hello"}


def test_cache_miss_returns_none(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    assert cache.get("nonexistent") is None


def test_cache_expired_returns_none(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    cache.set("expire_me", {"data": 1}, ttl=0)
    time.sleep(0.05)
    assert cache.get("expire_me") is None


def test_cache_cleanup_removes_expired(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    cache.set("old", {"data": 1}, ttl=0)
    cache.set("fresh", {"data": 2}, ttl=3600)
    time.sleep(0.05)
    removed = cache.cleanup()
    assert removed >= 1
    assert cache.get("old") is None
    assert cache.get("fresh") == {"data": 2}


def test_cache_key_hashing_is_deterministic(tmp_cache_dir):
    cache = CacheManager(tmp_cache_dir)
    cache.set("same_key", {"v": 1}, ttl=60)
    cache.set("same_key", {"v": 2}, ttl=60)
    assert cache.get("same_key") == {"v": 2}
