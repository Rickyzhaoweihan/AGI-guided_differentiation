"""Tests for token-bucket rate limiter."""
import time
from scripts.utils.rate_limiter import RateLimiter


def test_limiter_allows_burst_up_to_capacity():
    limiter = RateLimiter(rate=10.0, capacity=3)
    # Should allow 3 requests immediately (burst capacity)
    for _ in range(3):
        wait = limiter.wait()
        assert wait < 0.05  # near-instant


def test_limiter_throttles_after_burst():
    limiter = RateLimiter(rate=10.0, capacity=1)
    limiter.wait()  # consume the 1 token
    start = time.monotonic()
    limiter.wait()  # must wait ~0.1s for refill
    elapsed = time.monotonic() - start
    assert elapsed >= 0.08  # at least ~80ms


def test_limiter_default_from_env_no_key(monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    limiter = RateLimiter.default()
    assert limiter.rate == 3.0


def test_limiter_default_from_env_with_key(monkeypatch):
    monkeypatch.setenv("NCBI_API_KEY", "fake_key")
    limiter = RateLimiter.default()
    assert limiter.rate == 10.0
