"""Token-bucket rate limiter for NCBI API calls."""
from __future__ import annotations

import os
import threading
import time


class RateLimiter:
    """Thread-safe token-bucket rate limiter.

    Args:
        rate: Maximum sustained requests per second.
        capacity: Burst capacity (tokens available immediately).
    """

    def __init__(self, rate: float = 3.0, capacity: int = 1) -> None:
        self.rate = rate
        self.capacity = capacity
        self._tokens = float(capacity)
        self._last = time.monotonic()
        self._lock = threading.Lock()

    @classmethod
    def default(cls) -> RateLimiter:
        """Create a limiter configured from NCBI_API_KEY env var."""
        has_key = bool(os.environ.get("NCBI_API_KEY"))
        rate = 10.0 if has_key else 3.0
        return cls(rate=rate, capacity=1)

    def wait(self) -> float:
        """Block until a token is available. Returns seconds waited."""
        with self._lock:
            now = time.monotonic()
            elapsed = now - self._last
            self._tokens = min(self.capacity, self._tokens + elapsed * self.rate)
            self._last = now

            if self._tokens >= 1.0:
                self._tokens -= 1.0
                return 0.0

            deficit = 1.0 - self._tokens
            sleep_time = deficit / self.rate
            self._tokens = 0.0

        time.sleep(sleep_time)

        with self._lock:
            self._last = time.monotonic()

        return sleep_time
