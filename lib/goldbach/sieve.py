"""에라토스테네스 체 — numpy 가속 소수 생성."""

import numpy as np

_cache = {"limit": 0, "is_prime": None, "primes": None, "plimit": 0}


def sieve(limit: int) -> np.ndarray:
    """limit 이하 소수 여부 불리언 배열 반환 (인덱스 = 수)."""
    if limit <= _cache["limit"] and _cache["is_prime"] is not None:
        return _cache["is_prime"][: limit + 1]
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, int(limit**0.5) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    _cache.update(limit=limit, is_prime=is_prime, primes=None)
    return is_prime


def primes_upto(limit: int) -> np.ndarray:
    """limit 이하 모든 소수 배열."""
    # 주의: 소수 '목록'의 한계(plimit)는 체의 한계(limit)와 별개로 추적해야 한다.
    # (과거 버그: 작은 limit 호출 뒤 큰 limit 호출 시 잘린 목록 반환)
    if limit <= _cache["plimit"] and _cache["primes"] is not None:
        ps = _cache["primes"]
        return ps[ps <= limit]
    ps = np.flatnonzero(sieve(limit))
    _cache["primes"] = ps
    _cache["plimit"] = limit
    return ps


def is_prime(n: int) -> bool:
    """단일 수 소수 판정 (체 범위 밖이면 시험 나눗셈)."""
    if n < 2:
        return False
    if _cache["is_prime"] is not None and n <= _cache["limit"]:
        return bool(_cache["is_prime"][n])
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True
