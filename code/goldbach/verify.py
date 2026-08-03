"""골드바흐 분할 검증 엔진."""

from dataclasses import dataclass

import numpy as np

from .sieve import primes_upto, sieve


@dataclass
class Partition:
    n: int
    p: int
    q: int  # n = p + q, p <= q, 둘 다 소수


def find_partition(n: int) -> Partition | None:
    """짝수 n의 골드바흐 분할 하나 찾기 (가장 작은 p)."""
    if n < 4 or n % 2:
        return None
    is_p = sieve(n)
    for p in primes_upto(n // 2):
        if is_p[n - p]:
            return Partition(n, int(p), int(n - p))
    return None  # 발견되면 골드바흐 추측의 반례!


def all_partitions(n: int) -> list[Partition]:
    """짝수 n의 모든 골드바흐 분할."""
    if n < 4 or n % 2:
        return []
    is_p = sieve(n)
    return [
        Partition(n, int(p), int(n - p))
        for p in primes_upto(n // 2)
        if is_p[n - p]
    ]


def partition_counts(limit: int) -> np.ndarray:
    """4..limit 모든 짝수의 분할 개수 g(n)을 한꺼번에 계산.

    반환: counts[n] = g(n) (홀수 인덱스는 0). 골드바흐 혜성의 원료.
    """
    is_p = sieve(limit)
    ps = primes_upto(limit)
    counts = np.zeros(limit + 1, dtype=np.int64)
    for p in ps:
        # p + q = n, q >= p 인 소수 q → n = p+q 위치에 +1
        qs = ps[(ps >= p) & (ps <= limit - p)]
        np.add.at(counts, p + qs, 1)
    counts[1::2] = 0
    return counts


def verify_range(lo: int, hi: int) -> tuple[int, list[int]]:
    """[lo, hi] 짝수 전수 검증. (검증 개수, 반례 목록) 반환."""
    lo = max(4, lo + (lo % 2))
    is_p = sieve(hi)
    ps = primes_upto(hi)
    checked, failures = 0, []
    # 짝수 n마다 최소 소수 p를 찾음 — 작은 소수에서 대부분 즉시 해결됨
    small = ps[: min(len(ps), 100)]
    for n in range(lo, hi + 1, 2):
        checked += 1
        if any(is_p[n - p] for p in small if p <= n - 2):
            continue
        if not any(is_p[n - p] for p in ps[ps <= n // 2]):
            failures.append(n)
    return checked, failures
