"""골드바흐 분할 통계 — g(n) 분포와 Hardy-Littlewood 추정 비교.

핵심 사고: 골드바흐가 '왜 참일 수밖에 없어 보이는가'의 정량적 근거.
g(n)(분할 개수)은 n이 커질수록 오히려 증가한다 — 반례가 나오려면
g(n) = 0 이어야 하는데, 예측치는 ~ n / ln²n 으로 발산한다.
"""

from math import log

import numpy as np

from .sieve import primes_upto
from .verify import partition_counts

# 쌍둥이 소수 상수 C2 = prod (1 - 1/(p-1)^2) over odd primes
TWIN_PRIME_C2 = 0.6601618158468696


def hardy_littlewood_estimate(n: int) -> float:
    """확장 골드바흐 추측(Hardy-Littlewood)의 g(n) 점근 추정.

    g(n) ~ 2·C2 · [prod_{p|n, p>2} (p-1)/(p-2)] · n / ln²n
    (순서 무시 분할 기준으로 1/2 배)
    """
    if n < 4 or n % 2:
        return 0.0
    corr = 1.0
    m, p = n, 3
    while m % 2 == 0:
        m //= 2
    while p * p <= m:
        if m % p == 0:
            corr *= (p - 1) / (p - 2)
            while m % p == 0:
                m //= p
        p += 2
    if m > 1:  # 남은 소인수
        corr *= (m - 1) / (m - 2)
    return TWIN_PRIME_C2 * corr * n / log(n) ** 2


def comet_data(limit: int) -> tuple[np.ndarray, np.ndarray]:
    """골드바흐 혜성 원료: (짝수 n 배열, g(n) 배열)."""
    counts = partition_counts(limit)
    ns = np.arange(4, limit + 1, 2)
    return ns, counts[ns]


def minimum_g_growth(limit: int, buckets: int = 20) -> list[dict]:
    """구간별 g(n)의 최솟값 추이 — '0에 가까워지는가?'에 대한 답.

    반례가 존재하려면 어딘가에서 min g(n) = 0 이어야 한다.
    실제로는 구간 최솟값조차 증가 추세 → 추측이 참일 강력한 정황.
    """
    ns, gs = comet_data(limit)
    edges = np.linspace(4, limit, buckets + 1, dtype=int)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (ns >= lo) & (ns < hi)
        if not mask.any():
            continue
        i = np.argmin(gs[mask])
        rows.append({
            "구간": f"[{lo}, {hi})",
            "min g(n)": int(gs[mask][i]),
            "달성 n": int(ns[mask][i]),
            "HL 예측(구간초)": round(hardy_littlewood_estimate(int(lo) + int(lo) % 2), 1),
        })
    return rows


def hl_accuracy(sample_ns: list[int], limit: int | None = None) -> list[dict]:
    """표본 n들에서 실제 g(n) vs Hardy-Littlewood 예측 비교."""
    from .verify import all_partitions
    rows = []
    for n in sample_ns:
        actual = len(all_partitions(n))
        pred = hardy_littlewood_estimate(n)
        rows.append({"n": n, "실제 g(n)": actual, "HL 예측": round(pred, 1),
                     "비율": round(actual / pred, 3) if pred else None})
    return rows
