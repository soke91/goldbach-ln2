"""소수의 분포·증가 패턴 및 '랜덤인가?' 실험.

다루는 것:
1. 증가 패턴 — 소수 정리(PNT): pi(x) ~ x/ln x, n번째 소수 ~ n ln n
2. 간격 패턴 — 소수 간격 분포, 쌍둥이 소수, 크라메르 모델(지수 분포) 비교
3. 랜덤성 검사 — 소수는 결정론적이지만 통계적으로 '가짜 랜덤':
   - mod q 잉여류 균등분포 (디리클레) ... 랜덤처럼 보이는 증거
   - 체비쇼프 편향 (mod 4에서 3이 1보다 근소 우세) ... 랜덤이 아닌 증거
   - 연속 소수 끝자리 반발 (Lemke Oliver-Soundararajan 2016) ... 랜덤이 아닌 증거
"""

from math import log

import numpy as np

from .sieve import primes_upto


# ── 1. 증가 패턴 ──────────────────────────────────────────────

def pnt_table(xs: list[int]) -> list[dict]:
    """pi(x) 실제값 vs x/ln x vs 로그적분 근사 비교표."""
    rows = []
    ps = primes_upto(max(xs))
    for x in xs:
        actual = int(np.searchsorted(ps, x, side="right"))
        simple = x / log(x)
        li = _li(x)
        rows.append({
            "x": x, "pi(x)": actual,
            "x/ln x": round(simple, 1), "오차%": round((simple / actual - 1) * 100, 2),
            "Li(x)": round(li, 1), "Li오차%": round((li / actual - 1) * 100, 3),
        })
    return rows


def _li(x: float, steps: int = 10_000) -> float:
    """로그적분 Li(x) = ∫₂ˣ dt/ln t (사다리꼴 근사)."""
    ts = np.linspace(2, x, steps)
    return float(np.trapezoid(1 / np.log(ts), ts))


def nth_prime_growth(ns: list[int]) -> list[dict]:
    """n번째 소수 p_n vs n ln n 근사."""
    need = max(ns)
    limit = max(100, int(need * (log(need) + log(log(need)) + 1)) if need > 5 else 100)
    ps = primes_upto(limit)
    while len(ps) < need:
        limit *= 2
        ps = primes_upto(limit)
    return [
        {"n": n, "p_n": int(ps[n - 1]), "n·ln n": round(n * log(n), 1),
         "비율": round(ps[n - 1] / (n * log(n)), 4)}
        for n in ns
    ]


# ── 2. 간격 패턴 ──────────────────────────────────────────────

def gap_stats(limit: int) -> dict:
    """소수 간격 분포 — 최빈 간격, 쌍둥이 소수 개수, 최대 간격."""
    ps = primes_upto(limit)
    gaps = np.diff(ps)
    vals, cnts = np.unique(gaps, return_counts=True)
    return {
        "범위": limit,
        "소수 개수": len(ps),
        "평균 간격": round(float(gaps.mean()), 2),
        "이론 평균(ln x)": round(log(limit), 2),
        "최빈 간격": int(vals[cnts.argmax()]),
        "쌍둥이(간격2)": int(cnts[vals == 2][0]) if 2 in vals else 0,
        "최대 간격": int(gaps.max()),
        "간격 분포 상위": {int(v): int(c) for v, c in
                     sorted(zip(vals, cnts), key=lambda t: -t[1])[:8]},
    }


# ── 3. 랜덤성 실험 ────────────────────────────────────────────

def chebyshev_bias(limit: int) -> dict:
    """mod 4 경주: 소수가 4k+1 vs 4k+3 어느 쪽에 많은가.

    랜덤이라면 50:50이어야 하지만, 실제로는 4k+3이 거의 항상 근소 우세
    (체비쇼프 편향) — 소수가 '완전 랜덤이 아님'을 보여주는 고전적 현상.
    """
    ps = primes_upto(limit)
    ps = ps[ps > 2]
    r1 = int(np.sum(ps % 4 == 1))
    r3 = int(np.sum(ps % 4 == 3))
    return {"범위": limit, "4k+1": r1, "4k+3": r3,
            "차이": r3 - r1, "우세": "4k+3" if r3 > r1 else "4k+1"}


def last_digit_repulsion(limit: int, base: int = 10) -> dict:
    """연속 소수 끝자리 상관 — Lemke Oliver-Soundararajan (2016) 현상.

    소수가 랜덤이라면 연속 소수의 (끝자리, 다음 끝자리) 쌍은 균등해야
    하는데, 실제로는 같은 끝자리가 연달아 나오는 걸 '기피'한다.
    반환: 전이 행렬 {현재끝자리: {다음끝자리: 횟수}} + 같은끝자리 비율.
    """
    ps = primes_upto(limit)
    ps = ps[ps > base]  # 한 자리 소수 제외
    last = ps % base
    trans: dict[int, dict[int, int]] = {}
    for a, b in zip(last[:-1], last[1:]):
        trans.setdefault(int(a), {})[int(b)] = trans.get(int(a), {}).get(int(b), 0) + 1
    same = int(np.sum(last[:-1] == last[1:]))
    total = len(last) - 1
    digits = sorted(trans)
    expect = 1 / len(digits)  # 균등 가정 시 같은 끝자리 확률
    return {
        "범위": limit, "밑": base,
        "같은 끝자리 연속 비율": round(same / total, 4),
        "랜덤 기대치": round(expect, 4),
        "결론": "기대치보다 낮음 → 반발(랜덤 아님)" if same / total < expect else "기대치 이상",
        "전이행렬": {d: dict(sorted(trans[d].items())) for d in digits},
    }


def residue_equidistribution(limit: int, mod: int = 12) -> dict:
    """mod q 잉여류별 소수 분포 — 디리클레 정리의 수치 확인.

    gcd(r, q) = 1 인 잉여류에만 소수가 살고, 그 안에서는 거의 균등
    → '허용된 곳 안에서는 랜덤처럼' 행동한다는 증거.
    """
    ps = primes_upto(limit)
    counts = {r: int(np.sum(ps % mod == r)) for r in range(mod)}
    return {"범위": limit, "mod": mod, "분포": counts}
