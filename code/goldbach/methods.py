"""증명 방법론 시뮬레이션 — '어떻게 증명에 접근해 왔는가'를 코드로 사고.

역사적 접근법 4가지를 수치 실험으로 재현:

1. 원 방법 (Hardy-Littlewood-Vinogradov circle method)
   지수합 S(a) = sum_p e^{2πi·a·p} 의 |S|² 구조에서 분할 수를 읽는다.
   약한 골드바흐(홀수=소수 3개 합)는 이 방법으로 2013년 증명됨(Helfgott).

2. 슈니렐만 밀도 (Schnirelmann, 1930)
   '모든 자연수는 소수 최대 C개의 합' — C를 유한하게 잡는 데 성공한
   최초의 무조건적 결과. 현재 C=4까지 내려옴 (짝수는 C=2가 골드바흐).

3. 체 이론 (Brun, Chen)
   천징룬(1973): 충분히 큰 짝수 = 소수 + (소수 또는 반소수 P2).
   '거의 골드바흐'까지 온 최강 결과를 수치로 확인.

4. 확률적 휴리스틱 (Cramér 모델)
   n이 소수일 확률 ~ 1/ln n 로 보고 g(n) 기대값을 유도 —
   기대값이 발산하므로 반례 확률의 총합이 유한(보렐-칸텔리적 사고).
"""

from math import log

import numpy as np

from .sieve import is_prime, primes_upto, sieve


# ── 1. 원 방법 ────────────────────────────────────────────────

def circle_method_demo(n: int) -> dict:
    """지수합으로 g(n) 복원: g(n) = (1/N)·sum_a |S(a)|²·e^{-2πi·a·n} 의 이산판.

    FFT로 소수 지시함수의 자기합성곱을 구하는 것과 동치 —
    '해석적 방법이 조합적 개수를 세는' 원 방법의 본질을 보여준다.
    """
    N = n + 1
    ind = np.zeros(N)
    ind[primes_upto(n)] = 1.0
    S = np.fft.rfft(ind, 2 * N)          # 지수합 S(a)
    conv = np.fft.irfft(S * S, 2 * N)    # |S|² 역변환 = 표현 수 (순서 있음)
    ordered = int(round(conv[n]))
    # 순서 무시 + p=q 케이스 보정
    half = is_prime(n // 2) and n % 2 == 0
    g = (ordered + (1 if half else 0)) // 2
    return {"n": n, "지수합 기반 g(n)": g,
            "설명": "FFT(원 방법의 이산판)로 센 분할 수 — verify 모듈과 일치해야 함"}


# ── 2. 슈니렐만 밀도 ──────────────────────────────────────────

def schnirelmann_demo(limit: int) -> dict:
    """'소수+소수' 집합의 밀도 측정 — 슈니렐만 아이디어의 축소판.

    A = {p+q : p,q 소수} 가 짝수를 얼마나 덮는지, 덮이지 않는 짝수가
    있는지 확인. 골드바흐 = 'A가 4 이상 모든 짝수를 덮는다'.
    """
    is_p = sieve(limit)
    ps = primes_upto(limit)
    covered = np.zeros(limit + 1, dtype=bool)
    for p in ps:
        qs = ps[(ps >= p) & (ps <= limit - p)]
        covered[p + qs] = True
    evens = np.arange(4, limit + 1, 2)
    missed = evens[~covered[evens]]
    return {"범위": limit,
            "짝수 커버율": f"{covered[evens].mean():.6%}",
            "미커버 짝수": missed.tolist()[:10],
            "결론": "전부 커버" if len(missed) == 0 else f"반례 후보 {len(missed)}개?!"}


# ── 3. 천징룬 방향 (소수 + P2) ────────────────────────────────

def chen_demo(limit: int) -> dict:
    """짝수 n = p + m, m은 소수 또는 반소수(P2)로 쓰는 방법의 여유도.

    골드바흐(소수+소수)보다 조건을 살짝 푼 천의 정리는 증명됐다.
    각 짝수에서 '소수+소수' 해와 '소수+P2' 해의 개수를 비교해
    조건 완화가 얼마나 큰 여유를 주는지 본다.
    """
    is_p = sieve(limit)
    # 반소수(두 소수의 곱) 표시
    ps = primes_upto(limit)
    is_p2 = np.zeros(limit + 1, dtype=bool)
    for p in ps:
        if p * p > limit:
            break
        qs = ps[(ps >= p) & (ps <= limit // p)]
        is_p2[p * qs] = True
    ok = is_p | is_p2
    sample = range(max(4, limit - 20), limit + 1, 2)
    rows = []
    for n in sample:
        gp = sum(1 for p in ps[ps <= n // 2] if is_p[n - p])
        gc = sum(1 for p in ps[ps <= n - 2] if ok[n - p])
        rows.append({"n": n, "소수+소수": gp, "소수+(소수|P2)": gc})
    return {"표본": rows, "설명": "천의 조건 완화가 주는 해의 여유 — 증명 가능/불가능의 경계 감각"}


# ── 4. 확률적 휴리스틱 ────────────────────────────────────────

def cramer_heuristic(ns: list[int]) -> list[dict]:
    """Cramér 모델: '수 m이 소수일 확률 = 1/ln m'인 랜덤 세계에서
    g(n) 기대값과 '분할이 하나도 없을 확률'을 추정.

    P(반례) ~ exp(-c·n/ln²n) → n에 대해 초지수적으로 급감.
    모든 n에 대한 합이 유한 → 랜덤 세계에서 반례는 유한 개(사실상 0개).
    골드바흐가 '확률적으로 압도적'인 이유.
    """
    rows = []
    for n in ns:
        expected = n / (2 * log(n) ** 2)           # 기대 분할 수 (대략)
        p_fail = np.exp(-expected)                  # 포아송 근사: g(n)=0 확률
        rows.append({"n": n, "기대 g(n)": round(expected, 2),
                     "P(분할 없음)": f"{p_fail:.2e}"})
    return rows


# ── 5. 약한 골드바흐 (증명된 정리) 확인 ───────────────────────

def weak_goldbach_demo(lo: int, hi: int) -> dict:
    """홀수 n(≥7) = 소수 3개 합 (Helfgott 2013 증명 완료) 수치 확인.

    강한 골드바흐가 참이면 약한 골드바흐는 자동으로 따라온다
    (n-3 이 짝수이므로). 역은 성립하지 않는다 — 두 명제의 간극이
    바로 미해결의 벽이다.
    """
    is_p = sieve(hi)
    ps = primes_upto(hi)
    fails = []
    for n in range(lo | 1, hi + 1, 2):
        if n < 7:
            continue
        found = False
        for p in ps[ps <= n - 4]:
            m = n - p  # 짝수 — 두 소수 합인지
            if any(is_p[m - q] for q in ps[ps <= m // 2]):
                found = True
                break
        if not found:
            fails.append(n)
    return {"범위": f"[{lo},{hi}] 홀수", "반례": fails,
            "결론": "전부 성립 (정리이므로 당연)" if not fails else "구현 버그!"}
