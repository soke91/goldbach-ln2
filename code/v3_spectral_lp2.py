"""V3 실험 1 — 스펙트럼-제약 소수 덮개의 내재 비용 (LP).

min Σ_m w(m)  s.t.  w(m) ≥ 1 (m 소수), w(m) ≥ 0 (전 구간),
w의 스펙트럼 ⊂ {a/q: q ≤ Q} major-arc 주파수 집합 (+상수항).

w(m) = c₀ + Σ_j [A_j cos(2π f_j m) + B_j sin(2π f_j m)] — LP 변수 (c,A,B).
목적 Σ_m w(m) = N·c₀ + (진동항 합, 거의 0) ≈ N·c₀ → c₀ 최소화.

판정량: cost = (min Σw)/π(N) — 소수 1개당 덮개 비용. Λ² 실현 손실
(≈2.1)과 비교. Q-스케일링으로 "스펙트럼 자원 ↑ → 비용 ↓" 속도 측정.
"""

import numpy as np
from scipy.optimize import linprog

from goldbach.sieve import sieve

N = 60_000
is_p = sieve(N + 10)
ms = np.arange(3, N, 2)
pmask = is_p[ms]
print(f"N={N:,}  홀수 {len(ms):,}  소수 {pmask.sum():,}", flush=True)

for Q in [6, 12, 20, 30]:
    # major-arc 주파수 a/q, q ≤ Q, (a,q)=1, 0 < a/q < 1/2 (cos/sin 쌍)
    freqs = []
    from math import gcd
    for q in range(2, Q + 1):
        for a in range(1, q // 2 + 1):
            if gcd(a, q) == 1 and 0 < a / q <= 0.5:
                freqs.append(a / q)
    freqs = sorted(set(freqs))
    F = len(freqs)

    # 설계 행렬: 열 = [c0, A_1..A_F, B_1..B_F]
    ang = 2 * np.pi * np.outer(ms, freqs)
    X = np.hstack([np.ones((len(ms), 1)), np.cos(ang), np.sin(ang)])
    nvar = 1 + 2 * F

    # 목적: 총합 최소화 = X 열합 · 변수
    # 제약: 소수에서 w ≥ 1  (−X_p v ≤ −1) / 전체 w ≥ 0 (−X v ≤ 0)
    A_ub = -X
    b_ub = np.where(pmask, -1.0, 0.0)
    bounds = [(0.0, 10.0)] + [(-10.0, 10.0)] * (2 * F)
    keep = np.abs(X).sum(axis=0) > 1e-6
    X = X[:, keep]
    c_obj = X.sum(axis=0)
    A_ub = -X
    bounds = [b for b, k in zip(bounds, keep) if k]
    nvar = int(keep.sum())
    r = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=bounds,
                method="highs", options={"presolve": True})
    if not r.success:
        print(f"Q={Q}: LP 실패 {r.message}", flush=True)
        continue
    total = float(c_obj @ r.x)
    cost = total / pmask.sum()
    print(f"Q={Q:>3}  주파수 {F:>4}개  덮개비용/소수 = {cost:.3f}  "
          f"(Λ² 실현 벤치마크 ~2.1)", flush=True)
