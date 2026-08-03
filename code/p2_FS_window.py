"""F_S의 GM-창 가설 결정 시험 — 창 폭에 따른 F_S 이동.

가설: F_S ≈ 1 − ln(λ_eff)/ln x, λ_eff ∝ 창 폭 (GM 구간-분산비)
예측 (x = 10⁷): 창 1.15× → ~0.18 / 1.3× → 0.137(기준) / 1.6× → ~0.09
창-불변이면 가설 기각, GM-이동이면 F_S 유도 완성.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 10_000_000
NS = 900
print(f"{'창 비율':>7} {'F_S':>7} {'GM 예측':>8}")
for RATIO in [1.15, 1.3, 1.6, 2.2]:
    HI = int(X * RATIO)
    is_p = sieve(HI + 100_000)
    ps_all = primes_upto(HI // 2 + 50_000)
    y = int(round(X ** (1 / 3))) + 1
    qs = [int(q) for q in primes_upto(y) if q > 2]
    base = np.linspace(X, HI - 10, NS).astype(np.int64)
    cands = [int(v + (2 - v % 6) % 6) for v in base]
    S = []
    for n in cands:
        m = n - ps_all[ps_all <= n // 2]
        al = np.ones(len(m), dtype=bool)
        for q in qs:
            if n % q:
                al &= (m % q != 0)
        al &= m > 1
        S.append(int(al.sum()))
    S = np.array(S, float)
    W = np.ones(len(cands))
    for q in qs:
        mask = (np.array(cands) % q) != 0
        W[mask] *= 1 - 1 / (q - 1)
    x1 = np.linspace(0, 1, len(S))
    z = S / W
    eS = (z - np.polyval(np.polyfit(x1, z, 2), x1)) * W
    FS = float(eS.var() / S.mean())
    lam_eff = (RATIO - 1) * X / 3  # 2차 디트렌드의 유효 파장 ~ 창/3 (보정 상수)
    gm = max(1 - np.log(lam_eff) / np.log(X), 0.02)
    print(f"{RATIO:>7.2f} {FS:>7.4f} {gm:>8.4f}", flush=True)
