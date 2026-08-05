"""엔진 킬-테스트 — 분산 비대각 C_{k,k'}(M) = Σ_m μ(N−mk)μ(N−mk')의
실측 상쇄 (설계도 §4).

(k, k') 표본 500쌍, k,k' ~ K = N^{0.4}, M = N/K.  r = |C|/√M_eff.
제곱근-상쇄면 r ~ O(1) 반정규 — 엔진이 자연에서 작동함의 증거.
"""

import math

import numpy as np

X = 100_000_000

print("mu 계산...", flush=True)
mu = np.ones(X + 1, dtype=np.int8)
pm = np.ones(X + 1, dtype=bool)
pm[:2] = False
for p in range(2, int(X ** 0.5) + 1):
    if pm[p]:
        pm[p * p:: p] = False
        mu[p::p] *= -1
        mu[p * p:: p * p] = 0
val = np.arange(X + 1, dtype=np.int64)
for p in range(2, int(X ** 0.5) + 1):
    if pm[p]:
        val[p::p] //= p
        pp = p * p
        while pp <= X:
            val[pp::pp] //= p
            pp *= p
mu[val > 1] *= -1
mu[0] = 0
del val, pm
print("mu 완료됨", flush=True)

N = 99_999_998
K = int(X ** 0.4)          # ~ 1585
rng = np.random.default_rng(139)
rows = []
for t in range(500):
    k1 = int(rng.integers(K // 2, K))
    k2 = int(rng.integers(K // 2, K))
    if k1 == k2:
        continue
    M = (N - 1) // max(k1, k2)
    ms = np.arange(1, M + 1, dtype=np.int64)
    v1 = mu[N - k1 * ms]
    v2 = mu[N - k2 * ms]
    prod = (v1.astype(np.int16) * v2)
    C = int(prod.sum())
    M_eff = int(np.count_nonzero(prod))
    r = abs(C) / math.sqrt(max(M_eff, 1))
    rows.append((k1, k2, C, M_eff, r))
    if (t + 1) % 50 == 0:
        print(f"  {t+1}/500", flush=True)
        np.savez("dispersion_engine.npz", rows=np.array(rows))

d = np.array(rows)
r = d[:, 4]
print(f"\n[분산 엔진 비대각, 쌍 {len(d)}개 @ N=10⁸, K~{K}]")
print(f"r = |C|/√M_eff: 평균 {r.mean():.3f}  중앙값 {np.median(r):.3f}  "
      f"최악 {r.max():.3f}")
print(f"반정규 기준 0.798 — r ~ O(1)이면 엔진 작동 (제곱근 상쇄)")
big = d[r > 3]
print(f"r > 3 쌍: {len(big)}개 (반정규 기대 ~{0.0027*len(d):.1f})")
np.savez("dispersion_engine.npz", rows=d)
print("전체완료", flush=True)
