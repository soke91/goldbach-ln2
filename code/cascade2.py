"""E2 — 캐스케이드 2단: k-대역 전수 dual(k)² 집계 vs 독립 예측.

dual(k) = Σ_{√N<m<N/k} μ(m)μ(N−mk).  대역 [K, 2K] 전수:
  Q = Σ_k dual(k)²  vs  Q_ind = Σ_k V_k (V_k = 비영 항수 = 독립 기대)
비 Q/Q_ind ≈ 1 이면 듀얼 족이 독립-급 — 4점 평균 객체 자체 상쇄
(E1의 자연-측 근거). ≫1 이면 족 내 상관 = 캐스케이드 장애.
세 대역 (K = 2000, 8000, 20000).
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
SQ = int(N ** 0.5)

for K in [2000, 8000, 20000]:
    Q = 0.0
    Q_ind = 0.0
    duals = []
    for k in range(K, 2 * K):
        M0 = SQ + 1
        M1 = (N - 1) // k
        if M1 <= M0:
            continue
        ms = np.arange(M0, M1 + 1, dtype=np.int64)
        vals = mu[ms].astype(np.int16) * mu[N - k * ms]
        d = int(vals.sum(dtype=np.int64))
        v = int(np.count_nonzero(vals))
        Q += d * d
        Q_ind += v
        duals.append(d / math.sqrt(max(v, 1)))
    duals = np.array(duals)
    print(f"K=[{K},{2*K})  전수 {len(duals):,}  Q/Q_ind = {Q/Q_ind:.3f}  "
          f"(독립=1)  r분포: 평균|r| {np.abs(duals).mean():.3f} "
          f"최악 {np.abs(duals).max():.2f}", flush=True)
    np.savez(f"cascade2_{K}.npz", duals=duals, Q=Q, Q_ind=Q_ind)
print("전체완료", flush=True)
