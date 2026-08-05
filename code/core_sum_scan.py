"""가장 깊은 핵 — Σ_n Λ(n)μ(N−n)의 N-분포 (q=1 이항 상관 총합).

o(N)조차 무조건부 미해결 급 (이항 Chowla 계급). 200개 N @ 10⁸에서
T(N) = Σ_p log(p)·μ(N−p) 를 √-기준 (제곱근 상쇄 규모) 대비 측정.
r(N) = |T(N)| / √(Σ log²p·μ²(N−p)).  랜덤워크면 r ~ 반정규분포 (평균 0.8).
"""

import math

import numpy as np

from goldbach.sieve import primes_upto

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

ps = primes_upto(X - 100)
ps = ps[ps > 2].astype(np.int64)
logp = np.log(ps.astype(np.float64))
logp2 = logp ** 2

rng = np.random.default_rng(127)
base = X - 10_000_000
NS = sorted(set(int(base + 6 * k + (2 - base % 6) % 6)
                for k in rng.integers(0, 1_600_000, 220)))
rows = []
for i, N in enumerate(NS):
    muv = mu[N - ps].astype(np.float64)
    T = float((logp * muv).sum())
    V = float((logp2 * (muv != 0)).sum())
    r = abs(T) / math.sqrt(V)
    rows.append((N, T, r))
    if (i + 1) % 20 == 0:
        print(f"  {i+1}/{len(NS)}", flush=True)
        np.savez("core_sum_scan.npz", rows=np.array(rows))

d = np.array(rows)
r = d[:, 2]
print(f"\n[핵 총합 T(N) = Σ Λμ, N {len(d)}개 @ ~10⁸]")
print(f"r = |T|/√V: 평균 {r.mean():.3f}  중앙값 {np.median(r):.3f}  "
      f"최악 {r.max():.3f}")
print(f"반정규 기준: 평균 0.798, P(최악>{r.max():.1f}) 참조")
print(f"T 부호: + {int((d[:,1]>0).sum())} / − {int((d[:,1]<0).sum())}")
print(f"o(N) 위반 신호(r ≫ polylog) 여부: 최악 r = {r.max():.2f} — "
      f"√-상쇄 규모면 전무")
np.savez("core_sum_scan.npz", rows=d)
print("전체완료", flush=True)
