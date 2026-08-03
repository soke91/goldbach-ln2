"""S3-2 — 장벽 지도 저-theta 연장 (gamma -> 2 바닥 확인).

theta in [0.125, 0.24] (f(s) 유효역 2<=s=0.5/theta<=4 전체).
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

EG = 0.5772156649015329
X = 2_000_000
W = 100_000
M = 800
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)

cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
rng = np.random.default_rng(3)
cands = sorted(rng.choice(cands, M, replace=False).tolist())

print(f"X={X:,}  표본 {M}")
print(f"{'theta':>6} {'s_sieve':>7} {'s평균':>8} {'price':>7} {'L=f(s)':>8} "
      f"{'U/L':>7} {'gamma':>8}")

for th in [0.125, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24]:
    y = X ** th
    qs = [int(q) for q in primes_upto(int(y) + 1) if q > 2]
    s_vals = np.empty(M)
    for i, n in enumerate(cands):
        m = n - ps_all[ps_all <= n // 2]
        al = np.ones(len(m), dtype=bool)
        for q in qs:
            if n % q:
                al &= (m % q != 0)
        al &= m > 1
        surv = m[al]
        s_vals[i] = float(is_p[surv].mean())
    smean = s_vals.mean()
    price = 1 / (1 - smean)
    s_sieve = 0.5 / th
    L = 2 * np.exp(EG) * np.log(s_sieve - 1) / s_sieve if s_sieve > 2 else 0.0
    UL = 2 / L if L > 0 else float("inf")
    gam = UL / price
    print(f"{th:>6.3f} {s_sieve:>7.2f} {smean:>8.5f} {price:>7.4f} "
          f"{L:>8.4f} {UL:>7.3f} {gam:>8.3f}", flush=True)
