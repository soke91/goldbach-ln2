"""A5-3 / 킬-테스트 — "s-요동 = 생존자 집합의 순수 이항 샘플링 잡음" 가설.

s(n) = (survivors 중 소수 비율). n마다 생존자 수 N(n)이 있고, 만약 s(n)의
n-간 산포가 sqrt(s(1-s)/N) (이항 표준편차)와 일치하면 s 채널은 구조 없는
샘플링 잡음이다. 초과분(비 1.0 초과)은 미지의 구조.

여러 스케일에서 비 R = σ_obs / σ_bin 측정. R→1이면 가설 생존.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

for X, W, M in [(100_000, 60_000, 5000), (1_000_000, 180_000, 5000),
                (10_000_000, 600_000, 3000)]:
    is_p = sieve(X + W + 10)
    ps_all = primes_upto((X + W) // 2 + 10)
    y = int(round(X ** (1 / 3))) + 1
    qs = [int(q) for q in primes_upto(y) if q > 2]

    cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
    rng = np.random.default_rng(7)
    if len(cands) > M:
        cands = sorted(rng.choice(cands, M, replace=False).tolist())

    s_vals = np.empty(len(cands))
    N_vals = np.empty(len(cands))
    for i, n in enumerate(cands):
        m = n - ps_all[ps_all <= n // 2]
        al = np.ones(len(m), dtype=bool)
        for q in qs:
            if n % q:
                al &= (m % q != 0)
        al &= m > 1
        surv = m[al]
        N_vals[i] = len(surv)
        s_vals[i] = float(is_p[surv].mean())

    s = s_vals
    mu = s.mean()
    sig_obs = s.std()
    sig_bin = float(np.mean(np.sqrt(s * (1 - s) / N_vals)))
    R = sig_obs / sig_bin
    print(f"X={X:>12,}  표본 {len(s):>5}  N_surv 평균 {N_vals.mean():>9.0f}  "
          f"s {mu:.5f}  sig_obs {sig_obs:.5f}  sig_bin {sig_bin:.5f}  "
          f"R = {R:.4f}", flush=True)
