"""A5-3 창-계통 판정 — 10^7 고정, 창 폭만 변화 (W/X = 0.06 vs 0.3).

기존 10^7 측정(R2=0.708)은 W=600k (W/X=0.06). W=3M (0.3)으로 재측정.
R2가 W에 끌려가면 6점 곡선의 U자는 창 계통, 불변이면 스케일 실재.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 10_000_000
W = 3_000_000
M = 3000
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)
y = int(round(X ** (1 / 3))) + 1
qs = [int(q) for q in primes_upto(y) if q > 2]

cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
rng = np.random.default_rng(17)
cands = sorted(rng.choice(cands, M, replace=False).tolist())

s_vals = np.empty(M)
N_vals = np.empty(M)
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
    if (i + 1) % 300 == 0:
        print(f"  {i+1}/{M}", flush=True)

s = s_vals
sig_obs = s.std()
sig_bin = float(np.mean(np.sqrt(s * (1 - s) / N_vals)))
R2 = (sig_obs / sig_bin) ** 2
se = R2 * np.sqrt(2 / (M - 1))
print(f"\nX={X:,}  W/X=0.3  표본 {M}  s {s.mean():.5f}")
print(f"R2 = {R2:.4f} +- {se:.4f}   (기존 W/X=0.06: 0.7076)")
np.savez("s_binom_window.npz", cands=np.array(cands), s=s_vals, N=N_vals)
