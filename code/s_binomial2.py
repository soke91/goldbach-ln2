"""A5-3 확장 — 분산비 R² = var_obs/var_bin 의 스케일 적합 (ln2 극한 검증).

가설(신규): R²(X) → ln 2 = 0.69315. 스케일 4e6, 4e7, 1e8 추가 측정 후
기존 3점과 합쳐 1/lnX 및 1/ln²X 모형 적합.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

results = []
for X, W, M in [(4_000_000, 400_000, 3000), (40_000_000, 1_200_000, 1200),
                (100_000_000, 2_000_000, 600)]:
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
        if (i + 1) % 300 == 0:
            print(f"  X={X:,}: {i+1}/{len(cands)}", flush=True)

    s = s_vals
    sig_obs = s.std()
    sig_bin = float(np.mean(np.sqrt(s * (1 - s) / N_vals)))
    R2 = (sig_obs / sig_bin) ** 2
    results.append((X, R2))
    print(f"X={X:>13,}  표본 {len(s):>5}  s {s.mean():.5f}  "
          f"R2 = {R2:.4f}", flush=True)
    np.savez(f"s_binom_{X}.npz", cands=np.array(cands), s=s_vals, N=N_vals)

# 기존 3점 합류
prev = [(100_000, 0.8199), (1_000_000, 0.7276), (10_000_000, 0.7076)]
allpts = sorted(prev + results)
lx = np.array([np.log(x) for x, _ in allpts])
r2 = np.array([v for _, v in allpts])
print("\n[적합] R2 = a + b/lnX + c/lnX^2")
A = np.vstack([np.ones_like(lx), 1 / lx, 1 / lx ** 2]).T
coef, res_, *_ = np.linalg.lstsq(A, r2, rcond=None)
pred = A @ coef
print(f"  a(극한) = {coef[0]:.4f}  b = {coef[1]:+.3f}  c = {coef[2]:+.2f}")
print(f"  잔차 rms = {np.sqrt(np.mean((r2-pred)**2)):.5f}")
print(f"  ln2 = {np.log(2):.4f}  | 편차 = {coef[0] - np.log(2):+.4f}")
A1 = np.vstack([np.ones_like(lx), 1 / lx]).T
c1, *_ = np.linalg.lstsq(A1, r2, rcond=None)
print(f"[1차 모형] a = {c1[0]:.4f}  b = {c1[1]:+.3f}  "
      f"잔차 rms = {np.sqrt(np.mean((r2 - A1@c1)**2)):.5f}")
for (x, v), p in zip(allpts, pred):
    print(f"  X=10^{np.log10(x):.1f}  R2 {v:.4f}  적합 {p:.4f}")
