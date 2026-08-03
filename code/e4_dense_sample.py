"""고옥타브 표본 증강 — k = 25..28, 각 8000 표본 (잡음 ±0.008 → ±0.005)."""

from math import log

import numpy as np

from goldbach.sieve import primes_upto, sieve
from goldbach.stats import hardy_littlewood_estimate

LIMIT = 540_000_000
is_p = sieve(LIMIT)
ps = primes_upto(LIMIT // 2)
print(f"체 완료 {LIMIT:,}", flush=True)

hist = {12: 0.5062, 13: 0.4914, 14: 0.4961, 15: 0.5317, 16: 0.5296,
        17: 0.5372, 18: 0.5507, 19: 0.5573, 20: 0.5652, 21: 0.5691,
        22: 0.5708, 23: 0.5784, 24: 0.5900}

rows = [(np.sqrt(2.0 ** k * 2.0 ** (k + 1)), v, 0.008) for k, v in hist.items()]

for k in [25, 26, 27, 28]:
    lo, hi = 2 ** k, 2 ** (k + 1)
    base = np.linspace(lo, min(hi, LIMIT), 8000).astype(np.int64)
    nvals = base + (2 - base % 6) % 6
    zs, hls = [], []
    for n in nvals:
        n = int(n)
        sub = ps[: np.searchsorted(ps, n // 2, side="right")]
        g = int(is_p[n - sub].sum())
        hl = hardy_littlewood_estimate(n)
        zs.append(g / hl)
        hls.append(hl)
    zs, hls = np.array(zs), np.array(hls)
    x = np.log(nvals.astype(float))
    resid = zs - np.polyval(np.polyfit(x, zs, 1), x)
    r = float(resid.std() / np.mean(1 / np.sqrt(hls)))
    se = r / np.sqrt(2 * len(nvals))
    rows.append((np.sqrt(float(lo) * hi), r, se))
    print(f"2^{k}: {r:.4f} ± {se:.4f}", flush=True)

mid = np.array([q[0] for q in rows])
arr = np.array([q[1] for q in rows])
w = 1 / np.array([q[2] for q in rows]) ** 2
x = 1 / np.log(mid)
# 가중 최소제곱 (말단 창별)
for tail in [len(rows), 10, 8]:
    xs, ys, ws = x[-tail:], arr[-tail:], w[-tail:]
    A = np.vstack([np.ones(tail), xs]).T
    W = np.diag(ws)
    beta, cov = np.linalg.solve(A.T @ W @ A, A.T @ W @ ys), np.linalg.inv(A.T @ W @ A)
    a, se_a = beta[0], np.sqrt(cov[0, 0])
    print(f"말단 {tail:>2}옥타브 가중적합 극한 = {a:.4f} ± {se_a:.4f}")
print(f"ln 2 = {log(2):.4f}")
