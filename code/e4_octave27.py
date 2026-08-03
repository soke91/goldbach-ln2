"""정준(comet) 시리즈 확장 — 옥타브 2^26~2^27 (표본 3200개).

r = 순수요동 초과비 (HL1 기준, n≡2 mod 6, 옥타브 내 1차 드리프트 제거).
기존 16옥타브와 결합해 말단 적합 극한 재계산.
"""

from math import log

import numpy as np

from goldbach.sieve import primes_upto, sieve
from goldbach.stats import hardy_littlewood_estimate

LIMIT = 270_000_000
is_p = sieve(LIMIT)
ps = primes_upto(LIMIT // 2)
print(f"체 완료 {LIMIT:,} / 소수 {len(ps):,}", flush=True)

hist = {12: 0.5062, 13: 0.4914, 14: 0.4961, 15: 0.5317, 16: 0.5296,
        17: 0.5372, 18: 0.5507, 19: 0.5573, 20: 0.5652, 21: 0.5691,
        22: 0.5708, 23: 0.5784, 24: 0.5900, 25: 0.5839, 26: 0.5967}

rng = np.random.default_rng(27)
rows = [(np.sqrt(2.0 ** k * 2.0 ** (k + 1)), v) for k, v in hist.items()]

for k in [27]:
    lo, hi = 2 ** k, 2 ** (k + 1)
    base = np.linspace(lo, min(hi, LIMIT), 3200).astype(np.int64)
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
    rows.append((np.sqrt(float(lo) * hi), r))
    print(f"2^{k}: 순수 초과비 = {r:.4f} ± {se:.4f}", flush=True)

mid = np.array([q[0] for q in rows])
arr = np.array([q[1] for q in rows])
x = 1 / np.log(mid)
for tail in [len(rows), 10, 8, 6]:
    b, a = np.polyfit(x[-tail:], arr[-tail:], 1)
    print(f"말단 {tail:>2}옥타브 적합 극한 = {a:.4f}")
print(f"ln 2 = {log(2):.4f}")
