"""r_sb 교정 측정 — 상대 창(×1.3) + 일관된 2차 드리프트 제거 (함정 16 수정판)."""

import numpy as np

from goldbach.sieve import primes_upto, sieve

print(f"{'스케일':>13} {'표본':>5} {'r_sb(교정)':>10} {'±':>7}")
for SCALE, NS in [(3_000_000, 1200), (10_000_000, 1000), (30_000_000, 500)]:
    HI = int(SCALE * 1.3)
    is_p = sieve(HI + 100_000)
    ps_all = primes_upto(HI // 2 + 50_000)
    y = int(round(SCALE ** (1 / 3))) + 1
    qs = [int(q) for q in primes_upto(y) if q > 2]
    base = np.linspace(SCALE, HI - 10, NS).astype(np.int64)
    cands = [int(v + (2 - v % 6) % 6) for v in base]  # n ≡ 2 (mod 6)
    G, S = [], []
    for n in cands:
        m = n - ps_all[ps_all <= n // 2]
        al = np.ones(len(m), dtype=bool)
        for q in qs:
            if n % q:
                al &= (m % q != 0)
        al &= m > 1
        surv = m[al]
        G.append(int(is_p[surv].sum()))
        S.append(len(surv))
    G = np.array(G, float); S = np.array(S, float)
    s = G / S
    x = np.linspace(0, 1, len(s))
    sd = s - np.polyval(np.polyfit(x, s, 2), x)
    sbar = s.mean()
    r2 = float(sd.var() / (sbar * (1 - sbar) / S.mean()))
    r = np.sqrt(r2)
    se = r / np.sqrt(2 * len(s))
    print(f"{SCALE:>13,} {len(s):>5} {r:>10.4f} {se:>7.4f}", flush=True)
print(f"\n√ln2 = {np.sqrt(np.log(2)):.4f}")
