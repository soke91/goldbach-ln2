"""F_S (생존자 총수 Fano)의 스케일 안정성 스캔 — 상수 판정.

동일 설계(×1.3 창, n≡2(6), 2차 드리프트, W-구조 정규화)로
3×10⁶ / 10⁷ / 3×10⁷ 에서 F̃_S, F̃_g, F̃_P2 삼중 추적.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

def sing0(n):
    c, m_, p = 1.0, n, 3
    while m_ % 2 == 0: m_ //= 2
    while p * p <= m_:
        if m_ % p == 0:
            c *= (p - 1) / (p - 2)
            while m_ % p == 0: m_ //= p
        p += 2
    if m_ > 1: c *= (m_ - 1) / (m_ - 2)
    return c

print(f"{'스케일':>13} {'F̃_g':>7} {'F̃_S':>7} {'F̃_P2':>7} {'corr':>7}")
for SCALE, NS in [(3_000_000, 5000), (10_000_000, 2500)]:
    HI = int(SCALE * 1.3)
    is_p = sieve(HI + 100_000)
    ps_all = primes_upto(HI // 2 + 50_000)
    y = int(round(SCALE ** (1 / 3))) + 1
    qs = [int(q) for q in primes_upto(y) if q > 2]
    base = np.linspace(SCALE, HI - 10, NS).astype(np.int64)
    cands = [int(v + (2 - v % 6) % 6) for v in base]
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
    G = np.array(G, float); S = np.array(S, float); P2 = S - G
    Sg = np.array([sing0(int(n)) for n in cands])
    W = np.ones(len(cands))
    for q in qs:
        mask = (np.array(cands) % q) != 0
        W[mask] *= 1 - 1 / (q - 1)
    x = np.linspace(0, 1, len(G))
    def rc(cnt, st):
        z = cnt / st
        return (z - np.polyval(np.polyfit(x, z, 2), x)) * st
    eg, eS, e2 = rc(G, Sg), rc(S, W), rc(P2, W)
    Fg = float(eg.var() / G.mean())
    FS = float(eS.var() / S.mean())
    F2 = float(e2.var() / P2.mean())
    cr = float(np.mean(eg * e2) / np.sqrt(eg.var() * e2.var()))
    print(f"{SCALE:>13,} {Fg:>7.3f} {FS:>7.3f} {F2:>7.3f} {cr:>+7.3f}",
          flush=True)
