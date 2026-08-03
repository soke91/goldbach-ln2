"""② 프레임 상수 대캠페인 — 7스케일 다-옥타브, 극한 적합용 시리즈.

F_g / F_S / F_P2 를 1e6 ~ 6.4e7 (×2 간격 7스케일)에서 동일 설계로 측정.
각 스케일 표본수는 비용 역비례. 최종: 1/ln n 적합 극한 3종.
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

rows = []
for SCALE, NS in [(1_000_000, 4000), (2_000_000, 3000), (4_000_000, 2200),
                  (8_000_000, 1600), (16_000_000, 1000), (32_000_000, 600),
                  (64_000_000, 350)]:
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
    rows.append((SCALE, NS, Fg, FS, F2))
    print(f"{SCALE:>12,} ({NS}) F_g={Fg:.4f} F_S={FS:.4f} F_P2={F2:.4f}",
          flush=True)

mid = np.array([r[0] for r in rows], float)
x1 = 1 / np.log(mid)
print("\n[1/ln n 가중적합 극한]")
for name, idx in [("F_g", 2), ("F_S", 3), ("F_P2", 4)]:
    vals = np.array([r[idx] for r in rows])
    w = np.array([r[1] for r in rows], float)  # 표본수 가중
    A = np.vstack([np.ones(len(x1)), x1]).T
    Wm = np.diag(w)
    beta = np.linalg.solve(A.T @ Wm @ A, A.T @ Wm @ vals)
    print(f"  {name}: 극한 = {beta[0]:.4f} (기울기 {beta[1]:+.3f})")
print(f"\n(ln2)² = {np.log(2)**2:.4f} | ln2 = {np.log(2):.4f}")
