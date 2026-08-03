"""r_sb² 성분 분해 판정 — F_g, F_P2, Cov(g,P2) 직접 측정.

항등식: s = g/S (g = 소수 파트너 수 = 골드바흐 계수, S = 생존자 수)
분해: r_sb² = (1−s̄)F_g + s̄F_P2 − 2Cov(g,P2)/S̄
예측(추측 3 유도 경로): F_g ≈ F_P2 ≈ ln2, Cov ≈ 0 ⇒ r_sb² = ln2.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

SCALE = 10_000_000
NS = 1200
is_p = sieve(SCALE + 4_000_000)
ps_all = primes_upto(SCALE // 2 + 2_000_000)
y = int(round(SCALE ** (1 / 3))) + 1
qs = [int(q) for q in primes_upto(y) if q > 2]
cands = [int(v) for v in np.arange(SCALE, SCALE + 3_000_000, 2)
         if v % 6 != 0][:NS]
G, P2, Ss = [], [], []
for i, n in enumerate(cands):
    m = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(m), dtype=bool)
    for q in qs:
        if n % q:
            al &= (m % q != 0)
    al &= m > 1
    surv = m[al]
    g = int(is_p[surv].sum())
    G.append(g)
    P2.append(len(surv) - g)
    Ss.append(len(surv))
    if (i + 1) % 300 == 0:
        print(f"  {i+1}/{NS}", flush=True)
G = np.array(G, float); P2 = np.array(P2, float); S = np.array(Ss, float)

# 국소 드리프트 제거 (스케일 3M 구간의 완만한 추세)
x = np.linspace(0, 1, len(G))
def dt(v):
    return v - np.polyval(np.polyfit(x, v, 2), x)
Gd, P2d, Sd = dt(G), dt(P2), dt(S)

F_g = float(Gd.var() / G.mean())
F_2 = float(P2d.var() / P2.mean())
cov = float(np.mean(Gd * P2d))
sbar = G.mean() / S.mean()
r2_pred = (1 - sbar) * F_g + sbar * F_2 - 2 * cov / S.mean()
# 직접 r_sb²
s = G / S
sd = dt(s)
r2_direct = float(sd.var() / (sbar * (1 - sbar) / S.mean()))

print(f"\nn ~ 10^7, {NS}표본 (2차 드리프트 제거):")
print(f"s̄ = {sbar:.5f} | ḡ = {G.mean():,.0f} | P̄2 = {P2.mean():,.0f}")
print(f"F_g (골드바흐 Fano)  = {F_g:.4f}")
print(f"F_P2 (P2 Fano)      = {F_2:.4f}")
print(f"corr(g, P2)         = {cov/np.sqrt(Gd.var()*P2d.var()):+.4f}")
print(f"\n분해 예측 r_sb² = {r2_pred:.4f}")
print(f"직접 측정 r_sb² = {r2_direct:.4f}")
print(f"ln 2 = {np.log(2):.4f}")
