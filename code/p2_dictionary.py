"""프레임 사전 연결 — g/S/P2의 구조-제거 Fano 삼중 측정 (10⁷, 1200표본).

구조 모델: g는 𝔖₀(n) = ∏_{q|n,q>2}(q−1)/(q−2), S는 W(n) = ∏_{3≤q≤y, q∤n}(1−1/(q−1)).
각각 x(n)/구조(n)를 ln n 2차 다항으로 평활 적합 → 잔차를 카운트 단위로 환원 →
F̃ = Var(잔차)/평균카운트.
접속 확인: F̃_g ≈ comet r²(10⁷ 스케일 ≈ 0.58² ≈ 0.34)?
조합 검증: r_sb² =? (1−s̄)F̃_g + s̄F̃_P2 − 2Coṽ/S̄.
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

SCALE = 10_000_000
NS = 1200
HI = int(SCALE * 1.3)
is_p = sieve(HI + 100_000)
ps_all = primes_upto(HI // 2 + 50_000)
y = int(round(SCALE ** (1 / 3))) + 1
qs = [int(q) for q in primes_upto(y) if q > 2]
base = np.linspace(SCALE, HI - 10, NS).astype(np.int64)
cands = [int(v + (2 - v % 6) % 6) for v in base]

G, S = [], []
for i, n in enumerate(cands):
    m = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(m), dtype=bool)
    for q in qs:
        if n % q:
            al &= (m % q != 0)
    al &= m > 1
    surv = m[al]
    G.append(int(is_p[surv].sum()))
    S.append(len(surv))
    if (i + 1) % 300 == 0:
        print(f"  {i+1}/{NS}", flush=True)
G = np.array(G, float); S = np.array(S, float)
P2 = S - G
narr = np.array(cands, float)

# 구조 모델
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
Sg = np.array([sing0(int(n)) for n in cands])
W = np.ones(NS)
for q in qs:
    mask = (np.array(cands) % q) != 0
    W[mask] *= 1 - 1 / (q - 1)

x = np.log(narr); x = (x - x.mean()) / x.std()
def resid_counts(cnt, struct):
    z = cnt / struct
    fit = np.polyval(np.polyfit(x, z, 2), x)
    return (z - fit) * struct  # 카운트 단위 잔차

eg = resid_counts(G, Sg)
eS = resid_counts(S, W)
e2 = resid_counts(P2, W)   # P2 구조 ≈ S 구조 (1차 근사) — 잔차로 검증

Fg = float(eg.var() / G.mean())
FS = float(eS.var() / S.mean())
F2 = float(e2.var() / P2.mean())
sbar = G.mean() / S.mean()
cov_g2 = float(np.mean(eg * e2))
r2_comb = (1 - sbar) * Fg + sbar * F2 - 2 * cov_g2 / S.mean()

print(f"\nn ~ 10⁷ (×1.3 창), {NS}표본:")
print(f"F̃_g (𝔖-정규화)  = {Fg:.4f}   [comet r² 예상 ≈ 0.34]")
print(f"F̃_S (W-정규화)  = {FS:.4f}   [신규: 생존자-수 Fano]")
print(f"F̃_P2            = {F2:.4f}")
print(f"corr(ẽ_g, ẽ_P2) = {cov_g2/np.sqrt(eg.var()*e2.var()):+.4f}")
print(f"조합 r_sb² 예측  = {r2_comb:.4f}  [직접 측정 ≈ 0.64~0.78]")
print(f"참고: ln2 = 0.6931, (ln2)² = 0.4805")
