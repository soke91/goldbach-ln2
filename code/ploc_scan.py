"""P_loc 돌파-표적 전수 검증 — theta=1/3, 집계량 2개의 최악 편차.

g(n) = S(n) − P2(n) (정리 1 정확). 승리 조건: 개별 n에서
  |S/S_model − 1|·1 + |P2/P2_model − 1|·(P2*/S*) < g*/S* = 0.2226
즉 각 집계가 ~12% 이내면 충분. 전수 창에서 최악 편차 실측.

모형: S_model(n) = pi(n/2)·∏_{2<q≤y, q∤n}(q−2)/(q−1)·e^γ·ω(3)·보정 —
대신 **자기-보정 모형**: 창 전체 평균 형상 × n별 𝔖-인자 (법칙 완전
상쇄 구조 이용, 58차분). 편차 = 개별 n의 형상-모형 잔차.
주의: 자기-보정은 '창 평균을 앵커 1점으로 쓰는' P_loc 구조 그대로
(앵커+전이 아키텍처의 실측판).
"""

import numpy as np

from goldbach.sieve import primes_upto, sieve

X = 100_000_000
W = 3_000_000
is_p = sieve(X + W + 10)
ps_all = primes_upto((X + W) // 2 + 10)
y = int(round(X ** (1 / 3))) + 1
qs = [int(q) for q in primes_upto(y) if q > 2]

cands = [int(n) for n in range(X, X + W, 2) if n % 6 == 2]
M = len(cands)
print(f"전수 {M:,}개 @ 10⁸ (창 3e6)", flush=True)

def sing_ratio(n):
    # n별 𝔖-型 인자: ∏_{q≤y, q|n} (q−1)/(q−2)  (S·P2 공통 모드)
    c = 1.0
    for q in qs:
        if n % q == 0:
            c *= (q - 1) / (q - 2)
    return c

S_v = np.empty(M)
P2_v = np.empty(M)
Sg = np.empty(M)
for i, n in enumerate(cands):
    m = n - ps_all[ps_all <= n // 2]
    al = np.ones(len(m), dtype=bool)
    for q in qs:
        if n % q:
            al &= (m % q != 0)
    al &= m > 1
    surv = m[al]
    S_v[i] = len(surv)
    P2_v[i] = len(surv) - int(is_p[surv].sum())
    Sg[i] = sing_ratio(n)
    if (i + 1) % 50000 == 0:
        print(f"  {i+1}/{M}", flush=True)
    if (i + 1) % 100000 == 0:
        np.savez("ploc_scan_partial.npz", S=S_v[:i+1], P2=P2_v[:i+1],
                 Sg=Sg[:i+1], cands=np.array(cands[:i+1]))

# 자기-보정 모형: ln n 선형 드리프트 + 𝔖-인자
ln = np.log(np.array(cands, dtype=float))
lc = ln - ln.mean()
def devs(V):
    base = V / Sg                      # 공통 𝔖-모드 제거
    A = np.vstack([np.ones(M), lc]).T
    coef, *_ = np.linalg.lstsq(A, np.log(base), rcond=None)
    model = np.exp(A @ coef) * Sg
    return V / model - 1

dS = devs(S_v)
dP2 = devs(P2_v)
gstar_ratio = 1 - (P2_v / S_v).mean()   # g*/S* 실측 (~0.39 at 1e8 유한크기)
combo = np.abs(dS) + np.abs(dP2) * (P2_v / S_v).mean()
print(f"\n[P_loc 원장 @10⁸ 전수]  g*/S* = {gstar_ratio:.4f} (승리 문턱)")
print(f"S-편차:  σ {dS.std():.5f}  최악 {np.abs(dS).max():.5f}")
print(f"P2-편차: σ {dP2.std():.5f}  최악 {np.abs(dP2).max():.5f}")
print(f"결합 최악 {combo.max():.5f}  vs 문턱 {gstar_ratio:.4f}  "
      f"→ **슬랙 {gstar_ratio/combo.max():.0f}배**")
n_worst = cands[int(np.argmax(combo))]
print(f"최악 n = {n_worst}")
np.savez("ploc_scan.npz", S=S_v, P2=P2_v, Sg=Sg, cands=np.array(cands),
         dS=dS, dP2=dP2)
