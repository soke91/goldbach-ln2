"""캐스케이드 효율 측정 — 전 대역 L²에서 사다리-합의 설명력.

k ∈ [600, 1400] 전수. D(k), S(k) = Σ_{p∈P} A_p(k).
  (1) Σ|E|²/Σ|D|²  (E = D + S/L — 소박 전개)
  (2) 회귀: c* = argmin Σ|D + cS|², 잔차비 = 1 − R² — 캐스케이드
      1스텝의 실효 에너지 전달률 R².
P 창 2종 (좁은 [50,500] / 넓은 [3,3000])으로 L-의존 확인.
"""

import math

import numpy as np

X = 100_000_000

print("mu 계산...", flush=True)
mu = np.ones(X + 1, dtype=np.int8)
pm = np.ones(X + 1, dtype=bool)
pm[:2] = False
for p in range(2, int(X ** 0.5) + 1):
    if pm[p]:
        pm[p * p:: p] = False
        mu[p::p] *= -1
        mu[p * p:: p * p] = 0
val = np.arange(X + 1, dtype=np.int64)
for p in range(2, int(X ** 0.5) + 1):
    if pm[p]:
        val[p::p] //= p
        pp = p * p
        while pp <= X:
            val[pp::pp] //= p
            pp *= p
mu[val > 1] *= -1
mu[0] = 0
del val, pm
print("mu 완료됨", flush=True)

N = 99_999_998
SQ = int(N ** 0.5)

def primes_in(a, b):
    sv = np.ones(b + 1, dtype=bool)
    sv[:2] = False
    for q in range(2, int(b ** 0.5) + 1):
        if sv[q]:
            sv[q * q:: q] = False
    return [int(q) for q in np.nonzero(sv)[0] if q >= a]

ks = np.arange(600, 1400)
for P0, P1 in [(50, 500), (3, 3000)]:
    Pset = primes_in(P0, P1)
    L = sum(1.0 / q for q in Pset)
    Ds = np.empty(len(ks))
    Ss = np.empty(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        M0, M1 = SQ, (N - 1) // k
        ms = np.arange(M0 + 1, M1 + 1, dtype=np.int64)
        bv = mu[ms].astype(np.float64) * mu[N - k * ms]
        Ds[i] = bv.sum()
        S = 0.0
        for q in Pset:
            sel = ms[ms % q == 0]
            S += float((mu[sel].astype(np.float64)
                        * mu[N - k * sel]).sum())
        Ss[i] = S
        if (i + 1) % 100 == 0:
            print(f"  P[{P0},{P1}]: {i+1}/{len(ks)}", flush=True)
    E = Ds + Ss / L
    r2_naive = float((E ** 2).sum() / (Ds ** 2).sum())
    cstar = -float((Ds * Ss).sum() / (Ss ** 2).sum())
    resid = Ds + (-cstar) * (-Ss)
    R2 = 1 - float((resid ** 2).sum() / (Ds ** 2).sum())
    corr = float(np.corrcoef(Ds, Ss)[0, 1])
    print(f"\nP=[{P0},{P1}]  L={L:.3f}  |k대역|={len(ks)}")
    print(f"  소박 Σ|E|²/Σ|D|² = {r2_naive:.3f}  (1/L 계수)")
    print(f"  corr(D, S) = {corr:+.3f}   최적 c* = {cstar:+.4f} "
          f"(소박 1/L = {1/L:.4f})")
    print(f"  **캐스케이드 1스텝 에너지 전달률 R² = {R2:.3f}**", flush=True)
    np.savez(f"cascade_eff_{P0}_{P1}.npz", D=Ds, S=Ss, ks=ks)
print("전체완료", flush=True)
