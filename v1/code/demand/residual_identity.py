"""13번째 사냥 — 잔차장 정체 판정: R(k) ≈ −Σ_{m 소수} μ(N−mk)?

ω_P(m)-슬라이스 D_j(k) = Σ_{ω_P(m)=j} μ(m)μ(N−mk) (j=0..5),
R(k) = D + c1·S 재구성 후:
  corr(R, D_0) — D_0 = 소수-m 슬라이스 (대역에서 ω=0 ⟺ m 소수)
  슬라이스별 에너지·R-기여 분해.
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
P0, P1 = 3, 3000
sv = np.ones(P1 + 1, dtype=bool)
sv[:2] = False
for q in range(2, int(P1 ** 0.5) + 1):
    if sv[q]:
        sv[q * q:: q] = False
Pset = np.array([int(q) for q in np.nonzero(sv)[0] if q >= P0])

Mmax = (N - 1) // 600
omega = np.zeros(Mmax + 1, dtype=np.int8)
for q in Pset:
    omega[q::q] += 1

ks = np.arange(600, 1400, 2)
JMAX = 5
Dslice = np.zeros((len(ks), JMAX + 1))
for i, k in enumerate(ks):
    k = int(k)
    ms = np.arange(SQ + 1, (N - 1) // k + 1, dtype=np.int64)
    muv = mu[ms].astype(np.float64)
    om = omega[ms]
    tail = mu[N - k * ms]
    base = muv * tail
    for j in range(JMAX + 1):
        Dslice[i, j] = base[om == j].sum()
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{len(ks)}", flush=True)

D = Dslice.sum(axis=1)
S = (Dslice * np.arange(JMAX + 1)).sum(axis=1)
c1 = -float((D * S).sum() / (S * S).sum())
R = D + c1 * S
print(f"\nc1 = {c1:+.4f}")
print(f"{'j(ω_P)':>7} {'에너지비 Σ|D_j|²/Σ|D|²':>22} {'corr(R, D_j)':>14} "
      f"{'가중(1+c1·j)':>12}")
for j in range(JMAX + 1):
    e = float((Dslice[:, j] ** 2).sum() / (D ** 2).sum())
    c = float(np.corrcoef(R, Dslice[:, j])[0, 1])
    print(f"{j:>7} {e:>22.3f} {c:>14.3f} {1+c1*j:>12.3f}")
print(f"\n**corr(R, D_0) = {np.corrcoef(R, Dslice[:,0])[0,1]:+.4f}**")
print(f"corr(R, Σ_j w_j D_j 재구성 검산) = "
      f"{np.corrcoef(R, (Dslice*(1+c1*np.arange(JMAX+1))).sum(axis=1))[0,1]:.4f}")
print(f"(대역 m ∈ ({SQ:,}, {Mmax:,}]에서 ω_P=0 ⟺ m 소수 — 3000² 초과 확인: "
      f"{3000*3000 > Mmax})")
np.savez("residual_identity.npz", Dslice=Dslice, R=R, ks=ks)
print("전체완료", flush=True)
