"""TK/Ramaré 1-스텝 전개 검증 — 확대-캐스케이드의 오차항 실측.

전개(주장): D(k) ≈ −(1/L) Σ_{p∈P} A_p(k),  L = Σ_{p∈P} 1/p.
오차: E(k) = D(k) + (1/L)Σ_p A_p(k)  (TK 예측: |E| ≪ √M·(1/√(L·정규화)))
P = [P0, P1] 중간 소수. 측정: |E(k)|/√M vs |D(k)|/√M — 오차가 신호보다
충분히 작으면 캐스케이드 1스텝 성립 → 반복 가능성.
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

# 중간 소수 집합 P
P0, P1 = 50, 500
sv = np.ones(P1 + 1, dtype=bool)
sv[:2] = False
for q in range(2, int(P1 ** 0.5) + 1):
    if sv[q]:
        sv[q * q:: q] = False
Pset = [int(q) for q in np.nonzero(sv)[0] if q >= P0]
L = sum(1.0 / q for q in Pset)
print(f"P = [{P0}, {P1}]  |P| = {len(Pset)}  L = {L:.4f}", flush=True)

rng = np.random.default_rng(173)
ks = sorted(int(v) for v in rng.integers(600, 1400, 30))
rows = []
for k in ks:
    M0, M1 = SQ, (N - 1) // k
    ms = np.arange(M0 + 1, M1 + 1, dtype=np.int64)
    base_vals = mu[ms].astype(np.float64) * mu[N - k * ms]
    D = float(base_vals.sum())
    M_eff = float(np.count_nonzero(base_vals))
    S = 0.0
    for q in Pset:
        sel = ms[ms % q == 0]
        S += float((mu[sel].astype(np.float64) * mu[N - k * sel]).sum())
    E = D + S / L
    rows.append((k, D / math.sqrt(M_eff), E / math.sqrt(M_eff)))
d = np.array(rows)
print(f"\n[1-스텝 전개, k 30개]  (√M 정규화)")
print(f"|D|:  평균 {np.abs(d[:,1]).mean():.3f}  최악 {np.abs(d[:,1]).max():.3f}")
print(f"|E|:  평균 {np.abs(d[:,2]).mean():.3f}  최악 {np.abs(d[:,2]).max():.3f}")
print(f"오차/신호 비: {np.abs(d[:,2]).mean()/max(np.abs(d[:,1]).mean(),1e-9):.3f}")
print(f"TK 예상 오차 규모 ~ 1/√(L·보정) = {1/math.sqrt(L):.3f}")
np.savez("tk_expansion.npz", rows=d)
print("전체완료", flush=True)
