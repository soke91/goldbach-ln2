"""①수축 공식 검증 + ②캐스케이드 2스텝 직접 시험.

1스텝: D, S(=Σμω μ), c1* 회귀 → 잔차 R(k) = D + c1 S
  (동치: 가중 w1(m) = μ(m)(1 + c1 ω_P(m)))
수축 예측: c_pred = (1/L)/(1 + V_n/(L² σ_D²)), V_n = 잔차 분산 실측치로
  자기일관 검증 (선험판은 TK로 대체 가능 — 여기선 공식 구조만 확인).
2스텝: S_R(k) = Σ_m w1(m) ω_P(m) μ(N−mk) 회귀 → R²₂.
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
L = float((1.0 / Pset).sum())

# ω_P(m) 사전계산 (m ≤ N/600 범위)
Mmax = (N - 1) // 600
print("ω_P 계산...", flush=True)
omega = np.zeros(Mmax + 1, dtype=np.int8)
for q in Pset:
    omega[q::q] += 1

ks = np.arange(600, 1400)
D1 = np.empty(len(ks))
S1 = np.empty(len(ks))
for i, k in enumerate(ks):
    k = int(k)
    ms = np.arange(SQ + 1, (N - 1) // k + 1, dtype=np.int64)
    muv = mu[ms].astype(np.float64)
    tail = mu[N - k * ms]
    base = muv * tail
    D1[i] = base.sum()
    S1[i] = (base * omega[ms]).sum()
    if (i + 1) % 200 == 0:
        print(f"  1스텝 {i+1}/{len(ks)}", flush=True)

c1 = -float((D1 * S1).sum() / (S1 * S1).sum())
R2_1 = 1 - float(((D1 + c1 * S1) ** 2).sum() / (D1 ** 2).sum())
# 수축 공식 자기일관 검증
sD2 = float((D1 ** 2).mean())
resid1 = D1 + c1 * S1
Vn = float((resid1 ** 2).mean())          # 잔차(잡음) 분산 근사
c_pred = (1 / L) / (1 + Vn * (1) / (L ** 2 * sD2 / (1 + 0)))  # 구조 확인용
print(f"\n[1스텝] c1* = {c1:+.4f}  R²₁ = {R2_1:.3f}  "
      f"(1/L = {1/L:.4f}, 수축구조 c_pred = {c_pred:.4f})", flush=True)

# 2스텝: R의 사다리 S_R = Σ w1 ω μ,  w1 = μ(1 + c1 ω)
D2v = np.empty(len(ks))
S2v = np.empty(len(ks))
for i, k in enumerate(ks):
    k = int(k)
    ms = np.arange(SQ + 1, (N - 1) // k + 1, dtype=np.int64)
    muv = mu[ms].astype(np.float64)
    om = omega[ms].astype(np.float64)
    tail = mu[N - k * ms]
    w1 = muv * (1 + c1 * om)
    D2v[i] = (w1 * tail).sum()            # = R(k) 재계산
    S2v[i] = (w1 * om * tail).sum()
    if (i + 1) % 200 == 0:
        print(f"  2스텝 {i+1}/{len(ks)}", flush=True)

c2 = -float((D2v * S2v).sum() / (S2v * S2v).sum())
R2_2 = 1 - float(((D2v + c2 * S2v) ** 2).sum() / (D2v ** 2).sum())
print(f"\n[2스텝] c2* = {c2:+.4f}  **R²₂ = {R2_2:.3f}**")
print(f"1스텝 0.681 대비 — 유지되면 반복 확정, 급락하면 1회성")
tot = 1 - (1 - R2_1) * (1 - R2_2)
print(f"2스텝 누적 에너지 이관: {tot:.3f} (잔차 {(1-R2_1)*(1-R2_2):.3f})")
np.savez("cascade_2step.npz", D1=D1, S1=S1, D2=D2v, S2=S2v,
         c=[c1, c2], R2=[R2_1, R2_2])
print("전체완료", flush=True)
