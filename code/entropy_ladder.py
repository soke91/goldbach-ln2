"""확대-엔트로피 사다리 검증 — p-관계 항등식의 정밀 실측.

항등식(주장): A_p(k) := Σ_{m∈(M0,M1], p|m} μ(m)μ(N−mk)
            = −Σ_{m'∈(M0/p,M1/p], p∤m'} μ(m')μ(N−m'(pk))
(m = pm', μ(pm') = −μ(m') for p∤m'; p²|m 항은 μ=0으로 자동 소거)
→ 정확 항등식이어야 함 (오차 0). + 사다리 강도: |A_p|/|dual| 비율
(엔트로피 논증이 쓰는 정보량 = 1/p-급 부분합의 크기).
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

def dual_range(k, M0, M1, coprime_p=None):
    if M1 <= M0:
        return None, 0
    ms = np.arange(M0 + 1, M1 + 1, dtype=np.int64)
    if coprime_p:
        ms = ms[ms % coprime_p != 0]
    vals = mu[ms].astype(np.int16) * mu[N - k * ms]
    return int(vals.sum(dtype=np.int64)), int(np.count_nonzero(vals))

rng = np.random.default_rng(163)
ks = sorted(int(v) for v in rng.integers(500, 1500, 40))
print(f"{'p':>3} {'항등식 최대오차':>12} {'평균 |A_p|/|dual|':>16} "
      f"{'corr(A_p, dual(pk))':>20}")
for p in [2, 3, 5]:
    errs = []
    ratios = []
    As = []
    Ds = []
    for k in ks:
        M0, M1 = SQ, (N - 1) // k
        # A_p(k): p|m 부분합
        ms = np.arange(M0 + 1, M1 + 1, dtype=np.int64)
        sel = ms[ms % p == 0]
        A = int((mu[sel].astype(np.int16) * mu[N - k * sel]).sum(
            dtype=np.int64))
        # 우변: −Σ_{m'∈(M0/p, M1/p], p∤m'} μ(m')μ(N−m'pk)
        B, _ = dual_range(p * k, M0 // p, M1 // p, coprime_p=p)
        errs.append(abs(A - (-(B if B is not None else 0))))
        D, v = dual_range(k, M0, M1)
        if D is not None and abs(D) > 0:
            ratios.append(abs(A) / max(abs(D), 1))
        As.append(A)
        full_pk, _ = dual_range(p * k, M0 // p, M1 // p)
        Ds.append(full_pk if full_pk is not None else 0)
    c = float(np.corrcoef(As, Ds)[0, 1])
    print(f"{p:>3} {max(errs):>12} {np.mean(ratios):>16.3f} {c:>20.3f}",
          flush=True)
print("\n오차 0 = 항등식 정확 (사다리 가로대는 엄밀). 비율 = 논증 연료량.")
print("corr(A_p, dual(pk)) ≈ −1 근접일수록 스케일-전달이 결정론적.")
print("전체완료", flush=True)
