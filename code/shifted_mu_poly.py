"""(b)-보조정리 직접 측정 — 이동-μ 디리클레 다항식의 t-상쇄.

G(t) = Σ_{j∈[J0,J1]} μ(N−j) j^{−it}.  측정: |G(t)|/√J_eff 의
t-그리드 (0 ~ T) 통계 — 평균·최악. MR-형 보조정리의 예측: 전 t에서
√-상쇄 (랜덤 다항식과 동일). 비교 대조: 계수를 무작위 부호로 바꾼
다항식의 동일 통계.
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
J0, J1 = 50_000_000, 52_000_000
js = np.arange(J0, J1, dtype=np.int64)
coef = mu[N - js].astype(np.float64)
J_eff = float(np.count_nonzero(coef))
lnj = np.log(js.astype(np.float64))
print(f"J = [{J0:,}, {J1:,})  유효 계수 {J_eff:,.0f}", flush=True)

rng = np.random.default_rng(167)
coef_rand = np.where(coef != 0, rng.choice([-1.0, 1.0], len(coef)), 0.0)

ts = np.linspace(0.0, 2000.0, 1200)
r_mu = np.empty(len(ts))
r_rd = np.empty(len(ts))
for i, t in enumerate(ts):
    ph = np.exp(-1j * t * lnj)
    r_mu[i] = abs(np.dot(coef, ph)) / math.sqrt(J_eff)
    r_rd[i] = abs(np.dot(coef_rand, ph)) / math.sqrt(J_eff)
    if (i + 1) % 200 == 0:
        print(f"  {i+1}/{len(ts)}", flush=True)
        np.savez("shifted_mu_poly.npz", ts=ts[:i+1], r_mu=r_mu[:i+1],
                 r_rd=r_rd[:i+1])

print(f"\n[이동-μ 다항식 |G(t)|/√J, t ∈ [0, 2000], 1200점]")
print(f"μ(N−j) 계수:   평균 {r_mu.mean():.3f}  중앙값 "
      f"{np.median(r_mu):.3f}  최악 {r_mu.max():.3f}  (t={ts[np.argmax(r_mu)]:.1f})")
print(f"무작위 대조:   평균 {r_rd.mean():.3f}  중앙값 "
      f"{np.median(r_rd):.3f}  최악 {r_rd.max():.3f}")
print(f"레일리 기준: 평균 ≈ 0.886. μ-열이 무작위와 동행하면")
print(f"(b)-보조정리(이동-μ MR 평균값)의 자연-측 참 확정")
np.savez("shifted_mu_poly.npz", ts=ts, r_mu=r_mu, r_rd=r_rd)
print("전체완료", flush=True)
