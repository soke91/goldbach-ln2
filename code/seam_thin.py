"""SEAM 추측의 정확한 객체 — 얇은-수열 클래스별 μ-합 직접 실측.

이음새: 모듈러스 L = [k,k'] ~ K² (K ∈ (x^{0.3}, x^{1/3}]), 범위 y ~ x
클래스당 원소 ~ y/L (얇음: L > √y). 객체:
  S(L, a) = Σ_{w ≡ a (L), w ≤ y} μ(w)   (클래스별 μ-합)
실측: L ~ 10⁵ (모델), y = 10⁸ → 클래스당 ~10³ 원소.
r(L,a) = |S|/√(y/L·6/π²) — 얇은-수열에서도 반정규인가 (SEAM의 내용).
표본: L 40개 × 클래스 25개.
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

rng = np.random.default_rng(223)
rows = []
for i in range(40):
    L = int(rng.integers(80_000, 200_000))    # L ~ x^{0.6} 모델
    for _ in range(25):
        a = int(rng.integers(0, L))
        vals = mu[a::L].astype(np.int64)
        S = int(vals.sum())
        n_eff = int(np.count_nonzero(vals))
        r = abs(S) / math.sqrt(max(n_eff, 1))
        rows.append((L, a, S, n_eff, r))
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/40", flush=True)
        np.savez("seam_thin.npz", rows=np.array(rows))

d = np.array(rows)
r = d[:, 4]
print(f"\n[얇은-수열 μ-합, L~10⁵ (클래스당 ~600 유효원소), 1000클래스]")
print(f"r = |S|/√n_eff: 평균 {r.mean():.3f}  중앙값 {np.median(r):.3f}  "
      f"최악 {r.max():.3f}  (반정규 0.798)")
print(f"r > 3: {int((r > 3).sum())}개 (기대 ~{0.0027 * len(r):.1f})")
np.savez("seam_thin.npz", rows=d)
print("전체완료", flush=True)
