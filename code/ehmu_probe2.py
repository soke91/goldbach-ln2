"""EH_mu 활력징후 v3 — 소수 모듈러스 전용 뫼비우스 등분포 지형.

D_p(θ) = Σ_{소수 q ≤ x^θ} max_a |M(x;q,a)| / x, 랜덤워크 기준 대비.
소수 q 1229개(≤10⁴) × bincount(0.5s) ≈ 12분.
"""

import math

import numpy as np

X = 100_000_000

print("mu 계산 중...", flush=True)
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
del val
print("mu 완료됨. 검증 Σμ(≤1e6) =",
      int(mu[:1_000_001].astype(np.int64).sum()), "(참값 212)", flush=True)

Qmax = 10_000
sq = np.ones(Qmax + 1, dtype=bool)
sq[:2] = False
for p in range(2, 101):
    if sq[p]:
        sq[p * p:: p] = False
prime_qs = [int(q) for q in np.nonzero(sq)[0] if q >= 3]
print(f"소수 q {len(prime_qs)}개 (≤{Qmax:,})", flush=True)

mu_f = mu.astype(np.float32)
del mu, pm
idx = np.arange(X + 1, dtype=np.int64)

cps = sorted(set(int(round(X ** t)) for t in [0.20, 0.30, 0.40, 0.45, 0.50]))
cp_th = {int(round(X ** t)): t for t in [0.20, 0.30, 0.40, 0.45, 0.50]}
tot = 0.0
rw = 0.0
worst = (0, 0, 0.0)
res = []
ci = 0
for j, q in enumerate(prime_qs):
    while ci < len(cps) and q > cps[ci]:
        th = cp_th[cps[ci]]
        print(f"θ={th:.2f}  Q={cps[ci]:,}  D_p={tot/X:.6f}  "
              f"RW_p={rw/X:.6f}  비 {tot/rw:.3f}  최악 q={worst[0]} "
              f"({worst[2]:.1f}×RW)", flush=True)
        res.append((th, cps[ci], tot / X, rw / X, tot / rw))
        np.savez("ehmu_probe_p.npz", res=np.array(res))
        ci += 1
    sums = np.bincount(idx % q, weights=mu_f, minlength=q).astype(np.int64)
    best = int(np.abs(sums).max())
    tot += best
    rw += math.sqrt(X / q)
    r = best / math.sqrt(X / q)
    if r > worst[2]:
        worst = (q, best, r)
    if (j + 1) % 100 == 0:
        print(f"  ...{j+1}/{len(prime_qs)} (q={q:,})", flush=True)
while ci < len(cps):
    th = cp_th[cps[ci]]
    print(f"θ={th:.2f}  Q={cps[ci]:,}  D_p={tot/X:.6f}  RW_p={rw/X:.6f}  "
          f"비 {tot/rw:.3f}  최악 q={worst[0]} ({worst[2]:.1f}×RW)",
          flush=True)
    res.append((th, cps[ci], tot / X, rw / X, tot / rw))
    ci += 1
np.savez("ehmu_probe_p.npz", res=np.array(res))
print("전체완료", flush=True)
