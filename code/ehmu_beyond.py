"""EH_mu 미답 지대 (θ > 1/2) 표본 정찰 — 추측의 실질 영역.

로그-간격 소수 q ∈ (10⁴, 10⁷] 표본 240개에서 r(q) = max_a|M(x;q,a)| /
√(x/q) 분포 측정. θ = lnq/lnx 축으로 프로파일 — 1/2 경계를 넘어도
로그-팽창 랜덤워크가 유지되는가 (EH_μ의 미증명 코어 직접 관측).
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
del val, pm
mu_f = mu.astype(np.float32)
del mu
idx = np.arange(X + 1, dtype=np.int64)

def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

rng = np.random.default_rng(113)
targets = np.exp(np.linspace(math.log(1.2e4), math.log(1e7), 240))
rows = []
for i, t in enumerate(targets):
    q = int(t) | 1
    while not is_prime(q):
        q += 2
    sums = np.bincount(idx % q, weights=mu_f, minlength=q).astype(np.int64)
    r = float(np.abs(sums).max() / math.sqrt(X / q))
    th = math.log(q) / math.log(X)
    rows.append((q, th, r))
    if (i + 1) % 24 == 0:
        print(f"  {i+1}/240  q={q:,} θ={th:.3f} r={r:.2f}", flush=True)
        np.savez("ehmu_beyond.npz", rows=np.array(rows))

d = np.array(rows)
print("\n[θ-구간별 대-RW 비 r(q)]")
for lo, hi in [(0.5, 0.55), (0.55, 0.6), (0.6, 0.65), (0.65, 0.7),
               (0.7, 0.75), (0.75, 0.8), (0.8, 0.85), (0.85, 0.88)]:
    mk = (d[:, 1] >= lo) & (d[:, 1] < hi)
    if mk.sum():
        print(f"θ {lo:.2f}-{hi:.2f}: 평균 {d[mk,2].mean():.2f}  "
              f"최악 {d[mk,2].max():.2f}  (표본 {int(mk.sum())})")
np.savez("ehmu_beyond.npz", rows=d)
print("전체완료", flush=True)
