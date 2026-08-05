"""E1 전제 검증 — 확대-족의 엔트로피-이동 상관: corr_k(dual(k), dual(pk)).

Tao 엔트로피-감소의 확대판 전제: dual(k)와 dual(pk)가 μ(p)-관계로
정보를 공유. p = 2, 3, 5, 7에서 k-대역 상관 실측 + 무작위 대조.
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

def dual(k):
    M0 = SQ + 1
    M1 = (N - 1) // k
    if M1 <= M0:
        return None
    ms = np.arange(M0, M1 + 1, dtype=np.int64)
    vals = mu[ms].astype(np.int16) * mu[N - k * ms]
    v = int(np.count_nonzero(vals))
    return int(vals.sum(dtype=np.int64)) / math.sqrt(max(v, 1))

ks = list(range(500, 1900))
base = {}
for k in ks:
    d = dual(k)
    if d is not None:
        base[k] = d
print(f"기저 dual: {len(base)}개 (k 500~1900)", flush=True)

for p in [2, 3, 5, 7]:
    pairs = [(base[k], dual(p * k)) for k in ks
             if k in base and p * k <= SQ]
    pairs = [(a, b) for a, b in pairs if b is not None]
    a = np.array([x for x, _ in pairs])
    b = np.array([y for _, y in pairs])
    if len(a) > 30:
        c = float(np.corrcoef(a, b)[0, 1])
        # 대조: 무작위 재쌍
        rng = np.random.default_rng(157 + p)
        c0 = float(np.corrcoef(a, rng.permutation(b))[0, 1])
        print(f"p={p}: corr(dual(k), dual({p}k)) = {c:+.4f}  "
              f"(쌍 {len(a)}, 무작위 대조 {c0:+.4f})", flush=True)
print("전체완료", flush=True)
