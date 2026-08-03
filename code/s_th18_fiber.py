"""theta=1/8 섬유 균일성 — 과녁 지대의 chi2 (상금 질문).

대조군 n (중간 인수 없음) 25개 @ 10^10. 섬유 q in (19, 300]:
F_q(n) = #{합성 생존자: q | m} (SKIP 세그먼트 표본).
chi2/(K-1) vs 포아송 1.0 — theta=1/3의 sub-Poisson(0.52~0.56)이
과녁 거칠기에서도 성립하는지.
"""

import numpy as np

from goldbach.sieve import primes_upto

X = 10_000_000_000
SEG = 4_000_000
SKIP = 20
M = 25

def mr_bases(n):
    if n < 3_215_031_751:
        return (2, 3, 5, 7)
    if n < 341_550_071_728_321:
        return (2, 3, 5, 7, 11, 13, 17)
    return (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

def is_prime_u64(n):
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in mr_bases(n):
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

med_primes = [int(p) for p in primes_upto(2000) if p > 19]
qs_fiber = np.array([int(p) for p in primes_upto(300) if p > 19])
K = len(qs_fiber)

def no_med(n):
    return all(n % p for p in med_primes)

rng = np.random.default_rng(61)
base = X + (2 - X % 6) % 6
cands = []
while len(cands) < M:
    n = int(base + 6 * int(rng.integers(0, X // 200 // 6)))
    if no_med(n):
        cands.append(n)

y = int(round(X ** (1 / 8)))
small = [int(p) for p in primes_upto(y + 1) if p > 2]
base_primes = primes_upto(int((X // 2 + SEG) ** 0.5) + 1)

def primes_in_segment(lo, hi):
    if lo < 3:
        lo = 3
    odd_lo = lo if lo % 2 else lo + 1
    sz = (hi - odd_lo + 1) // 2
    seg = np.ones(sz, dtype=bool)
    for p in base_primes:
        if p == 2:
            continue
        start = max(p * p, ((odd_lo + p - 1) // p) * p)
        if start % 2 == 0:
            start += p
        if start >= hi:
            continue
        seg[(start - odd_lo) // 2:: p] = False
    return odd_lo + 2 * np.nonzero(seg)[0]

F = np.zeros((M, K))
for i, n in enumerate(cands):
    half = n // 2
    seg_idx = 0
    lo = 3
    while lo <= half:
        hi = min(lo + SEG, half + 1)
        if seg_idx % SKIP == 0:
            ps = primes_in_segment(lo, hi)
            m = n - ps
            al = m > 1
            for qq in small:
                if n % qq:
                    al &= (m % qq != 0)
            m = m[al]
            comp = np.array([v for v in m.tolist()
                             if not is_prime_u64(int(v))], dtype=np.int64)
            for j, q in enumerate(qs_fiber):
                F[i, j] += np.count_nonzero(comp % q == 0)
        seg_idx += 1
        lo = hi
    print(f"  {i+1}/{M}  (합성 섬유합 {F[i].sum():.0f})", flush=True)
    np.savez("s_th18_fiber.npz", F=F[: i + 1],
             cands=np.array(cands[: i + 1]))

a_q = F.mean(axis=0)
chi2 = np.empty(M)
for i in range(M):
    b = F[i].sum() / a_q.sum()
    E = a_q * b
    chi2[i] = float(((F[i] - E) ** 2 / E).sum())
r = chi2 / (K - 1)
print(f"\nchi2/(K-1): 평균 {r.mean():.4f}  sigma {r.std():.4f}  "
      f"최악 {r.max():.4f}")
print(f"theta=1/3 참조: 전형 0.52~0.56, 포아송 1.0")
