"""전선1 수치 — theta=1/8 합성 생존자의 spf-분포 실측 (지도 가중 검증).

n 8개 @ 10^10, SKIP=40. 합성 생존자를 x^{1/3}=2154까지 시행나눗셈 →
spf 히스토그램 (alpha = ln spf/ln n 구간별) vs Buchstab 지도 가중.
spf > x^{1/3} 은 P2 코어 (보조인자 자동 소수).
"""

import numpy as np

from goldbach.sieve import primes_upto

X = 10_000_000_000
SEG = 4_000_000
SKIP = 40
M = 8

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

rng = np.random.default_rng(71)
base = X + (2 - X % 6) % 6
cands = sorted(int(base + 6 * k) for k in
               rng.integers(0, X // 200 // 6, M))

y = int(round(X ** (1 / 8)))
small = [int(p) for p in primes_upto(y + 1) if p > 2]
x13 = int(round(X ** (1 / 3)))
med = np.array([int(p) for p in primes_upto(x13 + 10) if p > 19])
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

lnX = np.log(X)
edges = [0.125, 0.15, 1/6, 0.2, 0.25, 0.30, 1/3]
hist = np.zeros(len(edges) + 1)   # 마지막 = 코어 (spf > x^{1/3})
tot_comp = 0
tot_all = 0
for idx, n in enumerate(cands):
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
            tot_all += len(m)
            comp = np.array([v for v in m.tolist()
                             if not is_prime_u64(int(v))], dtype=np.int64)
            tot_comp += len(comp)
            # spf via 벡터화 시행나눗셈
            spf = np.full(len(comp), -1, dtype=np.int64)
            rem = comp.copy()
            todo = np.ones(len(comp), dtype=bool)
            for q in med:
                if not todo.any():
                    break
                hitq = todo & (rem % q == 0)
                spf[hitq] = q
                todo &= ~hitq
            alpha = np.where(spf > 0, np.log(np.maximum(spf, 2)) / np.log(n),
                             1.0)
            for j, e in enumerate(edges):
                lo_e = edges[j - 1] if j else 0
                hist[j] += np.count_nonzero((alpha >= (edges[j-1] if j else 0))
                                            & (alpha < e))
            hist[-1] += np.count_nonzero(spf < 0)
        seg_idx += 1
        lo = hi
    print(f"  {idx+1}/{M}", flush=True)
    np.savez("s_th18_spf.npz", hist=hist, tot=[tot_comp, tot_all])

print(f"\n합성 {tot_comp:,} / 전체 {tot_all:,} (합성 비율 "
      f"{tot_comp/tot_all:.4f}; Buchstab 예상 ~0.745@y19)")
labels = []
prev = 0.0
for e in edges:
    labels.append(f"{prev:.3f}-{e:.3f}")
    prev = e
labels.append(">1/3 (P2코어)")
print(f"{'alpha 구간':>14} {'질량(실측)':>10}")
for lb, h in zip(labels, hist):
    print(f"{lb:>14} {h/tot_comp:>10.4f}")
print("\n지도 가중 참조 (Buchstab, u=8): 중간 0.65 / 코어 0.35")
