"""theta=1/8 극단 최악 후보 — 완전 y-매끈 n (모든 소인수 <= 19).

n = 2^a 3^b 5^c 7^d 11^e 13^f 17^g 19^h ~ 10^10, n = 2 mod 6 (=> b=0),
부스트 인자곱 정확히 1. 이 극단 클래스의 s가 대조군 s0와 일치하는지
(법칙 예측: 일치) + 하방 이탈 없는지.
"""

import itertools

import numpy as np

from goldbach.sieve import primes_upto

X = 10_000_000_000
SEG = 4_000_000
SKIP = 20

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

# n = 2 mod 6: 짝수, 3 나누어지지 않음 -> b=0. 소인수 {2,5,7,11,13,17,19}
prs = [2, 5, 7, 11, 13, 17, 19]
lo_t, hi_t = 0.95 * X, 1.35 * X
found = []
def gen(idx, val):
    if val > hi_t:
        return
    if idx == len(prs):
        if val >= lo_t and val % 6 == 2:
            found.append(int(val))
        return
    p = prs[idx]
    v = 1
    while val * v <= hi_t:
        gen(idx + 1, val * v)
        v *= p
gen(0, 1)
found = sorted(set(found))
print(f"완전-매끈 후보 {len(found)}개: {found[:12]}", flush=True)
targets = found[:10]

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

rows = []
for n in targets:
    half = n // 2
    n_pr = n_srv = 0
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
            n_srv += len(m)
            n_pr += sum(1 for v in m.tolist() if is_prime_u64(int(v)))
        seg_idx += 1
        lo = hi
    s = n_pr / n_srv if n_srv else float("nan")
    rows.append((n, s, n_srv))
    print(f"  n={n}  s={s:.5f}  N_eff={n_srv:,}", flush=True)
    np.savez("s_th18_pure_smooth.npz", data=np.array(rows, dtype=np.float64))

s_arr = np.array([r[1] for r in rows])
print(f"\n완전-매끈 s: 평균 {s_arr.mean():.5f} ± {s_arr.std():.5f}  "
      f"min {s_arr.min():.5f}")
print(f"대조군(준-매끈) 참조: 0.25444 ± 0.00019, 40개 min 0.25388")
