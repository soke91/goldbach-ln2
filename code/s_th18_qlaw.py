"""theta=1/8 잔류 구조 법칙 — 통제된 중간-소수 클래스 실험.

클래스: q in {19, 23, 29, 37, 53, 101, 211, 499, 997} + 대조군(중간인수 없음).
각 클래스에서 n ~ 10^10, n = 2 mod 6, q | n, 다른 (18, 2000] 인수 없음.
클래스당 4개 n, s 측정 (SKIP=20) -> s(q) 프로파일.

이론 후보: 소박 곱 = s0 * (q-1)/(q-2) [트윈 부스트]
           + 죽은-류 보정 = /(1 - P(q|m, m rough)) 의 결합 형태 판정.
"""

import sys

import numpy as np

from goldbach.sieve import primes_upto

X = 10_000_000_000
SEG = 4_000_000
SKIP = 20
PER = 4
QLIST = [19, 23, 29, 37, 53, 101, 211, 499, 997, None]

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

med_primes = [int(p) for p in primes_upto(2000) if p > 18]

def no_other_med(n, q):
    for p in med_primes:
        if p != q and n % p == 0:
            return False
    return True

def pick_ns(q, count, seed):
    rng = np.random.default_rng(seed)
    out = []
    if q is None:
        base = X + (2 - X % 6) % 6
        while len(out) < count:
            n = int(base + 6 * int(rng.integers(0, X // 200 // 6)))
            if no_other_med(n, -1):
                out.append(n)
        return out
    step = 6 * q
    n0 = None
    for k in range(6 * q):
        cand = X + k
        if cand % 6 == 2 and cand % q == 0:
            n0 = cand
            break
    while len(out) < count:
        n = int(n0 + step * int(rng.integers(0, X // 300 // step)))
        if (n // q) % q != 0 and no_other_med(n, q):
            out.append(n)
    return out

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

def measure(n):
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
    return (n_pr / n_srv if n_srv else float("nan")), n_srv

all_rows = []
for qi, q in enumerate(QLIST):
    ns = pick_ns(q, PER, 53 + qi)
    svals = []
    for n in ns:
        s, N = measure(n)
        svals.append(s)
        all_rows.append((q if q else 0, n, s, N))
        np.savez("s_th18_qlaw.npz", data=np.array(all_rows, dtype=np.float64))
    sv = np.array(svals)
    tag = f"q={q}" if q else "대조군"
    boost = ""
    if q:
        boost = f"  트윈부스트 예측비 {(q-1)/(q-2):.4f}"
    print(f"{tag:>8}: s = {sv.mean():.5f} ± {sv.std():.5f}{boost}", flush=True)

d = np.array(all_rows, dtype=np.float64)
ctrl = d[d[:, 0] == 0][:, 2].mean()
print(f"\n대조군 기준 s0 = {ctrl:.5f}")
print(f"{'q':>5} {'s비':>8} {'트윈예측':>8} {'죽은류만':>8}")
for q in [v for v in QLIST if v]:
    sq = d[d[:, 0] == q][:, 2].mean()
    print(f"{q:>5} {sq/ctrl:>8.4f} {(q-1)/(q-2):>8.4f} "
          f"{1/(1-2.6/q):>8.4f}")
