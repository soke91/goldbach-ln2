"""theta=1/8 생존자 스캐너 v2 — 세그먼트 체 + p-공간 계통 샘플링.

한 n에 대해: p를 세그먼트 체로 순회(SKIP개 세그먼트당 1개 처리),
m = n-p가 y-거칠기(y = n^{1/8}) 통과 시 적응형 MR로 소수성 판정.
s(n) 추정 = 소수 비율 (표본 N_eff ~ 수백만 -> SE ~ 1e-4).

사용: python s_theta18_scan.py <X> <n개수> <SKIP>
기본: 1e10, 20, 20. 부분 저장 s_th18_<X>.npz.
"""

import sys

import numpy as np

from goldbach.sieve import primes_upto

X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10_000_000_000
NN = int(sys.argv[2]) if len(sys.argv) > 2 else 20
SKIP = int(sys.argv[3]) if len(sys.argv) > 3 else 20
SEG = 4_000_000

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

rng = np.random.default_rng(47)
base = X + (2 - X % 6) % 6
cands = sorted(set(int(base + 6 * k) for k in
                   rng.integers(0, X // 100 // 6, NN)))
y = int(round(X ** (1 / 8)))
small = [int(q) for q in primes_upto(y + 1) if q > 2]
print(f"X~{X:.2e}  n {len(cands)}개  y={y}  거칠기 홀소수 {small}  "
      f"SKIP={SKIP}", flush=True)

def primes_in_segment(lo, hi, base_primes):
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

results = []
base_primes = primes_upto(int((X // 2 + SEG) ** 0.5) + 1)
for n in cands:
    half = n // 2
    n_pr = 0
    n_srv = 0
    seg_idx = 0
    lo_seg = 3
    while lo_seg <= half:
        hi_seg = min(lo_seg + SEG, half + 1)
        if seg_idx % SKIP == 0:
            ps = primes_in_segment(lo_seg, hi_seg, base_primes)
            m = n - ps
            al = m > 1
            for q in small:
                if n % q:
                    al &= (m % q != 0)
            m = m[al]
            n_srv += len(m)
            n_pr += sum(1 for v in m.tolist() if is_prime_u64(int(v)))
        seg_idx += 1
        lo_seg = hi_seg
    s = n_pr / n_srv if n_srv else float("nan")
    results.append((n, s, n_srv))
    print(f"  n={n}  s={s:.5f}  N_eff={n_srv:,}", flush=True)
    arr = np.array(results, dtype=np.float64)
    np.savez(f"s_th18_{X}.npz", data=arr)

s_arr = np.array([r[1] for r in results])
print(f"\ns: 평균 {s_arr.mean():.5f} ± {s_arr.std():.5f}")
print(f"점근 예측 1/(8*omega(8)) = {1/(8*0.5615):.5f} (omega(8)~e^-gamma)")
