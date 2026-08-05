# -*- coding: utf-8 -*-
"""
Increment 157: does the factorization law close the SEAM-band reading?

The seam band (small k, prime-indexed, large moduli) has read
0.29-0.43 across the program's history -- the strongest apparent
sub-half-normal. By the law this must be mask accounting: N carries
7^2 (and 2), so q | k classes have forced-zero densities, and
annihilated pairs enter naive averages as zeros. Test at scale:
6000 pairs in k in [252, 464], prime index p in [P/2, P]:
  - r1_naive: |C|/sqrt(nonzero) including annihilated pairs as 0
    (the historical estimator),
  - clean bucket: gcd(k1 k2, N) = 1 -> prediction 0.798,
  - conditioning ladder by gcd(k1 k2, N),
  - annihilated / near-annihilated fraction accounting.
N = 99_999_998 = 2 x 7^2 x 1020409 (1020409 = 1013 x 1007.3? -- the
script prints the small-prime spectrum it uses).
"""
import numpy as np, time
from math import gcd

def mobius_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool); pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p*p::p] = False
            mu[p::p] *= -1
            mu[p*p::p*p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
    mu[val > 1] *= -1
    return mu

def primes_in(lo, hi):
    sieve = np.ones(hi - lo + 1, dtype=bool)
    if lo <= 1:
        sieve[:max(0, 2 - lo)] = False
    for p in range(2, int(hi ** 0.5) + 1):
        start = max(p * p, ((lo + p - 1) // p) * p)
        sieve[start - lo::p] = False
    return np.nonzero(sieve)[0] + lo

def factor_small(n, bound=1000):
    fs = {}
    m = n
    d = 2
    while d <= bound and d * d <= n:
        while m % d == 0:
            fs[d] = fs.get(d, 0) + 1
            m //= d
        d += 1
    return fs, m

def main():
    rng = np.random.default_rng(20260817)
    N = 99_999_998
    fs, cof = factor_small(N)
    print(f"N = {N} small factors {fs} cofactor {cof}", flush=True)
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)

    target = 6000
    C = np.zeros(target); sup = np.zeros(target)
    gN = np.zeros(target, dtype=np.int64)
    done = 0
    while done < target:
        k1 = int(rng.integers(252, 465)); k2 = int(rng.integers(252, 465))
        if k1 == k2:
            continue
        P1 = min(110_000, (N - 2) // max(k1, k2))
        P0 = P1 // 2
        pps = primes_in(P0, P1)
        w = N - pps.astype(np.int64) * k1
        wp = N - pps.astype(np.int64) * k2
        ok = (w > 1) & (wp > 1)
        t = mu[w[ok]].astype(np.int64) * mu[wp[ok]]
        C[done] = t.sum(); sup[done] = np.count_nonzero(t)
        gN[done] = gcd(k1 * k2, N)
        done += 1
        if done % 500 == 0:
            print(f"{done}/{target}  t={time.time()-t0:.0f}s", flush=True)

    r1_naive = np.abs(C) / np.sqrt(np.maximum(sup, 1))  # 0 when killed
    print(f"annihilated (sup=0): {np.mean(sup==0):.3f}   "
          f"near-annihilated (sup<100): {np.mean(sup<100):.3f}",
          flush=True)
    print(f"naive mean r1 (historical estimator) = "
          f"{r1_naive.mean():.3f}", flush=True)

    def rep(mask, tag):
        m = mask & (sup > 100)
        n = m.sum()
        if n < 50:
            print(f"  {tag:30s} n={n:5d}  (too few)", flush=True)
            return
        m2 = np.mean(C[m]**2 / sup[m])
        print(f"  {tag:30s} n={n:5d}  r1_eff={r1_naive[m].mean():.3f}  "
              f"m2_eff={m2:.3f}", flush=True)

    print("--- law conditioning (sup>100 only) ---", flush=True)
    rep(np.ones(target, dtype=bool), "all viable")
    rep(gN == 1, "gN = 1 (pred 0.798)")
    rep(gN == 2, "gN = 2")
    rep(gN % 7 == 0, "7 | gN")
    rep(gN % 14 == 0, "14 | gN")
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
