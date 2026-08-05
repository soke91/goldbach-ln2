# -*- coding: utf-8 -*-
"""
Increment 161: which cells break at v2(N) = 2?

Inc. 160 found viable-class m2_eff = 0.917 (-2.4 sigma) at
N = 999_999_996 (v2 = 2, v3 = 1). Since a mixture of unit cells
averages to 1, some exact cell must be genuinely below unit on its
support. Decompose: (v2(k), v2(k')) x (v3(k), v3(k')) cells, freeze
nothing else (record gcd with the small factorization of N beyond
2,3 -- here N = 2^2 * 3 * 83333333, so only 2 and 3 matter below
1000). Per cell: n, mean support fraction, m2_eff, kurt_eff.
"""
import numpy as np, time
from math import gcd

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i*i::i] = False
    return np.nonzero(s)[0]

def mobius_range(A, B, base_primes, seg=10_000_000):
    out = np.empty(B - A, dtype=np.int8)
    for lo in range(A, B, seg):
        hi = min(lo + seg, B)
        mu = np.ones(hi - lo, dtype=np.int8)
        val = np.arange(lo, hi, dtype=np.int64)
        for p in base_primes:
            p = int(p)
            start = ((lo + p - 1) // p) * p - lo
            mu[start::p] *= -1
            val[start::p] //= p
            pp = p * p
            if pp < hi:
                start2 = ((lo + pp - 1) // pp) * pp - lo
                mu[start2::pp] = 0
        mu[val > 1] *= -1
        out[lo - A:hi - A] = mu
    return out

def primes_in(lo, hi):
    sieve = np.ones(hi - lo + 1, dtype=bool)
    if lo <= 1:
        sieve[:max(0, 2 - lo)] = False
    for p in range(2, int(hi ** 0.5) + 1):
        start = max(p * p, ((lo + p - 1) // p) * p)
        sieve[start - lo::p] = False
    return np.nonzero(sieve)[0] + lo

def vq(n, q, cap=3):
    v = 0
    while n % q == 0 and v < cap:
        n //= q; v += 1
    return v

def main():
    rng = np.random.default_rng(20260821)
    N = 999_999_996
    K0, K1 = 3000, 4000
    P0 = N // (2 * K1); P1 = int(1.6 * P0)
    A = max(2, N - P1 * (K1 - 1) - 10); B = N - P0 * K0 + 10
    t0 = time.time()
    bp = primes_upto(int(B ** 0.5) + 1)
    print(f"building mu [{A},{B})...", flush=True)
    mu = mobius_range(A, B, bp)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    ps = primes_in(P0, P1)

    target = 12000
    recs = []
    done = 0
    while done < target:
        k = int(rng.integers(K0, K1)); kp = int(rng.integers(K0, K1))
        if k == kp:
            continue
        pmax = (N - A - 10) // max(k, kp)
        pp = ps[ps <= pmax]
        t = mu[N - pp * k - A].astype(np.int64) * mu[N - pp * kp - A]
        sup = int(np.count_nonzero(t))
        recs.append((float(t.sum()), sup, len(pp),
                     vq(k, 2), vq(kp, 2), vq(k, 3), vq(kp, 3)))
        done += 1
        if done % 2000 == 0:
            print(f"{done}/{target}  t={time.time()-t0:.0f}s", flush=True)

    C = np.array([r[0] for r in recs])
    sup = np.array([r[1] for r in recs], float)
    nt = np.array([r[2] for r in recs], float)
    v2 = np.array([tuple(sorted((r[3], r[4]))) for r in recs])
    v3 = np.array([tuple(sorted((r[5], r[6]))) for r in recs])

    def rep(mask, tag):
        m = mask & (sup > 100)
        n = m.sum()
        if n < 60:
            if mask.sum() >= 60:
                print(f"  {tag:22s} n={mask.sum():5d}  "
                      f"supfrac={np.mean(sup[mask]/nt[mask]):.3f}  "
                      f"(annihilated)", flush=True)
            return
        m2 = np.mean(C[m]**2 / sup[m])
        kt = np.mean(C[m]**4 / sup[m]**2) / max(m2, 1e-12)**2
        print(f"  {tag:22s} n={n:5d}  supfrac="
              f"{np.mean(sup[m]/nt[m]):.3f}  m2_eff={m2:.3f}  "
              f"kurt={kt:.2f}", flush=True)

    print("--- v2 cells (v3 = (0,0) frozen) ---", flush=True)
    v3free = (v3[:, 0] == 0) & (v3[:, 1] == 0)
    for a in range(4):
        for b in range(a, 4):
            m = v3free & (v2[:, 0] == a) & (v2[:, 1] == b)
            if m.sum() >= 60:
                rep(m, f"v2=({a},{b})")

    print("--- v3 cells (v2 = (0,0) frozen) ---", flush=True)
    v2free = (v2[:, 0] == 0) & (v2[:, 1] == 0)
    for a in range(4):
        for b in range(a, 4):
            m = v2free & (v3[:, 0] == a) & (v3[:, 1] == b)
            if m.sum() >= 60:
                rep(m, f"v3=({a},{b})")

    print("--- all viable ---", flush=True)
    rep(np.ones(len(recs), dtype=bool), "all")
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
