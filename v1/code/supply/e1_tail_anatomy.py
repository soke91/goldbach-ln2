# -*- coding: utf-8 -*-
"""
Increment 145: anatomy of the heavy tail of C_{k,k'}.

Increment 144 found kurt(C) = 5.3-5.8 (vs Gaussian 3) while the
variance is suppressed (m2 = 0.36-0.40): C is a MIXTURE. Candidate
drivers of the tail, tested by conditioning:
  (a) g = gcd(k, k')      -- shared arithmetic between the two dilates
  (b) gcd-with-N           -- k or k' sharing factors with N
  (c) near-ratio           -- k'/k close to a rational with tiny height
  (d) smoothness of k*k'   -- structure-law analog (many small factors)
For each bucket: n, mean |C|/sqrt(n_p), variance ratio, kurtosis.
If one conditioning flattens kurtosis to ~3 within buckets, the tail
is that arithmetic -- the mixture is then explainable (and the
fourth-moment route can price it).
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

def smooth_part(n, B=100):
    s = 1
    m = n
    for p in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
              73,79,83,89,97):
        while m % p == 0:
            m //= p; s *= p
    return s

def stats(tag, vals, nps):
    vals = np.asarray(vals, float); nps = np.asarray(nps, float)
    if len(vals) < 30:
        print(f"  {tag:28s} n={len(vals):5d}  (too few)", flush=True)
        return
    r1 = np.mean(np.abs(vals) / np.sqrt(nps))
    m2 = np.mean(vals**2 / nps)
    kurt = np.mean(vals**4) / np.mean(vals**2)**2
    print(f"  {tag:28s} n={len(vals):5d}  r1={r1:.3f}  m2={m2:.3f}  "
          f"kurt={kurt:.2f}", flush=True)

def main():
    rng = np.random.default_rng(20260806)
    N = 200_000_000 - 2  # even
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    K0, K1 = 2000, 4000
    P0 = N // (2 * K1); P1 = 2 * P0
    ps = primes_in(P0, P1)

    rows = []  # (C, n_p, g, gN, smooth(kk'), minheight)
    target = 6000
    done = 0
    while done < target:
        k = int(rng.integers(K0, K1)); kp = int(rng.integers(K0, K1))
        if k == kp:
            continue
        pmax = (N - 2) // max(k, kp)
        pp = ps[ps <= pmax]
        if len(pp) < 200:
            continue
        C = int(np.sum(mu[N - pp * k].astype(np.int64) *
                       mu[N - pp * kp].astype(np.int64)))
        g = gcd(k, kp)
        gN = gcd(k * kp, N)
        sm = smooth_part(k // g) * smooth_part(kp // g)
        # height of k'/k in lowest terms
        h = max(k // g, kp // g)
        rows.append((C, len(pp), g, gN, sm, h))
        done += 1
        if done % 500 == 0:
            print(f"{done}/{target}  t={time.time()-t0:.0f}s", flush=True)

    C = np.array([r[0] for r in rows], float)
    np_ = np.array([r[1] for r in rows], float)
    g = np.array([r[2] for r in rows])
    gN = np.array([r[3] for r in rows])
    sm = np.array([r[4] for r in rows], float)
    h = np.array([r[5] for r in rows], float)

    print("=== ALL ===", flush=True)
    stats("all", C, np_)

    print("=== (a) by gcd(k,k') ===", flush=True)
    stats("g = 1", C[g == 1], np_[g == 1])
    stats("g in [2,10]", C[(g > 1) & (g <= 10)], np_[(g > 1) & (g <= 10)])
    stats("g > 10", C[g > 10], np_[g > 10])

    print("=== (b) by gcd(kk', N) ===", flush=True)
    stats("gN = 1", C[gN == 1], np_[gN == 1])
    stats("gN in [2,10]", C[(gN > 1) & (gN <= 10)], np_[(gN > 1) & (gN <= 10)])
    stats("gN > 10", C[gN > 10], np_[gN > 10])

    print("=== (c) by reduced height max(k/g, k'/g) ===", flush=True)
    med = np.median(h)
    stats(f"height <= {med:.0f}", C[h <= med], np_[h <= med])
    stats(f"height > {med:.0f}", C[h > med], np_[h > med])
    q10 = np.quantile(h, 0.1)
    stats(f"height <= {q10:.0f} (10%)", C[h <= q10], np_[h <= q10])

    print("=== (d) by 100-smooth part of (k/g)(k'/g) ===", flush=True)
    meds = np.median(sm)
    stats(f"smooth <= {meds:.0f}", C[sm <= meds], np_[sm <= meds])
    stats(f"smooth > {meds:.0f}", C[sm > meds], np_[sm > meds])
    q90 = np.quantile(sm, 0.9)
    stats(f"smooth > {q90:.0f} (top10%)", C[sm > q90], np_[sm > q90])

    # tail membership: which conditioning over-represents the tail?
    print("=== tail (|C|/sqrt(n_p) > 1.0) composition ===", flush=True)
    tail = (np.abs(C) / np.sqrt(np_)) > 1.0
    print(f"  tail fraction: {tail.mean():.4f}", flush=True)
    for name, mask in (("g>1", g > 1), ("gN>1", gN > 1),
                       ("height<=q10", h <= q10), ("smooth>q90", sm > q90)):
        base = mask.mean()
        intail = mask[tail].mean() if tail.sum() else float('nan')
        print(f"  {name:12s} base={base:.3f}  in-tail={intail:.3f}  "
              f"lift={intail/base if base else float('nan'):.2f}", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
