# -*- coding: utf-8 -*-
"""
Increment 146: the local model of the dispersion suppression.

Increment 145 decomposed the mixture: pairs with gcd(kk', N) = 1 are
EXACTLY Gaussian (r1 0.740, m2 0.872, kurt 3.12); shared factors with
N suppress the variance. Hypothesis (structure-law analog): the
variance ratio is multiplicative over the primes q | N dividing kk',
    m2(class) = m2_free * prod_{q | N, q | kk'} f(q),
with f(q) < 1 a local factor. This script extracts f(q) per prime,
tests multiplicativity, checks within-class Gaussianity, and compares
f(q) against candidate local models 1/(q-1)^2-type / (q-2)/(q-1)-type.

N = 199999998 = 2 * 3^2 * 11 * 73 * 101 * 137  (rich small spectrum)
"""
import numpy as np, time
from math import gcd
from itertools import combinations

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

def factor_small(n):
    fs = []
    for p in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
              73,79,83,89,97,101,103,107,109,113,127,131,137,139,149):
        if n % p == 0:
            fs.append(p)
            while n % p == 0:
                n //= p
    return fs, n

def main():
    rng = np.random.default_rng(20260807)
    N = 199_999_998
    fs, rest = factor_small(N)
    print(f"N = {N} small primes {fs} cofactor {rest}", flush=True)
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    K0, K1 = 2000, 4000
    P0 = N // (2 * K1); P1 = 2 * P0
    ps = primes_in(P0, P1)

    target = 20000
    C = np.zeros(target); np_ = np.zeros(target)
    sig = []  # frozenset of N-primes dividing k*k'
    done = 0
    while done < target:
        k = int(rng.integers(K0, K1)); kp = int(rng.integers(K0, K1))
        if k == kp:
            continue
        pmax = (N - 2) // max(k, kp)
        pp = ps[ps <= pmax]
        if len(pp) < 200:
            continue
        c = int(np.sum(mu[N - pp * k].astype(np.int64) *
                       mu[N - pp * kp].astype(np.int64)))
        C[done] = c; np_[done] = len(pp)
        kk = k * kp
        sig.append(frozenset(q for q in fs if kk % q == 0))
        done += 1
        if done % 2000 == 0:
            print(f"{done}/{target}  t={time.time()-t0:.0f}s", flush=True)

    sig = np.array([hash(s) for s in sig]), sig
    hashes, sets = sig

    def cell(mask, tag):
        n = mask.sum()
        if n < 80:
            return None
        m2 = np.mean(C[mask]**2 / np_[mask])
        kurt = np.mean(C[mask]**4) / np.mean(C[mask]**2)**2
        print(f"  {tag:24s} n={n:6d}  m2={m2:.3f}  kurt={kurt:.2f}",
              flush=True)
        return m2

    empty = np.array([len(s) == 0 for s in sets])
    print("=== free class (gcd(kk',N)=1) ===", flush=True)
    m2_free = cell(empty, "free")

    print("=== single-prime marginal f(q) = m2(q|kk')/m2_free ===",
          flush=True)
    fq = {}
    for q in fs:
        mask = np.array([(q in s) for s in sets])
        m2q_only = None
        # marginal: q | kk' regardless of others -- and exact: {q} only
        only = np.array([s == frozenset([q]) for s in sets])
        m2o = cell(only, f"exactly {{{q}}}")
        if m2o and m2_free:
            fq[q] = m2o / m2_free
            print(f"    f({q}) = {fq[q]:.3f}   candidates: "
                  f"(q-2)/(q-1)={max(q-2,0)/(q-1):.3f}  "
                  f"1-2/(q-1)^2={1-2/(q-1)**2 if q>2 else float('nan'):.3f}  "
                  f"((q-1)/(q-2))^-2={((q-2)/(q-1))**2 if q>2 else float('nan'):.3f}",
                  flush=True)

    print("=== multiplicativity: 2-prime classes vs product ===",
          flush=True)
    for qa, qb in combinations([q for q in fs if q in fq], 2):
        pair = np.array([s == frozenset([qa, qb]) for s in sets])
        if pair.sum() < 80:
            continue
        m2p = np.mean(C[pair]**2 / np_[pair])
        pred = m2_free * fq[qa] * fq[qb]
        print(f"  {{{qa},{qb}}} n={pair.sum():5d}  m2={m2p:.3f}  "
              f"pred={pred:.3f}  ratio={m2p/pred:.3f}", flush=True)

    print("=== within-class Gaussianity (kurt ~ 3?) summary above ===",
          flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
