# -*- coding: utf-8 -*-
"""
Increment 148: exact-cell decomposition of the survivor correlation.

After zero-accounting (inc. 147): the free class is EXACTLY Gaussian
(m2_eff 0.97-1.02, kurt 3.0), but shared classes keep m2_eff < 1 on
their nonzero support -- either residual cell-mixing or a genuine
negative correlation among surviving terms. Test: freeze the exact
local data -- for q in {2, 3} (the primes of N with v_q(N) = 1, 2),
classify pairs by (v_q(k), v_q(k')) capped at 2, restrict to pairs
whose OTHER N-primes are absent (gcd with 11*73*101*137 = 1), and
report z, m2_eff, kurt_eff per exact cell. If m2_eff ~ 1 in every
exact cell: all suppression is local bookkeeping (mix + zeros).
If a cell keeps m2_eff < 1: a real survivor-correlation with that
exact local signature -- the second layer of the local model.
N = 199999998 = 2 * 3^2 * 11 * 73 * 101 * 137.
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

def vq(n, q):
    v = 0
    while n % q == 0:
        n //= q; v += 1
    return min(v, 2)

def main():
    rng = np.random.default_rng(20260809)
    N = 199_999_998
    OTHER = 11 * 73 * 101 * 137
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    K0, K1 = 2000, 4000
    P0 = N // (2 * K1); P1 = 2 * P0
    ps = primes_in(P0, P1)

    target = 30000
    recs = []  # (C, nz, v2k, v2kp, v3k, v3kp)
    done = 0
    while done < target:
        k = int(rng.integers(K0, K1)); kp = int(rng.integers(K0, K1))
        if k == kp:
            continue
        if gcd(k * kp, OTHER) != 1:
            continue  # freeze the other N-primes out
        pmax = (N - 2) // max(k, kp)
        pp = ps[ps <= pmax]
        if len(pp) < 200:
            continue
        t = (mu[N - pp * k].astype(np.int64) *
             mu[N - pp * kp].astype(np.int64))
        recs.append((t.sum(), np.count_nonzero(t),
                     vq(k, 2), vq(kp, 2), vq(k, 3), vq(kp, 3)))
        done += 1
        if done % 3000 == 0:
            print(f"{done}/{target}  t={time.time()-t0:.0f}s", flush=True)

    C = np.array([r[0] for r in recs], float)
    nz = np.array([r[1] for r in recs], float)
    v2 = np.array([(min(r[2], r[3]), max(r[2], r[3])) for r in recs])
    v3 = np.array([(min(r[4], r[5]), max(r[4], r[5])) for r in recs])

    def report(mask, tag):
        n = mask.sum()
        if n < 60:
            return
        good = mask & (nz > 0)
        m2e = np.mean(C[good]**2 / nz[good])
        r1e = np.mean(np.abs(C[good]) / np.sqrt(nz[good]))
        kte = np.mean(C[good]**4 / nz[good]**2) / max(m2e, 1e-12)**2
        z = 1 - (nz[mask].sum() / (mask.sum() * np.mean(nz[mask] / (1 - 0))))  # not used
        print(f"  {tag:26s} n={n:6d}  m2_eff={m2e:.3f}  r1_eff={r1e:.3f}  "
              f"kurt_eff={kte:.2f}", flush=True)

    print("--- exact (v2 unordered pair, v3 = (0,0)) ---", flush=True)
    v3free = (v3[:, 0] == 0) & (v3[:, 1] == 0)
    for a in range(3):
        for b in range(a, 3):
            m = v3free & (v2[:, 0] == a) & (v2[:, 1] == b)
            report(m, f"v2=({a},{b}) v3-free")

    print("--- exact (v3 unordered pair, v2 = (0,0)) ---", flush=True)
    v2free = (v2[:, 0] == 0) & (v2[:, 1] == 0)
    for a in range(3):
        for b in range(a, 3):
            m = v2free & (v3[:, 0] == a) & (v3[:, 1] == b)
            report(m, f"v3=({a},{b}) v2-free")

    print("--- joint worst cells ---", flush=True)
    m = (v2[:, 0] >= 1) & (v3[:, 1] >= 1)
    report(m, "both v2min>=1 & v3max>=1")
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
