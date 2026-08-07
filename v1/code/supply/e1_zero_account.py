# -*- coding: utf-8 -*-
"""
Increment 147: the zero-accounting test (decisive for the local model).

Hypothesis: ALL of the variance suppression in C_{k,k'} is forced
mu = 0 terms (shared factors with N force divisibility of N - pk,
killing squarefreeness at computable density) -- i.e., on its nonzero
support the field is exactly Gaussian. Test: per pair, count the zero
terms directly, set z = zero-density, and check
    m2_eff = E[C^2] / (n_p (1 - z))  ~  0.85-1.0   universally,
    r1_eff = |C| / sqrt(n_p (1 - z)) ~  0.798      universally.
If yes at every class: C = Gaussian x computable local density --
the dispersion field has NO arithmetic mystery beyond forced zeros.
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

def main():
    rng = np.random.default_rng(20260808)
    for N in (199_999_998, 100_000_000):  # rich spectrum / 2^8*5^8-type
        print(f"=== N = {N} ===", flush=True)
        t0 = time.time()
        mu = mobius_upto(N)
        print(f"mu ready {time.time()-t0:.0f}s", flush=True)
        K0, K1 = 2000, 4000
        P0 = N // (2 * K1); P1 = 2 * P0
        ps = primes_in(P0, P1)

        target = 12000
        C = np.zeros(target); npz = np.zeros(target)  # nonzero-support size
        np_ = np.zeros(target); gN = np.zeros(target, dtype=np.int64)
        done = 0
        while done < target:
            k = int(rng.integers(K0, K1)); kp = int(rng.integers(K0, K1))
            if k == kp:
                continue
            pmax = (N - 2) // max(k, kp)
            pp = ps[ps <= pmax]
            if len(pp) < 200:
                continue
            a = mu[N - pp * k].astype(np.int64)
            b = mu[N - pp * kp].astype(np.int64)
            t = a * b
            C[done] = t.sum()
            np_[done] = len(pp)
            npz[done] = np.count_nonzero(t)
            gN[done] = gcd(k * kp, N)
            done += 1
            if done % 2000 == 0:
                print(f"{done}/{target}  t={time.time()-t0:.0f}s",
                      flush=True)

        z = 1 - npz / np_
        m2_raw = C**2 / np_
        m2_eff = np.where(npz > 0, C**2 / np.maximum(npz, 1), 0.0)
        r1_eff = np.abs(C) / np.sqrt(np.maximum(npz, 1))

        def report(mask, tag):
            n = mask.sum()
            if n < 60:
                return
            print(f"  {tag:22s} n={n:6d}  z={z[mask].mean():.3f}  "
                  f"m2_raw={m2_raw[mask].mean():.3f}  "
                  f"m2_eff={m2_eff[mask].mean():.3f}  "
                  f"r1_eff={r1_eff[mask].mean():.3f}  "
                  f"kurt_eff={np.mean(C[mask]**4/np.maximum(npz[mask],1)**2)/max(np.mean(C[mask]**2/np.maximum(npz[mask],1)),1e-12)**2:.2f}",
                  flush=True)

        print("--- zero-accounting (target: m2_eff ~ 0.85-1.0, "
              "r1_eff ~ 0.74-0.80 everywhere) ---", flush=True)
        report(np.ones(target, dtype=bool), "all")
        report(gN == 1, "gN=1 (free)")
        report((gN > 1) & (gN <= 10), "gN 2-10")
        report((gN > 10) & (gN <= 100), "gN 11-100")
        report(gN > 100, "gN>100")
        hi = z > np.quantile(z, 0.9)
        report(hi, "z top-10%")
        # correlation: does z explain m2_raw?
        c = np.corrcoef(z, m2_raw)[0, 1]
        print(f"  corr(z, m2_raw) = {c:.3f}", flush=True)
        del mu
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
