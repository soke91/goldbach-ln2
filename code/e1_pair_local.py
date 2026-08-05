# -*- coding: utf-8 -*-
"""
Increment 156: the pair-level local channel of the integer field.

Inc. 155: single-k D(k) is already half-normal (0.784 raw, 0.800
leave-one-out; no mean field). V3's 0.688 deficit lives in the PAIR
object T = Sum_m mu(N-k1 m) mu(N-k2 m), which has a channel the
single object lacks: for q | (k1-k2), the two arguments run CONGRUENT
mod q (N-k1 m == N-k2 m mod q for all m) -- correlated local data.
Progressive conditioning (inc-145 style):
  bucket by D_small = product of q <= 30 dividing (k1-k2),
  and by gcd(k1 k2, N) as before. Prediction if this is the whole
  story: the clean bucket (D_small = 1, gcd(k1k2, N) = 1) reads
  exactly 0.798; deficit concentrates monotonically in the
  q | (k1-k2) buckets.
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

SMALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

def dsmall(n):
    d = 1
    for q in SMALL:
        if n % q == 0:
            d *= q
    return d

def main():
    rng = np.random.default_rng(20260816)
    N = 99_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)

    target = 6000
    T = np.zeros(target); sup = np.zeros(target)
    ds = np.zeros(target, dtype=np.int64)
    gN = np.zeros(target, dtype=np.int64)
    done = 0
    while done < target:
        k1 = int(rng.integers(790, 1590)); k2 = int(rng.integers(790, 1590))
        if k1 == k2:
            continue
        M = (N - 1) // max(k1, k2)
        ms = np.arange(1, M + 1, dtype=np.int64)
        t = mu[N - k1 * ms].astype(np.int64) * mu[N - k2 * ms]
        T[done] = t.sum(); sup[done] = np.count_nonzero(t)
        ds[done] = dsmall(abs(k1 - k2))
        gN[done] = gcd(k1 * k2, N)
        done += 1
        if done % 500 == 0:
            print(f"{done}/{target}  t={time.time()-t0:.0f}s", flush=True)

    r1 = np.abs(T) / np.sqrt(np.maximum(sup, 1))

    def rep(mask, tag):
        n = mask.sum()
        if n < 50:
            return
        m2 = np.mean(T[mask]**2 / np.maximum(sup[mask], 1))
        print(f"  {tag:34s} n={n:5d}  r1_eff={r1[mask].mean():.3f}  "
              f"m2_eff={m2:.3f}", flush=True)

    print("--- all ---", flush=True)
    rep(np.ones(target, dtype=bool), "all (V3-style benchmark 0.688)")
    print("--- progressive conditioning ---", flush=True)
    clean = (ds == 1) & (gN == 1)
    rep(clean, "CLEAN: dsmall=1 & gN=1  (pred 0.798)")
    rep((ds == 1) & (gN > 1), "dsmall=1, gN>1")
    rep((ds > 1) & (gN == 1), "dsmall>1, gN=1")
    rep((ds > 1) & (gN > 1), "dsmall>1, gN>1")
    print("--- by dsmall factor ---", flush=True)
    for q in (2, 3, 5, 7):
        rep(ds % q == 0, f"q={q} | (k1-k2)")
    rep(ds == 2, "dsmall = 2 exactly")
    rep(ds == 6, "dsmall = 6 exactly")
    big = ds > 30
    rep(big, "dsmall > 30")
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
