# -*- coding: utf-8 -*-
"""
Increment 158: the factorization law at 10^9 -- blind scale stamp.

Pre-registered predictions (before running):
  free-class (gcd(kk', N) = 1) prime-indexed pairs at N ~ 1e9:
    m2_eff = 1.00 +- 0.04,  kurtosis = 3.0 +- 0.3,
    r1_eff = 0.798 +- 0.02.
Segmented Mobius over the needed window keeps memory ~1 GB.
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
    """mu on [A, B) as int8 array of length B-A (segmented)."""
    out = np.empty(B - A, dtype=np.int8)
    for lo in range(A, B, seg):
        hi = min(lo + seg, B)
        n = hi - lo
        mu = np.ones(n, dtype=np.int8)
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

def main():
    rng = np.random.default_rng(20260819)
    N = 999_999_998
    K0, K1 = 3000, 4000
    P0 = N // (2 * K1)          # ~125000
    P1 = int(1.6 * P0)
    lo_w = N - P1 * (K1 - 1)    # smallest w needed
    hi_w = N - P0 * K0          # largest w
    A, B = max(2, lo_w - 10), hi_w + 10
    print(f"N={N}  mu window [{A}, {B})  size {(B-A)/1e6:.0f}M",
          flush=True)
    t0 = time.time()
    bp = primes_upto(int(B ** 0.5) + 1)
    mu = mobius_range(A, B, bp)
    print(f"segmented mu ready {time.time()-t0:.0f}s", flush=True)
    ps = primes_in(P0, P1)
    print(f"index primes: {len(ps)}", flush=True)

    target = 6000
    m2s, r1s, Cs = [], [], []
    done = 0
    while done < target:
        k = int(rng.integers(K0, K1)); kp = int(rng.integers(K0, K1))
        if k == kp or gcd(k * kp, N) != 1:
            continue
        pmax = (N - A - 10) // max(k, kp)
        pp = ps[ps <= pmax]
        if len(pp) < 300:
            continue
        w1 = N - pp * k - A
        w2 = N - pp * kp - A
        t = mu[w1].astype(np.int64) * mu[w2]
        nz = int(np.count_nonzero(t))
        if nz < 100:
            continue
        c = float(t.sum())
        m2s.append(c * c / nz)
        r1s.append(abs(c) / np.sqrt(nz))
        Cs.append(c / np.sqrt(nz))
        done += 1
        if done % 500 == 0:
            print(f"{done}/{target}  t={time.time()-t0:.0f}s", flush=True)

    m2 = float(np.mean(m2s))
    r1 = float(np.mean(r1s))
    Cs = np.array(Cs)
    kurt = float(np.mean(Cs**4) / np.mean(Cs**2)**2)
    print(f"RESULT (free class, N ~ 1e9, {done} pairs):", flush=True)
    print(f"  m2_eff = {m2:.3f}   (pred 1.00 +- 0.04)", flush=True)
    print(f"  r1_eff = {r1:.3f}   (pred 0.798 +- 0.02)", flush=True)
    print(f"  kurt   = {kurt:.2f}    (pred 3.0 +- 0.3)", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
