# -*- coding: utf-8 -*-
"""
Increment 160: the combined blind -- mask + Gaussian predicted BEFORE
any mu value is seen, at 1e9, at a fresh N with different 2-adic
structure (v2(N) = 2).

Phase 1 (no mu access): for 3000 sampled pairs, compute s_pred by the
mask algorithm alone (enumeration over units mod q^2, q <= 50 + tail),
and PRINT the pre-registered predictions:
  P1 = predicted annihilation fraction  (s_pred < 0.01),
  P2 = predicted naive mean r1          (0.798 * mean sqrt(s_pred)),
  P3 = predicted viable r1_eff          (0.798),
  P4 = predicted viable m2_eff          (1.00).
Phase 2: build segmented mu and measure the same quantities.
"""
import numpy as np, time
from math import gcd

QCUT = 50

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

def rho_q(q, N, k, kp):
    q2 = q * q
    cnt = tot = 0
    for p in range(q2):
        if p % q == 0:
            continue
        tot += 1
        if (N - p * k) % q2 == 0 or (N - p * kp) % q2 == 0:
            cnt += 1
    return cnt / tot

def main():
    rng = np.random.default_rng(20260820)
    N = 999_999_996  # v2 = 2, 3 | N
    K0, K1 = 3000, 4000
    SPS = [int(q) for q in primes_upto(QCUT)]
    ps_tail = primes_upto(100000)
    tf = 1.0
    for q in ps_tail:
        if q > QCUT:
            tf *= (1 - 2 / (int(q) * int(q)))
    print(f"N = {N}  tail factor = {tf:.5f}", flush=True)

    # sample pairs first (fixed list, used in both phases)
    pairs = []
    while len(pairs) < 3000:
        k = int(rng.integers(K0, K1)); kp = int(rng.integers(K0, K1))
        if k != kp:
            pairs.append((k, kp))

    t0 = time.time()
    s_pred = np.zeros(len(pairs))
    for i, (k, kp) in enumerate(pairs):
        pred = tf
        for q in SPS:
            pred *= (1 - rho_q(q, N, k, kp))
        s_pred[i] = pred
        if i % 500 == 499:
            print(f"mask {i+1}/3000  t={time.time()-t0:.0f}s", flush=True)

    ann_pred = float(np.mean(s_pred < 0.01))
    naive_pred = float(0.798 * np.mean(np.sqrt(s_pred)))
    print("=== PRE-REGISTERED (no mu accessed) ===", flush=True)
    print(f"  P1 annihilation fraction = {ann_pred:.3f}", flush=True)
    print(f"  P2 naive mean r1        = {naive_pred:.3f}", flush=True)
    print(f"  P3 viable r1_eff        = 0.798", flush=True)
    print(f"  P4 viable m2_eff        = 1.00", flush=True)

    # phase 2: measure
    P0 = N // (2 * K1); P1 = int(1.6 * P0)
    A = max(2, N - P1 * (K1 - 1) - 10); B = N - P0 * K0 + 10
    bp = primes_upto(int(B ** 0.5) + 1)
    print(f"building mu window [{A},{B}) ({(B-A)/1e6:.0f}M)...",
          flush=True)
    mu = mobius_range(A, B, bp)
    print(f"mu ready t={time.time()-t0:.0f}s", flush=True)
    ps = primes_in(P0, P1)

    C = np.zeros(len(pairs)); sup = np.zeros(len(pairs))
    nterm = np.zeros(len(pairs))
    for i, (k, kp) in enumerate(pairs):
        pmax = (N - A - 10) // max(k, kp)
        pp = ps[ps <= pmax]
        t = mu[N - pp * k - A].astype(np.int64) * mu[N - pp * kp - A]
        C[i] = t.sum(); sup[i] = np.count_nonzero(t); nterm[i] = len(pp)
        if i % 500 == 499:
            print(f"measure {i+1}/3000  t={time.time()-t0:.0f}s",
                  flush=True)

    ann_obs = float(np.mean(sup / nterm < 0.01))
    r1_naive = float(np.mean(np.abs(C) / np.sqrt(np.maximum(sup, 1))))
    viable = sup > 100
    r1_eff = float(np.mean(np.abs(C[viable]) / np.sqrt(sup[viable])))
    m2_eff = float(np.mean(C[viable]**2 / sup[viable]))
    print("=== OBSERVED ===", flush=True)
    print(f"  P1 annihilation fraction = {ann_obs:.3f}  "
          f"(pred {ann_pred:.3f})", flush=True)
    print(f"  P2 naive mean r1        = {r1_naive:.3f}  "
          f"(pred {naive_pred:.3f})", flush=True)
    print(f"  P3 viable r1_eff        = {r1_eff:.3f}  (pred 0.798)",
          flush=True)
    print(f"  P4 viable m2_eff        = {m2_eff:.3f}  (pred 1.00)",
          flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
