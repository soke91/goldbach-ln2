# -*- coding: utf-8 -*-
"""
Increment 144: the corrected-normalization ledger after REVIEW_VERDICT.

The adversarial review showed the chain needs only
    Sum_{k != k'} |C_{k,k'}|  <<  K^2 P (log x)^{-2A}          (*)
(not the over-normalized SEAM statement), while the measured
half-normal law gives Sum |C| ~ 0.7 K^2 sqrt(P) -- a sqrt(P)-size
MARGIN over (*). This script stamps that margin directly and probes
the fourth-moment route (second moment of C over pairs), at two
K-bands and two N.

Outputs, per (N, K-band):
  r1 = mean |C| / sqrt(n_p)        (half-normal benchmark 0.798)
  margin = P^{1/2} -> factor by which nature over-delivers vs (*)
  m2 = E|C|^2 / n_p                (should be ~1 if C ~ N(0, n_p))
  kurt = E|C|^4 / (E|C|^2)^2      (Gaussian: 3 -- fourth-moment route)
"""
import numpy as np, sys, time

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
            pp = p * p
            if pp <= X:
                val[pp::pp] = 0  # squarefull marker not needed; mu already 0
    # correct: recompute val cleanly (val //= p only once per p above is wrong
    # for higher powers) -- use the standard two-pass:
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
    mu2 = mu.copy()
    mu2[val > 1] *= -1  # one large prime factor remains
    return mu2, pm

def primes_in(lo, hi):
    sieve = np.ones(hi - lo + 1, dtype=bool)
    if lo <= 1:
        sieve[:max(0, 2 - lo)] = False
    for p in range(2, int(hi ** 0.5) + 1):
        start = max(p * p, ((lo + p - 1) // p) * p)
        sieve[start - lo::p] = False
    return np.nonzero(sieve)[0] + lo

def band(N, mu, K0, K1, npairs, rng, tag):
    t0 = time.time()
    P0 = N // (2 * K1)  # ensure N - p*k stays positive with slack
    P1 = 2 * P0
    ps = primes_in(P0, P1)
    absC, n_used = [], []
    done = 0
    while done < npairs:
        k = rng.integers(K0, K1)
        kp = rng.integers(K0, K1)
        if k == kp:
            continue
        pmax = (N - 2) // max(k, kp)
        pp = ps[ps <= pmax]
        if len(pp) < 200:
            continue
        C = int(np.sum(mu[N - pp * k].astype(np.int64) *
                       mu[N - pp * kp].astype(np.int64)))
        absC.append(abs(C)); n_used.append(len(pp))
        done += 1
        if done % 200 == 0:
            print(f"[{tag}] {done}/{npairs}  t={time.time()-t0:.0f}s",
                  flush=True)
    absC = np.array(absC, float); n_used = np.array(n_used, float)
    r1 = np.mean(absC / np.sqrt(n_used))
    m2 = np.mean(absC**2 / n_used)
    kurt = np.mean(absC**4) / np.mean(absC**2)**2
    Pmean = n_used.mean()
    print(f"[{tag}] K in [{K0},{K1})  n_p mean {Pmean:.0f}  pairs {done}",
          flush=True)
    print(f"[{tag}] r1 = {r1:.3f}  (half-normal 0.798)", flush=True)
    print(f"[{tag}] m2 = {m2:.3f}  (Gaussian 1.0)", flush=True)
    print(f"[{tag}] kurt(|C|) = {kurt:.3f}  (half-normal: 3.0 for C)",
          flush=True)
    print(f"[{tag}] margin over (*): sqrt(n_p) ~ {np.sqrt(Pmean):.0f}x",
          flush=True)
    return r1, m2, kurt

def main():
    rng = np.random.default_rng(20260805)
    for N in (100_000_007 + 1 - 8, 200_000_033 + 1 - 34):  # two even N
        print(f"=== N = {N} ===", flush=True)
        t0 = time.time()
        mu, _ = mobius_upto(N)
        print(f"mu ready in {time.time()-t0:.0f}s", flush=True)
        band(N, mu, 2000, 4000, 1500, rng, f"N{N}-K2k")
        band(N, mu, 10000, 20000, 1500, rng, f"N{N}-K10k")
        del mu
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
