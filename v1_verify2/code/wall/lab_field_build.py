# -*- coding: utf-8 -*-
"""
Build C(N), V(N), A(N), S(N) on every even N <= X, from scratch.

This is infrastructure for the v1_verify2 Phase 1 measurement arm. It is
written from the STATEMENTS in v1/paper/wall_v1.tex, not from v1's code:

  C(N) = sum_{n<N} Lambda(n) mu(N-n)                       (sec:wall)
  V(N) = sum_{v<N} mu^2(v) Lambda(N-v)^2                   (prop:V)
  W(N) = sum_{w<N} Lambda(w)^2                             (prop:V)
  A(N) = prod_{q not| N} (1 - 1/(q(q-1)))                  (prop:V)
  S(N) = 2 C_2 prod_{p|N, p>2} (p-1)/(p-2)                 (Hardy-Littlewood)
  Z(N) = C(N)/sqrt(V(N))                                   (sec:floor)

Cells are the 32 divisibility patterns of N by {3,5,7,11,13}; depth d is
the number of them dividing N (sec:floor: "cells are indexed by depth d,
the number of 3,5,7,11,13 dividing N").

Both convolutions are done by FFT. Accuracy is asserted against a direct
O(N) evaluation at a sample of N before anything is cached.

Usage:  python lab_field_build.py [X]        default X = 16_000_000
Cache:  v1_verify2_log/cache/field_<X>.npz   (gitignored tree)
"""

import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np

ARTIN = 0.37395581361920228805  # prod_p (1 - 1/(p(p-1)))
HL_C2 = 0.66016181584686957393  # twin-prime constant

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CACHE = os.path.join(ROOT, "v1_verify2_log", "cache")


def smallest_prime_factor(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    spf[2::2] = 2
    i = 3
    while i * i <= X:
        if spf[i] == 0:
            sl = spf[i * i:: 2 * i]
            sl[sl == 0] = i
        i += 2
    rest = np.nonzero(spf == 0)[0]
    rest = rest[rest >= 2]
    spf[rest] = rest
    spf[0] = spf[1] = 0
    return spf


def mobius(X, primes):
    """mu by the classical sieve: one sign flip per prime divisor, then
    zero on multiples of every p^2."""
    mu = np.ones(X + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes:
        p = int(p)
        mu[p:: p] *= -1
        pp = p * p
        if pp <= X:
            mu[pp:: pp] = 0
    return mu


def von_mangoldt(X, spf):
    lam = np.zeros(X + 1, dtype=np.float64)
    is_p = np.zeros(X + 1, dtype=bool)
    is_p[2:] = spf[2:] == np.arange(2, X + 1)
    primes = np.nonzero(is_p)[0]
    lam[primes] = np.log(primes.astype(np.float64))
    for p in primes[primes <= int(X ** 0.5) + 1]:
        q = int(p) * int(p)
        lp = float(np.log(p))
        while q <= X:
            lam[q] = lp
            q *= int(p)
    return lam, primes


def fftconv_prefix(a, b, out_len):
    """First out_len+1 entries of the linear convolution of a and b."""
    L = len(a) + len(b) - 1
    n = 1
    while n < L:
        n <<= 1
    fa = np.fft.rfft(a, n)
    fb = np.fft.rfft(b, n)
    fa *= fb
    del fb
    res = np.fft.irfft(fa, n)
    del fa
    return res[: out_len + 1].copy()


def local_factors(X, primes):
    """A(N) and S(N) for every N <= X, by sieving over prime divisors."""
    inv_a = np.ones(X + 1, dtype=np.float64)   # prod_{q|N} (1-1/(q(q-1)))
    ss = np.full(X + 1, 2.0 * HL_C2, dtype=np.float64)
    for p in primes:
        p = int(p)
        f = 1.0 - 1.0 / (p * (p - 1.0))
        inv_a[p:: p] *= f
        if p > 2:
            ss[p:: p] *= (p - 1.0) / (p - 2.0)
    A = ARTIN / inv_a
    return A, ss


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 16_000_000
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f"field_{X}.npz")
    print(f"lab_field_build: X = {X:,}")

    print("  sieving smallest prime factor ...")
    spf = smallest_prime_factor(X)
    print("  von mangoldt ...")
    lam, primes = von_mangoldt(X, spf)
    del spf
    print(f"  {len(primes):,} primes")
    print("  mobius ...")
    mu = mobius(X, primes)

    print("  C = Lambda * mu   (fft) ...")
    C = fftconv_prefix(lam, mu.astype(np.float64), X)
    print("  V = Lambda^2 * mu^2  (fft) ...")
    lam2 = lam * lam
    mu2 = (mu != 0).astype(np.float64)
    V = fftconv_prefix(lam2, mu2, X)
    del mu2

    # ---- accuracy assertion, direct evaluation at a sample of N -------
    rng = np.random.default_rng(20260807)
    sample = rng.integers(1000, X, size=12)
    worstC = worstV = 0.0
    muf = mu.astype(np.float64)
    for N in sample:
        N = int(N)
        c = float(np.dot(lam[1:N], muf[N - 1:0:-1]))
        v = float(np.dot(lam2[1:N], (muf[N - 1:0:-1] != 0)))
        worstC = max(worstC, abs(c - C[N]))
        worstV = max(worstV, abs(v - V[N]) / max(v, 1.0))
    print(f"  fft check: max |C_fft - C_direct| = {worstC:.3e}")
    print(f"  fft check: max rel err in V       = {worstV:.3e}")
    if worstC > 1e-3 or worstV > 1e-12:
        print("  FFT ACCURACY FAILURE")
        return 1
    del lam2, muf

    print("  W = cumulative sum of Lambda^2 ...")
    W = np.cumsum(lam * lam)

    print("  local factors A(N), S(N) ...")
    A, S = local_factors(X, primes)

    print("  cells: divisibility pattern by {3,5,7,11,13} ...")
    N_even = np.arange(2, X + 1, 2, dtype=np.int64)
    cell = np.zeros(N_even.shape, dtype=np.int8)
    depth = np.zeros(N_even.shape, dtype=np.int8)
    for i, p in enumerate((3, 5, 7, 11, 13)):
        hit = (N_even % p) == 0
        cell |= (hit.astype(np.int8) << i)
        depth += hit.astype(np.int8)

    print(f"  saving {path}")
    np.savez_compressed(
        path,
        X=np.int64(X),
        N=N_even,
        C=C[N_even],
        V=V[N_even],
        W=W[N_even - 1],
        A=A[N_even],
        S=S[N_even],
        cell=cell,
        depth=depth,
    )
    print(f"  done. {len(N_even):,} even N")
    print()
    print(f"  even N <= {X:,} : {len(N_even):,}")
    print(f"  (conj:wall item 1 quotes '6.3e6 values' for this field)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
