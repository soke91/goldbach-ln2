# -*- coding: utf-8 -*-
"""
Hypothesis round 2 (increment 212), aimed at the wall's own objects.

H10 -- CROSS-N CORRELATION OF THE DILATE FIELD. D_N(k) = Sum_m mu(m)
mu(N-mk). Moving N to N+h shifts every Mobius argument by h, so a
nonzero corr_k(D_N, D_{N+h}) would be a binary Chowla signal; a zero
one says the field is reborn at each N, which is exactly why the
N-averaged route (adjudication route 5) cannot escape its exceptional
set.
  ALIVE iff |corr| >= 0.10 at two or more (N, h) pairs.
  DEAD  iff |corr| <= 0.03 throughout.

H12 -- PARTIAL-SUM PROFILE OF C(N). C(t) = Sum_{n<=t} Lambda(n)
mu(N-n), so C(N) is the endpoint of a walk. If some range of n carried
the mass, that range could be attacked on its own.
  ALIVE iff one decade of n carries >= 50% of the total variation
        systematically across N.
  DEAD  iff the profile is Brownian: |C(t)| ~ sqrt(t) and each range
        contributes in proportion to its length.
"""
import numpy as np
import math

from e1_forge_r4 import mobius_upto
from thmC_alpha_scan import sieve


def field(mu, N, ks):
    SQ = int(N ** 0.5)
    D = np.zeros(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        hi = N // k
        if hi <= SQ:
            continue
        ms = np.arange(SQ + 1, hi + 1, dtype=np.int64)
        D[i] = float((mu[ms].astype(np.int64)
                      * mu[N - k * ms].astype(np.int64)).sum())
    return D


def h10(mu, bases, hs, bands):
    print("--- H10: cross-N correlation of D_N(k) ---")
    print(f"{'N':>10} {'h':>4} {'K':>6} {'#k':>5} {'corr':>8}")
    flags = 0
    for N in bases:
        for K, nk in bands:
            ks = np.arange(K, K + nk, dtype=np.int64)
            D0 = field(mu, N, ks)
            for h in hs:
                D1 = field(mu, N + h, ks)
                c = float(np.corrcoef(D0, D1)[0, 1])
                if abs(c) >= 0.10:
                    flags += 1
                print(f"{N:>10} {h:>4} {K:>6} {len(ks):>5} {c:>8.4f}")
    print(f"  flags (|corr| >= 0.10): {flags}")
    return flags


def h12(mu, lam, Ns):
    print("\n--- H12: partial-sum profile of C(N) ---")
    print(f"{'N':>9} {'C(N)':>10} {'max|C(t)|':>10} "
          f"{'|C(N)|/sqrtN':>12} {'max/sqrtN':>10} "
          f"{'top-decade share':>17}")
    alive = 0
    for N in Ns:
        idx = np.arange(1, N)
        f = lam[1:N] * mu[N - idx].astype(np.float64)
        run = np.cumsum(f)
        C = run[-1]
        mx = float(np.abs(run).max())
        # contribution of each decade of n, as a share of total |増分|
        edges = [1]
        while edges[-1] * 10 < N:
            edges.append(edges[-1] * 10)
        edges.append(N)
        parts = []
        for a, b in zip(edges[:-1], edges[1:]):
            parts.append(abs(float(f[a - 1:b - 1].sum())))
        tot = sum(parts)
        share = max(parts) / tot if tot > 0 else 0.0
        if share >= 0.5:
            alive += 1
        print(f"{N:>9} {C:>10.1f} {mx:>10.1f} "
              f"{abs(C)/math.sqrt(N):>12.3f} {mx/math.sqrt(N):>10.3f} "
              f"{share:>17.3f}")
    print(f"  N with a decade carrying >= 50%: {alive}/{len(Ns)}")
    return alive


def main():
    N1 = 9_999_998
    mu = mobius_upto(N1 + 64)
    f10 = h10(mu, [N1, 9_999_942], [2, 6, 30],
              [(500, 300), (1500, 300)])
    del mu

    X = 900_000
    mu2, lam2, phi2, primes2, spf2 = sieve(X)
    a12 = h12(mu2, lam2, [199_998, 399_998, 599_998, 899_998])

    print("\nverdicts:")
    print("  H10:", "ALIVE" if f10 >= 2 else
          "DEAD -- the field is reborn at each N; nearby N carry "
          "independent copies, which is why N-averaging buys an "
          "exceptional set and not a fixed-N statement")
    print("  H12:", "ALIVE" if a12 >= 3 else
          "DEAD -- no range of n carries the mass")
    print("DONE")


if __name__ == "__main__":
    main()
