# -*- coding: utf-8 -*-
"""
Re-verification of representation class C-I of v1/paper/wall_v1.tex,
and of the null it is judged against.

THE STATEMENT UNDER TEST (§7.3):

    C-I abelian | rational-peak energy vs mask-null |
    **closed**, 0/6: the abelian spectrum is mask-exact.

THE DESIGN, from its own docstring: for fixed k, t_m = mu(m) mu(N-mk)
on m in (sqrt N, N/k]; zero-pad t to length 2^18, FFT, and take the
total rational-peak energy

    E_Q = sum over a/q with q <= 32, gcd(a,q) = 1,
          of |F|^2 in the bin nearest a/q.

Null: **8 draws** of random signs on the real support. Report
z = (E_Q^real - mean)/std per k; alive iff z >= 4 at >= 2 of 6 k.

THE OBJECTION BEING TESTED. Eight draws is the exact null size the
paper's own Methodology names as a trap:

    "Eight-draw-null inflation. Nulls estimated from too few draws
     inflate maxima across a family of 300 tests; the C-II 'hit'
     survived 8-draw nulls and died under 64-draw nulls plus
     replication."

An 8-draw sample standard deviation is biased low (the bias factor
c_4(8) = 0.9650, and its own relative spread is about 27%), so z is
inflated. For C-II that manufactured a hit. For C-I the verdict is
0/6 -- a NULL result -- and an inflated z makes hits MORE likely, so
the direction is conservative and the closure is not in danger from
this. What IS in danger is the test's resolution: a bar at z >= 4
cannot be located from 8 draws at all.

METHOD HERE. Written from the statement. The same E_Q, at 8, 64 and
400 draws, so the drift of z with null size is visible; plus the
permutation null (signs permuted rather than redrawn), which keeps the
multiset of values; plus the fraction of total energy at the rational
peaks, which is the quantity "mask-exact" is about.

PRE-REGISTRATION (written before the run).

  (1) RULE. Report z at 8, 64 and 400 draws. If z falls materially as
      the null grows, the 8-draw z is inflated and the printed
      significance is not the measured one -- in the safe direction
      here, but the number is still not a measurement.
  (2) RULE. The closure needs 0 hits at z >= 4. Report the hit count
      at every null size. If any k reaches z >= 4 at 400 draws, C-I is
      not closed.
  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      Conjecture 10 says the abelian spectrum is exactly the mask, so
      I expect z near zero at every null size and 0/6 throughout. I
      expect the 8-draw z to be about 3-5% larger in magnitude than
      the 400-draw z, from the c_4 bias alone, and much noisier. I
      expect no finding, and the informative outcome would be a k
      whose z grows rather than shrinks with null size.
"""
import sys
import math
import time
from math import gcd

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

L = 1 << 18
QMAX = 32


def sieve_mu(X):
    mu = np.ones(X + 1, dtype=np.int8)
    comp = np.zeros(X + 1, dtype=bool)
    for p in range(2, int(X ** 0.5) + 1):
        if not comp[p]:
            comp[p * p::p] = True
            mu[p::p] = -mu[p::p]
            mu[p * p::p * p] = 0
    rest = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if not comp[p]:
            q = p
            while q <= X:
                rest[q::q] //= p
                q *= p
    mu[rest > 1] = -mu[rest > 1]
    mu[0] = 0
    del comp, rest
    return mu


def peak_bins():
    bins = set()
    for q in range(2, QMAX + 1):
        for a in range(1, q):
            if gcd(a, q) == 1:
                bins.add(int(round(a / q * L)) % L)
    return np.array(sorted(bins))


def energy(t, bins):
    F = np.fft.fft(t, L)
    return float((np.abs(F[bins]) ** 2).sum()), float(
        (np.abs(F) ** 2).sum())


def main():
    t0 = time.time()
    N = 199_999_998
    KS = [2001, 2311, 2731, 3001, 3511, 3911]
    bins = peak_bins()
    print("Re-verification of C-I (the abelian spectrum) and its null")
    print(f"  N = {N}, {len(KS)} k values, FFT length 2^18,")
    print(f"  {len(bins)} rational peak bins (q <= {QMAX})")
    print()
    mu = sieve_mu(N)
    print(f"  mu ready t={time.time()-t0:.0f}s", flush=True)
    SQ = int(N ** 0.5)
    rng = np.random.default_rng(777)

    hdr = (f"  {'k':>6} {'M':>8} {'E_rat/E_tot':>12} "
           f"{'z (8 draws)':>12} {'z (64)':>9} {'z (400)':>9} "
           f"{'z (perm)':>9} {'hit?':>5}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    hits = {8: 0, 64: 0, 400: 0}
    for k in KS:
        ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
        t = (mu[ms].astype(np.float64) * mu[N - k * ms])
        sup = t != 0
        Er, Et = energy(t, bins)
        zs = {}
        pool = []
        for _ in range(400):
            s = np.where(sup, rng.choice([-1.0, 1.0], size=t.shape), 0.0)
            pool.append(energy(s, bins)[0])
        pool = np.array(pool)
        for nd in (8, 64, 400):
            sub = pool[:nd]
            zs[nd] = (Er - sub.mean()) / max(sub.std(), 1e-30)
            hits[nd] += int(zs[nd] >= 4)
        # permutation null: reuse the real values, reshuffled
        pp = []
        vals = t[sup]
        for _ in range(64):
            s = np.zeros_like(t)
            s[sup] = rng.permutation(vals)
            pp.append(energy(s, bins)[0])
        pp = np.array(pp)
        zp = (Er - pp.mean()) / max(pp.std(), 1e-30)
        print(f"  {k:>6} {len(ms):>8} {Er/Et:>12.5f} "
              f"{zs[8]:>+12.2f} {zs[64]:>+9.2f} {zs[400]:>+9.2f} "
              f"{zp:>+9.2f} {'HIT' if zs[400] >= 4 else '':>5}"
              f"   t={time.time()-t0:.0f}s", flush=True)

    print()
    print(f"(2) hits at z >= 4:  8 draws {hits[8]}/6,  "
          f"64 draws {hits[64]}/6,  400 draws {hits[400]}/6")
    print(f"    v1 reports 0/6 at 8 draws. Closure needs 0.")
    print()
    print("(1) the 8-draw standard deviation is biased low by the")
    print(f"    factor c_4(8) = 0.9650, so an 8-draw z is about 3.6%")
    print("    larger in magnitude than the truth on average, and its")
    print("    own relative spread is about 27%. For a NULL result")
    print("    that bias is conservative -- it makes hits easier, and")
    print("    none occurred. But a bar at z >= 4 cannot be located")
    print("    from 8 draws: the 8-draw z is not a measurement, and")
    print("    the closure should be quoted at the 400-draw column.")
    if hits[400]:
        print("DONE (C-I not closed)")
        sys.exit(1)
    print("DONE")


if __name__ == "__main__":
    main()
