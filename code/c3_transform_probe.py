# -*- coding: utf-8 -*-
"""
C-III item (1): does a legitimate dilate -> shift transform exist?
(increment 201)

Re-audit round 2 re-opened C-III and stated what it needs. Item (1) is
a legitimate transform in place of the draft's refuted hybrid object.
This script pins down the geometry and tests the one empirical
question the design raises.

THE GEOMETRY. In the L^2 off-diagonal
    OD = Sum_{m != m'} mu(m) mu(m') Sum_{k~K} mu(N-mk) mu(N-m'k),
put a = N - mk and b = N - m'k. Then m' a - m b = (m'-m) N: as k runs,
the pair (a,b) traverses lattice points on a fixed LINE of rational
slope m/m'. A shift correlation Sum_a mu(a) mu(a+h) is the same object
on a line of slope 1. So:

    dilate problem  =  mu-mu correlations along lines of arbitrary
                       rational slope,
    shift problem   =  the slope-1 case.

Item (1) asks for a transform carrying the first to the second.

THE CIRCULARITY CHECK (algebra, verified below). The natural attempt
re-indexes by the difference d = m' - m and fixes (d, k). Then the
shift h = dk is FIXED and a = N - mk runs over the progression
a = N (mod k), so the inner factor really is a fixed-shift correlation
    Sum_{a = N (mod k)} mu(a) mu(a - dk).
But the shift dk is a MULTIPLE OF THE MODULUS k, and rescaling the
progression by k, a = N - mk, sends mu(a)mu(a-dk) back to
mu(N-mk) mu(N-(m-d)k) -- the dilate object we started from. The shift
reading exists and is exactly circular. This script verifies the
decomposition identity numerically so the algebra is not taken on
trust.

THE EMPIRICAL QUESTION. Is the slope family, in aggregate, any tamer
than the shift family at matched parameters? If slope correlations
were measurably smaller than shift correlations, a transform would be
buying something real and worth constructing. If they match, the two
families are the same difficulty and item (1) is a pure
provability question, not a structural one.

PRE-REGISTERED (fixed before the run):
  LEAD   iff mean |T(m,m')| / sqrt(K) over the slope family is <= 0.6,
         i.e. measurably below the half-normal value 0.798 that the
         shift family sits at.
  CLOSED iff the two means agree to within 10%: same difficulty, and
         item (1) cannot be motivated by the slope family being tamer.
"""
import numpy as np
import math

from e1_forge_r4 import mobius_upto

HALF_NORMAL = math.sqrt(2.0 / math.pi)   # 0.7979


def slope_family(mu, N, K, nk, ms, rng):
    """|T(m,m')| / sqrt(K) for random coprime-ish pairs (m,m')."""
    ks = np.arange(K, K + nk, dtype=np.int64)
    out = []
    for _ in range(len(ms)):
        m1, m2 = rng.choice(ms, size=2, replace=False)
        if m1 == m2:
            continue
        hi = min(N // int(m1), N // int(m2))
        kk = ks[ks < hi]
        if kk.size < 32:
            continue
        a = mu[N - int(m1) * kk].astype(np.int64)
        b = mu[N - int(m2) * kk].astype(np.int64)
        t = float((a * b).sum())
        out.append(abs(t) / math.sqrt(kk.size))
    return np.array(out)


def shift_family(mu, X0, L, hs):
    """|Sum_{a in [X0, X0+L)} mu(a)mu(a+h)| / sqrt(L) for each h."""
    a = np.arange(X0, X0 + L, dtype=np.int64)
    base = mu[a].astype(np.int64)
    out = []
    for h in hs:
        t = float((base * mu[a + int(h)].astype(np.int64)).sum())
        out.append(abs(t) / math.sqrt(L))
    return np.array(out)


def identity_check(mu, N, K, nk, ms):
    """OD computed two ways: directly, and via the (d,k) re-indexing."""
    ks = np.arange(K, K + nk, dtype=np.int64)
    direct = 0.0
    for m1 in ms:
        for m2 in ms:
            if m1 == m2:
                continue
            hi = min(N // int(m1), N // int(m2))
            kk = ks[ks < hi]
            if kk.size == 0:
                continue
            direct += (mu[m1] * mu[m2]
                       * float((mu[N - int(m1) * kk].astype(np.int64)
                                * mu[N - int(m2) * kk].astype(np.int64)
                                ).sum()))
    # via d = m2 - m1, shift h = d*k on the progression a = N mod k
    viad = 0.0
    for m1 in ms:
        for m2 in ms:
            if m1 == m2:
                continue
            d = int(m2) - int(m1)
            hi = min(N // int(m1), N // int(m2))
            kk = ks[ks < hi]
            if kk.size == 0:
                continue
            a = N - int(m1) * kk                      # a = N mod k
            viad += (mu[m1] * mu[m2]
                     * float((mu[a].astype(np.int64)
                              * mu[a - d * kk].astype(np.int64)).sum()))
    return direct, viad


def main():
    N = 99_999_998
    rng = np.random.default_rng(20260906)
    mu = mobius_upto(N)

    print("=== circularity check: the (d,k) shift reading is an "
          "identity, not a transform ===")
    ms_small = np.array([m for m in range(10001, 10061)
                         if mu[m] != 0], dtype=np.int64)
    d1, d2 = identity_check(mu, N, 1000, 200, ms_small[:24])
    print(f"  OD direct            = {d1:.1f}")
    print(f"  OD via (d,k) reading = {d2:.1f}")
    print(f"  difference           = {abs(d1-d2):.3e}   "
          f"{'IDENTICAL (as predicted: the shift reading is a '
             'relabelling)' if abs(d1-d2) < 1e-6 else 'MISMATCH'}")

    print("\n=== slope family vs shift family at matched parameters ===")
    print(f"{'K':>7} {'#pairs':>7} {'slope mean':>11} "
          f"{'#shifts':>8} {'shift mean':>11} {'half-normal':>12} "
          f"{'ratio':>7}")
    for K, nk in ((1000, 512), (3000, 512)):
        ms = np.array([m for m in range(int(N ** 0.5) + 1,
                                        int(N ** 0.5) + 4000)
                       if mu[m] != 0], dtype=np.int64)
        sl = slope_family(mu, N, K, nk, ms[:600], rng)
        hs = rng.integers(1, 4000, size=400)
        sh = shift_family(mu, 40_000_000, nk, hs)
        r = sl.mean() / sh.mean()
        print(f"{K:>7} {sl.size:>7} {sl.mean():>11.4f} "
              f"{sh.size:>8} {sh.mean():>11.4f} {HALF_NORMAL:>12.4f} "
              f"{r:>7.4f}")

    print("\nPre-registered: LEAD if slope mean <= 0.6; CLOSED if the "
          "two agree within 10%.")
    print("DONE")


if __name__ == "__main__":
    main()
