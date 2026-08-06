# -*- coding: utf-8 -*-
"""
Hypothesis round 2, redone against explicit noise floors
(increment 213).

Round 2 returned ALIVE on both hypotheses and both verdicts were
artefacts of badly chosen thresholds:

  H10 used |corr| >= 0.10 with 300 points, where the null SE of a
  Pearson correlation is 1/sqrt(300) = 0.058 -- so the threshold sat
  at 1.7 sigma, below the noise floor, and two hits out of twelve
  tests is what chance delivers (expected 1.1).

  H12 used "one decade of n carries >= 50% of the total", but the top
  decade is always the last one, which spans ~89% of the range. Under
  Brownian behaviour a range of length L contributes ~sqrt(L), so the
  last decade's expected share is ~0.95; the measured 0.53-0.74 is
  therefore LESS than random, not more.

Both are redone here with the null computed first and the statistic
chosen so that "no structure" has a known value.

H10' : report z = corr * sqrt(n) over many (N, h, K) triples and
       compare the collection to a standard normal.
       ALIVE iff mean|z| >= 2.5 or any |z| >= 4 replicated at a second
       band; DEAD if the z's look standard normal.

H12' : for equal-log-length ranges of n, report
       g = |Sum_range f| / sqrt(#terms), which is constant (= rms f)
       under Brownian behaviour.
       ALIVE iff some range has g >= 3x the median g across ranges,
       at the same range for a majority of N; DEAD if g is flat.
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
    print("--- H10': cross-N correlation, against the null SE ---")
    zs = []
    print(f"{'N':>10} {'h':>4} {'K':>6} {'#k':>6} {'corr':>8} {'z':>7}")
    for N in bases:
        for K, nk in bands:
            ks = np.arange(K, K + nk, dtype=np.int64)
            D0 = field(mu, N, ks)
            for h in hs:
                D1 = field(mu, N + h, ks)
                c = float(np.corrcoef(D0, D1)[0, 1])
                z = c * math.sqrt(len(ks))
                zs.append(z)
                print(f"{N:>10} {h:>4} {K:>6} {len(ks):>6} "
                      f"{c:>8.4f} {z:>7.2f}")
    z = np.array(zs)
    print(f"  n tests = {z.size}   mean z = {z.mean():+.3f}   "
          f"mean |z| = {np.abs(z).mean():.3f} (null 0.798)   "
          f"max |z| = {np.abs(z).max():.2f}")
    return float(np.abs(z).mean()), float(np.abs(z).max())


def h12(mu, lam, Ns, nr=8):
    print("\n--- H12': normalised contribution per range ---")
    print(f"{'N':>9} " + " ".join(f"{i:>7}" for i in range(nr))
          + f" {'max/med':>8}")
    ratios = []
    argmaxes = []
    for N in Ns:
        idx = np.arange(1, N)
        f = lam[1:N] * mu[N - idx].astype(np.float64)
        # equal-log-length ranges from N^(1/4) to N
        lo = int(N ** 0.25)
        edges = np.unique(np.geomspace(lo, N, nr + 1).astype(np.int64))
        gs = []
        for a, b in zip(edges[:-1], edges[1:]):
            seg = f[a - 1:b - 1]
            if seg.size < 10:
                gs.append(float('nan')); continue
            gs.append(abs(float(seg.sum())) / math.sqrt(seg.size))
        g = np.array(gs)
        med = float(np.nanmedian(g))
        r = float(np.nanmax(g) / med) if med > 0 else float('nan')
        ratios.append(r)
        argmaxes.append(int(np.nanargmax(g)))
        print(f"{N:>9} " + " ".join(f"{v:>7.3f}" for v in g)
              + f" {r:>8.2f}")
    print(f"  max/median by N: {[f'{r:.2f}' for r in ratios]}   "
          f"argmax ranges: {argmaxes}")
    same = max(argmaxes.count(a) for a in set(argmaxes))
    return ratios, same, len(Ns)


def main():
    N1 = 9_999_998
    mu = mobius_upto(N1 + 64)
    mz, mx = h10(mu, [N1, 9_999_942, 9_999_866], [2, 6, 30],
                 [(500, 1200), (1500, 1200)])
    del mu

    X = 900_000
    mu2, lam2, phi2, primes2, spf2 = sieve(X)
    ratios, same, nN = h12(mu2, lam2,
                           [199_998, 399_998, 599_998, 899_998])

    print("\nverdicts:")
    print("  H10':", "ALIVE" if (mz >= 2.5 or mx >= 4) else
          "DEAD -- the z's are standard normal; nearby N carry "
          "independent copies of the field")
    alive12 = (max(ratios) >= 3 and same >= 3)
    print("  H12':", "ALIVE" if alive12 else
          "DEAD -- the normalised contribution is flat across ranges; "
          "the walk is Brownian and no range carries the mass")
    print("DONE")


if __name__ == "__main__":
    main()
