# -*- coding: utf-8 -*-
"""
Re-measure prop:V, sec:margin, and conj:wall items 1 and 3.
(v1_verify2, Phase 1, blind -- written from the statements only.)

STATEMENTS UNDER TEST, from v1/paper/wall_v1.tex:

  prop:V   V(N) = W(N) A(N) (1+o(1)), A(N) = prod_{q not| N}(1-1/(q(q-1))).
           "the residual standard deviation over every even N <= 1.6e7 is
            0.000323 for A against 0.245235 for S -- a factor of 760"
           "the ratio V(N)/W(N) against A(N) has mean 1.000000 with
            standard deviation 0.000145 in the top octave, and agrees to
            five decimals in each of the six radical cells"

  conj:wall item 1
           "Excess kurtosis -0.0005 (z=-0.3) and E|G|/sd(G) short of
            sqrt(2/pi) by 0.00018 (z=-0.8), on 6.3e6 values, removing
            cell means alone. Under an S N-based scale the same data give
            excess kurtosis +0.1704 at z=98."

  conj:wall item 3
           "max|G| tracks the Gumbel law with mean deviation +0.54+-0.45
            from E[max] over eight octaves; aggregate tail counts against
            the Gaussian expectation give ratios 0.999 at t=3, 0.997 at
            t=4, 0.878 at t=5"

  sec:margin
           "measured, max|C|/N falls from 0.056 to 0.0082 over the range
            computed, and the margin at N=1e8 is a factor N^{0.454}"

PRE-REGISTRATION (fixed before this ran).

  Decision rule, per item: recompute under the reading the paper's words
  force. Where the words admit more than one reading, compute ALL of
  them and print each, because finding 4 of the first-pass rules ("define
  the statistic") says a number with two defensible readings is the
  failure mode here.
    REPRODUCED     : some stated reading returns the paper's figure to
                     the precision the paper quotes.
    NOT REPRODUCED : no reading does. Then report the spread across
                     readings, so the reader can see whether the
                     disagreement is a definition or a defect.

  Predictions written before running:
   (a) prop:V's A-vs-S contrast REPRODUCES. It is a large, robust effect
       and the local-factor derivation is elementary and I checked it by
       hand.
   (b) item 1's kurtosis REPRODUCES in sign and rough size but the
       sample size will come out 8.0e6, not the 6.3e6 quoted, because
       every even N <= 1.6e7 is 8.0e6 values. (This is already logged as
       blind finding B5 from the z=98 arithmetic; here it is measured.)
   (c) item 3's t=5 tail ratio will have a Poisson spread of order 50%,
       making the quoted 0.878 uninformative.
   (d) sec:margin's max|C|/N endpoints REPRODUCE but N^{0.454} does not
       (already logged as B7 from the paper's own formula).

  What would refute each: the corresponding figure coming out inside the
  paper's quoted precision under a reading stated in the paper.
"""

import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CACHE = os.path.join(ROOT, "v1_verify2_log", "cache")


def phi(z):
    from math import erfc, sqrt
    return 0.5 * erfc(-z / sqrt(2.0))


def gumbel_max_normal(n):
    """E[max of n iid standard normals], and the Gumbel scale."""
    import math
    ln = math.log(n)
    b = math.sqrt(2 * ln) - (math.log(ln) + math.log(4 * math.pi)) / (
        2 * math.sqrt(2 * ln))
    a = 1.0 / math.sqrt(2 * ln)
    return b + 0.5772156649015329 * a, a


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 16_000_000
    path = os.path.join(CACHE, f"field_{X}.npz")
    if not os.path.exists(path):
        print(f"missing cache {path}; run lab_field_build.py {X} first")
        return 2
    z = np.load(path)
    Ni_all = z["N"]
    n_all = len(Ni_all)

    # Rule 5 of this tree ("check whether the design was run"): count how
    # much of the design is non-degenerate before believing any average.
    good = (z["V"] > 0) & (z["W"] > 0)
    ndeg = n_all - int(good.sum())

    Ni = Ni_all[good]
    N = Ni.astype(np.float64)
    C = z["C"][good]
    V = z["V"][good]
    W = z["W"][good]
    A = z["A"][good]
    S = z["S"][good]
    cell = z["cell"][good]
    depth = z["depth"][good]
    n = len(N)

    print("audit_propV_and_wall  (v1_verify2 Phase 1, blind)")
    print(f"field: every even N <= {X:,}   ->   {n_all:,} values")
    print(f"       degenerate (V=0 or W=0, i.e. N <= 4): {ndeg}")
    print(f"       non-degenerate field: {n:,} values")
    print(f"       [conj:wall item 1 quotes '6.3e6 values' for this field;")
    print(f"        that is {n / 6.3e6:.3f}x smaller than what is here]")
    print("=" * 72)
    print()

    # ================================================== prop:V
    print("--- prop:V: is the local factor A, and how well ------------------")
    ratio = V / W                       # the paper's "V(N)/W(N) against A(N)"
    rA = ratio / A
    rS = ratio / S
    top = N > X / 2.0                   # the top octave
    print(f"  top octave: {top.sum():,} values (N > {X / 2:,.0f})")
    print(f"  V/W divided by A(N):  mean {rA[top].mean():.6f}   "
          f"sd {rA[top].std(ddof=1):.6f}   [paper: mean 1.000000, sd 0.000145]")
    print(f"  V/W divided by S(N):  mean {rS[top].mean():.6f}   "
          f"sd {rS[top].std(ddof=1):.6f}")
    print()

    print("  'rescaling each candidate to the measured mean, so that only")
    print("   its shape in N is judged' -- residual sd over every even N:")
    for nm, r in (("A", rA), ("S", rS)):
        rr = r / r.mean()
        print(f"    {nm}: residual sd of (V/W)/{nm}, rescaled to mean 1 "
              f"= {rr.std(ddof=1):.6f}")
    # second reading: compare V against candidate * N log N, shape only
    for nm, cand in (("A", A), ("S", S)):
        pred = cand * N * np.log(N)
        rr = V / pred
        rr = rr / rr.mean()
        print(f"    {nm}: residual sd of V/({nm}*N*logN), rescaled "
              f"= {rr.std(ddof=1):.6f}")
    sdA = (rA / rA.mean()).std(ddof=1)
    sdS = (rS / rS.mean()).std(ddof=1)
    print(f"  paper: 0.000323 for A, 0.245235 for S, factor 760")
    print(f"  here : {sdA:.6f} for A, {sdS:.6f} for S, "
          f"factor {sdS / sdA:.1f}")
    print()

    print("  the six radical cells (2|N, 2*3|N, ..., 2*3*5*7*11*13|N):")
    rad = 1
    for p in (1, 3, 5, 7, 11, 13):
        rad *= p
        m = (Ni % (2 * rad)) == 0
        print(f"    2*{rad:<6} | N :  {m.sum():>9,} values   "
              f"mean (V/W)/A = {rA[m].mean():.7f}")
    print()

    # ================================================== sec:margin
    print("--- sec:margin: max|C|/N over the range --------------------------")
    absC = np.abs(C)
    print(f"  {'octave':>22} {'count':>10} {'max|C|/N':>12} {'max|C|':>14}")
    oct_lo = X
    rows = []
    while oct_lo > 60000:
        lo, hi = oct_lo / 2.0, oct_lo
        m = (N > lo) & (N <= hi)
        if m.sum() == 0:
            break
        i = np.argmax(absC[m] / N[m])
        val = (absC[m] / N[m])[i]
        rows.append((lo, hi, m.sum(), val))
        print(f"  {lo:>10.3g}-{hi:<10.3g} {m.sum():>10,} {val:>12.5f} "
              f"{absC[m].max():>14.1f}")
        oct_lo /= 2.0
    print()
    print(f"  paper: 'max|C|/N falls from 0.056 to 0.0082 over the range "
          f"computed'")
    print(f"  here : {rows[-1][3]:.4f} at the bottom octave -> "
          f"{rows[0][3]:.4f} at the top")
    import math
    print()
    print("  other readings of 'max|C|/N', since the paper does not say")
    print("  over what window the max is taken:")
    for name, lo, hi in (("N <= 1e5", 0, 1e5), ("N <= 1.4e7", 0, 1.4e7),
                         ("N <= 1.6e7", 0, 1.6e7),
                         ("1e5 < N <= 1.4e7", 1e5, 1.4e7)):
        m = (N > lo) & (N <= hi)
        if m.sum() == 0:
            continue
        mx = float((absC[m] / N[m]).max())
        print(f"    max|C|/N over {name:<18} = {mx:.5f}")
    for name, lo, hi in (("N <= 1e5", 0, 1e5), ("N <= 1.4e7", 0, 1.4e7)):
        m = (N > lo) & (N <= hi)
        mx = float(absC[m].max())
        print(f"    max|C| over   {name:<18} = {mx:.1f}   "
              f"/hi = {mx / hi:.5f}")
    print()
    gmargin = 1.0 / rows[0][3]
    print(f"  margin N/max|C| at the top octave = {gmargin:.1f} = "
          f"N^{math.log(gmargin) / math.log(X):.4f}")
    print(f"  paper claims the margin at N=1e8 is N^0.454 = "
          f"{1e8 ** 0.454:.0f}")
    print("  fitting max|C|/N ~ N^-a across the octaves above:")
    lx = np.log(np.array([0.5 * (r[0] + r[1]) for r in rows]))
    ly = np.log(np.array([r[3] for r in rows]))
    a_fit = -np.polyfit(lx, ly, 1)[0]
    print(f"    a = {a_fit:.4f}  ->  margin ~ N^{a_fit:.4f}, "
          f"and at N=1e8 the extrapolated margin is "
          f"{(1.0 / rows[0][3]) * (1e8 / (0.75 * X)) ** a_fit:.0f}")
    print()

    # ================================================== conj:wall item 1
    print("--- conj:wall item 1: the Gaussian bulk --------------------------")
    Z = C / np.sqrt(V)

    def demean(field, labels, k):
        out = field.copy()
        for c in range(k):
            m = labels == c
            if m.sum():
                out[m] -= out[m].mean()
        return out

    # "cells" is not defined in the paper beyond "indexed by depth d, the
    # number of 3,5,7,11,13 dividing N". Both readings are computed.
    print("  two readings of 'cell', since the paper defines only the index:")
    for nm, lab, k in (("32 divisibility patterns", cell, 32),
                       ("6 depths", depth, 6)):
        g = demean(Z, lab, k)
        gs = g / g.std(ddof=1)
        e = (gs ** 4).mean() - 3.0
        print(f"    {nm:<26}: excess kurtosis {e:+.6f}  "
              f"z={e / np.sqrt(24.0 / n):+.1f}")
    print()

    G = demean(Z, cell, 32)
    sd = G.std(ddof=1)
    Gs = G / sd
    ek = (Gs ** 4).mean() - 3.0
    se_ek = np.sqrt(24.0 / n)
    absmean = np.abs(Gs).mean()
    target = np.sqrt(2.0 / np.pi)
    se_abs = np.sqrt(1 - 2 / np.pi) / np.sqrt(n)
    print(f"  n = {n:,}  (paper says 6.3e6)")
    print(f"  excess kurtosis            = {ek:+.6f}   "
          f"z = {ek / se_ek:+.2f}   [paper: -0.0005, z=-0.3]")
    print(f"  E|G|/sd(G)                 = {absmean:.6f}   "
          f"sqrt(2/pi) = {target:.6f}")
    print(f"  shortfall                  = {absmean - target:+.6f}   "
          f"z = {(absmean - target) / se_abs:+.2f}   [paper: -0.00018, z=-0.8]")
    print()
    # the S*N-based scale
    ZS = C / np.sqrt(S * N)
    GS = demean(ZS, cell, 32)
    GSs = GS / GS.std(ddof=1)
    ekS = (GSs ** 4).mean() - 3.0
    print(f"  under the S(N)*N scale:")
    print(f"  excess kurtosis            = {ekS:+.6f}   "
          f"z = {ekS / se_ek:+.1f}   [paper: +0.1704, z=98]")
    print()
    # and with no cell-mean removal at all, for contrast
    Zr = (Z - Z.mean()) / Z.std(ddof=1)
    print(f"  (no cell-mean removal: excess kurtosis "
          f"{(Zr ** 4).mean() - 3.0:+.6f})")
    print()

    # ================================================== conj:wall item 3
    print("--- conj:wall item 3: the tail -----------------------------------")
    print(f"  {'t':>4} {'observed':>10} {'expected':>12} {'ratio':>8} "
          f"{'Poisson se of ratio':>21}")
    for t in (3, 4, 5, 6):
        obs = int((np.abs(Gs) > t).sum())
        exp = 2 * phi(-t) * n
        r = obs / exp if exp > 0 else float("nan")
        se_r = np.sqrt(max(obs, 1)) / exp
        print(f"  {t:>4} {obs:>10,} {exp:>12.1f} {r:>8.3f} "
              f"{'+-' + format(se_r, '.3f'):>21}")
    print(f"  [paper: 0.999 at t=3, 0.997 at t=4, 0.878 at t=5, "
          f"no spread quoted]")
    print()

    print("  max|G| per octave against the Gumbel expectation:")
    devs = []
    oct_lo = X
    print(f"  {'octave':>22} {'n_oct':>10} {'max|G|':>9} {'E[max]':>9} "
          f"{'dev':>8} {'dev/scale':>10}")
    while oct_lo > 60000:
        lo, hi = oct_lo / 2.0, oct_lo
        m = (N > lo) & (N <= hi)
        if m.sum() == 0:
            break
        no = int(m.sum())
        emax, scale = gumbel_max_normal(2 * no)   # |G|: two-sided
        obs = float(np.abs(Gs[m]).max())
        devs.append(obs - emax)
        print(f"  {lo:>10.3g}-{hi:<10.3g} {no:>10,} {obs:>9.3f} "
              f"{emax:>9.3f} {obs - emax:>+8.3f} "
              f"{(obs - emax) / scale:>+10.2f}")
        oct_lo /= 2.0
    devs = np.array(devs)
    print(f"  mean deviation over {len(devs)} octaves = {devs.mean():+.3f} "
          f"+- {devs.std(ddof=1) / np.sqrt(len(devs)):.3f}   "
          f"[paper: +0.54 +- 0.45 over eight octaves]")
    print()

    print("  extremes: are they at deep radicals or generic N?")
    idx = np.argsort(-np.abs(Gs))[:20]
    dd = depth[idx]
    print(f"    depth of the 20 largest |G|: "
          f"{np.bincount(dd, minlength=6).tolist()} (by depth 0..5)")
    print(f"    depth distribution of the whole field: "
          f"{(np.bincount(depth, minlength=6) / n).round(4).tolist()}")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
