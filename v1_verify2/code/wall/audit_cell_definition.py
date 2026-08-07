# -*- coding: utf-8 -*-
"""
conj:wall items 1 and 3: does a finer "cell" rescue the quoted numbers?
(v1_verify2, Phase 1, blind.)

Established so far by audit_propV_readings.py: the paper's field is
10^5 <= N <= 1.6e7 (all three of prop:V's figures -- 0.000323, 0.245235,
760 -- reproduce to full quoted precision there, and nowhere else). So
the pipeline here is the paper's pipeline.

On that field, "removing cell means alone" with cells = the 32
divisibility patterns of {3,5,7,11,13} gives excess kurtosis z = +8.3 and
a t=5 tail ratio of 6.3, against the paper's z = -0.3 and 0.878.

The paper never defines "cell" beyond "indexed by depth d, the number of
3,5,7,11,13 dividing N" (sec:floor) and "indexed by which small primes
divide N" (conj:wall). This script asks whether a finer or wider index
closes the gap -- in particular whether the mask depends on the
VALUATIONS v_q rather than on divisibility alone, which the per-depth sd
table already hints at (depth 4 has sd 1.27x the pooled sd after its mean
is removed, so something inside that cell is still structured).

PRE-REGISTRATION.  Decision rule: for each candidate index, remove cell
means alone (never cell sds -- that is a different operation and is
reported separately), and record excess kurtosis, its z, and the t=3,4,5
tail ratios.
  RESCUED   : some candidate index returns z within +-2 AND all three
              tail ratios within their Poisson spread of the paper's
              0.999 / 0.997 / 0.878.
  NOT RESCUED : none does.

Prediction: PARTIALLY RESCUED. I expect valuation-indexed cells to cut
the kurtosis substantially, because the mask plainly varies inside a
divisibility cell, but I do NOT expect the t=5 tail to fall below 1,
because the extreme values sit at depth >= 3 where no mean-removal can
change a spread. If that holds, the operative finding survives: the
tail claim of item 3 is not reproducible by removing means under any
index.
"""

import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np
from math import erfc, sqrt


def phi(z):
    return 0.5 * erfc(-z / sqrt(2.0))


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CACHE = os.path.join(ROOT, "v1_verify2_log", "cache")


def valuation(n, p, cap):
    v = np.zeros(n.shape, dtype=np.int64)
    cur = n.copy()
    for _ in range(cap):
        d = (cur % p) == 0
        if not d.any():
            break
        v[d] += 1
        cur[d] //= p
    return np.minimum(v, cap)


def stats(g, n, label):
    gs = g / g.std(ddof=1)
    ek = float((gs ** 4).mean() - 3.0)
    zk = ek / np.sqrt(24.0 / n)
    row = [label, ek, zk]
    for t in (3, 4, 5):
        obs = int((np.abs(gs) > t).sum())
        exp = 2 * phi(-t) * n
        row.append(obs / exp)
        row.append(obs)
    return row


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 16_000_000
    z = np.load(os.path.join(CACHE, f"field_{X}.npz"))
    good = (z["V"] > 0) & (z["W"] > 0)
    Ni = z["N"][good]
    Nf = Ni.astype(np.float64)
    C = z["C"][good]
    V = z["V"][good]

    # the field established by prop:V
    fld = Nf > 1e5
    Ni = Ni[fld]
    Z = (C[fld] / np.sqrt(V[fld]))
    n = len(Z)
    print("audit_cell_definition   (v1_verify2 Phase 1, blind)")
    print(f"field: 1e5 < N <= {X:,}   ->  {n:,} values")
    print("       (the field on which prop:V's three figures reproduce)")
    print("=" * 78)
    print()

    # ---- candidate cell indices ------------------------------------
    cands = []

    div5 = np.zeros(n, dtype=np.int64)
    for i, p in enumerate((3, 5, 7, 11, 13)):
        div5 |= ((Ni % p == 0).astype(np.int64) << i)
    cands.append(("divisibility by {3,5,7,11,13}  (32)", div5))

    div7 = div5.copy()
    for i, p in enumerate((17, 19)):
        div7 |= ((Ni % p == 0).astype(np.int64) << (5 + i))
    cands.append(("divisibility by {3..19}       (128)", div7))

    # valuations, capped
    for cap in (2, 3):
        lab = np.zeros(n, dtype=np.int64)
        base = 1
        for p in (3, 5, 7, 11, 13):
            lab += valuation(Ni, p, cap) * base
            base *= (cap + 1)
        cands.append((f"valuations v_q<= {cap} of "
                      f"{{3,5,7,11,13}} ({base})", lab))

    # include v_2 as well
    for cap in (2, 3):
        lab = np.zeros(n, dtype=np.int64)
        base = 1
        for p in (2, 3, 5, 7, 11, 13):
            lab += valuation(Ni, p, cap) * base
            base *= (cap + 1)
        cands.append((f"valuations v_q<= {cap} of "
                      f"{{2,3,5,7,11,13}} ({base})", lab))

    # conj:wall's preamble says the mask is "removed by finite modular
    # enumeration", which reads as residue classes rather than
    # divisibility patterns. Both are tried.
    for Q in (2310, 30030):
        cands.append((f"residue class N mod {Q} "
                      f"({Q // 2})", (Ni % Q).astype(np.int64)))

    print(f"  {'cell index':<40}{'cells':>7}{'exc kurt':>11}{'z':>8}"
          f"{'t=3':>8}{'t=4':>8}{'t=5':>8}{'n(t=5)':>8}")
    print("  " + "-" * 96)
    for label, lab in cands:
        u, inv = np.unique(lab, return_inverse=True)
        g = Z.copy()
        sums = np.bincount(inv, weights=Z, minlength=len(u))
        cnts = np.bincount(inv, minlength=len(u))
        means = sums / np.maximum(cnts, 1)
        g = g - means[inv]
        r = stats(g, n, label)
        print(f"  {label:<40}{len(u):>7}{r[1]:>+11.4f}{r[2]:>+8.1f}"
              f"{r[3]:>8.3f}{r[5]:>8.3f}{r[7]:>8.3f}{r[8]:>8}")

    print()
    print("  paper (item 1 and item 3):        "
          f"{'':>7}{-0.0005:>+11.4f}{-0.3:>+8.1f}"
          f"{0.999:>8.3f}{0.997:>8.3f}{0.878:>8.3f}")
    print()

    # ---- the operation the paper does NOT describe -------------------
    print("  for contrast, the same indices with per-cell STANDARDISATION")
    print("  (means AND spreads removed) -- not what the paper describes:")
    print(f"  {'cell index':<40}{'cells':>7}{'exc kurt':>11}{'z':>8}"
          f"{'t=3':>8}{'t=4':>8}{'t=5':>8}{'n(t=5)':>8}")
    print("  " + "-" * 96)
    for label, lab in cands:
        u, inv = np.unique(lab, return_inverse=True)
        g = Z.copy()
        cnts = np.bincount(inv, minlength=len(u))
        sums = np.bincount(inv, weights=Z, minlength=len(u))
        means = sums / np.maximum(cnts, 1)
        g = g - means[inv]
        sq = np.bincount(inv, weights=g * g, minlength=len(u))
        sds = np.sqrt(sq / np.maximum(cnts - 1, 1))
        sds[sds == 0] = 1.0
        ok = cnts[inv] > 5
        g = np.where(ok, g / sds[inv], g)
        r = stats(g, n, label)
        print(f"  {label:<40}{len(u):>7}{r[1]:>+11.4f}{r[2]:>+8.1f}"
              f"{r[3]:>8.3f}{r[5]:>8.3f}{r[7]:>8.3f}{r[8]:>8}")
    print()
    exp5 = 2 * phi(-5) * n
    print(f"  expected count at t=5 on this field: {exp5:.1f}, so the "
          f"Poisson spread")
    print(f"  of any t=5 ratio is about +-{1 / np.sqrt(exp5):.0%}. The "
          f"paper's 0.878 at t=5 is")
    print(f"  {exp5 * 0.878:.1f} events against {exp5:.1f} expected -- "
          f"indistinguishable from 1.000.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
