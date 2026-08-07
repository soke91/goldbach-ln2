# -*- coding: utf-8 -*-
"""
sec:c3 item (2) -- the Heath-Brown weight's mass outside a <= M^{0.05}.
(v1_verify2, Phase 1, blind.)

STATEMENT UNDER TEST, verbatim:

  "Measuring the absolute Heath--Brown weight
   W(a) = sum_j binom(J,j) A_j(a) D_j(M/a) across a-sizes, the fraction
   with a > M^{0.05} is 0.939, 0.949, 0.960, 0.947 at J = 3, 4, 6, 8,
   rising to 0.961 and 0.969 at M = 10^6. So about 95% of the identity's
   weight lies outside the region the classification covers ...
   the identity with cut z needs z^J >= x, and its j-th term has
   a <= z^j, which at j = J is x for every admissible (z,J) --- while
   the weight concentrates in exactly those high-j terms (0.824 at j=3
   when J=3; 0.847 at j in {6,7,8} when J=8)."

THE READING USED, stated because the paper defines neither A_j nor D_j.
Heath-Brown's identity with cut z and depth J writes the Mobius side as
j factors m_1..m_j each <= z, so with a = m_1...m_j:

  A_j(a) = #{(m_1,..,m_j) : prod m_i = a, m_i <= z}
           = the j-fold divisor function restricted to factors <= z,
           computed as the j-th Dirichlet power of 1_{[1,z]}.
  D_j(y) = #{(n_1,..,n_j) : prod n_i <= y} = sum_{n<=y} d_j(n),
           the summatory j-fold divisor function.
  z = M^{1/J}, the smallest admissible cut (z^J = M).

"Absolute" means the mu-signs are dropped, which is what makes W a
weight rather than a sum. Both quantities are exact integers here.

PRE-REGISTRATION.

  Decision rule. Compute the share of W with a > M^{0.05}, for
  J = 3,4,6,8 at M = 1e4, 1e5, 1e6, and the per-j shares of W.
    REPRODUCED   : the four fractions land near 0.939/0.949/0.960/0.947
                   at one M, and near 0.961/0.969 at M = 1e6.
    NOT REPRODUCED: they do not, under this reading; then report the
                   spread across readings of D_j, since that is the
                   ambiguity.
  Also report, as the paper's own argument requires, the share of the
  weight carried by the top-j terms.

  Prediction written before running.  The QUALITATIVE claim will hold
  overwhelmingly -- M^{0.05} at M=1e6 is 1.995, i.e. the threshold is
  essentially "a > 1", and almost all of the weight has a >= 2 by
  construction. I predict the fractions come out ABOVE the quoted ones,
  closer to 0.98--1.00, and that the interesting content is not the
  fraction but that the threshold M^{0.05} is a number near 2 for every
  M this program can reach: at M=1e6 the classification's covered region
  is a in {1}, so "95% of the weight lies outside" understates it.
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np


def dirichlet_mul(f, g, M):
    """(f * g)(n) for n <= M, f and g supported on [1, M]."""
    out = np.zeros(M + 1)
    nzf = np.nonzero(f)[0]
    for d in nzf:
        d = int(d)
        if d < 1 or d > M:
            continue
        out[d::d] += f[d] * g[1: M // d + 1]
    return out


def main():
    print("c3_hb_weight   (v1_verify2 Phase 1, blind)")
    print("=" * 74)
    print("reading: A_j = j-th Dirichlet power of 1_{[1,z]}, z = M^{1/J};")
    print("         D_j(y) = sum_{n<=y} d_j(n);  W(a) = sum_j C(J,j)"
          " A_j(a) D_j(M/a)")
    print()

    for M in (10_000, 100_000, 1_000_000):
        thr = M ** 0.05
        print(f"--- M = {M:,}   threshold M^0.05 = {thr:.3f}"
              f"   (so 'covered' means a <= {int(thr)}) ---")
        # d_j summatory functions
        one = np.zeros(M + 1)
        one[1:] = 1.0
        dj = [None, one.copy()]
        for j in range(2, 9):
            dj.append(dirichlet_mul(dj[-1], one, M))
        Dj = [None] + [np.cumsum(x)[: M + 1] for x in dj[1:]]

        print(f"  {'J':>4}{'z':>9}{'frac a > M^0.05':>18}"
              f"{'top-j share':>14}{'per-j shares':>34}")
        for J in (3, 4, 6, 8):
            z = M ** (1.0 / J)
            ind = np.zeros(M + 1)
            ind[1: int(z) + 1] = 1.0
            A = [None, ind.copy()]
            for j in range(2, J + 1):
                A.append(dirichlet_mul(A[-1], ind, M))
            a = np.arange(M + 1)
            W = np.zeros(M + 1)
            perj = []
            from math import comb
            for j in range(1, J + 1):
                y = np.zeros(M + 1)
                y[1:] = Dj[j][(M // a[1:]).astype(np.int64)]
                term = comb(J, j) * A[j] * y
                perj.append(term.sum())
                W += term
            tot = W.sum()
            frac = W[a > thr].sum() / tot
            perj = np.array(perj) / tot
            top = perj[-1] if J == 3 else perj[5:].sum() if J == 8 \
                else perj[-1]
            s = " ".join(f"{v:.3f}" for v in perj)
            print(f"  {J:>4}{z:>9.2f}{frac:>18.4f}{top:>14.3f}   {s}")
        print()
    print("[paper: 0.939, 0.949, 0.960, 0.947 at J=3,4,6,8, rising to")
    print(" 0.961 and 0.969 at M=1e6; per-j 0.824 at j=3 when J=3,")
    print(" 0.847 at j in {6,7,8} when J=8]")
    print()
    print("note: M^{0.05} is 1.585 at M=1e4 and 1.995 at M=1e6, so the")
    print("'covered region' a <= M^{0.05} is exactly {a=1} at every M")
    print("reachable here. The fraction outside it is therefore 1 minus")
    print("the share of the single term a=1, and its closeness to 1 is a")
    print("statement about the threshold, not about the identity.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
