# -*- coding: utf-8 -*-
"""
C-III item (2): where does the Heath-Brown mass actually sit?
(increment 205)

C3_DRAFT §2 classifies the one-sided opening into exactly two leaves,
P-I and the Central Object, and justifies the CO modulus bound by

    "modulus q = ak with a <= y^{O(1)}, hence q <= x^{1/3 + O(1/J)}".

The review's coordinate 2 objected that the Heath-Brown identity allows
the mu-side a up to y^J ~ x, so a third leaf exists. That is a
quantitative claim about where the identity's mass sits, and it is
measured here rather than argued.

SETUP. Heath-Brown with J blocks and cut z = M^{1/J}:

    mu(m) = Sum_{j=1..J} (-1)^{j-1} C(J,j)
              Sum_{a_1..a_j b_1..b_j = m, a_i <= z} mu(a_1)...mu(a_j).

Write a = a_1...a_j (the mu-side, each block <= z) and b = b_1...b_j
(the free side). Bookkeeping must bound each piece, so the relevant
weight is the ABSOLUTE one:

    W(a) = Sum_j C(J,j) * A_j(a) * D_j(M/a),
    A_j(a) = #{(a_1..a_j): prod = a, a_i <= z, each squarefree},
    D_j(x) = Sum_{b<=x} d_j(b)   (the free side's ordered factorisations).

The question: what fraction of Sum_a W(a) has a <= M^{1/3}, the region
in which the draft's modulus bound q = ak <= x^{1/3+O(1/J)} holds?

PRE-REGISTERED (fixed before the run):
  DRAFT TENABLE  iff the weight fraction with a > M^{1/3} is small
                 (<= 5%): the CO corner then carries essentially all
                 the mass and the modulus bound is effectively valid.
  THIRD LEAF REAL iff that fraction is substantial (>= 25%): the
                 identity puts real mass at mu-side sizes the draft's
                 bookkeeping excluded, and those pieces have a rough
                 coefficient with no divisor structure, hence no
                 Voronoi entry.
"""
import numpy as np
import math


def mobius_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool)
    pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p * p::p] = False
            mu[p::p] *= -1
            mu[p * p::p * p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
    mu[val > 1] *= -1
    mu[0] = 0
    return mu


def dirichlet_mul(f, g, M):
    """(f*g)(n) for n <= M, f,g given as arrays indexed 0..M."""
    out = np.zeros(M + 1)
    for d in range(1, M + 1):
        if f[d] == 0.0:
            continue
        out[d::d] += f[d] * g[1: M // d + 1]
    return out


def main():
    import sys
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
    Js = [int(v) for v in sys.argv[2:]] or [3, 4, 5, 6, 8]
    print(f"M = {M}")
    mu = mobius_upto(M)

    for J in Js:
        z = int(M ** (1.0 / J))
        # base of the mu-side: |mu(a)| on a <= z, else 0
        base = np.zeros(M + 1)
        base[1:z + 1] = (mu[1:z + 1] != 0).astype(np.float64)
        # base of the free side: all-ones
        one = np.zeros(M + 1)
        one[1:] = 1.0

        # d_j prefix sums, and A_j, built by repeated convolution
        A = base.copy()          # A_1
        d = one.copy()           # d_1 = 1
        totals = {}
        for j in range(1, J + 1):
            if j > 1:
                A = dirichlet_mul(A, base, M)
                d = dirichlet_mul(d, one, M)
            D = np.cumsum(d)     # D_j(x) = sum_{b<=x} d_j(b)
            c = math.comb(J, j)
            idx = np.nonzero(A)[0]
            idx = idx[idx >= 1]
            w = c * A[idx] * D[(M // idx).astype(np.int64)]
            totals[j] = (idx, w)

        # The draft needs a <= y^{O(1)} = x^{O(1/J)}, i.e. a = M^{o(1)}.
        # Report the whole tail profile rather than one threshold.
        thetas = [0.05, 0.10, 0.20, 1.0 / 3]
        tot = sum(float(w.sum()) for _, w in totals.values())
        line = []
        for th in thetas:
            cut = M ** th
            big = sum(float(w[idx > cut].sum())
                      for idx, w in totals.values())
            line.append(big / tot)
        print(f"\nJ = {J},  z = M^(1/J) = {z},  max a = z^J = {z**J} "
              f"(vs M = {M})")
        print("  weight fraction with a > M^theta:")
        print("   " + "  ".join(f"th={th:.2f}: {f:.4f}"
                                for th, f in zip(thetas, line)))
        # where the mass sits in j
        shares = {j: float(w.sum()) / tot for j, (idx, w) in totals.items()}
        print("   mass by block count j: "
              + "  ".join(f"j={j}: {s:.3f}" for j, s in sorted(shares.items())))
    print("\nReading: the draft's CO bookkeeping needs a = M^{o(1)}, so the")
    print("th=0.05 and th=0.10 columns are the relevant ones. Compare the")
    print("same J across the two M to separate genuine J-decay from the")
    print("finite-size artefact that z = M^(1/J) becomes tiny.")
    print("DONE")


if __name__ == "__main__":
    main()
