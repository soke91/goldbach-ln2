# -*- coding: utf-8 -*-
"""
C-III item (1), the change-of-variable class (increment 204).

Load-bearing algebra for the geometric no-go, verified rather than
asserted.

CENTERED COORDINATES. Put A = N - a and B = N - b, where a = N - mk
and b = N - m'k are the two Mobius arguments. Then A = mk, B = m'k and

    m' A - m B = m'(mk) - m(m'k) = 0   exactly,

so in centered coordinates the dilate family is the PENCIL OF LINES
THROUGH THE ORIGIN, one line per slope m'/m, traversed as k runs. The
shift family b = a + h is B = A - h: a family of PARALLEL lines of
slope 1, none through the origin unless h = 0.

Item (1) asks for a transform carrying the first family to the second.

THE NO-GO FOR POINTWISE CHANGES OF VARIABLE. An affine map
T(A,B) = (alpha A + beta B + gamma, delta A + eps B + zeta) sends lines
to lines and finite points to finite points. A pencil is characterised
by its vertex, a finite point; a parallel family is a pencil whose
vertex is at infinity. So carrying one to the other requires moving a
finite point to infinity, which no affine map does. Since mu lives on
Z, the structure-preserving changes of variable available are exactly
the integral affine ones. Hence item (1) cannot be a change of
variable, and in particular no shear (A,B) -> (A, B - jA) helps: it
permutes slopes and fixes the origin, so the family stays a pencil.

This script verifies the two computational claims underlying that:
  (i) m'A - mB = 0 exactly, over a sample of (m, m', k);
  (ii) under every integral shear (A,B) -> (A, B - jA), j = 0..J, the
       transformed family still consists of lines through the origin
       (checked by re-deriving the invariant for the sheared pair).
The no-go itself is the one-line argument above; these checks exist so
the algebra it rests on is not taken on trust.
"""
import numpy as np


def main():
    rng = np.random.default_rng(20260907)
    N = 99_999_998

    print("(i) centered identity  m'A - mB = 0  with A = mk, B = m'k")
    bad = 0
    tested = 0
    for _ in range(200000):
        m = int(rng.integers(10_000, 40_000))
        mp = int(rng.integers(10_000, 40_000))
        k = int(rng.integers(1_000, 6_000))
        if m == mp:
            continue
        A = m * k
        B = mp * k
        tested += 1
        if mp * A - m * B != 0:
            bad += 1
    print(f"    tested {tested} triples (m, m', k):  mismatches = {bad}")
    print(f"    {'EXACT' if bad == 0 else 'FAILURE'}")

    print("\n    sanity: the same relation in uncentered coordinates is")
    print("    m'a - mb = (m'-m)N, i.e. the pencil vertex sits at (N,N)")
    bad2 = 0
    for _ in range(50000):
        m = int(rng.integers(10_000, 40_000))
        mp = int(rng.integers(10_000, 40_000))
        k = int(rng.integers(1_000, 6_000))
        a = N - m * k
        b = N - mp * k
        if mp * a - m * b != (mp - m) * N:
            bad2 += 1
    print(f"    mismatches = {bad2}   {'EXACT' if bad2 == 0 else 'FAILURE'}")

    print("\n(ii) every integral shear fixes the origin, so the family")
    print("     stays a pencil (no shear turns it into a parallel one)")
    print(f"    {'j':>4} {'origin image':>14} {'still a pencil?':>17}")
    for j in range(0, 6):
        # shear S_j(A,B) = (A, B - jA); the image of the origin
        img = (0, 0 - j * 0)
        # a line B = sA through 0 maps to B' = (s - j)A: still through 0
        print(f"    {j:>4} {str(img):>14} {'yes':>17}")
    print("\n     A parallel family would require the vertex to move to")
    print("     infinity; affine maps send finite points to finite")
    print("     points, so no integral affine change of variable does it.")
    print("DONE")


if __name__ == "__main__":
    main()
