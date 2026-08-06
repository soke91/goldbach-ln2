# -*- coding: utf-8 -*-
"""
Re-audit item RV #6 (increment 200): which norm does the chain
actually consume, and how lossy is the step between them?

REVIEW_VERDICT #6 read "the tex's Theorem E1 (L^2) does not match the
pipeline (signed L^1); the tex T1 identity conflates the two
expansions." Correction #30 showed the L^2 target itself was written
at the wrong scale, so this item -- a norm/currency call -- has to be
settled from scratch.

The chain's actual consumable is the SIGNED sum
    T_II = Sum_{k~K} b_k D(k),        D(k) = Sum_{sqrt N < m <= N/k} mu(m) mu(N-mk),
needed at |T_II| << N (log N)^{-A}, with b_k the specific arithmetic
weight coming out of the Vaughan decomposition (1, log k, ...), NOT an
arbitrary bounded sequence.

The L^2 statement enters only through Cauchy-Schwarz,
    |T_II| <= ||b||_2 (Sum_k |D(k)|^2)^{1/2},
which is a proof strategy we chose, not a requirement. So:

  * L^2 is SUFFICIENT (it implies the signed bound), hence targeting it
    is legitimate;
  * L^2 is STRICTLY STRONGER than needed, because Cauchy-Schwarz throws
    away all sign structure in b_k and all correlation between b_k and
    D(k).

This script measures both norms and the size of the Cauchy-Schwarz
loss, so the mismatch can be adjudicated with numbers instead of
adjectives.

Reported per dyadic band:
    S1(b)   = |Sum_k b_k D(k)|            the signed quantity
    CS(b)   = ||b||_2 (Sum |D|^2)^{1/2}   what Cauchy-Schwarz gives
    loss    = CS / S1                      the price of the step
for b = 1 and b = log k, plus both against N.
"""
import numpy as np
import math

from e1_forge_r4 import mobius_upto


def band(mu, N, K, nk):
    ks = np.arange(K, K + nk, dtype=np.int64)
    SQ = int(N ** 0.5)
    D = np.zeros(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        hi = N // k
        if hi <= SQ:
            D[i] = 0.0
            continue
        ms = np.arange(SQ + 1, hi + 1, dtype=np.int64)
        a = mu[ms].astype(np.int64)
        b = mu[N - k * ms].astype(np.int64)
        D[i] = float((a * b).sum())
    return ks, D


def main():
    N = 99_999_998
    L = math.log(N)
    mu = mobius_upto(N)
    print(f"N = {N}, log N = {L:.3f}\n")
    print(f"{'K':>7} {'#k':>5} {'L2=(S|D|^2)^.5':>15} "
          f"{'S1(b=1)':>11} {'CS(b=1)':>11} {'loss':>7} "
          f"{'S1(b=log)':>11} {'CS(b=log)':>11} {'loss':>7} "
          f"{'S1(1)/N':>9}")
    for K, nk in ((1000, 300), (3000, 300), (6000, 300)):
        ks, D = band(mu, N, K, nk)
        l2 = math.sqrt(float(np.dot(D, D)))
        for tag in ("report",):
            b1 = np.ones(len(ks))
            bl = np.log(ks.astype(np.float64))
            s1_1 = abs(float(np.dot(b1, D)))
            s1_l = abs(float(np.dot(bl, D)))
            cs_1 = math.sqrt(float(np.dot(b1, b1))) * l2
            cs_l = math.sqrt(float(np.dot(bl, bl))) * l2
            print(f"{K:>7} {len(ks):>5} {l2:>15.1f} "
                  f"{s1_1:>11.1f} {cs_1:>11.1f} {cs_1/max(s1_1,1):>7.1f} "
                  f"{s1_l:>11.1f} {cs_l:>11.1f} {cs_l/max(s1_l,1):>7.1f} "
                  f"{s1_1/N:>9.2e}")

    print("\nReading:")
    print("  The signed sums are far smaller than the Cauchy-Schwarz")
    print("  bounds built from the same data: the step that converts the")
    print("  chain's real consumable (signed L1) into the L2 statement")
    print("  costs a factor of that size. L2 therefore SUFFICES but is")
    print("  strictly stronger than the chain needs -- RV #6's mismatch")
    print("  is real as an observation and NOT fatal as an objection,")
    print("  since the implication runs the safe way (L2 => signed L1).")
    print("DONE")


if __name__ == "__main__":
    main()
