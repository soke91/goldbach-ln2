# -*- coding: utf-8 -*-
"""
Audit of the E1 target itself (increment 198).

AMPLITUDE_ADJUDICATION.md states the consumable as

    (as printed)   Sum_{k~K} |D(k)|^2  <<  (log N)^{-A} Sum_{k~K} M_k,

with D(k) = Sum_{sqrt N < m <= N/k} mu(m) mu(N-mk) and M_k the length
of the m-range. But this repository's own measurements read the ratio
Sum |D|^2 / Sum supp as 0.966 / 0.950 / 0.922 -- i.e. ~1, not
(log N)^{-A}. Both cannot be right. This script re-derives the target
from the chain and measures all three quantities.

Derivation. The wall is the type-II term
    T_II = Sum_{k~K} b_k D(k),   |b_k| << log N,
which the chain needs at  |T_II| << N (log N)^{-A}.  Cauchy-Schwarz in
k gives |T_II| <= ||b||_2 (Sum_k |D(k)|^2)^{1/2}, with
||b||_2^2 ~ K (log N)^2, so the requirement is

    (CORRECT)      Sum_{k~K} |D(k)|^2  <<  (log N)^{-2A-2} Sum_{k~K} M_k^2,

since Sum_{k~K} M_k^2 ~ K (N/K)^2 = N^2/K.  The printed form has M_k
where M_k^2 belongs.  The difference is a factor M ~ N/K, which is a
POWER of N -- so the printed target demands beating square-root
cancellation by a further log power, which the measurements show is
false, while the correct target is a log-power saving over the TRIVIAL
bound, which square-root cancellation clears with room to spare.

This script reports, per dyadic band:
    S2   = Sum_k |D(k)|^2                  (what nature does)
    T1   = Sum_k M_k                       (the printed normaliser)
    T2   = Sum_k M_k^2                      (the correct normaliser)
    need = (log N)^{-2A-2} * T2   at A = 4  (what the chain requires)
    slack = need / S2                       (>1 means nature clears it)

and the asymptotic form of the slack, N / (K (log N)^{2A+2}).
"""
import numpy as np
import math

from e1_forge_r4 import mobius_upto


def band(mu, N, K, nk):
    ks = np.arange(K, K + nk)
    SQ = int(N ** 0.5)
    S2 = 0.0
    T1 = 0.0
    T2 = 0.0
    supp = 0.0
    for k in ks:
        k = int(k)
        hi = N // k
        if hi <= SQ:
            continue
        ms = np.arange(SQ + 1, hi + 1, dtype=np.int64)
        a = mu[ms].astype(np.int64)
        b = mu[N - k * ms].astype(np.int64)
        d = float((a * b).sum())
        S2 += d * d
        Mk = float(ms.size)
        T1 += Mk
        T2 += Mk * Mk
        supp += float(np.count_nonzero(a * b))
    return S2, T1, T2, supp, len(ks)


def main():
    N = 99_999_998
    A = 4
    L = math.log(N)
    mu = mobius_upto(N)
    print(f"N = {N},  log N = {L:.3f},  A = {A}\n")
    print(f"{'K':>7} {'#k':>5} {'S2/T1':>9} {'S2/supp':>9} "
          f"{'S2/T2':>12} {'need/T2':>12} {'SLACK':>12} "
          f"{'asympt':>12}")
    for K, nk in ((1000, 300), (3000, 300), (6000, 300)):
        S2, T1, T2, supp, n = band(mu, N, K, nk)
        if S2 == 0.0:
            print(f"{K:>7} {n:>5}  (band empty: K must be < sqrt N)")
            continue
        need = T2 * L ** (-(2 * A + 2))
        slack = need / S2
        asym = N / (K * L ** (2 * A + 2))
        print(f"{K:>7} {n:>5} {S2/T1:>9.4f} {S2/supp:>9.4f} "
              f"{S2/T2:>12.3e} {need/T2:>12.3e} {slack:>12.3e} "
              f"{asym:>12.3e}")

    print("\nReading:")
    print("  S2/supp ~ 1        -> square-root cancellation, exactly as")
    print("                        Conjecture L predicts (this is what")
    print("                        the repository has been measuring).")
    print("  S2/T1 ~ 1 too      -> so the PRINTED target, S2 << (log N)^-A T1,")
    print("                        is FALSE by five orders of magnitude:")
    print("                        it demands more than nature delivers.")
    print("  SLACK              -> = (N/K) / (0.3 (log N)^{2A+2}), which")
    print("                        tends to infinity since K <= N^{1/3},")
    print("                        but is NOT visible at N = 10^8: there")
    print("                        (log N)^10 = 4.5e12 dwarfs N/K = 1e5.")
    print("                        The correct target is asymptotic; the")
    print("                        printed one is false at every N.")
    print("DONE")


if __name__ == "__main__":
    main()
