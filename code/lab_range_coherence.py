# -*- coding: utf-8 -*-
"""
Transform Lab, session 10 (increment 230): does grouping transform P by
dyadic RANGES of p beat grouping by single p?

Proposition P.4 asks for Sum_p log p |D_p| = o(N) -- absolute values
taken one prime at a time, so no cancellation across p is used. That is
a deliberately weak hypothesis, and the measured margin belongs to it.
But it is not the only sufficient condition available. Grouping the
same identity by dyadic ranges R = [P, 2P) gives

    C(N) - Lambda(N-1) = Sum_R C_R(N),
    C_R(N) := Sum_{p in R} log p D_p(N)
            = Sum_{2<=v<N} mu(v) Lambda(N-v) L_R(v) / log v,
    L_R(v) := Sum_{p | v, p in R} log p,

and Sum_R |C_R| = o(N) is a WEAKER hypothesis than P.4's, since
|C_R| <= Sum_{p in R} log p |D_p| termwise. The question is how much
weaker -- i.e. how much cancellation there is among the D_p inside a
range.

WHY IT MATTERS FOR WHERE WE ARE. Session 6's dyadic profile put all of
the surviving mass at large p, where D_p is a short sum: at p > N/2 it
is a SINGLE term, mu(p) Lambda(N-p)/log p, and mu(p) = -1 always, so
no cancellation exists inside D_p at all. Range grouping is the natural
repair: for p in [P, 2P) with P large, w = v/p runs over a long range
too, so C_R is a genuine bilinear form in (p, w) with both variables
prescribed -- the shape Type-II technology is built for. The design
rule from Proposition P.3 is respected throughout: mu(w) stays inside
the group.

BUT NOTE WHAT RANGE GROUPING CANNOT FIX. For p > N/2 the only v with
p | v is v = p itself, so C_R over the top range is
-Sum_{p in R} Lambda(N-p): every term negative, no cancellation even
after grouping. That part of the demand is irreducible within transform
P, and it has size ~ S(N) N / (2 log N), i.e. a fraction ~1/log N of
the trivial bound. Predicted before running: |C_R| / A_R = 1.0000
exactly in the top range. Anything else is a bug here.

NULLS, from the data (hazard 4, CLOSURE_REAUDIT). For each range,
    null_R := sqrt( Sum_v mu^2(v) (Lambda(N-v) L_R(v) / log v)^2 ),
the exact second moment of the same v-sum. Reported beside
    A_R := Sum_{p in R} log p |D_p|   (P.4's demand on the range)
    T_R := Sum_{p in R} log p M_p     (trivial mass of the range).

CRITERIA, stated with the nulls.
  * |C_R| / A_R near 1  => range grouping adds nothing there.
  * |C_R| / A_R falling with the number of primes in the range
    => the D_p cancel against each other, and the weaker hypothesis is
    materially weaker.
  * The aggregate Sum_R |C_R| against S_abs = Sum_R A_R is the verdict:
    if it is substantially smaller AND still decaying as a fraction of
    the trivial bound, range grouping is the better sufficient
    condition and the one to prove.
"""
import numpy as np
import math

from lab_prime_factor_split import sieve


def main():
    X = 400_000
    mu, lam, spf, primes = sieve(X)
    NS = (50_000, 100_000, 200_000, 400_000)

    agg = []
    for N in NS:
        v = np.arange(1, N)
        muv = mu[1:N].astype(np.float64)
        lamr = lam[N - v]
        logv = np.log(v.astype(np.float64)); logv[0] = 1.0
        base = lamr / logv                      # Lambda(N-v)/log v
        base[0] = 0.0
        triv = float((np.abs(muv) * lamr)[1:].sum())
        ps = primes[primes < N]
        lp = np.log(ps.astype(np.float64))

        rows = []
        b = 2
        while b < N:
            hi = min(2 * b, N)
            sel = (ps >= b) & (ps < hi)
            if sel.sum():
                L = np.zeros(N - 1)
                for p in ps[sel]:
                    L[int(p) - 1::int(p)] += math.log(int(p))
                x = muv * base * L
                ax = np.abs(muv) * base * L
                C_R = float(x.sum())
                T_R = float(ax.sum())
                null_R = float(math.sqrt(np.dot(x, x)))
                A_R = 0.0
                for p in ps[sel]:
                    idx = np.arange(int(p), N, int(p)) - 1
                    A_R += math.log(int(p)) * abs(
                        float((muv[idx] * base[idx]).sum()))
                rows.append((b, hi, int(sel.sum()), C_R, A_R, T_R,
                             null_R))
            b *= 2
        agg.append((N, triv, rows))

    print("(A) dyadic profile at N = 400000 -- coherent vs absolute")
    N, triv, rows = agg[-1]
    print(f"{'p range':>16} {'#p':>6} {'T_R/triv':>9} {'A_R/triv':>9} "
          f"{'|C_R|/triv':>11} {'|C_R|/A_R':>10} {'|C_R|/null':>11}")
    for (b, hi, np_, C_R, A_R, T_R, null_R) in rows:
        print(f"{b:>7}-{hi:>8} {np_:>6} {T_R/triv:>9.4f} "
              f"{A_R/triv:>9.4f} {abs(C_R)/triv:>11.6f} "
              f"{abs(C_R)/max(A_R,1e-12):>10.4f} "
              f"{abs(C_R)/max(null_R,1e-12):>11.4f}")
    print("    A_R is P.4's demand on the range; |C_R| is the weaker")
    print("    demand that range grouping asks instead")

    print("\n(B) the aggregate -- which sufficient condition is easier")
    print(f"{'N':>8} {'Sum|C_R|':>11} {'S_abs':>11} {'triv':>11} "
          f"{'coh/triv':>9} {'abs/triv':>9} {'coh/abs':>8} {'coh/N':>8}")
    for (N, triv, rows) in agg:
        coh = sum(abs(r[3]) for r in rows)
        ab = sum(r[4] for r in rows)
        print(f"{N:>8} {coh:>11.1f} {ab:>11.1f} {triv:>11.1f} "
              f"{coh/triv:>9.4f} {ab/triv:>9.4f} {coh/ab:>8.4f} "
              f"{coh/N:>8.4f}")
    print("    the wall needs either of these to be o(N); coh/abs is")
    print("    how much range grouping buys over P.4")

    print("\n(C) the irreducible part -- ranges above N/2, where v = p")
    print(f"{'N':>8} {'top |C_R|':>11} {'/triv':>9} {'S(N)N/(2 logN)':>15} "
          f"{'|C|/A there':>12}")
    for (N, triv, rows) in agg:
        top = [r for r in rows if r[0] >= N / 2]
        c = sum(abs(r[3]) for r in top)
        a = sum(r[4] for r in top)
        pred = 1.3203 * N / (2 * math.log(N))
        print(f"{N:>8} {c:>11.1f} {c/triv:>9.4f} {pred:>15.1f} "
              f"{c/max(a,1e-12):>12.4f}")
    print("    these v are prime, mu(v) = -1 with no sign variation, so")
    print("    no grouping recovers cancellation; the share falls like")
    print("    1/log N and that is the floor transform P cannot cross")
    print("DONE")


if __name__ == "__main__":
    main()
