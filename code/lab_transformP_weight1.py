# -*- coding: utf-8 -*-
"""
Transform Lab, session 9 (increment 229): transform P stated and
measured on C(N) itself, with no log weight.

Sessions 6-8 worked with C_log(N) = Sum_v mu(v) Lambda(N-v) log v,
because log v = Sum_{p|v} log p is what makes the split exact. That
detour has a cost: recovering C(N) from C_log(N) by partial summation
loses a factor, since the trivial bound on the partial sums contributes
O(N log log N / log N) and caps any conclusion at C(N) = O(N loglog N /
log N). Qualitatively enough for the wall, quantitatively not.

The detour is avoidable. For squarefree v >= 2,

    Sum_{p | v} log p / log v = 1     exactly,

so inserting the SAME identity divided by log v gives, with no weight
at all,

    C(N) - Lambda(N-1) = Sum_{p<N} log p * Dt_p(N),
    Dt_p(N) := Sum_{v<N, p | v} mu(v) Lambda(N-v) / log v.

(The v = 1 term is mu(1) Lambda(N-1) = Lambda(N-1) = O(log N), which no
p can reach and which is written out rather than absorbed.) This is the
form worth formalising, and it is the form this script measures, so
that what is proved and what is measured are the same object.

MEASUREMENTS, as in session 6 but on C(N):
 (A) exactness of the identity;
 (B) losslessness: Sum_p log p Mt_p = Sum_{v>=2} mu^2(v) Lambda(N-v),
     which is the trivial bound for C(N) -- so the split redistributes
     the trivial bound and does not enlarge it;
 (C) the absolute aggregate St_abs = Sum_p log p |Dt_p|, which uses no
     cancellation across p, against the trivial bound and against the
     null.

NULL, from the data. Dt_p is a signed sum of t_v = mu(v) Lambda(N-v) /
log v over v divisible by p, so its scale is the square root of its own
second moment Vt_p = Sum_{v<N, p|v} mu^2(v) (Lambda(N-v)/log v)^2,
computed in the same loop (hazard 4, CLOSURE_REAUDIT).

WHAT THE WALL NEEDS HERE. C(N) = o(N) with the trivial bound
Sum_p log p Mt_p ~ A(N) N, so the demand is St_abs = o(N): the same
question as in session 6, now asked of the object the chain actually
consumes.
"""
import numpy as np
import math

from lab_prime_factor_split import sieve


def split1(N, mu, lam, primes):
    v = np.arange(1, N)
    muv = mu[1:N].astype(np.float64)
    lamr = lam[N - v]
    logv = np.log(v.astype(np.float64))
    logv[0] = 1.0                      # v = 1 handled separately
    t = muv * lamr / logv
    t[0] = 0.0
    at = np.abs(muv) * lamr / logv
    at[0] = 0.0
    C = float(np.dot(muv, lamr))
    C_head = float(lamr[0])            # the v = 1 term, Lambda(N-1)
    triv = float((np.abs(muv) * lamr)[1:].sum())

    ps = primes[primes < N]
    D = np.empty(len(ps)); M = np.empty(len(ps)); V = np.empty(len(ps))
    for i, p in enumerate(ps):
        idx = np.arange(int(p), N, int(p)) - 1        # into t[]
        tt = t[idx]
        D[i] = tt.sum()
        M[i] = at[idx].sum()
        V[i] = np.dot(tt, tt)
    return C, C_head, triv, ps, D, M, V


def main():
    X = 400_000
    mu, lam, spf, primes = sieve(X)
    NS = (50_000, 100_000, 200_000, 400_000)

    print("(A) the identity  C(N) - Lambda(N-1) = Sum_p log p Dt_p(N)")
    print(f"{'N':>8} {'C(N)':>12} {'Lam(N-1)':>10} {'Sum lp*Dt':>13} "
          f"{'abs.diff':>10}")
    rows = []
    for N in NS:
        C, C_head, triv, ps, D, M, V = split1(N, mu, lam, primes)
        lp = np.log(ps.astype(np.float64))
        lhs = float(np.dot(lp, D))
        print(f"{N:>8} {C:>12.4f} {C_head:>10.4f} {lhs:>13.4f} "
              f"{abs(lhs - (C - C_head)):>10.2e}")
        rows.append((N, C, C_head, triv, ps, D, M, V, lp))
    print("    exact: Sum_{p|v} log p / log v = 1 on squarefree v >= 2,")
    print("    and no p reaches v = 1, so that term stands outside")

    print("\n(B) losslessness -- Sum_p log p Mt_p vs the trivial bound")
    print(f"{'N':>8} {'Sum lp*Mt_p':>14} {'triv = Sum mu^2 Lam':>20} "
          f"{'rel.diff':>10} {'triv/N':>8}")
    for (N, C, C_head, triv, ps, D, M, V, lp) in rows:
        s = float(np.dot(lp, M))
        print(f"{N:>8} {s:>14.2f} {triv:>20.2f} "
              f"{abs(s-triv)/triv:>10.2e} {triv/N:>8.4f}")
    print("    triv/N sits at A(N) = 0.787, the Lambda-weighted")
    print("    squarefree density -- the same constant Huang-Li's")
    print("    Theorem 1 carries (DEPENDENCY_AUDIT, third pass)")

    print("\n(C) the absolute aggregate -- no cancellation across p")
    print(f"{'N':>8} {'St_abs':>11} {'St_null':>11} {'triv':>11} "
          f"{'abs/triv':>9} {'null/triv':>10} {'abs/null':>9} "
          f"{'abs/N':>8}")
    for (N, C, C_head, triv, ps, D, M, V, lp) in rows:
        S_abs = float(np.dot(lp, np.abs(D)))
        S_null = float(np.dot(lp, np.sqrt(V)))
        print(f"{N:>8} {S_abs:>11.1f} {S_null:>11.1f} {triv:>11.1f} "
              f"{S_abs/triv:>9.4f} {S_null/triv:>10.4f} "
              f"{S_abs/S_null:>9.4f} {S_abs/N:>8.4f}")
    print("    the wall needs St_abs = o(N); abs/triv is the fraction")
    print("    of the trivial bound that survives with no p-cancellation")

    print("\n(D) dyadic profile at N = 400000")
    (N, C, C_head, triv, ps, D, M, V, lp) = rows[-1]
    print(f"{'p range':>16} {'#p':>7} {'mass frac':>10} {'abs frac':>9} "
          f"{'mean rho':>9} {'null rho':>9} {'ratio':>7}")
    b = 2
    while b < N:
        sel = (ps >= b) & (ps < min(2 * b, N))
        if sel.sum():
            m = M[sel]; d = D[sel]; l = lp[sel]; vv = V[sel]
            wm = float(np.dot(l, m))
            wa = float(np.dot(l, np.abs(d)))
            safe = np.maximum(m, 1e-12)
            rho = float(np.dot(l * m, np.abs(d) / safe)) / max(wm, 1e-12)
            rn = float(np.dot(l * m, np.sqrt(vv) / safe)) / max(wm, 1e-12)
            print(f"{b:>7}-{min(2*b, N):>8} {int(sel.sum()):>7} "
                  f"{wm/triv:>10.4f} {wa/triv:>9.4f} {rho:>9.4f} "
                  f"{rn:>9.4f} {rho/max(rn,1e-12):>7.3f}")
        b *= 2
    print("DONE")


if __name__ == "__main__":
    main()
