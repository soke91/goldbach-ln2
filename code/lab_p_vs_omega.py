# -*- coding: utf-8 -*-
"""
The p-decomposition and the omega-decomposition are the same sum
(increment 272).

Proposition P.7 writes, for the p that carry most of the mass,

    D_p = -Lambda(N-p)/log p + Sum_{r<=q<=(N-1)/p} Lambda(N-qp)/log(qp),

so that

    Sum_p log p D_p = -Sum_p Lambda(N-p)
                      + Sum_{p,q} log p Lambda(N-qp)/log(pq) + ...

The first sum is the omega = 1 class of the wall and the second is the
omega = 2 class. That is session 8's alternating series in omega(v),
reached from the other direction: the p-grouping and the
omega-grouping are the SAME sum grouped two ways, and at the p that
carry the mass each D_p has one or two terms, so the p-grouping is
nearly trivial there.

WHAT IS CHECKED.
 (A) The identity, exactly: Sum_p log p D_p against Sum_j (-1)^j T_j
     with T_j = Sum_{omega(v)=j} mu^2(v) Lambda(N-v), both computed
     directly. They must agree to rounding.
 (B) The sizes: T_1 and T_2 separately, their difference, and the
     final C(N), so the depth of the cancellation is visible.
 (C) The share of Sum_p log p |D_p| contributed by p whose D_p has one
     or two terms -- the range where P.7 makes the p-grouping trivial.

There is no null here: (A) is an identity and its check is arithmetic.
(B) and (C) are descriptive, and are reported as such.
"""
import numpy as np
import math
import time

JMAX = 9


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8); mu[1] = 1
    om = np.zeros(X + 1, dtype=np.int8)
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
        om[i] = om[j] + (0 if j % p == 0 else 1)
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p); lg = math.log(int(p))
        while q <= X:
            lam[q] = lg; q *= int(p)
    del spf
    return mu, om, lam, primes


def main():
    X = 900_000
    t0 = time.time()
    mu, om, lam, primes = sieve(X)
    ps = primes[primes < X]
    lp = np.log(ps.astype(np.float64))
    print(f"sieve  t={time.time()-t0:.0f}s", flush=True)

    for N in (30030 * 20, 900000 - 2):
        if N > X:
            continue
        v = np.arange(1, N, dtype=np.int64)
        base = lam[N - v] / np.log(np.maximum(v, 2).astype(np.float64))
        mv = mu[1:N].astype(np.float64)
        t = mv * base; t[0] = 0.0
        a = np.abs(mv) * base; a[0] = 0.0
        jj = om[1:N]

        Sp = 0.0; Sabs = 0.0; few = 0.0
        for i in range(len(ps)):
            p = int(ps[i])
            if p >= N:
                break
            sl = slice(p - 1, None, p)
            d = float(t[sl].sum())
            Sp += lp[i] * d
            Sabs += lp[i] * abs(d)
            nt = int(np.count_nonzero(t[sl]))
            if nt <= 2:
                few += lp[i] * abs(d)

        # Sum_p log p D_p = Sum_v t_v * log v = Sum_v mu(v) Lambda(N-v),
        # i.e. WITHOUT the 1/log v. The omega classes must be built on
        # the same weight or the comparison is of two different objects
        # -- which the first version of this script did.
        aw = np.abs(mv) * lam[N - v]; aw[0] = 0.0
        tw = mv * lam[N - v]; tw[0] = 0.0
        T = np.array([float(aw[jj == j].sum()) for j in range(JMAX + 1)])
        alt = float(sum((-1) ** j * T[j] for j in range(1, JMAX + 1)))
        Cw = float(tw.sum())

        print(f"\n=== N = {N}   (rad-depth "
              f"{len([q for q in (3,5,7,11,13,17) if N % q == 0])})")
        print(f"(A) identity")
        print(f"    Sum_p log p D_p      = {Sp:+.6f}")
        print(f"    Sum_j (-1)^j T_j     = {alt:+.6f}")
        # computed independently as Sum_{v>=2} mu(v) Lambda(N-v), not
        # rearranged from either of the two above; the first version
        # printed an expression that reduced to `alt` algebraically and
        # so could not fail.
        print(f"    Sum_{{v>=2}} mu Lam    = {Cw:+.6f}")
        print(f"    |difference|         = {abs(Sp - alt):.3e}")
        print(f"(B) sizes")
        for j in range(1, 5):
            print(f"    T_{j} = {T[j]:>14.1f}")
        print(f"    T_1 - T_2 + T_3 - T_4 ... = {alt:+.1f}"
              f"    ({abs(alt)/max(T[1],1):.2e} of T_1)")
        print(f"(C) share of Sum_p log p |D_p| from p with <= 2 terms")
        print(f"    {few/Sabs:.4f}   (total {Sabs:.1f})")
    print("\nDONE")


if __name__ == "__main__":
    main()
