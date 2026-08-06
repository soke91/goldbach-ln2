# -*- coding: utf-8 -*-
"""
The wall's second moment unfolds into shift-averaged Chowla
(increment 216).

The sweep gave C(N) = sqrt(S(N) N) x unit Gaussian. That is a
statement about the VARIANCE of C over N, and the variance unfolds:

  Sum_{N in W} |C(N)|^2
     = Sum_{n,n'} Lambda(n) Lambda(n') Sum_{N in W, N>max(n,n')}
           mu(N-n) mu(N-n')
     = Sum_h r_W(h) S_W(h),
  r_W(h) = Sum_{n-n'=h} Lambda(n)Lambda(n'),
  S_W(h) = Sum_v mu(v) mu(v+h) over the induced v-range.

So the variance of the wall's scalar is a PRIME-PAIR-WEIGHTED
SHIFT-AVERAGED CHOWLA SUM -- and shift-averaged Chowla is exactly what
Matomaki-Radziwill-Tao control. This script (a) verifies the identity
by brute force, and (b) prices the gap between what the identity needs
and what MRT delivers.

Measured outcome, including a prediction of ours that failed:
  - the identity holds to machine precision (2e-8, 1.2e-7);
  - the off-diagonal is POSITIVE and comparable to the diagonal
    (0.545 and 0.475 of the total at X = 2000, 4000). The docstring
    originally predicted it would have to CANCEL a log X; it does the
    opposite. So roughly half the variance of the wall's scalar is
    genuine shifted-Mobius correlation, not diagonal mass;
  - MRT gives Sum_{|h|<=H} |S(h)| = o(HX), an average saving of o(1)
    per shift over the trivial |S(h)| <= X, while pinning the identity
    to the precision a fixed-N statement needs would require order
    1/X per shift. The shortfall is a factor of order X, not a log
    power -- and a second moment gives almost-all-N in any case, which
    is the known exceptional-set result.
"""
import numpy as np
import math


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8); mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1)
    for p in primes:
        q = int(p); lp = math.log(int(p))
        while q <= X:
            lam[q] = lp; q *= int(p)
    return mu, lam


def main():
    for X, W0 in ((2000, 1000), (4000, 2000)):
        mu, lam = sieve(X)
        Ws = list(range(W0, X + 1, 2))
        # LHS
        lhs = 0.0
        for N in Ws:
            idx = np.arange(1, N)
            lhs += float(np.dot(lam[1:N],
                                mu[N - idx].astype(np.float64))) ** 2
        # RHS by the unfolded double sum
        rhs = 0.0
        diag = 0.0
        ns = np.nonzero(lam[:X + 1])[0]
        for n in ns:
            for np_ in ns:
                lo = max(int(n), int(np_)) + 1
                s = 0.0
                for N in Ws:
                    if N <= lo - 1:
                        continue
                    s += float(mu[N - int(n)]) * float(mu[N - int(np_)])
                t = lam[n] * lam[np_] * s
                rhs += t
                if n == np_:
                    diag += t
        print(f"X = {X}, |W| = {len(Ws)}")
        print(f"  LHS  Sum_N |C(N)|^2      = {lhs:.4f}")
        print(f"  RHS  Sum_h r(h) S(h)     = {rhs:.4f}")
        print(f"  difference               = {abs(lhs-rhs):.3e}"
              f"   {'IDENTITY OK' if abs(lhs-rhs) < 1e-6*max(1,abs(lhs)) else 'MISMATCH'}")
        print(f"  diagonal h = 0           = {diag:.4f}"
              f"   ({diag/lhs:.3f} of total)")
        print(f"  off-diagonal             = {rhs-diag:.4f}"
              f"   ({(rhs-diag)/lhs:+.3f} of total)\n")

    print("Gap: MRT gives an average saving of o(1) per shift; the")
    print("identity needs order 1/X per shift. Shortfall ~ X, not a")
    print("log power. The second-moment route lands on the known")
    print("exceptional-set result, not on a fixed-N statement.")
    print("DONE")


if __name__ == "__main__":
    main()
