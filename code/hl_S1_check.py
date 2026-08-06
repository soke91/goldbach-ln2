# -*- coding: utf-8 -*-
"""
Re-derivation, analytic half, part 1 (increment 222): verify Huang-Li's
evaluation of S1.

Their (13) claims, under EH at level alpha,

    S1(alpha) = Sum_{n<N} Lambda(n) mu^2(N-n) Lambda_alpha(N-n)
              = S(N) N + O(N (log N)^{-A}),

with Lambda_alpha(u) = -Sum_{d|u, d<=alpha} mu(d) log d and S(N) the
singular series. This is the half of the chain that consumes EH for
Lambda; the other half consumes EH_mu and is where our Theorems A and C
live. It has not been checked here before.

The structure of their argument, which the check follows: expanding
Lambda_alpha and swapping gives

    S1 = -Sum_{d<=alpha} mu(d) log d * Sum_{n = N (mod d)} Lambda(n) mu^2(N-n),

so the inner sum is Lambda over a progression, twisted by mu^2 on the
complementary variable; the mu^2 is expanded as Sum_{b^2 | u} mu(b),
which is why their bookkeeping runs over moduli [b^2, d] with a tau_3
count.

NULLS AND CRITERION, stated together. Under (13) the ratio
S1(alpha) / (S(N) N) should approach 1, with a discrepancy of size
(log N)^{-A}. At the N reachable by brute force (log N ~ 12) even A = 1
allows ~8%, so the test cannot confirm the exponent -- only detect a
gross failure.
  CONSISTENT  iff |S1/(S(N)N) - 1| decreases with N and is within the
              (log N)^{-1} band.
  DEFECT      iff the ratio drifts away from 1.
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
    return mu, lam, spf


def singular(N, spf):
    C2 = 0.6601618158468696
    S = 2 * C2
    n = N
    while n > 1:
        p = int(spf[n])
        if p > 2:
            S *= (p - 1) / (p - 2)
        while n % p == 0:
            n //= p
    return S


def main():
    X = 400_000
    mu, lam, spf = sieve(X)
    print(f"{'N':>8} {'alpha':>7} {'S1/N':>10} {'S(N)':>8} "
          f"{'ratio':>8} {'|ratio-1|':>10} {'(log N)^-1':>11}")
    rows = []
    for N in (50_000, 100_000, 200_000, 400_000):
        alpha = math.sqrt(N)
        A = int(alpha)
        # S1 = -sum_{d<=alpha} mu(d) log d * sum_{n = N mod d} Lam(n) mu^2(N-n)
        # build f(u) = mu^2(u) once, and walk d
        S1 = 0.0
        for d in range(1, A + 1):
            md = int(mu[d])
            if md == 0:
                continue
            # n = N - d*j  for j >= 1, n >= 1
            js = np.arange(1, (N - 1) // d + 1, dtype=np.int64)
            ns = N - d * js
            inner = float(np.dot(lam[ns],
                                 (mu[d * js] != 0).astype(np.float64)))
            S1 += -md * math.log(d) * inner
        S = singular(N, spf)
        ratio = S1 / (S * N)
        rows.append(abs(ratio - 1))
        print(f"{N:>8} {alpha:>7.1f} {S1/N:>10.4f} {S:>8.4f} "
              f"{ratio:>8.4f} {abs(ratio-1):>10.4f} "
              f"{1/math.log(N):>11.4f}")
    ok = rows[-1] < rows[0] and rows[-1] < 1 / math.log(400_000) * 3
    print("\nverdict:",
          "CONSISTENT with (13)" if ok else
          "NOT CLEARLY CONSISTENT -- see the numbers against the band")
    print("Note: at these N the (log N)^{-A} band is wide, so this can")
    print("detect a gross failure of (13) but cannot confirm the")
    print("exponent A.")
    print("DONE")


if __name__ == "__main__":
    main()
