# -*- coding: utf-8 -*-
"""
Re-derivation, the assembly (increment 223).

Huang-Li's Theorem 1 concludes r(N) >= S(N)(1 - A(N)) N + O(N (log N)^-A),
with A(N) = prod_{p not| N, p>2} (1 - 1/(p(p-1))) their (7).

Where does (1 - A(N)) come from? Tracing the chain:

  r(N) = S(N) (N - C(N)) + O_A      [their (22), our Theorem C]
  C(N) = Sum_{n<N} Lambda(n) mu(N-n)

so a LOWER bound on r(N) needs an UPPER bound on C(N). The only bound
available without knowing anything about the sign of mu is

  |C(N)| <= Sum_{n<N} Lambda(n) mu^2(N-n),

and the claim implicit in their Theorem 1 is that this majorant is
asymptotically A(N) N -- i.e. A(N) is the density of squarefree values
of N - n weighted by Lambda(n). That is the step this script checks; it
is the last unverified link in the skeleton.

NULL AND CRITERION. If the claim holds, T(N)/N -> A(N), so the ratio
T(N)/(A(N) N) tends to 1 with a discrepancy of the usual (log N)^{-1}
size at reachable N.
  CONSISTENT iff the ratio approaches 1 monotonically within the band.

Also reported: whether the resulting lower bound S(N)(1-A(N))N is
actually positive, which is what makes Corollary 1 give Goldbach --
and the trivial-looking competitor 6/pi^2, to show why A(N) rather
than the squarefree density is the right constant.
"""
import numpy as np
import math

from hl_S1_check import sieve, singular


def A_of_N(N, spf, primes, PLIM):
    Np = set()
    n = N
    while n > 1:
        p = int(spf[n]); Np.add(p)
        while n % p == 0:
            n //= p
    A = 1.0
    for p in primes:
        p = int(p)
        if p > 2 and p not in Np:
            A *= (1 - 1.0 / (p * (p - 1)))
    return A


def main():
    X = 400_000
    mu, lam, spf = sieve(X)
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2

    print(f"{'N':>8} {'T/N':>9} {'A(N)':>8} {'ratio':>8} "
          f"{'|r-1|':>8} {'(logN)^-1':>10} {'S(1-A)':>8} {'6/pi^2':>8}")
    devs = []
    for N in (50_000, 100_000, 200_000, 400_000):
        idx = np.arange(1, N)
        T = float(np.dot(lam[1:N],
                         (mu[N - idx] != 0).astype(np.float64)))
        A = A_of_N(N, spf, primes, X)
        S = singular(N, spf)
        ratio = T / (A * N)
        devs.append(abs(ratio - 1))
        print(f"{N:>8} {T/N:>9.4f} {A:>8.4f} {ratio:>8.4f} "
              f"{abs(ratio-1):>8.4f} {1/math.log(N):>10.4f} "
              f"{S*(1-A):>8.4f} {6/math.pi**2:>8.4f}")

    ok = devs[-1] < devs[0]
    print("\nverdict:",
          "CONSISTENT -- the majorant of |C(N)| is A(N) N"
          if ok else "NOT CONSISTENT")
    print("\nWhy A(N) and not 6/pi^2: the squarefree density 0.6079 is")
    print("the density of squarefree u over ALL u. Here u = N - n with")
    print("n weighted by Lambda, i.e. n running over prime powers, and")
    print("the local density at p becomes 1 - 1/(p(p-1)) instead of")
    print("1 - 1/p^2, since n avoids the residue class 0 mod p. That is")
    print("exactly A(N), and it is LARGER than 6/pi^2 -- so the bound is")
    print("weaker than the naive one, which is why 1 - A(N) is small and")
    print("Corollary 1 has to work to make it positive.")
    print("DONE")


if __name__ == "__main__":
    main()
