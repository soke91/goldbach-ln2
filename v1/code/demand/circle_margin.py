# -*- coding: utf-8 -*-
"""
The circle-method margin for C(N) = Sum_{n<N} Lambda(n) mu(N-n)
(increment 196).

C(N) is the wall: Huang-Li's Theorem 1 says r(N) ~ S(N) N is
EQUIVALENT to C(N) = o(N). The divisor-switch route to it is closed
over its entire weight space (Theorem D). This script measures whether
the OTHER classical mechanism -- the circle method -- has any room at
all, by computing the exact norms its estimates would have to beat.

    C(N) = Integral_0^1 S_Lambda(a) S_mu(-a) e(-N a) da,
    S_Lambda(a) = Sum_{n<=N} Lambda(n) e(na),
    S_mu(a)     = Sum_{m<=N} mu(m) e(ma).

Every circle-method treatment of a binary problem pays one of two
bills:

  (i)  Cauchy-Schwarz:  |C| <= ||S_Lambda||_2 ||S_mu||_2
       = (N log N)^{1/2} (6N/pi^2)^{1/2} ~ 0.78 N (log N)^{1/2},
       which is a factor (log N)^{1/2} ABOVE the trivial bound
       |C| <= psi(N) ~ N.

  (ii) a pointwise bound on one factor against the L^1 norm of the
       other:  |C| <= sup_a |S_mu(a)| * ||S_Lambda||_1.

Route (ii) is where any hope lives, and Parseval already caps it:
sup_a |S_mu| >= ||S_mu||_2 ~ 0.78 N^{1/2} for free, and
||S_Lambda||_1 is of order N^{1/2} up to logs, so the product is of
order N -- the trivial bound, with no margin. The purpose of the
measurement is to see whether the constant works out below 1 and
whether the margin

    MARGIN(N) := N / ( sup_a |S_mu(a)| * ||S_Lambda||_1 )

grows, is flat, or decays with N.

PRE-REGISTERED (fixed before the run):
  ROOM      iff MARGIN(N) grows at least like a power of log N across
            the scales tested (then a pointwise mu-bound at the
            Parseval floor would give C(N) = o(N) and the route is
            worth developing).
  NO ROOM   iff MARGIN(N) is bounded or decaying (then the route is
            closed unconditionally, because sup >= L2 is free and
            cannot be improved).

Also reported: the true C(N)/N, to see where the actual object sits
relative to both bills.
"""
import numpy as np
import math


def sieve_mu_lambda(X):
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
    # Lambda
    lam = np.zeros(X + 1)
    sieve = np.ones(X + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(X ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    for p in np.nonzero(sieve)[0]:
        q = int(p)
        lp = math.log(p)
        while q <= X:
            lam[q] = lp
            q *= int(p)
    return mu, lam


def run(N, over=4):
    """Norms of S_mu and S_Lambda on a grid of Q = over*N points."""
    mu, lam = sieve_mu_lambda(N)
    Q = over * N
    a = np.zeros(Q)
    a[:N + 1] = mu[:N + 1]
    Smu = np.fft.fft(a)
    b = np.zeros(Q)
    b[:N + 1] = lam[:N + 1]
    Slam = np.fft.fft(b)

    absmu = np.abs(Smu)
    abslam = np.abs(Slam)
    sup_mu = float(absmu.max())
    l2_mu = float(np.sqrt(absmu.dot(absmu) / Q))
    l1_lam = float(abslam.sum() / Q)
    l2_lam = float(np.sqrt(abslam.dot(abslam) / Q))
    sup_lam = float(abslam.max())

    # true C(N) for the even N nearest below
    Ne = N if N % 2 == 0 else N - 1
    idx = np.arange(1, Ne)
    C = float(np.dot(lam[1:Ne], mu[Ne - idx].astype(np.float64)))
    return dict(N=N, Ne=Ne, sup_mu=sup_mu, l2_mu=l2_mu, l1_lam=l1_lam,
                l2_lam=l2_lam, sup_lam=sup_lam, C=C)


def main():
    print("Circle-method margin for C(N) = Sum Lambda(n) mu(N-n)")
    print("Parseval floor: sup|S_mu| >= ||S_mu||_2 ~ sqrt(6N/pi^2) "
          "= 0.7797 sqrt(N)\n")
    hdr = (f"{'N':>9} {'sup|Smu|/rtN':>13} {'||Smu||2/rtN':>13} "
           f"{'||Slam||1/rtN':>14} {'CS bound/N':>11} "
           f"{'MARGIN':>8} {'C(N)/N':>9}")
    print(hdr)
    rows = []
    for N in (2 ** 14, 2 ** 16, 2 ** 18, 2 ** 20):
        r = run(N)
        rt = math.sqrt(N)
        cs = r['l2_lam'] * r['l2_mu'] / N
        margin = N / (r['sup_mu'] * r['l1_lam'])
        rows.append((N, margin))
        print(f"{N:>9} {r['sup_mu']/rt:>13.4f} {r['l2_mu']/rt:>13.4f} "
              f"{r['l1_lam']/rt:>14.4f} {cs:>11.4f} "
              f"{margin:>8.4f} {r['C']/r['Ne']:>9.4f}")

    print("\n=== PRE-REGISTERED READING ===")
    for N, m in rows:
        print(f"  N = {N:>9}: MARGIN = {m:.4f}   (log N)^(1/2) = "
              f"{math.log(N)**0.5:.3f}")
    growing = rows[-1][1] > 2 * rows[0][1]
    print("verdict:",
          "ROOM (margin grows)" if growing else
          "NO ROOM -- the binary circle method sits at or above the "
          "trivial bound; sup >= L2 is free and cannot be improved, "
          "so no pointwise Mobius bound can rescue it")
    print("DONE")


if __name__ == "__main__":
    main()
