# -*- coding: utf-8 -*-
"""
Theorem D load-bearing check (increment 195): the weight-space
trade-off on the demand side.

Setting. The divisor switch turns any weighted demand functional into
    B_w * C(N) = [complete part] - [residual],   C(N) = Sum Lambda(n) mu(N-n),
    B_w = Sum_{k<K, (k,N)=1} mu(k) w_k / phi(k).
Write w = 1 * b (i.e. w_k = Sum_{d|k} b_d, b = mu * w -- every weight
has this form). Two exact facts then oppose each other:

 (1) EXTRACTION.  B_w = Sum_{d<K,(d,N)=1} b_d * mu(d)/phi(d) * rho_{dN}(K/d),
     where rho_{dN}(x) = Sum_{j<x, (j,dN)=1} mu(j)/phi(j).
     Huang-Li's Lemma 1 gives |rho| << exp(-c sqrt(log x)), so B_w is
     controlled by the part of b sitting AT the truncation point K:
     mass at d <= K^{1-eps} is damped by exp(-c sqrt(eps log K)).

 (2) BV-ACCESSIBILITY.  The residual expands, via w_k = Sum_{d|k} b_d,
     into Lambda-sums over progressions to moduli m*d with
     m < N^{1-theta'}.  Bombieri-Vinogradov needs m*d <= N^{1/2-delta},
     i.e. b supported on d <= N^{theta'-1/2-delta}.

The two requirements are separated by a factor N^{1/2}: that is the
sqrt(N) barrier, appearing on the demand side as a gap between where a
weight must live to see C(N) and where it may live for BV to close its
residual.

This script verifies the load-bearing quantity of (1) directly:
  - rho(x) is computed exactly and compared with exp(-c sqrt(log x));
  - B_w is computed exactly, by brute force over k < K, for the
    single-divisor weights w_k = [d0 | k] (b = delta_{d0}), and
    checked against the factorised formula;
  - the resulting damping is tabulated against d0, showing that
    |B_w| is of order 1 only when d0 is within exp(O((log log)^2)) of
    K, and is exponentially small for d0 below N^{theta'-1/2}.
"""
import numpy as np
import math


def mobius_upto(X):
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
    return mu


def totient_upto(X):
    phi = np.arange(X + 1, dtype=np.int64)
    for p in range(2, X + 1):
        if phi[p] == p:                      # p prime
            phi[p::p] -= phi[p::p] // p
    return phi


def coprime_mask(X, n):
    """boolean mask of j <= X with gcd(j, n) = 1"""
    m = np.ones(X + 1, dtype=bool)
    q = n
    p = 2
    while p * p <= q:
        if q % p == 0:
            m[p::p] = False
            while q % p == 0:
                q //= p
        p += 1
    if q > 1:
        m[q::q] = False
    return m


def main():
    N = 99_999_998
    theta = 0.56
    K = int(N ** theta)
    print(f"N = {N},  theta' = {theta},  K = N^theta' = {K}")
    print(f"BV admits weight-support only up to "
          f"N^(theta'-1/2) = {int(N**(theta-0.5))}\n")

    mu = mobius_upto(K + 1)
    phi = totient_upto(K + 1)

    # --- rho_N(x) = Sum_{j<x, (j,N)=1} mu(j)/phi(j), exact running sums
    cop = coprime_mask(K, N)
    term = np.zeros(K + 1)
    idx = np.arange(1, K + 1)
    term[1:] = (mu[1:K + 1].astype(np.float64)
                / phi[1:K + 1].astype(np.float64)) * cop[1:K + 1]
    rho = np.cumsum(term)

    print("rho_N(x) = Sum_{j<=x,(j,N)=1} mu(j)/phi(j)   "
          "vs Huang-Li Lemma 1 bound exp(-c sqrt(log x))")
    print(f"{'x':>10} {'rho(x)':>12} {'exp(-sqrt(log x))':>18} "
          f"{'ratio':>8}")
    for e in range(0, 8):
        x = min(10 ** e, K)
        b = math.exp(-math.sqrt(math.log(max(x, 2))))
        print(f"{x:>10} {rho[x]:>12.6f} {b:>18.6f} "
              f"{abs(rho[x])/b:>8.3f}")

    # --- B_w for single-divisor weights w_k = [d0 | k]
    print("\nB_w for w_k = [d0 | k]  (b = delta_{d0}):  exact brute "
          "force over k < K, vs the factorised formula")
    print(f"{'d0':>10} {'K/d0':>10} {'B_w (brute)':>14} "
          f"{'formula':>12} {'|B_w|*phi(d0)':>14}")
    ks = np.arange(1, K)
    base = (mu[1:K].astype(np.float64) / phi[1:K].astype(np.float64)
            ) * cop[1:K]
    for d0 in (1, 2, 3, 5, 30, 210, 2310, 30030,
               int(N ** (theta - 0.5)), K // 8, K // 2, K - 1):
        if d0 < 1 or d0 >= K:
            continue
        if math.gcd(d0, N) != 1 or mu[d0] == 0:
            # keep the table to squarefree d0 coprime to N
            d0 += 1
            while d0 < K and (math.gcd(d0, N) != 1 or mu[d0] == 0):
                d0 += 1
            if d0 >= K:
                continue
        sel = (ks % d0 == 0)
        brute = float(base[sel].sum())
        # formula: mu(d0)/phi(d0) * rho_{d0 N}(K/d0)
        x = K // d0
        cop2 = coprime_mask(max(x, 1), d0 * (N // math.gcd(N, d0)))
        t = np.zeros(max(x, 1) + 1)
        if x >= 1:
            t[1:] = (mu[1:x + 1].astype(np.float64)
                     / phi[1:x + 1].astype(np.float64)) * cop2[1:x + 1]
        form = (mu[d0] / phi[d0]) * float(t.sum())
        print(f"{d0:>10} {K//d0:>10} {brute:>14.8f} {form:>12.8f} "
              f"{abs(brute)*phi[d0]:>14.6f}")

    print("\nReading: |B_w| * phi(d0) is exactly |rho(K/d0)|.  It is of "
          "order 1 only when K/d0 is O(1) -- i.e. the weight's Mobius "
          "transform must sit AT the truncation point K = N^theta'.  "
          "For d0 <= N^(theta'-1/2), the range BV admits, it is "
          "exponentially damped.")
    print("DONE")


if __name__ == "__main__":
    main()
