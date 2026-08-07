# -*- coding: utf-8 -*-
"""
Re-verification of Lemma 13 (`lem:MP`), the second-moment identity, as
STATED in v1/paper/wall_v1.tex.

THE STATEMENT UNDER TEST (wall_v1.tex, Lemma `lem:MP`), verbatim:

    For any X,
        sum_{N<=X} C(N)^2  =  sum_{|h|<X} M(h) P(h),
        M(h) = sum_v mu(v) mu(v+h),   P(h) = sum_w Lambda(w) Lambda(w+h),
    the inner sums being over the ranges that keep both arguments
    below X.

and the paper's proof: "expanding C(N)^2 as a double sum over (n,v) and
(n',v') with n+v = n'+v' = N and summing over N leaves the single
constraint n-n' = v'-v =: h. Collecting the Lambda-pairs and the
mu-pairs separately gives the two factors."

METHOD HERE. Written from the statement only. Both sides are computed
by direct enumeration (O(X^2) double loops in numpy, no FFT, no
Parseval), at sizes small enough that the enumeration is the
definition. C(N) = sum_{n<N} Lambda(n) mu(N-n) is built by an explicit
convolution over n.

PRE-REGISTRATION (written before the run).

  (1) RULE: |LHS/RHS - 1| < 1e-10 at every X tested. The lemma says the
      two sides are equal exactly; anything else refutes it as stated.

  (2) PREDICTION, recorded so it cannot be reported as a surprise.
      The paper's proof drops a range coupling. Summing over N <= X
      leaves the three free variables (n, n', v) subject to n + v <= X,
      a SIMPLEX; the product M(h)P(h) sums them over a BOX. The box
      contains, in addition, every pair with X < n+v <= 2X. Those extra
      terms are the squares of truncated convolutions, all nonnegative
      in aggregate, of total size comparable to the LHS itself. So I
      predict RHS > LHS by a factor near 2 -- specifically I predict
      the ratio LHS/RHS near 0.5-0.6, drifting slowly, and NOT ->1.

  (3) The repair I expect to hold instead, tested as (C) below:
        sum_{N} C_X(N)^2  =  sum_{|h|<X} M(h) P(h)
      where C_X(N) = sum_{n+v=N, n<=X, v<=X} Lambda(n) mu(v) is the
      TRUNCATED convolution and N runs over its full support up to 2X.
      This is Parseval for the two truncated sequences and should be
      exact. If (C) passes while (1) fails, the defect is located
      precisely: the left-hand side of the published lemma is not the
      quantity the right-hand side computes.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def mobius_and_lambda(X):
    """mu and Lambda on [0, X], by plain sieving. Independent of v1."""
    mu = np.ones(X + 1, dtype=np.int64)
    primes_flag = np.zeros(X + 1, dtype=bool)
    remaining = np.arange(X + 1, dtype=np.int64)
    mu[0] = 0
    for p in range(2, X + 1):
        if remaining[p] == p:          # p untouched so far => prime
            primes_flag[p] = True
            mu[p::p] *= -1
            remaining[p::p] //= p
            p2 = p * p
            if p2 <= X:
                mu[p2::p2] = 0
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in range(2, X + 1):
        if primes_flag[p]:
            lg = math.log(p)
            q = p
            while q <= X:
                lam[q] = lg
                q *= p
    return mu.astype(np.float64), lam


def C_full(mu, lam, X):
    """C(N) = sum_{n<N} Lambda(n) mu(N-n) for N = 0..X, by direct
    enumeration over n. No FFT."""
    C = np.zeros(X + 1, dtype=np.float64)
    for n in range(2, X + 1):
        if lam[n] == 0.0:
            continue
        # N from n+1 to X ; v = N-n from 1 to X-n
        C[n + 1:X + 1] += lam[n] * mu[1:X - n + 1]
    return C


def C_truncated_full_support(mu, lam, X):
    """C_X(N) = sum_{n+v=N, 1<=n<=X, 1<=v<=X} Lambda(n) mu(v), for
    N = 0..2X. Direct enumeration."""
    C = np.zeros(2 * X + 1, dtype=np.float64)
    for n in range(2, X + 1):
        if lam[n] == 0.0:
            continue
        C[n + 1:n + X + 1] += lam[n] * mu[1:X + 1]
    return C


def autocorr_box(a, X):
    """M(h) = sum_v a(v) a(v+h) with BOTH arguments in [1, X], for
    h = 0..X-1. Direct enumeration."""
    out = np.zeros(X, dtype=np.float64)
    for h in range(X):
        out[h] = float(np.dot(a[1:X + 1 - h], a[1 + h:X + 1]))
    return out


def main():
    print("Re-verification of Lemma 13 (lem:MP) of v1/paper/wall_v1.tex")
    print("Both sides computed by direct enumeration, from the statement.")
    print()
    print("  LHS_paper = sum_{N<=X} C(N)^2        (the lemma's left side)")
    print("  RHS_paper = sum_{|h|<X} M(h) P(h)    (the lemma's right side)")
    print("  LHS_trunc = sum_{all N} C_X(N)^2     (truncated convolution,")
    print("                                        full support up to 2X)")
    print()
    hdr = (f"{'X':>7} {'LHS_paper':>16} {'RHS_paper':>16} "
           f"{'LHS/RHS':>9} {'LHS_trunc':>16} {'trunc/RHS - 1':>15}")
    print(hdr)
    print("-" * len(hdr))

    ok_paper = True
    ok_repair = True
    ratios = []
    for X in (500, 1000, 2000, 4000, 8000, 16000):
        mu, lam = mobius_and_lambda(X)

        C = C_full(mu, lam, X)
        lhs_paper = float((C[: X + 1] ** 2).sum())

        M = autocorr_box(mu, X)
        P = autocorr_box(lam, X)
        # |h| < X, so h = 0 once and every h != 0 twice (h and -h;
        # M(-h)=M(h), P(-h)=P(h) by relabelling).
        rhs_paper = float(M[0] * P[0] + 2.0 * np.dot(M[1:], P[1:]))

        Ct = C_truncated_full_support(mu, lam, X)
        lhs_trunc = float((Ct ** 2).sum())

        r = lhs_paper / rhs_paper
        ratios.append((X, r))
        ok_paper &= abs(r - 1.0) < 1e-10
        ok_repair &= abs(lhs_trunc / rhs_paper - 1.0) < 1e-10

        print(f"{X:>7} {lhs_paper:>16.8e} {rhs_paper:>16.8e} "
              f"{r:>9.5f} {lhs_trunc:>16.8e} "
              f"{lhs_trunc / rhs_paper - 1.0:>15.2e}")

    print()
    print(f"(1) the lemma as stated holds to 1e-10 : "
          f"{'PASS' if ok_paper else 'FAIL'}")
    print(f"(3) the repaired identity holds to 1e-10 : "
          f"{'PASS' if ok_repair else 'FAIL'}")
    print()
    print("(2) the pre-registered prediction was LHS/RHS near 0.5-0.6,")
    print("    not tending to 1. Measured:")
    for X, r in ratios:
        print(f"      X = {X:>6}   LHS/RHS = {r:.5f}")

    # ---- (D) locate the defect exactly: the missing block ----
    X = 4000
    mu, lam = mobius_and_lambda(X)
    Ct = C_truncated_full_support(mu, lam, X)
    C = C_full(mu, lam, X)
    below = float((Ct[: X + 1] ** 2).sum())
    above = float((Ct[X + 1:] ** 2).sum())
    agree = float(np.abs(Ct[: X + 1] - C[: X + 1]).max())
    print()
    print(f"(D) where the difference sits, at X = {X}:")
    print(f"      max_N<=X |C_X(N) - C(N)|            = {agree:.3e}"
          f"   (0 => the two agree on N <= X)")
    print(f"      sum_{{N<=X}}   C_X(N)^2               = {below:.8e}")
    print(f"      sum_{{X<N<=2X}} C_X(N)^2              = {above:.8e}")
    print(f"      share of the RHS carried by N > X   = "
          f"{above / (below + above):.4f}")
    print()
    print("      The right-hand side of the published lemma counts the")
    print("      block X < N <= 2X, which the left-hand side excludes.")
    print("      That block is not an error term: it carries the share")
    print("      printed above.")
    print("DONE")


if __name__ == "__main__":
    main()
