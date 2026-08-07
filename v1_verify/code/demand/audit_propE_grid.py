# -*- coding: utf-8 -*-
"""
Re-verification of Proposition 7 (`prop:E`, "the circle method has zero
margin") of v1/paper/wall_v1.tex.

THE STATEMENT UNDER TEST, verbatim (wall_v1.tex and theorem_A.tex
Section `sec:circle`):

  (i)  ||S_Lambda||_2 ||S_mu||_2 ~ (6/pi^2)^{1/2} N (log N)^{1/2},
       exceeding the trivial bound by a growing factor;
  (ii) any bound sup_alpha |S_mu| * ||S_Lambda||_1 is at least
       ||S_mu||_2 ||S_Lambda||_1 >> N^{1/2} * N^{1/2} = N.
  "The measured margin N/(sup|S_mu| ||S_Lambda||_1) is
   0.168, 0.175, 0.158, 0.152 at N = 2^14,...,2^20 -- below 1 and
   decaying."

with the quoted table

     N              2^14    2^16    2^18    2^20
     ||S_mu||_2/rtN 0.7798  0.7797  0.7797  0.7797
     sup|S_mu|/rtN  3.058   2.742   2.853   2.801
     ||S_L||_1/rtN  1.946   2.084   2.219   2.346
     margin         0.168   0.175   0.158   0.152

computed, per theorem_A.tex, by "exact FFT on a 4N-point grid".

THE OBJECTION BEING TESTED -- the grid. Two of the four quantities are
NOT grid-invariant:

  * ||S_mu||_2 is exact on any grid of >= N points (Parseval), so it is
    a control, not a measurement -- which is why it reproduces
    sqrt(6/pi^2) to four places.
  * sup_alpha |S_mu| over a finite grid is a LOWER bound for the true
    supremum; a coarser grid can only understate it.
  * ||S_Lambda||_1 = int_0^1 |S_Lambda| is a Riemann sum. S_Lambda has
    major-arc peaks of width ~1/N and height ~N/phi(q); a 4N-point grid
    puts about 4 samples across each peak. An L^1 Riemann sum over such
    a spiky integrand converges slowly.

Since the margin is sup|S_mu| * ||S_Lambda||_1, both grid errors push
the margin the SAME way -- upward -- so the reported margin is an
overstatement of the route's viability, and the "decaying" reading is
a claim about the third significant figure of a Riemann sum.

METHOD HERE. The same four quantities on grids of g*N points for
g = 4, 8, 16, 32, 64, and the margin recomputed at each. Written from
the statement; the FFT is the definition of the trigonometric sum on
the grid, so there is nothing to re-implement except the grid.

PRE-REGISTRATION (written before the run).

  (1) RULE. The quoted table is a measurement only if it is converged:
      require every entry to move by less than 1% between g = 4 and
      g = 64. Anything larger and the third digit is a grid artefact.
  (2) PREDICTION, recorded so it cannot be reported as a surprise.
      ||S_mu||_2 will be exactly grid-invariant. sup|S_mu| will RISE
      with g. ||S_Lambda||_1 will RISE with g, by several percent,
      because the L^1 norm of a spiky function is underestimated by
      coarse sampling. So the margin will FALL with g, and I predict
      the g=4 values 0.168...0.152 are too HIGH by at least a few
      percent.
  (3) The verdict of Proposition 7 -- margin below 1 -- is NOT in
      doubt: it sits at 0.16, a factor 6 from the boundary, and both
      grid errors move it further from 1. What is tested is whether
      the quoted digits, and the word "decaying", are measurements.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def sieve_mu_lambda(X):
    mu = np.ones(X + 1, dtype=np.int64)
    is_p = np.zeros(X + 1, dtype=bool)
    rem = np.arange(X + 1, dtype=np.int64)
    mu[0] = 0
    for p in range(2, X + 1):
        if rem[p] == p:
            is_p[p] = True
            mu[p::p] *= -1
            rem[p::p] //= p
            if p * p <= X:
                mu[p * p::p * p] = 0
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in range(2, X + 1):
        if is_p[p]:
            lg = math.log(p)
            q = p
            while q <= X:
                lam[q] = lg
                q *= p
    return mu.astype(np.float64), lam


def norms(a, N, M):
    """|a-hat| on an M-point grid: L2 (normalised), L1, sup."""
    F = np.fft.rfft(np.pad(a[: N + 1], (0, M - N - 1)), M)
    A = np.abs(F)
    # rfft gives half the circle; the other half is the mirror image,
    # so both L1 and L2 means are the mean over the full circle.
    if M % 2 == 0:
        w = np.concatenate(([1.0], 2.0 * np.ones(len(A) - 2), [1.0]))
    else:
        w = np.concatenate(([1.0], 2.0 * np.ones(len(A) - 1)))
    l1 = float((w * A).sum() / M)
    l2 = math.sqrt(float((w * A ** 2).sum() / M))
    return l1, l2, float(A.max())


def main():
    print("Re-verification of Proposition 7: grid convergence")
    print()
    EXPS = (14, 16, 18, 20)
    GS = (4, 8, 16, 32, 64)

    XMAX = 1 << max(EXPS)
    mu, lam = sieve_mu_lambda(XMAX)

    for e in EXPS:
        N = 1 << e
        print(f"N = 2^{e} = {N}")
        hdr = (f"    {'grid':>7} {'||S_mu||_2/rtN':>15} "
               f"{'sup|S_mu|/rtN':>14} {'||S_L||_1/rtN':>14} "
               f"{'(i)/N':>8} {'margin':>9}")
        print(hdr)
        print("    " + "-" * (len(hdr) - 4))
        vals = []
        for g in GS:
            M = g * N
            _, l2mu, supmu = norms(mu, N, M)
            l1lam, l2lam, _ = norms(lam, N, M)
            rt = math.sqrt(N)
            margin = N / (supmu * l1lam)
            vals.append((l2mu / rt, supmu / rt, l1lam / rt,
                         l2mu * l2lam / N, margin))
            print(f"    {g:>5}N {l2mu/rt:>15.4f} {supmu/rt:>14.4f} "
                  f"{l1lam/rt:>14.4f} {l2mu*l2lam/N:>8.4f} "
                  f"{margin:>9.4f}")
        a, b = vals[0], vals[-1]
        names = ("||S_mu||_2", "sup|S_mu|", "||S_L||_1", "(i)", "margin")
        moves = [abs(y / x - 1) for x, y in zip(a, b)]
        print(f"    move from 4N to 64N: " +
              ", ".join(f"{nm} {mv:+.2%}" for nm, mv in zip(names, moves)))
        print(f"    v1 quotes at 4N: margin "
              f"{ {14:0.168, 16:0.175, 18:0.158, 20:0.152}[e] }")
        print()

    print("(1) RULE was: converged if every entry moves < 1% from 4N to")
    print("    64N. See the 'move' lines.")
    print("(3) the verdict 'margin < 1' is unaffected either way.")
    print("DONE")


if __name__ == "__main__":
    main()
