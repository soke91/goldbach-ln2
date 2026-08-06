# -*- coding: utf-8 -*-
"""
Transform Lab, session 11 (increment 231): does the range-grouped
bilinear form survive the only manoeuvre available to it?

Session 10 left the question in analysis rather than measurement: the
range-grouped form

    C_R(N) = -Sum_{w} mu(w) H_{R,w},
    H_{R,w} := Sum_{p in R, p not| w} log p * Lambda(N - pw) / log(pw)

is a bilinear form in (p, w) with mu on one side and Lambda on the
other. Is it provable? The standard route for such a form is Cauchy-
Schwarz in the rough variable -- here w, carrying mu -- followed by
expanding the square, which is the dispersion method:

    |C_R|^2 <= ( Sum_{w<=W} mu^2(w) ) * ( Sum_{w<=W} H_{R,w}^2 ),

and the second factor expands into Sum_{p,p'} Sum_w Lambda(N-pw)
Lambda(N-p'w), a correlation of primes along two linear forms.

WHY THIS IS THE DECIDING TEST. The off-diagonal terms of that expansion
have POSITIVE main terms -- Sum_w Lambda(N-pw) Lambda(N-p'w) is a count
of representations and is of size ~ c W -- so Cauchy-Schwarz discards
the sign information that made C_R small in the first place and
reconstructs something of trivial size. That is the identical disease
Proposition E diagnosed for the circle method on C(N), and the identical
disease measured for the Cauchy-Schwarz step in REVIEW_VERDICT #6
(a loss of order sqrt(K)). The question is whether it recurs here, and
by how much.

BACK-OF-ENVELOPE, BEFORE MEASURING. With W ~ N/P, Sum_w Lambda(N-pw) ~
N/p and hence Sum_w H^2 ~ (log^2 P / log^2 N) * (number of pairs) * cW,
giving |C_R|^2 <~ c N^2 / log^2 N, i.e. |C_R| <~ sqrt(c) N / log N.
Since there are ~log N ranges and the wall needs Sum_R |C_R| = o(N),
the per-range demand is o(N / log N). So the Cauchy-Schwarz bound lands
EXACTLY at the demand, with no margin -- the same zero-margin verdict
as Proposition E. This script checks that against the arithmetic.

MEASUREMENTS AND CRITERIA, with the references on the same line.
 (A) per range: |C_R|, the Cauchy-Schwarz bound CS_R, and the loss
     CS_R / |C_R|.
 (B) the aggregate Sum_R CS_R against N, which is the scale the wall
     needs to beat. NULL / REFERENCE: the wall needs Sum_R |C_R| = o(N),
     so Sum_R CS_R / N >= 1 means the manoeuvre has already spent the
     entire budget before any arithmetic is done.
        MARGIN     iff Sum_R CS_R / N is below 1 and falling.
        NO MARGIN  iff it sits at or above 1, or is flat.
 (C) the diagonal alone, Sum_R sqrt( Sum_w (diagonal part) ), which is
     the part no cancellation can remove. If even the diagonal exceeds
     N the route is dead outright; if the diagonal is comfortable and
     only the full square is not, the damage is specifically the
     positive off-diagonal main terms, which is a sharper diagnosis.
"""
import numpy as np
import math

from lab_prime_factor_split import sieve


def main():
    X = 200_000
    mu, lam, spf, primes = sieve(X)
    NS = (50_000, 100_000, 200_000)

    agg = []
    for N in NS:
        v = np.arange(1, N)
        muv = mu[1:N].astype(np.float64)
        triv = float((np.abs(muv) * lam[N - v]).sum())
        ps = primes[primes < N]
        rows = []
        b = 2
        while b < N:
            hi = min(2 * b, N)
            sel = ps[(ps >= b) & (ps < hi)]
            if len(sel) == 0:
                b *= 2
                continue
            W = (N - 1) // b                      # w < N/p for p >= b
            H = np.zeros(W + 1)
            Hd = np.zeros(W + 1)                  # diagonal accumulator
            for p in sel:
                p = int(p)
                wmax = min(W, (N - 1) // p)
                if wmax < 1:
                    continue
                ww = np.arange(1, wmax + 1)
                t = (math.log(p) * lam[N - p * ww]
                     / np.log((p * ww).astype(np.float64)))
                ok = (ww % p != 0)                # enforce p not| w
                t = np.where(ok, t, 0.0)
                H[1:wmax + 1] += t
                Hd[1:wmax + 1] += t * t
            sqf = (mu[1:W + 1] != 0).astype(np.float64)
            C_R = -float(np.dot(mu[1:W + 1].astype(np.float64),
                                H[1:W + 1]))
            n_sqf = float(sqf.sum())
            L2 = float(np.dot(H[1:W + 1] * sqf, H[1:W + 1]))
            L2d = float(np.dot(Hd[1:W + 1], sqf))
            CS = math.sqrt(n_sqf * L2)
            CSd = math.sqrt(n_sqf * L2d)
            rows.append((b, hi, len(sel), C_R, CS, CSd))
            b *= 2
        agg.append((N, triv, rows))

    print("(A) per range at N = 200000 -- what Cauchy-Schwarz costs")
    N, triv, rows = agg[-1]
    print(f"{'p range':>16} {'#p':>6} {'|C_R|':>11} {'CS_R':>12} "
          f"{'CS/|C_R|':>10} {'CS_R/N':>9} {'diag CS':>12} {'d/CS':>7}")
    for (b, hi, npz, C_R, CS, CSd) in rows:
        print(f"{b:>7}-{hi:>8} {npz:>6} {abs(C_R):>11.1f} {CS:>12.1f} "
              f"{CS/max(abs(C_R),1e-9):>10.1f} {CS/N:>9.4f} "
              f"{CSd:>12.1f} {CSd/CS:>7.4f}")

    print("\n(B) the aggregate -- has the manoeuvre spent the budget?")
    print(f"{'N':>8} {'Sum|C_R|':>11} {'Sum CS_R':>12} {'triv':>11} "
          f"{'CS/N':>8} {'CS/triv':>9} {'CS/Sum|C_R|':>12}")
    for (N, triv, rows) in agg:
        c = sum(abs(r[3]) for r in rows)
        s = sum(r[4] for r in rows)
        print(f"{N:>8} {c:>11.1f} {s:>12.1f} {triv:>11.1f} "
              f"{s/N:>8.3f} {s/triv:>9.3f} {s/c:>12.2f}")
    print("    the wall needs Sum_R |C_R| = o(N); CS/N >= 1 means the")
    print("    Cauchy-Schwarz step alone has already spent the budget")

    print("\n(C) the diagonal alone -- the part no cancellation removes")
    print(f"{'N':>8} {'Sum diag CS':>13} {'/N':>8} {'Sum CS':>12} "
          f"{'diag share':>11}")
    for (N, triv, rows) in agg:
        d = sum(r[5] for r in rows)
        s = sum(r[4] for r in rows)
        print(f"{N:>8} {d:>13.1f} {d/N:>8.3f} {s:>12.1f} "
              f"{d/s:>11.4f}")
    print("    a small diagonal share means the damage is specifically")
    print("    the positive off-diagonal main terms of the expansion")
    print("DONE")


if __name__ == "__main__":
    main()
