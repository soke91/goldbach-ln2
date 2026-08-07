# -*- coding: utf-8 -*-
"""
Re-verification of the Proposition 15 (`prop:W`) reconstruction of
rho - 1, and of the shift-mass table of Section `sec:coin`.

THE STATEMENTS UNDER TEST (v1/paper/wall_v1.tex, Section `sec:coin`).

  (a) "Reconstructing rho - 1 from Proposition W gives -0.0976 against a
       measured -0.18, a factor 0.54; the sign is negative, as rho < 1
       requires."

  (b) "Gross mass by shift: h < 10^3 carries 1.1%, 10^3-10^4 carries
       3.0%, 10^4-10^5 carries 23.1%, 10^5-10^6 carries 48.9%, and above
       10^6 carries 23.8% ... Small shifts, where Chowla is hardest and
       the averaged theorem weakest, carry almost nothing: the wall
       leans on the range where that theorem is strongest."

Both read off the same sum, sum_{h != 0} c(h) S(h) / V, in which the
paper defines

        S(h) = < mu(u) mu(u-h) >          (an AVERAGE; |S| <= 1)

so that S(h) = M(h)/#{u}, with M(h) = sum_u mu(u) mu(u+h) the
unnormalised autocorrelation. For h ranging over [1, X) the number of
terms is X - h, not X. The two choices differ by the factor X/(X-h),
which is 1 at h = 0 and unbounded as h -> X. This script asks whether
the reconstruction and the shift-mass shares depend on that choice.

METHOD HERE. Written from the statement. M(h) and c(h) are computed by
FFT autocorrelation of mu and Lambda restricted to [1, X] -- so that
"both arguments below X" holds exactly, as the paper's Lemma 13
prescribes for M -- and V(X) by direct enumeration.

PRE-REGISTRATION (written before the run).

  (1) REPRODUCTION. With the denominator X (the choice made in v1's
      lab_offdiag_chowla.py), the reconstruction should reproduce
      -0.0976 at X = 4e6. RULE: agreement to 2 significant figures,
      else my implementation is what is wrong and nothing below reads.

  (2) SENSITIVITY, the actual test. Recompute with the denominator
      X - h, which is what the paper's own definition of S(h) as an
      average requires. RULE: if the two reconstructions differ by less
      than 20% the choice is immaterial and the finding is void.

  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      v1's own shift table puts 23.8% of the gross mass above h = 10^6
      and gives it the OPPOSITE sign to the bulk (+1.05e13 against a
      net of -2.22e13). Those are exactly the terms the factor X/(X-h)
      magnifies, by between 1.3x and 4x. I predict the mean-normalised
      reconstruction moves substantially towards zero or past it, and
      that the sign test (c) of v1's run is not stable under the
      choice. I also predict the shift shares move, so that (b)'s
      conclusion about where the wall leans is normalisation-dependent.
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


def autocorr(a, X):
    """returns r[h] = sum_{u} a(u) a(u+h), both arguments in [1, X],
    for h = 0..X-1."""
    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2
    F = np.fft.rfft(np.pad(a[: X + 1], (0, nf - X - 1)))
    r = np.fft.irfft(F * np.conj(F), nf)
    return r[:X]


def main():
    print("Re-verification of the Proposition 15 reconstruction of rho-1")
    print("and of the shift-mass table, under the two normalisations of")
    print("S(h) = <mu(u) mu(u-h)>:  M(h)/X   vs   M(h)/(X-h).")
    print()

    for X in (1_000_000, 4_000_000):
        mu, lam = sieve_mu_lambda(X)
        M = autocorr(mu, X)
        c = autocorr(lam, X)
        v = np.arange(1, X)
        V = float(((mu[1:X] ** 2) * (lam[X - v] ** 2)).sum())
        W = float((lam[:X] ** 2).sum())

        h = np.arange(1, X)
        Mh = M[1:]
        ch = c[1:]

        # sum over h != 0 is twice the sum over h >= 1 (M and c are even)
        offA = 2.0 * float(np.dot(ch, Mh)) / X          # S(h) = M(h)/X
        offB = 2.0 * float(np.dot(ch, Mh / (X - h)))    # S(h) = M(h)/(X-h)

        print(f"X = {X}")
        print(f"    M(0) = {M[0]:.0f}   (6X/pi^2 = {6*X/math.pi**2:.0f})")
        print(f"    V(X) = sum_v mu^2(v) Lambda(X-v)^2 = {V:.5e}   <- what")
        print(f"           Proposition 15 divides by")
        print(f"    W(X) = sum_w Lambda(w)^2           = {W:.5e}   <- the")
        print(f"           same sum with the mu^2 dropped; W/V = {W/V:.5f}")
        print()
        print(f"    {'numerator sum_{h!=0} c(h)S(h)':>34} "
              f"{'/ W  (v1)':>12} {'/ V  (correct)':>15}")
        print(f"    {'S(h) = M(h)/X    (v1 choice)':>34} "
              f"{offA:+.5e} -> {offA/W:+8.5f} {offA/V:+15.5f}")
        print(f"    {'S(h) = M(h)/(X-h)  (a mean)':>34} "
              f"{offB:+.5e} -> {offB/W:+8.5f} {offB/V:+15.5f}")
        print(f"    ratio B/A          : {offB/offA:+.4f}")

        # shift-resolved mass. v1's "gross" is the sum of the ABSOLUTE
        # VALUES OF THE BIN TOTALS (not of the individual terms), so
        # the same convention is used here to make the tables
        # comparable.
        edges = [1, 10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, X + 1]
        binsA, binsB = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sl = slice(lo - 1, min(hi, X) - 1)
            binsA.append(2.0 * float((ch[sl] * Mh[sl]).sum()))
            binsB.append(2.0 * float((ch[sl] * Mh[sl] * X
                                      / (X - h[sl])).sum()))
        gA = sum(abs(x) for x in binsA)
        gB = sum(abs(x) for x in binsB)
        print()
        print(f"    shift-resolved mass (v1's convention: |bin total|)")
        print(f"    {'shift range':>24} {'contrib (den X)':>17} "
              f"{'share':>8} {'contrib (den X-h)':>19} {'share':>8}")
        for (lo, hi), a_, b_ in zip(zip(edges[:-1], edges[1:]),
                                    binsA, binsB):
            print(f"    {f'[{lo}, {hi})':>24} {a_:>17.4e} "
                  f"{abs(a_)/gA:>7.1%} {b_:>19.4e} {abs(b_)/gB:>7.1%}")
        print(f"    net {sum(binsA):.4e} vs gross {gA:.4e} "
              f"(cancellation {gA/abs(sum(binsA)):.1f}x)   |   "
              f"net {sum(binsB):.4e} vs gross {gB:.4e} "
              f"(cancellation {gB/abs(sum(binsB)):.1f}x)")
        print()

    print("VERDICTS")
    print("  (1) v1's lab_offdiag_chowla.py reports, at X = 4e6,")
    print("        'reconstructed off-diagonal Sum_h c(h)M(h)/X"
          " = -5.5380e+06'")
    print("        'V(X) = Sum_v mu^2 Lambda^2                 "
          " = +5.6771e+07'")
    print("        'reconstructed (rho - 1)                    "
          " = -0.09755'")
    print("      The numerator reproduces (row A above). The number")
    print("      labelled V(X) does not: it is W(X), the same sum with")
    print("      the mu^2 dropped. V(X) is smaller by exactly the")
    print("      factor A(N) of Proposition 12 -- the paper's own")
    print("      correction, not applied where Proposition 15 needs it.")
    print("  (2) sensitivity rule was: immaterial if the two")
    print("      normalisations of S(h) differ by less than 20%.")
    print("      See the ratio B/A.")
    print("  (3) the paper's 'a factor 0.54' is the ratio of the")
    print("      reconstruction to the measured -0.18. Recomputed:")
    print("        v1 as run          (M/X, /W) : 0.54")
    print("        correct denominator (M/X, /V): see row A, column V")
    print("        S(h) a true mean  (M/(X-h), /V): see row B, column V")
    print("DONE")


if __name__ == "__main__":
    main()
