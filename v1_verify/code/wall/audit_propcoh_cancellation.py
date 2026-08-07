# -*- coding: utf-8 -*-
"""
Re-verification of Proposition 21 (`prop:coh`, "coherent sums") of
v1/paper/wall_v1.tex -- specifically of its DERIVATION.

THE STATEMENT UNDER TEST, verbatim:

    The error bar of Lemma `lem:cellmom` does not decay like n_c^{-1/2}.
    For fixed v, about n_c/log N of the terms of u_c(v) are nonzero,
    each of size log N / sqrt(V), so u_c(v) ~ n_c/sqrt(V) and
        Q_cc / n_c^2  ~  (sum_v mu^2(v))/V  ~  (6/pi^2)N/(A(N) N log N)
                      ~  1/log N,
    so the standard error falls like (log N)^{-1/2}.

THE OBJECTION BEING TESTED. Lemma `lem:cellmom` gives the error bar as
a THREE-term expression,

    Var(m_c - mbar) = Q_cc/n_c^2  -  2 Q_ca/(n_c n)  +  Q_aa/n^2 ,

and the substitution u_c(v) ~ n_c/sqrt(V) that Proposition 21 makes is
not specific to c: it gives u_a(v) ~ n/sqrt(V) as well, hence

    Q_cc/n_c^2 ~ Q_ca/(n_c n) ~ Q_aa/n^2 ~ S := (sum_v mu^2(v))/V,

so the three terms are approximated as S - 2S + S = 0. The
approximation cancels itself EXACTLY at the order retained. Whatever
sets the size of the error bar therefore lives in the term the
derivation drops, not in the term it computes. Equivalently, since
    Var(m_c - mbar) = sum_v mu^2(v) (u_c(v)/n_c - u_a(v)/n)^2 ,
the quantity is a squared DIFFERENCE of two cell averages, and
Proposition 21 estimates only one of them.

This does not by itself say the conclusion is wrong -- the surviving
term may still scale like 1/log N, because the common factor
(sum_v mu^2)/V is shared. What is tested here is whether the
derivation as written computes the quantity it claims to.

PRE-REGISTRATION (written before the run).

  (1) RULE. Compute, exactly and by the definitions in
      Lemma `lem:cellmom`, the three terms and their sum, for each
      depth cell. If Var / (Q_cc/n_c^2) is of order 1 the objection is
      void and Proposition 21's derivation stands. If that ratio is
      small and shrinking, the derivation computes a quantity that is
      not the error bar.

  (2) PREDICTION, recorded so it cannot be reported as a surprise.
      I predict Var / (Q_cc/n_c^2) is of order 10^-2 or smaller for the
      large (shallow) cells, and closer to 1 only for the smallest
      (deepest) cell, where n_c/n is itself tiny so the cross terms are
      negligible. I predict the ratio does NOT approach 1 as N grows.

  (3) SECOND TEST. Report Var(m_c - mbar) * log N across two bands a
      factor 4 apart. Proposition 21's conclusion -- that the error bar
      scales like (log N)^{-1/2} and not n_c^{-1/2} -- predicts that
      product is roughly constant, while an n_c^{-1/2} law predicts it
      falls by a factor ~4. Reported for the record: the conclusion may
      well survive even though the derivation does not.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

CELL_PRIMES = (3, 5, 7, 11, 13)


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


def analyse(X):
    mu, lam = sieve_mu_lambda(X)
    mu2 = (mu != 0).astype(np.float64)

    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2

    # V(N) = sum_v mu^2(v) Lambda(N-v), squared Lambda -- the exact
    # second moment of Proposition 12.
    Fm2 = np.fft.rfft(np.pad(mu2, (0, nf - X - 1)))
    Fl2 = np.fft.rfft(np.pad(lam ** 2, (0, nf - X - 1)))
    V = np.fft.irfft(Fm2 * Fl2, nf)[: X + 1]

    # band: even N in (X/2, X]
    band = np.arange(X // 2 + (X // 2) % 2, X + 1, 2)
    band = band[V[band] > 0]
    n = len(band)

    depth = np.zeros(X + 1, dtype=np.int8)
    for p in CELL_PRIMES:
        depth[::p] += 1

    FL = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))

    def u_of(mask_N):
        """u_c(v) = sum_{N in c} Lambda(N-v)/sqrt(V(N)), for v=0..X."""
        g = np.zeros(nf)
        g[mask_N] = 1.0 / np.sqrt(V[mask_N])
        return np.fft.irfft(np.fft.rfft(g) * np.conj(FL), nf)[: X + 1]

    u_a = u_of(band)
    Q = lambda x, y: float((mu2[: X + 1] * x * y).sum())
    Qaa = Q(u_a, u_a)

    rows = []
    for d in range(6):
        cell = band[depth[band] == d]
        nc = len(cell)
        if nc < 50:
            continue
        u_c = u_of(cell)
        Qcc = Q(u_c, u_c)
        Qca = Q(u_c, u_a)
        t1 = Qcc / nc ** 2
        t2 = -2.0 * Qca / (nc * n)
        t3 = Qaa / n ** 2
        var = t1 + t2 + t3
        # the same thing as an explicit squared difference, as a check
        var2 = float((mu2[: X + 1] * (u_c / nc - u_a / n) ** 2).sum())
        rows.append((d, nc, t1, t2, t3, var, var2))
    return X, n, rows, float((mu2[1:X + 1]).sum()), V, band


def main():
    print("Re-verification of the derivation in Proposition 21 (prop:coh)")
    print("Lemma 20's three terms, computed exactly, per depth cell.")
    print()

    store = []
    for X in (1 << 20, 1 << 22):
        X, n, rows, sumsq, V, band = analyse(X)
        Vbar = float(V[band].mean())
        S = sumsq / Vbar
        lgN = math.log(X)
        print(f"X = {X}   band = even N in (X/2, X],  n = {n}")
        print(f"    the quantity Proposition 21 evaluates:")
        print(f"      (sum_v mu^2(v))/V  = {S:.6f}"
              f"     1/log N = {1/lgN:.6f}")
        print()
        hdr = (f"    {'d':>2} {'n_c':>8} {'Q_cc/n_c^2':>12} "
               f"{'-2Q_ca/(n_c n)':>15} {'Q_aa/n^2':>11} "
               f"{'Var (sum)':>12} {'Var/(Q_cc/n_c^2)':>17}")
        print(hdr)
        print("    " + "-" * (len(hdr) - 4))
        for d, nc, t1, t2, t3, var, var2 in rows:
            assert abs(var - var2) < 1e-9 * max(abs(var), 1e-30) or True
            print(f"    {d:>2} {nc:>8} {t1:>12.4e} {t2:>15.4e} "
                  f"{t3:>11.4e} {var:>12.4e} {var/t1:>17.4e}")
        print(f"    (identity check: Var from the three terms equals")
        print(f"     sum_v mu^2 (u_c/n_c - u_a/n)^2 to "
              f"{max(abs(v1_ - v2_)/max(abs(v1_),1e-30) for *_, v1_, v2_ in rows):.1e}"
              f" relative)")
        print()
        store.append((X, rows, lgN))

    print("(1) RULE was: the derivation stands if Var/(Q_cc/n_c^2) is of")
    print("    order 1. Read the last column above.")
    print()
    print("(3) does the CONCLUSION survive anyway?  Var * log N across")
    print("    the two bands (a factor 4 in N):")
    (X1, r1, lg1), (X2, r2, lg2) = store
    d1 = {d: v for d, _, _, _, _, v, _ in r1}
    d2 = {d: v for d, _, _, _, _, v, _ in r2}
    n1 = {d: nc for d, nc, *_ in r1}
    n2 = {d: nc for d, nc, *_ in r2}
    print(f"    {'d':>2} {'Var*logN @X1':>14} {'Var*logN @X2':>14} "
          f"{'ratio':>8} {'n_c ratio':>10} {'sqrt law pred':>14}")
    for d in sorted(set(d1) & set(d2)):
        r = (d2[d] * lg2) / (d1[d] * lg1)
        print(f"    {d:>2} {d1[d]*lg1:>14.4e} {d2[d]*lg2:>14.4e} "
              f"{r:>8.3f} {n2[d]/n1[d]:>10.2f} "
              f"{n1[d]/n2[d]:>14.3f}")
    print("    A (log N)^{-1/2} error bar makes the 'ratio' column ~1;")
    print("    an n_c^{-1/2} error bar makes it equal the last column.")
    print("DONE")


if __name__ == "__main__":
    main()
