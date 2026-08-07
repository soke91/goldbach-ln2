# -*- coding: utf-8 -*-
"""
conj:wall item 4 -- the zeta-ordinate regression, with the coin control.
(v1_verify2, Phase 1, blind.)

STATEMENT UNDER TEST, verbatim:

  "G is Gaussian in distribution but not phase-random in log N.
   Regressing G on cos(gamma log N), sin(gamma log N) for the first ten
   ordinates of the zeta zeros gives R^2 = 3.90e-3 against a
   200-surrogate maximum of 5.09e-6, and every ordinate individually at
   z >= 23. ... The 0.39% is a floor, not the share --- the zeros do not
   stop at the tenth."

WHY A COIN CONTROL IS MANDATORY HERE.  lem:coin, stated by this same
paper, says: replace mu(v) by arbitrary signs eps(v) on the squarefree
support; V(N) is unchanged, so Z_eps = C_eps/sqrt(V) is the same kind of
object, and "any estimator whose output is reproduced when mu is
replaced by eps is not measuring mu."  Item 4 reports no such control.
It is the one item of conj:wall that makes a POSITIVE claim, so it is
the one where the control matters most.

There is a concrete reason to expect leakage.  C(N) = sum_v mu(v)
Lambda(N-v) contains Lambda, and the explicit formula puts oscillations
e(gamma log N / 2pi) into any smooth average of Lambda.  A regression of
C/sqrt(V) on cos(gamma log N) can therefore pick up the PRIMES' zeros
rather than any property of mu.  The coin control separates the two: it
keeps Lambda and V exactly, and destroys only mu.

PRE-REGISTRATION (fixed before this ran).

  Decision rule.
    (a) Reproduce R^2 for the first ten ordinates on the field where
        prop:V reproduces (1e5 < N <= 1.6e7), and the 200-surrogate
        maximum, surrogate frequencies drawn uniformly from the range
        the ten ordinates span.
    (b) Run the identical regression on Z_eps built from one random sign
        vector on the squarefree support, and on a second, independent
        one.
    (c) Report each ordinate individually against the surrogate spread.

    PASSES H7 : R^2(mu) exceeds R^2(coin) by a factor >= 10, i.e. the
                statistic is substantially about mu.
    FAILS H7  : R^2(coin) is within a factor 3 of R^2(mu). Then the
                measurement is reproduced with mu destroyed, and by the
                paper's own lem:coin it is not measuring mu.
    AMBIGUOUS : between the two.

  Prediction written before running.  I expect (a) to REPRODUCE -- the
  effect is enormous (1500x the null) and enormous effects usually
  survive re-measurement. I predict H7 FAILS or is AMBIGUOUS: the
  oscillation is carried by Lambda through the explicit formula, so the
  coin field should show a comparable R^2. If instead the coin R^2 sits
  at the surrogate null, item 4 is clean and the finding is only that
  the control was not reported.

  Falsifier for the prediction: R^2(coin) at the surrogate level.
"""

import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CACHE = os.path.join(ROOT, "v1_verify2_log", "cache")
sys.path.insert(0, HERE)

# first ten ordinates of the nontrivial zeros of zeta
GAMMAS = np.array([14.134725142, 21.022039639, 25.010857580,
                   30.424876126, 32.935061588, 37.586178159,
                   40.918719012, 43.327073281, 48.005150881,
                   49.773832478])


def r2_of(y, L, freqs):
    """R^2 of the least-squares fit of y on {cos(f L), sin(f L)}.

    Built through normal equations so no n x 2k design matrix is held.
    y must already be centred.
    """
    k = len(freqs)
    XtX = np.empty((2 * k, 2 * k))
    Xty = np.empty(2 * k)
    cols = []
    for f in freqs:
        a = f * L
        cols.append(np.cos(a))
        cols.append(np.sin(a))
    for i in range(2 * k):
        Xty[i] = np.dot(cols[i], y)
        for j in range(i, 2 * k):
            v = np.dot(cols[i], cols[j])
            XtX[i, j] = XtX[j, i] = v
    del cols
    try:
        beta = np.linalg.solve(XtX, Xty)
    except np.linalg.LinAlgError:
        beta = np.linalg.lstsq(XtX, Xty, rcond=None)[0]
    return float(np.dot(beta, Xty) / np.dot(y, y))


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 16_000_000
    nsur = int(sys.argv[2]) if len(sys.argv) > 2 else 200

    z = np.load(os.path.join(CACHE, f"field_{X}.npz"))
    goodm = (z["V"] > 0) & (z["W"] > 0)
    Ni = z["N"][goodm]
    Nf = Ni.astype(np.float64)
    fld = Nf > 1e5
    Ni, Nf = Ni[fld], Nf[fld]
    C = z["C"][goodm][fld]
    V = z["V"][goodm][fld]
    cell = z["cell"][goodm][fld]
    n = len(Nf)
    L = np.log(Nf)

    print("audit_zero_spectrum   (v1_verify2 Phase 1, blind)")
    print(f"field: 1e5 < N <= {X:,}  ->  {n:,} values")
    print("=" * 74)

    def demean(field):
        out = field.copy()
        for c in range(32):
            m = cell == c
            if m.sum() > 2:
                out[m] -= out[m].mean()
        return out - out.mean()

    G = demean(C / np.sqrt(V))

    # ---------------------------------------------------- (a)
    print()
    print("--- (a) the ten ordinates, and the surrogate null ---------------")
    r2_mu = r2_of(G, L, GAMMAS)
    print(f"    R^2 (first ten ordinates, mu)     = {r2_mu:.4e}"
          f"   [paper: 3.90e-3]")
    print(f"    naive null E[R^2] = 2k/n          = {2 * len(GAMMAS) / n:.4e}")

    rng = np.random.default_rng(20260807)
    lo, hi = GAMMAS.min(), GAMMAS.max()
    sur = np.empty(nsur)
    for s in range(nsur):
        f = rng.uniform(lo, hi, size=len(GAMMAS))
        sur[s] = r2_of(G, L, f)
    print(f"    {nsur} surrogates, frequencies ~ U({lo:.1f}, {hi:.1f}):")
    print(f"      max    = {sur.max():.4e}   [paper: 5.09e-6]")
    print(f"      mean   = {sur.mean():.4e}    sd = {sur.std(ddof=1):.4e}")
    print(f"      R^2(mu) / surrogate max = {r2_mu / sur.max():.1f}")
    print()

    # ---------------------------------------------------- (b) H7
    print("--- (b) the coin control that lem:coin requires ------------------")
    print("    replacing mu by random signs on the squarefree support;")
    print("    V(N), Lambda and the cells are all untouched.")
    from lab_field_build import (smallest_prime_factor, von_mangoldt,
                                 mobius, fftconv_prefix)
    spf = smallest_prime_factor(X)
    lam, primes = von_mangoldt(X, spf)
    del spf
    mu = mobius(X, primes)
    supp = (mu != 0)
    print(f"    squarefree support: {supp.sum():,}")

    coin_r2 = []
    for trial in range(2):
        rg = np.random.default_rng(1234 + trial)
        eps = np.zeros(X + 1)
        eps[supp] = rg.integers(0, 2, size=int(supp.sum())) * 2.0 - 1.0
        eps[0] = 0.0
        Ce = fftconv_prefix(lam, eps, X)
        Ge = demean(Ce[Ni] / np.sqrt(V))
        r2e = r2_of(Ge, L, GAMMAS)
        coin_r2.append(r2e)
        print(f"    coin trial {trial + 1}: R^2 = {r2e:.4e}"
              f"   ratio R^2(mu)/R^2(coin) = {r2_mu / r2e:7.2f}")
        del Ce, Ge, eps
    print()
    ratio = r2_mu / max(coin_r2)
    if ratio >= 10:
        verdict = "PASSES H7 (the statistic is about mu)"
    elif ratio <= 3:
        verdict = "FAILS H7 (reproduced with mu destroyed)"
    else:
        verdict = "AMBIGUOUS"
    print(f"    VERDICT: {verdict}   (best coin ratio {ratio:.2f})")
    print()

    # ---------------------------------------------------- (c)
    print("--- (c) each ordinate individually -------------------------------")
    sur1 = np.empty(nsur)
    for s in range(nsur):
        sur1[s] = r2_of(G, L, np.array([rng.uniform(lo, hi)]))
    m1, s1 = sur1.mean(), sur1.std(ddof=1)
    print(f"    single-frequency surrogate null: mean {m1:.3e}, "
          f"sd {s1:.3e}, max {sur1.max():.3e}")
    print(f"    {'gamma':>12}{'R^2':>12}{'z':>9}{'coin R^2':>12}{'coin z':>9}")
    rg = np.random.default_rng(1234)
    eps = np.zeros(X + 1)
    eps[supp] = rg.integers(0, 2, size=int(supp.sum())) * 2.0 - 1.0
    Ce = fftconv_prefix(lam, eps, X)
    Ge = demean(Ce[Ni] / np.sqrt(V))
    for g in GAMMAS:
        a = r2_of(G, L, np.array([g]))
        b = r2_of(Ge, L, np.array([g]))
        print(f"    {g:>12.6f}{a:>12.3e}{(a - m1) / s1:>9.1f}"
              f"{b:>12.3e}{(b - m1) / s1:>9.1f}")
    print("    [paper: 'every ordinate individually at z >= 23']")
    print()
    print("    note: the paper does not state the surrogate frequency")
    print("    range. A null drawn from far higher frequencies would be")
    print("    tighter and would inflate every z above.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
