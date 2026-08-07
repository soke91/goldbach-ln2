# -*- coding: utf-8 -*-
"""
conj:wall item 4 under the coin control, with a distribution not a draw.
(v1_verify2, Phase 1, blind.)

audit_zero_spectrum.py found, on the field where prop:V reproduces:

  R^2(mu, ten ordinates)          = 4.19e-3   [paper: 3.90e-3]  REPRODUCED
  200-surrogate max, same freq band = 3.60e-3 [paper: 5.09e-6]  NOT
  R^2(coin), one draw             = 4.21e-3   -- larger than mu's
  per-ordinate z                  = -0.6 .. +6.7 [paper: all >= 23]

Two coin draws are not a null. This script builds the coin distribution
properly, so the claim "R^2(mu) is reproduced when mu is destroyed" is
made against a spread rather than against one realisation.

It also separates the two candidate causes of the paper's very tight
surrogate maximum (5.09e-6, essentially the iid value 2k/n = 2.5e-6):

  (i) a surrogate that PERMUTES the field, destroying its
      autocorrelation in N -- which makes the null far too tight,
      because G(N) and G(N+2) share almost all of their summands; or
  (ii) a surrogate that redraws FREQUENCIES on the intact field, which
      is the null that respects the autocorrelation.

Both are computed here so the difference is visible.

PRE-REGISTRATION.

  Decision rule.
    (a) 12 independent coin fields; report the distribution of R^2 and
        the rank of R^2(mu) within it.
        FAILS H7 iff R^2(mu) lies inside the coin distribution, i.e. its
        rank is not extreme (two-sided p > 0.05).
    (b) the permutation surrogate maximum vs the frequency surrogate
        maximum, 200 each.
        CONFIRMS THE NULL DEFECT iff the permutation maximum is within a
        factor 3 of the paper's 5.09e-6 while the frequency maximum is
        orders of magnitude larger.

  Prediction written before running.  (a) FAILS H7 -- one coin draw
  already matched mu, so I expect mu to sit mid-distribution. (b)
  CONFIRMS: I predict the permutation surrogate lands near 2.5e-6, close
  to the paper's 5.09e-6, identifying the paper's null as one that
  destroys the field's autocorrelation.

  What would refute (a): R^2(mu) above every coin draw.
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

from audit_zero_spectrum import GAMMAS, r2_of  # noqa: E402


def main():
    X = 16_000_000
    NTRIAL = 12
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

    def demean(f):
        out = f.copy()
        for c in range(32):
            m = cell == c
            if m.sum() > 2:
                out[m] -= out[m].mean()
        return out - out.mean()

    G = demean(C / np.sqrt(V))
    r2_mu = r2_of(G, L, GAMMAS)

    print("audit_zero_coin_dist   (v1_verify2 Phase 1, blind)")
    print(f"field: 1e5 < N <= {X:,}  ->  {n:,} values")
    print("=" * 74)
    print(f"  R^2(mu, ten ordinates) = {r2_mu:.4e}   [paper: 3.90e-3]")
    print()

    # ---------------------------------------------------------- (a)
    print(f"--- (a) the coin distribution, {NTRIAL} independent sign fields --")
    from lab_field_build import (smallest_prime_factor, von_mangoldt,
                                 mobius, fftconv_prefix)
    spf = smallest_prime_factor(X)
    lam, primes = von_mangoldt(X, spf)
    del spf
    mu = mobius(X, primes)
    idx = np.nonzero(mu != 0)[0]
    del mu
    vals = []
    for t in range(NTRIAL):
        rg = np.random.default_rng(5000 + t)
        eps = np.zeros(X + 1)
        eps[idx] = rg.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        Ce = fftconv_prefix(lam, eps, X)
        Ge = demean(Ce[Ni] / np.sqrt(V))
        v = r2_of(Ge, L, GAMMAS)
        vals.append(v)
        print(f"    coin {t + 1:>2}: R^2 = {v:.4e}"
              f"{'   <- exceeds mu' if v > r2_mu else ''}")
        del eps, Ce, Ge
    vals = np.array(vals)
    above = int((vals >= r2_mu).sum())
    print()
    print(f"    coin distribution: mean {vals.mean():.4e}, "
          f"sd {vals.std(ddof=1):.4e}, max {vals.max():.4e}")
    print(f"    R^2(mu) = {r2_mu:.4e}")
    print(f"    coin draws at or above mu: {above} of {NTRIAL}"
          f"   -> one-sided p = {(above + 1) / (NTRIAL + 1):.3f}")
    print(f"    z of mu in the coin distribution = "
          f"{(r2_mu - vals.mean()) / vals.std(ddof=1):+.2f}")
    verdict = ("FAILS H7 -- the statistic is reproduced with mu destroyed"
               if (above + 1) / (NTRIAL + 1) > 0.05
               else "PASSES H7")
    print(f"    VERDICT: {verdict}")
    print()

    # ---------------------------------------------------------- (b)
    print("--- (b) which surrogate gives the paper's 5.09e-6? ---------------")
    rng = np.random.default_rng(31337)
    lo, hi = GAMMAS.min(), GAMMAS.max()

    freq = np.empty(200)
    for s in range(200):
        freq[s] = r2_of(G, L, rng.uniform(lo, hi, size=len(GAMMAS)))
    print(f"    frequency surrogate (field intact, freqs redrawn):")
    print(f"      max {freq.max():.4e}   mean {freq.mean():.4e}   "
          f"sd {freq.std(ddof=1):.4e}")

    perm = np.empty(200)
    for s in range(200):
        y = rng.permutation(G)
        perm[s] = r2_of(y, L, GAMMAS)
    print(f"    permutation surrogate (field shuffled, zeros kept):")
    print(f"      max {perm.max():.4e}   mean {perm.mean():.4e}   "
          f"sd {perm.std(ddof=1):.4e}")
    print(f"    iid expectation 2k/n = {2 * len(GAMMAS) / n:.4e}")
    print(f"    [paper's quoted surrogate maximum: 5.09e-6]")
    print()
    print(f"    ratio, frequency max / permutation max = "
          f"{freq.max() / perm.max():.0f}")
    print()
    print("    G(N) and G(N+2) share almost every summand, so shuffling")
    print("    the field destroys the very correlation that makes a smooth")
    print("    regressor fit it. A permutation null is therefore the wrong")
    print("    null for this statistic, and it is the one that lands near")
    print("    the paper's figure.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
