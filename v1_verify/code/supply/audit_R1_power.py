# -*- coding: utf-8 -*-
"""
Re-verification of kill-test R1 of v1/paper/wall_v1.tex, and of the
precision it claims.

THE STATEMENT UNDER TEST (§7.2):

    R1 | zero-spectrum visibility (explicit formula) | **dead**:
    R^2_zeros = 0.2152 vs random-frequency 0.2196 +- 0.0055. In
    standard errors of its own null the measurement sits -0.80 below,
    and the data exclude a 7.5% enhancement at three standard errors.

THE DESIGN, from its own docstring: least-squares projection of
{D(k)} onto span{Re T_gamma, Im T_gamma : first 30 zeros} over 300 k,
against the same-dimension projection onto random-frequency templates
(6 draws of 30 frequencies uniform in [10,105]), with
T_gamma(k) = sum_m m^{i gamma} mu(N-mk). ALIVE iff
R^2_zeros >= 2 x mean(R^2_random).

THE OBJECTION BEING TESTED -- what the 7.5% is a percentage of.
Thirty zeros give sixty regressors on three hundred points, so the
chance R^2 is 60/300 = 0.20; the design's own code prints that
number. Both the measurement (0.2152) and the null (0.2196) therefore
sit within about 0.02 of a floor that is pure free parameters. The
GENUINE explanatory power of either template set is R^2 - 0.20, i.e.
about 0.015 and 0.020.

Three standard errors is 3 x 0.0055 = 0.0165 in absolute R^2. Against
0.2196 that is 7.5%, which is what the paper quotes. Against the
genuine component 0.020 it is 82%. So the test excludes the zeros
contributing more than roughly TWICE what random frequencies
contribute, and nothing finer. The comparison itself is sound -- both
template sets carry the same 60 parameters, so the floor cancels --
but the quoted precision is a percentage of a quantity that is 91%
free parameters.

Two further points the design does not carry: the null is estimated
from SIX draws, which is smaller than the eight the paper's own
Methodology names as a trap; and the adjusted R^2, which prices the
free parameters explicitly, is never reported.

METHOD HERE. The same design at a size where many null draws are
affordable, with the templates built as one matrix product per k.
Reported: raw R^2, adjusted R^2, the chance floor, the null at 6 and
at 40 draws, and the excludable effect expressed both ways.

PRE-REGISTRATION (written before the run).

  (1) RULE. The DEAD verdict needs R^2_zeros < 2 x mean(R^2_random),
      and with both sitting near the chance floor that is not in
      doubt. What is tested is the precision: report 3 sd as a
      fraction of R^2_random and as a fraction of
      R^2_random - p/n. If the two differ by an order of magnitude,
      the quoted 7.5% is not the test's resolution.
  (2) RULE. Report the null at 6 draws and at 40. If the 6-draw sd
      differs from the 40-draw sd by more than 30%, the quoted
      +-0.0055 is not a measurement of the spread.
  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      Both R^2 will sit just above p/n and their difference will be
      within a standard error of zero, so DEAD stands. I predict the
      3-sd bound is of the order of the whole genuine signal, so the
      test excludes only a doubling, and I predict the 6-draw sd is
      noticeably off the 40-draw one. I expect no change of verdict.
"""
import sys
import math
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
         37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
         52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
         67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
         79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
         92.491899, 94.651344, 95.870634, 98.831194, 101.317851]


def sieve_mu(X):
    mu = np.ones(X + 1, dtype=np.int8)
    comp = np.zeros(X + 1, dtype=bool)
    for p in range(2, int(X ** 0.5) + 1):
        if not comp[p]:
            comp[p * p::p] = True
            mu[p::p] = -mu[p::p]
            mu[p * p::p * p] = 0
    rest = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if not comp[p]:
            q = p
            while q <= X:
                rest[q::q] //= p
                q *= p
    mu[rest > 1] = -mu[rest > 1]
    mu[0] = 0
    del comp, rest
    return mu


def main():
    t0 = time.time()
    N = 40_000_000
    K0, NK = 2000, 300
    ks = np.arange(K0, K0 + NK, dtype=np.int64)
    mu = sieve_mu(N)
    SQ = int(N ** 0.5)
    print(f"mu ready t={time.time()-t0:.0f}s", flush=True)

    # cache per k: log m and the weight w = mu(N - m k)
    cache = []
    y = np.zeros(NK)
    for i, k in enumerate(ks):
        ms = np.arange(SQ + 1, N // int(k) + 1, dtype=np.int64)
        w = mu[N - int(k) * ms].astype(np.float64)
        y[i] = float(np.dot(mu[ms].astype(np.float64), w))
        cache.append((np.log(ms.astype(np.float64)), w))
    M = len(cache[0][0])
    print(f"y and cache ready, M = {M} terms per k, "
          f"t={time.time()-t0:.0f}s", flush=True)

    def design(freqs):
        X = np.empty((NK, 2 * len(freqs)))
        f = np.asarray(freqs)[:, None]
        for i, (lm, w) in enumerate(cache):
            ph = f * lm[None, :]
            X[i, 0::2] = np.cos(ph) @ w
            X[i, 1::2] = np.sin(ph) @ w
        return X

    def r2_of(X):
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ c
        return 1.0 - float(r @ r) / float(((y - y.mean()) ** 2).sum())

    p = 2 * len(ZEROS)
    Rz = r2_of(design(ZEROS))
    print(f"R^2_zeros = {Rz:.4f}   (v1 quotes 0.2152)", flush=True)

    rng = np.random.default_rng(20260902)
    rr = []
    for d in range(40):
        rr.append(r2_of(design(rng.uniform(10, 105, size=len(ZEROS)))))
        if d % 10 == 9:
            print(f"  null draw {d+1}/40  t={time.time()-t0:.0f}s",
                  flush=True)
    rr = np.array(rr)
    m6, s6 = float(rr[:6].mean()), float(rr[:6].std())
    m40, s40 = float(rr.mean()), float(rr.std())

    floor = p / NK
    adj = lambda r: 1 - (1 - r) * (NK - 1) / (NK - p - 1)
    print()
    print(f"  regressors p = {p}, points n = {NK}, "
          f"chance floor p/n = {floor:.4f}")
    print(f"  R^2_zeros            {Rz:.4f}   adjusted {adj(Rz):+.4f}")
    print(f"  R^2_random,  6 draws {m6:.4f} +- {s6:.4f}"
          f"   (v1 quotes 0.2196 +- 0.0055)")
    print(f"  R^2_random, 40 draws {m40:.4f} +- {s40:.4f}")
    print(f"  6-draw sd / 40-draw sd = {s6/max(s40,1e-12):.2f}")
    print()
    print(f"(1) the verdict: ALIVE needs R^2_zeros >= 2 x "
          f"{m40:.4f} = {2*m40:.4f}; measured {Rz:.4f} -> "
          f"{'ALIVE' if Rz >= 2*m40 else 'DEAD'}")
    print(f"    z of the measurement against the 40-draw null: "
          f"{(Rz-m40)/s40:+.2f}   (v1 quotes -0.80)")
    print()
    print("(2) what three standard errors excludes:")
    print(f"      3 sd                          = {3*s40:.4f} in R^2")
    print(f"      as a fraction of R^2_random   = "
          f"{3*s40/m40:.1%}   <- the form the paper quotes")
    print(f"      as a fraction of the genuine")
    print(f"        component R^2_random - p/n  = "
          f"{3*s40/max(m40-floor,1e-9):.1%}")
    print()
    print("    Both template sets carry the same 60 free parameters,")
    print("    so the floor cancels from the comparison and the DEAD")
    print("    verdict is sound. It does not cancel from the quoted")
    print("    precision: a percentage of R^2 is a percentage of a")
    print("    quantity that is mostly free parameters.")
    print("DONE")


if __name__ == "__main__":
    main()
