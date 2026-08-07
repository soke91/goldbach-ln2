# -*- coding: utf-8 -*-
"""
conj:wall items 1 and 3: what exactly is Gaussian, and on what field?
(v1_verify2, Phase 1, blind.)

audit_propV_and_wall.py reproduced the S(N)N-scaled excess kurtosis
(+0.182 here vs +0.1704 quoted) and the S-candidate residual sd
(0.245236 here vs 0.245235 quoted, six figures) -- so the pipeline
agrees with the paper's on the quantities that do reproduce. It did NOT
reproduce:

  excess kurtosis of G      quoted -0.0005 (z=-0.3), measured +0.0252
  tail ratio at t=4         quoted  0.997,           measured  1.60
  tail ratio at t=5         quoted  0.878,           measured 12.9
  extremes at generic N     quoted,                  measured all 20
                                                     largest |G| at
                                                     depth >= 3

This script asks which of the three under-specified choices accounts for
the gap: (i) the field's lower cutoff, (ii) what "cell" means, (iii)
what "removing cell means alone" removes.

PRE-REGISTRATION (fixed before this ran).

  Decision rule.  Enumerate the cross product of the three choices and
  report excess kurtosis and the t=3,4,5 tail ratios for each.
    EXPLAINED   : some cell in the cross product returns the paper's
                  -0.0005 and 0.999/0.997/0.878 together. Then the
                  finding is under-specification, not error, and the
                  reading that works is named.
    NOT EXPLAINED : no cell does. Then the measurement is not
                  reproducible from the paper's description.

  Prediction written before running.  EXPLAINED, by per-cell
  standardisation: I expect that dividing by each cell's own sd (rather
  than removing its mean alone) collapses the kurtosis to ~0 and the
  tails to ~1. If so, the operative finding is that the paper's stated
  operation ("removing cell means alone", under the scale sqrt(V) "and
  no other") is not the operation that produces its numbers, and that a
  second, per-cell scale is doing the work.

  Falsifier for that prediction: per-cell standardisation leaves the
  kurtosis well away from zero, in which case the cutoff is the cause.

  Cross-check that the pipeline is sound: rho = mean(C^2)/mean(V) must
  come out near the paper's 0.82 at N ~ 1.4e7 (sec:coin, rem:rho). If it
  does not, everything above is my bug and the script says so.
"""

import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np
from math import erfc, sqrt


def phi(z):
    return 0.5 * erfc(-z / sqrt(2.0))


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CACHE = os.path.join(ROOT, "v1_verify2_log", "cache")


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 16_000_000
    z = np.load(os.path.join(CACHE, f"field_{X}.npz"))
    good = (z["V"] > 0) & (z["W"] > 0)
    Ni = z["N"][good]
    N = Ni.astype(np.float64)
    C = z["C"][good]
    V = z["V"][good]
    cell = z["cell"][good]
    depth = z["depth"][good]

    print("audit_wall_gaussianity  (v1_verify2 Phase 1, blind)")
    print("=" * 74)

    # ---- pipeline cross-check: rho ---------------------------------
    print()
    print("--- cross-check: rho = mean(C^2)/mean(V), the paper's first "
          "convention ---")
    print("    (sec:coin measures rho-1 = -0.18; rem:rho quotes 0.841 in "
          "these units")
    print("     at N~1e8, and says the three conventions agree to 0.75% at "
          "N~1.4e7)")
    for lo, hi in ((1e5, 1.4e7), (7e6, 1.4e7), (8e6, 1.6e7), (1e6, 2e6)):
        m = (N > lo) & (N <= hi)
        r1 = (C[m] ** 2).mean() / V[m].mean()
        r2 = ((C[m] ** 2) / V[m]).mean()
        r3 = (np.pi / 2) * (np.abs(C[m]) / np.sqrt(V[m])).mean() ** 2
        print(f"    N in ({lo:.3g}, {hi:.3g}]  n={m.sum():>9,}   "
              f"mean(C^2)/mean(V)={r1:.4f}   mean(C^2/V)={r2:.4f}   "
              f"(pi/2)mean(|C|/sqrtV)^2={r3:.4f}")
        print(f"        spread across the three conventions: "
              f"{(max(r1, r2, r3) / min(r1, r2, r3) - 1) * 100:.2f}%")
    print()

    Z = C / np.sqrt(V)

    def variants(mask, label):
        sub = Z[mask]
        cl = cell[mask]
        dp = depth[mask]
        n = len(sub)
        out = []
        for nm, lab, k, scale in (
            ("raw Z (nothing removed)", None, 0, False),
            ("minus cell mean (32 patterns)", cl, 32, False),
            ("minus cell mean (6 depths)", dp, 6, False),
            ("per-cell standardised (32)", cl, 32, True),
            ("per-cell standardised (6 depths)", dp, 6, True),
        ):
            g = sub.copy()
            if lab is None:
                g = g - g.mean()
            else:
                for c in range(k):
                    mm = lab == c
                    if mm.sum() > 2:
                        g[mm] -= g[mm].mean()
                        if scale:
                            g[mm] /= g[mm].std(ddof=1)
            gs = g / g.std(ddof=1)
            ek = (gs ** 4).mean() - 3.0
            am = np.abs(gs).mean()
            row = [nm, n, ek, ek / np.sqrt(24.0 / n), am - sqrt(2 / np.pi)]
            for t in (3, 4, 5):
                obs = int((np.abs(gs) > t).sum())
                exp = 2 * phi(-t) * n
                row.append(obs / exp)
            out.append(row)
        print(f"--- field: {label}   (n = {n:,}) " + "-" * 12)
        print(f"    {'treatment':<34}{'exc kurt':>10}{'z':>9}"
              f"{'E|G| dev':>11}{'t=3':>8}{'t=4':>8}{'t=5':>9}")
        for r in out:
            print(f"    {r[0]:<34}{r[2]:>+10.4f}{r[3]:>+9.1f}"
                  f"{r[4]:>+11.5f}{r[5]:>8.3f}{r[6]:>8.3f}{r[7]:>9.3f}")
        print(f"    [paper: -0.0005, z=-0.3, -0.00018, 0.999, 0.997, 0.878]")
        print()

    for lo, hi, lab in (
        (0, X, "every even N <= 1.6e7 (the paper's words)"),
        (1e5, 1.4e7, "1e5 <= N <= 1.4e7 (the range prop:coh uses)"),
        (1e6, X, "1e6 <= N <= 1.6e7"),
        (3.4e6, X, "3.4e6 <= N <= 1.6e7 (the cutoff giving 6.3e6 values)"),
        (X / 2, X, "top octave"),
    ):
        m = (N > lo) & (N <= hi)
        variants(m, lab)

    # ---- where does the variance live? -----------------------------
    print("--- the variance of Z is not constant across cells ---------------")
    print("    (item 1 removes cell MEANS; if the cells also differ in")
    print("     SPREAD, the pooled field is a mixture and its kurtosis is")
    print("     inflated by the mixture alone)")
    m = N > X / 2
    print(f"    top octave, {m.sum():,} values")
    print(f"    {'depth':>6}{'count':>10}{'mean Z':>10}{'sd Z':>10}"
          f"{'sd/pooled':>11}")
    pooled = Z[m].std(ddof=1)
    for d in range(6):
        mm = m & (depth == d)
        if mm.sum() > 2:
            print(f"    {d:>6}{mm.sum():>10,}{Z[mm].mean():>+10.4f}"
                  f"{Z[mm].std(ddof=1):>10.4f}"
                  f"{Z[mm].std(ddof=1) / pooled:>11.4f}")
    print()
    print("    a two-component normal mixture with sd ratio r and weights")
    print("    (w, 1-w) has excess kurtosis > 0 for any r != 1; that is the")
    print("    mechanism, and removing means does not touch it.")
    print()

    # ---- who are the extremes? -------------------------------------
    print("--- item 3: 'the extremes are attained at generic N, not at deep "
          "radicals' ---")
    for lo, hi, lab in ((0, X, "every even N"),
                        (1e5, 1.4e7, "1e5..1.4e7"),
                        (X / 2, X, "top octave")):
        m = (N > lo) & (N <= hi)
        g = Z[m].copy()
        cl = cell[m]
        for c in range(32):
            mm = cl == c
            if mm.sum() > 2:
                g[mm] -= g[mm].mean()
        gs = g / g.std(ddof=1)
        dsub = depth[m]
        idx = np.argsort(-np.abs(gs))[:50]
        cnt = np.bincount(dsub[idx], minlength=6)
        base = np.bincount(dsub, minlength=6) / len(dsub)
        exp = (base * 50).round(1)
        print(f"    {lab:<14} depth of 50 largest |G|: {cnt.tolist()}")
        print(f"    {'':<14} expected if generic     : {exp.tolist()}")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
