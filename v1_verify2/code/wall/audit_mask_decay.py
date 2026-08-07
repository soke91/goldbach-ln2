# -*- coding: utf-8 -*-
"""
sec:floor's decay table: are the six exponents measurable?
(v1_verify2, Phase 1, blind.)

audit_cell_floor.py established, using the exact floor of lem:cellmom
(itself confirmed against 60 independent-sign draws):

  - max_c |z_c| = 9.1 in the top band, Bonferroni cleared in all eight
    octaves   [paper: 8.4, cleared]                        REPRODUCED
  - prop:coh's b = 0.0394..0.0397 at the three shallowest depths
    [paper: 0.0379, 0.0378, 0.0379]                        REPRODUCED
  - but the decay exponents' standard errors came out 2x to 100x LARGER
    than the paper's, and in the top band delta_c at depths 0, 1, 2 is
    +0.0, +0.7, -0.1 standard errors -- i.e. no detectable mask at all.

The paper's table gives depth 0 an exponent of 0.6289 +- 0.0121, a
significance of 52.0, and it is one of the four steps that carry the
"rises monotonically as the cell gets shallower" claim. If depth 0's
amplitude is consistent with zero at every scale, that exponent is
fitted to noise.

This script prints the whole delta_c +- se_c grid so the question can be
answered by looking, and refits the exponents under BOTH cell readings,
since sec:floor says "cells are indexed by depth" while conj:wall says
"indexed by which small primes divide N":

   (i) six depth cells
  (ii) the 32 divisibility patterns, pooled by depth in the fit
       (more points per depth, hence smaller fitted errors)

PRE-REGISTRATION.

  Decision rule. For each depth report, at every octave, delta_c, the
  exact se from lem:cellmom, and z. Then fit log|delta| on log N,
  weighted by the exact errors, with the exponent's se from the fit
  covariance -- the procedure sec:floor states.
    MEASURABLE   : |a_d / s.e.| >= 5 under at least one cell reading.
    NOT MEASURABLE : it is not.
  Report how many octaves have |z_c| >= 3 for each depth; an exponent
  fitted to amplitudes that are everywhere within noise is not a
  measurement whatever its nominal error.

  Prediction written before running. I predict depths 3, 4, 5 are
  measurable under both readings and depths 0, 1, 2 are not under
  either, because their amplitudes sit at or below the exact floor in
  the top band and the floor falls only like (log N)^{-1/2} while the
  paper's own fitted a_d would have the amplitude falling like N^{-0.4}
  or faster -- so the signal-to-floor ratio gets WORSE, not better,
  toward large N, and the shallow cells cannot have been detectable at
  the top of the range.
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


def next_pow2(x):
    n = 1
    while n < x:
        n <<= 1
    return n


def wls_slope(lx, y, sy):
    w = 1.0 / sy ** 2
    Sw, Sx = w.sum(), (w * lx).sum()
    Sxx, Sy, Sxy = (w * lx * lx).sum(), (w * y).sum(), (w * lx * y).sum()
    den = Sw * Sxx - Sx * Sx
    return (Sw * Sxy - Sx * Sy) / den, np.sqrt(Sw / den)


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 16_000_000
    z = np.load(os.path.join(CACHE, f"field_{X}.npz"))
    goodm = (z["V"] > 0) & (z["W"] > 0)
    Ni = z["N"][goodm]
    Vv = z["V"][goodm]
    Z = z["C"][goodm] / np.sqrt(Vv)
    depth = z["depth"][goodm]
    patt = z["cell"][goodm]

    from lab_field_build import smallest_prime_factor, von_mangoldt, mobius
    print("audit_mask_decay   (v1_verify2 Phase 1, blind)")
    print("  sieving ...")
    spf = smallest_prime_factor(X)
    lam, primes = von_mangoldt(X, spf)
    del spf
    mu = mobius(X, primes)
    mu2 = (mu != 0).astype(np.float64)
    del mu
    print("=" * 78)

    bands = []
    b = X
    while b > 6e4:
        bands.append((b / 2.0, b))
        b /= 2.0
    bands = bands[::-1]

    def analyse(lo, hi, labels, k):
        sel = (Ni > lo) & (Ni <= hi)
        Nb, Zb, Vb, lb = Ni[sel], Z[sel], Vv[sel], labels[sel]
        n = len(Nb)
        B = int(hi)
        nfft = next_pow2(2 * B + 4)
        L = np.zeros(B + 1)
        L[1:] = lam[1: B + 1]
        fl = np.conjugate(np.fft.rfft(L, nfft))
        del L
        m2 = mu2[: B + 1]

        def corr(w):
            fw = np.fft.rfft(w, nfft)
            fw *= fl
            r = np.fft.irfft(fw, nfft)[: B + 1].copy()
            del fw
            return r

        wa = np.zeros(B + 1)
        wa[Nb] = 1.0 / np.sqrt(Vb)
        ua = corr(wa)
        Qaa = float(np.dot(m2 * ua, ua))
        res = {}
        for c in range(k):
            m = lb == c
            nc = int(m.sum())
            if nc < 3:
                continue
            wc = np.zeros(B + 1)
            wc[Nb[m]] = 1.0 / np.sqrt(Vb[m])
            uc = corr(wc)
            Qcc = float(np.dot(m2 * uc, uc))
            Qca = float(np.dot(m2 * uc, ua))
            var = Qcc / nc ** 2 - 2 * Qca / (nc * n) + Qaa / n ** 2
            res[c] = (nc, float(Zb[m].mean() - Zb.mean()),
                      np.sqrt(max(var, 1e-300)))
            del wc, uc
        del wa, ua, fl
        return res

    # ---------------------------------------------- (i) depth cells
    print()
    print("--- the delta_c +- se_c grid, cells = the six depths -------------")
    grid = {}
    for lo, hi in bands:
        grid[hi] = analyse(lo, hi, depth, 6)
    print(f"    {'octave top':>12}", end="")
    for d in range(6):
        print(f"{'d=' + str(d):>16}", end="")
    print()
    for lo, hi in bands:
        print(f"    {hi:>12,.0f}", end="")
        for d in range(6):
            r = grid[hi].get(d)
            if r is None:
                print(f"{'--':>16}", end="")
            else:
                print(f"{r[1]:>+8.4f}/{r[1] / r[2]:>6.1f}", end="")
        print()
    print("    (each entry is delta_c / z_c, z against the exact floor)")
    print()
    print(f"    {'depth':>6}{'octaves with |z|>=3':>22}{'of':>4}")
    for d in range(6):
        cnt = sum(1 for _, hi in bands
                  if d in grid[hi] and abs(grid[hi][d][1] / grid[hi][d][2]) >= 3)
        print(f"    {d:>6}{cnt:>22}{len(bands):>4}")

    # ---------------------------------------------- (ii) 32 patterns
    print()
    print("--- refit with the 32 divisibility patterns, pooled by depth ----")
    grid32 = {}
    for lo, hi in bands:
        grid32[hi] = analyse(lo, hi, patt, 32)
    pd = {c: bin(c).count("1") for c in range(32)}

    pa = {5: (0.1434, 0.0155), 4: (0.2152, 0.0065), 3: (0.2713, 0.0040),
          2: (0.3686, 0.0052), 1: (0.0437, 0.0556), 0: (0.6289, 0.0121)}
    print()
    print(f"    {'depth':>6}{'a_d (6 cells)':>15}{'s.e.':>9}{'a/se':>7}"
          f"{'a_d (32 pooled)':>17}{'s.e.':>9}{'a/se':>7}"
          f"{'paper a_d':>11}{'paper s.e.':>11}{'paper a/se':>11}")
    for d in (5, 4, 3, 2, 1, 0):
        row = [d]
        for src, k in ((grid, 6), (grid32, 32)):
            lx, y, sy = [], [], []
            for lo, hi in bands:
                if k == 6:
                    items = [(d, src[hi][d])] if d in src[hi] else []
                else:
                    items = [(c, src[hi][c]) for c in src[hi] if pd[c] == d]
                for _, (nc, dl, se) in items:
                    if abs(dl) <= 0:
                        continue
                    lx.append(np.log(np.sqrt(lo * hi)))
                    y.append(np.log(abs(dl)))
                    sy.append(se / abs(dl))
            if len(lx) < 3:
                row += [np.nan, np.nan]
                continue
            s, ss = wls_slope(np.array(lx), np.array(y), np.array(sy))
            row += [-s, ss]
        p = pa[d]
        print(f"    {d:>6}{row[1]:>15.4f}{row[2]:>9.4f}"
              f"{abs(row[1] / row[2]):>7.1f}"
              f"{row[3]:>17.4f}{row[4]:>9.4f}{abs(row[3] / row[4]):>7.1f}"
              f"{p[0]:>11.4f}{p[1]:>11.4f}{p[0] / p[1]:>11.1f}")
    print()
    print("    MEASURABLE was pre-registered as |a_d/s.e.| >= 5 under at")
    print("    least one reading.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
