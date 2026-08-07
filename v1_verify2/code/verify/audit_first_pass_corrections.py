# -*- coding: utf-8 -*-
"""
Phase 2: the first pass's own corrections, independently recomputed.
(v1_verify2.)

`v1_verify/paper/wall_v1_corrected.tex` rewrites nine passages of the
paper and puts new numbers in them. The first pass is the only witness
to its own repairs, which is exactly the position `v1` was in before
`v1_verify` existed. This script re-derives the ones that are
measurable here, from the statements, using the machinery this tree
already built and validated (`prop:V` reproduces to six figures).

CHECKED HERE

  (a) sec:coin's replacement reconstruction: -0.124 and -0.052 under the
      two normalisations of S(h), and the claim that v1's -0.0976 comes
      from dividing by W(X) instead of V(X).
  (b) Proposition 15's replacement: Gamma = (sum_{h!=0} c(h))/V
      ~ N/(A(N) log N), quoted as 1.549e3 / 1.852e4 / 3.580e5 at
      N = 1e4 / 1.6e5 / 4e6 with Gamma*logN/N = 1.4266 / 1.3868 / 1.3605.
  (c) conj:wall item 3's replacement tail counts 21441/21463, 502/503.6,
      4/4.6.  The expected values pin the field size at 7.95e6, which is
      exactly this tree's field (1e5 < N <= 1.6e7). So the OBSERVED
      counts are directly comparable, and this tree measured a heavy
      tail on that field under the paper's stated cell index. This is
      the one place the two passes contradict each other, so the cell
      index is swept until the first pass's counts are reproduced.

NOT CHECKED HERE: K1's 0.904, R2's G_1 = 1.513, R1's six-draw spread
(all outside this tree's coverage), sec:R4's sign-randomisation z, and
sec:c3's 0.933..0.977 (the first pass imposes mu(m_i) != 0 inside A_j,
which this tree did not; noted in FINDINGS.md rather than re-run).

PRE-REGISTRATION.

  Decision rule: recompute each figure independently; CONFIRMED if it
  lands inside the precision quoted, REFUTED otherwise. For (c),
  additionally report which cell index reproduces the first pass's
  counts, since the paper's stated index does not.

  Prediction: (a) and (b) CONFIRMED -- both are deterministic sums with
  no threshold and no null, and the first pass's diagnosis of the
  W-vs-V swap is checkable arithmetic. (c) the counts will be
  reproduced only by a cell index wider than the paper's stated
  {3,5,7,11,13}, which would make the correction inherit the same
  undefined-statistic problem it was written to fix.
"""

import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np
from math import erfc, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CACHE = os.path.join(ROOT, "v1_verify2_log", "cache")
sys.path.insert(0, os.path.join(ROOT, "v1_verify2", "code", "wall"))


def phi(z):
    return 0.5 * erfc(-z / sqrt(2.0))


def autocorr(a, X):
    n = 1
    while n < 2 * len(a):
        n <<= 1
    f = np.fft.rfft(a, n)
    f *= np.conjugate(f)
    return np.fft.irfft(f, n)[:X].copy()


def main():
    from lab_field_build import smallest_prime_factor, von_mangoldt, mobius

    print("audit_first_pass_corrections   (v1_verify2 Phase 2)")
    print("=" * 76)

    # =============================================== (a) and (b)
    print()
    print("--- (a) sec:coin's replacement reconstruction, X = 4e6 ----------")
    X = 4_000_000
    spf = smallest_prime_factor(X)
    lam, primes = von_mangoldt(X, spf)
    del spf
    mu = mobius(X, primes)
    mu2f = (mu != 0).astype(np.float64)

    M = autocorr(mu.astype(np.float64), X)          # M(h)
    P = autocorr(lam, X)                            # c(h)
    W = float(P[0])
    # V(X) = sum_v mu^2(v) Lambda(X-v)^2, the paper's prop:V object
    lam2 = lam * lam
    V = float(np.dot(mu2f[1:X], lam2[X - 1:0:-1]))
    print(f"    W(X) = sum_w Lambda(w)^2                = {W:.6e}"
          f"   [first pass: 5.67847e7]")
    print(f"    V(X) = sum_v mu^2(v) Lambda(X-v)^2      = {V:.6e}"
          f"   [first pass: 4.46842e7]")
    # A(N) = prod_{q not| N} (1-1/(q(q-1))); X = 4e6 = 2^8 * 5^6, so the
    # factors at BOTH q=2 and q=5 come out of the Artin constant.
    A_of_X = 0.3739558136
    for q in (2, 5):
        A_of_X /= (1.0 - 1.0 / (q * (q - 1.0)))
    print(f"    W/V = {W / V:.5f}    1/A(N) predicted = {1 / A_of_X:.5f}"
          f"   [first pass: 1.27080 vs 1.27021]")
    print()

    h = np.arange(1, X)
    cS = P[1:X]
    quoted = {("M/X", "W"): -0.09768, ("M/X", "V"): -0.12413,
              ("M/(X-h)", "W"): -0.04080, ("M/(X-h)", "V"): -0.05185}
    print(f"    {'S(h) reading':>12}{'denominator':>13}{'rho-1':>12}"
          f"{'first pass':>13}{'ratio to -0.18':>16}")
    for sname, Sh in (("M/X", M[1:X] / X),
                      ("M/(X-h)", M[1:X] / (X - h))):
        for dname, D in (("W", W), ("V", V)):
            val = 2.0 * float(np.dot(cS, Sh)) / D
            q = quoted[(sname, dname)]
            flag = "  <== CONFIRMED" if abs(val - q) < 5e-4 else "  <== DIFFERS"
            print(f"    {sname:>12}{dname:>13}{val:>12.5f}{q:>13.5f}"
                  f"{val / -0.18:>16.3f}{flag}")
    print()
    print("    [the paper's own figure is -0.0976, 'a factor 0.54']")

    print()
    print("--- (b) Proposition 15's replacement: the amplification Gamma ---")
    print(f"    Gamma = (sum_{{h!=0}} c(h)) / V = (psi(X)^2 - W)/V")
    print(f"    {'N':>10}{'Gamma':>14}{'Gamma logN/N':>15}"
          f"{'first pass':>26}")
    fp = {10_000: (1.549e3, 1.4266), 160_000: (1.852e4, 1.3868),
          4_000_000: (3.580e5, 1.3605)}
    for Xi in (10_000, 160_000, 4_000_000):
        psi = float(lam[: Xi + 1].sum())
        Wi = float((lam[: Xi + 1] ** 2).sum())
        Vi = float(np.dot(mu2f[1:Xi], lam2[Xi - 1:0:-1]))
        G = (psi * psi - Wi) / Vi
        g, r = fp[Xi]
        print(f"    {Xi:>10,}{G:>14.4e}{G * np.log(Xi) / Xi:>15.4f}"
              f"{f'{g:.4g} / {r}':>26}")
    print("    (the consequence: S(h)=o(1) buys nothing; the requirement")
    print("     is S(h) = o(log N / N), a gap of a factor N/log N)")
    del M, P, lam2, mu, mu2f, lam

    # =============================================== (c)
    print()
    print("--- (c) conj:wall item 3's replacement tail counts ---------------")
    Xw = 16_000_000
    z = np.load(os.path.join(CACHE, f"field_{Xw}.npz"))
    good = (z["V"] > 0) & (z["W"] > 0)
    Ni = z["N"][good]
    Nf = Ni.astype(np.float64)
    fld = Nf > 1e5
    Ni = Ni[fld]
    Z = z["C"][good][fld] / np.sqrt(z["V"][good][fld])
    n = len(Z)
    print(f"    field 1e5 < N <= 1.6e7: n = {n:,}")
    for t, e in ((3, 21463), (4, 503.6), (5, 4.6)):
        print(f"      expected at t={t}: mine {2 * phi(-t) * n:>9.1f}   "
              f"first pass {e}")
    print("    -> the first pass's expected values pin the same field.")
    print()

    def counts(labels, k, standardise):
        u, inv = np.unique(labels, return_inverse=True)
        g = Z.copy()
        c = np.bincount(inv, minlength=len(u))
        s = np.bincount(inv, weights=Z, minlength=len(u))
        g -= (s / np.maximum(c, 1))[inv]
        if standardise:
            sq = np.bincount(inv, weights=g * g, minlength=len(u))
            sd = np.sqrt(sq / np.maximum(c - 1, 1))
            sd[sd == 0] = 1.0
            g = np.where(c[inv] > 5, g / sd[inv], g)
        gs = g / g.std(ddof=1)
        return [int((np.abs(gs) > t).sum()) for t in (3, 4, 5)], len(u)

    print(f"    {'cell index':<44}{'cells':>7}{'t=3':>8}{'t=4':>7}{'t=5':>6}")
    print(f"    {'first pass (corrected paper)':<44}{'?':>7}"
          f"{21441:>8}{502:>7}{4:>6}")
    print("    " + "-" * 65)
    sets = [((3, 5, 7, 11, 13), "as sec:floor states"),
            ((3, 5, 7, 11, 13, 17), ""),
            ((3, 5, 7, 11, 13, 17, 19), ""),
            ((3, 5, 7, 11, 13, 17, 19, 23), "'depth 8'"),
            ((2, 3, 5, 7, 11, 13, 17, 19), "'depth 8' incl. 2")]
    for ps, note in sets:
        lab = np.zeros(n, dtype=np.int64)
        for i, p in enumerate(ps):
            lab |= ((Ni % p == 0).astype(np.int64) << i)
        for st in (False, True):
            c, k = counts(lab, len(ps), st)
            tag = f"{{{','.join(map(str, ps))}}}"
            op = "standardised" if st else "means only  "
            print(f"    {tag:<28}{op:<16}{k:>7}"
                  f"{c[0]:>8}{c[1]:>7}{c[2]:>6}"
                  f"{('   ' + note) if note and not st else ''}")
    # The remaining candidate: the first pass warns that pooling octaves
    # is a scale mixture, so it may have normalised within each octave
    # before pooling. Both forms of that are tried.
    print()
    print("    treatments that also normalise by octave (the first pass")
    print("    warns that pooling octaves is a scale mixture):")
    octv = np.floor(np.log2(Ni / 1e5)).astype(np.int64)
    for ps in ((3, 5, 7, 11, 13), (3, 5, 7, 11, 13, 17, 19, 23)):
        base = np.zeros(n, dtype=np.int64)
        for i, p in enumerate(ps):
            base |= ((Ni % p == 0).astype(np.int64) << i)
        lab = base * 64 + octv
        for st in (False, True):
            c, k = counts(lab, 0, st)
            tag = f"{{{','.join(map(str, ps))}}} x octave"
            op = "standardised" if st else "means only  "
            print(f"    {tag:<28}{op:<16}{k:>7}"
                  f"{c[0]:>8}{c[1]:>7}{c[2]:>6}")
    print()
    print("    the paper states cells indexed by {3,5,7,11,13} and removal")
    print("    of cell MEANS alone; that row is the first one above.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
