# -*- coding: utf-8 -*-
"""
Re-verification of kill-test R2 of v1/paper/wall_v1.tex, and of whether
its reported numbers determine its verdict.

THE STATEMENT UNDER TEST (§7.2):

    R2 | determinant / Kloosterman phase | **dead**: regression
    R^2 = -0.0001/+0.0004 against controls +-0.0002 over 48,000
    coprime pairs. The verdict rests on the measurement (-0.38 s.e.),
    not on the criterion.

THE PRE-REGISTRATION, from the design's own docstring:

    "Pre-registered: ALIVE iff mean G >= 2 or regression capture >= 2x
     random control; else DEAD."

THE IMPLEMENTED CRITERION, from the same file:

    alive = (Gs.mean() >= 2) or (Gs_m.mean() >= 2) or \\
            (mr > 0 and R2d >= 2*mr)

with `mr` the mean random-phase control R^2. The guard `mr > 0` is not
in the pre-registration, and it matters in both directions:

  * if mr <= 0 the regression arm CANNOT fire, whatever R2d is;
  * if mr > 0 but near zero -- and the paper reports the controls at
    +-0.0002, i.e. near zero -- then "2 x mr" is not a bar, and
    R2d = +0.0004 clears it for any mr < 0.0002.

The paper reports R2d = +0.0004 on one of the two runs and the control
spread as +-0.0002. **Those numbers do not determine the verdict**:
whether the regression arm fired depends on the sign and size of mr,
which is not reported. This script computes it.

METHOD HERE. Written from the statement plus the design's own
description: for 40 base k and 600 coprime partners k' each,

    C~(k,k') = sum_{p ~ P} mu(N-pk) mu(N-pk') / sqrt(#nonzero),
    phi_1 = 2 pi (N kbar' mod k)/k,  phi_2 = 2 pi (N kbar mod k')/k',

regress C~ on cos/sin of both phases. The random-phase control is run
with 50 draws rather than 8, and two further nulls are added that the
design does not carry: the adjusted R^2, and a permutation of C~
against a fixed design.

PRE-REGISTRATION (written before the run).

  (1) RULE. Report mr with its sign. If mr <= 0, the regression arm of
      R2's criterion was unreachable and the DEAD verdict rests on the
      coherent-gain arm alone. If 0 < mr < R2d/2, the criterion as
      IMPLEMENTED says ALIVE, and the reported verdict disagrees with
      the reported criterion.
  (2) RULE, on the substance rather than the bookkeeping. The honest
      test is R2d against the control distribution in standard errors.
      Report that z. |z| < 3 confirms DEAD on the measurement, which
      is what the paper says the verdict rests on.
  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      The phases are determined by N mod k and modular inverses, and
      Conjecture 10 says the field is a mask times a featureless
      fluctuation, so I expect R2d within a couple of standard errors
      of the control and DEAD to stand on the measurement. I expect mr
      to be small and of either sign, and therefore the CRITERION to
      be undetermined by the published numbers -- which is the
      finding, independent of the verdict.
  (4) The paper says "48,000 coprime pairs". The design is 40 base k
      times 600 partners = 24,000. Reported.
"""
import sys
import time
from math import gcd

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


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


def primes_in(a, b):
    s = np.ones(b + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(b ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    idx = np.nonzero(s)[0]
    return idx[(idx >= a) & (idx <= b)].astype(np.int64)


def main():
    t0 = time.time()
    N = 199_999_998
    K0, K1 = 2000, 4000
    NK, NPART = 40, 600
    rng = np.random.default_rng(20260903)

    mu = sieve_mu(N)
    print(f"mu ready t={time.time()-t0:.0f}s", flush=True)
    P0 = N // (2 * K1)
    ps = primes_in(P0, 2 * P0)
    print(f"primes p in [{P0}, {2*P0}] : {len(ps)}", flush=True)

    kvals = rng.choice(np.arange(K0, K1), size=NK, replace=False)
    allC, allp1, allp2 = [], [], []
    G1s, G2s, Grs = [], [], []
    for ik, k in enumerate(kvals):
        k = int(k)
        kps = [x for x in range(K0, K1) if x != k and gcd(x, k) == 1]
        kps = kps[:NPART]
        C = np.full(len(kps), np.nan)
        p1 = np.zeros(len(kps))
        p2 = np.zeros(len(kps))
        for i, kp in enumerate(kps):
            pmax = (N - 2) // max(k, kp)
            pp = ps[ps <= pmax]
            t = (mu[N - pp * k].astype(np.int64)
                 * mu[N - pp * kp].astype(np.int64))
            nz = int(np.count_nonzero(t))
            if nz > 100:
                C[i] = t.sum() / np.sqrt(nz)
            p1[i] = 2 * np.pi * ((N % k) * pow(kp, -1, k) % k) / k
            p2[i] = 2 * np.pi * ((N % kp) * pow(k, -1, kp) % kp) / kp
        ok = ~np.isnan(C)
        v, q1, q2 = C[ok], p1[ok], p2[ok]
        base = float((v ** 2).sum())
        if len(v) and base > 0:
            G1s.append(abs(np.sum(v * np.exp(-1j * q1))) ** 2 / base)
            G2s.append(abs(np.sum(v * np.exp(-1j * q2))) ** 2 / base)
            rp = rng.uniform(0, 2 * np.pi, size=len(v))
            Grs.append(abs(np.sum(v * np.exp(-1j * rp))) ** 2 / base)
        allC.append(v)
        allp1.append(q1)
        allp2.append(q2)
        if ik % 10 == 9:
            print(f"  k {ik+1}/{NK}  t={time.time()-t0:.0f}s",
                  flush=True)

    v = np.concatenate(allC)
    q1 = np.concatenate(allp1)
    q2 = np.concatenate(allp2)
    n = len(v)
    print()
    print(f"pairs used: {n}   (design is {NK} x {NPART} = "
          f"{NK*NPART}; the paper says 48,000)")
    print()

    def r2_of(X, y):
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ c
        return 1.0 - float(r @ r) / float(((y - y.mean()) ** 2).sum())

    Xd = np.column_stack([np.cos(q1), np.sin(q1),
                          np.cos(q2), np.sin(q2)])
    R2d = r2_of(Xd, v)
    ctrl = []
    for _ in range(50):
        rp = rng.uniform(0, 2 * np.pi, size=(n, 2))
        Xr = np.column_stack([np.cos(rp[:, 0]), np.sin(rp[:, 0]),
                              np.cos(rp[:, 1]), np.sin(rp[:, 1])])
        ctrl.append(r2_of(Xr, v))
    ctrl = np.array(ctrl)
    mr, sr = float(ctrl.mean()), float(ctrl.std())
    perm = np.array([r2_of(Xd, v[rng.permutation(n)])
                     for _ in range(50)])

    G1s, G2s, Grs = map(np.array, (G1s, G2s, Grs))
    nk = len(G1s)
    print(f"(a) coherent gain, over {nk} usable k of {NK}")
    print(f"    the design's bar is the ABSOLUTE value 2, stated")
    print(f"    against a theoretical null of 1. The null's own spread")
    print(f"    is not computed there; it is computed here.")
    # a further null: fresh random phases, many draws, same k-structure
    grand = []
    for _ in range(400):
        acc = []
        for vv in allC:
            if len(vv) and (vv ** 2).sum() > 0:
                ph = rng.uniform(0, 2 * np.pi, size=len(vv))
                acc.append(abs(np.sum(vv * np.exp(-1j * ph))) ** 2
                           / float((vv ** 2).sum()))
        grand.append(np.mean(acc))
    grand = np.array(grand)
    gm, gs = float(grand.mean()), float(grand.std())
    print(f"    random-phase null over {nk} k, 400 draws: "
          f"mean {gm:.3f}, sd {gs:.3f}")
    for nm, val in (("G1 (phase e(-N kbar'/k))", float(G1s.mean())),
                    ("G2 (mirror phase)", float(G2s.mean())),
                    ("the run's own random draw", float(Grs.mean()))):
        print(f"      {nm:>28} = {val:.3f}"
              f"   z = {(val-gm)/gs:+.2f}"
              f"   vs bar 2: {'clears' if val >= 2 else 'below'}")
    print(f"    the bar 2 sits at z = {(2-gm)/gs:+.2f} of this null.")
    print()
    print(f"(b) regression      R2_det   = {R2d:+.6f}")
    print(f"    random-phase control, 50 draws: mean mr = {mr:+.6f}, "
          f"sd = {sr:.6f}")
    print(f"    permutation control, 50 draws : mean "
          f"{perm.mean():+.6f}, sd {perm.std():.6f}")
    print(f"    adjusted R^2 (4 predictors, n = {n}) = "
          f"{1-(1-R2d)*(n-1)/(n-5):+.6f}")
    print(f"    pure-noise R^2 = 4/n = {4/n:.6f}")
    print()
    print(f"(1) the criterion, as IMPLEMENTED:")
    print(f"      mr > 0                : {mr > 0}")
    print(f"      R2d >= 2*mr           : {R2d >= 2*mr} "
          f"(2*mr = {2*mr:+.6f})")
    print(f"      regression arm fires  : {bool(mr > 0 and R2d >= 2*mr)}")
    print(f"    the criterion, as PRE-REGISTERED "
          f"('capture >= 2x random control', no sign guard):")
    print(f"      R2d >= 2*mr           : {R2d >= 2*mr}")
    print(f"    coherent-gain arms fire : "
          f"{bool(np.mean(G1s) >= 2 or np.mean(G2s) >= 2)}")
    print()
    print(f"(2) the measurement, which is what the paper says the")
    print(f"    verdict rests on: R2_det against the control")
    print(f"      z = {(R2d - mr)/sr:+.2f}"
          f"   (paper quotes -0.38 s.e.)")
    print()
    print("    A criterion of the form 'X >= 2 * control' is not a bar")
    print("    when the control is a quantity whose own mean is zero to")
    print("    within its spread. The verdict here stands on the z, as")
    print("    the paper says; the criterion it is written against does")
    print("    not decide anything.")
    print("DONE")


if __name__ == "__main__":
    main()
