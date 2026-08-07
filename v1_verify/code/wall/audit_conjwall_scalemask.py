# -*- coding: utf-8 -*-
"""
Re-verification of measurement 1 of Conjecture 14 (`conj:wall`) in
v1/paper/wall_v1.tex, and of the form of the conjecture itself.

THE STATEMENT UNDER TEST, verbatim:

    1. "The bulk is Gaussian, under this scale and no other.
        Excess kurtosis -0.0005 (z = -0.3) and E|G|/sd(G) short of
        sqrt(2/pi) by 0.00018 (z = -0.8), on 6.3e6 values, REMOVING
        CELL MEANS ALONE. Under an S N-based scale the same data give
        excess kurtosis +0.1704 at z = 98: a normaliser that
        manufactures a heavy tail."

and the conjecture the measurement supports,

    Conjecture 14:  C(N) = m(N) + sqrt(V(N)) G(N),
    with m(N) a deterministic LOCATION mask and G(N) Gaussian.

THE OBJECTION BEING TESTED. Under Conjecture 14, once the location
mask m(N) is removed and the exact scale sqrt(V(N)) divided out, what
is left is Gaussian. Removing cell MEANS is exactly "removing m(N)".
So item 1's protocol -- "removing cell means alone" -- is the protocol
the conjecture prescribes, and the near-zero kurtosis is the
conjecture's own prediction being met.

The question here is whether the two things item 1 attributes to
mean-removal and to a rival normaliser are what actually produce those
numbers. Specifically:

  (Q1) Does removing cell MEANS alone give excess kurtosis ~ 0, or
       does it leave a positive excess?
  (Q2) If a positive excess survives mean-removal, does it come from
       the cells having different VARIANCES -- a scale mask, which
       Conjecture 14 does not contain -- rather than different means?
  (Q3) Is the +0.1704 a property of an "S N-based scale", or of the
       V-based scale with means-only removal?

METHOD HERE. Written from the statement. Z = C/sqrt(V) with V the
exact second moment of Proposition 12. Cells keyed by which of the
first d odd primes divide N. Three conditions per octave band:

    RAW    Z itself
    MEAN   Z minus its cell mean          <- item 1's stated protocol
    STD    (Z - cell mean) / cell sd      <- also removes a scale mask

A mixture of Gaussians with a spread of variances has excess kurtosis
3 * Var(s^2) / (E s^2)^2 exactly. That number is computed from the
measured per-cell variances and compared with the measured excess, so
(Q2) is answered quantitatively and not by inference from a collapse.

Cells with few members give noisy sd estimates, and dividing by a
noisy sd shrinks tails by construction. A minimum-occupancy filter is
therefore applied and reported, so that the STD column cannot flatter
itself.

PRE-REGISTRATION (written before the run).

  (1) RULE for (Q1). If MEAN leaves |excess kurtosis| < 3 sigma it
      confirms item 1 as written and there is no finding.
  (2) RULE for (Q2). The scale-mask reading is declared only if the
      mixture prediction 3 Var(s^2)/(E s^2)^2 accounts for at least
      70% of the excess left by MEAN, at fixed occupancy floor.
  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      v1's own `lab_gaussian_half_audit.py` result file reports, at
      depth 8 with centring only, excess kurtosis +0.1704 at z = 98,
      and reports that per-cell STANDARDISATION removes 101.1% of it
      while the depth ladder removes -3.3%. Its own pre-registered
      decision rule then declares "(ii) A VARIANCE MASK", which its
      docstring calls "a NEW mask the conjecture does not name". So I
      predict MEAN leaves ~ +0.17, STD collapses it, the mixture
      formula accounts for most of it, and item 1's "removing cell
      means alone" describes the STD column and not the MEAN column.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ODD_PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)


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


def conv(a, b, X):
    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2
    return np.fft.irfft(np.fft.rfft(np.pad(a, (0, nf - len(a))))
                        * np.fft.rfft(np.pad(b, (0, nf - len(b)))), nf)[: X + 1]


def exkurt(x):
    c = x - x.mean()
    v = float((c * c).mean())
    return float((c ** 4).mean()) / v ** 2 - 3.0


def hn(x):
    c = x - x.mean()
    return float(np.abs(c).mean()) / math.sqrt(float((c * c).mean()))


def main():
    X = 16_000_000
    LO = 100_000
    DEPTH = 8
    MINOCC = 200

    mu, lam = sieve_mu_lambda(X)
    C = conv(mu, lam, X)
    V = conv((mu != 0).astype(np.float64), lam ** 2, X)
    Ns = np.arange(LO + LO % 2, X + 1, 2)
    Z = C[Ns] / np.sqrt(V[Ns])

    key = np.zeros(len(Ns), dtype=np.int64)
    for i, q in enumerate(ODD_PRIMES[:DEPTH]):
        key |= ((Ns % q) == 0).astype(np.int64) << i

    print("Re-verification of Conjecture 14, measurement 1")
    print(f"Z = C(N)/sqrt(V(N)); cells keyed by "
          f"{ODD_PRIMES[:DEPTH]}; min cell occupancy {MINOCC}")
    print()
    print("  RAW  = Z            MEAN = Z - cell mean  (item 1's stated")
    print("  protocol)          STD  = (Z - cell mean)/cell sd")
    print()
    hdr = (f"{'band':>22} {'n':>9} {'cells':>6} {'RAW kurt':>9} "
           f"{'MEAN kurt':>10} {'z':>7} {'STD kurt':>9} {'z':>7} "
           f"{'mixture pred':>13} {'accounts':>9}")
    print(hdr)
    print("-" * len(hdr))

    b = LO
    rows = []
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) < 20000:
            b = hi
            continue
        z = Z[sel]
        k = key[sel]
        uq, inv = np.unique(k, return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uq))
        tot = np.bincount(inv, weights=z, minlength=len(uq))
        m = tot / cnt
        cen = z - m[inv]
        tot2 = np.bincount(inv, weights=cen ** 2, minlength=len(uq))
        var = tot2 / np.maximum(cnt - 1, 1)
        keep = cnt[inv] >= MINOCC
        nk = int(keep.sum())
        zc = cen[keep]
        zs = cen[keep] / np.sqrt(var[inv][keep])

        big = cnt >= MINOCC
        w = cnt[big].astype(np.float64)
        s2 = var[big]
        Es2 = float((w * s2).sum() / w.sum())
        Vs2 = float((w * (s2 - Es2) ** 2).sum() / w.sum())
        mix = 3.0 * Vs2 / Es2 ** 2

        kr, km, ks = exkurt(z), exkurt(zc), exkurt(zs)
        se = math.sqrt(24.0 / nk)
        rows.append((b, hi, nk, km, ks, mix))
        print(f"{f'{b}-{hi}':>22} {nk:>9} {int(big.sum()):>6} "
              f"{kr:>9.4f} {km:>10.4f} {km/se:>7.1f} {ks:>9.4f} "
              f"{ks/se:>7.1f} {mix:>13.4f} "
              f"{mix/km if km else float('nan'):>8.1%}")
        b = hi

    print()
    print("(1) RULE was: item 1 stands as written if MEAN leaves")
    print("    |excess kurtosis| < 3 sigma. Read the 'MEAN kurt' z.")
    print("(2) RULE was: the scale-mask reading is declared if the")
    print("    mixture prediction accounts for >= 70% of the MEAN")
    print("    excess. Read the last column.")
    print()
    print("    v1's own lab_gaussian_half_audit.txt, depth 8, centring")
    print("    only: excess kurtosis +0.1704 at z = 98; per-cell")
    print("    standardisation removes 101.1%, the depth ladder -3.3%;")
    print("    its pre-registered rule then declares a VARIANCE MASK.")
    print("    The paper attributes +0.1704 to 'an S N-based scale' and")
    print("    attributes the near-zero value to 'removing cell means")
    print("    alone'.")
    print()

    # occupancy sensitivity: does STD collapse only because the sd is
    # estimated from few points?
    print("(3) occupancy sensitivity of the STD column (top octave)")
    sel = (Ns >= X // 2) & (Ns < X)
    z = Z[sel]
    k = key[sel]
    uq, inv = np.unique(k, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uq))
    tot = np.bincount(inv, weights=z, minlength=len(uq))
    cen = z - (tot / cnt)[inv]
    var = np.bincount(inv, weights=cen ** 2,
                      minlength=len(uq)) / np.maximum(cnt - 1, 1)
    print(f"    {'min occ':>9} {'n kept':>9} {'cells':>6} "
          f"{'MEAN kurt':>10} {'STD kurt':>9}")
    for mo in (50, 200, 1000, 5000):
        keep = cnt[inv] >= mo
        if keep.sum() < 10000:
            continue
        print(f"    {mo:>9} {int(keep.sum()):>9} "
              f"{int((cnt>=mo).sum()):>6} {exkurt(cen[keep]):>10.4f} "
              f"{exkurt(cen[keep]/np.sqrt(var[inv][keep])):>9.4f}")
    print()
    print("    a STD column that keeps collapsing as the occupancy floor")
    print("    RISES is removing real scale structure; one that only")
    print("    collapses at low occupancy is fitting noise.")
    print("DONE")


if __name__ == "__main__":
    main()
