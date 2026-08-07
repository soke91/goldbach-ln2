# -*- coding: utf-8 -*-
"""
Does the mask's exponent really vary with depth? (increment 334)

A3 SAYS THERE IS NO SINGLE MASK EXPONENT. Increment 305b fitted the
decay per depth and reported

    depth 5  a = 0.1434      depth 2  a = 0.3686
    depth 4  a = 0.2152      depth 1  a = 0.0437
    depth 3  a = 0.2713      depth 0  a = 0.6289

and concluded, from the roughly fourfold spread, that "there is no
single mask exponent" and that #69's one g was fitting a mixture.

NOT ONE OF THOSE SIX NUMBERS CARRIES A STANDARD ERROR. That is hazard 8
-- a spread quoted without the uncertainty of the things spread -- in
the same document that named it. Six exponents fitted to fifteen points
each, several of those points resting on cells with two or three
members, can easily scatter by a factor four without any of them
differing.

This run supplies the missing errors and asks whether a common exponent
is rejected.

PRE-REGISTRATION (fixed before the run).

  (E1) THE STANDARD ERROR ON EACH EXPONENT, from the weighted fit's own
       covariance, using the exact per-cell errors of increment 305 as
       weights. Reported per depth beside the exponent. RULE: none --
       this is the number that was missing.

  (E2) IS A COMMON EXPONENT REJECTED? Form the chi-square of the six
       exponents about their inverse-variance-weighted mean.
       RULE: a common exponent is rejected iff chi^2/dof > 2. If it is
       not rejected, A3's "no single exponent" overstates its evidence
       and the mixture reading of #69 needs restating.

  (E3) THE COMMON EXPONENT, if one is admissible: the
       inverse-variance-weighted mean and its own error. Reported
       whether or not (E2) rejects, since a rejected common value is
       still the right summary of what the data say about the average.

  WHAT WOULD REFUTE. (E2) rejecting confirms A3 as written, and the
  spread is real. That is a perfectly good outcome and the one A3
  currently asserts without evidence.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

QS = [3, 5, 7, 11, 13]
DEPTHS = [5, 4, 3, 2, 1, 0]


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]
            sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8)
    mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i])
        j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p)
        lg = math.log(int(p))
        while q <= X:
            lam[q] = lg
            q *= int(p)
    return mu, lam


def wfit_se(x, y, w):
    """Weighted least squares y = c0 + c1 x, with the slope's own SE."""
    S = w.sum()
    Sx = (w * x).sum()
    Sxx = (w * x * x).sum()
    Sy = (w * y).sum()
    Sxy = (w * x * y).sum()
    D = S * Sxx - Sx * Sx
    c1 = (S * Sxy - Sx * Sy) / D
    c0 = (Sxx * Sy - Sx * Sxy) / D
    r = y - (c0 + c1 * x)
    dof = max(len(x) - 2, 1)
    s2 = (w * r * r).sum() / dof
    var_c1 = s2 * S / D
    return c1, math.sqrt(max(var_c1, 0.0))


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam = sieve(X)
    nf = 1
    while nf < 2 * (X + 1):
        nf *= 2
    supp = (mu != 0).astype(np.float64)
    F_supp = np.fft.rfft(np.pad(supp, (0, nf - X - 1)))
    F_lam = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))
    Fl_c = np.conj(F_lam)
    V = np.fft.irfft(F_supp * np.fft.rfft(
        np.pad(lam ** 2, (0, nf - X - 1))), nf)[: X + 1]
    Ns = np.arange(lo, X + 1, 2)
    invV = 1.0 / np.sqrt(V[Ns])
    Creal = np.fft.irfft(np.fft.rfft(
        np.pad(mu.astype(np.float64), (0, nf - X - 1))) * F_lam,
        nf)[: X + 1]
    Z = Creal[Ns] * invV
    muw = supp[: X + 1]
    print(f"sieve + V + C  t={time.time()-t0:.0f}s", flush=True)

    div = [(Ns % q) == 0 for q in QS]
    cell = {}
    for d in DEPTHS:
        m = np.ones(len(Ns), dtype=bool)
        for j in range(len(QS)):
            m &= div[j] if j < d else ~div[j]
        cell[d] = m

    bands = []
    b = float(lo)
    while b < X:
        hi = min(b * math.sqrt(2.0), X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 500:
            bands.append((b, hi, sel))
        b = hi
    print(f"{len(bands)} half-octave bands", flush=True)

    def ucorr(vals):
        w = np.zeros(nf)
        w[Ns] = vals
        return np.fft.irfft(Fl_c * np.fft.rfft(w), nf)[: X + 1]

    rows = {d: [] for d in DEPTHS}
    for bi, (b0, hi, sel) in enumerate(bands):
        n = int(sel.sum())
        u_all = ucorr(np.where(sel, invV, 0.0))
        mu_all = muw * u_all
        Qaa = float(np.dot(mu_all, u_all))
        gm = float(Z[sel].mean())
        Nmid = math.sqrt(b0 * hi)
        for d in DEPTHS:
            m = sel & cell[d]
            nc = int(m.sum())
            if nc < 2:
                continue
            u = ucorr(np.where(m, invV, 0.0))
            var = (float(np.dot(muw * u, u)) / nc ** 2
                   - 2 * float(np.dot(mu_all, u)) / (nc * n)
                   + Qaa / n ** 2)
            if var <= 0:
                continue
            dm = float(Z[m].mean()) - gm
            if abs(dm) < 1e-12:
                continue
            rows[d].append((Nmid, dm, math.sqrt(var)))
        if (bi + 1) % 5 == 0:
            print(f"  band {bi+1}/{len(bands)}  "
                  f"t={time.time()-t0:.0f}s", flush=True)

    print(f"\n(E1) the exponent per depth, with the error it never had")
    print(f"{'depth':>6} {'pts':>4} {'a':>9} {'SE(a)':>9} "
          f"{'a/SE':>7} {'305b quoted':>12}")
    quoted = {5: 0.1434, 4: 0.2152, 3: 0.2713, 2: 0.3686,
              1: 0.0437, 0: 0.6289}
    est, ses, ds = [], [], []
    for d in DEPTHS:
        r = rows[d]
        if len(r) < 4:
            print(f"{d:>6} {len(r):>4}   too few points")
            continue
        Nm = np.array([x[0] for x in r])
        dm = np.array([x[1] for x in r])
        se = np.array([x[2] for x in r])
        y = np.log(np.abs(dm))
        w = (np.abs(dm) / se) ** 2
        a, sa = wfit_se(np.log(Nm), y, w)
        est.append(-a)
        ses.append(sa)
        ds.append(d)
        print(f"{d:>6} {len(r):>4} {-a:>9.4f} {sa:>9.4f} "
              f"{abs(a)/max(sa,1e-12):>7.1f} {quoted[d]:>12.4f}")

    est = np.array(est); ses = np.array(ses)
    wgt = 1.0 / ses ** 2
    amean = float((est * wgt).sum() / wgt.sum())
    amean_se = float(1.0 / math.sqrt(wgt.sum()))
    chi2 = float((wgt * (est - amean) ** 2).sum())
    dof = max(len(est) - 1, 1)
    okE2 = (chi2 / dof) > 2.0
    print(f"\n    (E2) a common exponent is REJECTED (chi2/dof > 2): "
          f"{'YES' if okE2 else 'NO'}  "
          f"(chi2 = {chi2:.2f} on {dof} dof, chi2/dof = "
          f"{chi2/dof:.2f})")
    print(f"    (E3) inverse-variance-weighted common exponent: "
          f"a = {amean:.4f} +/- {amean_se:.4f}")
    print(f"         spread of the six point estimates: "
          f"{est.max()/max(est.min(),1e-9):.1f}x, which is what A3 "
          f"quotes")

    if okE2:
        v = (f"the spread is real: chi2/dof = {chi2/dof:.1f} rejects a "
             f"common exponent, so A3 stands as written and #69's "
             f"single g was indeed fitting a mixture")
    else:
        v = (f"A3 overstates its evidence. The six exponents scatter by "
             f"{est.max()/max(est.min(),1e-9):.1f}x but carry errors "
             f"large enough that a common value is NOT rejected "
             f"(chi2/dof = {chi2/dof:.2f}), and that common value is "
             f"a = {amean:.4f} +/- {amean_se:.4f}. 'There is no single "
             f"mask exponent' was a spread quoted without the "
             f"uncertainty of the things spread -- hazard 8, in the "
             f"document that named it")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
