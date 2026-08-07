# -*- coding: utf-8 -*-
"""
Are the six mask exponents PREDICTED by A2's closed form? (inc. 335)

WHAT INCREMENT 334 LEFT. Supplying the errors A3 never had confirmed
that the six exponents differ (chi2/dof = 251) and, unasked, produced a
structure: excluding depth 1, which is not measurable, the exponent
rises MONOTONICALLY as the cell gets shallower,

    depth 5  0.1434 +/- 0.0155      depth 2  0.3686 +/- 0.0052
    depth 4  0.2152 +/- 0.0065      depth 0  0.6289 +/- 0.0121
    depth 3  0.2713 +/- 0.0040

with the steps at 5 to 30 standard errors. That ordering was measured,
not explained.

WHY IT SHOULD BE PREDICTABLE NOW. Increment 333 closed A2: the cell
floor's SIZE is the excess of the singular series over same-cell pairs,

    dm_c ~ a * ( E_same,c[S_2] - E_all[S_2] )   =:  a * D_c ,

to within 6% once the cells are weighted by n_c/n. If that is the
mechanism, it must also govern the floor's DECAY, and D_c's decay is
computable without running the wall at all.

AND THE UNKNOWN COEFFICIENT CANCELS. a is fitted per band and is common
to all six depths, so it shifts every predicted exponent by the same
amount. THE DIFFERENCES BETWEEN DEPTHS ARE THEREFORE PREDICTED WITH NO
FREE PARAMETER AT ALL. This run compares centred exponents, which is
the comparison that owes nothing to a.

D_c IS COMPUTED EXACTLY HERE, not sampled. Increment 312 drew pairs;
that sampling is what turned out to be the fault A2 died of. The mean
of S_2 over same-cell pairs is an autocorrelation,

    sum_{N,N' in c, N != N'} S_2(|N-N'|)  =  2 sum_{k>0} A_c(k) S_2(2k),

with A_c the autocorrelation of the cell's indicator inside the band,
so one FFT per cell per band gives it with no Monte-Carlo error. That
is hazard 8 answered for the prediction as well as for the measurement.

PRE-REGISTRATION (fixed before the run).

  (P1) BOTH SETS OF EXPONENTS, CENTRED. Report the measured exponent
       with its error (recomputed here, and it must reproduce 334) and
       the predicted exponent from D_c, each minus its
       inverse-variance-weighted mean. RULE: the measured column must
       reproduce 334's six values to 0.01, or this is not the same fit
       and nothing below reads.

  (P2) DOES THE MECHANISM ORDER THE DEPTHS? Spearman rank correlation
       between the centred measured and centred predicted exponents,
       over the five measurable depths -- depth 1 is excluded, because
       334 established its exponent is not measurable (a/SE = 0.8), and
       that exclusion is fixed here before the run.
       RULE: the mechanism orders the depths iff the correlation is
       at least +0.9.

  (P3) DOES IT EXPLAIN THE SIZE OF THE VARIATION? Compare the rms of
       the centred measured exponents against the rms of the residual
       (centred measured minus centred predicted).
       RULE: the mechanism explains the variation iff the residual rms
       is at most one third of the raw rms.

  WHAT WOULD REFUTE. (P2) failing means the singular-series excess does
  not even order the depths: A2's closed form would explain the floor's
  size within a band and say nothing about its decay, and A3's law
  would have to come from somewhere else. (P2) passing with (P3)
  failing is the interesting middle -- the mechanism is right in
  direction and wrong in magnitude, which is where A2 itself sat for
  twenty-nine increments.
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
C2 = 0.66016181584686957392


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
    return mu, lam, primes


def singular_series(X, primes):
    """S_2(h) for every h <= X, vectorised over the odd primes."""
    s = np.full(X + 1, 2.0 * C2)
    for p in primes:
        q = int(p)
        if q == 2:
            continue
        s[q::q] *= (q - 1.0) / (q - 2.0)
    s[0] = 0.0
    return s


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
    return c1, math.sqrt(max(s2 * S / D, 0.0))


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    return float((ra * rb).sum()
                 / math.sqrt((ra * ra).sum() * (rb * rb).sum()))


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    S2 = singular_series(X, primes)
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
    print(f"sieve + S_2 + V + C  t={time.time()-t0:.0f}s", flush=True)

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

    def esame(ind, off, step_pts):
        """Exact mean of S_2 over distinct pairs inside a point set.

        ind        boolean over the band's even N, in position order
        off        the N value of position 0
        step_pts   spacing of consecutive positions in N (here 2)
        """
        m = len(ind)
        nfb = 1
        while nfb < 2 * m:
            nfb *= 2
        x = ind.astype(np.float64)
        F = np.fft.rfft(x, nfb)
        A = np.fft.irfft(F * np.conj(F), nfb)[:m]
        k = np.arange(1, m)
        h = k * step_pts
        num = 2.0 * float(np.dot(A[1:], S2[h]))
        den = 2.0 * float(A[1:].sum())
        if den <= 0:
            return None
        return num / den

    meas = {d: [] for d in DEPTHS}
    pred = {d: [] for d in DEPTHS}
    for bi, (b0, hi, sel) in enumerate(bands):
        n = int(sel.sum())
        Nmid = math.sqrt(b0 * hi)
        u_all = ucorr(np.where(sel, invV, 0.0))
        mu_all = muw * u_all
        Qaa = float(np.dot(mu_all, u_all))
        gm = float(Z[sel].mean())

        # exact singular-series mean over ALL pairs in the band
        Nb = Ns[sel]
        off = int(Nb[0])
        pos = (Nb - off) // 2
        L = int(pos[-1]) + 1
        base = np.zeros(L, dtype=bool)
        base[pos] = True
        e_all = esame(base, off, 2)

        for d in DEPTHS:
            m = sel & cell[d]
            nc = int(m.sum())
            if nc < 2 or e_all is None:
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
            meas[d].append((Nmid, dm, math.sqrt(var)))

            ic = np.zeros(L, dtype=bool)
            ic[(Ns[m] - off) // 2] = True
            e_c = esame(ic, off, 2)
            if e_c is None:
                continue
            D = e_c - e_all
            if abs(D) < 1e-15:
                continue
            pred[d].append((Nmid, D))
        if (bi + 1) % 5 == 0:
            print(f"  band {bi+1}/{len(bands)}  "
                  f"t={time.time()-t0:.0f}s", flush=True)

    quoted = {5: 0.1434, 4: 0.2152, 3: 0.2713, 2: 0.3686,
              1: 0.0437, 0: 0.6289}
    print(f"\n(P1) measured against predicted, both raw and centred")
    print(f"{'depth':>6} {'a_meas':>9} {'SE':>8} {'334':>8} "
          f"{'a_pred':>9} {'SE':>8} {'D_c sign':>9}")
    am, sm, ap, sp, ds, bad = [], [], [], [], [], 0.0
    for d in DEPTHS:
        rm, rp = meas[d], pred[d]
        if len(rm) < 4 or len(rp) < 4:
            print(f"{d:>6}   too few points")
            continue
        Nm = np.array([x[0] for x in rm])
        dm = np.array([x[1] for x in rm])
        se = np.array([x[2] for x in rm])
        a, sa = wfit_se(np.log(Nm), np.log(np.abs(dm)),
                        (np.abs(dm) / se) ** 2)
        Np = np.array([x[0] for x in rp])
        Dv = np.array([x[1] for x in rp])
        p, spv = wfit_se(np.log(Np), np.log(np.abs(Dv)),
                         np.ones(len(Np)))
        bad = max(bad, abs(-a - quoted[d]))
        am.append(-a); sm.append(sa)
        ap.append(-p); sp.append(spv)
        ds.append(d)
        sgn = "+" if Dv.mean() > 0 else "-"
        print(f"{d:>6} {-a:>9.4f} {sa:>8.4f} {quoted[d]:>8.4f} "
              f"{-p:>9.4f} {spv:>8.4f} {sgn:>9}")

    am = np.array(am); sm = np.array(sm)
    ap = np.array(ap); sp = np.array(sp)
    ds = np.array(ds)
    okP1 = bad <= 0.01
    print(f"\n    (P1) the measured column reproduces 334 to "
          f"{bad:.4f}: {'PASS' if okP1 else 'FAIL'}")

    keep = ds != 1
    wm = 1.0 / sm[keep] ** 2
    cm = am[keep] - (am[keep] * wm).sum() / wm.sum()
    cp = ap[keep] - ap[keep].mean()
    print(f"\n    centred, over the five measurable depths "
          f"(depth 1 excluded before the run)")
    print(f"{'depth':>6} {'measured':>10} {'predicted':>11} "
          f"{'residual':>10}")
    for i, d in enumerate(ds[keep]):
        print(f"{d:>6} {cm[i]:>10.4f} {cp[i]:>11.4f} "
              f"{cm[i]-cp[i]:>10.4f}")

    rho = spearman(cm, cp)
    raw = float(np.sqrt((cm ** 2).mean()))
    res = float(np.sqrt(((cm - cp) ** 2).mean()))
    okP2 = rho >= 0.9
    okP3 = res <= raw / 3.0
    print(f"\n    (P2) the mechanism orders the depths "
          f"(Spearman >= +0.9): {'PASS' if okP2 else 'FAIL'}  "
          f"(rho = {rho:+.3f})")
    print(f"    (P3) it explains the size (residual rms <= raw/3): "
          f"{'PASS' if okP3 else 'FAIL'}  "
          f"(raw {raw:.4f}, residual {res:.4f}, "
          f"ratio {res/max(raw,1e-12):.2f})")

    if not okP1:
        v = ("the measured exponents do not reproduce increment 334, "
             "so this is not the same fit and no comparison reads")
    elif okP2 and okP3:
        v = (f"A3's law is A2's mechanism. The singular-series excess "
             f"over same-cell pairs both orders the six depths "
             f"(Spearman {rho:+.2f}) and accounts for "
             f"{100*(1-res/raw):.0f}% of the spread in the exponents, "
             f"with NO free parameter -- the fitted a is common to all "
             f"depths and cancels from the centred comparison. The "
             f"mask decays faster in shallow cells because the "
             f"same-cell singular-series excess decays faster there")
    elif okP2:
        v = (f"the mechanism has the right direction and the wrong "
             f"magnitude: it orders the depths (Spearman {rho:+.2f}) "
             f"but leaves {100*res/raw:.0f}% of the exponent spread as "
             f"residual. A2's closed form governs the floor's size "
             f"within a band and only partly its decay")
    else:
        v = (f"the singular-series excess does not order the depths "
             f"(Spearman {rho:+.2f}). A2's closed form explains the "
             f"floor's size within a band and says nothing about its "
             f"decay, so A3's law comes from somewhere else and this "
             f"natural extension is refuted")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
