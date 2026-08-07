# -*- coding: utf-8 -*-
"""
The location mask's decay law, per depth, with exact error bars
(increment 305b).

WHAT 305a ESTABLISHED. With the coin floor in closed form and the
per-cell variance exact, the mask is resolved in EVERY band to
1.6e7 -- max |z| running 11.1, 10.9, 12.9, 12.6, 11.2, 10.3, 9.4, 8.4 --
while its aggregate between-cell variance B is unmeasurable above
N ~ 8e5. The two disagree because B weights each cell by n_c/n and the
mask lives in the RARE deep cells: the largest-|z| cell holds about
1/15015 of a band. So #69's statistic could not have measured the mask
at large N whatever the floor.

WHAT IS STILL OPEN. #69's exponent stays withdrawn and nothing has
replaced it. The per-cell means do fall -- the cell 3*5*7*11*13 reads
-7.09, -5.91, -6.40, -5.82, -4.98, -4.40, -3.89, -3.58 across the eight
octaves -- but eight points read off a table is not a law, and this
program has twice fitted a power of N where a power of log N was right
(#36, the missing sqrt(log N) in the variance law; #47, an exponent
inflated 45% by fitting through a pre-asymptotic turn).

So this run fits the law properly, and the thing that makes it possible
is that every point now carries an EXACT standard error:

    sd(m_c - gm)^2 = Q_cc/n_c^2 - 2 Q_ca/(n_c n) + Q_aa/n^2,
    Q_cd = sum_v mu^2(v) u_c(v) u_d(v),
    u_c(v) = sum_{N in c} Lambda(N-v)/sqrt(V(N)).

Half-octave bands give 15 points per depth instead of 8, and only the
six cells being tested need a transform, so the cost is 7 FFTs a band.

PRE-REGISTRATION (fixed before the run).

  Track six cells by DEPTH -- the exact divisor pattern, not "at least":
    depth 5  3,5,7,11,13 all divide N
    depth 4  3,5,7,11 divide, 13 does not
    depth 3  3,5,7 divide, 11 and 13 do not
    depth 2  3,5 divide, 7,11,13 do not
    depth 1  3 divides, 5,7,11,13 do not
    depth 0  none of them divides
  Points with n_c < 2 are dropped; every other point is weighted by
  1/sd^2 from the exact variance.

  THREE MODELS for |m| against N, fitted by weighted least squares on
  log|m|:
    K   constant           (no decay)
    P   A N^-a             (power of N)
    L   A (log N)^-b       (power of log N)

  RULES.
    (1) DECAY: the mask decays at depth d iff both P and L beat K by a
        factor 2 in weighted RSS. If they do not, the falling numbers
        in 305a's table are within their own errors and no decay is
        established.
    (2) WHICH LAW: P is preferred iff its weighted RSS is below L's by
        a factor 2, and L iff the reverse. Otherwise the two are NOT
        DISTINGUISHABLE over this range -- which, over a factor 160 in
        N, is the outcome I expect and the one this program has twice
        failed to report when it was true.
    (3) DEPTH: report a and b per depth. No rule is attached; the
        question of whether the exponent depends on depth has never
        been asked and this run is the first look.

  WHAT WOULD REFUTE. (1) failing at every depth would mean the mask's
  amplitude is flat in N on this range and the decay seen in 305a is an
  artefact of reading a table without errors.
"""
import math
import time

import numpy as np

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


def wfit(x, y, w):
    """Weighted least squares y = c0 + c1 x; returns coeffs and RSS."""
    A = np.stack([np.ones_like(x), x], axis=1)
    W = w[:, None]
    c = np.linalg.lstsq(A * np.sqrt(W), y * np.sqrt(w), rcond=None)[0]
    r = y - A @ c
    return c, float((w * r * r).sum())


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam = sieve(X)
    nfft = 1
    while nfft < 2 * (X + 1):
        nfft *= 2
    supp = (mu != 0)
    suppf = supp.astype(np.float64)
    F_supp = np.fft.rfft(np.pad(suppf, (0, nfft - X - 1)))
    F_lam = np.fft.rfft(np.pad(lam, (0, nfft - X - 1)))
    Fl_c = np.conj(F_lam)
    V = np.fft.irfft(F_supp * np.fft.rfft(
        np.pad(lam ** 2, (0, nfft - X - 1))), nfft)[: X + 1]
    print(f"sieve + V  t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(lo, X + 1, 2)
    invV = 1.0 / np.sqrt(V[Ns])
    Creal = np.fft.irfft(np.fft.rfft(
        np.pad(mu.astype(np.float64), (0, nfft - X - 1))) * F_lam,
        nfft)[: X + 1]
    Zreal = Creal[Ns] * invV
    muw = suppf[: X + 1]

    div = [(Ns % q) == 0 for q in QS]
    cellmask = {}
    for d in DEPTHS:
        m = np.ones(len(Ns), dtype=bool)
        for j in range(len(QS)):
            m &= div[j] if j < d else ~div[j]
        cellmask[d] = m

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
        w = np.zeros(nfft)
        w[Ns] = vals
        return np.fft.irfft(Fl_c * np.fft.rfft(w), nfft)[: X + 1]

    rows = {d: [] for d in DEPTHS}
    for bi, (b0, hi, sel) in enumerate(bands):
        n = int(sel.sum())
        u_all = ucorr(np.where(sel, invV, 0.0))
        mu_all = muw * u_all
        Qaa = float(np.dot(mu_all, u_all))
        gm = float(Zreal[sel].mean())
        Nmid = math.sqrt(b0 * hi)
        for d in DEPTHS:
            m = sel & cellmask[d]
            nc = int(m.sum())
            if nc < 2:
                continue
            u = ucorr(np.where(m, invV, 0.0))
            Qcc = float(np.dot(muw * u, u))
            Qca = float(np.dot(mu_all, u))
            var = Qcc / nc ** 2 - 2 * Qca / (nc * n) + Qaa / n ** 2
            if var <= 0:
                continue
            dm = float(Zreal[m].mean()) - gm
            rows[d].append((Nmid, nc, dm, math.sqrt(var)))
        print(f"  band {bi+1}/{len(bands)}  N~{Nmid:.3e}  "
              f"t={time.time()-t0:.0f}s", flush=True)

    print(f"\nthe mask by depth: cell mean minus band mean, exact se")
    for d in DEPTHS:
        r = rows[d]
        print(f"\n  depth {d}  ({'*'.join(str(q) for q in QS[:d]) or 'none'}"
              f" divide N, the rest do not)   {len(r)} points")
        print(f"{'N':>12} {'n_c':>8} {'m_c - gm':>10} {'se':>8} {'z':>8}")
        for Nmid, nc, dm, se in r:
            print(f"{Nmid:>12.3e} {nc:>8} {dm:>10.4f} {se:>8.4f} "
                  f"{dm/se:>8.2f}")

    print(f"\nthe decay law, weighted by the exact errors")
    print(f"{'depth':>6} {'pts':>4} {'RSS const':>10} {'RSS N^-a':>10} "
          f"{'a':>8} {'RSS log^-b':>11} {'b':>8} {'verdict':>22}")
    verdicts = {}
    for d in DEPTHS:
        r = [x for x in rows[d] if abs(x[2]) > 0]
        if len(r) < 4:
            print(f"{d:>6} {len(r):>4}   too few points to fit")
            verdicts[d] = "too few points"
            continue
        Nm = np.array([x[0] for x in r])
        dm = np.array([x[2] for x in r])
        se = np.array([x[3] for x in r])
        y = np.log(np.abs(dm))
        # delta method: sd(log|m|) = se/|m|
        w = (np.abs(dm) / se) ** 2
        cK = float((w * y).sum() / w.sum())
        rK = float((w * (y - cK) ** 2).sum())
        cP, rP = wfit(np.log(Nm), y, w)
        cL, rL = wfit(np.log(np.log(Nm)), y, w)
        decay = (rK / max(rP, 1e-300) > 2.0) and (rK / max(rL, 1e-300) > 2.0)
        if not decay:
            vd = "no decay established"
        elif rP * 2.0 < rL:
            vd = "power of N"
        elif rL * 2.0 < rP:
            vd = "power of log N"
        else:
            vd = "not distinguishable"
        verdicts[d] = vd
        print(f"{d:>6} {len(r):>4} {rK:>10.2f} {rP:>10.2f} "
              f"{-cP[1]:>8.4f} {rL:>11.2f} {-cL[1]:>8.4f} {vd:>22}")

    ndec = sum(1 for d in DEPTHS if verdicts.get(d, "") not in
               ("no decay established", "too few points"))
    npow = sum(1 for d in DEPTHS if verdicts.get(d) == "power of N")
    nlog = sum(1 for d in DEPTHS if verdicts.get(d) == "power of log N")
    nund = sum(1 for d in DEPTHS if verdicts.get(d) == "not distinguishable")
    print(f"\n    (1) depths where a decay beats a constant by 2x: "
          f"{ndec}/{len(DEPTHS)}")
    print(f"    (2) of those: power of N {npow}, power of log N {nlog}, "
          f"not distinguishable {nund}")
    if ndec == 0:
        v = ("the mask's amplitude is flat in N within its own errors on "
             "this range; 305a's falling column is not a decay")
    elif nund >= max(npow, nlog):
        v = ("the mask decays, and over a factor 160 in N the data do "
             "NOT separate a power of N from a power of log N. Neither "
             "may be quoted as the law; #69's exponent stays withdrawn "
             "and is not replaced")
    elif npow > nlog:
        v = ("the mask decays as a power of N at the majority of depths")
    else:
        v = ("the mask decays as a power of log N at the majority of "
             "depths -- the same species of correction as #36 and #47")
    print(f"    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
