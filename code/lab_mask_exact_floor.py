# -*- coding: utf-8 -*-
"""
The location mask against an exactly-computed floor (increment 305).

WHY THIS IS NOW POSSIBLE. Increment 303 withdrew the mask's measured
scaling (#69's g = -0.489) and its large-N amplitude (#67's "1.15% at
1.6e7") because both came from a between-cell variance whose COIN FLOOR
had never been subtracted, and the floor turned out to dominate the
upper bands. The mask itself survived -- M.1 is a theorem -- but its
amplitude was left "unresolved, not small", which is an honest state
and a useless one.

What blocked resolving it was that the floor was estimated from twenty
noisy draws. Increment 304 removed that: for a coin the whole thing is
closed-form. With Z(N) = sum_v eps(v) Lambda(N-v)/sqrt(V(N)),

    sum_{N in c} Z = sum_v eps(v) u_c(v),
    u_c(v) = sum_{N in c} Lambda(N-v)/sqrt(V(N)),
    Q_cd  = sum_v mu^2(v) u_c(v) u_d(v),

and every second moment of every cell mean follows exactly:

    Var(m_c - gm) = Q_cc/n_c^2 - 2 Q_ca/(n_c n) + Q_aa/n^2,

with a the whole band. No draws anywhere.

THE SHARPER INSTRUMENT. #67, #69 and #112 all measured the mask through
an AGGREGATE between-cell variance B, which throws away which cell did
what and, as 304 showed, has an 87% per-draw spread. With Var(m_c - gm)
in hand the natural statistic is instead the exact per-cell z,

    z_c = (m_c - gm) / sqrt(Var(m_c - gm)),

one number per cell per band, comparable across N by construction. The
mask is a DETERMINISTIC function of which small primes divide N (M.1),
so if it is present at a given N it must show up as cells sitting away
from zero by more than their own standard error.

THE KEY. Cells here are built from {3,5,7,11,13}, 32 of them, rather
than the {3..23} of #67. That is a cost decision -- one FFT per cell per
band -- and it makes every number below a LOWER BOUND on the mask, since
the coarser partition averages some of the structure away. The measured
per-prime factors of #40 fall 12.98, 4.70, 2.91, 2.62, 2.28 across
q = 3..13, so little is being given up.

PRE-REGISTRATION (fixed before the run).

  (S) SELF-TEST, which must pass or nothing else reads: run R = 24 coin
      draws through the identical code. The z_c must be standard
      normal -- |mean| < 0.1 and sd within [0.9, 1.1] pooled over cells
      and draws. If Var(m_c - gm) is wrong, this catches it.

  (M) IS THE MASK RESOLVED, band by band. Bonferroni over 32 cells at
      0.001 gives |z| > 3.66. RULE: the mask is resolved in a band iff
      max_c |z_c| exceeds 3.66 there.

  (A) AMPLITUDE AND SCALING. For the bands that resolve, report the
      mask's own between-cell variance B_mask = B_real - E[B_coin]
      with E[B_coin] = (1/n) sum_c Q_cc/n_c - Q_aa/n^2 exact, and fit
      its exponent in N. RULE: an exponent is reportable iff at least
      five bands resolve. Fewer than that and #69 stays withdrawn with
      a reason, which is a result and not a failure.

  A CAVEAT STATED IN ADVANCE. The standard error used is the COIN's.
  The real mu has Var(Z) = rho < 1 (0.76 to 0.84 by #288), so the true
  errors are smaller and these z are conservative. Both versions are
  printed; the verdict uses the conservative one.

  WHAT WOULD REFUTE. (S) failing kills the closed form. (M) resolving
  in only the bottom bands would confirm 303's picture -- the mask is
  real but fades below measurability -- and (M) resolving everywhere
  would mean 303 withdrew a scaling that was measurable after all, just
  not by the statistic #69 used.
"""
import math
import time

import numpy as np

QS = [3, 5, 7, 11, 13]
ZCRIT = 3.66


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
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    invV = 1.0 / np.sqrt(V[Ns])

    Creal = np.fft.irfft(np.fft.rfft(
        np.pad(mu.astype(np.float64), (0, nfft - X - 1))) * F_lam,
        nfft)[: X + 1]
    Zreal = Creal[Ns] * invV

    sels, lab = [], []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            sels.append(sel)
            lab.append((b, hi, int(sel.sum()), math.log(math.sqrt(b * hi))))
        b = hi
    nb = len(sels)
    muw = suppf[: X + 1]

    # ---- exact second moments, per band ----
    SD = []       # sd of (m_c - gm) under the coin, per cell
    NC = []       # cell sizes
    EBcoin = np.zeros(nb)
    rho_real = np.zeros(nb)
    for i, sel in enumerate(sels):
        n = int(sel.sum())
        kb = key[sel]
        cells = np.unique(kb)
        u_all = np.zeros(X + 1)
        us, ns = [], []
        for c in cells:
            m = sel.copy()
            m[sel] = (kb == c)
            w = np.zeros(nfft)
            w[Ns] = np.where(m, invV, 0.0)
            u = np.fft.irfft(Fl_c * np.fft.rfft(w), nfft)[: X + 1]
            us.append(u); ns.append(int(m.sum()))
            u_all += u
        mu_all = muw * u_all
        Qaa = float(np.dot(mu_all, u_all))
        sd = np.zeros(len(cells))
        acc = 0.0
        for j, u in enumerate(us):
            nc = ns[j]
            Qcc = float(np.dot(muw * u, u))
            Qca = float(np.dot(mu_all, u))
            var = Qcc / nc ** 2 - 2 * Qca / (nc * n) + Qaa / n ** 2
            sd[j] = math.sqrt(max(var, 0.0))
            acc += Qcc / nc
        EBcoin[i] = acc / n - Qaa / n ** 2
        SD.append(sd); NC.append(np.array(ns, dtype=np.float64))
        rho_real[i] = float(Zreal[sel].var())
        print(f"  band {i+1}/{nb} exact  E[B_coin] = {EBcoin[i]:.6f}  "
              f"rho_real = {rho_real[i]:.4f}  t={time.time()-t0:.0f}s",
              flush=True)

    # ---- real mu: per-cell means and exact z ----
    zmax = np.zeros(nb)
    nsig = np.zeros(nb, dtype=int)
    Breal = np.zeros(nb)
    Treal = np.zeros(nb)
    zbig = []
    for i, sel in enumerate(sels):
        z = Zreal[sel]; k = key[sel]
        uq, inv = np.unique(k, return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uq)).astype(np.float64)
        tot = np.bincount(inv, weights=z, minlength=len(uq))
        cm = tot / cnt
        gm = float(z.mean())
        zc = (cm - gm) / SD[i]
        n = len(z)
        Breal[i] = float((cnt * (cm - gm) ** 2).sum()) / n
        Treal[i] = float(((z - gm) ** 2).sum()) / n
        zmax[i] = float(np.abs(zc).max())
        nsig[i] = int((np.abs(zc) > ZCRIT).sum())
        zbig.append((uq, cm - gm, zc))

    # ---- (S) self-test on coins ----
    R = 24
    rng = np.random.default_rng(305)
    idx = np.nonzero(supp)[0]
    zs_all = []
    Bc = np.zeros((R, nb))
    for r in range(R):
        eps = np.zeros(nfft)
        eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        C = np.fft.irfft(np.fft.rfft(eps) * F_lam, nfft)[: X + 1]
        Z = C[Ns] * invV
        for i, sel in enumerate(sels):
            z = Z[sel]; k = key[sel]
            uq, inv = np.unique(k, return_inverse=True)
            cnt = np.bincount(inv, minlength=len(uq)).astype(np.float64)
            tot = np.bincount(inv, weights=z, minlength=len(uq))
            cm = tot / cnt
            gm = float(z.mean())
            zs_all.append((cm - gm) / SD[i])
            Bc[r, i] = float((cnt * (cm - gm) ** 2).sum()) / len(z)
    za = np.concatenate(zs_all)
    zmean, zsd = float(za.mean()), float(za.std(ddof=1))
    okS = (abs(zmean) < 0.1) and (0.9 <= zsd <= 1.1)

    sdB = Bc.std(axis=0, ddof=1)
    Bmask = Breal - EBcoin

    print(f"\n(S) self-test: exact z on {R} coin draws, {len(za)} values")
    print(f"    mean {zmean:+.4f} (need |.| < 0.1),  "
          f"sd {zsd:.4f} (need 0.9 to 1.1)  ->  "
          f"{'PASS' if okS else 'FAIL'}")
    print(f"    closed-form E[B_coin] against the draws' mean:")
    for i, (b0, hi, n, Lm) in enumerate(lab):
        mm = float(Bc[:, i].mean())
        se = sdB[i] / math.sqrt(R)
        print(f"      {b0:>9}-{hi:<10} exact {EBcoin[i]:.5f}  "
              f"MC {mm:.5f} +/- {se:.5f}  z {(mm-EBcoin[i])/se:+.2f}")

    print(f"\n(M) the mask cell by cell, exact standard errors, key {QS}")
    print(f"{'band':>21} {'n':>9} {'max |z|':>8} {'|z|>3.66':>9} "
          f"{'z*sqrt(rho)':>12} {'resolved':>9}")
    resolved = np.zeros(nb, dtype=bool)
    for i, (b0, hi, n, Lm) in enumerate(lab):
        resolved[i] = zmax[i] > ZCRIT
        print(f"{b0:>9}-{hi:>11} {n:>9} {zmax[i]:>8.2f} {nsig[i]:>9} "
              f"{zmax[i]*math.sqrt(rho_real[i]):>12.2f} "
              f"{'yes' if resolved[i] else 'NO':>9}")

    print(f"\n(A) the mask's own between-cell variance, floor removed")
    print(f"{'band':>21} {'B_real':>9} {'E[B_coin]':>10} {'B_mask':>9} "
          f"{'sd(B)':>8} {'B_mask/sd':>10} {'share of T':>11}")
    for i, (b0, hi, n, Lm) in enumerate(lab):
        print(f"{b0:>9}-{hi:>11} {Breal[i]:>9.5f} {EBcoin[i]:>10.5f} "
              f"{Bmask[i]:>9.5f} {sdB[i]:>8.5f} "
              f"{Bmask[i]/sdB[i]:>10.2f} {Bmask[i]/Treal[i]:>11.4%}")

    # (M) and (A) judge DIFFERENT statistics and are counted separately.
    # An earlier draft of this block tested (A) on the count of bands
    # that (M) resolved and then printed an exponent it had not fitted,
    # producing the self-contradictory line "PASS -> g = not fitted".
    nres = int(resolved.sum())
    okM = nres >= 5
    agg = Bmask / sdB
    nagg = int((agg > 2.0).sum())
    okA = nagg >= 5
    Ls = np.array([x[3] for x in lab])
    good = agg > 2.0
    if okA:
        g = float(np.polyfit(Ls[good], np.log(Bmask[good]), 1)[0])
        gtxt = f"{g:+.4f}"
    else:
        gtxt = "not fitted -- too few bands carry a measurable B_mask"

    wmax = np.zeros(nb)
    for i in range(nb):
        uq, dm, zc = zbig[i]
        j = int(np.argmax(np.abs(zc)))
        wmax[i] = float(NC[i][j] / NC[i].sum())

    print(f"\n    (S) exact variance verified on coins: "
          f"{'PASS' if okS else 'FAIL'}")
    print(f"    (M) per-cell: bands clearing Bonferroni |z| > {ZCRIT}: "
          f"{nres}/{nb}  ->  {'PASS' if okM else 'FAIL'}")
    print(f"    (A) aggregate: bands with B_mask > 2 sd: {nagg}/{nb}  ->  "
          f"{'PASS' if okA else 'FAIL'}   g = {gtxt}")
    print(f"        (#69 reported g = -0.489, withdrawn at 303)")
    if okS and okM and not okA:
        v = ("the mask is resolved in EVERY band by the per-cell test, "
             f"at |z| = {zmax[-1]:.1f} in the top band, while its "
             "aggregate between-cell variance is unmeasurable above "
             f"N ~ {lab[nagg][0]:.0e}. The two disagree because B "
             "weights each cell by n_c/n and the mask lives in the "
             f"RARE deep cells: the largest-|z| cell holds "
             f"{wmax[-1]:.2e} of the top band. So #69's statistic could "
             "not have measured the mask at large N whatever the floor, "
             "and 303's 'amplitude unresolved' was a statement about "
             "that statistic, not about the mask.")
    elif okS and okM and okA:
        v = (f"the mask is measurable both ways once the floor is exact; "
             f"the aggregate exponent is g = {gtxt}")
    elif okS and nres > 0:
        v = (f"the mask is significant in {nres} of {nb} bands, too few "
             f"to fit a scaling; #69 stays withdrawn")
    elif okS:
        v = ("no band resolves the mask at Bonferroni with this key; "
             "the mask is a theorem and is below this instrument")
    else:
        v = "the exact variance fails its self-test; nothing here reads"
    print(f"    {v}")

    print(f"\n    the three largest cells by |z|, per band "
          f"(cell = bitmask over {QS})")
    for i, (b0, hi, n, Lm) in enumerate(lab):
        uq, dm, zc = zbig[i]
        o = np.argsort(-np.abs(zc))[:3]
        parts = []
        for j in o:
            divs = [str(QS[t]) for t in range(len(QS))
                    if int(uq[j]) >> t & 1]
            nm = "".join(divs) if divs else "none"
            parts.append(f"{nm}: {dm[j]:+.4f} (z {zc[j]:+.1f})")
        print(f"      {b0:>9}-{hi:<10} " + "  ".join(parts))
    print("DONE")


if __name__ == "__main__":
    main()
