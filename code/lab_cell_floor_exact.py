# -*- coding: utf-8 -*-
"""
The cell floor in closed form, with no Monte Carlo (increment 304d).

WHY THIS RUN EXISTS. Increments 304b and 304c judged the mechanism for
#113 against a MEASURED coin B/T from eight draws. Comparing the two
runs afterwards shows that statistic is not stable: band 6.4e6-1.28e7
gave 0.0324 in one and 0.0781 in the other, a factor 2.4 apart. So
304c's "RULE B passes 7 of 8 bands, within 1.5x" was evaluated against
a target whose own run-to-run spread I had never measured, and it is
not trustworthy at that tolerance. That is a fault in my test design,
recorded as its own correction.

The fix is not more draws. For a coin the whole quantity is available
exactly. Write Z(N) = sum_v eps(v) Lambda(N-v)/sqrt(V(N)). Then for any
set c of N,

    sum_{N in c} Z(N) = sum_v eps(v) u_c(v),
    u_c(v) = sum_{N in c} Lambda(N-v)/sqrt(V(N)),

and since eps is independent +/-1 on {mu != 0},

    E[(sum_{N in c} Z)^2] = sum_v mu^2(v) u_c(v)^2

with no approximation whatever. Every term of

    E[B] = sum_c (n_c/n) E[(m_c - gm)^2],   E[T] = 1 - E[gm^2]

then follows from the u_c, which are cross-correlations of Lambda with
the 1/sqrt(V)-weighted indicator of each cell -- one FFT per cell.
E[Z^2] = V/V = 1 exactly, which is what makes E[T] that simple.

To keep the number of cells to eight the key here is built from
{3, 5, 7} rather than {3,...,23}. That coarsens the cells and lowers
the floor; it does not change what is being asked, and the singular-
series prediction is recomputed for the same coarse key.

PRE-REGISTRATION (fixed before the run).

  (S) SELF-TEST, which must pass or nothing else is readable: the mean
      B/T over R = 12 actual coin draws must sit within 3 standard
      errors of the closed form. If the algebra above is wrong this
      catches it, and the same run reports the per-draw spread that
      304c should have had.

  (P) PREDICTION: closed-form B/T against a(sep) * (E_same[S_2] -
      E_all[S_2]), with a taken at the band's own mean pair separation
      exactly as in 304c. RULE: within 1.5x in a majority of bands.
      Now the target carries no sampling error at all, so the tolerance
      means what it says.

  WHAT WOULD REFUTE. (S) failing kills the closed form. (P) failing
  leaves the singular series of the shift as the floor's shape and
  cause -- both already established at 0.9997 correlation and by the
  placebo key -- but not its size.
"""
import math
import time

import numpy as np

QS = [3, 5, 7]
C2 = 0.66016181584686957392
NPER = 10


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
    return mu, lam, spf


def sing(h, spf):
    s = 2.0 * C2
    m = int(h)
    while m % 2 == 0:
        m //= 2
    while m > 1:
        p = int(spf[m])
        s *= (p - 1.0) / (p - 2.0)
        while m % p == 0:
            m //= p
    return s


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam, spf = sieve(X)
    nfft = 1
    while nfft < 2 * (X + 1):
        nfft *= 2
    supp = (mu != 0)
    suppf = supp.astype(np.float64)
    F_supp = np.fft.rfft(np.pad(suppf, (0, nfft - X - 1)))
    lamp = np.pad(lam, (0, nfft - X - 1))
    F_lam = np.fft.rfft(lamp)
    V = np.fft.irfft(F_supp * np.fft.rfft(
        np.pad(lam ** 2, (0, nfft - X - 1))), nfft)[: X + 1]
    print(f"sieve + V  t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    invV = 1.0 / np.sqrt(V[Ns])

    sels, lab = [], []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            sels.append(sel)
            lab.append((b, hi, int(sel.sum())))
        b = hi
    nb = len(sels)
    muw = suppf[: X + 1]
    Fl_c = np.conj(F_lam)

    def ucorr(w_on_Ns):
        """u(v) = sum_N w(N) Lambda(N-v), for v in [0, X]."""
        w = np.zeros(nfft)
        w[Ns] = w_on_Ns
        return np.fft.irfft(Fl_c * np.fft.rfft(w), nfft)[: X + 1]

    exact = np.zeros(nb)
    for i, sel in enumerate(sels):
        n = int(sel.sum())
        kb = key[sel]
        cells = np.unique(kb)
        u_all = np.zeros(X + 1)
        parts = []
        for c in cells:
            m = sel.copy()
            m[sel] = (kb == c)
            nc = int(m.sum())
            u = ucorr(np.where(m, invV, 0.0))
            parts.append((nc, u))
            u_all += u
        q_all = float(np.dot(muw * u_all, u_all))
        Egm2 = q_all / n ** 2
        EB = 0.0
        for nc, u in parts:
            q_cc = float(np.dot(muw * u, u))
            q_ca = float(np.dot(muw * u, u_all))
            EB += (nc / n) * (q_cc / nc ** 2 - 2 * q_ca / (nc * n) + Egm2)
        exact[i] = EB / (1.0 - Egm2)
        print(f"  band {i+1}/{nb} exact B/T = {exact[i]:.6f}  "
              f"t={time.time()-t0:.0f}s", flush=True)

    # ---- (S) self-test against real coin draws ----
    R = 12
    rng = np.random.default_rng(3042)
    idx = np.nonzero(supp)[0]
    mc = np.zeros((R, nb))
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
            gm = float(z.mean())
            B = float((cnt * (tot / cnt - gm) ** 2).sum()) / len(z)
            T = float(((z - gm) ** 2).sum()) / len(z)
            mc[r, i] = B / T
    mmean = mc.mean(axis=0)
    msd = mc.std(axis=0, ddof=1)
    mse = msd / math.sqrt(R)

    # ---- a(h) by scale, and the singular-series excess ----
    scales = [2 * 10 ** e for e in range(1, 7)]
    a_sc = np.full((len(scales), nb), np.nan)
    pad = np.zeros(nfft)
    for si, sc in enumerate(scales):
        hs = [sc + 2 * t for t in range(NPER)]
        S = np.array([sing(h, spf) for h in hs])
        rr = np.full((NPER, nb), np.nan)
        for j, h in enumerate(hs):
            g = pad.copy()
            g[: X + 1 - h] = lam[: X + 1 - h] * lam[h:]
            num = np.fft.irfft(F_supp * np.fft.rfft(g), nfft)[: X + 1]
            keep = Ns + h <= X
            Nk = Ns[keep]
            rv = num[Nk] / np.sqrt(V[Nk] * V[Nk + h])
            for i, sel in enumerate(sels):
                s = sel[keep]
                if int(s.sum()) > 1000:
                    rr[j, i] = float(rv[s].mean())
        A = np.stack([S, np.ones_like(S)], axis=1)
        for i in range(nb):
            if not np.isnan(rr[:, i]).any():
                a_sc[si, i] = float(
                    np.linalg.lstsq(A, rr[:, i], rcond=None)[0][0])
        print(f"  scale h~{sc}  t={time.time()-t0:.0f}s", flush=True)

    exc = np.zeros(nb)
    for i, sel in enumerate(sels):
        Nb = Ns[sel]; kb = key[sel]
        i1 = rng.integers(0, len(Nb), 80_000)
        i2 = rng.integers(0, len(Nb), 80_000)
        ok = i1 != i2
        ha = np.abs(Nb[i1[ok]] - Nb[i2[ok]])
        order = np.argsort(kb, kind="stable")
        ks = kb[order]
        bnd = np.append(np.searchsorted(ks, np.unique(ks), side="left"),
                        len(ks))
        j1, j2 = [], []
        for c in range(len(bnd) - 1):
            a0, b0 = bnd[c], bnd[c + 1]
            if b0 - a0 < 2:
                continue
            t = min((b0 - a0) * 4, 40_000)
            u = order[rng.integers(a0, b0, t)]
            w = order[rng.integers(a0, b0, t)]
            g = u != w
            j1.append(u[g]); j2.append(w[g])
        hsm = np.abs(Nb[np.concatenate(j1)] - Nb[np.concatenate(j2)])
        sa = float(np.mean([sing(int(h), spf) for h in
                            ha[rng.integers(0, len(ha), 50_000)]]))
        ss = float(np.mean([sing(int(h), spf) for h in
                            hsm[rng.integers(0, len(hsm), 50_000)]]))
        exc[i] = ss - sa

    # ---- report ----
    print(f"\n(S) closed form against {R} coin draws, key from {QS}")
    print(f"{'band':>21} {'exact':>9} {'MC mean':>9} {'MC sd':>9} "
          f"{'se':>9} {'z':>7}")
    zs = np.zeros(nb)
    for i, (b0, hi, n) in enumerate(lab):
        zs[i] = (mmean[i] - exact[i]) / mse[i]
        print(f"{b0:>9}-{hi:>11} {exact[i]:>9.5f} {mmean[i]:>9.5f} "
              f"{msd[i]:>9.5f} {mse[i]:>9.5f} {zs[i]:>7.2f}")
    okS = bool((np.abs(zs) < 3).all())

    print(f"\n(P) closed form against the singular-series prediction")
    print(f"{'band':>21} {'sep':>9} {'a used':>9} {'E_s-E_a':>9} "
          f"{'pred':>9} {'exact':>9} {'ratio':>7}")
    rat = np.zeros(nb)
    for i, (b0, hi, n) in enumerate(lab):
        sep = b0 / 3.0
        si = int(np.argmin([abs(math.log(s / sep)) for s in scales]))
        col = a_sc[:, i]
        while np.isnan(col[si]) and si > 0:
            si -= 1
        pr = col[si] * exc[i]
        rat[i] = pr / exact[i]
        print(f"{b0:>9}-{hi:>11} {sep:>9.0f} {col[si]:>9.5f} "
              f"{exc[i]:>9.5f} {pr:>9.5f} {exact[i]:>9.5f} {rat[i]:>7.2f}")
    nP = int(((rat > 1 / 1.5) & (rat < 1.5)).sum())
    okP = nP > nb / 2

    print(f"\n    (S) every band within 3 se of the closed form: "
          f"{'PASS' if okS else 'FAIL'}  (max |z| = {np.abs(zs).max():.2f})")
    print(f"    per-draw relative spread of B/T: "
          f"{float(np.mean(msd/mmean)):.1%} on average -- this is what "
          f"304c's ratio test was resting on")
    print(f"    (P) prediction within 1.5x of the closed form in a "
          f"majority: {'PASS' if okP else 'FAIL'}  ({nP}/{nb} bands)")
    if okS and okP:
        v = ("#113 is answered. The cell floor is the singular series of "
             "the shift: same-cell pairs have h divisible by small "
             "primes more often, S_2(h) is larger there, and the excess "
             "is a property of PAIRS, so no sample size removes it. "
             "Shape (corr 0.9997), cause (placebo key) and size (this "
             "run) all hold.")
    elif okS:
        v = ("the closed form is verified but the singular-series "
             "prediction does not reproduce its size; shape and cause "
             "stand, size stays open")
    else:
        v = "the closed form fails its own self-test; nothing here reads"
    print(f"    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
