# -*- coding: utf-8 -*-
"""
Does the shift's correlation coefficient survive to large h?
(increment 304c)

THE GAP THIS CLOSES. lab_cell_floor_shift.py established two of three
clauses about the cell floor of #113:

  (1) PASS  rho(h) = a S_2(h) + b, with corr(rho, S_2) = 0.9997 to
            1.0000 in all eight bands and |b| < 0.001 -- the shape is
            the singular series of the shift, essentially exactly.
  (3) PASS  permuting the cell labels across N, which preserves every
            cell size and leaves Z byte-identical, collapses the coin's
            B/T from 0.053 to 0.0034 at the bottom band and from 0.038
            to 0.00015 at the top -- i.e. to the (k-1)/n an independent
            sample would give. The floor is caused by the cell-to-
            divisibility correspondence, not by estimation noise.
  (2) FAIL  the predicted B/T = a (E_same[S_2] - E_all[S_2]) came out
            1.39x to 2.51x the measured value, outside the
            pre-registered 1.5x band in seven of eight bands.

The suspect for (2) is named and is mine: a was fitted on h in [2, 80]
and then applied to pair separations of order 10^6. Nothing in that run
tested whether a is the same at h = 10^6 as at h = 10. Hardy-Littlewood
says the SHAPE in h is S_2(h) at any scale, but the coefficient here is
a ratio of a shifted prime-pair count to sqrt(V(N)V(N+h)), and for h
comparable to N the two windows overlap less.

So this run measures a as a function of the scale of h instead of
assuming it.

PRE-REGISTRATION (fixed before the run).

  Compute rho(h) exactly -- again with no Monte Carlo, from
  rho(h) = (mu^2 * g_h)(N) / sqrt(V(N)V(N+h)) averaged over N in band,
  g_h(w) = Lambda(w)Lambda(w+h) -- for ten even h at each of six
  scales h ~ 2*10^1 ... 2*10^6, and fit a separately at each scale.

  RULE A (is there an h-dependence at all).
    a is scale-free iff max(a)/min(a) over the six scales is below 1.15
    in a majority of bands.

  RULE B (does it account for (2)). Recompute the prediction using the
    a measured at the scale nearest the band's mean pair separation,
    which for a band [b, 2b) sampled uniformly is b/3. Predicted B/T is
    then within 1.5x of measured in a majority of bands.

  THE FORK, both outcomes meaningful and stated now.
    A passes  -> a really is scale-free, my extrapolation was sound,
                 and (2)'s failure is a genuine size mismatch: the
                 singular series gives the floor's shape and cause but
                 overpredicts its size by about a factor two, which
                 stays on the books as unexplained.
    A fails and B passes -> the coefficient decays with h, the
                 extrapolation was the fault, and the mechanism
                 accounts for the floor's shape, cause AND size.
    A fails and B fails -> the decay is real but does not close the
                 gap; the mechanism keeps shape and cause only.

  I am not predicting which. The point of writing the fork down is that
  the third outcome is the one I would be tempted to round away.
"""
import math
import time

import numpy as np

QS = [3, 5, 7, 11, 13, 17, 19, 23]
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
    supp = (mu != 0).astype(np.float64)
    F_supp = np.fft.rfft(np.pad(supp, (0, nfft - X - 1)))
    V = np.fft.irfft(F_supp * np.fft.rfft(
        np.pad(lam ** 2, (0, nfft - X - 1))), nfft)[: X + 1]
    print(f"sieve + V  t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i

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
            r = num[Nk] / np.sqrt(V[Nk] * V[Nk + h])
            for i, sel in enumerate(sels):
                s = sel[keep]
                if int(s.sum()) > 1000:
                    rr[j, i] = float(r[s].mean())
        for i in range(nb):
            if not np.isnan(rr[:, i]).any():
                A = np.stack([S, np.ones_like(S)], axis=1)
                a_sc[si, i] = float(
                    np.linalg.lstsq(A, rr[:, i], rcond=None)[0][0])
        print(f"  scale h~{sc}  t={time.time()-t0:.0f}s", flush=True)

    # measured coin B/T and the singular-series excess, as in 304b
    rng = np.random.default_rng(3041)
    idx = np.nonzero(supp)[0]
    R = 8
    meas = np.zeros((R, nb))
    for r in range(R):
        eps = np.zeros(nfft)
        eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        C = np.fft.irfft(np.fft.rfft(eps) * np.fft.rfft(
            np.pad(lam, (0, nfft - X - 1))), nfft)[: X + 1]
        Z = C[Ns] / np.sqrt(V[Ns])
        for i, sel in enumerate(sels):
            z = Z[sel]; k = key[sel]
            uq, inv = np.unique(k, return_inverse=True)
            cnt = np.bincount(inv, minlength=len(uq)).astype(np.float64)
            tot = np.bincount(inv, weights=z, minlength=len(uq))
            gm = float(z.mean())
            B = float((cnt * (tot / cnt - gm) ** 2).sum()) / len(z)
            T = float(((z - gm) ** 2).sum()) / len(z)
            meas[r, i] = B / T
    meas = meas.mean(axis=0)

    exc = np.zeros(nb)
    for i, sel in enumerate(sels):
        Nb = Ns[sel]; kb = key[sel]
        i1 = rng.integers(0, len(Nb), 60_000)
        i2 = rng.integers(0, len(Nb), 60_000)
        ok = i1 != i2
        ha = np.abs(Nb[i1[ok]] - Nb[i2[ok]])
        order = np.argsort(kb, kind="stable")
        ks = kb[order]
        bnd = np.append(np.searchsorted(ks, np.unique(ks), side="left"),
                        len(ks))
        j1, j2 = [], []
        for c in range(len(bnd) - 1):
            lo_c, hi_c = bnd[c], bnd[c + 1]
            if hi_c - lo_c < 2:
                continue
            t = min((hi_c - lo_c) * 4, 20_000)
            u = order[rng.integers(lo_c, hi_c, t)]
            w = order[rng.integers(lo_c, hi_c, t)]
            g = u != w
            j1.append(u[g]); j2.append(w[g])
        hsme = np.abs(Nb[np.concatenate(j1)] - Nb[np.concatenate(j2)])
        sa = float(np.mean([sing(int(h), spf) for h in
                            ha[rng.integers(0, len(ha), 40_000)]]))
        ss = float(np.mean([sing(int(h), spf) for h in
                            hsme[rng.integers(0, len(hsme), 40_000)]]))
        exc[i] = ss - sa

    print(f"\nthe coefficient a in rho(h) = a S_2(h) + b, by scale of h")
    hdr = "".join(f"{'h~'+str(s):>11}" for s in scales)
    print(f"{'band':>21}{hdr}{'max/min':>9}")
    spread = np.zeros(nb)
    for i, (b0, hi, n) in enumerate(lab):
        col = a_sc[:, i]
        fin = col[~np.isnan(col)]
        spread[i] = float(fin.max() / fin.min()) if len(fin) > 1 else np.nan
        cells = "".join(("{:>11.5f}".format(v) if not np.isnan(v)
                         else "{:>11}".format("-")) for v in col)
        print(f"{b0:>9}-{hi:>11}{cells}{spread[i]:>9.2f}")

    print(f"\nprediction at the band's own pair separation (b/3)")
    print(f"{'band':>21} {'sep':>9} {'a used':>9} {'pred B/T':>9} "
          f"{'meas B/T':>9} {'ratio':>7}")
    rat = np.zeros(nb)
    for i, (b0, hi, n) in enumerate(lab):
        sep = b0 / 3.0
        si = int(np.argmin([abs(math.log(s / sep)) for s in scales]))
        col = a_sc[:, i]
        while np.isnan(col[si]) and si > 0:
            si -= 1
        au = col[si]
        pr = au * exc[i]
        rat[i] = pr / meas[i]
        print(f"{b0:>9}-{hi:>11} {sep:>9.0f} {au:>9.5f} {pr:>9.5f} "
              f"{meas[i]:>9.5f} {rat[i]:>7.2f}")

    nA = int((spread[~np.isnan(spread)] < 1.15).sum())
    okA = nA > nb / 2
    nB = int(((rat > 1 / 1.5) & (rat < 1.5)).sum())
    okB = nB > nb / 2
    print(f"\n    RULE A  a scale-free (max/min < 1.15) in a majority: "
          f"{'PASS' if okA else 'FAIL'}  ({nA}/{nb} bands)")
    print(f"    RULE B  prediction at the band's own separation within "
          f"1.5x: {'PASS' if okB else 'FAIL'}  ({nB}/{nb} bands)")
    if okA and okB:
        v = ("a is scale-free and the prediction lands anyway -- the "
             "earlier gap was the pair-separation sampling, not the "
             "coefficient")
    elif okA:
        v = ("a is scale-free, so the extrapolation was sound and the "
             "size gap of clause (2) is real and unexplained; the "
             "singular series of the shift gives the floor's shape and "
             "cause but overpredicts its size")
    elif okB:
        v = ("a decays with h; that decay was the whole of clause (2)'s "
             "failure, and the shift's singular series accounts for the "
             "cell floor's shape, cause and size")
    else:
        v = ("a decays with h but the decay does not close the size gap; "
             "shape and cause stand, size does not")
    print(f"    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
