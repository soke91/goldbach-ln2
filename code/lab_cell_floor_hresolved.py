# -*- coding: utf-8 -*-
"""
The cell floor's size, with the single-a approximation removed
(increment 304e).

WHAT IS ALREADY SETTLED ABOUT #113, AND WHAT IS NOT.

  SHAPE   -- settled (304b). For a coin, rho(h) = a S_2(h) + b exactly
             in the sense that corr(rho, S_2) = 0.9997 to 1.0000 across
             all eight bands with |b| < 0.001.
  CAUSE   -- settled (304b). Permuting the cell labels across N, which
             preserves every cell size and leaves Z byte-identical,
             collapses B/T from 0.053 to 0.0034 at the bottom band and
             from 0.038 to 0.00015 at the top: the (k-1)/n of an
             independent sample. The floor is the cell-to-divisibility
             correspondence, not estimation noise.
  SIZE    -- NOT settled. 304d put the target in closed form (verified
             against 12 coin draws, max |z| = 1.04) and the prediction
             a(sep) * (E_same[S_2] - E_all[S_2]) overshot it by 1.13x
             to 1.87x, failing in five of eight bands.

AND ONE THING THAT HAS TO BE SAID ABOUT 304c. Its RULE B "passed 7 of 8
bands within 1.5x". 304d measured the per-draw spread of the statistic
that rule was judging: 86.9%. With eight draws that is a standard error
near 31%, so a 1.5x tolerance could not have discriminated anything.
The pass was noise. It is recorded as a correction, and it is why 304d
replaced the measured target with the closed form.

WHAT THIS RUN CHANGES. Exactly one thing: the prediction used a single
coefficient a, evaluated at the band's mean pair separation b/3. But
304c also showed a is NOT scale-free -- it is flat for h up to about
N/10 and falls away above that. Pairs inside a band have separations
spread over the whole band, so collapsing that to one a is an
approximation, and since a decays it is an approximation that
overstates. The model's actual content is

    B/T  ~  E_same[rho(h)] - E_all[rho(h)],
    rho(h) = a(h) S_2(h) + b(h),

with a and b interpolated in log h from seven measured scales. Nothing
else moves: same closed-form target, same sampled pairs, same key.

HONESTY ABOUT ORDER. This refinement is being run AFTER seeing (P)
fail, which is the classic way to manufacture agreement. Three things
guard it. The tolerance stays exactly where it was, 1.5x, fixed before
304d. The refinement was named in 304c's own output as the untested
step, before this result existed. And it removes an approximation
rather than adding a parameter -- a(h) is measured, not fitted to the
target.

PRE-REGISTRATION.
  (S) self-test: the closed form must again sit within 3 se of R = 12
      coin draws.
  (R) RULE: the h-resolved prediction is within 1.5x of the closed form
      in a majority of bands, AND its ratio does not drift with N by
      more than the single-a version's drift (1.31 -> 1.87). The second
      clause matters: a refinement that fixes the level but keeps the
      trend has not explained the trend.

  WHAT WOULD REFUTE. Either clause failing leaves #113 with shape and
  cause established and size open, which is a perfectly reportable
  state and the one I will report if that is what comes out.
"""
import math
import time

import numpy as np

QS = [3, 5, 7]
C2 = 0.66016181584686957392
NPER = 10
SCALES = [20, 200, 2000, 20000, 200000, 2000000, 6000000]


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


def sing_vec(hs, spf):
    out = np.empty(len(hs))
    for i, h in enumerate(hs):
        s = 2.0 * C2
        m = int(h)
        while m % 2 == 0:
            m //= 2
        while m > 1:
            p = int(spf[m])
            s *= (p - 1.0) / (p - 2.0)
            while m % p == 0:
                m //= p
        out[i] = s
    return out


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
    F_lam = np.fft.rfft(np.pad(lam, (0, nfft - X - 1)))
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

    # ---- closed-form target ----
    exact = np.zeros(nb)
    for i, sel in enumerate(sels):
        n = int(sel.sum())
        kb = key[sel]
        u_all = np.zeros(X + 1)
        parts = []
        for c in np.unique(kb):
            m = sel.copy()
            m[sel] = (kb == c)
            w = np.zeros(nfft)
            w[Ns] = np.where(m, invV, 0.0)
            u = np.fft.irfft(Fl_c * np.fft.rfft(w), nfft)[: X + 1]
            parts.append((int(m.sum()), u))
            u_all += u
        Egm2 = float(np.dot(muw * u_all, u_all)) / n ** 2
        EB = 0.0
        for nc, u in parts:
            EB += (nc / n) * (float(np.dot(muw * u, u)) / nc ** 2
                              - 2 * float(np.dot(muw * u, u_all)) / (nc * n)
                              + Egm2)
        exact[i] = EB / (1.0 - Egm2)
        print(f"  band {i+1}/{nb} exact B/T = {exact[i]:.6f}  "
              f"t={time.time()-t0:.0f}s", flush=True)

    # ---- a(h), b(h) at seven scales ----
    a_sc = np.full((len(SCALES), nb), np.nan)
    b_sc = np.full((len(SCALES), nb), np.nan)
    pad = np.zeros(nfft)
    for si, sc in enumerate(SCALES):
        hs = [sc + 2 * t for t in range(NPER)]
        S = sing_vec(hs, spf)
        A = np.stack([S, np.ones_like(S)], axis=1)
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
        for i in range(nb):
            if not np.isnan(rr[:, i]).any():
                cf = np.linalg.lstsq(A, rr[:, i], rcond=None)[0]
                a_sc[si, i], b_sc[si, i] = float(cf[0]), float(cf[1])
        print(f"  scale h~{sc}  t={time.time()-t0:.0f}s", flush=True)

    logS = np.log(np.array(SCALES, dtype=np.float64))

    # ---- h-resolved prediction ----
    rng = np.random.default_rng(3043)
    pred_h = np.zeros(nb)
    pred_1 = np.zeros(nb)
    for i, sel in enumerate(sels):
        Nb = Ns[sel]
        kb = key[sel]
        fin = ~np.isnan(a_sc[:, i])
        xs, ays, bys = logS[fin], a_sc[fin, i], b_sc[fin, i]

        def rho_of(h):
            lh = np.log(np.maximum(h, 2.0))
            a = np.interp(lh, xs, ays)
            bb = np.interp(lh, xs, bys)
            return a * sing_vec(h, spf) + bb

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
            a0, b0 = bnd[c], bnd[c + 1]
            if b0 - a0 < 2:
                continue
            t = min((b0 - a0) * 4, 60_000)
            u = order[rng.integers(a0, b0, t)]
            w = order[rng.integers(a0, b0, t)]
            g = u != w
            j1.append(u[g]); j2.append(w[g])
        hs_ = np.abs(Nb[np.concatenate(j1)] - Nb[np.concatenate(j2)])
        ha = ha[rng.integers(0, len(ha), 50_000)]
        hs_ = hs_[rng.integers(0, len(hs_), 50_000)]
        pred_h[i] = float(rho_of(hs_).mean() - rho_of(ha).mean())
        sep = lab[i][0] / 3.0
        a1 = float(np.interp(math.log(sep), xs, ays))
        pred_1[i] = a1 * float(sing_vec(hs_, spf).mean()
                               - sing_vec(ha, spf).mean())

    # ---- (S) self-test ----
    R = 12
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
    mmean, msd = mc.mean(axis=0), mc.std(axis=0, ddof=1)
    mse = msd / math.sqrt(R)
    zs = (mmean - exact) / mse
    okS = bool((np.abs(zs) < 3).all())

    # ---- report ----
    print(f"\na(h) by scale, key from {QS}")
    print(f"{'band':>21}" + "".join(f"{'h~'+str(s):>10}" for s in SCALES))
    for i, (b0, hi, n) in enumerate(lab):
        cells = "".join(("{:>10.5f}".format(v) if not np.isnan(v)
                         else "{:>10}".format("-")) for v in a_sc[:, i])
        print(f"{b0:>9}-{hi:>11}{cells}")

    print(f"\nsize of the cell floor: closed form against both predictions")
    print(f"{'band':>21} {'exact':>9} {'single-a':>9} {'ratio':>7} "
          f"{'h-resolved':>11} {'ratio':>7} {'MC z':>6}")
    r1 = pred_1 / exact
    rh = pred_h / exact
    for i, (b0, hi, n) in enumerate(lab):
        print(f"{b0:>9}-{hi:>11} {exact[i]:>9.5f} {pred_1[i]:>9.5f} "
              f"{r1[i]:>7.2f} {pred_h[i]:>11.5f} {rh[i]:>7.2f} "
              f"{zs[i]:>6.2f}")

    nR = int(((rh > 1 / 1.5) & (rh < 1.5)).sum())
    ok1 = nR > nb / 2
    drift_h = float(rh.max() / rh.min())
    drift_1 = float(r1.max() / r1.min())
    ok2 = drift_h <= drift_1
    print(f"\n    (S) closed form within 3 se of {R} coin draws: "
          f"{'PASS' if okS else 'FAIL'}  (max |z| = {np.abs(zs).max():.2f})")
    print(f"    (R1) h-resolved within 1.5x in a majority: "
          f"{'PASS' if ok1 else 'FAIL'}  ({nR}/{nb} bands)")
    print(f"    (R2) ratio drift not worse than single-a: "
          f"{'PASS' if ok2 else 'FAIL'}  "
          f"(h-resolved {drift_h:.2f}x, single-a {drift_1:.2f}x)")
    if okS and ok1 and ok2:
        v = ("#113 is answered in full. The cell floor is the singular "
             "series of the shift -- shape, cause and size. Same-cell "
             "pairs have h divisible by small primes more often, "
             "S_2(h) is larger there, and the excess is a property of "
             "PAIRS, so no sample size removes it.")
    elif okS and (ok1 or ok2):
        v = ("the h-resolved model improves the size but does not meet "
             "both clauses; shape and cause stand, size is indicated "
             "and not established")
    elif okS:
        v = ("removing the single-a approximation does not close the "
             "size gap; shape and cause stand, size stays open")
    else:
        v = "the closed form fails its self-test; nothing here reads"
    print(f"    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
