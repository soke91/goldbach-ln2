# -*- coding: utf-8 -*-
"""
The mask's share, per cell instead of pooled (increment 318)

WHY. `OPEN_QUESTIONS.md` A4 is the mask's share of the variance, and
what it says is not "the number is wrong" but "**a share statistic that
is not size-weighted. None is defined yet.**" #67's figures -- 14.1% at
1e5, 1.15% at 1.6e7 -- were withdrawn because they came from a
between-cell variance `B = sum_c (n_c/n)(m_c - m)^2`, which weights
each cell by how many N are in it, and the mask lives in the RARE
cells: the largest-|z| cell holds 6.6e-5 of the top band (#118).

The two statistics answer different questions and only one has ever
been computed.

  size-weighted   how much of the fluctuation of a RANDOM even N is
                  the mask? Small at large N, and correctly so --
                  a random N is in no deep cell.
  per cell        how much of the fluctuation of an N IN THAT CELL is
                  the mask? Never computed.

The second is the one every use of the mask in this program actually
needs, since the mask is invoked precisely when N has many small prime
factors.

THE STATISTIC. In units of sqrt(V), write d_c = m_c - gm for the cell's
mean offset and rho = Var(Z) for the band. The fluctuation of an N in
cell c has a deterministic part d_c and a random part of variance rho,
so

    share_c = d_c^2 / (d_c^2 + rho).

d_c is measured with error, so d_c^2 is biased upward by exactly the
variance of the estimator -- which increment 305 made available in
closed form, sd_c^2 = Q_cc/n_c^2 - 2Q_ca/(n_c n) + Q_aa/n^2. The
debiased square is d_c^2 - sd_c^2, and rule (T1) is what shows the
debiasing works.

PRE-REGISTRATION (fixed before the run).

  (T1) SELF-TEST ON A COIN, which must pass or nothing else reads.
       A coin has no mask, so its debiased d_c^2 must scatter about
       zero. RULE: over 16 coin draws and all cells, the mean debiased
       square is within 3 standard errors of zero, and the implied
       share is below 0.02.

  (T2) THE PER-CELL SHARE for the real mu, reported per band for the
       deepest cell and as a distribution over cells. No threshold --
       this is the number A4 says does not exist.

  (T3) THE SIZE-WEIGHTED SHARE computed alongside, to show the two
       are different questions rather than one of them being wrong.
       RULE: the size-weighted share falls with N while the deepest
       cell's per-cell share does not fall as fast. If they behave the
       same way there is only one statistic here and A4 is empty.

  WHAT WOULD REFUTE. (T1) failing means the debiasing is wrong and
  every share below is inflated. (T3) failing would mean the pooled
  figure was an adequate summary after all.
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
R = 16


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
    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    invV = 1.0 / np.sqrt(V[Ns])
    Creal = np.fft.irfft(np.fft.rfft(
        np.pad(mu.astype(np.float64), (0, nfft - X - 1))) * F_lam,
        nfft)[: X + 1]
    Zreal = Creal[Ns] * invV
    muw = suppf[: X + 1]
    print(f"sieve + V + C  t={time.time()-t0:.0f}s", flush=True)

    sels, lab = [], []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            sels.append(sel)
            lab.append((b, hi, int(sel.sum())))
        b = hi

    SD, NC, KEYS = [], [], []
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
            us.append(u)
            ns.append(int(m.sum()))
            u_all += u
        mu_all = muw * u_all
        Qaa = float(np.dot(mu_all, u_all))
        sd = np.zeros(len(cells))
        for j, u in enumerate(us):
            nc = ns[j]
            var = (float(np.dot(muw * u, u)) / nc ** 2
                   - 2 * float(np.dot(mu_all, u)) / (nc * n)
                   + Qaa / n ** 2)
            sd[j] = math.sqrt(max(var, 0.0))
        SD.append(sd)
        NC.append(np.array(ns, dtype=np.float64))
        KEYS.append(cells)
        print(f"  band {i+1}/{len(sels)} exact  t={time.time()-t0:.0f}s",
              flush=True)

    def cellstats(Z, i):
        sel = sels[i]
        z = Z[sel]
        k = key[sel]
        uq, inv = np.unique(k, return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uq)).astype(np.float64)
        tot = np.bincount(inv, weights=z, minlength=len(uq))
        gm = float(z.mean())
        d = tot / cnt - gm
        rho = float(z.var())
        return d, rho, cnt, gm

    # ---- (T1) coin self-test ----
    idx = np.nonzero(supp)[0]
    rng = np.random.default_rng(318)
    deb = [[] for _ in sels]
    for r in range(R):
        eps = np.zeros(nfft)
        eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        C = np.fft.irfft(np.fft.rfft(eps) * F_lam, nfft)[: X + 1]
        Z = C[Ns] * invV
        for i in range(len(sels)):
            d, rho, cnt, gm = cellstats(Z, i)
            deb[i].append(d ** 2 - SD[i] ** 2)
    allc = np.concatenate([np.concatenate(x) for x in deb])
    mdeb = float(allc.mean())
    sdeb = float(allc.std(ddof=1)) / math.sqrt(len(allc))
    zc = mdeb / sdeb
    share_coin = max(mdeb, 0.0) / (max(mdeb, 0.0) + 1.0)
    okT1 = abs(zc) < 3.0 and share_coin < 0.02
    print(f"\n(T1) coin: mean debiased d^2 over {len(allc)} cell-draws "
          f"= {mdeb:+.5f} +/- {sdeb:.5f}  (z = {zc:+.2f}),")
    print(f"     implied share {share_coin:.4f}  ->  "
          f"{'PASS' if okT1 else 'FAIL'}")

    # ---- (T2)(T3) the real mask ----
    print(f"\n(T2)(T3) the mask's share, per cell and pooled")
    print(f"{'band':>21} {'rho':>6} {'deepest d':>10} {'share':>8} "
          f"{'cells>50%':>10} {'N in them':>11} {'pooled':>9}")
    deep_key = (1 << len(QS)) - 1
    rows = []
    for i, (b0, hi, n) in enumerate(lab):
        d, rho, cnt, gm = cellstats(Zreal, i)
        db = d ** 2 - SD[i] ** 2
        sh = np.maximum(db, 0.0) / (np.maximum(db, 0.0) + rho)
        j = int(np.nonzero(KEYS[i] == deep_key)[0][0])
        big = sh > 0.5
        frac = float(cnt[big].sum() / cnt.sum())
        pooled = float((cnt * d ** 2).sum() / cnt.sum())
        pooled = pooled / (pooled + rho)
        rows.append((b0, hi, rho, d[j], sh[j], int(big.sum()), frac,
                     pooled))
        print(f"{b0:>9}-{hi:>11} {rho:>6.3f} {d[j]:>+10.3f} "
              f"{sh[j]:>8.3f} {int(big.sum()):>10} {frac:>11.2e} "
              f"{pooled:>9.4f}")

    ds = np.array([r[4] for r in rows])
    ps = np.array([r[7] for r in rows])
    fall_deep = ds[0] / ds[-1]
    fall_pool = ps[0] / ps[-1]
    okT3 = fall_pool > fall_deep
    print(f"\n    (T3) the pooled share falls faster than the deepest "
          f"cell's: {'PASS' if okT3 else 'FAIL'}")
    print(f"         pooled  {ps[0]:.4f} -> {ps[-1]:.4f}  "
          f"(factor {fall_pool:.1f})")
    print(f"         deepest {ds[0]:.4f} -> {ds[-1]:.4f}  "
          f"(factor {fall_deep:.2f})")
    if okT1 and okT3:
        v = (f"A4 is answered. For a RANDOM even N the mask is "
             f"{ps[-1]:.2%} of the fluctuation at the top band, which "
             f"is what #67 was measuring and why it looked negligible. "
             f"For an N in the deepest cell it is {ds[-1]:.1%} -- the "
             f"mask is not small, it is RARE, and the pooled statistic "
             f"cannot tell those apart")
    elif okT1:
        v = ("the two statistics behave alike, so the pooled figure was "
             "an adequate summary and A4 is empty")
    else:
        v = ("the debiasing does not zero out on a coin; every share "
             "here is inflated and none of it reads")
    print(f"    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
