# -*- coding: utf-8 -*-
"""
The cell floor is the singular series of the shift (increment 304b).

WHERE THIS STANDS. #113 asked why a coin's between-cell variance is
B/T ~ 0.05 and does not fall with sample size. The first attempt at an
answer (lab_cell_floor_mechanism.py, this same increment) FAILED both
of its pre-registered rules, and the failure was informative rather
than fatal: it measured

    rho_bar = [(Sum Z)^2 - Sum Z^2] / (n(n-1) Var Z)  ~  0.2 to 0.46

against B/T ~ 0.06. Those two can never agree, and the fault is in the
test I wrote, not in the data. B is taken about the GRAND mean,

    B = sum_c (n_c/n)(m_c - gm)^2,

so a component common to every N cancels exactly, while rho_bar is
almost entirely that common component. The right identity is

    B/T  ~  rho_same - rho_all,

the EXCESS correlation of same-cell pairs over generic pairs. That is
recorded as its own correction; what follows tests the mechanism the
first run's numbers actually point at.

THE MECHANISM. For a coin, Var(Z(N)) = V(N)/V(N) = 1 exactly, so
covariance is correlation, and averaging over the coin gives it in
CLOSED FORM -- no Monte Carlo at all:

    rho(h) = E_eps[Z(N) Z(N+h)]
           = Sum_v mu^2(v) Lambda(N-v) Lambda(N+h-v)
             / sqrt(V(N) V(N+h))
           = (mu^2 * g_h)(N) / sqrt(V(N)V(N+h)),
      g_h(w) = Lambda(w) Lambda(w+h).

The numerator is a prime-pair count at shift h, so Hardy-Littlewood
puts the singular series S_2(h) = 2 C_2 prod_{p|h, p>2} (p-1)/(p-2)
inside it. And S_2(h) is LARGER when h has small prime factors.

Now the cells: two N in the same cell agree on which q in {3..23}
divide them. For q dividing both, q | h outright. For q dividing
neither, both residues are nonzero and equal with probability
1/(q-1) > 1/q. Either way same-cell pairs have h divisible by small
primes MORE often than generic pairs, hence a larger S_2(h), hence a
larger correlation. That is a floor with no sample-size dependence
whatever, because it is a property of pairs and not of counts.

Note what this says: the coin has no arithmetic in eps, but it still
has arithmetic -- in Lambda, through the shift. That is hazard 7 from
the other side.

PRE-REGISTRATION (fixed before the run).

  (1) SHAPE IN h. Compute rho(h) exactly for every even h in [2, 80],
      per octave band. Fit rho(h) = a S_2(h) + b by least squares.
      RULE: the correlation between rho(h) and S_2(h) across h exceeds
      0.8 in the majority of bands. If rho(h) does not track S_2(h),
      the mechanism is wrong.

  (2) QUANTITATIVE PREDICTION, no free parameter left. Sample pairs
      within a band, same-cell and generic, and form
          predicted B/T = a * (E_same[S_2] - E_all[S_2]).
      The intercept b cancels in the difference, so a alone fixes it.
      RULE: predicted B/T is within a factor 1.5 of the measured coin
      B/T in the majority of bands. This is the real test: a is fitted
      to rho(h) at small h and then used to predict a completely
      different statistic.

  (3) PLACEBO KEY. Re-run the coin's B/T with the cell labels randomly
      permuted across N. That preserves every cell size exactly and
      leaves Z byte-identical; it destroys only the correspondence
      between a cell and the divisibility of N, which is the one thing
      the mechanism says matters.
      RULE: placebo B/T falls to at most a quarter of the real-key
      value. If the floor survives the placebo it is estimation noise
      after all and the mechanism is wrong.

  WHAT WOULD REFUTE. Any of the three failing. They are close to
  independent: (1) is the shape, (2) the size, (3) the cause.
"""
import math
import time

import numpy as np

QS = [3, 5, 7, 11, 13, 17, 19, 23]
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
    return mu, lam, spf


def sing(h, spf):
    """S_2(h) = 2 C_2 prod_{p|h, p>2} (p-1)/(p-2), for even h."""
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


def sing_arr(hs, spf):
    return np.array([sing(int(h), spf) for h in hs])


def main():
    X = 16_000_000
    lo = 100_000
    HMAX = 80
    t0 = time.time()
    mu, lam, spf = sieve(X)
    nfft = 1
    while nfft < 2 * (X + 1):
        nfft *= 2
    supp = (mu != 0).astype(np.float64)

    F_supp = np.fft.rfft(np.pad(supp, (0, nfft - X - 1)))
    lam2 = np.pad(lam ** 2, (0, nfft - X - 1))
    V = np.fft.irfft(F_supp * np.fft.rfft(lam2), nfft)[: X + 1]
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
            lab.append((b, hi, int(sel.sum()), math.log(math.sqrt(b * hi))))
        b = hi
    nb = len(sels)

    # ---- (1) exact rho(h) per band, no Monte Carlo ----
    hs = np.arange(2, HMAX + 1, 2)
    rho = np.zeros((len(hs), nb))
    pad = np.zeros(nfft)
    for j, h in enumerate(hs):
        g = pad.copy()
        g[: X + 1 - h] = lam[: X + 1 - h] * lam[h:]
        num = np.fft.irfft(F_supp * np.fft.rfft(g), nfft)[: X + 1]
        Nk = Ns[Ns + h <= X]
        r = num[Nk] / np.sqrt(V[Nk] * V[Nk + h])
        for i, sel in enumerate(sels):
            s = sel[: len(Nk)]
            rho[j, i] = float(r[s].mean())
        if (j + 1) % 10 == 0:
            print(f"  h={h}  t={time.time()-t0:.0f}s", flush=True)

    S = sing_arr(hs, spf)
    A = np.stack([S, np.ones_like(S)], axis=1)
    coef = np.linalg.lstsq(A, rho, rcond=None)[0]
    a_fit, b_fit = coef[0], coef[1]
    corr = np.array([float(np.corrcoef(S, rho[:, i])[0, 1]) for i in range(nb)])

    # ---- (2) singular-series excess of same-cell pairs ----
    rng = np.random.default_rng(3040)
    exc = np.zeros(nb)
    for i, sel in enumerate(sels):
        Nb = Ns[sel]
        kb = key[sel]
        m = min(len(Nb), 200_000)
        i1 = rng.integers(0, len(Nb), m)
        i2 = rng.integers(0, len(Nb), m)
        ok = i1 != i2
        h_all = np.abs(Nb[i1[ok]] - Nb[i2[ok]])
        order = np.argsort(kb, kind="stable")
        ks = kb[order]
        bnd = np.searchsorted(ks, np.unique(ks), side="left")
        bnd = np.append(bnd, len(ks))
        j1, j2 = [], []
        for c in range(len(bnd) - 1):
            lo_c, hi_c = bnd[c], bnd[c + 1]
            sz = hi_c - lo_c
            if sz < 2:
                continue
            t = min(sz * 4, 40_000)
            u = order[rng.integers(lo_c, hi_c, t)]
            w = order[rng.integers(lo_c, hi_c, t)]
            g = u != w
            j1.append(u[g]); j2.append(w[g])
        h_same = np.abs(Nb[np.concatenate(j1)] - Nb[np.concatenate(j2)])
        sa = float(np.mean([sing(int(h), spf) for h in
                            h_all[rng.integers(0, len(h_all), 60_000)]]))
        ss = float(np.mean([sing(int(h), spf) for h in
                            h_same[rng.integers(0, len(h_same), 60_000)]]))
        exc[i] = ss - sa
    pred = a_fit * exc

    # ---- (3) measured coin B/T, real key and placebo key ----
    def bt(Z, kk, sel):
        z = Z[sel]; k = kk[sel]
        uq, inv = np.unique(k, return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uq)).astype(np.float64)
        tot = np.bincount(inv, weights=z, minlength=len(uq))
        gm = float(z.mean())
        B = float((cnt * (tot / cnt - gm) ** 2).sum()) / len(z)
        T = float(((z - gm) ** 2).sum()) / len(z)
        return B / T

    R = 8
    idx = np.nonzero(supp)[0]
    real_bt = np.zeros((R, nb))
    plac_bt = np.zeros((R, nb))
    for r in range(R):
        eps = np.zeros(nfft)
        eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        C = np.fft.irfft(np.fft.rfft(eps) * np.fft.rfft(
            np.pad(lam, (0, nfft - X - 1))), nfft)[: X + 1]
        Z = C[Ns] / np.sqrt(V[Ns])
        kp = key.copy()
        rng.shuffle(kp)
        for i, sel in enumerate(sels):
            real_bt[r, i] = bt(Z, key, sel)
            plac_bt[r, i] = bt(Z, kp, sel)
        print(f"  coin {r+1}/{R}  t={time.time()-t0:.0f}s", flush=True)
    meas = real_bt.mean(axis=0)
    plac = plac_bt.mean(axis=0)

    # ---- report ----
    print(f"\nexact coin correlation rho(h) against the singular series")
    print(f"{'band':>21} {'corr(rho,S2)':>13} {'a':>10} {'b':>10} "
          f"{'a*S2bar':>9} {'1/logN':>8}")
    Sbar = float(S.mean())
    for i, (b0, hi, n, Lm) in enumerate(lab):
        print(f"{b0:>9}-{hi:>11} {corr[i]:>13.4f} {a_fit[i]:>10.5f} "
              f"{b_fit[i]:>10.5f} {a_fit[i]*Sbar:>9.5f} {1.0/Lm:>8.5f}")

    print(f"\npredicted vs measured cell floor")
    print(f"{'band':>21} {'E_same-E_all':>13} {'pred B/T':>9} "
          f"{'meas B/T':>9} {'ratio':>7} {'placebo':>9} {'plac/meas':>10}")
    for i, (b0, hi, n, Lm) in enumerate(lab):
        print(f"{b0:>9}-{hi:>11} {exc[i]:>13.5f} {pred[i]:>9.5f} "
              f"{meas[i]:>9.5f} {pred[i]/meas[i]:>7.2f} {plac[i]:>9.5f} "
              f"{plac[i]/meas[i]:>10.3f}")

    n1 = int((corr > 0.8).sum())
    ok1 = n1 > nb / 2
    rat = pred / meas
    n2 = int(((rat > 1 / 1.5) & (rat < 1.5)).sum())
    ok2 = n2 > nb / 2
    n3 = int((plac / meas <= 0.25).sum())
    ok3 = n3 == nb

    print(f"\n    (1) corr(rho, S_2) > 0.8 in a majority: "
          f"{'PASS' if ok1 else 'FAIL'}  ({n1}/{nb} bands)")
    print(f"    (2) predicted B/T within 1.5x of measured in a majority: "
          f"{'PASS' if ok2 else 'FAIL'}  ({n2}/{nb} bands)")
    print(f"    (3) placebo key at most a quarter of the real key: "
          f"{'PASS' if ok3 else 'FAIL'}  ({n3}/{nb} bands)")
    good = sum([ok1, ok2, ok3])
    if good == 3:
        v = ("#113 is answered: the cell floor is the singular series "
             "of the shift. Same-cell pairs have h divisible by small "
             "primes more often, S_2(h) is larger there, and the excess "
             "correlation is a property of PAIRS, so no sample size "
             "removes it.")
    elif good == 2:
        v = ("two of three clauses hold; the mechanism is indicated but "
             "not established, and #113 stays open with a named "
             "candidate")
    else:
        v = ("the shift's singular series does not explain the floor; "
             "#113 remains open")
    print(f"    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
