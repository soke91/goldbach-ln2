# -*- coding: utf-8 -*-
"""
Does a coin reproduce the zeta LINES, as opposed to the aggregate R^2?

WHY THIS RUN EXISTS. `audit_coin_control_v1claims.py` reproduced
correction #110's aggregate comparison closely (real R^2 = 4.64e-3,
coin mean 3.27e-3, 4 of 20 coin draws at or above real, ratio 1.42x
against #110's 1.30x), and #110 concluded from it that "the lines are
Lambda's" and withdrew the claim as "empty, restating that Lambda is
in the convolution".

But the same run found, on the LOCAL-background test, real 6/10
ordinates above their own local 99th percentile against 0/10 for a
single coin draw. Those two readings point opposite ways, and the
reason they can is that the aggregate R^2 at these frequencies is
dominated by the BROADBAND low-frequency power that any autocorrelated
series has -- real and coin alike. A statistic dominated by a
background both fields share cannot decide whether one of them has a
line on top of it. The local test subtracts that background; the
aggregate one does not.

One coin draw is not evidence. This script runs the local test on the
real field and on eight independent coins.

PRE-REGISTRATION (written before the run).

  (1) RULE. If the coin hit counts are distributed around the real
      one, #110 is right and the lines are Lambda's. If the coin hit
      counts cluster near the chance level (0-1 of 10 at the 1% level)
      while the real field gives 6, then a coin does NOT reproduce the
      lines, #110's withdrawal rests on a statistic that cannot see
      them, and the paper's item 4 is -- as a qualitative statement --
      correct after all.
  (2) Reported either way: the aggregate R^2 of each coin, so the two
      statistics can be read side by side on the same draws.
  (3) PREDICTION. My round-5 run and the single coin above both point
      the same way, so I predict the coins cluster at 0-2 of 10. I
      note that this contradicts my own report of #110 as "confirmed",
      and that the aggregate half of #110 IS confirmed -- both can be
      true, and the distinction is exactly which statistic is used.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

QS = (3, 5, 7, 11, 13)
GAMMAS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
          37.586178, 40.918719, 43.327073, 48.005151, 49.773832]
NCTRL = 40
NCOIN = 8


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


def main():
    X = 4_000_000
    LO = 100_000
    mu, lam = sieve_mu_lambda(X)
    supp = (mu != 0).astype(np.float64)
    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2
    FL = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))
    V = np.fft.irfft(np.fft.rfft(np.pad(supp, (0, nf - X - 1)))
                     * np.fft.rfft(np.pad(lam ** 2, (0, nf - X - 1))),
                     nf)[: X + 1]

    def wall(s):
        return np.fft.irfft(
            np.fft.rfft(np.pad(s, (0, nf - X - 1))) * FL, nf)[: X + 1]

    Ns = np.arange(LO + LO % 2, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    uniq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
    bands = []
    b = LO
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 20000:
            bands.append(sel)
        b = hi

    def pipeline(Cfull):
        Cv = Cfull[Ns]
        tot = np.bincount(inv, weights=Cv, minlength=len(uniq))
        Z = (Cv - (tot / cnt)[inv]) / np.sqrt(V[Ns])
        out = Z.copy()
        for sel in bands:
            out[sel] = (Z[sel] - Z[sel].mean()) / Z[sel].std(ddof=1)
        return out - out.mean()

    L = np.log(Ns.astype(np.float64))
    # precompute the control frequencies once, so every field is judged
    # against the same set
    rng2 = np.random.default_rng(4242)
    ctrls = []
    for t in GAMMAS:
        c = []
        while len(c) < NCTRL:
            x = float(rng2.uniform(t - 4.0, t + 4.0))
            if min(abs(x - u) for u in GAMMAS) < 1.0:
                continue
            c.append(x)
        ctrls.append(c)
    freqs = [t for t in GAMMAS] + [x for c in ctrls for x in c]
    CS = np.empty((len(freqs), 2, len(L)))
    for i, f in enumerate(freqs):
        c = np.cos(f * L)
        s = np.sin(f * L)
        CS[i, 0] = c - c.mean()
        CS[i, 1] = s - s.mean()

    def amps(z):
        out = np.empty(len(freqs))
        for i in range(len(freqs)):
            c, s = CS[i]
            G = np.array([[c @ c, c @ s], [c @ s, s @ s]])
            co = np.linalg.solve(G, np.array([c @ z, s @ z]))
            out[i] = math.hypot(co[0], co[1])
        return out

    cols = []
    for t in GAMMAS:
        cols.append(np.cos(t * L))
        cols.append(np.sin(t * L))
    B = np.stack(cols, axis=1)
    B -= B.mean(axis=0)
    Ginv = np.linalg.inv(B.T @ B)

    def r2_of(z):
        bty = B.T @ z
        return float((Ginv @ bty) @ bty / (z @ z))

    def hits_of(z):
        a = amps(z)
        h = 0
        for j in range(len(GAMMAS)):
            lo = len(GAMMAS) + j * NCTRL
            p99 = float(np.quantile(a[lo:lo + NCTRL], 0.99))
            h += a[j] > p99
        return h

    Zr = pipeline(wall(mu))
    print("Local-background test: does a coin reproduce the LINES?")
    print(f"X = {X}, {NCTRL} local controls per ordinate, "
          f"{NCOIN} coin draws")
    print()
    print(f"{'field':>14} {'aggregate R^2':>15} {'ordinates above':>17}")
    print(f"{'':>14} {'':>15} {'local p99, of 10':>17}")
    print("-" * 48)
    print(f"{'real mu':>14} {r2_of(Zr):>15.4e} {hits_of(Zr):>17d}")

    rng = np.random.default_rng(90210)
    idx = np.nonzero(supp)[0]
    hs, rs = [], []
    for k in range(NCOIN):
        eps = np.zeros(X + 1)
        eps[idx] = rng.choice([-1.0, 1.0], size=len(idx))
        z = pipeline(wall(eps))
        h, r = hits_of(z), r2_of(z)
        hs.append(h)
        rs.append(r)
        print(f"{'coin ' + str(k + 1):>14} {r:>15.4e} {h:>17d}")
    print()
    print(f"    coin aggregate R^2: mean {np.mean(rs):.4e}, "
          f"max {np.max(rs):.4e}   (real {r2_of(Zr):.4e})")
    print(f"    coin hits: {hs}, mean {np.mean(hs):.2f}"
          f"   (real {hits_of(Zr)})")
    print(f"    chance expectation at the 1% level: 0.1 of 10")
    print()
    print("(1) RULE was: if the coin hit counts sit around the real")
    print("    one, #110 is right and the lines are Lambda's. If they")
    print("    cluster at chance while the real field does not, the")
    print("    aggregate statistic #110 used cannot see the lines.")
    print("DONE")


if __name__ == "__main__":
    main()
