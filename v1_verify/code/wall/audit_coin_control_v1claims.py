# -*- coding: utf-8 -*-
"""
The coin control of Lemma 18 (`lem:coin`), applied to the two figures
of `v1/paper/wall_v1.tex` that its own correction record withdrew.

WHAT LEMMA 18 SAYS, verbatim from the paper:

    Let eps(v) = +-1 be arbitrary signs on {v : mu(v) != 0} and zero
    elsewhere. Then eps^2 = mu^2 pointwise, so V(N) is unchanged, and
    the field C_eps(N) = sum_v eps(v) Lambda(N-v) has the same exact
    second moment as C(N) for every N. Consequently ANY ESTIMATOR
    WHOSE OUTPUT IS REPRODUCED WHEN mu IS REPLACED BY eps IS NOT
    MEASURING mu.

WHAT `v1_log/docs/CLOSURE_REAUDIT.md` RECORDS.

  #106: "every measurement of rho = Var C/V has gone through one
        pipeline ... it does, and it accounts for the whole effect.
        Substituting eps(v) = +-1 at random on {mu != 0} leaves V
        identical and every other step byte-identical, so a coin must
        give rho = 1. It gives 0.761, 0.780, 0.787, 0.792, 0.813,
        0.828, 0.851, 0.860 -- reproducing the real curve with z
        between -0.5 and +0.4 in every band. THE CENTRED ESTIMATOR
        CANNOT TELL mu FROM A COIN. Withdrawn: ... and the
        QUANTITATIVE half of #86, which compared a reconstruction
        against 'the measured -0.18'."

  #110: "#94 and #96 read the zeta ordinates in C(N) as a property of
        THE WALL, against a permutation null | THE LINES ARE
        Lambda's. Replacing mu by a random +-1 on the same support,
        through the same Lambda and the identical pipeline: real
        R^2 = 3.896e-3 against a coin mean of 2.994e-3, coin max
        5.515e-3, and 6 of 20 coin draws at or above the real value.
        The ratio is 1.30x, not the 1566x reported. ... Withdrawn:
        #94's 'the wall's fluctuation is Gaussian in distribution but
        not phase-random in log N' as a property of the wall -- it
        remains literally true and is EMPTY, restating that Lambda is
        in the convolution; and #96's 0.39% share."

WHAT THE PAPER NEVERTHELESS SAYS.

  Section `sec:coin`: "Reconstructing rho-1 from Proposition W gives
  -0.0976 against a measured -0.18, a factor 0.54" -- the withdrawn
  quantitative half of #86, with its withdrawn target.

  `conj:wall` item 4: "G is Gaussian in distribution but not
  phase-random in log N. Regressing G on cos(gamma log N),
  sin(gamma log N) ... gives R^2 = 3.90e-3 against a 200-surrogate
  maximum of 5.09e-6, and every ordinate individually at z >= 23. ...
  The 0.39% is a floor, not the share" -- #94's withdrawn sentence and
  #96's withdrawn number, verbatim.

THIS SCRIPT re-runs the control from the statement of Lemma 18, on
both figures, without reading v1's coin-control code.

PRE-REGISTRATION (written before the run).

  (1) RULE for rho. If a coin reproduces the measured rho band by band
      -- within, say, 3 standard errors of the real value in every
      band -- then by Lemma 18 the "measured -0.18" is not a
      measurement of mu, and comparing Proposition 15's reconstruction
      against it establishes nothing about mu.
  (2) RULE for the zeta regression. If coin draws reach the real R^2,
      the spectral claim is a statement about Lambda and not about the
      wall's mu-fluctuation.
  (3) SECOND TEST on the zeta claim, because my own round-5 run found
      6 of 10 ordinates clearing a LOCAL background null and concluded
      "the effect is real". That test cannot attribute the effect. So
      the local-background test is re-run ON A COIN. If the coin also
      clears 6-ish of 10 locally, my round-5 conclusion stands only as
      a statement about Lambda, and item 4 is empty exactly as #110
      says.
  (4) PREDICTION, recorded so it cannot be reported as a surprise.
      C = mu * Lambda, and Lambda carries the zeta zeros through the
      explicit formula, so I predict the coin reproduces both figures
      and that #106 and #110 are confirmed. I expect my round-5
      "the effect is real" to survive only in the empty form.
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
    L2 = np.fft.rfft(np.pad(lam ** 2, (0, nf - X - 1)))
    V = np.fft.irfft(np.fft.rfft(np.pad(supp, (0, nf - X - 1))) * L2,
                     nf)[: X + 1]

    def wall(signs):
        return np.fft.irfft(
            np.fft.rfft(np.pad(signs, (0, nf - X - 1))) * FL, nf)[: X + 1]

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
            bands.append((b, hi, sel))
        b = hi

    def pipeline(Cfull):
        """the paper's pipeline: de-mask by cells, standardise by
        sqrt(V), band-standardise per octave."""
        Cv = Cfull[Ns]
        tot = np.bincount(inv, weights=Cv, minlength=len(uniq))
        Z = (Cv - (tot / cnt)[inv]) / np.sqrt(V[Ns])
        out = Z.copy()
        for _lo, _hi, sel in bands:
            out[sel] = (Z[sel] - Z[sel].mean()) / Z[sel].std(ddof=1)
        return Z, out

    def rho_bands(Cfull):
        """rho = mean(C^2)/mean(V) per octave, with the cell means
        removed -- the centred estimator of #106."""
        Cv = Cfull[Ns]
        tot = np.bincount(inv, weights=Cv, minlength=len(uniq))
        Cc = Cv - (tot / cnt)[inv]
        return [float((Cc[sel] ** 2).mean() / V[Ns][sel].mean())
                for _lo, _hi, sel in bands]

    L = np.log(Ns.astype(np.float64))
    cols = []
    for t in GAMMAS:
        cols.append(np.cos(t * L))
        cols.append(np.sin(t * L))
    B = np.stack(cols, axis=1)
    B -= B.mean(axis=0)
    Ginv = np.linalg.inv(B.T @ B)

    def r2_of(z):
        zc = z - z.mean()
        bty = B.T @ zc
        return float((Ginv @ bty) @ bty / (zc @ zc))

    C_real = wall(mu)
    _, Zr = pipeline(C_real)
    r2_real = r2_of(Zr)
    rho_real = rho_bands(C_real)

    rng = np.random.default_rng(2026)
    idx = np.nonzero(supp)[0]
    NDRAW = 20
    coin_r2, coin_rho = [], []
    for _ in range(NDRAW):
        eps = np.zeros(X + 1)
        eps[idx] = rng.choice([-1.0, 1.0], size=len(idx))
        Cc = wall(eps)
        _, Zc = pipeline(Cc)
        coin_r2.append(r2_of(Zc))
        coin_rho.append(rho_bands(Cc))
    coin_r2 = np.array(coin_r2)
    coin_rho = np.array(coin_rho)

    print("The coin control of Lemma 18, applied to two figures the")
    print("correction record withdrew.")
    print()
    print("(1) rho = mean(C^2)/mean(V) per octave, cell means removed")
    hdr = (f"    {'band':>22} {'real rho':>9} {'coin mean':>10} "
           f"{'coin sd':>8} {'z of real':>10}")
    print(hdr)
    print("    " + "-" * (len(hdr) - 4))
    zs = []
    for i, (lo_, hi_, _s) in enumerate(bands):
        m, s = coin_rho[:, i].mean(), coin_rho[:, i].std()
        z = (rho_real[i] - m) / s
        zs.append(z)
        print(f"    {f'{lo_}-{hi_}':>22} {rho_real[i]:>9.4f} "
              f"{m:>10.4f} {s:>8.4f} {z:>10.2f}")
    print(f"    v1's #106 coin curve: 0.761 0.780 0.787 0.792 0.813 "
          f"0.828 0.851 0.860")
    print(f"    |z| of the real value against the coin: max "
          f"{max(abs(z) for z in zs):.2f}")
    ok1 = max(abs(z) for z in zs) < 3.0
    print(f"    (1) a coin reproduces the measured rho: "
          f"{'YES -- the estimator cannot tell mu from a coin' if ok1 else 'NO'}")
    print()

    print("(2) the zeta-ordinate regression, same pipeline")
    print(f"    real R^2            {r2_real:.6e}")
    print(f"    coin mean           {coin_r2.mean():.6e}")
    print(f"    coin max of {NDRAW}      {coin_r2.max():.6e}")
    print(f"    coin draws >= real  {int((coin_r2 >= r2_real).sum())}"
          f" of {NDRAW}")
    print(f"    ratio real/coin mean {r2_real/coin_r2.mean():.2f}x")
    print(f"    v1's #110: real 3.896e-3, coin mean 2.994e-3, coin max")
    print(f"    5.515e-3, 6 of 20 coin draws at or above real, "
          f"ratio 1.30x")
    print(f"    the paper quotes a 200-surrogate maximum of 5.09e-6, "
          f"i.e. a")
    print(f"    ratio near 770x, against a null that is not a wall.")
    print()

    # (3) the local-background test, on real and on one coin
    def amp(z, f):
        c = np.cos(f * L)
        s = np.sin(f * L)
        c -= c.mean()
        s -= s.mean()
        zc = z - z.mean()
        G = np.array([[c @ c, c @ s], [c @ s, s @ s]])
        co = np.linalg.solve(G, np.array([c @ zc, s @ zc]))
        return math.hypot(co[0], co[1])

    print("(3) the local-background test of round 5, re-run on a coin")
    eps = np.zeros(X + 1)
    eps[idx] = rng.choice([-1.0, 1.0], size=len(idx))
    _, Zcoin = pipeline(wall(eps))
    rng2 = np.random.default_rng(4242)
    print(f"    {'gamma':>10} {'real above p99':>15} {'coin above p99':>15}")
    hits_r = hits_c = 0
    for t in GAMMAS:
        ctr = []
        while len(ctr) < 60:
            x = float(rng2.uniform(t - 4.0, t + 4.0))
            if min(abs(x - u) for u in GAMMAS) < 1.0:
                continue
            ctr.append(x)
        pr = float(np.quantile([amp(Zr, x) for x in ctr], 0.99))
        pc = float(np.quantile([amp(Zcoin, x) for x in ctr], 0.99))
        ar, ac = amp(Zr, t), amp(Zcoin, t)
        hits_r += ar > pr
        hits_c += ac > pc
        print(f"    {t:>10.4f} {'YES' if ar > pr else 'no':>15} "
              f"{'YES' if ac > pc else 'no':>15}")
    print(f"    real: {hits_r}/10 above their local 99th percentile")
    print(f"    coin: {hits_c}/10 above their local 99th percentile")
    print()
    print("    If the coin clears comparably, then round 5's finding")
    print("    'the effect is real' is a statement about Lambda, which")
    print("    carries the zeros by the explicit formula, and NOT about")
    print("    the wall's mu-fluctuation. That is #110's verdict.")
    print("DONE")


if __name__ == "__main__":
    main()
