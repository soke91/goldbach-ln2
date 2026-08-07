# -*- coding: utf-8 -*-
"""
Re-verification of measurement 4 of Conjecture 14 (`conj:wall`) in
v1/paper/wall_v1.tex: the claim that G carries a zeta-zero signature.

THE STATEMENT UNDER TEST, verbatim:

    "G is Gaussian in distribution but not phase-random in log N.
     Regressing G on cos(gamma log N), sin(gamma log N) for the first
     ten ordinates of the zeta zeros gives R^2 = 3.90e-3 against a
     200-surrogate maximum of 5.09e-6, and every ordinate individually
     at z >= 23."

THE OBJECTION BEING TESTED -- the null. The verdict is "R^2 above the
surrogate maximum", so everything rests on what the surrogates are.
C(N) for nearby N share most of their terms mu(v)Lambda(N-v); the paper
says so itself, in the discussion of Lemma `lem:coin`: "nearby N share
the same mu(v)Lambda(N-v) terms and their C(N) are positively
correlated". So Z is a SERIALLY CORRELATED (red) series.

The regression frequencies are very low. Over 1e5 <= N <= 1.6e7 the
log-N span is 5.08, so the fundamental is 2*pi/5.08 = 1.237 and the
first ordinate gamma = 14.13 is the ELEVENTH Fourier bin. Essentially
all of a red series' power sits in bins that low.

A null that destroys the serial correlation therefore has far less
low-frequency power than the data, and any red series whatever will
clear it -- zeta zeros or no zeta zeros. The right null preserves the
autocorrelation and randomises only the phases.

METHOD HERE. Z is rebuilt from the paper's description (mask removed by
finite modular enumeration over 3,5,7,11,13; divided by sqrt(V);
standardised per octave) and the same 20-column design is regressed.
Three nulls are then compared:

  NULL P  value permutation of Z across N. This is what v1's
          lab_wall_spectral_share.py actually does (`rng.permutation`),
          although its own docstring calls the surrogates
          "phase-randomised".
  NULL F  phase randomisation: keep |FFT(Z)| exactly, replace every
          phase by a uniform draw. This preserves the autocorrelation
          of Z exactly and destroys any phase alignment -- it is the
          null the docstring describes.
  NULL C  circular rotation of Z by a random offset. Preserves the
          autocorrelation exactly and by construction contains no
          alignment between Z's features and log N.

PRE-REGISTRATION (written before the run).

  (1) RULE. The claim stands if R^2 exceeds the maximum over the
      autocorrelation-preserving nulls F and C, not merely over P.
  (2) PREDICTION, recorded so it cannot be reported as a surprise.
      I predict null P reproduces v1's ~5e-6 (it is white noise, so
      its R^2 is about 20/n), while nulls F and C sit orders of
      magnitude higher, and that the measured R^2 does NOT clear them.
      If instead R^2 clears F and C too, the finding is void and the
      paper's measurement is stronger than its own null design.
"""
import sys
import math
import time

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


def conv(a, b, X):
    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2
    return np.fft.irfft(np.fft.rfft(np.pad(a, (0, nf - len(a))))
                        * np.fft.rfft(np.pad(b, (0, nf - len(b)))),
                        nf)[: X + 1]


def main():
    X = 4_000_000
    LO = 100_000
    t0 = time.time()
    mu, lam = sieve_mu_lambda(X)
    C = conv(mu, lam, X)
    V = conv((mu != 0).astype(np.float64), lam ** 2, X)
    print(f"sieve + convolutions  t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(LO + LO % 2, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i

    Cv = C[Ns]
    uniq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
    tot = np.bincount(inv, weights=Cv, minlength=len(uniq))
    Z = (Cv - (tot / cnt)[inv]) / np.sqrt(V[Ns])
    b = LO
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            Z[sel] = (Z[sel] - Z[sel].mean()) / Z[sel].std(ddof=1)
        b = hi

    n = len(Z)
    L = np.log(Ns.astype(np.float64))
    cols = []
    for t in GAMMAS:
        cols.append(np.cos(t * L))
        cols.append(np.sin(t * L))
    B = np.stack(cols, axis=1)
    B = B - B.mean(axis=0)
    G = B.T @ B
    Ginv = np.linalg.inv(G)

    def fit(y):
        yc = y - y.mean()
        bty = B.T @ yc
        coef = Ginv @ bty
        r2 = float(coef @ bty / (yc @ yc))
        amp = np.sqrt(coef[0::2] ** 2 + coef[1::2] ** 2)
        return r2, amp

    def r2_of(y):
        return fit(y)[0]

    span = L[-1] - L[0]
    print()
    print(f"n = {n}   log-N span = {span:.3f}   fundamental "
          f"= {2*math.pi/span:.4f}")
    print(f"gamma_1 = {GAMMAS[0]:.4f} is Fourier bin "
          f"{GAMMAS[0]/(2*math.pi/span):.1f} of {n//2}")
    print(f"white-noise chance R^2 = 20/n = {20/n:.3e}")
    # how red is Z?
    for lag in (1, 2, 5, 20, 100):
        r = float(np.corrcoef(Z[:-lag], Z[lag:])[0, 1])
        print(f"    lag-{lag:<4d} autocorrelation of Z = {r:+.4f}")

    r2, amp = fit(Z)
    print()
    print(f"measured R^2 = {r2:.6e}   ({100*r2:.4f}% of Var Z)")
    print()

    rng = np.random.default_rng(295)
    NS = 200
    F = np.fft.rfft(Z - Z.mean())
    mag = np.abs(F)

    def draw(kind):
        if kind == "P":
            return Z[rng.permutation(n)]
        if kind == "F":
            ph = rng.uniform(0, 2 * math.pi, len(F))
            ph[0] = 0.0
            return np.fft.irfft(mag * np.exp(1j * ph), n)
        return np.roll(Z, int(rng.integers(1, n)))

    names = {"P": "P value permutation (v1)",
             "F": "F phase randomisation",
             "C": "C circular rotation"}
    r2n, ampn = {}, {}
    for k in ("P", "F", "C"):
        rr, aa = zip(*(fit(draw(k)) for _ in range(NS)))
        r2n[k] = np.array(rr)
        ampn[k] = np.array(aa)

    hdr = (f"{names['P']:>26} {'mean R^2':>12} {'max of 200':>12} "
           f"{'measured/max':>13} {'clears?':>8}")
    print(f"{'null':>26} {'mean R^2':>12} {'max of 200':>12} "
          f"{'measured/max':>13} {'clears?':>8}")
    print("-" * len(hdr))
    for k in ("P", "F", "C"):
        v = r2n[k]
        print(f"{names[k]:>26} {v.mean():>12.4e} {v.max():>12.4e} "
              f"{r2/v.max():>13.2f} "
              f"{'YES' if r2 > v.max() else 'NO':>8}")

    print()
    print('the paper also claims "every ordinate individually at z >= 23".')
    print("per-ordinate amplitude z, against each null's own spread:")
    print(f"{'gamma':>10} {'amplitude':>11} {'z (null P, v1)':>15} "
          f"{'z (null F)':>11} {'z (null C)':>11}")
    zs = {k: [] for k in ("P", "F", "C")}
    for i, t in enumerate(GAMMAS):
        row = [f"{t:>10.4f}", f"{amp[i]:>11.5f}"]
        for k in ("P", "F", "C"):
            m, s = float(ampn[k][:, i].mean()), float(ampn[k][:, i].std())
            z = (amp[i] - m) / s
            zs[k].append(z)
            row.append(f"{z:>{15 if k=='P' else 11}.2f}")
        print(" ".join(row))
    for k in ("P", "F", "C"):
        print(f"    min z over the ten ordinates, null {k}: "
              f"{min(zs[k]):+.2f}")

    print()
    print("(1) RULE was: the claim stands only if R^2 clears the")
    print("    autocorrelation-preserving nulls F and C, not merely the")
    print("    value-permutation null P.")
    print()
    print("    NOTE on conservatism. Nulls F and C preserve the whole")
    print("    N-domain power spectrum, hence also whatever power sits")
    print("    at the tested frequencies; they are therefore if")
    print("    anything too strict a comparison for a real line. They")
    print("    are nevertheless the right order of magnitude for the")
    print("    null, and null P is not: P is white, and the regressors")
    print("    live in the lowest few dozen Fourier bins.")
    print("DONE")


if __name__ == "__main__":
    main()
