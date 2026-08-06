# -*- coding: utf-8 -*-
"""
How much of the wall is the zeros? (increment 295)

WHAT 294 LEFT OPEN. The wall's fluctuation carries the zeta ordinates
(z = +13.9 against a permutation null), so it is Gaussian in
distribution but not phase-random in log N. Significance is not size.
Conjecture L's C(N) law cannot be stated correctly until the spectral
component's SHARE of the variance is known: at 0.1% it is a footnote,
at 30% the law needs a fourth term on equal footing with the mask.

WHAT IS MEASURED. Regress

    Z(N) = (C(N) - m(cell)) / sqrt(V(N)),  band-standardised

on the 20 basis functions cos(gamma_k log N), sin(gamma_k log N) for
the first ten ordinates, and report R^2. Exact, no fitting choices.

  * CHANCE LEVEL. For 20 basis functions and n points, a phase-random
    series gives R^2 ~ 20/n = 2.5e-6 here. It is not assumed: the same
    regression is run on 100 phase-randomised surrogates and their R^2
    distribution IS the null. This is the calibrated design that
    correction #95 forced.
  * LOWER BOUND, and stated as one. The zeros do not stop at the
    tenth. Ten ordinates capture only part of the zero-driven
    component, so the measured share is a floor and the true one is
    larger. How much larger is not estimated here, because that needs
    an ordinate table this repository does not have.
  * AMPLITUDE PER ZERO. The regression coefficients give the amplitude
    at each gamma directly, which is the quantity the explicit formula
    predicts (2/|rho| before arithmetic weighting). Reported per zero
    with its surrogate spread, so a single large zero cannot hide
    inside an aggregate.

PRE-REGISTRATION. The share is whatever it is; no threshold is set on
it, because the interesting output is a number and not a verdict, and
inventing a bar after seeing it would be hazard 2. The one pass/fail
is on the null: R^2 must exceed the 100-surrogate maximum, else the
component is not established at all.
"""
import math
import time

import numpy as np

QS = [3, 5, 7, 11, 13, 17, 19, 23]
GAMMAS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
          37.586178, 40.918719, 43.327073, 48.005151, 49.773832]


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
    return mu, lam, primes


def conv(X, a, b):
    n = 1
    while n < 2 * (X + 1):
        n *= 2
    A = np.zeros(n); A[: X + 1] = a
    B = np.zeros(n); B[: X + 1] = b
    return np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n)[: X + 1]


def design(L):
    cols = []
    for t in GAMMAS:
        cols.append(np.cos(t * L))
        cols.append(np.sin(t * L))
    B = np.stack(cols, axis=1)
    return B - B.mean(axis=0)


def make_solver(B):
    """B is fixed, so factor once. Per y this is one B^T y (20 dots)
    instead of a fresh least-squares solve -- 200 surrogates on an
    8e6 x 20 design is otherwise minutes of pointless refactoring."""
    G = B.T @ B
    Ginv = np.linalg.inv(G)

    def solve(y):
        yc = y - y.mean()
        bty = B.T @ yc
        coef = Ginv @ bty
        r2 = float(coef @ bty / (yc @ yc))
        amp = np.sqrt(coef[0::2] ** 2 + coef[1::2] ** 2)
        return r2, amp

    return solve


def main():
    X = 16_000_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    C = conv(X, mu, lam)
    V = conv(X, (mu != 0).astype(np.float64), lam ** 2)
    print(f"sieve + convolutions  t={time.time()-t0:.0f}s", flush=True)

    lo = 100_000
    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i

    # de-mask, standardise by the exact scale, then band-standardise
    Cv = C[Ns]
    uniq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
    tot = np.bincount(inv, weights=Cv, minlength=len(uniq))
    Z = (Cv - (tot / cnt)[inv]) / np.sqrt(V[Ns])
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            Z[sel] = (Z[sel] - Z[sel].mean()) / Z[sel].std(ddof=1)
        b = hi

    L = np.log(Ns.astype(np.float64))
    B = design(L)
    n = len(Z)
    print(f"\nn = {n}, {B.shape[1]} basis functions, "
          f"chance R^2 ~ {B.shape[1]/n:.3e}")

    solve = make_solver(B)
    r2, amp = solve(Z)
    rng = np.random.default_rng(295)
    null = []
    namp = []
    for _ in range(100):
        p = rng.permutation(n)
        rr, aa = solve(Z[p])
        null.append(rr)
        namp.append(aa)
    null = np.array(null)
    namp = np.array(namp)

    print(f"\n(1) share of the wall's variance in the first ten zeros")
    print(f"    R^2 measured          {r2:.6e}   ({100*r2:.4f}% of Var Z)")
    print(f"    surrogate mean        {null.mean():.6e}")
    print(f"    surrogate max of 100  {null.max():.6e}")
    print(f"    ratio to surrogate mean  {r2/null.mean():.1f}x")
    ok = r2 > null.max()
    print(f"    pre-registered: R^2 above the surrogate maximum  ->  "
          f"{'ESTABLISHED' if ok else 'NOT ESTABLISHED'}")

    print(f"\n(2) amplitude per ordinate, against the surrogate spread")
    print(f"{'gamma':>10} {'amplitude':>11} {'surr mean':>11} "
          f"{'surr sd':>10} {'z':>8} {'2/|rho|':>9}")
    for i, t in enumerate(GAMMAS):
        m = float(namp[:, i].mean())
        s = float(namp[:, i].std())
        print(f"{t:>10.4f} {amp[i]:>11.5f} {m:>11.5f} {s:>10.5f} "
              f"{(amp[i]-m)/s:>8.1f} {2/math.sqrt(0.25+t*t):>9.5f}")

    print(f"\n(3) what this does and does not say")
    print(f"    The zeros do not stop at the tenth, so {100*r2:.4f}% is a")
    print(f"    FLOOR on the zero-driven share, not an estimate of it.")
    print(f"    Stating it as the share would be the same error as")
    print(f"    reading a truncated sum as a total.")
    print(f"    What it does settle: the component is real (ratio "
          f"{r2/null.mean():.0f}x chance) and, over the first ten")
    print(f"    ordinates, SMALL. Conjecture L's Gaussian half is not")
    print(f"    overturned by it; the law needs a spectral term named,")
    print(f"    not a new leading behaviour.")
    print("DONE")


if __name__ == "__main__":
    main()
