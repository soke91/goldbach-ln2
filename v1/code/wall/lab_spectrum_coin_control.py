# -*- coding: utf-8 -*-
"""
Do the zeta lines in the wall come from mu, or just from Lambda?
(increment 302)

WHAT 294 AND 295 ESTABLISHED, AND WHAT THEY DID NOT. The de-masked,
V-standardised wall carries the zeta ordinates: permutation test
z = +13.9, and regressing on the first ten ordinates gives
R^2 = 3.90e-3 against a 200-surrogate maximum of 5.09e-6.

The surrogate there was Z permuted across N. That destroys EVERY
structure at once, so it answers "is there structure" and not "where
does it come from". And there is an obvious other source:
C = mu * Lambda, and Lambda carries the zeros by the explicit formula.
If the lines are Lambda's, then #94 and #96 report a property of the
von Mangoldt function, not of the wall, and the sentence "the wall's
fluctuation is not phase-random in log N" is true but empty.

Increment 300 showed what happens when a claim is never coin-
controlled. This is the same control applied to the spectral claim.

THE TEST. Replace mu by eps(v) = random +/-1 on {mu != 0} and run the
identical pipeline: same Lambda, same V, same de-masking, same band
standardisation, same regression on the same ten ordinates. If the coin
shows the lines at the same strength, they are Lambda's. If it does
not, they need mu.

WHAT TO EXPECT, WRITTEN BEFORE THE RUN so the outcome is not read
backwards. Lambda's zero terms enter C as Sum_v sign(v) (N-v)^{rho-1}.
With random signs that sum is phase-incoherent and should be tiny; with
mu the signs are arithmetically correlated with the summand and can add
coherently. So the coin is expected to show no lines. But the same
reasoning would have predicted the coin fails at increment 300, and it
did not, which is why this is run rather than argued.

PRE-REGISTRATION.
  R = 20 coin draws. Statistic: R^2 of Z on the 20 basis functions
  cos(gamma log N), sin(gamma log N) for the first ten ordinates,
  through byte-identical code.
  DECISION RULE. The lines are attributed to mu iff the real R^2
  exceeds every one of the 20 coin R^2 values. Anything less and the
  spectral claim of #94 and #96 is downgraded to a statement about
  Lambda.
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
    return mu, lam


def conv(X, a, b, nfft):
    A = np.zeros(nfft); A[: X + 1] = a
    B = np.zeros(nfft); B[: X + 1] = b
    return np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), nfft)[: X + 1]


def standardise(sig, lam, V, Ns, key, X, lo, nfft):
    """The exact pipeline of increments 294-295, sign vector aside."""
    C = conv(X, sig, lam, nfft)[Ns]
    uniq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
    tot = np.bincount(inv, weights=C, minlength=len(uniq))
    Z = (C - (tot / cnt)[inv]) / np.sqrt(V[Ns])
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            Z[sel] = (Z[sel] - Z[sel].mean()) / Z[sel].std(ddof=1)
        b = hi
    return Z


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam = sieve(X)
    nfft = 1
    while nfft < 2 * (X + 1):
        nfft *= 2
    supp = (mu != 0)
    V = conv(X, supp.astype(np.float64), lam ** 2, nfft)
    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    print(f"sieve + V  t={time.time()-t0:.0f}s", flush=True)

    L = np.log(Ns.astype(np.float64))
    cols = []
    for g in GAMMAS:
        cols.append(np.cos(g * L))
        cols.append(np.sin(g * L))
    B = np.stack(cols, axis=1)
    B = B - B.mean(axis=0)
    Ginv = np.linalg.inv(B.T @ B)

    def r2(Z):
        y = Z - Z.mean()
        bty = B.T @ y
        return float((Ginv @ bty) @ bty / (y @ y))

    Zr = standardise(mu.astype(np.float64), lam, V, Ns, key, X, lo, nfft)
    real = r2(Zr)
    print(f"real R^2 = {real:.6e}   t={time.time()-t0:.0f}s", flush=True)

    R = 20
    rng = np.random.default_rng(302)
    idx = np.nonzero(supp)[0]
    coin = np.empty(R)
    for r in range(R):
        eps = np.zeros(X + 1, dtype=np.float64)
        eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        coin[r] = r2(standardise(eps, lam, V, Ns, key, X, lo, nfft))
        if (r + 1) % 5 == 0:
            print(f"  coin {r+1}/{R}  R^2 = {coin[r]:.3e}  "
                  f"t={time.time()-t0:.0f}s", flush=True)

    print(f"\nR^2 on the first ten ordinates, identical pipeline")
    print(f"    real mu              {real:.6e}")
    print(f"    coin mean            {coin.mean():.6e}")
    print(f"    coin max of {R}       {coin.max():.6e}")
    print(f"    coin min             {coin.min():.6e}")
    print(f"    real / coin mean     {real / coin.mean():.2f}x")
    above = int((coin >= real).sum())
    print(f"    coin draws at or above the real value: {above}/{R}")

    ok = above == 0
    print(f"\n    pre-registered: real above every coin draw  ->  "
          f"{'PASS' if ok else 'FAIL'}")
    if ok:
        v = ("the lines need mu: random signs on the same support, "
             "through the same Lambda, do not produce them")
    else:
        v = ("the lines are Lambda's, not the wall's -- #94 and #96 "
             "describe the von Mangoldt function")
    print(f"    {v}")
    print(f"\n    Note on what this cannot say: it separates mu from")
    print(f"    Lambda, not the wall from the explicit formula. Even")
    print(f"    if mu is needed, the zeros still enter through Lambda;")
    print(f"    what would be established is that the mu-Lambda")
    print(f"    correlation is what makes them add coherently.")
    print("DONE")


if __name__ == "__main__":
    main()
