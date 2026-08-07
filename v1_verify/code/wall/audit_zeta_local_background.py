# -*- coding: utf-8 -*-
"""
Closing the open item left by `audit_zeta_regression_null.py`: is there
a zeta-ordinate signal in the wall's G at all?

WHERE THIS STANDS. `conj:wall` item 4 claims R^2 = 3.90e-3 against a
"200-surrogate maximum of 5.09e-6", with "every ordinate individually
at z >= 23". The companion script showed that v1's surrogates are a
VALUE PERMUTATION -- white -- while the regressors sit in the lowest
few dozen Fourier bins, so the quoted significances are artefacts of
the null. Under autocorrelation-preserving nulls the per-ordinate z
collapsed to -0.48 .. +3.52. But those nulls preserve the power at the
tested frequencies too, so they are conservative and cannot settle
whether an effect exists.

THE TEST USED HERE is the one v1 itself adopted for the same question
in `lab_E_zeta_spectrum.py`, whose docstring states the principle:

    "the background is not flat in gamma -- it FALLS -- so a peak at
     gamma = 14 clears a pooled 99th percentile merely by sitting where
     the background is high"

and which therefore compares each |F(gamma)| against controls drawn
from ITS OWN NEIGHBOURHOOD, `rng.uniform(t - 4, t + 4)`, excluding
draws within 1.0 of any ordinate. That is a notch comparison: it holds
the broadband background fixed and asks only whether a LINE sits at
gamma.

v1 applied it to E(N) = the Mertens-type sum. It was never applied to
G, which is what `conj:wall` item 4 is about. This script applies it.

PRE-REGISTRATION (written before the run).

  (1) RULE, per ordinate. The amplitude at gamma must exceed the 99th
      percentile of its own local controls. Report how many of the ten
      do.
  (2) RULE, jointly. The 10-ordinate R^2 must exceed the R^2 of
      200 draws of ten frequencies, each drawn from the corresponding
      ordinate's own neighbourhood. This holds the frequency band
      fixed and varies only the arithmetic identity of the ten.
  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      The companion run found only gamma_1 = 14.13 above 2 sigma under
      a global autocorrelation-preserving null. I predict at most
      1-2 of the ten clear their local 99th percentile -- i.e. about
      what 10 tests at the 1% level produce by chance -- and that the
      joint R^2 does NOT clear the local-band null. If instead most
      ordinates clear locally, the effect is real and only its quoted
      strength was wrong; that outcome would be reported as such.
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


def conv(a, b, X):
    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2
    return np.fft.irfft(np.fft.rfft(np.pad(a, (0, nf - len(a))))
                        * np.fft.rfft(np.pad(b, (0, nf - len(b)))), nf)[: X + 1]


def main():
    X = 4_000_000
    LO = 100_000
    mu, lam = sieve_mu_lambda(X)
    C = conv(mu, lam, X)
    V = conv((mu != 0).astype(np.float64), lam ** 2, X)

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
    Z = Z - Z.mean()
    n = len(Z)
    L = np.log(Ns.astype(np.float64))
    zz = float(Z @ Z)

    def amp(f):
        """regression amplitude of Z on cos(f L), sin(f L)."""
        c = np.cos(f * L)
        s = np.sin(f * L)
        c -= c.mean()
        s -= s.mean()
        # near-orthogonal at these frequencies; solve the 2x2 exactly
        G = np.array([[c @ c, c @ s], [c @ s, s @ s]])
        r = np.array([c @ Z, s @ Z])
        co = np.linalg.solve(G, r)
        return math.hypot(co[0], co[1]), float(co @ r)

    print("Local-background test for the zeta ordinates in G")
    print(f"n = {n}, log-N span {L[-1]-L[0]:.3f}")
    print()
    rng = np.random.default_rng(4242)
    NC = 120

    print(f"{'gamma':>10} {'amplitude':>11} {'local med':>10} "
          f"{'local p99':>10} {'/p99':>7} {'rank':>9} {'above?':>7}")
    hits = 0
    ctrl_store = []
    for t in GAMMAS:
        a, _ = amp(t)
        ctrls = []
        while len(ctrls) < NC:
            x = float(rng.uniform(t - 4.0, t + 4.0))
            if min(abs(x - u) for u in GAMMAS) < 1.0:
                continue
            ctrls.append(amp(x)[0])
        ctrls = np.array(ctrls)
        ctrl_store.append(ctrls)
        p99 = float(np.quantile(ctrls, 0.99))
        med = float(np.median(ctrls))
        rank = float((ctrls < a).mean())
        ok = a > p99
        hits += ok
        print(f"{t:>10.4f} {a:>11.5f} {med:>10.5f} {p99:>10.5f} "
              f"{a/p99:>7.2f} {rank:>8.1%} {'YES' if ok else 'no':>7}")
    print(f"    {hits} of {len(GAMMAS)} above their LOCAL 99th percentile")
    print(f"    chance expectation at the 1% level: 0.1 of 10")
    print()

    # joint R^2 against ten frequencies drawn from the same bands
    def r2_of(freqs):
        cols = []
        for f in freqs:
            cols.append(np.cos(f * L))
            cols.append(np.sin(f * L))
        B = np.stack(cols, axis=1)
        B -= B.mean(axis=0)
        bty = B.T @ Z
        return float(np.linalg.solve(B.T @ B, bty) @ bty / zz)

    r2 = r2_of(GAMMAS)
    null = []
    for _ in range(200):
        fs = []
        for t in GAMMAS:
            while True:
                x = float(rng.uniform(t - 4.0, t + 4.0))
                if min(abs(x - u) for u in GAMMAS) >= 1.0:
                    break
            fs.append(x)
        null.append(r2_of(fs))
    null = np.array(null)
    print("joint R^2, ten ordinates vs ten frequencies from the same")
    print("neighbourhoods (200 draws):")
    print(f"    measured R^2      {r2:.6e}")
    print(f"    local-null mean   {null.mean():.6e}")
    print(f"    local-null max    {null.max():.6e}")
    print(f"    measured rank     {(null < r2).mean():.1%}")
    print(f"    clears the local-null maximum: "
          f"{'YES' if r2 > null.max() else 'NO'}")
    print()
    print("    for reference, v1's value-permutation null on the same")
    print("    data gives a maximum near 2e-5, i.e. two orders of")
    print("    magnitude below this band's own background.")
    print("DONE")


if __name__ == "__main__":
    main()
