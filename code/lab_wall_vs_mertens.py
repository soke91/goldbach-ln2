# -*- coding: utf-8 -*-
"""
Is the wall essentially the Mertens function? (increment 294)

WHERE THIS COMES FROM. Increment 293's discarded control turned out to
be a real signal: C(N) = Sum_v mu(v) Lambda(N-v) contains Lambda, whose
explicit formula carries the zeta zeros. Pushing that one step further,
Lambda has average 1, so the leading behaviour of C is

    C(N) = Sum_{v<N} mu(v) Lambda(N-v)
         ~ Sum_{v<N} mu(v)  +  (zero corrections)
         =  M(N-1)          +  ...,

with M the MERTENS FUNCTION. And the sizes agree: |M(N)| is conjectured
to be N^{1/2+o(1)}, and increment 281 measured |C| ~ N^{0.5457}, i.e.
sqrt(N) times a small log power. If C were essentially M, the whole
"wall" would be the Mertens function wearing a Goldbach costume -- and
five days of this program have never asked.

The heuristic is crude (Lambda's mean is 1 only after averaging, and
mu(v) is not smooth), so this is a question, not a claim.

PRE-REGISTRATION (fixed before the run).

  (A) CORRELATION. corr(C(N), M(N)) over even N, per octave band, both
      raw and with the location mask removed from C. DECISION RULE:
      "C is essentially M" requires |corr| > 0.5 and not falling with
      N. A small correlation refutes the identification; a large and
      rising one supports it.

  (B) THE RATIO C/M, its mean and spread. If C = M + lower order the
      ratio concentrates; if C and M are merely the same SIZE the
      ratio scatters around zero mean.

  (C) A SIZE CHECK THAT IS NOT A CORRELATION. The exponents of |C| and
      |M| fitted on the same bands by the same estimator. Equal
      exponents are necessary but nowhere near sufficient -- stated so
      the size agreement is not mistaken for evidence of identity,
      which is the error increment 293 corrected in itself.

  (D) THE SPECTRUM. Does the de-masked, V-standardised wall carry the
      zeta ordinates, as E does? Same local-null machinery as
      increment 293, same phase-randomised surrogate. Conjecture L
      says the fluctuation is "exactly Gaussian ... no class
      structure"; a line spectrum is structure of a kind the
      distribution tests of increment 283 cannot see, since an
      oscillatory component of small amplitude leaves kurtosis alone.

  A NOTE ON WHAT WOULD BE OVERCLAIMING. Finding lines in C does NOT
  contradict increment 283. Gaussian in distribution and phase-random
  in log N are different statements, and this program has only tested
  the first.
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


def spectrum(Y, Ns, nbins=1 << 19):
    L = np.log(Ns.astype(np.float64))
    span = L[-1] - L[0]
    idx = ((L - L[0]) / span * (nbins - 1)).astype(np.int64)
    acc = np.bincount(idx, weights=Y * (2.0 / Ns), minlength=nbins)
    acc = acc - acc.mean()
    F = np.fft.rfft(acc)
    g = 2.0 * math.pi * np.arange(len(F)) / span
    return g, np.abs(F) / span * 2.0


def line_strength(Y, Ns):
    """Sum of |F| at the ten ordinates, in units of the band median.

    The local-percentile null used at increment 293 is NOT calibrated:
    the spectral resolution here is 2*pi/span = 1.238, so a +/-4
    window holds only about 6.5 INDEPENDENT frequencies, and a 99th
    percentile estimated from six independent values is meaningless --
    interpolating 4000 draws inside them does not add information. The
    surrogate duly scored 4/10 where a calibrated null would give
    about 0.1/10. So the null is replaced by a PERMUTATION TEST: the
    same statistic on many phase-randomised surrogates, which is
    calibrated by construction and needs no percentile model.
    """
    g, A = spectrum(Y, Ns)
    band = (g >= 5.0) & (g <= 70.0)
    gb, Ab = g[band], A[band]
    med = float(np.median(Ab))
    return float(sum(np.interp(t, gb, Ab) for t in GAMMAS)) / med


def perm_test(tag, Y, Ns, rng, nperm=200):
    real = line_strength(Y, Ns)
    null = np.array([line_strength(Y[rng.permutation(len(Y))], Ns)
                     for _ in range(nperm)])
    pv = float((null >= real).sum() + 1) / (nperm + 1)
    z = (real - null.mean()) / null.std()
    print(f"    {tag:<38} S={real:8.2f}  null {null.mean():6.2f}+/-{null.std():5.2f}  z={z:+7.1f}  p={pv:.4f}")
    return real, pv, z


def fit(x, y):
    return float(np.polyfit(x, y, 1)[0])


def main():
    X = 16_000_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    C = conv(X, mu, lam)
    Mert = np.cumsum(mu.astype(np.int64))
    V = conv(X, (mu != 0).astype(np.float64), lam ** 2)
    print(f"sieve + convolutions  t={time.time()-t0:.0f}s", flush=True)

    lo = 100_000
    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    Cv = C[Ns]
    Mv = Mert[Ns - 1].astype(np.float64)

    print("\n(A)+(B) is C essentially M?")
    print(f"{'band':>21} {'corr(C,M)':>10} {'corr demask':>12} "
          f"{'mean C/M':>10} {'sd C/M':>9} {'logNmid':>8}")
    rows = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) < 1000:
            b = hi
            continue
        c, m = Cv[sel], Mv[sel]
        uniq, inv = np.unique(key[sel], return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
        tot = np.bincount(inv, weights=c, minlength=len(uniq))
        cd = c - (tot / cnt)[inv]
        ok = np.abs(m) > 1e-9
        rat = c[ok] / m[ok]
        r1 = float(np.corrcoef(c, m)[0, 1])
        r2 = float(np.corrcoef(cd, m)[0, 1])
        L = math.log(math.sqrt(b * hi))
        rows.append((r1, r2, float(np.abs(c).mean()),
                     float(np.abs(m).mean()), L))
        print(f"{b:>9}-{hi:>11} {r1:>10.4f} {r2:>12.4f} "
              f"{float(np.median(rat)):>10.3f} {float(rat.std()):>9.1f} "
              f"{L:>8.3f}")
        b = hi

    r1s = [x[0] for x in rows]
    print(f"\n    corr(C,M): first {r1s[0]:+.4f}, last {r1s[-1]:+.4f}")
    ok_A = abs(r1s[-1]) > 0.5 and abs(r1s[-1]) >= abs(r1s[0]) - 0.05
    print(f"    pre-registered |corr| > 0.5 and not falling  ->  "
          f"{'C IS essentially M' if ok_A else 'REFUTED'}")

    L = np.array([x[4] for x in rows])
    bC = fit(L, np.log([x[2] for x in rows]))
    bM = fit(L, np.log([x[3] for x in rows]))
    print(f"\n(C) exponents on the same bands, same estimator")
    print(f"    beta_|C| = {bC:.4f}    beta_|M| = {bM:.4f}")
    print("    Equal exponents are necessary, not sufficient. Increment")
    print("    293 had to correct exactly this confusion in itself.")

    print("\n(D) does the wall carry the zeta ordinates?")
    rng = np.random.default_rng(294)
    Vv = np.sqrt(V[Ns])
    uniq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
    tot = np.bincount(inv, weights=Cv, minlength=len(uniq))
    Z = (Cv - (tot / cnt)[inv]) / Vv
    _, pC, zC = perm_test("C de-masked / sqrt(V)", Z, Ns, rng)
    _, pM, zM = perm_test("M(N)/sqrt(N)  (a known zero sum)",
                          Mv / np.sqrt(Ns), Ns, rng)
    # Increment 293's E result rests on the same uncalibrated local
    # null that is being replaced here, so it is re-run with the
    # permutation statistic rather than left standing on a null now
    # known to be wrong.
    C2b = 0.6601618158468696
    Sb = np.full(X + 1, 2 * C2b)
    for pp in primes:
        pp = int(pp)
        if pp > 2:
            Sb[pp::pp] *= (pp - 1) / (pp - 2)
    rr = conv(X, lam, lam)
    Eb = (rr[Ns] - Sb[Ns] * Ns) / np.sqrt(Ns)
    _, pE, zE = perm_test("E(N)/sqrt(N)  [increment 293, recalibrated]",
                          Eb, Ns, rng)
    print("")
    print("    permutation test, 200 surrogates each; the null is the")
    print("    same values with the log-N phase alignment destroyed")
    if pC <= 0.01:
        v = ("the wall's fluctuation carries the zeros -- Gaussian in "
             "distribution but NOT phase-random in log N")
    else:
        v = "no line structure in the wall at this power"
    print(f"    {v}")
    print("    This does not contradict increment 283: Gaussian in")
    print("    distribution and phase-random in log N are different")
    print("    statements, and only the first was ever tested.")
    print("DONE")


if __name__ == "__main__":
    main()
