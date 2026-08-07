# -*- coding: utf-8 -*-
"""
Conjecture L, the half that matters (increment 237): the TAIL and the
EXTREMES of the wall's own fluctuation.

Conjecture L asserts that every mu-family factorises as a deterministic
local mask times an exactly Gaussian fluctuation. The mask half is
blind-verified (corr 1.0000). The Gaussian half has been checked at
pair, cell, matrix and E1-ratio level -- all BULK statistics. But
EH_mu, the hypothesis the whole chain consumes, is a MAXIMUM:

    Sum_q max_y max_a | ... |  <<_A  N (log N)^{-A}.

A max statement is decided by the tail, not the bulk. A distribution can
match a Gaussian to four moments and still have a max that is off by a
constant factor -- and a constant factor in the max, summed over q, is
exactly the kind of loss that decides whether a chain closes. So the
tail is the half of Conjecture L that the chain actually needs, and it
has not been tested.

THE OBJECT. The wall itself, over its whole family:

    C(N) = Sum_{v+n=N} mu(v) Lambda(n),   G(N) := C(N) / sqrt(S(N) N),

computed for EVERY even N <= X at once by additive convolution (FFT),
which is what makes a tail test possible at all -- 2 x 10^6 samples
rather than the handful a direct sum would allow.

NULLS, ALL COMPUTED BEFORE THE THRESHOLDS.
  * scale: Conjecture L's law says sd(G) = 1. Reported. The SHAPE tests
    below are then run on G/sd(G), so that a wrong constant in the
    normalisation cannot be mistaken for a wrong tail.
  * tails: P(|G| > t) against 2*Phi(-t) for t = 1,2,3,4,5, with the
    binomial standard error sqrt(p(1-p)/n) on the same line.
  * extremes: for n iid N(0,1),
        E[max] = b_n + gamma/a_n,  sd[max] = pi/(sqrt(6) a_n),
        a_n = sqrt(2 log n),
        b_n = a_n - (log log n + log(4 pi)) / (2 a_n).
    Measured max|G| is reported as a z-score against that.
  * top-k: the k-th largest of n iid N(0,1) has expectation
    approximately Phi^{-1}(1 - (k-0.5)/n); reported for k = 1..10.

CRITERION.
  GAUSSIAN TAIL iff every tail probability lies within 3 binomial SE of
                2*Phi(-t), and |z| < 3 for the max.
  HEAVIER       iff tails exceed the Gaussian systematically -- which
                would mean the max-form of EH_mu costs more than the
                bulk statistics suggest.
  LIGHTER       iff they fall short systematically -- which would be
                good news for the chain and bad news for Conjecture L
                as stated, since it claims EXACTLY Gaussian.

WHY EITHER DEVIATION MATTERS. Conjecture L says "exactly Gaussian". A
lighter tail refutes it as stated while helping the chain; a heavier
tail refutes it and hurts. Only agreement leaves it standing.
"""
import numpy as np
import math
import time
from math import erf, sqrt, log, pi


def Phi(t):
    return 0.5 * (1.0 + erf(t / sqrt(2.0)))


def Phi_inv(p):
    # Acklam-style rational approximation, adequate at these precisions
    a = [-3.969683028665376e+01, 2.209460984245205e+02,
         -2.759285104469687e+02, 1.383577518672690e+02,
         -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02,
         -1.556989798598866e+02, 6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01,
         2.445134137142996e+00, 3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if p < pl:
        q = sqrt(-2 * log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > ph:
        q = sqrt(-2 * log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5; r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def main():
    X = 4_000_000
    t0 = time.time()

    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8); mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1)
    for p in primes:
        q = int(p); lp = math.log(int(p))
        while q <= X:
            lam[q] = lp; q *= int(p)
    print(f"sieve  t={time.time()-t0:.0f}s", flush=True)

    # C(N) for ALL N at once, by additive convolution
    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    A = np.zeros(n_fft); A[: X + 1] = mu
    B = np.zeros(n_fft); B[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]
    del A, B
    print(f"convolution  t={time.time()-t0:.0f}s", flush=True)

    # singular series by sieve
    C2 = 0.6601618158468696
    S = np.full(X + 1, 2 * C2)
    for p in primes:
        p = int(p)
        if p > 2:
            S[p::p] *= (p - 1) / (p - 2)

    lo = 100_000
    Ns = np.arange(lo + (lo % 2), X + 1, 2)
    G = C[Ns] / np.sqrt(S[Ns] * Ns)
    n = len(G)
    print(f"samples n = {n}   t={time.time()-t0:.0f}s", flush=True)

    m = float(G.mean()); sd = float(G.std())
    g3 = float(((G - m) ** 3).mean() / sd ** 3)
    g4 = float(((G - m) ** 4).mean() / sd ** 4)
    print(f"\n(A) scale and bulk")
    print(f"  mean {m:+.4f} (null 0, SE {1/math.sqrt(n):.4f})")
    print(f"  sd   {sd:.4f} (Conjecture L's law says 1)")
    print(f"  skew {g3:+.4f} (null 0, SE {math.sqrt(6/n):.4f})")
    print(f"  kurt {g4:.4f} (null 3, SE {math.sqrt(24/n):.4f})")

    # ---------------------------------------------------------------
    # BEFORE any tail claim: is the normalisation N-independent?
    # A scale MIXTURE is leptokurtic and heavy-tailed automatically, so
    # pooling N over a wide range would manufacture exactly the
    # deviation this script is looking for. Hazard 1, scale-
    # normalisation drift (CLOSURE_REAUDIT). Split by dyadic band.
    print(f"\n(A2) is sd(G) constant in N? (if not, pooling mixes")
    print(f"     scales and manufactures a heavy tail by itself)")
    print(f"{'N band':>22} {'count':>9} {'mean':>9} {'sd':>9} "
          f"{'kurt':>8}")
    bands = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if sel.sum() > 1000:
            g = G[sel]
            mm = float(g.mean()); ss = float(g.std())
            kk = float(((g - mm) ** 4).mean() / ss ** 4)
            bands.append((b, hi, int(sel.sum()), mm, ss, kk, sel))
            print(f"{b:>10}-{hi:>11} {int(sel.sum()):>9} {mm:>9.4f} "
                  f"{ss:>9.4f} {kk:>8.4f}")
        b *= 2
    sds = [r[4] for r in bands]
    drift = max(sds) / min(sds)
    print(f"     sd range [{min(sds):.4f}, {max(sds):.4f}], "
          f"ratio {drift:.4f}")
    print(f"     within-band kurtosis is the honest shape statistic;")
    print(f"     the pooled one is contaminated whenever ratio > 1")

    # band-standardised sample: each band by its OWN sd
    Zb = np.empty_like(G)
    for (b, hi, cnt, mm, ss, kk, sel) in bands:
        Zb[sel] = (G[sel] - mm) / ss
    covered = np.zeros(len(G), bool)
    for r in bands:
        covered |= r[6]
    Zb = Zb[covered]
    nb = len(Zb)
    print(f"     band-standardised sample: n = {nb}")

    Z = (G - m) / sd          # shape only; scale already reported
    print(f"\n(B) tails of G/sd -- POOLED (see (A2) before reading)")
    print(f"{'t':>4} {'P(|Z|>t) meas':>15} {'2Phi(-t)':>11} "
          f"{'binom SE':>10} {'z':>8}")
    tail_ok = True
    for t in (1, 2, 3, 4, 5):
        obs = float((np.abs(Z) > t).mean())
        exp = 2 * Phi(-t)
        se = math.sqrt(exp * (1 - exp) / n)
        z = (obs - exp) / se
        if abs(z) >= 3:
            tail_ok = False
        print(f"{t:>4} {obs:>15.6f} {exp:>11.6f} {se:>10.6f} {z:>8.2f}")

    print(f"\n(B2) the same tails, BAND-STANDARDISED -- the honest one")
    print(f"{'t':>4} {'P(|Zb|>t) meas':>16} {'2Phi(-t)':>11} "
          f"{'binom SE':>10} {'z':>8}")
    tail_ok_b = True
    for t in (1, 2, 3, 4, 5):
        obs = float((np.abs(Zb) > t).mean())
        exp = 2 * Phi(-t)
        se = math.sqrt(exp * (1 - exp) / nb)
        z = (obs - exp) / se
        if abs(z) >= 3:
            tail_ok_b = False
        print(f"{t:>4} {obs:>16.6f} {exp:>11.6f} {se:>10.6f} {z:>8.2f}")

    n = nb
    Z = Zb
    a_n = math.sqrt(2 * math.log(n))
    b_n = a_n - (math.log(math.log(n)) + math.log(4 * pi)) / (2 * a_n)
    gam = 0.5772156649
    e_max = b_n + gam / a_n
    sd_max = pi / (math.sqrt(6) * a_n)
    obs_max = float(np.abs(Z).max())
    zmax = (obs_max - e_max) / sd_max
    print(f"\n(C) the extreme -- max|Z| against the Gaussian-max law")
    print(f"  observed max|Z| = {obs_max:.4f}")
    print(f"  E[max] = {e_max:.4f}  sd[max] = {sd_max:.4f}"
          f"   z = {zmax:+.2f}")
    print(f"  (sqrt(2 log n) = {a_n:.4f})")

    print(f"\n(D) the top ten, against their order statistics")
    top = np.sort(np.abs(Z))[-10:][::-1]
    print(f"{'k':>3} {'|Z|_(k)':>10} {'expected':>10} {'diff':>8}")
    for k in range(1, 11):
        e = Phi_inv(1 - (k - 0.5) / (2 * n))    # |Z|: two-sided
        print(f"{k:>3} {top[k-1]:>10.4f} {e:>10.4f} "
              f"{top[k-1]-e:>+8.4f}")

    if drift > 1.02:
        head = (f"NORMALISATION DRIFTS (sd ratio {drift:.3f}) -- the "
                f"pooled tail in (B) is contaminated; read (B2)")
    else:
        head = f"normalisation stable (sd ratio {drift:.3f})"
    print(f"\nverdict: {head}")
    print("        ",
          "GAUSSIAN TAIL -- band-standardised tails within 3 SE "
          "and |z_max| < 3"
          if tail_ok_b and abs(zmax) < 3 else
          "DEVIATION SURVIVES BAND-STANDARDISATION -- see (B2)")
    print("what this decides: EH_mu is a max statement, so the tail is")
    print("the half of Conjecture L the chain actually consumes.")
    print("DONE")


if __name__ == "__main__":
    main()
