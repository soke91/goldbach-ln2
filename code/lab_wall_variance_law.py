# -*- coding: utf-8 -*-
"""
The wall's variance law, corrected (increment 238).

This program has carried the law

    C(N) = sqrt( S(N) * N ) * G(N),   G ~ N(0,1)

since the early measurements. Increment 237's tail test computed C(N)
for EVERY even N <= 4*10^6 by additive convolution and found
sd(G) = 2.24, not 1 -- and, worse, sd(G) DRIFTING monotonically upward
with N, from 2.069 to 2.306 across the range. A normalisation whose
residual scale drifts is a normalisation with a missing factor, and a
drifting scale pooled across N manufactures a heavy tail by itself.
This is hazard 1 (scale-normalisation drift) for the third time.

WHAT THIS SCRIPT TESTS. If the true law is

    Var C(N)  =  kappa * S(N) * N * (log N)^alpha,

then sd(G)^2 / (log N)^alpha is constant in N for the right alpha. The
script fits alpha freely and reports the constancy of the ratio for
alpha = 0 (the recorded law), alpha = 1, and the fitted value.

NULLS AND CRITERION, on the same line.
  * alpha = 0 is the recorded law. It predicts sd(G) constant; the
    measured drift is 2.069 -> 2.306, a ratio of 1.114 across a factor
    25 in N, which is 34 binomial-free standard errors of nothing --
    the band sds are each computed from 5*10^4 to 8*10^5 samples, so
    their sampling error is below 0.005.
  * alpha = 1 predicts sd(G)^2 / log N constant.
  * CONFIRMED iff the coefficient of variation of sd^2/(log N)^alpha
    across bands is below 2% for the fitted alpha AND the fitted alpha
    is within 0.15 of 1.
  * The trivial-scale sanity check is printed too: under a random-sign
    model Var C(N) would be Sum_v mu^2(v) Lambda(N-v)^2, computed here
    exactly, and the ratio of the measured variance to it is reported.
    That is the honest normaliser, and kappa is whatever it is.
"""
import numpy as np
import math
import time


def sieve(X):
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
    return mu, lam, spf, primes


def main():
    X = 4_000_000
    t0 = time.time()
    mu, lam, spf, primes = sieve(X)

    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    A = np.zeros(n_fft); A[: X + 1] = mu
    B = np.zeros(n_fft); B[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]
    # the random-sign normaliser, exactly: Sum_v mu^2(v) Lambda(N-v)^2
    A2 = np.zeros(n_fft); A2[: X + 1] = np.abs(mu)
    B2 = np.zeros(n_fft); B2[: X + 1] = lam ** 2
    V = np.fft.irfft(np.fft.rfft(A2) * np.fft.rfft(B2), n_fft)[: X + 1]
    del A, B, A2, B2
    print(f"convolutions  t={time.time()-t0:.0f}s", flush=True)

    C2 = 0.6601618158468696
    S = np.full(X + 1, 2 * C2)
    for p in primes:
        p = int(p)
        if p > 2:
            S[p::p] *= (p - 1) / (p - 2)

    lo = 100_000
    Ns = np.arange(lo, X + 1, 2)
    G = C[Ns] / np.sqrt(S[Ns] * Ns)

    print(f"\n{'band':>22} {'count':>8} {'sd(G)':>8} {'sd^2':>8} "
          f"{'logNmid':>8} {'sd^2/logN':>10} {'VarC/SumLam^2':>14}")
    xs, ys = [], []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if sel.sum() > 1000:
            g = G[sel]
            sd = float(g.std())
            mid = math.sqrt(b * hi); L = math.log(mid)
            varC = float(((C[Ns[sel]] - C[Ns[sel]].mean()) ** 2).mean())
            rs = varC / float(V[Ns[sel]].mean())
            xs.append(L); ys.append(sd * sd)
            print(f"{b:>10}-{hi:>11} {int(sel.sum()):>8} {sd:>8.4f} "
                  f"{sd*sd:>8.4f} {L:>8.3f} {sd*sd/L:>10.4f} "
                  f"{rs:>14.4f}")
        b *= 2

    xs = np.array(xs); ys = np.array(ys)
    alpha = float(np.polyfit(np.log(xs), np.log(ys), 1)[0])
    kappa = float(np.mean(ys / xs ** alpha))

    def cv(a):
        r = ys / xs ** a
        return float(r.std() / r.mean())

    print(f"\nfits")
    print(f"  alpha = 0 (the recorded law) : CV of sd^2 "
          f"= {cv(0.0)*100:.2f}%")
    print(f"  alpha = 1                    : CV of sd^2/logN "
          f"= {cv(1.0)*100:.2f}%")
    print(f"  fitted alpha = {alpha:.4f}      : CV "
          f"= {cv(alpha)*100:.2f}%   kappa = {kappa:.4f}")

    ok = cv(alpha) < 0.02 and abs(alpha - 1.0) < 0.15
    print("\nverdict:",
          f"CONFIRMED -- Var C(N) = {kappa:.3f} * S(N) * N * log N, "
          f"and the recorded law is missing a factor sqrt(log N)"
          if ok else
          "NOT CONFIRMED -- see the CV column")
    print("consequence: G(N) = C(N)/sqrt(S(N) N) is NOT the standardised")
    print("fluctuation; C(N)/sqrt(kappa S(N) N log N) is. Any statement")
    print("normalised the old way drifts by sqrt(log N) across a range.")
    print("DONE")


if __name__ == "__main__":
    main()
