# -*- coding: utf-8 -*-
"""
The mask as a second-order term in the Goldbach count (increment 252).

LOCATION_MASK.md establishes that C(N) = Sum_n Lambda(n) mu(N-n) is not
mean-zero: it carries a deterministic term of sqrt(N) scale whose
coefficient is set by which small primes divide N, reaching about
-20 sqrt(N) for N built from the first seven primes.

Theorem C of this campaign is the unconditional identity

    r~(N) = Sum_{n<N} Lambda(n) Lambda(N-n)
          = S(N) ( N - C(N) )  +  E_3(N)  +  O_A(N (log N)^{-A}).

So a deterministic term in C(N) is a deterministic term in the GOLDBACH
COUNT ITSELF, at second order:

    R(N) := r~(N) - S(N) N  ~=  - S(N) C(N)  +  E_3(N).

For deep N that predicts a POSITIVE excess of about S(N) * 20 sqrt(N)
~ 100 sqrt(N) -- more Goldbach representations than S(N) N, by a
computable amount depending on rad(N). This connects the mask to an
observable that has nothing to do with mu, and it is the sharpest
external check available on the whole thread.

WHAT IS TESTED, WITH NULLS ON THE SAME LINE.
 (A) The regression of R(N) on -S(N) C(N), over every even N in the
     range. PREDICTION: slope 1, since Theorem C is an identity and
     the only other term is E_3. NULL: slope 0, i.e. R does not track
     the mask at all. A permutation control (N-labels shuffled) is run
     on the identical statistic and printed whether or not it flatters.
 (B) The same by depth, since the mask is largest where the most small
     primes divide N and that is where the prediction is sharpest.
 (C) The residual E_3(N) = R(N) + S(N) C(N), reported as a power of N,
     against MEASUREMENTS section 9's recorded N^0.599. It must stay
     smaller than the term being tested, or (A) is measuring E_3 rather
     than the mask; at N = 4e6, N^0.6 is about 9e3 while S C is about
     2e5 for the deep cells, so there is room -- but it is checked and
     not assumed.

LIMIT, STATED. Theorem C carries an O_A(N (log N)^{-A}) which at these
N is not small in absolute terms. The test is therefore of the
CORRELATION and SLOPE between R and -S C, which that error term cannot
manufacture: it depends on N smoothly and not on rad(N).
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
    return mu, lam, primes


def conv(a, b, n_fft, X):
    A = np.zeros(n_fft); A[: X + 1] = a
    B = np.zeros(n_fft); B[: X + 1] = b
    return np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]


def main():
    X = 4_000_000
    lo = 100_000
    QS = [3, 5, 7, 11, 13, 17, 19, 23]
    t0 = time.time()
    mu, lam, primes = sieve(X)
    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    rt = conv(lam, lam, n_fft, X)            # r~(N) = Lambda * Lambda
    C = conv(mu.astype(np.float64), lam, n_fft, X)
    print(f"convolutions t={time.time()-t0:.0f}s", flush=True)

    C2 = 0.6601618158468696
    S = np.full(X + 1, 2 * C2)
    for p in primes:
        p = int(p)
        if p > 2:
            S[p::p] *= (p - 1) / (p - 2)

    Ns = np.arange(lo, X + 1, 2)
    R = rt[Ns] - S[Ns] * Ns
    P = -S[Ns] * C[Ns]                       # the predicted second order
    E3 = R - P

    def slope(y, x):
        xm, ym = x.mean(), y.mean()
        return float(np.dot(x - xm, y - ym) / np.dot(x - xm, x - xm))

    rng = np.random.default_rng(20260806)
    Pp = P[rng.permutation(len(P))]

    print(f"\n(A) does the Goldbach count track the mask?")
    print(f"  corr(R, -S C)  = {float(np.corrcoef(R, P)[0,1]):+.4f}")
    print(f"  slope          = {slope(R, P):+.4f}"
          f"     (predicted 1, null 0)")
    print(f"  permuted corr  = {float(np.corrcoef(R, Pp)[0,1]):+.4f}"
          f"     (null 0, SE {1/math.sqrt(len(R)):.5f})")
    print(f"  permuted slope = {slope(R, Pp):+.4f}")
    print(f"  n = {len(R)}")

    print(f"\n(B) by depth -- how many of {{3..23}} divide N")
    nsm = np.zeros(len(Ns), dtype=np.int8)
    for q in QS:
        nsm += (Ns % q == 0).astype(np.int8)
    print(f"{'depth':>6} {'n':>8} {'mean R':>12} {'mean -S C':>12} "
          f"{'ratio':>8} {'slope':>8} {'corr':>8}")
    for k in range(0, 6):
        sel = nsm == k
        if sel.sum() < 300:
            continue
        r_, p_ = R[sel], P[sel]
        print(f"{k:>6} {int(sel.sum()):>8} {r_.mean():>12.1f} "
              f"{p_.mean():>12.1f} "
              f"{r_.mean()/p_.mean() if p_.mean() else float('nan'):>8.3f} "
              f"{slope(r_, p_):>8.3f} "
              f"{float(np.corrcoef(r_, p_)[0,1]):>8.4f}")

    print(f"\n(C) the residual E_3 = R + S C, against its recorded size")
    print(f"{'N band':>20} {'n':>8} {'rms E_3':>12} {'rms S C':>12} "
          f"{'E3/SC':>8} {'log E3/log N':>13}")
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if sel.sum() > 1000:
            e = float(np.sqrt((E3[sel] ** 2).mean()))
            s = float(np.sqrt((P[sel] ** 2).mean()))
            mid = math.sqrt(b * hi)
            print(f"{b:>9}-{hi:>10} {int(sel.sum()):>8} {e:>12.1f} "
                  f"{s:>12.1f} {e/s:>8.3f} "
                  f"{math.log(e)/math.log(mid):>13.4f}")
        b *= 2
    print("    MEASUREMENTS section 9 recorded the discrepancy as")
    print("    ~ N^0.599; the last column is the same exponent measured")
    print("    here. If E_3 dominates S C, (A) is measuring E_3.")
    print("DONE")


if __name__ == "__main__":
    main()
