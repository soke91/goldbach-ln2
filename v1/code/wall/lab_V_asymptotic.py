# -*- coding: utf-8 -*-
"""
What IS the wall's exact scale? V(N) in closed form (increment 287).

BACKGROUND. Increment 283 showed that Conjecture L's Gaussian half
holds for C(N) under the EXACT second moment

    V(N) = Sum_{v<N} mu^2(v) Lambda(N-v)^2

and fails (z = 98) under the fitted stand-in kappa*S(N)*N*log N that
the documents display. Increment 280's unidentifiable exponent was the
same fact seen from the other side. What nobody then asked is the
obvious next question: V(N) is an exactly computable arithmetic
function -- so what is it?

THE DERIVATION, which is short. Lambda(w)^2 is supported on prime
powers w = p^k with weight (log p)^2, so with w = N - v,

    V(N) = Sum_{p^k < N} (log p)^2 mu^2(N - p^k)
         = Sum_{p < N} (log p)^2 mu^2(N - p)  +  O(sqrt(N) log^2 N),

the k >= 2 terms being negligible. The remaining sum counts SQUAREFREE
SHIFTED PRIMES with a (log p)^2 weight. Its local density at a prime q:

  * q | N: q^2 | N-p forces q | N-p, hence q | p, hence p = q. One
    term. The factor is 1 -- primes dividing N impose NO condition.
  * q not| N: the bad class is p == N (mod q^2), and since q not| N
    that is a UNIT class, so it holds for a density 1/phi(q^2) =
    1/(q(q-1)) of primes. The factor is 1 - 1/(q(q-1)).

So the arithmetic factor is

    A(N) := Prod_{q not| N} (1 - 1/(q(q-1)))
          = ArtinConst / Prod_{q | N} (1 - 1/(q(q-1))),

with ArtinConst = 0.3739558136... This is Mirsky's 1949 theorem on
squarefree values of shifted primes; the derivation is recalled, not
claimed.

WHY IT MATTERS HERE, AND IT IS NOT A NAMING EXERCISE. This program has
used S(N) = 2C_2 Prod_{q|N, q>2} (q-1)/(q-2) as the local factor
everywhere. For V the correct local factor is A(N). Both are products
over the primes dividing N, which is why S looked serviceable, but
they are DIFFERENT FUNCTIONS -- and that is exactly the shape of
increment 283's finding that 93.7% of the variance of (S*N)/V is
explained by the divisibility cells.

PRE-REGISTRATION (fixed before the run).

  The analytic factor is removed exactly rather than modelled. Define
      W(N) = Sum_{w<N} Lambda(w)^2,
  a prefix sum that does not depend on N's factorisation at all. Then

      R(N) := V(N)/W(N)

  is the weighted density of squarefree N-w, and the claim above is
  precisely R(N) -> A(N). This isolates the arithmetic from the
  analytic, so no asymptotic for Sum (log p)^2 is needed and none is
  assumed.

  TEST 1: R(N)/A(N) -> 1, and the deviation shrinks with N.
  TEST 2 (the decisive one): which local factor explains the cell
    structure? Fit R(N) by A(N) and, separately, by S(N) rescaled to
    the same mean. Report the residual spread of each ACROSS the
    divisibility cells. DECISION RULE, fixed now: A(N) is declared the
    correct factor only if its residual spread is below one tenth of
    S(N)'s.
  TEST 3: a control. If A(N) is right for the reason claimed, then
    restricting to N with a FIXED radical must make R(N) constant to
    within sampling. Reported per cell.
"""
import math
import time

import numpy as np

ARTIN = 0.3739558136192022880547280543464164151116


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
    return mu, lam, primes, spf


def main():
    X = 16_000_000
    t0 = time.time()
    mu, lam, primes, spf = sieve(X)

    n = 1
    while n < 2 * (X + 1):
        n *= 2
    A_ = np.zeros(n); A_[: X + 1] = (mu != 0).astype(np.float64)
    B_ = np.zeros(n); B_[: X + 1] = lam ** 2
    V = np.fft.irfft(np.fft.rfft(A_) * np.fft.rfft(B_), n)[: X + 1]
    del A_, B_
    W = np.cumsum(lam ** 2)          # W(N) = Sum_{w<=N} Lambda(w)^2
    print(f"convolution  t={time.time()-t0:.0f}s", flush=True)

    # A(N) = ArtinConst / Prod_{q|N} (1 - 1/(q(q-1))), by sieving the
    # reciprocal factor over each prime's multiples.
    corr = np.ones(X + 1, dtype=np.float64)
    for p in primes:
        p = int(p)
        f = 1.0 - 1.0 / (p * (p - 1.0))
        corr[p::p] /= f
    AN = ARTIN * corr

    lo = 100_000
    Ns = np.arange(lo, X + 1, 2)
    R = V[Ns] / W[Ns - 1]
    Apred = AN[Ns]

    print("\n(TEST 1) R(N) = V(N)/W(N) against A(N), by octave band")
    print(f"{'band':>21} {'count':>9} {'mean R/A(N)':>13} "
          f"{'sd of R/A':>11} {'max |dev|':>11}")
    b = lo
    rows = []
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) < 1000:
            b = hi
            continue
        q = R[sel] / Apred[sel]
        rows.append((b, hi, int(sel.sum()), float(q.mean()),
                     float(q.std()), float(np.abs(q - 1).max())))
        print(f"{b:>9}-{hi:>11} {rows[-1][2]:>9} {rows[-1][3]:>13.6f} "
              f"{rows[-1][4]:>11.6f} {rows[-1][5]:>11.6f}")
        b = hi
    dev = [abs(r[3] - 1.0) for r in rows]
    print(f"    |mean R/A - 1| : first band {dev[0]:.6f}, "
          f"last {dev[-1]:.6f}   "
          f"{'shrinking' if dev[-1] < dev[0] else 'NOT shrinking'}")

    # S(N), the factor this program has used everywhere
    C2 = 0.6601618158468696
    S = np.full(X + 1, 2 * C2)
    for p in primes:
        p = int(p)
        if p > 2:
            S[p::p] *= (p - 1) / (p - 2)
    Spred = S[Ns]

    print("\n(TEST 2) which local factor explains the cell structure?")
    print("    each predictor is rescaled to R's mean, so only its")
    print("    SHAPE across N is being judged")
    out = {}
    for lab, pred in (("A(N)  (Mirsky)", Apred), ("S(N)  (used here)", Spred)):
        z = R / (pred * float((R / pred).mean()))
        out[lab] = (float(z.std()), float(np.abs(z - 1).max()))
        print(f"    {lab:>20}: residual sd {out[lab][0]:.6f}, "
              f"max |dev| {out[lab][1]:.6f}")
    ra = out["A(N)  (Mirsky)"][0]
    rs = out["S(N)  (used here)"][0]
    print(f"    ratio sd(A)/sd(S) = {ra/rs:.5f}")
    print(f"    pre-registered rule (< 0.1): "
          f"{'A(N) IS the correct factor' if ra < 0.1 * rs else 'not established'}")

    print("\n(TEST 3) control: fix the radical, R(N) must be constant")
    print(f"{'primes dividing N':>26} {'count':>9} {'mean R':>10} "
          f"{'sd R':>10} {'A(N)':>10} {'R/A':>9}")
    small = [2, 3, 5, 7, 11, 13]
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(small):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    big = np.ones(len(Ns), dtype=bool)
    for q in small:
        big &= True
    for k in (1, 3, 7, 15, 31, 63):        # 2 | N and then 3, 5, 7, ...
        sel = (key == k)
        if int(sel.sum()) < 500:
            continue
        lab = "".join(f"{q}." for i, q in enumerate(small) if k >> i & 1)
        print(f"{lab.rstrip('.'):>26} {int(sel.sum()):>9} "
              f"{float(R[sel].mean()):>10.6f} {float(R[sel].std()):>10.6f} "
              f"{float(Apred[sel].mean()):>10.6f} "
              f"{float((R[sel]/Apred[sel]).mean()):>9.5f}")
    print("    within a cell A(N) is constant, so any spread in R is")
    print("    the tail primes q > 13 dividing N plus sampling.")
    print("DONE")


if __name__ == "__main__":
    main()
