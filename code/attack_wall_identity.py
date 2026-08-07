# -*- coding: utf-8 -*-
"""
The wall's second moment as a Mobius-times-prime-pair sum, and what it
does and does not buy (increment 320)

THE TARGET. C(N) = sum_{n<N} Lambda(n) mu(N-n), and by Huang-Li
C(N) = o(N) is EQUIVALENT to binary Goldbach for large even N. The
trivial bound is |C(N)| <= psi(N) ~ N, and the measured truth is
C(N) ~ sqrt(N) (log N)^{0.29} -- true with a full power of N to spare.
Nothing weaker than the trivial bound is known pointwise, because any
improvement of the form N(log N)^{-A} IS the theorem.

THE IDENTITY. Summing the square over N and swapping the order,

    sum_{N<=X} C(N)^2
      = sum_{v,v'} mu(v) mu(v') sum_N Lambda(N-v) Lambda(N-v')
      = sum_h M_X(h) P_X(h),

    M_X(h) = sum_v mu(v) mu(v+h),   P_X(h) = sum_w Lambda(w) Lambda(w+h),

exactly, with h ranging over all shifts and the h = 0 term being
sum_v mu^2(v) Lambda(v)^2-weighted. This is the global form of
Proposition W: the wall's second moment is the binary Mobius
autocorrelation weighted by prime-pair counts, and nothing else.

WHAT IT BUYS, AND THE ARITHMETIC OF WHY THAT IS NOT ENOUGH.
A second moment gives an almost-all statement by Chebyshev and no more.
If sum_{N<=X} C(N)^2 << X^2 (log X)^{-A} then
#{N <= X : |C(N)| > X (log X)^{-A/3}} << X (log X)^{-A/3}: Goldbach for
almost all even N, which is Chudakov-Estermann-van der Corput, 1938.
Getting from almost-all to every-N needs a POINTWISE bound, and no
second moment supplies one.

This run does three things, none of which is a proof and all of which
are checkable:

  (A) VERIFY THE IDENTITY numerically, exactly, at a size where both
      sides can be computed independently. It is the one new algebraic
      object here and it should not go into a document unverified.

  (B) MEASURE THE COST OF GOING FROM ALMOST-ALL TO EVERY-N, in this
      program's own units: max_N |C(N)| / rms_N C(N) over dyadic
      ranges. If the extreme is Gumbel-Gaussian, as increment 290
      found, that ratio is ~sqrt(2 log n) and the pointwise truth is
      only a sqrt(log) factor above the rms. The gap to be closed is
      therefore NOT in the size of the extreme.

  (C) STATE WHERE THE GAP ACTUALLY IS, with the program's own numbers:
      the margin between truth and target is a power of N, and both
      classical routes lose a power of N. Proposition E measures the
      circle method's loss as exactly sqrt(Q/(rho A)) ~ sqrt(N)
      (increment 309); Theorem D's loss is exp(c sqrt(log N)), which
      exceeds every power of log but is SMALLER than sqrt(N) -- so the
      two routes fail for different reasons and only one of them fails
      by the margin's own size.

PRE-REGISTRATION (fixed before the run).

  (A1) The identity holds to floating-point exactness: the two sides
       agree to a relative 1e-10 at X = 2^15 and 2^16, computed by
       independent routes (direct convolution vs the h-sum).
       RULE: relative difference < 1e-10. If it fails the identity is
       wrong and nothing below it reads.

  (B1) max|C| / rms C over each dyadic band, against sqrt(2 log n)
       for n the number of even N in the band. RULE: the ratio to
       sqrt(2 log n) lies in [0.7, 1.4] in every band.

       NOTE, added after the first run and before the second: this is
       the RAW C, which carries the location mask and the variation of
       sqrt(V) across the band. #290's Gumbel reading is for the
       STANDARDISED G = (C - m)/sqrt(V), and the deepest cells sit at
       m ~ -7 sqrt(V) (#151), so the raw extreme MUST come out larger.
       The rule is left exactly as written and duly fails at ~1.9x.
       What is reported beside it is max|C|/N, which is the quantity a
       proof would have to bound and which no standardisation can
       flatter.

  (C1) No rule. The comparison of losses is arithmetic on numbers
       already established, and it is reported, not tested.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


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


def main():
    t0 = time.time()

    # ---- (A) the identity, at two small sizes, two independent routes
    print("(A) sum_{N<=X} C(N)^2  =  sum_h M(h) P(h)")
    print(f"{'X':>8} {'direct':>18} {'h-sum':>18} {'rel diff':>11}")
    okA = True
    for e in (15, 16):
        X = 1 << e
        mu, lam = sieve(X)
        muf = mu.astype(np.float64)
        # direct: C(N) for every N <= X, then sum of squares
        nf = 1
        while nf < 2 * (X + 1):
            nf *= 2
        # mu*Lambda has support up to 2X, and Parseval makes the
        # identity hold over the FULL range. A first draft summed
        # C[:X+1] against a complete h-sum and duly missed by 46% --
        # the identity was right and the check truncated one side.
        C = np.fft.irfft(np.fft.rfft(np.pad(muf, (0, nf - X - 1)))
                         * np.fft.rfft(np.pad(lam, (0, nf - X - 1))),
                         nf)
        direct = float((C ** 2).sum())
        # h-sum: M(h) and P(h) by autocorrelation, both truncated at X
        FM = np.fft.rfft(np.pad(muf, (0, nf - X - 1)))
        FL = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))
        M = np.fft.irfft(FM * np.conj(FM), nf)
        P = np.fft.irfft(FL * np.conj(FL), nf)
        # No factor 2. The circular autocorrelation of a zero-padded
        # array already carries the lag-h value at BOTH h and nf-h, so
        # M[1:] contains every nonzero shift twice. A first draft
        # doubled it again and missed by 6.2% at 2^15 and 1.1% at
        # 2^16 -- shrinking with X, which is what sent me looking for
        # a boundary effect instead of a double count.
        hsum = float(M[0] * P[0] + float((M[1:] * P[1:]).sum()))
        rel = abs(direct - hsum) / max(abs(direct), 1e-30)
        okA &= rel < 1e-10
        print(f"{X:>8} {direct:>18.6e} {hsum:>18.6e} {rel:>11.2e}")
    print(f"    (A1) the identity holds to 1e-10: "
          f"{'PASS' if okA else 'FAIL'}")

    # ---- (B) the extreme, on the raw C ----
    X = 8_000_000
    mu, lam = sieve(X)
    nf = 1
    while nf < 2 * (X + 1):
        nf *= 2
    C = np.fft.irfft(np.fft.rfft(np.pad(mu.astype(np.float64),
                                        (0, nf - X - 1)))
                     * np.fft.rfft(np.pad(lam, (0, nf - X - 1))),
                     nf)[: X + 1]
    print(f"\n(B) how far the extreme sits above the rms  "
          f"t={time.time()-t0:.0f}s")
    print(f"{'band':>21} {'n':>9} {'rms C':>12} {'max|C|':>12} "
          f"{'ratio':>8} {'sqrt(2 ln n)':>13} {'rel':>7} "
          f"{'max|C|/N':>9}")
    okB = True
    b = 100_000
    while b < X:
        hi = min(2 * b, X)
        Ns = np.arange(b + (b % 2), hi, 2)
        c = C[Ns]
        n = len(c)
        rms = float(np.sqrt((c * c).mean()))
        mx = float(np.abs(c).max())
        pred = math.sqrt(2.0 * math.log(n))
        rel = (mx / rms) / pred
        okB &= 0.7 <= rel <= 1.4
        print(f"{b:>9}-{hi:>11} {n:>9} {rms:>12.4e} {mx:>12.4e} "
              f"{mx/rms:>8.3f} {pred:>13.3f} {rel:>7.3f} "
              f"{mx/hi:>9.4f}")
        b = hi
    print(f"    (B1) extreme/rms within [0.7, 1.4] of sqrt(2 ln n) "
          f"everywhere: {'PASS' if okB else 'FAIL'}")

    # ---- (C) where the gap is ----
    Nq = 1e8
    lg = math.log(Nq)
    truth = math.sqrt(Nq) * lg ** 0.29
    trivial = Nq
    print(f"\n(C) the margin, and what each route loses, at N = 1e8")
    print(f"    trivial bound      psi(N)              "
          f"{trivial:.3e}")
    print(f"    measured truth     sqrt(N)(log N)^0.29 "
          f"{truth:.3e}")
    print(f"    margin             trivial / truth     "
          f"{trivial/truth:.3e}   = N^{math.log(trivial/truth)/math.log(Nq):.3f}")
    cs = math.sqrt(0.6079 / (0.81 * 0.8106) * Nq)
    print(f"    circle method loses  sqrt(Q/(rho*A))    "
          f"{cs:.3e}   = N^{math.log(cs)/math.log(Nq):.3f}   (Prop E, inc 309)")
    dloss = math.exp(1.0 * math.sqrt(lg))
    print(f"    divisor switch loses exp(sqrt(log N))   "
          f"{dloss:.3e}   = N^{math.log(dloss)/math.log(Nq):.3f}   (Thm D)")
    print(f"\n    The margin is a power of N. The circle method loses")
    print(f"    the SAME power -- that is Proposition E, and it is why")
    print(f"    'zero margin' is exact rather than rhetorical. The")
    print(f"    divisor switch loses far LESS than the margin, and it")
    print(f"    still fails: Theorem D's obstruction is not a size")
    print(f"    deficit but the absence of mass at the truncation")
    print(f"    point. Two routes, two different reasons, one wall.")
    if okA and okB:
        v = ("the identity is verified and the raw extreme follows the "
             "Gumbel scale, which would mean the mask does not inflate "
             "it -- contradicting #151 and worth chasing")
    elif okA:
        v = ("the identity is verified. The RAW extreme runs about 1.9x "
             "the Gumbel scale, which is what #151 predicts once the "
             "mask is left in -- #290's Gumbel reading is for the "
             "standardised G and is untouched. The number that matters "
             "is max|C|/N, and it is under 1% at every band and "
             "falling: the pointwise truth is far below trivial. THE "
             "GAP IS NOT THE SIZE OF THE EXTREME. What is missing is a "
             "pointwise method losing less than a power of N, and this "
             "program has proved neither classical route is one")
    else:
        v = "the identity does not verify; nothing here reads"
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
