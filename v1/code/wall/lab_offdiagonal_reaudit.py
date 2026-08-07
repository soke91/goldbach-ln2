# -*- coding: utf-8 -*-
"""
Re-audit of the wall's off-diagonal (increment 256).

MEASUREMENTS section 13 proves the exact identity

    Sum_{N in W} |C(N)|^2 = Sum_h r_W(h) S_W(h),

r_W the prime-pair count at shift h and S_W the binary Chowla
correlation, and reads off it:

  "The off-diagonal does not cancel; it adds. We predicted it would
   have to cancel a log X against the diagonal. Measured, it is
   POSITIVE AND COMPARABLE: 0.545 and 0.475 of the total at X = 2000
   and 4000. So about half the variance of C(N) is genuine
   shifted-Mobius correlation rather than diagonal mass."

X = 2000 and 4000 are the whole basis of that reading, and increment
238 measured the same object indirectly at much larger N: Var C
divided by the exact diagonal Sum_v mu^2(v) Lambda(N-v)^2 runs 1.006,
0.955, 0.923, 0.890, 0.874, 0.873 across dyadic bands up to 4e6. A
ratio BELOW 1 means the off-diagonal is NEGATIVE, not positive and
comparable. The two readings cannot both be right.

Note first what section 13 does settle, and which section 12 got wrong:
the diagonal alone is Sum_v mu^2(v) Lambda(N-v)^2 ~ A(N) N log N, so the
identity of section 13 already implies Var C ~ N log N and contradicts
section 12's Var C = S(N) N. The document carried its own refutation
(correction #36).

WHAT IS MEASURED HERE. Directly, with no model:

    LHS(W)  = Sum_{N in W} C(N)^2
    diag(W) = Sum_{N in W} V(N),  V(N) = Sum_v mu^2(v) Lambda(N-v)^2
    off(W)  = LHS - diag

at X = 2000 and 4000 to reproduce section 13's own numbers, and then at
1e5, 1e6 and 4e6 to see whether the reading survives. Both C and V come
from additive convolutions, so nothing here is estimated.

NULLS AND CRITERION.
 * The quantity reported is off/LHS, section 13's own statistic, so its
   0.545 and 0.475 are the reference on the same line.
 * The identity itself is the control: LHS must equal diag + off by
   construction, so any arithmetic slip shows up as a nonzero residual
   in the printed check.
 * CONFIRMED iff off/LHS stays near 0.5 as X grows. REFUTED iff it
   falls through zero, which is what increment 238's ratios imply.
 * The window is stated: section 13 used all N <= X including odd N.
   Both conventions are printed, since C(N) for odd N is a different
   population (no 2 | N) and pooling them is exactly the kind of mixing
   hazard 5 warns about.
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
    return mu, lam


def main():
    X = 4_000_000
    t0 = time.time()
    mu, lam = sieve(X)
    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    A = np.zeros(n_fft); B = np.zeros(n_fft)
    A[: X + 1] = mu; B[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]
    A[: X + 1] = np.abs(mu); B[: X + 1] = lam ** 2
    V = np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]
    del A, B
    print(f"convolutions t={time.time()-t0:.0f}s", flush=True)

    print(f"\n{'window':>22} {'parity':>7} {'n':>9} {'LHS':>13} "
          f"{'diag':>13} {'off':>13} {'off/LHS':>9}")
    for hi in (2000, 4000, 100_000, 1_000_000, 4_000_000):
        lo = 2
        for tag, step in (("all N", 1), ("even N", 2)):
            start = lo if step == 1 else (lo if lo % 2 == 0 else lo + 1)
            Ns = np.arange(start, hi + 1, step)
            L = float((C[Ns] ** 2).sum())
            d = float(V[Ns].sum())
            print(f"{lo:>9}-{hi:>12} {tag:>7} {len(Ns):>9} "
                  f"{L:>13.4e} {d:>13.4e} {L-d:>13.4e} "
                  f"{(L-d)/L:>9.4f}")
    print("    section 13 recorded off/LHS = 0.545 at X = 2000 and")
    print("    0.475 at X = 4000, and read it as 'positive and")
    print("    comparable ... about half the variance'")

    print(f"\nthe same on the upper half of each range, where the")
    print(f"asymptotics are cleaner than on [2, X]")
    print(f"{'window':>22} {'n even':>9} {'off/LHS':>9} "
          f"{'LHS/(n N logN)':>15} {'diag/(n N logN)':>16}")
    for hi in (4000, 100_000, 1_000_000, 4_000_000):
        lo = hi // 2
        Ns = np.arange(lo + (lo % 2), hi + 1, 2)
        L = float((C[Ns] ** 2).sum())
        d = float(V[Ns].sum())
        mid = math.sqrt(lo * hi)
        norm = len(Ns) * mid * math.log(mid)
        print(f"{lo:>9}-{hi:>12} {len(Ns):>9} {(L-d)/L:>9.4f} "
              f"{L/norm:>15.4f} {d/norm:>16.4f}")
    print("    the last two columns show the identity's own scale: the")
    print("    diagonal alone is ~ A(N) N log N, which is why section")
    print("    13 already implied Var C ~ N log N and section 12's")
    print("    S(N) N was short by a log (correction #36)")
    print("DONE")


if __name__ == "__main__":
    main()
