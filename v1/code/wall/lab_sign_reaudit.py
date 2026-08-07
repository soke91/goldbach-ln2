# -*- coding: utf-8 -*-
"""
Re-audit of the sign-balance test against the location mask
(increment 248).

sweep_B item B4 measured the sign balance of C(N) against a null of
1/2 and read P(C>0) = 0.4800 at z = -1.55, i.e. "no signal". That null
assumed C(N) has zero mean. Increments 239-247 established that it does
not: C(N) carries a deterministic term set by which small primes divide
N. So B4's reading was UNDERPOWERED rather than null, and it was taken
on about 1500 values of N.

Here it is redone on 1.95 x 10^6 values, where the standard error of a
proportion is 1/(2 sqrt n) = 3.6 x 10^-4 rather than 1.3 x 10^-2 -- a
factor 36 more resolution.

AND IT CARRIES A SHARP PREDICTION. Increment 247 found that the mask
does not vanish off 3 | N but CHANGES SIGN: mean R_A = -6.10e-4 for
3 | N against +2.45e-4 for 3 not| N. A deterministic mean shifts a
sign balance in its own direction, so:

  PREDICTION 1  P(C > 0) < 1/2 when 3 | N
  PREDICTION 2  P(C > 0) > 1/2 when 3 not| N
  PREDICTION 3  the gap between the two is larger than either
                deviation from 1/2 on its own

Both predictions are directional and were fixed before the measurement.
The null for each is 1/2 with SE 1/(2 sqrt n) on that subgroup, printed
on the same line.

 (D) COMPLETENESS CHECK. After subtracting the enumerated cell mean --
     the mask of increment 240 -- the sign balance should return to 1/2
     within its error. If it does not, the mask is not the whole of the
     deterministic part. NULL: 1/2 with the same SE.
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


def report(tag, x):
    n = len(x)
    p = float((x > 0).mean())
    se = 1.0 / (2 * math.sqrt(n))
    print(f"{tag:>34} {n:>9} {p:>9.5f} {se:>9.5f} "
          f"{(p - 0.5)/se:>9.2f}")
    return p, (p - 0.5) / se


def main():
    X = 4_000_000
    lo = 100_000
    QS = [3, 5, 7, 11, 13, 17, 19, 23]
    t0 = time.time()
    mu, lam, primes = sieve(X)
    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    F = np.zeros(n_fft); G = np.zeros(n_fft)
    F[: X + 1] = mu; G[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(F) * np.fft.rfft(G), n_fft)[: X + 1]
    del F, G
    print(f"convolution t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(lo, X + 1, 2)
    Cv = C[Ns]

    print(f"\n(A) B4 redone at full power")
    print(f"{'group':>34} {'n':>9} {'P(C>0)':>9} {'SE':>9} {'z':>9}")
    report("all even N  [B4's statistic]", Cv)
    p3, z3 = report("3 | N        [predicted < 1/2]", Cv[Ns % 3 == 0])
    p3n, z3n = report("3 not| N     [predicted > 1/2]",
                      Cv[Ns % 3 != 0])
    print(f"  sweep_B's original reading was P = 0.4800, z = -1.55 on")
    print(f"  about 1500 values -- an SE of 1.3e-2 against 3.6e-4 here")

    print(f"\n(B) the gap between the two, and its own null")
    dn = min(int((Ns % 3 == 0).sum()), int((Ns % 3 != 0).sum()))
    se_gap = math.sqrt(0.25 / (Ns % 3 == 0).sum()
                       + 0.25 / (Ns % 3 != 0).sum())
    print(f"  P(3|N) - P(3 not|N) = {p3 - p3n:+.5f}"
          f"   SE {se_gap:.5f}   z = {(p3-p3n)/se_gap:+.2f}")
    print(f"  null 0 (no difference between the two populations)")

    print(f"\n(C) by depth -- how many small primes divide N")
    print(f"{'group':>34} {'n':>9} {'P(C>0)':>9} {'SE':>9} {'z':>9}")
    nsm = np.zeros(len(Ns), dtype=np.int8)
    for q in QS:
        nsm += (Ns % q == 0).astype(np.int8)
    for k in range(0, 6):
        sel = nsm == k
        if sel.sum() > 500:
            report(f"{k} of {{3..23}} divide N", Cv[sel])

    print(f"\n(D) completeness -- after removing the enumerated mask")
    cell = np.zeros(len(Ns), dtype=np.int64)
    for i, q in enumerate(QS):
        cell |= ((Ns % q == 0).astype(np.int64) << i)
    ncell = 1 << len(QS)
    Cc = Cv.copy()
    for c in range(ncell):
        idx = np.nonzero(cell == c)[0]
        if len(idx) >= 30:
            Cc[idx] = Cv[idx] - Cv[idx].mean()
    print(f"{'group':>34} {'n':>9} {'P(C>0)':>9} {'SE':>9} {'z':>9}")
    report("all even N, mask removed", Cc)
    report("3 | N,      mask removed", Cc[Ns % 3 == 0])
    report("3 not| N,   mask removed", Cc[Ns % 3 != 0])
    print("  if these return to 1/2 the enumerated mask is the whole")
    print("  of the deterministic part; if not, something is left")
    print("DONE")


if __name__ == "__main__":
    main()
