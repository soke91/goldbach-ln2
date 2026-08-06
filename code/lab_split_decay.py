# -*- coding: utf-8 -*-
"""
How fast does the 3 | N mean split actually decay? (increment 254)

MEASUREMENTS section 10, hypothesis H2, found the means of C(N) split
by 3 | N at about 6.7 sigma, chased it in three disjoint windows, and
concluded:

  "Real but decaying. ... the split measured in sqrt(N) units is 1.007,
   0.719, 0.536, scaling as N^-0.41 -- i.e. the absolute split grows
   only like N^0.09. It is a lower-order term, NOT A MASK on
   C(N)/sqrt(N), and at the largest window normalising by S(N) already
   takes it below 3 sigma."

Increments 239-253 establish that it IS the mask, and correction #45
shows that normalising by S(N) does not remove it. What has not been
checked is the DECAY RATE, which is the whole basis of the "lower-order
term" reading -- and it was fitted to three points carrying z = -4.38,
-3.17 and -2.29, i.e. relative errors of roughly 23, 32 and 44 percent.
A two-parameter power law fitted to three points of that quality
determines very little.

THE TWO HYPOTHESES, both stated before measuring.
  RECORDED   split/sqrt(N) ~ N^-0.41. Over N from 3e5 to 1.4e6 that is
             a fall by a factor 0.53.
  MASK       LOCATION_MASK.md measures C/sqrt(N) for the deep family
             falling like (log N)^-1.5, which over the same range is a
             fall by a factor 0.84 -- nearly flat by comparison.
  The three recorded points cannot separate 0.53 from 0.84; full power
  can.

WHAT IS MEASURED. The split mean(C | 3|N) - mean(C | 3 not|N) over
every even N in [1e5, 4e6], by dyadic band, in three normalisations:
raw, /sqrt(N), and /sqrt(N log N). Each with the standard error of a
difference of means computed from the same data, so the z is honest.
Then the decay exponent is fitted in N and in log N and both are
reported against the recorded -0.41.

NULLS ON THE SAME LINE. The null for each split is 0 with SE
sqrt(var1/n1 + var2/n2). A permutation control -- 3 | N labels shuffled
within each band -- is run on the identical statistic and printed for
every band, so the design's own floor is visible and cannot be skipped.
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


def zof(x, m):
    a, c = x[m], x[~m]
    se = math.sqrt(a.var(ddof=1)/len(a) + c.var(ddof=1)/len(c))
    return float(a.mean() - c.mean()) / se


def main():
    X = 4_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    A = np.zeros(n_fft); B = np.zeros(n_fft)
    A[: X + 1] = mu; B[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]
    del A, B
    print(f"convolution t={time.time()-t0:.0f}s", flush=True)

    C2c = 0.6601618158468696
    S = np.full(X + 1, 2 * C2c)
    for p_ in primes:
        p_ = int(p_)
        if p_ > 2:
            S[p_::p_] *= (p_ - 1) / (p_ - 2)
    Ns = np.arange(lo, X + 1, 2)
    Cv = C[Ns]
    Sv = S[Ns]
    d3 = Ns % 3 == 0
    rng = np.random.default_rng(20260806)

    print(f"\n{'band':>21} {'n':>9} {'split':>10} {'SE':>9} {'z':>9} "
          f"{'perm z':>8} {'/sqrtN':>9} {'z of C/sqrtS':>13} "
          f"{'z of C/S':>9}")
    mids, spl = [], []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if sel.sum() > 5000:
            a = Cv[sel & d3]; c = Cv[sel & ~d3]
            s = float(a.mean() - c.mean())
            se = math.sqrt(a.var(ddof=1) / len(a) + c.var(ddof=1) / len(c))
            lab = rng.permutation(np.concatenate(
                [np.ones(len(a), bool), np.zeros(len(c), bool)]))
            allv = np.concatenate([a, c])
            sp = float(allv[lab].mean() - allv[~lab].mean())
            mid = math.sqrt(b * hi)
            mids.append(mid); spl.append(s)
            print(f"{b:>9}-{hi:>10} {int(sel.sum()):>9} {s:>10.2f} "
                  f"{se:>9.2f} {s/se:>9.2f} {sp/se:>8.2f} "
                  f"{s/math.sqrt(mid):>9.4f} "
                  f"{zof(Cv[sel]/np.sqrt(Sv[sel]), d3[sel]):>13.2f} "
                  f"{zof(Cv[sel]/Sv[sel], d3[sel]):>9.2f}")
        b *= 2

    mids = np.array(mids); spl = np.array(spl)
    y = np.log(-spl / np.sqrt(mids))
    aN = float(np.polyfit(np.log(mids), y, 1)[0])
    aL = float(np.polyfit(np.log(np.log(mids)), y, 1)[0])
    print(f"\n  fitted decay of split/sqrt(N):")
    print(f"    as N^a       a = {aN:+.4f}    "
          f"[MEASUREMENTS section 10 recorded -0.41]")
    print(f"    as (log N)^a a = {aL:+.4f}    "
          f"[LOCATION_MASK's deep-family exponent is -1.5]")
    r0 = -spl[0] / math.sqrt(mids[0])
    r1 = -spl[-1] / math.sqrt(mids[-1])
    print(f"  measured fall across the range: {r1/r0:.3f}  over a")
    print(f"  factor {mids[-1]/mids[0]:.1f} in N")
    print(f"    N^-0.41 would give {(mids[-1]/mids[0])**-0.41:.3f}")
    print(f"    (log N)^-1.5 would give "
          f"{(math.log(mids[-1])/math.log(mids[0]))**-1.5:.3f}")
    print("DONE")


if __name__ == "__main__":
    main()
