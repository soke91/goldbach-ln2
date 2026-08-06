# -*- coding: utf-8 -*-
"""
What is left after the mask (increment 249): truncation, or skewness,
or something else?

Increment 248 removed the enumerated cell mean over q in {3,...,23} and
the sign balance of C(N) returned from z = -94.3 and +61.6 to +4.83 and
-5.17. That residual has three candidate explanations and they are
distinguishable.

 (i)  TRUNCATION. The enumeration stops at q = 23; primes 29, 31, ...
      dividing N are not centred out. If so the residual should shrink
      as the truncation is raised.
 (ii) SKEWNESS. Subtracting a cell MEAN forces the mean to zero, not
      the median. A skewed distribution with zero mean has
      P(X > 0) != 1/2, so the residual could be pure skewness and
      nothing to do with a missing mask at all. This alternative was
      not considered in increment 248 and is the more likely of the
      two, since it needs no missing structure.
 (iii) SOMETHING ELSE, which is what is left if neither of the above
      accounts for it.

DESIGN.
 (A) Sweep the truncation: enumerate over QS[:k] for k = 1..9, i.e.
     {3}, {3,5}, ..., {3,...,37}, subtract cell means with the sparse
     cell fallback of increment 240 (a cell with fewer than 30 members
     takes the mean of the same pattern minus its rarest prime, so that
     noise is not laundered into the prediction), and report the
     residual sign-balance z at each k.
     PREDICTION for (i): |z| decreases monotonically toward 0.
     PREDICTION for (iii): |z| plateaus.
 (B) Measure the within-cell skewness directly, and compare the
     residual sign imbalance against what that skewness alone predicts.
     For a distribution with mean 0 the Edgeworth expansion gives
     F(0) ~= 1/2 + gamma_1 phi(0)/6, hence

         P(X > 0) - 1/2  ~=  - gamma_1 / (6 sqrt(2 pi)),

     with a MINUS sign: negative skewness pushes the median above the
     mean, so more values land above 0, not fewer. (Written with a plus
     on the first pass, which inverted every verdict; the corrected
     form is used below and the error is recorded rather than quietly
     repaired.) A measured gamma_1 therefore predicts a specific
     imbalance with no free parameter.
     PREDICTION for (ii): the observed residual matches that.

NULLS. After perfect centring the sign-balance z is standard normal, so
|z| < 3 is "explained". Each skewness carries SE sqrt(6/n) on its cell.
The Edgeworth prediction has no fitted constant.
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


def centred(Cv, Ns, QS):
    cell = np.zeros(len(Ns), dtype=np.int64)
    for i, q in enumerate(QS):
        cell |= ((Ns % q == 0).astype(np.int64) << i)
    ncell = 1 << len(QS)
    cnt = np.bincount(cell, minlength=ncell)
    mean = np.zeros(ncell)
    for c in range(ncell):
        if cnt[c] >= 30:
            mean[c] = float(Cv[cell == c].mean())
    for c in range(ncell):                    # sparse-cell fallback
        if cnt[c] < 30 and cnt[c] > 0:
            cc = c
            for i in range(len(QS) - 1, -1, -1):
                if cc >> i & 1:
                    cc &= ~(1 << i)
                    if cnt[cc] >= 30:
                        break
            mean[c] = mean[cc]
    return Cv - mean[cell], cell


def zbal(x):
    n = len(x)
    p = float((x > 0).mean())
    return (p - 0.5) * 2 * math.sqrt(n), p, n


def main():
    X = 4_000_000
    lo = 100_000
    ALL = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    t0 = time.time()
    mu, lam, primes = sieve(X)
    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    F = np.zeros(n_fft); G = np.zeros(n_fft)
    F[: X + 1] = mu; G[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(F) * np.fft.rfft(G), n_fft)[: X + 1]
    del F, G
    Ns = np.arange(lo, X + 1, 2)
    Cv = C[Ns]
    d3 = Ns % 3 == 0
    print(f"convolution t={time.time()-t0:.0f}s", flush=True)

    print(f"\n(A) raising the truncation")
    print(f"{'QS':>34} {'cells':>7} {'z all':>9} {'z 3|N':>9} "
          f"{'z 3 not|N':>11}")
    z0, _, _ = zbal(Cv)
    za, _, _ = zbal(Cv[d3]); zb, _, _ = zbal(Cv[~d3])
    print(f"{'(none, raw C)':>34} {'-':>7} {z0:>9.2f} {za:>9.2f} "
          f"{zb:>11.2f}")
    for k in range(1, len(ALL) + 1):
        QS = ALL[:k]
        Cc, _ = centred(Cv, Ns, QS)
        z1, _, _ = zbal(Cc)
        z2, _, _ = zbal(Cc[d3]); z3, _, _ = zbal(Cc[~d3])
        lab = ",".join(str(q) for q in QS)
        print(f"{lab:>34} {1<<k:>7} {z1:>9.2f} {z2:>9.2f} "
              f"{z3:>11.2f}   t={time.time()-t0:.0f}s", flush=True)
    print("    truncation => |z| falls toward 0; something else => it")
    print("    plateaus. Read the 3|N and 3 not|N columns, not the")
    print("    pooled one, which cancels (increment 248)")

    print(f"\n(B) is the residual just skewness?")
    QS = ALL[:8]
    Cc, cell = centred(Cv, Ns, QS)
    print(f"{'group':>18} {'n':>9} {'skew g1':>9} {'SE':>7} "
          f"{'P-1/2 meas':>11} {'Edgeworth':>10} {'ratio':>7}")
    for tag, sel in (("all", np.ones(len(Ns), bool)),
                     ("3 | N", d3), ("3 not| N", ~d3)):
        x = Cc[sel]
        n = len(x)
        sd = float(x.std())
        g1 = float(((x - x.mean()) ** 3).mean() / sd ** 3)
        pm = float((x > 0).mean()) - 0.5
        edge = -g1 / (6 * math.sqrt(2 * math.pi))
        print(f"{tag:>18} {n:>9} {g1:>+9.4f} {math.sqrt(6/n):>7.4f} "
              f"{pm:>+11.5f} {edge:>+10.5f} "
              f"{pm/edge if edge else float('nan'):>7.2f}")
    print("    Edgeworth has NO fitted constant: for mean-zero data")
    print("    P(X>0) - 1/2 = -g1/(6 sqrt(2 pi)), the MINUS mattering")
    print("    because negative skew lifts the median above the mean.")
    print("    A ratio near 1 means the residual IS skewness and no")
    print("    mask is missing. The pooled row mixes two populations")
    print("    with opposite signs and is not readable (increment 248).")
    print("DONE")


if __name__ == "__main__":
    main()
