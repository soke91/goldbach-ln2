# -*- coding: utf-8 -*-
"""
Hypothesis H (increment 208): the final-boss scalar IS the Goldbach
deficit that this program's earliest corpus already measured.

Theorem C gives, unconditionally,
    E_3(alpha) = r(N) - S(N) (N - C(N)) + O_A(N (log N)^{-A}),
with r(N) = Sum_{n<N} Lambda(n)Lambda(N-n),
     C(N) = Sum_{n<N} Lambda(n) mu(N-n),
     S(N) the singular series.
Rearranged:

    C(N) = N - r(N)/S(N) + E_3/S(N) + O_A(...),

so up to the E_3 term, C(N)/N is exactly the RELATIVE DEFICIT of the
Goldbach count against its Hardy-Littlewood prediction. That deficit
is what the ln 2 phase of this program measured at length (the comet,
the dispersion constant 0.693). The two halves of the campaign have
never been connected; if they are the same object, the wall's own
scalar already has a large empirical dossier.

TESTS
  (a) does C(N) track the deficit D(N) := N - r(N)/S(N)?
  (b) how big is C(N)/N, and does it shrink? (the chain needs o(N))
  (c) the residual C(N) - D(N) is E_3/S(N) up to O_A -- measure it.

PRE-REGISTERED (fixed before the run):
  CONFIRMED iff corr(C, D) >= 0.9 over the sample AND the mean of
            |C - D| / N is <= 0.10 and does not grow across octaves.
  REFUTED   iff corr(C, D) < 0.5: then C(N) is not the deficit and the
            two corpora are about different objects.
"""
import numpy as np
import math


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
    lam = np.zeros(X + 1)
    for p in primes:
        q = int(p)
        lp = math.log(int(p))
        while q <= X:
            lam[q] = lp
            q *= int(p)
    return mu, lam, primes, spf


def singular(N, spf):
    C2 = 0.6601618158468696          # twin-prime constant
    S = 2 * C2
    n = N
    while n > 1:
        p = int(spf[n])
        if p > 2:
            S *= (p - 1) / (p - 2)
        while n % p == 0:
            n //= p
    return S


def main():
    X = 2_000_000
    print(f"sieving to {X} ...", flush=True)
    mu, lam, primes, spf = sieve(X)

    groups = [(120_000, 80), (240_000, 80), (480_000, 80),
              (960_000, 80), (1_900_000, 80)]
    print(f"\n{'N0':>10} {'#N':>4} {'corr(C,D)':>10} {'mean|C-D|/N':>12} "
          f"{'mean C/N':>10} {'sd C/N':>9} {'mean D/N':>10}")
    allC, allD = [], []
    rows = []
    for N0, cnt in groups:
        Cs, Ds = [], []
        for t in range(cnt):
            N = N0 + 2 * t
            if N % 2:
                N += 1
            idx = np.arange(1, N)
            lamv = lam[1:N]
            C = float(np.dot(lamv, mu[N - idx].astype(np.float64)))
            r = float(np.dot(lamv, lam[N - idx]))
            S = singular(N, spf)
            Cs.append(C / N)
            Ds.append((N - r / S) / N)
        Cs = np.array(Cs); Ds = np.array(Ds)
        cc = float(np.corrcoef(Cs, Ds)[0, 1])
        md = float(np.mean(np.abs(Cs - Ds)))
        rows.append((N0, cc, md))
        allC.append(Cs); allD.append(Ds)
        print(f"{N0:>10} {cnt:>4} {cc:>10.4f} {md:>12.5f} "
              f"{Cs.mean():>10.5f} {Cs.std():>9.5f} {Ds.mean():>10.5f}")

    C = np.concatenate(allC); D = np.concatenate(allD)
    cc = float(np.corrcoef(C, D)[0, 1])
    md = float(np.mean(np.abs(C - D)))
    print(f"\npooled: corr = {cc:.4f}   mean|C-D|/N = {md:.5f}   "
          f"({C.size} values of N)")
    # scaling law of the wall's own scalar: sd(C)/N ~ N^{-alpha}
    Ns = np.array([g[0] for g in groups], dtype=np.float64)
    sds = np.array([c.std() for c in allC])
    diffs = np.array([r[2] for r in rows])
    aC = np.polyfit(np.log(Ns), np.log(sds), 1)[0]
    aE = np.polyfit(np.log(Ns), np.log(diffs), 1)[0]
    print(f"\nscaling (5 octave groups):")
    print(f"  sd(C/N)   ~ N^({aC:+.3f})   =>  |C(N)| ~ N^({1+aC:.3f})")
    print(f"  |C-D|/N   ~ N^({aE:+.3f})   =>  |E_3/S| ~ N^({1+aE:.3f})")
    print("  (the chain needs only C(N) = o(N), i.e. exponent < 1)")

    grows = rows[-1][2] > rows[0][2]
    print("verdict:",
          "CONFIRMED -- C(N) is the Goldbach deficit; the two corpora "
          "are about the same object"
          if (cc >= 0.9 and md <= 0.10 and not grows) else
          ("REFUTED -- different objects" if cc < 0.5 else
           "PARTIAL -- see numbers against the pre-registered rule"))
    print("DONE")


if __name__ == "__main__":
    main()
