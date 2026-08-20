# -*- coding: utf-8 -*-
"""
a8_delta_trend.py  --  pass4 addendum.  Follow-up to a7.

PRE-REGISTRATION (fixed before the run).

WHAT IS MEASURED.  Only the correction term itself, over a longer
range of N than a7 could reach (a7 also recomputes S_2 by divisor
enumeration, which is quadratic; Delta alone is linear).

  Delta(N, theta) = sum_{2<=m<=alpha} mu(m) log m
                       sum_{k < (N-1)/alpha, (k,m)=1} mu^2(k) Lambda(N-mk),
  alpha = N**theta.

  Reported as Delta/N at N = 2e4, 8e4, 3.2e5, 1.28e6, 5.12e6 and
  theta in {0.30, 0.40, 0.45}, together with the trivial bound
  (log N)^2 that P1 section 7 quotes for it.

WHAT WOULD FALSIFY WHAT.
  P1 section 7 says only that Delta's trivial bound N (log N)^2
  "exceeds the target O(N (log N)^{-A}), so Delta is not negligible
  and needs its own treatment".  It does NOT claim Delta is large
  asymptotically -- its repair argument predicts a main term
  O(N e^{-c sqrt(log N)}) plus a Bombieri-Vinogradov error.
  So: CONSISTENT with P1 if |Delta|/N drifts downward or stays
  small and bounded; INCONSISTENT with P1's repair (not with the
  defect report) if |Delta|/N grows toward the trivial (log N)^2.

  This is a size measurement, not a detection test; no threshold
  is set from the data.

NULL.  None applies.  Delta is a deterministic finite sum with mu
supplied by arithmetic, not sampled; a sign control would replace
the object under measurement rather than control it.  The relevant
reference is the trivial bound, which is printed alongside.
"""
import numpy as np
import math

def sieve_all(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    primes = np.nonzero(s)[0]
    mu = np.ones(n + 1, dtype=np.int64)
    mu[0] = 0
    lam = np.zeros(n + 1)
    for p in primes:
        p = int(p)
        lp = math.log(p)
        mu[p::p] *= -1
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            lam[q] = lp
            q *= p
    return primes, mu, lam


print(__doc__.strip())
print()
LIM = 5120000
primes, mu, lam = sieve_all(LIM)
mu2 = (mu != 0)

print("  %-9s %-6s %-10s %-8s %-13s %-11s %-12s"
      % ("N", "theta", "alpha", "#m", "Delta", "Delta/N", "(log N)^2"))
print("  " + "-" * 80)
for N in (20000, 80000, 320000, 1280000, 5120000):
    for th in (0.30, 0.40, 0.45):
        alpha = N ** th
        Kp = (N - 1) / alpha
        kmax = int(math.ceil(Kp)) - 1
        while kmax >= Kp:
            kmax -= 1
        D = 0.0
        nm = 0
        for m in range(2, int(math.floor(alpha)) + 1):
            if mu[m] == 0:
                continue
            nm += 1
            top = min(kmax, (N - 1) // m)
            if top < 1:
                continue
            ks = np.arange(1, top + 1)
            keep = mu2[1:top + 1].copy()
            mm = m
            for p in primes:
                p = int(p)
                if p * p > mm:
                    break
                if mm % p == 0:
                    keep[p - 1::p] = False
                    while mm % p == 0:
                        mm //= p
            if mm > 1 and mm <= top:
                keep[mm - 1::mm] = False
            D += mu[m] * math.log(m) * float(np.sum(lam[N - ks[keep] * m]))
        print("  %-9d %-6.2f %-10.1f %-8d %-13.5e %-11.5f %-12.1f"
              % (N, th, alpha, nm, D, D / N, math.log(N) ** 2))
