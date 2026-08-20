# -*- coding: utf-8 -*-
"""
a7_hl18_delta.py  --  pass4 addendum, after the blind was lifted.

PRE-REGISTRATION (fixed before the run).

WHAT IS MEASURED.  Huang-Li arXiv:2005.03811v2, section 3.2, page 16.
The line immediately above their (18) carries the n-DEPENDENT bound
k < (N-n)/alpha; their (18) carries the n-FREE bound k < (N-1)/alpha.
This script computes, exactly, at N in {2000, 4000, 8000, 16000} and
theta in {0.30, 0.40, 0.45}, alpha = N**theta:

  S2   = sum_{n<N} Lambda(n) mu^2(N-n) sum_{d|(N-n), d>alpha} mu(d) log(1/d)
         (Huang-Li's definition of S_2(alpha), section 3, page 12)

  R18  = sum_{k < (N-1)/alpha} mu(k) sum_{n<N, n = N (mod k)}
              Lambda(n) mu(N-n) log(k/(N-n))
         (the right-hand side of their (18) as printed)

  Delta_P1 = sum_{2<=m<=alpha} mu(m) log m
                sum_{k < (N-1)/alpha, (k,m)=1} mu^2(k) Lambda(N-mk)
         (the correction term of P1 section 7, eq (25))

WHAT WOULD FALSIFY WHAT.
  (a) "the two ranges differ": FALSIFIED if S2 == R18 to machine
      precision at every (N, theta).
  (b) "P1's Delta is exactly that difference": FALSIFIED if
      |S2 - (R18 + Delta_P1)| / N exceeds 1e-12 at any (N, theta).
  (c) "the missing term is not negligible": the trivial bound P1
      quotes is Delta << N (log N)^2.  Reported here as Delta/N so
      that its actual size at accessible N is on the record.  No
      threshold is attached to (c); it is a size, not a test.

NULL.  None applies and none would mean anything.  All three
quantities are deterministic finite sums over the same integers with
no sign input of their own; (a) and (b) are exact identities between
rearrangements, so a sign control would move every side alike.
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
LIM = 16000
primes, mu, lam = sieve_all(LIM)
mu2 = (mu != 0).astype(np.float64)

print("  %-7s %-6s %-9s %-14s %-14s %-14s %-12s %-11s %-10s"
      % ("N", "theta", "alpha", "S2", "R18", "Delta_P1",
         "|S2-R18-D|/N", "(R18-S2)/N", "Delta/N"))
print("  " + "-" * 116)

for N in (2000, 4000, 8000, 16000):
    for th in (0.30, 0.40, 0.45):
        alpha = N ** th
        Kp = (N - 1) / alpha

        # ---- S2, straight from the definition
        S2 = 0.0
        for n in range(1, N):
            if lam[n] == 0.0:
                continue
            u = N - n
            if mu[u] == 0:
                continue
            t = 0.0
            for d in range(1, u + 1):
                if u % d == 0 and d > alpha:
                    t += mu[d] * math.log(1.0 / d)
            S2 += lam[n] * t

        # ---- R18, the right-hand side of (18) as printed
        R18 = 0.0
        for k in range(1, int(math.ceil(Kp))):
            if k >= Kp or mu[k] == 0:
                continue
            # n < N with n = N mod k, i.e. N-n = mk
            acc = 0.0
            for m in range(1, (N - 1) // k + 1):
                u = m * k
                n = N - u
                if n < 1 or lam[n] == 0.0 or mu[u] == 0:
                    continue
                acc += lam[n] * mu[u] * math.log(k / float(u))
            R18 += mu[k] * acc

        # ---- Delta as P1 section 7 defines it
        D = 0.0
        for m in range(2, int(math.floor(alpha)) + 1):
            if mu[m] == 0:
                continue
            lm = math.log(m)
            acc = 0.0
            for k in range(1, int(math.ceil(Kp))):
                if k >= Kp or mu[k] == 0 or math.gcd(k, m) != 1:
                    continue
                n = N - m * k
                if n < 1:
                    continue
                acc += lam[n]
            D += mu[m] * lm * acc

        print("  %-7d %-6.2f %-9.2f %-14.6e %-14.6e %-14.6e %-12.2e %-11.5f %-10.5f"
              % (N, th, alpha, S2, R18, D,
                 abs(S2 - R18 - D) / N, (R18 - S2) / N, D / N))
print()
print("  (R18-S2)/N is the size of what (18) silently gained; Delta/N is minus that.")
