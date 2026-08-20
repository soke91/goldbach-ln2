# -*- coding: utf-8 -*-
"""
a2_identities.py  --  pass4, blind mathematical re-verification.

PRE-REGISTRATION (fixed before the run).

WHAT IS MEASURED, one block per claim, each re-implemented here from the
statement in the paper and not from the papers' own scripts:

 (I1) P1 Lemma 11 (density).  For every squarefree m < 400, the exact
      rational residual  sum_{g|m} mu(g)/(phi(m/g) g phi(g)) - 1/m.
      FALSIFIED by any nonzero residual.
 (I2) P1 eq:(30).  A(N)*Gtilde(x) with
      Gtilde(x)=sum_{m<=x} mu(m)lambda(m)1_{(m,N)=1} log m / m,
      lambda(m)=prod_{p|m}(1-1/(p(p-1)))^{-1}, at N=4e6, x=4e6,
      against -S(N).  Paper prints -1.760250 against S(N)=1.760432.
      FALSIFIED if the recomputation differs from -1.760250 in the
      sixth decimal, or if the sign is not negative.
 (I3) P3 Lemma 5 (second-moment identity), at X = 400, 800, 1600:
      sum_N Chat(N)^2  against  sum_{|h|<X} M(h)P(h), exactly.
      FALSIFIED by any relative discrepancy above 1e-12.
 (I4) P3 Note 8 (the truncation is load-bearing): the same right side
      against the UNtruncated left side restricted to N<=X, whose ratio
      the paper says is "near 1.57 at X from 800 to 3200" and "does not
      tend to 1".  FALSIFIED if the ratio is within 0.01 of 1, or if it
      moves monotonically toward 1 across the four X.
 (I5) P5 Proposition 6: sum_{k>=1} sum_{m: mk<=N-1} mu(m)mu(N-mk)
      against mu(N-1), at N = 200, 1000, 5000, 20000.
      FALSIFIED by any mismatch.
 (I6) P3 Measurement 12's "about 430 standard deviations": the actual
      standard deviation of T_eps(x) over 200 independent sign draws at
      x=2e6, against |T(x)|, and against the max-of-20 figure 0.000935
      that the results file prints.
      FALSIFIED if |T|/sd lands within 10% of 430, i.e. if the printed
      figure is a standard-deviation count after all.

NULL.  I1-I5 are identities with no sign input; a sign control would
change both sides alike and settle nothing, so none is run.  I6 IS a
null measurement: the control is the i.i.d. sign ensemble on
supp(mu^2), 200 draws, seed 20260808.  Its threshold (the 10% band)
is set from the printed claim and not from the measured effect.
"""
import numpy as np
import math
from fractions import Fraction


def sieve_all(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    primes = np.nonzero(s)[0]
    mu = np.ones(n + 1, dtype=np.int64)
    mu[0] = 0
    phi = np.arange(n + 1, dtype=np.int64)
    lam = np.zeros(n + 1)
    for p in primes:
        p = int(p)
        lp = math.log(p)
        mu[p::p] *= -1
        phi[p::p] -= phi[p::p] // p
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            lam[q] = lp
            q *= p
    return primes, mu, phi, lam


print(__doc__.strip())
print()
LIM = 4000000
primes, mu, phi, lam = sieve_all(LIM)

# ---------------------------------------------------------------- I1
print("I1  P1 Lemma 11, exact rationals, every squarefree m < 400")


def phi_int(n):
    r = n
    m = n
    for p in primes:
        p = int(p)
        if p * p > m:
            break
        if m % p == 0:
            r -= r // p
            while m % p == 0:
                m //= p
    if m > 1:
        r -= r // m
    return r


bad = 0
tested = 0
for m in range(1, 400):
    if mu[m] == 0:
        continue
    tested += 1
    tot = Fraction(0)
    divs = [d for d in range(1, m + 1) if m % d == 0]
    for g in divs:
        tot += Fraction(int(mu[g]), phi_int(m // g) * g * phi_int(g))
    if tot != Fraction(1, m):
        bad += 1
print("    squarefree m tested: %d   mismatches: %d   -> %s"
      % (tested, bad, "HOLDS" if bad == 0 else "FAILS"))

# ---------------------------------------------------------------- I2
print()
print("I2  P1 eq:(30), A(N)*Gtilde(x) against -S(N)")
N = 4000000
pf = []
m = N
for p in primes:
    p = int(p)
    if p * p > m:
        break
    if m % p == 0:
        pf.append(p)
        while m % p == 0:
            m //= p
if m > 1:
    pf.append(m)
A = 1.0
S = 2.0
for q in primes:
    q = int(q)
    if q > 2:
        S *= (1.0 - 1.0 / float((q - 1) ** 2))
    if N % q != 0:
        A *= (1.0 - 1.0 / float(q * (q - 1)))
for p in pf:
    if p > 2:
        S *= (1.0 + 1.0 / float(p - 2))
lamb = np.ones(LIM + 1)
for p in primes:
    p = int(p)
    lamb[p::p] *= 1.0 / (1.0 - 1.0 / (p * (p - 1)))
ms = np.arange(1, LIM + 1)
copr = np.ones(LIM + 1, dtype=bool)
for p in pf:
    copr[p::p] = False
f = mu[1:LIM + 1] * lamb[1:LIM + 1] * copr[1:LIM + 1]
G = float(np.sum(f * np.log(ms) / ms))
print("    N=%d  prime support %s  S(N)=%.6f  A(N)=%.6f" % (N, pf, S, A))
print("    A(N)*Gtilde(4e6) = %+.6f   -S(N) = %+.6f   paper prints -1.760250"
      % (A * G, -S))
print("    |A*Gt + S| = %.6f   sign negative: %s" % (abs(A * G + S), A * G < 0))

# ---------------------------------------------------------------- I3
print()
print("I3  P3 Lemma 5, truncated convolution, exact")


def MPsides(X):
    L = lam[:X + 1].copy()
    Mu = mu[:X + 1].astype(np.float64)
    C = np.convolve(L, Mu)
    lhs = float(np.dot(C, C))
    rhs = 0.0
    for h in range(-(X - 1), X):
        a = abs(h)
        Mh = float(np.dot(Mu[1:X + 1 - a], Mu[1 + a:X + 1]))
        Ph = float(np.dot(L[1:X + 1 - a], L[1 + a:X + 1]))
        rhs += Mh * Ph
    return lhs, rhs, C


for X in (400, 800, 1600):
    l, r, _ = MPsides(X)
    print("    X=%-5d lhs=%.10e  rhs=%.10e  rel=%.2e"
          % (X, l, r, abs(l - r) / abs(l)))

# ---------------------------------------------------------------- I4
print()
print("I4  P3 Note 8, the untruncated left side on N<=X against the same box")
for X in (800, 1600, 3200, 6400):
    _, rhs, _ = MPsides(X)
    L = lam[:X + 1].copy()
    Mu = mu[:X + 1].astype(np.float64)
    Cfull = np.convolve(L, Mu)
    lhs = float(np.dot(Cfull[:X + 1], Cfull[:X + 1]))
    print("    X=%-6d simplex=%.6e  box=%.6e  box/simplex=%.4f"
          % (X, lhs, rhs, rhs / lhs))

# ---------------------------------------------------------------- I5
print()
print("I5  P5 Proposition 6, the untruncated dilate double sum")
for Nv in (200, 1000, 5000, 20000):
    tot = 0
    for k in range(1, Nv):
        mmax = (Nv - 1) // k
        if mmax == 0:
            break
        mm = np.arange(1, mmax + 1)
        tot += int(np.dot(mu[mm], mu[Nv - mm * k]))
    print("    N=%-6d double sum = %-4d   mu(N-1) = %-4d   %s"
          % (Nv, tot, int(mu[Nv - 1]),
             "ok" if tot == int(mu[Nv - 1]) else "MISMATCH"))

# ---------------------------------------------------------------- I6
print()
print("I6  P3 Measurement 12, what 'about 430 standard deviations' counts")
x = 2000000
v = np.arange(1, x + 1)
odd_sqf = (mu[1:x + 1] != 0) & (v % 2 == 1)
idx = v[odd_sqf]
T = float(np.sum(mu[idx] * mu[2 * idx])) / x
rng = np.random.default_rng(20260808)
supp = np.nonzero(mu[1:2 * x + 1] != 0)[0] + 1
pos = np.zeros(2 * x + 1, dtype=np.int64)
pos[supp] = np.arange(len(supp))
draws = 200
vals = np.empty(draws)
i1 = pos[idx]
i2 = pos[2 * idx]
for d in range(draws):
    eps = rng.integers(0, 2, size=len(supp)).astype(np.int8) * 2 - 1
    vals[d] = float(np.sum(eps[i1].astype(np.int64) * eps[i2])) / x
sd = float(vals.std(ddof=1))
mx = float(np.abs(vals).max())
print("    T(x)                = %+.6f   (paper -0.405295; -4/pi^2 = %+.6f)"
      % (T, -4 / math.pi ** 2))
print("    coin sd (200 draws) = %.6e   max|T_eps| over 200 = %.6e" % (sd, mx))
print("    |T|/sd = %.1f     |T|/0.000935 (the max-of-20 in the results file) = %.1f"
      % (abs(T) / sd, abs(T) / 0.000935))
print("    is |T|/sd within 10%% of 430?  %s"
      % ("yes" if abs(abs(T) / sd - 430) <= 43 else "no"))
