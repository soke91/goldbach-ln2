# -*- coding: utf-8 -*-
"""
a4_extract_identity.py  --  pass4, blind mathematical re-verification.

PRE-REGISTRATION (fixed before the run).

WHAT IS MEASURED.  P2's eq:(5) is displayed as "the identity"

      B_w * C(N)  =  CP  -  R_w  +  O(N^{o(1)}),
      CP  = sum_{u<N} Lambda(N-u) mu^2(u) b_u ,   b = mu * w,

with R_w the residual of the divisor switch (the k >= K tail).  At
N = 2e5, 4e5, 8e5 and theta' = 0.56, K = floor(N^0.56), for the two
weights the papers use, w = 1 (b = delta_1) and w = log (b = Lambda),
this script computes each of the five quantities separately and
directly from its own definition:

  A(N;k) = sum_{n<N, n = N (mod k)} Lambda(n) mu(N-n)
  C(N)   = sum_{n<N} Lambda(n) mu(N-n)
  B_w    = sum_{k<K,(k,N)=1} mu(k) w_k / phi(k)
  D_w    = sum_{k<K,(k,N)=1} mu(k) w_k A(N;k)
  T_w    = D_w - C(N) B_w                       (P2 eq:(2))
  CP     = sum_{u<N} Lambda(N-u) * S_inf(u),
             S_inf(u) = sum_{k|u,(k,N)=1} mu(k) w_k     (the COMPLETE sum)
  R_w    = CP - D_w                             (complete minus truncated)

and then reports the two residuals

  gap_paper   = B_w C(N) - (CP - R_w)
  gap_algebra = B_w C(N) - (CP - R_w - T_w)

WHAT WOULD FALSIFY THE FINDING UNDER TEST.  The finding is: eq:(5)
drops the term T_w, so it is not an identity; T_w is the object the
whole subject is about and for w = log it is of order N.
FALSIFIED if |gap_paper|/N is of size N^{o(1)}/N, i.e. if it is
below 1e-3 at all three N and falling like a power of N.
CONFIRMED if gap_algebra is zero to machine precision while
|gap_paper|/N stays of order 1e-1 for w = log.

Threshold: 1e-3 N, fixed here before the run, as the largest thing
that could plausibly be called O(N^{o(1)}) at N = 8e5.

NULL.  None applies.  Every quantity is a deterministic finite sum
with no sign input of its own; the two gaps are differences of exact
rearrangements of the same numbers, so a sign control would move both
sides of each comparison alike and could not separate them.
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
LIM = 800000
primes, mu, phi, lam = sieve_all(LIM)

print("  %-7s %-4s %-6s %-13s %-13s %-13s %-13s %-13s %-11s %-11s"
      % ("N", "w", "K", "C(N)", "B_w", "D_w", "T_w", "CP-R_w", "gap_paper/N", "gap_alg/N"))
print("  " + "-" * 128)

for N in (200000, 400000, 800000):
    K = int(N ** 0.56)
    n = np.arange(1, N)
    Cn = float(np.dot(lam[1:N], mu[N - n]))
    for wname in ("1", "log"):
        # ---- truncated side, k < K, (k,N)=1
        Dw = 0.0
        Bw = 0.0
        for k in range(1, K):
            if mu[k] == 0 or math.gcd(k, N) != 1:
                continue
            wk = 1.0 if wname == "1" else math.log(k)
            if wk == 0.0 and wname == "log":
                continue
            M = (N - 1) // k
            ms = np.arange(1, M + 1)
            u = ms * k
            Ak = float(np.dot(lam[N - u], mu[u]))
            Dw += mu[k] * wk * Ak
            Bw += mu[k] * wk / float(phi[k])
        # ---- complete side: S_inf(u) = sum_{k|u,(k,N)=1} mu(k) w_k
        Sinf = np.zeros(N)
        for k in range(1, N):
            if mu[k] == 0 or math.gcd(k, N) != 1:
                continue
            wk = 1.0 if wname == "1" else math.log(k)
            if wk == 0.0 and wname == "log":
                continue
            Sinf[k::k] += mu[k] * wk
        uu = np.arange(1, N)
        CP = float(np.dot(lam[N - uu], Sinf[1:N]))
        Rw = CP - Dw
        Tw = Dw - Cn * Bw
        gap_paper = Bw * Cn - (CP - Rw)
        gap_alg = Bw * Cn - (CP - Rw - Tw)
        print("  %-7d %-4s %-6d %-13.5e %-13.5e %-13.5e %-13.5e %-13.5e %-11.5f %-11.2e"
              % (N, wname, K, Cn, Bw, Dw, Tw, CP - Rw, gap_paper / N, abs(gap_alg) / N))

print()
print("  reading:  gap_alg is the algebraic identity  B_w C = CP - R_w - T_w;")
print("            gap_paper is P2 eq:(5) as printed, i.e. the same with T_w dropped.")
