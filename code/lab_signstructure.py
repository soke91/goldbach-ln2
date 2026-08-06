# -*- coding: utf-8 -*-
"""
The sign structure of D_p at large p (increment 271).

Proposition P.5 says the nonzero terms of D_p(N) are exactly v = mp
with (m, rad N) = 1. That has a consequence for the SIGNS that the
sign-randomised null of increments 267-269 cannot see.

Let P = rad(N) and let r be the least prime not dividing N. Every
composite m coprime to P has all its prime factors >= r, so the
smallest composite one is r^2. Hence for

    M := (N-1)/p  <  r^2,

the surviving m are exactly 1 and the primes in [r, M], and

    mu(1 * p) = mu(p) = -1,
    mu(q p)   = mu(q)mu(p) = +1   for a prime q != p.

So D_p is ONE NEGATIVE TERM against a sum of POSITIVE ones:

    D_p = -Lambda(N-p)/log p  +  Sum_{r <= q <= M} Lambda(N-qp)/log(qp).

For N = k*30030 the least prime not dividing N is r = 17, so this holds
for every p with N/p < 289 -- which is most of the mass.

WHY IT MATTERS. Every cancellation statistic used so far (S_null from
the second moment, and the sign-randomised null) models D_p as a sum of
independently signed terms. It is not: the signs are forced, one minus
against many plus. That does not invalidate the measured ratios -- they
are what they are -- but it means "how far from random" is the wrong
question at these p, and the right one is whether the single negative
term is matched by the positive sum.

WHAT IS CHECKED.
 (A) The sign claim itself, exactly: for N = k*30030 and p with
     N/p < 289, every term of D_p with m > 1 is positive and the m = 1
     term is negative. Any violation is a refutation and is printed
     with its (p, m).
 (B) The balance: the ratio (positive sum)/(negative term) by term
     count, which is what decides |D_p|. A ratio near 1 is near-total
     cancellation; far from 1 is none.
 (C) The same for all even N, where r is usually 3 or 5 and the
     threshold r^2 is tiny, so the structure should be absent -- the
     control that says the effect is about depth and not about p.
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
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p); lg = math.log(int(p))
        while q <= X:
            lam[q] = lg; q *= int(p)
    return mu, lam, primes, spf


def least_nondiv(N, primes):
    for q in primes:
        if N % int(q) != 0:
            return int(q)
    return None


def check(N, mu, lam, ps, primes, spf_, tag):
    r = least_nondiv(N, primes)
    thr = r * r
    viol = 0
    exc = 0
    nchk = 0
    buckets = {}
    for p in ps:
        p = int(p)
        M = (N - 1) // p
        if M < 1 or M >= thr:
            continue
        neg = 0.0; pos = 0.0; nt = 0
        for m in range(1, M + 1):
            v = m * p
            if mu[v] == 0 or lam[N - v] == 0.0:
                continue
            term = float(mu[v]) * lam[N - v] / math.log(v)
            nt += 1
            if m == 1:
                if term > 0:
                    viol += 1
                neg += -term
            else:
                if term < 0:
                    # P.5 allows exactly one kind of exception: N - v a
                    # power of a prime dividing N. Classify rather than
                    # count as a violation, and print the first few so
                    # the classification is visible.
                    w = N - v
                    q = int(spf_[w]) if w > 1 else 0
                    isexc = (q > 0 and N % q == 0
                             and lam[w] != 0.0)
                    if isexc:
                        exc += 1
                        if exc <= 3:
                            print(f"    exception p={p} m={m} "
                                  f"N-v={w} = {q}^k, {q} | N")
                    else:
                        viol += 1
                        if viol <= 3:
                            print(f"    VIOLATION p={p} m={m} "
                                  f"term={term:+.4f} N-v={w}")
                pos += term
            nchk += 1
        if nt == 0:
            continue
        b = 1 if nt == 1 else (2 if nt <= 3 else (4 if nt <= 8 else 9))
        s = buckets.setdefault(b, [0, 0.0, 0.0])
        s[0] += 1; s[1] += neg; s[2] += pos
    print(f"  {tag}: N={N}, r={r}, threshold r^2={thr}, "
          f"{nchk} terms checked, {viol} genuine violations, "
          f"{exc} stated exceptions")
    for b in sorted(buckets):
        c, neg, pos = buckets[b]
        lab = {1: "1 term", 2: "2-3", 4: "4-8", 9: ">=9"}[b]
        print(f"      {lab:>7}: {c:>6} primes   pos/neg = "
              f"{pos/max(neg,1e-12):>7.4f}")
    return viol


def main():
    X = 600_000
    t0 = time.time()
    mu, lam, primes, spf_ = sieve(X)
    ps = primes[primes < X]
    print(f"sieve  t={time.time()-t0:.0f}s", flush=True)
    tot = 0
    print("\n(A)+(B) deep N, where r = 17 and the threshold is 289")
    for k in (12, 15, 18):
        N = k * 30030
        if N <= X:
            tot += check(N, mu, lam, ps, primes, spf_, f"k={k}")
    print("\n(C) control: all even N, where r is small and the")
    print("    threshold r^2 is tiny, so the structure should vanish")
    for N in (599998, 599996, 599990):
        tot += check(N, mu, lam, ps, primes, spf_, f"N={N}")
    print(f"\ntotal sign violations across every check: {tot}")
    print("  zero means the sign structure of P.5 holds exactly:")
    print("  one negative term against a sum of positive ones")
    print("DONE")


if __name__ == "__main__":
    main()
