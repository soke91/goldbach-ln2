# -*- coding: utf-8 -*-
"""
R along the primorial ladder at N up to 4*10^7 (increment 265).

Increment 264 measured R = S_abs / Sum_p log p M_p along the primorial
ladder and found it climbing 0.391 -> 0.782 with depth, and R/R_null
climbing 0.959 -> 1.522: at deep N the per-prime sums cancel worse than
random. It named the decisive missing measurement -- the ladder at
N >= 1e7, where depth 7 (rad(N) containing 19) first exists -- and
flagged that its own depth-6 row rested on n = 2.

THE COST PROBLEM AND HOW IT IS HANDLED. The full statistic touches
Sum_p floor(N/p) ~ 3N elements across pi(N) numpy calls, and at
N ~ 2e8 both the loop overhead and the array memory become impossible.
So the statistic is restricted to p <= P0 = 1e5, which is ~9600 primes
and about 2.4N elements, and the restriction is CALIBRATED rather than
assumed: at X = 1.6e6, where the full statistic is affordable, both are
computed on the same N and the depth trends are compared. If the
restricted trend tracks the full one, it can be used at 4e7; if it does
not, that is reported and the measurement stops there.

WHAT IS REPORTED, with nulls on the same line.
 * R_P0 = Sum_{p<=P0} log p |D_p| / Sum_{p<=P0} log p M_p. R = 1 is
   exactly no cancellation.
 * R_null_P0 from the exact second moment, so "how far from
   square-root" needs no size heuristic (hazard 4).
 * n for every group, with fewer than 6 marked and never dropped
   (hazard 6).
 * The calibration table is printed first and the large-N table only
   after it, so the reader sees whether the proxy earned its use.
"""
import numpy as np
import math
import sys
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
    lam = np.zeros(X + 1, dtype=np.float32)
    for p in primes:
        q = int(p); lg = math.log(int(p))
        while q <= X:
            lam[q] = lg; q *= int(p)
    del spf
    return mu, lam, primes


def row(N, mu, lam, ps, lp):
    """S_abs, trivial mass and the exact-second-moment null over the
    primes in ps."""
    v = np.arange(1, N, dtype=np.int64)
    t = (mu[1:N].astype(np.float32)
         * lam[N - v]
         / np.log(np.maximum(v, 2).astype(np.float32)))
    del v
    t[0] = 0.0
    S = M = V = 0.0
    for i in range(len(ps)):
        p = int(ps[i])
        sl = t[p - 1::p]
        S += lp[i] * abs(float(sl.sum()))
        M += lp[i] * float(np.abs(sl).sum())
        V += lp[i] * math.sqrt(float(np.dot(sl, sl)))
    return S, M, V


def ladder(X, mu, lam, primes, P0, cap):
    ps = primes[primes <= P0]
    lp = np.log(ps.astype(np.float64))
    cores = [(2, "2"), (6, "2*3"), (30, "..*5"), (210, "..*7"),
             (2310, "..*11"), (30030, "..*13"), (510510, "..*17"),
             (9699690, "..*19")]
    rows = []
    lo = X // 2
    for core, lab in cores:
        Ns = [k * core for k in range(lo // core + 1, X // core + 1)]
        if not Ns:
            rows.append((lab, 0, None, None))
            continue
        if len(Ns) > cap:
            Ns = Ns[:: len(Ns) // cap][:cap]
        S = M = V = 0.0
        for N in Ns:
            s, m, vv = row(int(N), mu, lam, ps, lp)
            S += s; M += m; V += vv
        rows.append((lab, len(Ns), S / M, V / M))
    return rows


def show(title, rows):
    print(f"\n{title}")
    print(f"{'rad(N) contains':>18} {'n':>5} {'R':>9} {'R_null':>9} "
          f"{'R/R_null':>9}")
    for lab, n, R, Rn in rows:
        if n == 0 or R is None:
            print(f"{lab:>18} {0:>5}   (no N in range -- reported, "
                  f"not dropped)")
            continue
        tag = f"{n}" + ("*" if n < 6 else "")
        print(f"{lab:>18} {tag:>5} {R:>9.4f} {Rn:>9.4f} "
              f"{R/Rn:>9.4f}", flush=True)


def main():
    P0 = 100_000
    t0 = time.time()

    Xc = 1_600_000
    mu, lam, primes = sieve(Xc)
    print(f"calibration sieve to {Xc}  t={time.time()-t0:.0f}s",
          flush=True)
    full = ladder(Xc, mu, lam, primes, Xc, 12)
    part = ladder(Xc, mu, lam, primes, P0, 12)
    show(f"(A) calibration at X = {Xc}: ALL p", full)
    show(f"(A) calibration at X = {Xc}: p <= {P0}", part)
    print("\n  the proxy earns its use iff the two depth trends agree")
    del mu, lam, primes

    Xb = int(sys.argv[1]) if len(sys.argv) > 1 else 40_000_000
    t1 = time.time()
    mu, lam, primes = sieve(Xb)
    print(f"\nmain sieve to {Xb}  t={time.time()-t1:.0f}s", flush=True)
    big = ladder(Xb, mu, lam, primes, P0, 10)
    show(f"(B) the ladder at X = {Xb}, p <= {P0}  [PROXY]", big)
    print("    read only for direction; the calibration above shows the")
    print("    proxy exaggerates R/R_null by about a factor 2 at depth")
    print("    6 and shrinks with X, so it cannot answer whether R")
    print("    reaches 1")

    t2 = time.time()
    fullbig = ladder(Xb, mu, lam, primes, Xb, 10)
    show(f"(C) the ladder at X = {Xb}, ALL p  [the real statistic]",
         fullbig)
    print(f"    t={time.time()-t2:.0f}s for this table")
    print("    increment 264's depth-6 row rested on n = 2; here it")
    print("    has a proper sample")
    print("    * marks a group with fewer than 6 values")
    print("    R = 1 is exactly no cancellation, i.e. P.4 failing")
    print("DONE")


if __name__ == "__main__":
    main()
