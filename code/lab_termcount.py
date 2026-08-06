# -*- coding: utf-8 -*-
"""
P.4's demand by the number of terms per prime (increment 269).

Increment 268 found that for p > N/2 the sum D_p has a single term, so
no per-p cancellation statistic can say anything there -- the same fact
as the ceiling of TRANSFORM_P section 5. The natural generalisation has
not been measured: D_p has exactly

    j(p) = # { v < N : p | v }  =  floor((N-1)/p)

terms, so j is the right coordinate for this statistic and the demand
splits by it. Cancellation can buy at most a factor sqrt(j) at each p,
and nothing at all at j = 1.

WHAT IS COMPUTED, over all p < N and grouped by j:

    M_j = Sum_{p : j(p) = j} log p M_p          trivial mass
    S_j = Sum_{p : j(p) = j} log p |D_p|        P.4's demand
    Z_j = the same on sign-randomised data      what cancellation buys

and their cumulative sums from j = 1 upward, so one can read off how
much of the demand sits where cancellation is nearly powerless.

THE POINT. Section 5's ceiling used only j = 1, giving
Sum_p log p |D_p| >= Sum_{N/2<p<N} Lambda(N-p) ~ S(N) N / (2 log N).
If a large share of the mass sits at small j, the ceiling extends: at
j terms the best possible is a factor ~1/sqrt(j) of the mass, so

    Sum_p log p |D_p|  >~  Sum_j M_j / sqrt(j)

is a floor that no per-p cancellation can beat. Whether that floor is
o(N) is then a computable question about the distribution of M_j, and
it is the sharpest form of the ceiling available.

NULLS. Z_j is the sign-randomised level for that j-class, drawn from
the data, so "what cancellation buys" is measured and not modelled
(hazard 4). Classes with fewer than 3 primes are printed with that fact
(hazard 6). Both a deep-N group and all even N are run, since increment
268 showed the two behave differently at small p.
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
    del spf
    return mu, lam, primes


JS = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64, 128, 512, 4096, 10 ** 9]


def one(N, mu, lam, ps, lp, rng):
    v = np.arange(1, N, dtype=np.int64)
    base = lam[N - v] / np.log(np.maximum(v, 2).astype(np.float64))
    mv = mu[1:N].astype(np.float64)
    t = mv * base; t[0] = 0.0
    a = np.abs(mv) * base; a[0] = 0.0
    sup = a > 0
    ts = np.zeros_like(t)
    ts[sup] = a[sup] * rng.choice(np.array([-1.0, 1.0]),
                                  size=int(sup.sum()))
    nb = len(JS) - 1
    M = np.zeros(nb); S = np.zeros(nb); Z = np.zeros(nb)
    cnt = np.zeros(nb, dtype=np.int64)
    for i in range(len(ps)):
        p = int(ps[i])
        j = (N - 1) // p
        b = int(np.searchsorted(JS, j, side='right')) - 1
        if b < 0 or b >= nb:
            continue
        sl = slice(p - 1, None, p)
        M[b] += lp[i] * float(a[sl].sum())
        S[b] += lp[i] * abs(float(t[sl].sum()))
        Z[b] += lp[i] * abs(float(ts[sl].sum()))
        cnt[b] += 1
    return M, S, Z, cnt


def main():
    X = 1_200_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    ps = primes[primes < X]
    lp = np.log(ps.astype(np.float64))
    rng = np.random.default_rng(20260806)
    print(f"sieve, {len(ps)} primes  t={time.time()-t0:.0f}s", flush=True)

    lo = X // 2
    groups = [("deep k*30030",
               [k * 30030 for k in range(lo // 30030 + 1,
                                         X // 30030 + 1)][:5]),
              ("all even N",
               list(range(lo + 2, X + 1, 2 * ((X - lo) // 10)))[:5])]
    for name, Ns in groups:
        M = S = Z = C = None
        for N in Ns:
            m, s, z, c = one(int(N), mu, lam, ps, lp, rng)
            M = m if M is None else M + m
            S = s if S is None else S + s
            Z = z if Z is None else Z + z
            C = c if C is None else C + c
        tot = M.sum()
        print(f"\n=== {name}   n = {len(Ns)}   "
              f"t={time.time()-t0:.0f}s")
        print(f"{'terms j':>14} {'#p':>8} {'mass':>8} {'cum mass':>9} "
              f"{'S/M':>7} {'Z/M':>7} {'S/Z':>7} {'cum S/M':>8}")
        cm = cs = 0.0
        for b in range(len(JS) - 1):
            if C[b] == 0:
                continue
            cm += M[b]; cs += S[b]
            tag = f"{C[b]}" + ("*" if C[b] < 3 else "")
            hi = JS[b + 1] - 1
            lab = (f"{JS[b]}" if hi == JS[b]
                   else f"{JS[b]}-{hi if hi < 10**8 else 'inf'}")
            print(f"{lab:>14} {tag:>8} {M[b]/tot:>8.4f} "
                  f"{cm/tot:>9.4f} {S[b]/M[b]:>7.4f} {Z[b]/M[b]:>7.4f} "
                  f"{S[b]/max(Z[b],1e-12):>7.3f} {cs/tot:>8.4f}")
        floor = sum(M[b] / math.sqrt(0.5 * (JS[b] + min(JS[b+1]-1, 4096)))
                    for b in range(len(JS) - 1) if C[b] > 0)
        print(f"  totals: S/M = {S.sum()/tot:.4f}, "
              f"Z/M = {Z.sum()/tot:.4f}, S/Z = {S.sum()/Z.sum():.4f}")
        print(f"  crude 1/sqrt(j) floor on S/M = {floor/tot:.4f}"
              f"   (what perfect per-p cancellation could give)")
    print("\n  * marks a j-class with fewer than 3 primes")
    print("DONE")


if __name__ == "__main__":
    main()
