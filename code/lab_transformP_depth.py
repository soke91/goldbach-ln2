# -*- coding: utf-8 -*-
"""
Does transform P's demand grow with depth until it fails?
(increment 264)

Increment 263 found that P.4's demand S_abs(N) = Sum_p log p |D_p(N)|
is three times larger at deep N (S_abs/N = 0.7773 for N = k*30030) than
at shallow N (0.2576 for N = 2q), and that every earlier measurement
averaged over all even N and so reported the typical case rather than
the hard one. P.4 is a statement about every N, so the question is
whether the demand keeps growing with depth.

IT MUST BE READ AGAINST THE TRIVIAL BOUND, NOT AGAINST N. The trivial
bound for C(N) is Sum_{v>=2} mu^2(v) Lambda(N-v), which is itself
larger for N with many small prime factors, so S_abs/N conflates two
effects. The scale-free statistic is

    R(N) = S_abs(N) / Sum_p log p M_p(N),

which is 1 when no cancellation happens at all and is what increment
225 reported as 0.38 averaged over all N. P.4 needs R(N) -> 0.

THE LADDER. Cores 2, 6, 30, 210, 2310, 30030 and 510510 -- the
primorials -- give N with 0, 1, 2, 3, 4, 5, 6 small odd primes forced
into rad(N), at matched size. If R climbs toward 1 as the core grows,
P.4 fails for deep N and transform P's margin is an artefact of
averaging; if R levels off below 1, the margin survives and only its
value was misreported.

NULLS AND CRITERIA, on the same line.
 * R = 1 is exactly no cancellation, the failure point for P.4. It is
   printed as the reference on every row.
 * The random-sign null for each row is also computed from the data:
   R_null = Sum_p log p sqrt(V_p) / Sum_p log p M_p with
   V_p the exact second moment, so a row can be read as "how far from
   square-root" without a size heuristic (hazard 4).
 * CRITERION: transform P's margin survives depth iff R stays clearly
   below 1 and its gap to R_null does not close. It fails iff R climbs
   to within noise of 1 at any depth reachable here.
 * Every group is reported with its n, and any group with fewer than
   6 values is printed with that fact rather than dropped (hazard 6).
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


def row(N, mu, lam, ps, lp):
    v = np.arange(1, N)
    base = lam[N - v] / np.log(np.maximum(v, 2).astype(np.float64))
    mv = mu[1:N].astype(np.float64)
    t = mv * base; t[0] = 0.0
    a = np.abs(mv) * base; a[0] = 0.0
    q = t * t
    S = M = V = 0.0
    for i, p in enumerate(ps):
        p = int(p)
        sl = slice(p - 1, None, p)
        S += lp[i] * abs(t[sl].sum())
        M += lp[i] * a[sl].sum()
        V += lp[i] * math.sqrt(q[sl].sum())
    return S, M, V


def main():
    X = 1_600_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    ps = primes[primes < X]
    lp = np.log(ps.astype(np.float64))
    print(f"sieve, {len(ps)} primes  t={time.time()-t0:.0f}s", flush=True)

    cores = [(2, "2"), (6, "2*3"), (30, "2*3*5"), (210, "2*3*5*7"),
             (2310, "..*11"), (30030, "..*13"), (510510, "..*17")]
    lo = X // 2
    print(f"\n{'core':>12} {'#odd p':>7} {'n':>4} {'S/N':>8} "
          f"{'R = S/M':>9} {'R_null':>8} {'R/R_null':>9} {'M/N':>7}")
    for core, lab in cores:
        Ns = [k * core for k in range(lo // core + 1, X // core + 1)]
        if len(Ns) > 24:
            Ns = Ns[:: len(Ns) // 24][:24]
        if not Ns:
            print(f"{lab:>12} {'-':>7} {0:>4}   (no N in range, "
                  f"reported not dropped)")
            continue
        S = M = V = 0.0; sn = 0.0
        for N in Ns:
            s, m, v = row(int(N), mu, lam, ps, lp)
            S += s; M += m; V += v; sn += s / N
        nodd = len([q for q in (3, 5, 7, 11, 13, 17) if core % q == 0])
        tag = f"{len(Ns)}" + ("*" if len(Ns) < 6 else "")
        print(f"{lab:>12} {nodd:>7} {tag:>4} {sn/len(Ns):>8.4f} "
              f"{S/M:>9.4f} {V/M:>8.4f} {(S/M)/(V/M):>9.4f} "
              f"{M/sum(Ns):>7.4f}   t={time.time()-t0:.0f}s", flush=True)
    print("    * marks a group with fewer than 6 values")
    print("    R = 1 is exactly no cancellation, i.e. P.4 failing;")
    print("    R_null is the exact-second-moment random-sign level")
    print("DONE")


if __name__ == "__main__":
    main()
