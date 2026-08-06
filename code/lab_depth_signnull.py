# -*- coding: utf-8 -*-
"""
The depth ladder against a sign-randomised null (increment 267).

Increments 264-266 read the ladder against R_null = Sum_p log p
sqrt(V_p), with V_p the exact second moment of D_p, and called
R/R_null > 1 "worse than random". That null is wrong, and wrong in a
p-dependent way.

For a mean-zero fluctuation with MANY terms the expected absolute value
is E|D_p| = sqrt(2/pi) sqrt(V_p) = 0.798 sqrt(V_p), not sqrt(V_p). But
for large p the sum has FEW terms -- at p > N/2 exactly one -- and there
E|D_p| = sqrt(V_p) exactly. So the correct random-sign level runs from
0.798 sqrt(V_p) at small p to sqrt(V_p) at large p, and since the mass
sits at large p the printed ratio was neither the one nor the other.

THE FIX NEEDS NO FORMULA. Replace mu(v) by random signs on the same
support -- same |mu|, same Lambda weights, same p-structure -- and
recompute S_abs. That is the random-sign level for every p at once,
with the right term counts built in, and it is drawn from the data
rather than assumed (hazard 4).

WHAT IS REPORTED.
 * R = S_abs / Sum_p log p M_p, unchanged, with R = 1 no cancellation.
 * R_sign = the same statistic on sign-randomised data, averaged over
   draws, with the spread across draws printed so the reader sees the
   null's own error bar.
 * R / R_sign, which is the honest "how far from random" and replaces
   the R/R_null column of increments 264-266.
 * The old analytic column is printed beside it so the size of the
   correction is visible rather than silently absorbed.

CRITERION. If R/R_sign is near 1 at every depth, the ladder's climb is
entirely the growth of the trivial mass and there is no anomaly. If it
climbs with depth, the per-prime sums genuinely cancel worse than
random at deep N and the conclusion of increment 264 survives with a
corrected baseline.
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


def row(N, mu, lam, ps, lp, rng, ndraw):
    v = np.arange(1, N, dtype=np.int64)
    base = lam[N - v] / np.log(np.maximum(v, 2).astype(np.float64))
    mv = mu[1:N].astype(np.float64)
    t = mv * base; t[0] = 0.0
    a = np.abs(mv) * base; a[0] = 0.0
    q = t * t
    sup = a > 0
    S = M = V = 0.0
    draws = np.zeros(ndraw)
    ts = np.empty_like(t)
    for d in range(ndraw):
        ts[:] = a
        sg = rng.choice(np.array([-1.0, 1.0]), size=int(sup.sum()))
        ts[sup] = a[sup] * sg
        ts[~sup] = 0.0
        acc = 0.0
        for i in range(len(ps)):
            p = int(ps[i])
            acc += lp[i] * abs(float(ts[p - 1::p].sum()))
        draws[d] = acc
    for i in range(len(ps)):
        p = int(ps[i])
        sl = slice(p - 1, None, p)
        S += lp[i] * abs(float(t[sl].sum()))
        M += lp[i] * float(a[sl].sum())
        V += lp[i] * math.sqrt(float(q[sl].sum()))
    return S, M, V, draws


def main():
    X = 1_600_000
    NDRAW = 3
    t0 = time.time()
    mu, lam, primes = sieve(X)
    ps = primes[primes < X]
    lp = np.log(ps.astype(np.float64))
    rng = np.random.default_rng(20260806)
    print(f"sieve, {len(ps)} primes  t={time.time()-t0:.0f}s", flush=True)

    cores = [(2, "2 (all even N)"), (30, "..*5"), (210, "..*7"),
             (2310, "..*11"), (30030, "..*13"), (510510, "..*17")]
    lo = X // 2
    print(f"\n{'N divisible by':>16} {'n':>4} {'R':>8} {'R_sign':>9} "
          f"{'sd/draw':>8} {'R/R_sign':>9} {'old R/R_null':>13}")
    for core, lab in cores:
        Ns = [k * core for k in range(lo // core + 1, X // core + 1)]
        if not Ns:
            print(f"{lab:>16} {0:>4}   (none in range, reported)")
            continue
        if len(Ns) > 6:
            Ns = Ns[:: len(Ns) // 6][:6]
        S = M = V = 0.0
        D = np.zeros(NDRAW)
        for N in Ns:
            s, m, v, d = row(int(N), mu, lam, ps, lp, rng, NDRAW)
            S += s; M += m; V += v; D += d
        rs = D / M
        tag = f"{len(Ns)}" + ("*" if len(Ns) < 6 else "")
        print(f"{lab:>16} {tag:>4} {S/M:>8.4f} {rs.mean():>9.4f} "
              f"{rs.std():>8.4f} {(S/M)/rs.mean():>9.4f} "
              f"{(S/M)/(V/M):>13.4f}   t={time.time()-t0:.0f}s",
              flush=True)
    print("    rows are NESTED: 'N divisible by 30030' is a subset of")
    print("    'N divisible by 210', and the first row is every even N,")
    print("    not a depth-0 group -- the labels are lower bounds on")
    print("    depth, which increments 264-266 did not say")
    print("    * marks a group with fewer than 6 values")
    print("DONE")


if __name__ == "__main__":
    main()
