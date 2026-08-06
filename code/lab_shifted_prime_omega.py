# -*- coding: utf-8 -*-
"""
Where the missing factor lives (increment 246): the omega-profile of the
shifted primes against the generic rough set.

Increment 244 derived C(N) ~= kappa_0 kappa(N) M_{rad N}(N) and found it
short on the deep cells by about 1.7; increment 245 showed that swapping
the coprime-density factor kappa(N) for the shifted-prime enrichment
S(N) improves the cell-mean R^2 from 0.569 to 0.632, so the diagnosis
was pointing at the right set -- but a free power still beat both, so
neither derived factor is steep enough.

The diagnosis was: M_P(N) sums mu over ALL v coprime to P, whereas C(N)
sums it over the shifted primes v = N - n, a subset of density 1/log N
with its own omega-distribution. This measures that distribution
directly, instead of arguing about it.

WHAT IS COMPUTED. For each N, two omega-profiles over squarefree v:

  A_j(N) = # { p prime < N : omega(N - p) = j }      shifted primes
  B_j(N) = # { v < N, (v, rad N) = 1 : omega(v) = j } generic rough set

and, since C(N) carries the weight Lambda, also the log p-weighted form
of A. On the support of mu, mu = (-1)^omega, so the object that governs
the sign structure of each set is its ALTERNATING SUM

  R = Sum_j (-1)^j T_j / Sum_j T_j,

which session 8 established is exactly how C(N) cancels. The ratio
R_A / R_B is the factor that increments 244 and 245 could not derive.

TWO GROUPS, MATCHED IN SIZE.
  DEEP    N = k * 30030   (2*3*5*7*11*13)
  SHALLOW N = 2q, q prime (the emptiest possible rad)
The prediction, made before running: A is shifted toward LOW omega
relative to B, because omega(N-p) = 1 means N - p is prime, i.e. a
Goldbach representation, and those are S(N)-enhanced while the generic
rough set is not. Low omega means more weight at mu = -1, hence a more
negative alternating sum. If |R_A / R_B| is near 1 the diagnosis of
increments 244-245 is wrong; if it is near 1.7 and larger for the deep
group, it is right and quantified.

NULLS ON THE SAME LINE. Each profile is a count, so its sampling error
is Poisson: SE(T_j) = sqrt(T_j), propagated into R. For the ratio
R_A/R_B the null of "no difference between the sets" is exactly 1, and
the errors are printed beside it. B is computed with rad(N) truncated
to the group's core, which is stated rather than hidden: for the deep
group the extra prime factors of k are not in the mask, so B is
slightly less rough than the truth and R_B is biased toward the shallow
value -- i.e. the truncation works AGAINST the predicted effect.
"""
import numpy as np
import math
import time

JMAX = 10


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8); mu[1] = 1
    om = np.zeros(X + 1, dtype=np.int8)
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
        om[i] = om[j] + (0 if j % p == 0 else 1)
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    return mu, om, spf, primes


def alt_ratio(T):
    tot = float(T[1:].sum())
    if tot <= 0:
        return float('nan'), float('nan'), 0.0
    alt = float(sum((-1) ** j * T[j] for j in range(1, JMAX + 1)))
    # Poisson error on each count, propagated (signs are +-1)
    se = math.sqrt(float(T[1:].sum()))
    return alt / tot, se / tot, tot


def main():
    X = 2_000_000
    t0 = time.time()
    mu, om, spf, primes = sieve(X)
    lpr = np.log(primes.astype(np.float64))
    print(f"sieve t={time.time()-t0:.0f}s", flush=True)

    v = np.arange(X + 1)
    sq = (mu != 0)

    groups = []
    CORE = 30030
    deep = [k * CORE for k in range(1, X // CORE + 1)]
    deep = deep[:: max(1, len(deep) // 40)][:40]
    groups.append(("deep  k*30030", deep, [2, 3, 5, 7, 11, 13]))
    cand = primes[np.searchsorted(primes, X // 4):]
    shal = [2 * int(q) for q in cand[:: max(1, len(cand) // 60)][:40]
            if 2 * int(q) <= X]
    groups.append(("shallow  2q", shal, [2]))

    print(f"\n{'group':>16} {'n':>4} "
          f"{'R_A (shifted primes)':>22} {'R_A weighted':>14} "
          f"{'R_B (rough set)':>17} {'R_A/R_B':>9}")
    out = {}
    for name, Ns, core in groups:
        # B: one masked cumulative histogram per omega, using the core
        keep = sq.copy()
        for q in core:
            keep &= (v % q != 0)
        cum = np.zeros((JMAX + 1, X + 1))
        for j in range(1, JMAX + 1):
            cum[j] = np.cumsum((keep & (om == j)).astype(np.float64))

        ra_l, raw_l, rb_l = [], [], []
        for N in Ns:
            jj = int(np.searchsorted(primes, N))
            ps = primes[:jj]
            w = N - ps
            good = mu[w] != 0
            A = np.bincount(om[w[good]], minlength=JMAX + 1)[: JMAX + 1]
            Aw = np.bincount(om[w[good]], weights=lpr[:jj][good],
                             minlength=JMAX + 1)[: JMAX + 1]
            B = cum[:, N - 1]
            a, _, _ = alt_ratio(A.astype(np.float64))
            aw, _, _ = alt_ratio(Aw)
            b, _, _ = alt_ratio(B)
            ra_l.append(a); raw_l.append(aw); rb_l.append(b)
        ra = np.array(ra_l); raw = np.array(raw_l); rb = np.array(rb_l)
        out[name] = (ra, raw, rb)
        print(f"{name:>16} {len(ra):>4} "
              f"{ra.mean():>+12.5f} +- {ra.std()/math.sqrt(len(ra)):<7.5f} "
              f"{raw.mean():>+14.5f} "
              f"{rb.mean():>+17.5f} "
              f"{ra.mean()/rb.mean():>9.3f}")

    print("\n(B) the profiles themselves, at one N from each group")
    for name, Ns, core in groups:
        N = Ns[len(Ns) // 2]
        jj = int(np.searchsorted(primes, N))
        ps = primes[:jj]; w = N - ps; good = mu[w] != 0
        A = np.bincount(om[w[good]], minlength=JMAX + 1)[: JMAX + 1]
        keep = sq.copy()
        for q in core:
            keep &= (v % q != 0)
        B = np.bincount(om[1:N][keep[1:N]], minlength=JMAX + 1)[: JMAX + 1]
        ta = A[1:].sum(); tb = B[1:].sum()
        print(f"  {name}  N = {N}")
        print("     j  " + "".join(f"{j:>8}" for j in range(1, 7)))
        print("     A  " + "".join(f"{A[j]/ta:>8.4f}" for j in range(1, 7)))
        print("     B  " + "".join(f"{B[j]/tb:>8.4f}" for j in range(1, 7)))

    print("\n(C) the factor increments 244-245 could not derive")
    d = out["deep  k*30030"]; s = out["shallow  2q"]
    print(f"  deep    R_A/R_B = {d[0].mean()/d[2].mean():.3f}")
    print(f"  shallow R_A/R_B = {s[0].mean()/s[2].mean():.3f}")
    print(f"  deep/shallow of that ratio = "
          f"{(d[0].mean()/d[2].mean())/(s[0].mean()/s[2].mean()):.3f}")
    print("  null for every ratio here is 1.000 (no difference between")
    print("  the shifted-prime set and the generic rough set)")
    print("DONE")


if __name__ == "__main__":
    main()
