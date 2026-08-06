# -*- coding: utf-8 -*-
"""
Does the supply side carry a location mask too? (increment 250)

Increments 239-249 found that the wall C(N) is not mean-zero: it carries
a deterministic term set by which small primes divide N, arising because
n prime forces N - n coprime to rad(N). The campaign's own code states
the opposite assumption for the supply side -- e1_forge_r4.py line 19:
"Under Conjecture L the field {D(k)} is mean-zero and INDEPENDENT" --
and several E1 kill-tests are built on it. That assumption has never
been audited, and the same mechanism applies.

THE MECHANISM, TRANSPOSED. With

    D(k) = Sum_{sqrt(N) < m <= N/k} mu(m) mu(N - mk)

(the campaign's own definition, code/hyp_round.py), take q | k. Then
N - mk = N mod q, so

    q | N - mk   <=>   q | N.

So when q divides k but NOT N, the variable v = N - mk is forced
coprime to q -- for every m at once. That is exactly the forcing that
produced the wall's mask, now indexed by the primes dividing k rather
than N. If the mask mechanism is general, D(k) has a nonzero mean
depending on which small primes divide k and do not divide N.

PREDICTION, made before measuring. Write d(k) = # { q <= 23 : q | k,
q not| N }, the number of small primes that force coprimality. The
prediction is that mean D(k) becomes systematically nonzero as d(k)
grows, in the same direction as the wall's mask -- negative, since
removing small primes from the pool over-represents the primes, where
mu = -1.

NULLS, ON THE SAME LINE.
 * For each group of k, the null is mean D = 0 with SE = sd(D)/sqrt(n),
   computed from the same data.
 * A PERMUTATION CONTROL that cannot be skipped: the same statistic
   after permuting the k-labels, which destroys any dependence on the
   factorisation of k while preserving the distribution of D. Printed
   for every group whether or not it is flattering, since increment 243
   lost a control by a silent guard.
 * The sign balance P(D > 0) is reported beside the mean, because
   increment 248 showed a pooled mean can cancel while the structure is
   large -- and because a mean and a sign balance fail differently.

WHY IT MATTERS. If D(k) has a mask, then every E1 test that scored it
against a mean-zero null was measuring the mask plus the fluctuation and
attributing all of it to the fluctuation. That does not by itself
overturn those tests -- most were "is there any signal at all", which is
threshold-free -- but it has to be checked rather than assumed.
"""
import numpy as np
import math
import time


def mobius_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    isp = np.ones(X + 1, dtype=bool); isp[:2] = False
    r = int(X ** 0.5)
    for p in range(2, r + 1):
        if isp[p]:
            isp[p * p::p] = False
            mu[p::p] = -mu[p::p]
            mu[p * p::p * p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    isp2 = np.ones(X + 1, dtype=bool); isp2[:2] = False
    for p in range(2, r + 1):
        if isp2[p]:
            isp2[p * p::p] = False
    for p in np.nonzero(isp2[: r + 1])[0]:
        val[int(p)::int(p)] //= int(p)
    mu[val > 1] = -mu[val > 1]
    mu[0] = 0
    return mu


def field(mu, N, ks):
    SQ = int(N ** 0.5)
    D = np.zeros(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        hi = N // k
        if hi <= SQ:
            continue
        ms = np.arange(SQ + 1, hi + 1, dtype=np.int64)
        D[i] = float((mu[ms].astype(np.int64)
                      * mu[N - k * ms].astype(np.int64)).sum())
    return D


def main():
    N = 20_000_000
    QS = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"sieve to {N}  t={time.time()-t0:.0f}s", flush=True)

    ks = np.arange(600, 4600, dtype=np.int64)
    D = field(mu, N, ks)
    print(f"D(k) for {len(ks)} values of k  t={time.time()-t0:.0f}s",
          flush=True)

    Nf = set(q for q in QS if N % q == 0)
    d = np.zeros(len(ks), dtype=np.int8)
    for q in QS:
        if q not in Nf:
            d += (ks % q == 0).astype(np.int8)
    print(f"  N = {N}, its small prime factors among QS: "
          f"{sorted(Nf)}")

    rng = np.random.default_rng(20260806)
    Dp = D[rng.permutation(len(D))]

    # Before reading any sign statistic: what does the distribution
    # actually look like? A sign imbalance of 0.34 against a recorded
    # skewness of +0.04 (sweep_A A6) is not consistent for a smooth
    # unimodal law, so one of the two is measuring something else.
    # Exact zeros are the first suspect: k > N/sqrt(N) leaves the sum
    # empty and D = 0 by construction, which counts as "not > 0".
    print(f"\n(A0) the distribution of D(k), before any test of it")
    nz = int((D == 0).sum())
    empty = int((N // ks <= int(N ** 0.5)).sum())
    nzD = D[D != 0]
    print(f"  n = {len(D)},  exact zeros = {nz} ({nz/len(D):.4f})")
    print(f"    of which structurally empty (N/k <= sqrt N) = {empty}")
    print(f"  mean {D.mean():+.2f}  median {np.median(D):+.2f}  "
          f"sd {D.std():.2f}")
    g1 = float(((D - D.mean()) ** 3).mean() / D.std() ** 3)
    print(f"  skewness {g1:+.4f} (SE {math.sqrt(6/len(D)):.4f})"
          f"   [sweep_A A6 recorded +0.0395, z = +0.86]")
    qs = np.percentile(D, [1, 5, 25, 50, 75, 95, 99])
    print("  percentiles 1/5/25/50/75/95/99: "
          + " ".join(f"{q:+.0f}" for q in qs))
    p_all = float((D > 0).mean())
    p_nz = float((nzD > 0).mean())
    print(f"  P(D > 0) counting zeros as not-positive = {p_all:.4f}")
    print(f"  P(D > 0 | D != 0)                       = {p_nz:.4f}"
          f"   (null 0.5, SE {1/(2*math.sqrt(len(nzD))):.4f},"
          f" z = {(p_nz-0.5)*2*math.sqrt(len(nzD)):+.2f})")
    print("  the second line is the one that tests symmetry; the first")
    print("  is contaminated by the atom at zero")

    # The atom is not noise, it is derivable. If q^2 | gcd(k, N) then
    # q^2 | N - mk for EVERY m, so mu(N - mk) = 0 identically and
    # D(k) = 0 exactly. That is a deterministic local mask on the
    # SUPPORT -- which is precisely what Conjecture L asserts a
    # mu-family carries. Counted rather than asserted:
    sqdiv = np.zeros(len(ks), dtype=bool)
    q = 2
    while q * q <= 64:
        if N % (q * q) == 0:
            sqdiv |= (ks % (q * q) == 0)
        q += 1
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        if N % (q * q) == 0:
            sqdiv |= (ks % (q * q) == 0)
    emptysel = (N // ks) <= int(N ** 0.5)
    pred = sqdiv | emptysel
    print(f"\n  the atom, derived: D(k) = 0 identically when "
          f"q^2 | gcd(k, N)")
    print(f"    N = {N} = 2^8 * 5^7, so 4 | k or 25 | k forces it")
    print(f"    predicted zeros (q^2 | gcd) or (empty range) = "
          f"{int(pred.sum())}")
    print(f"    observed zeros                               = {nz}")
    print(f"    all predicted are observed: "
          f"{bool((D[pred] == 0).all())}")
    extra = int(((D == 0) & ~pred).sum())
    print(f"    unpredicted zeros = {extra} "
          f"({extra/max(int((~pred).sum()),1):.4f} of the rest; a "
          f"smooth law with sd {D.std():.0f} would give "
          f"{1/(D.std()*math.sqrt(2*math.pi)):.4f})")
    print("    so the zeros are a support mask, not a symmetry failure")

    print(f"\n(A) mean D(k) by d(k) = #{{q<=23 : q|k, q not|N}}")
    print(f"{'d(k)':>6} {'n':>6} {'mean D':>11} {'SE':>10} {'z':>8} "
          f"{'perm z':>8} {'P(D>0)':>8} {'z sign':>8}")
    for j in range(0, int(d.max()) + 1):
        sel = d == j
        n = int(sel.sum())
        if n < 30:
            print(f"{j:>6} {n:>6}    (fewer than 30, reported not "
                  f"skipped)")
            continue
        x = D[sel]; xp = Dp[sel]
        se = float(x.std(ddof=1) / math.sqrt(n))
        sep = float(xp.std(ddof=1) / math.sqrt(n))
        p = float((x > 0).mean())
        print(f"{j:>6} {n:>6} {x.mean():>11.1f} {se:>10.1f} "
              f"{x.mean()/se:>8.2f} {xp.mean()/sep:>8.2f} "
              f"{p:>8.4f} {(p-0.5)*2*math.sqrt(n):>8.2f}")

    print(f"\n(B) the same, prime by prime (q | k against q not| k)")
    print(f"{'q':>4} {'q|N?':>5} {'n(q|k)':>7} {'mean|q|k':>11} "
          f"{'mean|q not|k':>13} {'z of diff':>10} {'perm z':>8}")
    for q in QS:
        a = ks % q == 0
        if a.sum() < 30:
            continue
        x, y = D[a], D[~a]
        sd = math.sqrt(x.var(ddof=1) / len(x) + y.var(ddof=1) / len(y))
        xp, yp = Dp[a], Dp[~a]
        sdp = math.sqrt(xp.var(ddof=1) / len(xp)
                        + yp.var(ddof=1) / len(yp))
        print(f"{q:>4} {'yes' if q in Nf else 'no':>5} {int(a.sum()):>7} "
              f"{x.mean():>11.1f} {y.mean():>13.1f} "
              f"{(x.mean()-y.mean())/sd:>10.2f} "
              f"{(xp.mean()-yp.mean())/sdp:>8.2f}")
    print("    predicted: a shift for q | k with q not| N, none for")
    print("    q | k with q | N, and nothing in the permuted column")
    print("DONE")


if __name__ == "__main__":
    main()
