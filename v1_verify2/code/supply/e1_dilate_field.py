# -*- coding: utf-8 -*-
"""
E1's normalisation figures and sec:R4's block diagnostics.
(v1_verify2, Phase 1, blind.)

STATEMENTS UNDER TEST, verbatim from v1/paper/wall_v1.tex:

  sec:supply  D(k) = sum_{sqrt N < m <= N/k} mu(m) mu(N-mk),
              "at the square-root normalisation the demand would be
               sum|D|^2 / sum M_k << 8.8e-6, which our own measurements
               refute directly (0.305, 0.310, 0.319 at N=1e8)"

  sec:R4      "On 8000 values of k --- the entire band on which the
               type-II field is non-empty, since m > sqrt N forces
               k < sqrt N --- the block ratios are flat (B=8: 0.958 and
               1.023 at two N, against B=1 baselines 0.980 and 0.979),
               and the sharper diagnostic, the lag-1 autocorrelation of
               D(k)/sqrt(supp(k)), reads +0.0104 and +0.0127 against a
               standard error of 0.0112: dead zero"

              and the exact identity
               sum_k sum_{m: mk<=N-1} mu(m)mu(N-mk) = mu(N-1),
              "verified by brute force at N = 5000 and N = 20000"

PRE-REGISTRATION (fixed before this ran).

  Decision rule.
   (1) The exact identity: recompute at N=5000, 20000 and at four more N.
       PASS iff it equals mu(N-1) exactly at every N.
   (2) sum|D|^2 / sum M_k on dyadic bands at N=1e8. REPRODUCED iff three
       band values land in 0.30--0.32.
   (3) The block ratio, defined as sum_j S_B(j)^2 / sum_j (sum_{k in
       block j} supp(k)) -- the only reading under which a B=1 baseline
       can differ from 1, as the paper's 0.980/0.979 requires.
       REPRODUCED iff B=1 gives ~0.98 and B=8 stays flat.
   (4) The lag-1 autocorrelation of D(k)/sqrt(supp(k)) over the whole
       non-empty band, against its own standard error.

  Predictions written before running.
   (2) REPRODUCES near 0.3226 = prod_p (1-2/p^2), which is the density of
       m with mu(m) and mu(N-mk) both nonzero. If so the measured 0.305
       -- 0.319 is EXACT square-root cancellation on the surviving
       support, and the paper's phrasing ("refute directly") is about the
       wrong normalisation, not about a failure of cancellation.
   (3) REPRODUCES; I predict the B=1 baseline is 0.98 for the same
       reason -- it is sum D^2 / sum supp, which is 1 under exact
       square-root cancellation.
   (4) REPRODUCES near zero. I predict the quoted standard error 0.0112
       is 1/sqrt(n_k) with n_k ~ 8000, i.e. a COUNT-based error bar --
       the very thing prop:coh says is wrong for correlated summands.
       Whether it is wrong HERE is a separate question and is tested by
       a permutation null, which is the honest null for a lag statistic.
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np


def mobius_upto(X):
    """mu(0..X) as int8, memory-lean."""
    sieve = np.ones(X + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(X ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = False
    primes = np.nonzero(sieve)[0]
    del sieve
    mu = np.ones(X + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes:
        p = int(p)
        mu[p:: p] *= -1
        pp = p * p
        if pp <= X:
            mu[pp:: pp] = 0
    return mu, primes


def D_band(N, mu, kmax, chunk=4_000_000):
    """D(k) and supp(k) for k = 1..kmax, D(k)=sum_{sqrtN<m<=N/k}
    mu(m)mu(N-mk)."""
    root = int(N ** 0.5)
    D = np.zeros(kmax + 1)
    SUP = np.zeros(kmax + 1)
    for k in range(1, kmax + 1):
        mhi = N // k
        if mhi <= root:
            continue
        lo = root + 1
        acc = 0.0
        sup = 0
        while lo <= mhi:
            hi = min(mhi, lo + chunk - 1)
            m = np.arange(lo, hi + 1, dtype=np.int64)
            a = mu[m]
            idx = N - m * k
            b = mu[idx]
            acc += float(np.dot(a.astype(np.float64), b.astype(np.float64)))
            sup += int(np.count_nonzero(a.astype(np.int16) * b))
            lo = hi + 1
        D[k] = acc
        SUP[k] = sup
    return D, SUP


def main():
    print("e1_dilate_field   (v1_verify2 Phase 1, blind)")
    print("=" * 74)

    # ---------------------------------------------------------- (1)
    print()
    print("--- (1) the exact divisor-switch identity ------------------------")
    print("    sum_{k>=1} sum_{m: mk<=N-1} mu(m)mu(N-mk) =?= mu(N-1)")
    muS, _ = mobius_upto(60000)
    ok = True
    for N in (5000, 20000, 7919, 30030, 40001, 50000):
        tot = 0.0
        for k in range(1, N):
            mhi = (N - 1) // k
            if mhi < 1:
                break
            m = np.arange(1, mhi + 1)
            tot += float(np.dot(muS[m].astype(np.float64),
                                muS[N - m * k].astype(np.float64)))
        good = abs(tot - muS[N - 1]) < 1e-9
        ok &= good
        print(f"    N={N:>6}:  sum = {tot:>8.1f}   mu(N-1) = "
              f"{muS[N - 1]:>3}   {'ok' if good else 'MISMATCH'}")
    print(f"    identity: {'CONFIRMED' if ok else 'REFUTED'}")
    del muS

    # ---------------------------------------------------------- (2)-(4)
    N = 100_000_000
    print()
    print(f"--- sieving mu to {N:,} ---")
    mu, _ = mobius_upto(N)
    root = int(N ** 0.5)
    print(f"    sqrt N = {root:,};  the type-II field is non-empty exactly "
          f"for k < {root:,}")
    print(f"    [paper: 'on 8000 values of k --- the entire band']")

    print()
    print("--- (2) E1 at the two normalisations, dyadic bands, N=1e8 -------")
    print(f"    {'band k~K':>14}{'count':>8}{'sum|D|^2/sum M_k':>19}"
          f"{'sum|D|^2/sum supp':>20}{'sum|D|^2/sum M_k^2':>21}")
    for K in (100, 200, 400):
        ks = np.arange(K, 2 * K)
        D, S = D_band(N, mu, int(ks[-1]))
        d = D[ks]
        s = S[ks]
        M = N / ks - root
        print(f"    {K:>7}-{2 * K - 1:<6}{len(ks):>8}"
              f"{(d ** 2).sum() / M.sum():>19.4f}"
              f"{(d ** 2).sum() / s.sum():>20.4f}"
              f"{(d ** 2).sum() / (M ** 2).sum():>21.3e}")
    print(f"    [paper: 0.305, 0.310, 0.319 at N=1e8]")
    print(f"    prod_p (1-2/p^2) = 0.3226340989 -- the density of m with")
    print(f"    mu(m) and mu(N-mk) both nonzero")
    print(f"    [paper's target at the trivial normalisation: 8.8e-6]")

    print()
    print("--- (3)-(4) sec:R4 over the whole non-empty band -----------------")
    kmax = root - 1
    D, S = D_band(N, mu, kmax)
    k = np.arange(1, kmax + 1)
    d = D[1:]
    s = S[1:]
    live = s > 0
    print(f"    k = 1..{kmax:,};  non-degenerate (supp>0): "
          f"{int(live.sum()):,}")

    import os
    cdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "..", "v1_verify2_log", "cache")
    cdir = os.path.abspath(cdir)
    os.makedirs(cdir, exist_ok=True)
    np.savez_compressed(os.path.join(cdir, f"dilate_{N}.npz"),
                        D=D, SUP=S, N=np.int64(N))

    print()
    print("    the block ratio: the paper does not state its weight, and")
    print("    the two natural weightings disagree at B=1, where a")
    print("    baseline of 1 is what square-root cancellation predicts.")
    print(f"    {'B':>6}{'blocks':>9}{'ratio of sums':>16}"
          f"{'mean of ratios':>17}{'live-k only':>14}")
    for B in (1, 2, 4, 8, 16, 64, 512):
        nb = len(d) // B
        sb = d[: nb * B].reshape(nb, B).sum(axis=1)
        ss = s[: nb * B].reshape(nb, B).sum(axis=1)
        ratio_of_sums = (sb ** 2).sum() / ss.sum()
        ok = ss > 0
        mean_of_ratios = float(np.mean(sb[ok] ** 2 / ss[ok]))
        dl, sl = d[live], s[live]
        nb2 = len(dl) // B
        sb2 = dl[: nb2 * B].reshape(nb2, B).sum(axis=1)
        ss2 = sl[: nb2 * B].reshape(nb2, B).sum(axis=1)
        live_only = float(np.mean(sb2 ** 2 / ss2))
        print(f"    {B:>6}{nb:>9}{ratio_of_sums:>16.4f}"
              f"{mean_of_ratios:>17.4f}{live_only:>14.4f}")
    print(f"    [paper: B=8 gives 0.958 and 1.023 at two N, against B=1")
    print(f"     baselines 0.980 and 0.979]")

    print()
    x = d[live] / np.sqrt(s[live])
    x = x - x.mean()
    lag1 = float(np.dot(x[:-1], x[1:]) / np.dot(x, x))
    nn = len(x)
    print(f"    lag-1 autocorrelation of D(k)/sqrt(supp(k)) = {lag1:+.4f}")
    print(f"    count-based s.e. 1/sqrt(n) = {1 / np.sqrt(nn):.4f}"
          f"   (n = {nn:,})")
    print(f"    [paper: +0.0104 and +0.0127 against a s.e. of 0.0112]")
    rng = np.random.default_rng(7)
    perm = np.empty(400)
    for i in range(400):
        y = rng.permutation(x)
        perm[i] = float(np.dot(y[:-1], y[1:]) / np.dot(y, y))
    print(f"    permutation null (400 draws): mean {perm.mean():+.5f}, "
          f"sd {perm.std(ddof=1):.5f}")
    print(f"    z against the permutation null = "
          f"{(lag1 - perm.mean()) / perm.std(ddof=1):+.2f}")
    print(f"    1/sqrt(8000) = {1 / np.sqrt(8000):.4f}  <- the paper's "
          f"quoted 0.0112 is a COUNT-based bar")
    return 0


if __name__ == "__main__":
    sys.exit(main())
