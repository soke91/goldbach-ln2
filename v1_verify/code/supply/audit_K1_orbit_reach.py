# -*- coding: utf-8 -*-
"""
Is K1's "full 63-divisor orbit" full?

THE STATEMENT UNDER TEST (§7.2 of v1/paper/wall_v1.tex):

    K1 | multiplicative Fejer kernel on the exact dilation ladder |
    **dead**: full 63-divisor orbit R^2 = 0.466. Half the field's
    energy is invisible to its entire multiplicative orbit.

THE ARITHMETIC. E1's dilate field carries the type-II cut, so

    D(v) = sum_{sqrt N < m <= N/v} mu(m) mu(N - m v)

is EMPTY unless v < sqrt N. The orbit of a base point k is
{D(sk) : s | 30030 squarefree}, so the column at s is a column of
zeros unless

    s * k < sqrt N.

The largest s is 30030. A base point therefore contributes its whole
orbit only when N > (30030 k)^2, which at k = 2 is 3.6e9 and grows
quadratically. **The 63-divisor orbit of the type-II field is not
reachable at any N this program computes**, and at the parameters K1
used it is truncated hard.

This script measures the truncation exactly, at v1's own configuration
and at others, and then asks what R^2 = 0.466 is the R^2 of.

PRE-REGISTRATION (written before the run).

  (1) RULE. Count, at v1's configuration (N = 199,999,998,
      k in [500, 900], s | 30030), how many of the 63 orbit columns
      are non-degenerate. If that number is close to 63 the objection
      is void. If it is a small fraction, then "full 63-divisor orbit"
      names a design that was not run, and the pure-noise baseline
      the measured R^2 must be read against is p/n for the LIVE p,
      not 63/400.
  (2) RULE. Rerun the regression on the live columns only. If R^2 is
      unchanged, the dead columns cost nothing and only the wording is
      wrong. Report the live-column count as the p in the adjusted
      R^2.
  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      sqrt(199999998) = 14142, and k >= 500, so s must be under 28.
      The squarefree divisors of 30030 below 28 are
      1,2,3,5,6,7,10,11,13,14,15,21,22,26 -- fourteen, of which
      thirteen are non-trivial, and fewer as k grows. So I predict
      about 9 to 13 live columns of 63, i.e. roughly 80% of the orbit
      never evaluated, and I predict R^2 on the live columns
      reproduces 0.466 closely because the dead columns contribute
      nothing to a least-squares fit.
  (4) The companion run `audit_K1_orbit.py` puts the same design on
      the FULL field (no type-II cut), where all 63 columns are live,
      and gets R^2 = 0.904 against a coin null of 0.44 and a
      permutation null of 0.17. Under K1's own pre-registered
      thresholds (>= 0.95 alive, <= 0.9 dead) that is MARGINAL, not
      dead. Both readings are reported here; they are different
      objects and the paper does not say which one K1 used.
"""
import sys
import math
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

PR = [2, 3, 5, 7, 11, 13]


def sqfree_divisors():
    out = [1]
    for p in PR:
        out += [x * p for x in out]
    return sorted(out)


def sieve_mu(X):
    mu = np.ones(X + 1, dtype=np.int8)
    comp = np.zeros(X + 1, dtype=bool)
    for p in range(2, int(X ** 0.5) + 1):
        if not comp[p]:
            comp[p * p::p] = True
            mu[p::p] = -mu[p::p]
            mu[p * p::p * p] = 0
    # primes above sqrt(X) each divide their multiples exactly once
    rest = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if not comp[p]:
            q = p
            while q <= X:
                rest[q::q] //= p
                q *= p
    mu[rest > 1] = -mu[rest > 1]
    mu[0] = 0
    del comp, rest
    return mu


def main():
    S = [s for s in sqfree_divisors() if s > 1]
    print("Is K1's 'full 63-divisor orbit' full?")
    print()
    print("(1) live columns: D(sk) is empty unless s*k < sqrt(N)")
    print(f"    {'N':>14} {'sqrt N':>9} {'k range':>13} "
           f"{'live cols (min..max of 63)':>28}")
    for N, klo, khi in ((199_999_998, 500, 900),
                        (199_999_998, 2, 10),
                        (10 ** 10, 500, 900),
                        (10 ** 12, 500, 900)):
        rt = math.isqrt(N)
        lo = min(sum(1 for s in S if s * k < rt) for k in (klo, khi))
        hi = max(sum(1 for s in S if s * k < rt) for k in (klo, khi))
        print(f"    {N:>14} {rt:>9} {f'[{klo}, {khi}]':>13} "
              f"{f'{lo} .. {hi}':>28}")
    need = (30030 * 2) ** 2
    print(f"    the whole orbit is live only when N > (30030 k)^2,")
    print(f"    i.e. N > {need:.3e} even at k = 2.")
    print()

    # ---- reproduce v1's configuration -------------------------------
    t0 = time.time()
    N = 199_999_998
    rng = np.random.default_rng(20260830)
    ks = rng.choice(np.arange(500, 900), size=400, replace=False)
    rt = int(N ** 0.5)
    live = np.array([[1 if s * int(k) < rt else 0 for s in S]
                     for k in ks])
    percol = live.sum(axis=0)
    print(f"(2) v1's configuration, N = {N}, sqrt N = {rt}")
    print(f"    columns with at least one non-empty row: "
          f"{int((percol > 0).sum())} of {len(S)}")
    print(f"    columns non-empty for every row: "
          f"{int((percol == len(ks)).sum())} of {len(S)}")
    print(f"    live entries in the design matrix: "
          f"{100*live.mean():.1f}%")
    dead = [S[j] for j in range(len(S)) if percol[j] == 0]
    print(f"    the {len(dead)} columns that are identically zero:")
    print(f"      {dead}")
    print()

    print("    building the field (this is the slow part)", flush=True)
    mu = sieve_mu(N)
    print(f"    mu ready t={time.time()-t0:.0f}s", flush=True)

    def D_of(v):
        m = np.arange(rt + 1, N // v + 1, dtype=np.int64)
        if not len(m):
            return 0.0
        return float(np.dot(mu[m].astype(np.float64),
                            mu[N - v * m].astype(np.float64)))

    y = np.array([D_of(int(k)) for k in ks])
    X = np.zeros((len(ks), len(S)))
    for i, k in enumerate(ks):
        for j, s in enumerate(S):
            X[i, j] = D_of(s * int(k))
    print(f"    field ready t={time.time()-t0:.0f}s", flush=True)

    def r2(Xm):
        c, *_ = np.linalg.lstsq(Xm, y, rcond=None)
        r = y - Xm @ c
        return 1.0 - float(r @ r) / float(((y - y.mean()) ** 2).sum())

    keep = percol > 0
    r_all = r2(X)
    r_live = r2(X[:, keep])
    p_live = int(keep.sum())
    n = len(ks)
    print(f"    R^2 on all {len(S)} columns        {r_all:.4f}"
          f"   (v1 quotes 0.466)")
    print(f"    R^2 on the {p_live} live columns   {r_live:.4f}")
    print(f"    pure-noise R^2 = p/n: quoted basis {len(S)}/{n} = "
          f"{len(S)/n:.4f}, actual {p_live}/{n} = {p_live/n:.4f}")
    adj = 1 - (1 - r_live) * (n - 1) / (n - p_live - 1)
    print(f"    adjusted R^2 on the live columns  {adj:.4f}"
          f"   -> blind share {1-adj:.4f}")
    print()
    print("(4) the same design on the FULL field, all 63 columns live")
    print("    (audit_K1_orbit.py, N = 1e8, k coprime to 30030):")
    print("      R^2 = 0.9040, adjusted 0.8860")
    print("      permutation null 0.165 +- 0.051, coin null 0.44")
    print("      K1's own thresholds: >= 0.95 alive, <= 0.9 dead")
    print("      -> MARGINAL, and the coin reaches 0.44 of it")
    print("DONE")


if __name__ == "__main__":
    main()
