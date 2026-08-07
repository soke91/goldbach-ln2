# -*- coding: utf-8 -*-
"""
Re-verification of kill-test K1 of v1/paper/wall_v1.tex.

THE STATEMENT UNDER TEST, verbatim (§7.2):

    K1 | multiplicative Fejer kernel on the exact dilation ladder |
    **dead**: full 63-divisor orbit R^2 = 0.466. Half the field's
    energy is invisible to its entire multiplicative orbit.

The design, from the closure's own record: represent D(k) as an
optimally weighted combination of its orbit {D(sk) : s squarefree,
s | 30030}; if the representation is near-exact, orbit overlap across
k linearises the k-average. The test is the full-orbit least-squares
R^2 of D(k) on {D(sk)}, with the pre-registered thresholds
R^2 >= 0.95 alive, R^2 <= 0.9 dead.

WHAT IS MISSING. `v1`'s own `audit_killtest_nulls.py` classifies this
file as carrying **no null at all**, and `R^2` on 63 predictors is not
a quantity that can be read against zero: pure noise gives
`R^2 = 63/n`. The verdict (0.466 against a threshold of 0.9) does not
depend on the null, but the SENTENCE does -- "half the field's energy
is invisible to its entire multiplicative orbit" reads 1 - R^2 as a
share of energy, and 1 - R^2 is not that until the free parameters are
paid for.

METHOD HERE. Written from the statement. D is the dilate field over
its full m-range, D(k) = sum_{m <= (N-1)/k} mu(m) mu(N - mk), which is
the only version defined for every k in the orbit (the type-II cut
m > sqrt N empties the band for k > sqrt N, and the orbit multiplies k
by up to 30030). Three nulls, none of which v1 ran:

  ADJ   the adjusted R^2, 1 - (1-R^2)(n-1)/(n-p-1): what 1 - R^2 is
        once 63 free parameters are paid for.
  PERM  regress D(k) on the orbit of a DIFFERENT k -- the columns
        rebuilt from a fixed permutation of the base points. This
        destroys the ladder relation and keeps everything else, so it
        is the null for "does the multiplicative orbit matter".
  COIN  Lemma 17: rebuild the whole field with mu replaced by signs on
        the same support, and rerun. This is the null for "is the
        orbit's reach a property of mu".

PRE-REGISTRATION (written before the run).

  (1) RULE. The DEAD verdict stands or falls on R^2 <= 0.9 and is not
      in doubt. What is tested is the sentence: report the adjusted
      R^2, and the PERM and COIN nulls. If R^2 does not clear PERM,
      the orbit carries no information at all and the closure is
      stronger than stated. If the COIN reaches the real R^2, the
      orbit's reach is not mu's.
  (2) PREDICTION, recorded so it cannot be reported as a surprise.
      The orbit shares terms with the base sum by construction --
      D(sk) sums mu(m)mu(N-msk) and D(k) sums mu(m)mu(N-mk), and the
      m-ranges nest -- so I expect real R^2 well above PERM, and I
      expect the COIN to reach it too, because the sharing is
      arithmetic and not about mu. I expect the adjusted R^2 to come
      out materially below the raw one, so that "half" overstates the
      orbit's reach.
  (3) The orbit here is s | 30030 as stated, so 63 non-trivial
      divisors, and n = 400 base points as in v1. Pure-noise R^2 is
      then 63/400 = 0.157, which is the number "0.466" has to be read
      against.
"""
import sys
import math
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

RAD = 30030          # 2*3*5*7*11*13
NBASE = 400


def sieve_mu(X):
    mu = np.ones(X + 1, dtype=np.int8)
    comp = np.zeros(X + 1, dtype=bool)
    for p in range(2, int(X ** 0.5) + 1):
        if not comp[p]:
            comp[p * p::p] = True
    primes = np.nonzero(~comp)[0]
    for p in primes[primes >= 2]:
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= X:
            mu[p * p::p * p] = 0
    mu[0] = 0
    return mu


def sqfree_divisors(n):
    ps = []
    m = n
    d = 2
    while d * d <= m:
        if m % d == 0:
            ps.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        ps.append(m)
    out = [1]
    for p in ps:
        out += [x * p for x in out]
    return sorted(out)


def field(mu, N, vals):
    """D(v) = sum_{m <= (N-1)/v} mu(m) mu(N - m v), for v in vals."""
    out = {}
    for v in vals:
        m = np.arange(1, (N - 1) // v + 1, dtype=np.int64)
        out[v] = float(np.dot(mu[m].astype(np.float64),
                              mu[N - m * v].astype(np.float64)))
    return out


def r2_of(y, Xc):
    Xc = Xc - Xc.mean(axis=0)
    yc = y - y.mean()
    coef, *_ = np.linalg.lstsq(Xc, yc, rcond=None)
    res = yc - Xc @ coef
    return 1.0 - float(res @ res) / float(yc @ yc)


def main():
    t0 = time.time()
    N = 100_000_000
    S = [s for s in sqfree_divisors(RAD)]
    Snz = [s for s in S if s > 1]
    ks = []
    k = 1
    while len(ks) < NBASE:
        k += 1
        if math.gcd(k, RAD) == 1:
            ks.append(k)
    kmaxprod = max(ks) * RAD
    print("Re-verification of kill-test K1 (the multiplicative orbit)")
    print()
    print(f"  N = {N}, orbit s | {RAD} ({len(Snz)} non-trivial),")
    print(f"  {len(ks)} base points k coprime to {RAD}, k <= {max(ks)},")
    print(f"  largest s*k = {kmaxprod} <= N: "
          f"{'yes' if kmaxprod < N else 'NO'}")
    print(f"  pure-noise R^2 = {len(Snz)}/{len(ks)} = "
          f"{len(Snz)/len(ks):.4f}")
    print()

    mu = sieve_mu(N)
    print(f"  sieve t={time.time()-t0:.0f}s", flush=True)
    need = sorted({k * s for k in ks for s in S})
    D = field(mu, N, need)
    print(f"  field on {len(need)} points  t={time.time()-t0:.0f}s",
          flush=True)

    y = np.array([D[k] for k in ks])
    X = np.array([[D[k * s] for s in Snz] for k in ks])
    r2 = r2_of(y, X)
    n, p = len(ks), len(Snz)
    adj = 1.0 - (1.0 - r2) * (n - 1) / (n - p - 1)

    print(f"  real R^2                       {r2:.4f}"
          f"   (v1 quotes 0.466)")
    print(f"  adjusted R^2                   {adj:.4f}")
    print(f"  1 - R^2  ('half is invisible') {1-r2:.4f}")
    print(f"  1 - adjusted R^2               {1-adj:.4f}")
    print()

    rng = np.random.default_rng(1729)
    perm = []
    for _ in range(200):
        sig = rng.permutation(n)
        perm.append(r2_of(y, X[sig]))
    perm = np.array(perm)
    print(f"  PERM null (orbit of a different k, 200 draws)")
    print(f"    mean {perm.mean():.4f}, sd {perm.std():.4f}, "
          f"max {perm.max():.4f}")
    print(f"    real R^2 clears it: "
          f"{'YES' if r2 > perm.max() else 'NO'}"
          f"   z = {(r2-perm.mean())/perm.std():+.1f}")
    print()

    supp = np.nonzero(mu != 0)[0]
    coins = []
    for c in range(3):
        eps = np.zeros(N + 1, dtype=np.int8)
        eps[supp] = rng.choice([-1, 1], size=len(supp)).astype(np.int8)
        Dc = field(eps, N, need)
        yc = np.array([Dc[k] for k in ks])
        Xc = np.array([[Dc[k * s] for s in Snz] for k in ks])
        coins.append(r2_of(yc, Xc))
        print(f"  COIN draw {c+1}: R^2 = {coins[-1]:.4f}"
              f"   t={time.time()-t0:.0f}s", flush=True)
    coins = np.array(coins)
    print(f"  COIN null: mean {coins.mean():.4f}, "
          f"range [{coins.min():.4f}, {coins.max():.4f}]")
    print()
    print("(1) the DEAD verdict: R^2 <= 0.9 -> "
          f"{'DEAD' if r2 <= 0.9 else 'NOT DEAD'}")
    print("    the sentence: '1 - R^2' is the orbit's blind share only")
    print(f"    after the {p} free parameters are paid for, which makes")
    print(f"    it {1-adj:.3f} and not {1-r2:.3f}.")
    print("    Whether the reach is mu's is the COIN line above.")
    print("DONE")


if __name__ == "__main__":
    main()
