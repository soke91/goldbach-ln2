# -*- coding: utf-8 -*-
r"""
The control lab_level_forecast.py points at, checked against what it
actually uses.

WHAT IS AT STAKE

lab_level_forecast.py declines a control by pointing elsewhere: "gamma
carries whatever the sign pattern contributed, which
lab_dilate_extrapolation.py measured for mu and for a coin separately.
The control that matters was run there."

The pointer is real but it misses.  lab_dilate_extrapolation.py
measures c'(N), the MEDIAN of |A(N;k)|/sqrt(N/k), for mu and for
coins.  lab_level_forecast.py's own diagnostic then rejects the median
-- "gamma was calibrated on the MEDIAN, but B is a SUM and needs the
MEAN" -- and rebuilds the forecast on a MEAN-based gamma of 0.9803.
The coin was never measured for the mean.  That is the failure mode
{#rem:splitnull} found in lab_predictable_part.py: the named control
tests a different statistic from the claim.

By the criterion of {#rem:weightgapnull} and {#rem:extendnull} the
control is also usable here: a mean of absolute values is bounded away
from zero, so a coin cannot make it ill-conditioned.  So it is run.

BACKS: Remark {#rem:forecastnull} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  X1  Direction: the coin's mean-based gamma is below mu's at every
      draw, since a coin gets square-root cancellation inside A and
      mu does not.
  X2  The statistic is well conditioned under the control, as the
      criterion requires: the draws' spread of gamma is under 0.15 of
      their median.
  X3  The calibration is mu's: mu's gamma lies outside the range
      spanned by the 8 draws.
  X4  And so is the forecast: the N at which the model's K* reaches
      sqrt(N) -- the level-one-half crossing that {#rem:forecast}
      quotes -- lies, for mu, outside the range spanned by the draws.

REFUTATION RULE (fixed before the run)

  X1  REFUTED if any draw reaches mu's gamma. A failure would
      contradict {#rem:whycoinwins} and mean the coin is wrong.
  X2  REFUTED if the spread reaches 0.15 of the median, in which case
      the control is not usable and the decline was right after all.
  X3  REFUTED if mu's gamma falls inside the draws' range. That is the
      one that decides whether the calibration is mu's.
  X4  REFUTED if mu's crossing falls inside the draws' range.

  All four gate.

  NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere, with
  the field, the k-range and theta' identical, so the sign pattern on
  the long variable is the only difference. Eight draws, and gamma
  fitted from the MEAN, which is the statistic the forecast uses.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_forecast_null.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
FITN = 3                      # gamma fitted on the three smallest N
THETA = 0.56
CLIM = 4_000_000
SKMAX = 200_000               # S(K) built exactly to here
COINS = 8
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n):
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return pr, lam, mu


def factor_set(n):
    v, out, d = n, set(), 2
    while d * d <= v:
        if v % d == 0:
            out.add(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        out.add(v)
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    A_, S_ = artin, twin
    for q in (2, 5):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))
    THR = S_ * (1.0 - A_)
    say("  threshold S(N)(1-A(N)) for the family {2,5} = %.6f" % THR)

    # S(K) = sum over admissible k < K of (log k) k^{-1/2}
    kk = np.arange(2, SKMAX + 1, dtype=np.int64)
    adm = sqf[kk] & (kk % 2 != 0) & (kk % 5 != 0)
    kk = kk[adm]
    SK = np.cumsum(np.log(kk.astype(float)) / np.sqrt(kk.astype(float)))

    def Sof(K):
        j = int(np.searchsorted(kk, K))
        return float(SK[min(j, SK.size) - 1]) if j else 0.0

    rng = np.random.default_rng(SEED)
    sup = mu != 0
    signs = [mu.astype(np.float64)]
    for _ in range(COINS):
        c = np.zeros(NMAX + 1, dtype=np.float64)
        c[sup] = rng.choice([-1.0, 1.0], size=int(sup.sum()))
        signs.append(c)
    say("  1 mu and %d coins on supp(mu^2), field and k-range identical"
        % COINS)

    means = [[] for _ in signs]
    for N in NS:
        PN = factor_set(N)
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if sqf[k] and all(k % q for q in PN)])
        scale = np.sqrt(N / ks.astype(float))
        idx = np.arange(1, N, dtype=np.int64)
        for j, sg in enumerate(signs):
            f0 = np.zeros(N, dtype=np.float64)
            f0[1:] = lam[1:N] * sg[N - idx]
            A = np.empty(ks.size)
            for i, k in enumerate(ks):
                r = N % int(k)
                A[i] = (f0[r::int(k)].sum() if r
                        else f0[int(k)::int(k)].sum())
            del f0
            means[j].append(float((np.abs(A) / scale).mean()))
        say("  N = %-10d  #k = %d" % (N, ks.size))

    say()
    say("X1/X2  the mean of |A(N;k)|/sqrt(N/k)")
    say("  N            mu        coin min   coin median   coin max")
    x1 = True
    for i, N in enumerate(NS):
        c = [means[j][i] for j in range(1, len(signs))]
        if max(c) >= means[0][i]:
            x1 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, means[0][i], min(c), float(np.median(c)), max(c)))
    say("  X1 %s" % ("hold" if x1 else "REFUTED"))

    # gamma from |A| = gamma sqrt(N/k) sqrt(log N), fitted out of sample
    gam = [float(np.mean([means[j][i] / math.sqrt(math.log(NS[i]))
                          for i in range(FITN)]))
           for j in range(len(signs))]
    cg = gam[1:]
    sp = (max(cg) - min(cg)) / float(np.median(cg))
    x2 = sp < 0.15
    say()
    say("X2  the draws' spread of gamma, as a fraction of their median:"
        " %.4f  (cap 0.15)   %s" % (sp, "hold" if x2 else "REFUTED"))
    x3 = not (min(cg) <= gam[0] <= max(cg))
    say("X3  mu's gamma %.4f against the draws' range [%.4f, %.4f]   %s"
        % (gam[0], min(cg), max(cg), "hold" if x3 else "REFUTED"))

    def cross(g):
        """N at which the model's K* reaches sqrt(N).

        K* solves S(K*) = THR sqrt(N/log N)/g, and S is increasing, so
        K* < sqrt(N) exactly when S(sqrt(N)) > that target. Below the
        crossing the target is the smaller of the two; the bisection
        must move UP there, which an earlier version had backwards and
        which sent every answer to the top of the bracket.
        """
        lo, hi = math.log(1e4), math.log(1e12)
        for _ in range(300):
            mid = 0.5 * (lo + hi)
            Nn = math.exp(mid)
            want = THR * math.sqrt(Nn / mid) / g
            if Sof(math.sqrt(Nn)) > want:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi) / math.log(10)

    say()
    say("X4  where the model's K* reaches sqrt(N)")
    say("  who       gamma     crossing at 10^x")
    xs = [cross(g) for g in gam]
    say("  %-9s %-9.4f %.2f" % ("mu", gam[0], xs[0]))
    for j in range(1, len(signs)):
        say("  coin %-4d %-9.4f %.2f" % (j - 1, gam[j], xs[j]))
    cc = xs[1:]
    x4 = not (min(cc) <= xs[0] <= max(cc))
    say("  mu's %.2f against the draws' range [%.2f, %.2f]   %s"
        % (xs[0], min(cc), max(cc), "hold" if x4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). What the pointer missed. The control")
    say("  lab_level_forecast.py points at measures the MEDIAN of")
    say("  |A|/sqrt(N/k); the forecast uses the MEAN. Both, for mu:")
    say("  N            mean      median    mean/median")
    for i, N in enumerate(NS):
        PN = factor_set(N)
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if sqf[k] and all(k % q for q in PN)])
        scale = np.sqrt(N / ks.astype(float))
        idx = np.arange(1, N, dtype=np.int64)
        f0 = np.zeros(N, dtype=np.float64)
        f0[1:] = lam[1:N] * mu[N - idx].astype(np.float64)
        A = np.empty(ks.size)
        for i2, k in enumerate(ks):
            r = N % int(k)
            A[i2] = (f0[r::int(k)].sum() if r
                     else f0[int(k)::int(k)].sum())
        del f0
        v = np.abs(A) / scale
        say("  %-12d %-9.4f %-9.4f %.4f"
            % (N, float(v.mean()), float(np.median(v)),
               float(v.mean() / np.median(v))))
    say("  The two differ by about half again, which is exactly why")
    say("  lab_level_forecast.py rebuilt gamma -- and why the control it")
    say("  pointed at does not reach the quantity it then used.")

    say()
    say("=" * 70)
    ok = x1 and x2 and x3 and x4
    say("the forecast's calibration and crossing are mu's, and the "
        "pointer it declined behind did not cover them"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the mean of |A(N;k)|/sqrt(N/k) over the admissible k,",
        "           for mu and for 8 coins; the mean-based gamma fitted",
        "           out of sample on the three smallest N; and the N at",
        "           which the model K* solved from",
        "           B/N = gamma sqrt(log N/N) S(K) against the threshold",
        "           S(N)(1-A(N)) first reaches sqrt(N).",
        "NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere,",
        "      field, k-range and theta' identical. This is the control",
        "      lab_level_forecast.py declined by pointing at",
        "      lab_dilate_extrapolation.py, which measures the MEDIAN of",
        "      the same ratio; the forecast uses the MEAN, so the pointer",
        "      does not reach the quantity in question.",
        "FIELD: N = 2e5 through 3.2e6 by doubling, theta' = 0.56; k over",
        "       the squarefree k < N^0.56 coprime to N; S(K) enumerated",
        "       exactly over the squarefree k coprime to 10 up to 2e5;",
        "       the threshold is that of the family {2,5}; seed 20260808.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not ok:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
