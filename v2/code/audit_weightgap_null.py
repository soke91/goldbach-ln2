# -*- coding: utf-8 -*-
r"""
The control lab_weight_gap.py declined, run.

WHAT IS AT STAKE

lab_weight_gap.py backs Proposition {#prop:flatsum} and Remark
{#rem:weightgap}, and its result file declines a null: "Z1 is an
identity and the rest are comparisons between two weightings of the
SAME numbers, so a sign control would change both sides alike and
settle nothing."

That reason is the same shape as the one lab_predictable_part.py gave
before {#rem:splitnull} refuted it.  A coin would indeed move both
sums, but the claims are about their RATIO and about the shape of the
profile in the weight, and a coin gives both of those a reference:

  * Z2 says |sum H| / |sum (log k)H| is small, 0.1188 at the top N.
    If a coin gives the same, the smallness belongs to the weighting
    and not to mu.
  * Z4's diagnostic reads a single effective modulus off the near
    geometric profile in j, spread 1.0018 over four ratios.  A coin
    has no reason to be that clean, and if it is, the reading is an
    artefact of the parameterisation.

So the coin is run here on exactly lab_weight_gap.py's statistics.

BACKS: Remark {#rem:weightgapnull} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The ratio is mu's: |sum H|/|sum (log k)H| for mu is below every
      one of 8 coin draws, at every N.
  V2  The profile is mu's: the spread of the four consecutive ratios
      of |sum (log k)^j H| in j is smaller for mu than for every coin,
      at every N.
  V3  The decay ordering is mu's: mu's flat sum decays faster than its
      |E_3| (the exponent comparison of lab_weight_gap.py's Z3), and
      at most 4 of the 8 coins reproduce that ordering.
  V4  Sanity: the identity T_1 = sum H - C B_1 holds for the coin too,
      to better than 1e-12 relative -- it is algebra, so a failure
      means the coin construction is wrong.

REFUTATION RULE (fixed before the run)

  V1  REFUTED if any coin matches or beats mu at any N. That is the
      one that decides whether {#rem:weightgap}'s ratio means
      anything.
  V2  REFUTED if any coin is as clean as mu at any N, in which case
      the effective-modulus reading is an artefact.
  V3  REFUTED if 5 or more coins reproduce the ordering.
  V4  REFUTED at 1e-12 relative at any N or any draw.

  All four gate.

  NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere, with
  the field, the weights, the k-range and the truncation identical, so
  that the sign pattern on the long variable is the only difference.
  Eight draws.
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
OUT = os.path.join(ROOT, "results", "audit_weightgap_null.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
JS = [0.0, 0.25, 0.5, 0.75, 1.0]
THETA = 0.56
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


def phi_of(k):
    v, phi, d = k, 1, 2
    while d * d <= v:
        if v % d == 0:
            phi *= (d - 1)
            v //= d
            while v % d == 0:
                phi *= d
                v //= d
        d += 1
    if v > 1:
        phi *= (v - 1)
    return phi


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


def stats(N, ks, lg, iph, lam, sign):
    """lab_weight_gap.py's statistics for one sign array."""
    f0 = np.zeros(N, dtype=np.float64)
    idx = np.arange(1, N, dtype=np.int64)
    f0[1:] = lam[1:N] * sign[N - idx]
    C = float(f0.sum())
    A = np.empty(ks.size)
    for i, k in enumerate(ks):
        r = N % int(k)
        A[i] = f0[r::int(k)].sum() if r else f0[int(k)::int(k)].sum()
    del f0
    sg = sign[ks]
    H = sg * A
    T1 = float((sg * (A - C / iph)).sum())
    rhs = float(H.sum()) - C * float((sg / iph).sum())
    rel = abs(T1 - rhs) / max(abs(T1), 1e-300)
    S0 = abs(float(H.sum()))
    S1 = abs(float((lg * H).sum()))
    E3 = abs(float((sg * lg * (A - C / iph)).sum()))
    row = [abs(float(((lg ** j) * H).sum())) for j in JS]
    rr = [row[i + 1] / row[i] for i in range(len(row) - 1)
          if row[i] > 0]
    spread = (max(rr) / min(rr)) if len(rr) == len(JS) - 1 else float("inf")
    return dict(rel=rel, ratio=S0 / max(S1, 1e-300), spread=spread,
                flat=S0 / N, e3=E3 / N)


def loo(x, y, name, say):
    """Refit dropping each end in turn, and report the spread."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    f = [-float(np.polyfit(x[t], y[t], 1)[0])
         for t in (slice(None), slice(1, None), slice(0, -1))]
    sp = max(f) - min(f)
    say("  leave-one-out on %s: full %.4f, without the smallest N "
        "%.4f," % (name, f[0], f[1]))
    say("  without the largest %.4f -- spread %.4f" % (f[2], sp))
    say("SWEPT %s N-range %.4f" % (name, sp))
    return sp


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)

    rng = np.random.default_rng(SEED)
    sup = mu != 0
    coins = []
    for _ in range(COINS):
        c = np.zeros(NMAX + 1, dtype=np.float64)
        c[sup] = rng.choice([-1.0, 1.0], size=int(sup.sum()))
        coins.append(c)
    say("  %d coins on supp(mu^2), field and weights identical" % COINS)

    res = []
    for N in NS:
        PN = factor_set(N)
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if mu[k] != 0 and all(k % q for q in PN)])
        lg = np.log(ks.astype(float))
        iph = np.array([phi_of(int(k)) for k in ks], dtype=np.float64)
        m = stats(N, ks, lg, iph, lam, mu.astype(np.float64))
        cs = [stats(N, ks, lg, iph, lam, c) for c in coins]
        res.append((N, m, cs))
        say("  N = %-10d  #k = %d" % (N, ks.size))

    say()
    say("V4  the identity, for mu and for every coin")
    v4 = True
    worst = 0.0
    for N, m, cs in res:
        worst = max(worst, m["rel"], max(c["rel"] for c in cs))
    v4 = worst < 1e-12
    say("  worst relative error over all N and all draws: %.3e   %s"
        % (worst, "hold" if v4 else "REFUTED"))

    say()
    say("V1  |sum H| / |sum (log k)H|")
    say("  N            mu        coin min   coin median   coin max")
    v1 = True
    for N, m, cs in res:
        r = [c["ratio"] for c in cs]
        if m["ratio"] >= min(r):
            v1 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, m["ratio"], min(r), float(np.median(r)), max(r)))
    say("  V1 %s" % ("hold" if v1 else "REFUTED"))

    say()
    say("V2  spread of the four consecutive ratios in j")
    say("  N            mu        coin min   coin median   coin max")
    v2 = True
    for N, m, cs in res:
        r = [c["spread"] for c in cs]
        if m["spread"] >= min(r):
            v2 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, m["spread"], min(r), float(np.median(r)), max(r)))
    say("  V2 %s" % ("hold" if v2 else "REFUTED"))

    say()
    say("V3  does the flat sum decay faster than |E_3|?")
    x = np.log(np.array(NS, dtype=float))

    def expo(vals):
        return -float(np.polyfit(x, np.log(np.array(vals)), 1)[0])

    bf = expo([m["flat"] for N, m, cs in res])
    be = expo([m["e3"] for N, m, cs in res])
    say("  mu:   flat N^{-%.4f} against |E_3| N^{-%.4f}   ordering %s"
        % (bf, be, "holds" if bf > be else "fails"))
    loo(x, np.log(np.array([m["flat"] for N, m, cs in res])),
        "weightgap_flat_decay", say)
    loo(x, np.log(np.array([m["e3"] for N, m, cs in res])),
        "weightgap_E3_decay", say)
    nrep = 0
    say("  coin  flat exponent   |E_3| exponent   ordering")
    for j in range(COINS):
        cf = expo([cs[j]["flat"] for N, m, cs in res])
        ce = expo([cs[j]["e3"] for N, m, cs in res])
        ok = cf > ce
        nrep += int(ok)
        say("  %-5d %-15.4f %-16.4f %s" % (j, cf, ce, "holds" if ok
                                           else "fails"))
    v3 = (bf > be) and nrep <= 4
    say("  %d of %d coins reproduce the ordering   (cap 4)   %s"
        % (nrep, COINS, "hold" if v3 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). What the coin's ratio says the")
    say("  declined reason got wrong. A coin does move both sums, but")
    say("  not by the same factor: the ratio is a statistic in its own")
    say("  right and the coin gives it a reference. mu's value against")
    say("  the coin band, in units of the band's width:")
    say("  N            mu        band          (mu - min)/width")
    for N, m, cs in res:
        r = [c["ratio"] for c in cs]
        w = max(r) - min(r)
        say("  %-12d %-9.4f [%.4f, %.4f]  %+.2f"
            % (N, m["ratio"], min(r), max(r),
               (m["ratio"] - min(r)) / w if w > 0 else float("nan")))

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("the weight-gap statistics are mu's and the declined null was "
        "worth running" if ok else "REFUTED")

    head = [
        "STATISTIC: lab_weight_gap.py's own statistics recomputed for mu",
        "           and for 8 coins -- the identity residual, the ratio",
        "           |sum H|/|sum (log k)H|, the spread of the four",
        "           consecutive ratios of |sum (log k)^j H| in j, and the",
        "           fitted decay exponents of the flat sum and of |E_3|.",
        "NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere,",
        "      with the field, the weights, the k-range and the",
        "      truncation identical, so the sign pattern on the long",
        "      variable is the only difference. Eight draws. This is the",
        "      control lab_weight_gap.py declined on the ground that a",
        "      sign control would move both sums alike; it does move",
        "      both, but the claims are about their ratio and about the",
        "      shape of the profile, and those it gives a reference for.",
        "FIELD: N = 2e5 through 3.2e6 by doubling, theta' = 0.56; k over",
        "       the squarefree k < N^0.56 coprime to N; seed 20260808.",
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
