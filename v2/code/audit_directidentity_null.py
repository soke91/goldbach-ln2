# -*- coding: utf-8 -*-
r"""
The last declined control, taken rule by rule.

WHAT IS AT STAKE

lab_direct_identity.py declines a control with an argument rather than
a pointer: "Z1 is an exact arithmetic identity whose reference is
itself, Z2 compares against S(N)N which is the reference, and Z3 and
Z4 compare two sums over the SAME terms, so a sign control would move
both sides together."

Remark {#rem:weightgapnull} fixed the criterion for judging such an
argument: a control is worth running exactly when the statistic stays
well conditioned under it.  Applied rule by rule the argument is
right in three places and wrong in one.

  Z1 is an identity -- nothing to control.
  Z3 and Z4 divide by the TOTAL, and for a coin the total is a
     square-root fluctuation near zero, so those ratios are ill
     conditioned and the decline is correct.
  Z2 divides by S(N)N, a fixed nonzero constant.  That is well
     conditioned, and the control was owed.

The same criterion also repairs Z3: measured against N instead of
against the total, its partial sums are well conditioned and can be
controlled.  Both are done here.

BACKS: Remark {#rem:identitynull} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The identity is mu's: |T - R|/R is under 1e-12 for mu, where
      T is the untruncated sum and R = sum_{p<N} Lambda(N-p) log p,
      and |T_eps|/R is under 0.2 for every draw.
  V2  Z2's ratio is mu's: T/(S(N)N) is within 0.02 of 1 for mu and
      outside [0.5, 1.5] for every draw.
  V3  Z3 repaired: measured against N, mu's |partial sum over
      k < N^0.90|/N exceeds every draw's at every N.
  V4  And its sign is mu's: mu's partial is negative at every N,
      while at most 4 of the 8 draws are negative at every N.

REFUTATION RULE (fixed before the run)

  V1  REFUTED if mu misses 1e-12, or if any draw reaches 0.2. The
      first would mean the identity is wrong; the second that the
      untruncated sum is not mu's.
  V2  REFUTED if mu leaves [0.98, 1.02] or any draw enters
      [0.5, 1.5]. This is the rule the decline got wrong.
  V3  REFUTED if any draw reaches mu at any N.
  V4  REFUTED if 5 or more draws are negative throughout, which would
      make the one-signed mass an artefact of the range and not of mu.

  All four gate.

  NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere, with
  the field, the weight log k and the k-range identical. Eight draws.
  Z3 and Z4 as lab_direct_identity.py states them are NOT controlled
  here and should not be: dividing by a total that a coin sends to
  zero is exactly the ill-conditioned case of [rem:weightgapnull].
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
OUT = os.path.join(ROOT, "results", "audit_directidentity_null.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
THETAS = [0.56, 0.70, 0.90]
CLIM = 4_000_000
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
    isp = np.zeros(NMAX + 1, dtype=bool)
    isp[pr] = True
    sqf = mu != 0

    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    rng = np.random.default_rng(SEED)
    signs = [mu.astype(np.float64)]
    for _ in range(COINS):
        c = np.zeros(NMAX + 1, dtype=np.float64)
        c[sqf] = rng.choice([-1.0, 1.0], size=int(sqf.sum()))
        signs.append(c)
    say("  1 mu and %d coins on supp(mu^2), field and k-range identical"
        % COINS)

    res = []
    for N in NS:
        S = twin
        for q in sorted(factor_set(N)):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        ks = np.flatnonzero(sqf[2:N]).astype(np.int64) + 2
        lg = np.log(ks.astype(np.float64))
        R = float((lam[1:N] * lam[N - 1:0:-1] * isp[N - 1:0:-1]).sum())
        idx = np.arange(1, N, dtype=np.int64)
        tot, part = [], []
        for sg in signs:
            f0 = np.zeros(N, dtype=np.float64)
            f0[1:] = lam[1:N] * sg[N - idx]
            contrib = np.empty(ks.size, dtype=np.float64)
            for i, k in enumerate(ks):
                k = int(k)
                r = N % k
                a = f0[r::k].sum() if r else f0[k::k].sum()
                contrib[i] = lg[i] * sg[k] * a
            del f0
            tot.append(float(contrib.sum()))
            row = []
            for th in THETAS:
                j = int(np.searchsorted(ks, int(N ** th)))
                row.append(float(contrib[:j].sum()))
            part.append(row)
        res.append((N, S, R, tot, part))
        say("  N = %-10d  #k = %-9d target = %.4f N" % (N, ks.size, R / N))

    say()
    say("V1  the untruncated sum against the Goldbach count")
    say("  N            mu |T-R|/R   coin |T|/R min   median      max")
    v1 = True
    for N, S, R, tot, part in res:
        d = abs(tot[0] - R) / R
        c = [abs(t) / R for t in tot[1:]]
        if d >= 1e-12 or max(c) >= 0.2:
            v1 = False
        say("  %-12d %-12.3e %-16.4f %-11.4f %.4f"
            % (N, d, min(c), float(np.median(c)), max(c)))
    say("  V1 %s" % ("hold" if v1 else "REFUTED"))

    say()
    say("V2  T / (S(N)N), the ratio Z2 reports")
    say("  N            mu        coin min   coin median   coin max")
    v2 = True
    for N, S, R, tot, part in res:
        m = tot[0] / (S * N)
        c = [t / (S * N) for t in tot[1:]]
        if not (0.98 <= m <= 1.02):
            v2 = False
        if any(0.5 <= v <= 1.5 for v in c):
            v2 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, m, min(c), float(np.median(c)), max(c)))
    say("  V2 %s" % ("hold" if v2 else "REFUTED"))

    say()
    say("V3/V4  the partial sum over k < N^0.90, measured against N")
    say("  N            mu        coin min   coin median   coin max")
    v3 = True
    for N, S, R, tot, part in res:
        m = part[0][2] / N
        c = [p[2] / N for p in part[1:]]
        if max(abs(v) for v in c) >= abs(m):
            v3 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, m, min(c), float(np.median(c)), max(c)))
    say("  V3 mu's magnitude exceeds every draw   %s"
        % ("hold" if v3 else "REFUTED"))
    muneg = all(r[4][0][2] < 0 for r in res)
    nneg = sum(1 for j in range(1, COINS + 1)
               if all(r[4][j][2] < 0 for r in res))
    v4 = muneg and nneg <= 4
    say("  V4 mu negative throughout: %s; draws that are: %d of %d   %s"
        % (muneg, nneg, COINS, "hold" if v4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). Why Z3 and Z4 are left uncontrolled.")
    say("  As lab_direct_identity.py states them they divide by the")
    say("  TOTAL, and for a coin the total is what column two of V1")
    say("  shows -- a few percent of the count, wandering in sign. A")
    say("  ratio to that is the ill-conditioned case, so the decline is")
    say("  right for those two and the repair is to measure against N,")
    say("  which V3 does. The three truncations against N, for mu:")
    say("  N            k<N^0.56    k<N^0.70    k<N^0.90")
    for N, S, R, tot, part in res:
        say("  %-12d %-11.4f %-11.4f %.4f"
            % (N, part[0][0] / N, part[0][1] / N, part[0][2] / N))

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("the identity and its truncations are mu's, and the decline was "
        "right for two rules of four" if ok else "REFUTED")

    head = [
        "STATISTIC: the untruncated sum T = sum_k (log k) s(k) A_s(N;k)",
        "           for s = mu and 8 coins, against",
        "           R = sum_{p<N} Lambda(N-p) log p; the ratio T/(S(N)N);",
        "           and the partial sums over k < N^0.56, N^0.70, N^0.90",
        "           measured against N rather than against T.",
        "NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere,",
        "      field, weight log k and k-range identical, eight draws.",
        "      lab_direct_identity.py's Z3 and Z4 are deliberately NOT",
        "      controlled: as stated they divide by a total that a coin",
        "      sends near zero, which is the ill-conditioned case of",
        "      [rem:weightgapnull]. Measuring against N instead repairs",
        "      them, and that repaired form is what V3 controls.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k over every squarefree",
        "       2 <= k < N; S(N) from an Euler product at the fixed bound",
        "       4e6; seed 20260808.",
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
