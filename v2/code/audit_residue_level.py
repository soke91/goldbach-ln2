# -*- coding: utf-8 -*-
r"""
How far the level goes if only the residue has to be bounded.

WHAT IS AT STAKE

Remark {#rem:provablehalf} found that the elementary half is the half a
classical unconditional estimate reaches in shape: every condition
defining P's summation set is multiplicative or a residue condition to
a bounded modulus, so sum_{k<K}(log k)|P| would be o(N) for every fixed
theta' < 1 if the uniformity in k holds. Nothing of the kind exists for
R past k = N^{1/2}.

Grant that. Then |H| <= beta|P| + |R| gives B(N) <= B_R(N) + o(N), and
the route's condition [eq:nolog] becomes a condition on the residue
alone. The level it then permits is K*_R, and Remark {#rem:splitbudget}
already prints one: 9191 to 63399, at exponents 0.7477 down to 0.7382.

**But that K*_R is against the wrong budget.** {#rem:splitbudget}
crosses each half against S(N)N. The route needs S(N)(1-A(N))N, which
audit_model_transfer.py measured as a factor 4.7009 smaller and worth
about 0.1677 in the exponent. Subtracting that by hand from 0.74 lands
near 0.57 -- close enough to theta' = 0.56, and to 1/2, that the
difference between clearing the barrier and missing it is inside the
arithmetic. It has to be computed, not inferred.

BACKS: Remark {#rem:residuelevel} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  The control: crossing B_R against S(N)N reproduces the published
      K*_R -- 9191, 14289, 23131, 37439, 63399, read from
      results/lab_split_budget.txt -- to within 2 per cent at every N.
  U2  Against the operative budget S(N)(1-A(N))N the residue alone
      still clears the square-root barrier: log K*_R / log N exceeds
      0.5 at every N.
  U3  And it clears the theta' = 0.56 the papers use, at every N.
  U4  The margin is not closing: the exponent's least-squares slope
      against log N over the five N is not negative.

REFUTATION RULE (fixed before the run)

  U1  REFUTED at 2 per cent at any N, which would mean this is not the
      same split.
  U2  REFUTED if the exponent reaches 0.5 from above at any N. That is
      the one that matters: it would say the conditional reduction of
      {#rem:provablehalf} buys nothing at all, because the residue by
      itself cannot be carried past the barrier even at accessible N.
  U3  REFUTED if the exponent falls to 0.56 or below at any N.
  U4  REFUTED if the slope is negative, which would say the margin
      shrinks with N and the accessible range flatters the reduction.

  All four gate.

  NO NULL IS RUN and none applies. Nothing here is a detection against
  a background: a measured sum is crossed against a computed threshold
  and the crossing is located. The sign control for this field is
  lab_split_budget.py's size permutation, which established that the
  budget is spent at the bottom of the k-range, and
  lab_residue_cancellation.py's coin arm, which established that R's
  size is bought by cancellation at all.
"""

import io
import math
import os
import re
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_residue_level.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000
THETA = 0.56


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
    return lam, mu


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


def read_gap():
    """the K*_H budget gap in the exponent -- read, not copied"""
    p = os.path.join(ROOT, "results", "audit_model_transfer.txt")
    src = io.open(p, encoding="utf-8").read()
    return float(re.search(r"mean gap ([\d.]+)", src).group(1))


def read_published():
    """K*_R against S(N)N, and the exponents -- read, not copied"""
    p = os.path.join(ROOT, "results", "lab_split_budget.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("N            K*_H      K*_P      K*_R")
    ks = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        ks[int(f[0])] = int(f[3])
    j = src.index("log K*_H/log N")
    ex = {}
    for ln in src[j:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        ex[int(f[0])] = float(f[3])
    return ks, ex


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubk, pubex = read_published()
    gapH = read_gap()
    say("read %d published K*_R and %d exponents from "
        "results/lab_split_budget.txt, and the K*_H budget gap"
        % (len(pubk), len(pubex)))
    say("  %.4f from results/audit_model_transfer.txt" % gapH)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = sieves(NMAX)
    sqf = mu != 0
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("  sieve weight over the odd primes %s"
        % ", ".join(map(str, QS)))

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    res = []
    for N in NS:
        PN = factor_set(N)
        A_, S_ = artin, twin
        for q in sorted(PN):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))

        ks, Hs, Ps = [], [], []
        for k in range(2, KCAP):
            if not sqf[k] or any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
            if M < 2:
                continue
            ms = np.arange(1, M + 1, 2, dtype=np.int64)
            ms = ms[sqf[ms]]
            for q in factor_set(k):
                if q > 2:
                    ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            vals = N - ms * k
            g = mu[ms].astype(np.float64)
            w = np.ones(ms.size, dtype=np.float64)
            for q in QS:
                if k % q == 0:
                    continue
                w *= np.where(vals % q == 0, 0.0, q / (q - 1.0))
            ks.append(k)
            Hs.append(float((lam[vals] * g).sum()))
            Ps.append(float((g * w).sum()))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        P = np.array(Ps)
        beta = float((H * P).sum() / (P * P).sum())
        R = H - beta * P
        cum = np.cumsum(np.log(ks.astype(float)) * np.abs(R))
        res.append((N, ks, cum, S_, A_, beta))
        say("  N = %-10d  #k = %-7d beta = %.6f  S = %.4f  A = %.4f"
            % (N, ks.size, beta, S_, A_))

    def cross(cum, ks, thr):
        j = int(np.searchsorted(cum, thr))
        if j >= ks.size:
            return None
        return int(ks[j])

    # ------------------------------------------------------------- U1
    say()
    say("U1  the control: B_R crossed against S(N)N")
    say("  N            K*_R (here)   published    ratio")
    u1 = True
    for N, ks, cum, S_, A_, beta in res:
        k1 = cross(cum, ks, S_ * N)
        r = k1 / pubk[N] if k1 else float("nan")
        if not (abs(r - 1.0) < 0.02):
            u1 = False
        say("  %-12d %-13s %-12d %.4f"
            % (N, str(k1), pubk[N], r))
    say("  U1 %s   (cap 2 per cent)" % ("hold" if u1 else "REFUTED"))

    # ---------------------------------------------------------- U2/U3
    say()
    say("U2/U3  against the budget the route actually needs,")
    say("  S(N)(1-A(N))N, which is smaller by 1/(1-A) :")
    say("  N            budget factor   K*_R      log K*_R/log N  "
        "clears .56")
    u2 = u3 = True
    ex = []
    for N, ks, cum, S_, A_, beta in res:
        k2 = cross(cum, ks, S_ * (1.0 - A_) * N)
        if k2 is None:
            u2 = u3 = False
            say("  %-12d %-15.4f %-9s -          -"
                % (N, 1.0 / (1.0 - A_), "none"))
            continue
        e = math.log(k2) / math.log(N)
        ex.append(e)
        if e <= 0.5:
            u2 = False
        if e <= THETA:
            u3 = False
        say("  %-12d %-15.4f %-9d %-10.4f %s"
            % (N, 1.0 / (1.0 - A_), k2, e, "yes" if e > THETA else "NO"))
    say("  U2 exponent above 0.5 at every N   %s"
        % ("hold" if u2 else "REFUTED"))
    say("  U3 exponent above %.2f at every N   %s"
        % (THETA, "hold" if u3 else "REFUTED"))
    say()
    say("  the arithmetic this sweep covers, which gate check G34")
    say("  reads. Every N here is 2^a 5^b, so the sweep has ONE odd")
    say("  radical and the exponents above are a statement about that")
    say("  radical. audit_residue_arithmetic.py repeats the measurement")
    say("  across seven types and finds it below 1/2 at the")
    say("  primorial-like ones:")
    rads = set()
    for N, ks, cum, S_, A_, beta in res:
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  %d N, %d distinct odd radical%s: %s"
        % (len(res), len(rads), "" if len(rads) == 1 else "s",
           ", ".join(str(r) for r in sorted(rads))))
    say("RADICALS %d" % len(rads))

    say()
    say("  the two budget constants, in the form the gate reads, so")
    say("  that no exponent from one is ever read against the other:")
    for N, ks, cum, S_, A_, beta in res:
        say("BUDGET kstar_R_SN_N%d %.6f" % (N, S_))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, S_ * (1.0 - A_)))

    # ------------------------------------------------------------- U4
    say()
    say("U4  is the margin closing?")
    if len(ex) == len(NS):
        x = np.log(np.array(NS, dtype=float))
        sl = float(np.polyfit(x, np.array(ex), 1)[0])
        u4 = sl >= 0.0
        say("  exponents %s" % ", ".join("%.4f" % v for v in ex))
        say("  least-squares slope against log N = %+.6f   %s"
            % (sl, "hold" if u4 else "REFUTED"))
        f = [float(np.polyfit(x[s], np.array(ex)[s], 1)[0])
             for s in (slice(None), slice(1, None), slice(0, -1))]
        say("  leave-one-out: %+.6f, %+.6f, %+.6f -- spread %.6f"
            % (f[0], f[1], f[2], max(f) - min(f)))
    else:
        u4 = False
        say("  not enough crossings to fit   REFUTED")
    say("  U4 %s" % ("hold" if u4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). The two budgets side by side, and")
    say("  what the published exponent would have suggested. The gap")
    say("  audit_model_transfer.py measured between the two budgets")
    say("  for K*_H was %.4f in the exponent; here it is:" % gapH)
    say("  N            at S N     at S(1-A)N   gap    published K*_R exp")
    for i, (N, ks, cum, S_, A_, beta) in enumerate(res):
        k1 = cross(cum, ks, S_ * N)
        if k1 is None or i >= len(ex):
            continue
        e1 = math.log(k1) / math.log(N)
        say("  %-12d %-10.4f %-12.4f %-6.4f %.4f"
            % (N, e1, ex[i], e1 - ex[i], pubex[N]))
    say("  So the headline 0.74 is against a budget four and a half")
    say("  times too generous, and the number the route needs is the")
    say("  second column.")

    say()
    say("=" * 70)
    ok = u1 and u2 and u3 and u4
    say("the residue alone carries the level past theta' = 0.56"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           B_R(N;K) = sum_{k<K}(log k)|R(N;k)| first reaches a",
        "           budget, for the budget S(N)N used by",
        "           lab_split_budget.py and for the budget",
        "           S(N)(1-A(N))N that Proposition [prop:nolog] needs;",
        "           the exponent log K*_R / log N of the second; and its",
        "           least-squares slope against log N.",
        "NULL: none is run and none applies. A measured sum is crossed",
        "      against a computed threshold and the crossing located;",
        "      there is no background to detect against. The sign",
        "      controls for this field were run in lab_split_budget.py,",
        "      whose size permutation established where the budget is",
        "      spent, and in lab_residue_cancellation.py, whose coin arm",
        "      established that R's size is bought by cancellation.",
        "FIELD: N = 2e5 to 3.2e6 by doubling; k squarefree and coprime",
        "       to N with 2 <= k < 100000; m odd, squarefree, coprime to",
        "       the odd part of k, m < N/k; the sieve weight uses the",
        "       odd primes up to 30; beta refitted here as",
        "       sum(H P)/sum(P^2) on the same k-range; S(N) and A(N)",
        "       from Euler products at the fixed bound 4000000; the",
        "       published K*_R and exponents are read from",
        "       results/lab_split_budget.txt.",
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
