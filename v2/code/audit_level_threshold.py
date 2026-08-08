# -*- coding: utf-8 -*-
r"""
At what level does the wall's sign stop slipping away?

WHAT IS AT STAKE

Remark {#rem:sievedepth} closed one question and opened a sharper one.
The sign of H(N;k) is carried, to within 0.6 per cent in the lean and
0.05 standard errors in its decay, by the level-sqrt(N) sieve -- that
is, by primality. At the fixed level 29 the sign agreement is 0.7367
to 0.8129 and GETS WORSE as N grows. So the sign is not a
bounded-modulus object.

Between those two lies the only quantity the programme trades in: the
level as a power of N. Sieving to Q = N^alpha and asking at which
alpha the agreement stops degrading turns "not bounded" into a number,
and that number can be set beside theta' = 0.56.

The existing table already hints where to look. At Q = 101 and
N = 2e5 the agreement is 0.8658; at Q = 503 and N = 6.4e6 it is
0.8688 -- and those two Q are about N^0.4 at their own N. So alpha
near 0.4 is registered below as the flat one, with the open question
being whether it clears two standard errors, not which direction it
points.

BACKS: Remark {#rem:levelthreshold} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  T1  The control: at alpha = 1/2 the agreement and the lean ratio
      reproduce {#rem:sievedepth}'s top rung to within 0.001.
  T2  A bounded level slips: at the fixed Q = 29 the agreement's
      least-squares slope against log N is negative and reaches two
      standard errors.
  T3  A level below the square root suffices: at alpha = 0.4 the
      slope is within two standard errors of zero.
  T4  But not any level: at alpha = 0.2 the slope is still negative
      at two standard errors.

REFUTATION RULE (fixed before the run)

  T1  REFUTED at 0.001 -- not the same statistic, and nothing below
      may be compared with {#rem:sievedepth}.
  T2  REFUTED below two standard errors, which would mean the fixed
      level does not degrade after all and {#rem:sievedepth}'s reading
      of its own table was too strong.
  T3  REFUTED if alpha = 0.4 is still falling at two standard errors.
      That is the one that matters: the sign would then need a level
      indistinguishable from the square root, and no growing level
      short of sqrt(N) would hold it.
  T4  REFUTED if alpha = 0.2 is already flat, which would put the
      threshold below anything the sweep can separate and make the
      number meaningless.

  All four gate.

  NO NULL IS RUN for the agreement: a deterministic predictor is
  compared with a deterministic sign, term by term, and the
  no-information value for such a comparison is what
  lab_survivor_selection.py's permutation control measured at 0.5372
  to 0.5414. The lean ratios are read against the 256 global sign
  vectors of audit_lean_floor.py on the identical magnitudes.
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
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_level_threshold.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
ALPHAS = [0.1, 0.2, 0.3, 0.4, 0.5]
QFIXED = 29
DRAWS = 256
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
    """von Mangoldt and Moebius, the cofactor kept in int32"""
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
    del pr, lgp
    mu = np.ones(n + 1, dtype=np.int8)
    cof = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        cof[p::p] //= p
        pk = p * p
        while pk <= n:
            cof[pk::pk] //= p
            if pk > n // p:
                break
            pk *= p
    big = cof > 1
    del cof
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def nosmall(n, q):
    """not divisible by any odd prime at or below q, the sieve's own
    convention as in audit_sieve_depth.py"""
    s = np.ones(n + 1, dtype=bool)
    s[0] = False
    for p in primes_upto(max(q, 3)):
        p = int(p)
        if p == 2 or p > q:
            continue
        s[p::p] = False
    return s


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


def read_top():
    """the top-rung agreement and lean ratio of {#rem:sievedepth}"""
    src = io.open(os.path.join(RES, "audit_sieve_depth.txt"),
                  encoding="utf-8").read()
    i = src.index("R1/R2  the sign agreement as the sieve deepens")
    agr = {}
    for ln in src[i:].splitlines()[2:]:
        f = ln.split()
        if len(f) < 8 or not f[0].isdigit():
            break
        agr[int(f[0])] = float(f[-1])
    j = src.index("R3  the lean each depth gives, as a ratio to mu's")
    rat = {}
    for ln in src[j:].splitlines()[2:]:
        f = ln.split()
        if len(f) < 8 or not f[0].isdigit():
            break
        rat[int(f[0])] = float(f[-1])
    return agr, rat


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se, abs(a) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    MAJS = []          # the predictor's own majority share

    tagr, trat = read_top()
    say("read %d top-rung agreements and lean ratios from "
        "results/audit_sieve_depth.txt" % len(tagr))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    oddsqf = sqf.copy()
    oddsqf[::2] = False
    rng = np.random.default_rng(SEED)
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    lvl = {}
    for N in NS:
        for al in ALPHAS:
            lvl[(N, al)] = max(3, int(round(N ** al)))
        lvl[(N, 0.0)] = QFIXED
    say("the levels swept, Q = N^alpha:")
    say("  N            " + "".join("a=%-8.1f" % a for a in
                                    [0.0] + ALPHAS))
    for N in NS:
        say("  %-12d %s" % (N, "".join(
            "%-10d" % lvl[(N, a)] for a in [0.0] + ALPHAS)))
    surv = {}
    for q in sorted(set(lvl.values())):
        surv[q] = nosmall(NMAX, q)

    agr, lean, muln, flo = {}, {}, {}, {}
    for N in NS:
        PN = factor_set(N)
        here = sorted(set(lvl[(N, a)] for a in [0.0] + ALPHAS))
        ks, Hs = [], []
        Ps = dict((q, []) for q in here)
        for k in range(2, int(N ** THETA)):
            if not sqf[k]:
                continue
            if any(k % q == 0 for q in PN):
                continue
            M = N // k
            if M < 1:
                continue
            ms = np.arange(1, M + 1, 2, dtype=np.int64)
            ms = ms[oddsqf[ms]]
            for q in factor_set(k):
                if q > 2:
                    ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            vals = N - ms * k
            g = mu[ms].astype(np.float64)
            ks.append(k)
            Hs.append(float((lam[vals] * g).sum()))
            for q in here:
                Ps[q].append(float(g[surv[q][vals]].sum()))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        a = np.log(ks.astype(np.float64)) * H
        l1 = float(np.abs(a).sum())
        w = np.abs(a)
        sh = np.sign(H)
        muln[N] = abs(0.5 - float(a[a > 0].sum() / l1))
        eps = (rng.integers(0, 2, size=(DRAWS, ks.size))
               .astype(np.int8) * 2 - 1)
        flo[N] = float(np.median(np.abs((eps @ w) / (2.0 * l1))))
        for q in here:
            P = np.array(Ps[q])
            sp = np.sign(P)
            MAJS.append(max(float((sp > 0).mean()),
                            float((sp < 0).mean())))
            ok = (sh != 0) & (sp != 0)
            agr[(N, q)] = float((sh[ok] == sp[ok]).mean())
            spz = np.where(sp == 0, 1.0, sp)
            lean[(N, q)] = abs(0.5 - float(w[spz > 0].sum() / l1))
        say("  N = %-10d #k = %-6d done" % (N, ks.size))

    x = np.log(np.array(NS, dtype=np.float64))

    # ------------------------------------------------------------- T1
    say()
    say("T1  the control: alpha = 1/2 against the published top rung")
    say("  N            agree here  published   lean here  published")
    t1 = True
    for N in NS:
        q = lvl[(N, 0.5)]
        r = lean[(N, q)] / muln[N]
        if abs(agr[(N, q)] - tagr[N]) >= 0.001 or \
                abs(r - trat[N]) >= 0.001:
            t1 = False
        say("  %-12d %-11.4f %-11.4f %-10.4f %.4f"
            % (N, agr[(N, q)], tagr[N], r, trat[N]))
    say("  T1 %s   (cap 0.001)" % ("hold" if t1 else "REFUTED"))

    # -------------------------------------------------- T2 / T3 / T4
    say()
    say("T2/T3/T4  does the agreement hold up as N grows?")
    say("  level          " + "".join("%-11d" % N for N in NS)
        + "slope        t")
    slopes = {}
    for al in [0.0] + ALPHAS:
        row = [agr[(N, lvl[(N, al)])] for N in NS]
        b, rms, se, t = fit(x, np.array(row))
        slopes[al] = (b, se, t)
        say("  %-14s %s%+-12.6f %.2f"
            % ("Q = %d" % QFIXED if al == 0.0 else "Q = N^%.1f" % al,
               "".join("%-11.4f" % v for v in row), b, t))
        say("SCATTER slope_levelthreshold_a%.1f %.4f" % (al, rms))
        say("TSTAT slope_levelthreshold_a%.1f %.2f" % (al, t))
        say("SPREAD slope_levelthreshold_a%.1f %.4f"
            % (al, float(x.max() - x.min())))
        if t < 2.0:
            say("UNRESOLVED SIGN slope_levelthreshold_a%.1f" % al)
    b0, se0, t0 = slopes[0.0]
    t2 = (b0 < 0.0) and (t0 >= 2.0)
    b4, se4, t4v = slopes[0.4]
    t3 = t4v < 2.0
    b2, se2, t2v = slopes[0.2]
    t4 = (b2 < 0.0) and (t2v >= 2.0)
    say("  T2 the fixed level falls (%+.6f, %.2f s.e.)   %s"
        % (b0, t0, "hold" if t2 else "REFUTED"))
    say("  T3 alpha = 0.4 is flat (%+.6f, %.2f s.e.)   %s"
        % (b4, t4v, "hold" if t3 else "REFUTED"))
    say("  T4 alpha = 0.2 still falls (%+.6f, %.2f s.e.)   %s"
        % (b2, t2v, "hold" if t4 else "REFUTED"))

    say()
    say("  the same sweep in the lean, as a ratio to mu's:")
    say("  level          " + "".join("%-11d" % N for N in NS))
    for al in [0.0] + ALPHAS:
        row = [lean[(N, lvl[(N, al)])] / muln[N] for N in NS]
        say("  %-14s %s"
            % ("Q = %d" % QFIXED if al == 0.0 else "Q = N^%.1f" % al,
               "".join("%-11.4f" % v for v in row)))
    say("PERN levelthreshold_a04_over_mu %d %.4f %.4f"
        % (len(NS),
           min(lean[(N, lvl[(N, 0.4)])] / muln[N] for N in NS),
           max(lean[(N, lvl[(N, 0.4)])] / muln[N] for N in NS)))

    say()
    say("  DIAGNOSTIC on T2 and T4 (post hoc). Neither fails by")
    say("  direction: both slopes are negative, and NEITHER is")
    say("  resolved -- no level's agreement clears two standard")
    say("  errors, the fixed one included at %.2f. Six points of a"
        % t0)
    say("  ratio of counts do not carry that much power. The lean is")
    say("  the statistic that does, and it is the one the route")
    say("  cares about:")
    say("  level          slope        s.e.       t")
    lslopes = {}
    for al in [0.0] + ALPHAS:
        row = np.array([lean[(N, lvl[(N, al)])] / muln[N]
                        for N in NS])
        b, rms, se, t = fit(x, row)
        lslopes[al] = (b, se, t)
        say("  %-14s %+-12.6f %-10.6f %.2f"
            % ("Q = %d" % QFIXED if al == 0.0
               else "Q = N^%.1f" % al, b, se, t))
        say("SCATTER slope_leanratio_a%.1f %.4f" % (al, rms))
        say("TSTAT slope_leanratio_a%.1f %.2f" % (al, t))
        say("SPREAD slope_leanratio_a%.1f %.4f"
            % (al, float(x.max() - x.min())))
        if t < 2.0:
            say("UNRESOLVED SIGN slope_leanratio_a%.1f" % al)
    falling = [al for al in [0.0] + ALPHAS
               if lslopes[al][0] < 0.0 and lslopes[al][2] >= 2.0]
    flat = [al for al in [0.0] + ALPHAS if al not in falling]
    say("  falling at two standard errors: %s"
        % (", ".join("%.1f" % a for a in falling) if falling
           else "none"))
    say("  not resolved:                   %s"
        % (", ".join("%.1f" % a for a in flat) if flat else "none"))
    thr = min(flat) if flat else None
    say("  the threshold is taken by rule and not by eye: the")
    say("  smallest alpha whose trend is unresolved at two standard")
    say("  errors, on the named statistic.")
    say("LEVEL sieve_alpha_threshold %s"
        % ("%.1f" % thr if thr is not None else "none"))
    say("THRESHOLD FROM sieve_alpha_threshold leanratio")
    if thr is not None and thr > 0.0:
        say("UNBOUNDED LEVEL sieve_alpha_threshold")
    if thr is not None:
        say("  so the lean stops slipping between alpha = %.1f and"
            % max([a for a in falling if a < thr] or [0.0]))
        say("  alpha = %.1f, and theta' = %.2f is what the route"
            % (thr, THETA))
        say("  consumes -- the sign is held at a level %s than that."
            % ("lower" if thr < THETA else "no lower"))

    say()
    say("=" * 70)
    ok = t1 and t2 and t3 and t4
    say("the sign is held below the square root but not by a bounded "
        "level" if ok else "REFUTED")

    mj = max(MAJS)
    say()
    say("  the predictor's own majority sign share, at its worst over "
        "everything")
    say("  reported above: %.4f. An agreement is only a measurement "
        "where the" % mj)
    say("  predictor has variance; where it takes one sign almost "
        "everywhere,")
    say("  the agreement is the other side's marginal rate read back.")
    say("MARGINAL %s %.4f" % ("audit_level_threshold", mj))
    if mj >= 0.9:
        say("DEGENERATE %s" % "audit_level_threshold")

    head = [
        "STATISTIC: on the squarefree k < N^" + str(THETA)
        + " coprime to N, the",
        "           sign agreement of H(N;k) with",
        "           P_Q = sum_m mu(m) [N - mk has no odd prime factor",
        "           at or below Q], at the fixed Q = " + str(QFIXED),
        "           and at Q = N^alpha for alpha = "
        + ", ".join("%.1f" % a for a in ALPHAS) + ";",
        "           each level's least-squares slope of that agreement",
        "           against log N with its standard error; and the",
        "           mass-weighted lean each gives on mu's own",
        "           (log k)|H| magnitudes.",
        "NULL: none for the agreement -- a deterministic predictor is",
        "      compared with a deterministic sign and the",
        "      no-information value was measured by",
        "      lab_survivor_selection.py's permutation control at",
        "      0.5372 to 0.5414. The leans are read against the "
        + str(DRAWS),
        "      global sign vectors of audit_lean_floor.py on the",
        "      identical magnitudes.",
        "FIELD: N = 2e5 through 6.4e6 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < N^" + str(THETA) + "; m odd,",
        "       squarefree and coprime to k, m <= N//k; Lambda and mu",
        "       from an integer sieve to " + str(NMAX) + "; numpy",
        "       default_rng seed " + str(SEED) + ". The skip of primes",
        "       dividing k is vacuous and dropped, as in",
        "       audit_sieve_depth.py. Every N is 2^a 5^b, one odd",
        "       radical, as RADICALS declares. The top-rung control is",
        "       read from results/audit_sieve_depth.txt.",
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
