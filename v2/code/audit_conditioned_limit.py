# -*- coding: utf-8 -*-
r"""
A conditioning floor that is derived, and the field that survives it

WHAT IS AT STAKE

rem:limitsqueeze refuted four counts at once and located the fault in
the ladder rather than the reading: for prime p the inner variable runs
m in [K, N/p) with K = floor(N^theta), so the largest primes at the
smallest N have almost no m, and PSET's rule -- "contributes at all
eleven N" -- is exactly the rule that puts them there.  Its first rung
scattered g from 0.607612 to 1.243394 while its last is smooth and
monotone in p.

A floor for that is not chosen, it is demanded.  The sum over m has
about n terms of which about n/log N are prime hits, so its relative
fluctuation is about sqrt(log N / n), and resolving a difference of
RESOL in g needs

    n  >=  log N / RESOL^2 .

At RESOL = 0.01 that is 10^4 log N, and **it is severe**: it kills
every rung rem:limitsqueeze used below the millions.  So this run
climbs instead of widening -- six rungs from 6.4e6 to 204.8e6 -- and
asks the floor which rectangle of (p, N) it leaves.

rem:limitsqueeze's LL3 was refuted and stays refuted.  Its eleven
medians turned near 0.8858 and rose over the last four, an observation
made after the rule broke and therefore worth nothing until it is
registered.  **MM3 below is that registration, and it can break.**

BACKS: Remark {#rem:conditionedlimit} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  MM1 THE GATE.  At N = 25600000, which both ladders share, g(3),
      g(11), g(17) and g(19) reproduce rem:derivedlimit's POINT
      gtop_3 and gtop_11 and rem:limitdirection's POINT glast_17 and
      POINT glast_19 to six decimals.
  MM2 **THE FLOOR LEAVES A FIELD.**  The largest rectangle every cell
      of which clears n >= log N / RESOL^2 has at least 4 primes and
      at least 4 rungs.
  MM3 **THE CENTRE CLIMBS.**  On that rectangle the median of g rises
      at every rung step.
  MM4 **THE SPAN CLOSES.**  On that rectangle the span max(g) - min(g)
      falls at every rung step.
  MM5 THE p-SHAPE IS NOT THE CONDITIONING.  At the top rung, g is
      still strictly decreasing in p across the rectangle's primes.

REFUTATION RULE (fixed before the run)

  MM1 REFUTED outside six decimals on any of the four; nothing below
      is reported.
  MM2 REFUTED below 4 by 4.  **Then this question is not answerable
      on this machine** and the run says so: the derived floor, not
      taste, is what removes the field, and the honest report is that
      eq:derivedlimit cannot be tested at the resolution it needs
      until N is far larger.  Nothing below is then read as evidence.
  MM3 REFUTED if the median falls at any rung step.  With at most 5
      steps a two-sided coin tail of 5/5 is 0.062500 and of 4/5 is
      0.375000, so **only a clean sweep is a count here** and the
      rule demands one; a single fall means rem:limitsqueeze's
      post-hoc turn was noise and eq:derivedlimit loses its last
      support from motion.
  MM4 REFUTED if the span rises at any rung step, on the same
      arithmetic.
  MM5 REFUTED at 2 or more inversions.  Then the descent of g in p
      is partly the conditioning that rem:limitsqueeze exposed, and
      the residual eq:derivedlimit fails to explain is smaller than
      rem:derivedlimit's twenty-prime table made it look.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  Five rung steps is
  the whole budget: a clean sweep has a two-sided tail of 0.062500,
  which does not clear 0.05, and four rules are asked of the same six
  rungs.  **So even MM3 and MM4 both sweeping is a lean and not a
  count**, and must be written as one.  The rectangle is fixed by the
  floor and a stated tie-break -- larger area, then more rungs -- so
  that it is not chosen after g is seen; if two rectangles tie on both
  the run prints both and reports the smaller-prime one, and that
  tie is a defect of the rule, not a choice.

  WHAT THIS CANNOT DO.  Climbing to 204800000 buys five doublings and
  rem:derivedlimit measured g(3) moving 0.089063 over ten; the centre
  sits about 0.10 below 1.  **Nothing here can reach the limit** --
  it can only say whether the motion toward it survives a floor that
  the last run showed was missing.  And a correction to a model that
  rem:residuemodel refuted as a whole is still not a route to the
  exponent.
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
OUT = os.path.join(ROOT, "results", "audit_conditioned_limit.txt")
SRCD = os.path.join(ROOT, "results", "audit_derived_limit.txt")
SRCL = os.path.join(ROOT, "results", "audit_limit_direction.txt")

THETA = 0.56
NS = [6_400_000, 12_800_000, 25_600_000, 51_200_000,
      102_400_000, 204_800_000]
CAND = (3, 7, 11, 13, 17, 19, 23, 29, 31, 37)
NGATE = 25_600_000
RESOL = 0.01
DEC = 6
MINP = 4
MINN = 4
INVMAX = 2


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
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


def phi(n):
    r = n
    for p in factor_set(n):
        r = r // p * (p - 1)
    return r


def gvals(N, lam, sqf):
    """g(p, N) for CAND, with the surviving term count of each.

    The m-ranges are contiguous, so every mask is a slice of the
    sieve and a strided assignment -- no index array the size of N
    is ever built.
    """
    PN = factor_set(N)
    K = int(N ** THETA)
    D = (N - 1) // K

    def survivors(lo, hi, drop):
        keep = sqf[lo:hi].copy()
        for q in drop:
            keep[((-lo) % int(q))::int(q)] = False
        return np.flatnonzero(keep).astype(np.int64) + lo

    keep = sqf[K:N].copy()
    keep &= lam[K:N] == 0.0
    for q in PN:
        keep[((-K) % int(q))::int(q)] = False
    m1 = np.flatnonzero(keep).astype(np.int64) + K
    del keep
    a = float((lam[N - m1] * np.log(m1.astype(np.float64))).sum())
    del m1

    g, cnt = {}, {}
    for p in CAND:
        if p > D or p in PN:
            g[p], cnt[p] = float("nan"), 0
            continue
        mp = survivors(K, (N - 1) // p + 1, {p} | PN)
        cnt[p] = int(mp.size)
        v = float((lam[N - p * mp]
                   * np.log(mp.astype(np.float64))).sum())
        w = (1.0 / phi(p)) * (1.0 - p / D) / (1.0 - 1.0 / D)
        g[p] = (v / a) / w * (p + 1.0) / p
        del mp
    return g, cnt, D, K


def read_pub():
    d = io.open(SRCD, encoding="utf-8").read()
    l = io.open(SRCL, encoding="utf-8").read()
    out = {}
    for src, key in ((d, "gtop_3"), (d, "gtop_11"),
                     (l, "glast_17"), (l, "glast_19")):
        m = re.search(r"^POINT %s (\S+)\s*$" % key, src, re.M)
        if not m:
            raise SystemExit("a published value is missing: " + key)
        out[key] = m.group(1)
    return out


def tail(k, tot):
    hi = max(k, tot - k)
    return 2.0 * sum(math.comb(tot, j)
                     for j in range(hi, tot + 1)) / 2.0 ** tot


HEAD = [
    "STATISTIC: the median and the span of g(p, N) = c(p, N)(p+1)/p",
    "           over the largest rectangle of (p, N) every cell of",
    "           which clears a term-count floor derived from the",
    "           resolution demanded, and the order of g in p there.",
    "FIELD: N = %s;" % NS,
    "       p over CAND = %s," % (CAND,),
    "       K = floor(N^%.2f), D = floor((N-1)/K); a cell is kept"
    % THETA,
    "       when its surviving term count n satisfies",
    "       n >= log N / %.2f^2. Four values at N = %d are READ from"
    % (RESOL, NGATE),
    "       results/audit_derived_limit.txt and",
    "       results/audit_limit_direction.txt.",
    "DERIVED: the sum over m has about n terms of which about n/log N",
    "         are prime hits, so its relative fluctuation is about",
    "         sqrt(log N / n); resolving RESOL in g needs n >=",
    "         log N / RESOL^2. The floor is demanded, not chosen.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    say("READ audit_derived_limit.txt POINT gtop_3 %s" % pub["gtop_3"])
    say("READ audit_derived_limit.txt POINT gtop_11 %s"
        % pub["gtop_11"])
    say("READ audit_limit_direction.txt POINT glast_17 %s"
        % pub["glast_17"])
    say("READ audit_limit_direction.txt POINT glast_19 %s"
        % pub["glast_19"])
    say("PRINTBOUND audit_conditioned_limit %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, resolution %.2f, rectangle at least %d primes"
        % (THETA, RESOL, MINP))
    say("  by %d rungs, at most %d inversions allowed" % (MINN, INVMAX))
    say("  candidates %s" % (CAND,))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    del mu

    G, CNT, FLOOR = {}, {}, {}
    for N in NS:
        g, cnt, D, K = gvals(N, lam, sqf)
        G[N], CNT[N] = g, cnt
        FLOOR[N] = math.log(N) / RESOL ** 2
        say("  N = %-10d K = %-6d D = %-5d floor %d"
            % (N, K, D, int(FLOOR[N])))
    say("SCALES %d" % len(NS))

    # ------------------------------------------------------------- MM1
    say()
    say("MM1  the gate")
    checks = [("g(3)", G[NGATE][3], pub["gtop_3"]),
              ("g(11)", G[NGATE][11], pub["gtop_11"]),
              ("g(17)", G[NGATE][17], pub["glast_17"]),
              ("g(19)", G[NGATE][19], pub["glast_19"])]
    mm1 = True
    for lab, got, want in checks:
        ok = abs(round(got, DEC) - float(want)) < 10.0 ** (-DEC)
        mm1 &= ok
        say("  %-6s %.6f against its %s   %s"
            % (lab, got, want, "ok" if ok else "MISMATCH"))
    say("  MM1 %s" % ("hold" if mm1 else "REFUTED"))
    if not mm1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------------------- MM2
    say()
    say("MM2  what does the derived floor leave?")
    say("      p " + "".join("%12d" % N for N in NS))
    for p in CAND:
        say("  %5d " % p
            + "".join("%11d%s" % (CNT[N][p],
                                  "*" if CNT[N][p] >= FLOOR[N] else " ")
                      for N in NS))
    say("  (* = clears the floor for that N)")

    best = None
    for i in range(len(NS)):
        rungs = NS[i:]
        for j in range(1, len(CAND) + 1):
            ps = CAND[:j]
            if all(CNT[N][p] >= FLOOR[N] and CNT[N][p] > 0
                   for N in rungs for p in ps):
                cand = (len(ps) * len(rungs), len(rungs), ps, rungs)
                if best is None or cand[:2] > best[:2]:
                    best = cand
    if best is None:
        area, ps, rungs = 0, (), ()
    else:
        area, _, ps, rungs = best
    mm2 = len(ps) >= MINP and len(rungs) >= MINN
    say("  rectangle %d primes by %d rungs, area %d"
        % (len(ps), len(rungs), area))
    say("  primes %s" % (ps,))
    say("  rungs  %s" % (rungs,))
    say("POINT rectprimes %d" % len(ps))
    say("POINT rectrungs %d" % len(rungs))
    say("  MM2 %s   (floor: %d by %d)"
        % ("hold" if mm2 else "REFUTED", MINP, MINN))
    if not mm2:
        say()
        say("=" * 70)
        say("MM1 hold  MM2 REFUTED")
        say()
        say("the derived floor removes the field. this question is "
            "not answerable")
        say("on this machine at the resolution it needs, and nothing "
            "below the")
        say("floor is read as evidence -- which is what the rule said "
            "in advance.")
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        print("\nwrote %s" % os.path.normpath(OUT))
        return 0

    med = [float(np.median([G[N][p] for p in ps])) for N in rungs]
    spn = [float(max(G[N][p] for p in ps) - min(G[N][p] for p in ps))
           for N in rungs]

    # ------------------------------------------------------------- MM3
    say()
    say("MM3  does the centre climb on the rectangle?")
    up = sum(1 for i in range(len(med) - 1) if med[i + 1] > med[i])
    tot = len(med) - 1
    mm3 = up == tot
    say("  median %s" % " ".join("%.6f" % v for v in med))
    say("  rises %d/%d" % (up, tot))
    say("POINT medbottom %.6f" % med[0])
    say("POINT medtop %.6f" % med[-1])
    say("POINT medrises %d" % up)
    say("POINT medmotion %.6f" % (med[-1] - med[0]))
    say("POINT medgap %.6f" % (1.0 - med[-1]))
    say("  MM3 %s" % ("hold" if mm3 else "REFUTED"))

    # ------------------------------------------------------------- MM4
    say()
    say("MM4  does the span close on the rectangle?")
    dn = sum(1 for i in range(len(spn) - 1) if spn[i + 1] < spn[i])
    mm4 = dn == tot
    say("  span   %s" % " ".join("%.6f" % v for v in spn))
    say("  falls %d/%d" % (dn, tot))
    say("POINT spanbottom %.6f" % spn[0])
    say("POINT spantop %.6f" % spn[-1])
    say("POINT spanfallsc %d" % dn)
    say("POINT spanmotion %.6f" % (spn[0] - spn[-1]))
    say("  MM4 %s" % ("hold" if mm4 else "REFUTED"))

    # ------------------------------------------------------------- MM5
    say()
    say("MM5  is the descent in p still there at the top rung?")
    NT = rungs[-1]
    vs = [G[NT][p] for p in ps]
    inv = sum(1 for i in range(len(vs) - 1) if vs[i + 1] >= vs[i])
    mm5 = inv < INVMAX
    say("      p     g at the top rung   p/(p+1)   g - 1")
    for p in ps:
        say("  %5d   %.6f            %.6f  %+.6f"
            % (p, G[NT][p], p / (p + 1.0), G[NT][p] - 1.0))
    say("  inversions %d of %d steps" % (inv, len(vs) - 1))
    say("POINT inversions %d" % inv)
    say("  MM5 %s   (cap: %d)"
        % ("hold" if mm5 else "REFUTED", INVMAX))

    say()
    say("=" * 70)
    say("MM1 %s  MM2 %s  MM3 %s  MM4 %s  MM5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (mm1, mm2, mm3, mm4, mm5)))
    say()
    if mm3 and mm4:
        say("above a floor that is demanded rather than chosen, the "
            "centre climbs")
        say("and the span closes at every rung. five steps is the "
            "whole budget and")
        say("a clean sweep tails at %.6f, so this is a lean and not "
            "a count --" % tail(tot, tot))
        say("the rule said so before the run.")
    else:
        say("above the derived floor the motion does not survive. "
            "rem:limitsqueeze's")
        say("post-hoc turn was noise, and eq:derivedlimit keeps its "
            "derivation and")
        say("loses its evidence from motion.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
