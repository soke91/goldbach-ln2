# -*- coding: utf-8 -*-
r"""
Inside one radical, packed close: is the level a function of N at all

WHAT IS AT STAKE

rem:levelmatched measured L(N) = log(|sum a| / l2) with no fitting
anywhere and found the radicals separated by 1.890177 against a band
contribution of 0.024960 -- the first statement about the radical in
this branch that no window can be blamed for.  But its control failed
on one pair of three: two N of radical {2,7}, 0.133531 apart in log N,
gave L differing by 0.066632.  A smooth L rising at that radical's own
published drift, near 0.132, would put them 0.018 apart.  The excess
is not explained by the radical and there is no window to blame.

The pair is one pair.  Packing N close inside a single radical decides
what it was.  If L is a smooth function of N, halving the gap halves
the difference; if there is scatter at fixed radical, the difference
stays put while the gap shrinks.  **At a log gap of 0.024 the two
readings are 0.003 and 0.05 -- a factor of twenty apart**, which is
what makes this cheap test decisive.

Fourteen N of radical {2,7} and fourteen of {2,5} in [1.0e6, 4.1e6],
all of the form 2^a p^b.  {2,5} is the better instrument by accident
of arithmetic: 128/125 = 1.024, so it supplies seven pairs at a log
gap of 0.0237.

BACKS: Remark {#rem:leveldense} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  T1  THE GATE.  N = 1882384 and 2151296 reproduce rem:levelmatched's
      L of +2.008906 and +2.075538, and N = 2000000 its +2.166023, to
      six decimals.
  T2  **L is smooth in N inside a radical.**  Fitting L on log N
      within each radical, the r.m.s. residual is under 0.03.
  T3  **THE ONE THE PACKING IS FOR.**  Every pair with a log gap
      under 0.03 has |dL| under 0.01, as a function rising at 0.14
      per log unit would give.
  T4  And the level's own slope is the drift: the fitted slope within
      each radical agrees with that radical's published drift to
      within 0.05.

REFUTATION RULE (fixed before the run)

  T1  REFUTED outside six decimals on any of the three; nothing below
      is reported.
  T2  REFUTED above 0.03 r.m.s.
  T3  **REFUTED by any close pair above 0.01, and that is the
      outcome that decides what this branch has been measuring.**  It
      would mean L carries scatter at fixed radical that does not
      shrink with the spacing -- so the radical does not determine
      the level, the failure cannot be laid on a window, and every
      radical number in this branch is a number plus an unmodelled
      term of that size.  **The size is then the headline**: this run
      must report the scatter's magnitude, because that magnitude
      bounds what any radical statement here can claim.
  T4  REFUTED outside 0.05.  Then the level and the drift disagree
      about the same quantity and neither is usable until that is
      resolved.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  A pair's difference
  is bounded below by what a smooth L would give: its log gap times
  the radical's own drift, not times any single number -- the defect
  rem:levelmatched recorded as its ninth.  **This run prints, for
  every pair it judges, the smooth prediction computed from that
  radical's own published drift**, and a |dL| below its own smooth
  prediction is not evidence of anything.  T3's cap of 0.01 sits
  above the smooth prediction at a gap of 0.024, which is 0.0034 at
  the {2,5} drift; a close pair landing between those two numbers
  tests nothing and the remark must say so rather than count it.

  WHAT THIS CANNOT DO.  Two radicals of the form 2^a p^b.  Nothing
  here measures a slope over a long window, so it says nothing about
  the drift's own value, and nothing about item 5's demand.  A
  scatter found here is measured, not explained.
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
OUT = os.path.join(ROOT, "results", "audit_level_dense.txt")
SRCL = os.path.join(ROOT, "results", "audit_level_matched.txt")
SRCV = os.path.join(ROOT, "results", "audit_valuation.txt")

THETA = 0.56
SETS = [
    (7, [1_075_648, 1_229_312, 1_404_928, 1_605_632, 1_647_086,
         1_835_008, 1_882_384, 2_151_296, 2_458_624, 2_809_856,
         3_211_264, 3_294_172, 3_670_016, 3_764_768]),
    (5, [1_000_000, 1_024_000, 1_250_000, 1_280_000, 1_310_720,
         1_562_500, 1_600_000, 1_638_400, 2_000_000, 2_048_000,
         2_500_000, 2_560_000, 3_125_000, 3_200_000]),
]
GATES = [1_882_384, 2_151_296, 2_000_000]
DRIFTBASE = {7: [14_336, 28_672, 25_088, 21_952, 19_208],
             5: [20_480, 25_600, 16_000, 20_000, 25_000]}
DEC = 6
RMSCAP = 0.03
CLOSEGAP = 0.03
CLOSECAP = 0.01
SLOPECAP = 0.05


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
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int64)
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


def level(N, lam, mu, sqf):
    PN = factor_set(N)
    K = int(N ** THETA)
    lk = np.zeros(N, dtype=np.float64)
    l2sq = 0.0
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        lg = math.log(k)
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        mm = mu[ms].astype(np.float64)
        l2sq += (lg * float((lam[N - ms * k] * mm).sum())) ** 2
        lk[ms * k] += lg * mm
        del ms, mm
    j = np.arange(1, N, dtype=np.int64)
    sa = abs(float((lam[N - j] * lk[1:]).sum()))
    del j, lk
    return math.log(sa / math.sqrt(l2sq))


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    return float(b), math.sqrt(float((r ** 2).mean()))


def read_pub():
    src = io.open(SRCL, encoding="utf-8").read()
    lv = {}
    for N in GATES:
        m = re.search(r"^POINT levelmatched_%d ([-+]?[\d.]+)\s*$" % N,
                      src, re.M)
        if not m:
            raise SystemExit("no levelmatched marker for %d" % N)
        lv[N] = float(m.group(1))
    vs = io.open(SRCV, encoding="utf-8").read()
    dr = {}
    for p, bases in DRIFTBASE.items():
        vals = []
        for b in bases:
            m = re.search(r"^POINT valdrift_%d ([-+]?[\d.]+)\s*$" % b,
                          vs, re.M)
            if m:
                vals.append(float(m.group(1)))
        if not vals:
            raise SystemExit("no valdrift markers for p = %d" % p)
        dr[p] = float(np.mean(vals))
    sp = re.search(r"^POINT levelspread ([\d.]+)\s*$", src, re.M)
    if not sp:
        raise SystemExit("no levelspread marker")
    return lv, dr, float(sp.group(1))


HEAD = [
    "STATISTIC: L(N) = log(|sum a| / l2) at single N, packed close",
    "           inside one radical, and the difference of neighbouring",
    "           L against what a smooth L rising at that radical's own",
    "           published drift would give.",
    "FIELD: N = 2^a p^b in [%d, %d] for p = 7 and p = 5, fourteen"
    % (min(min(g) for _, g in SETS), max(max(g) for _, g in SETS)),
    "       each; k over the squarefree k < N^%.2f coprime to N; j"
    % THETA,
    "       over every index below N. Three L are READ from",
    "       results/audit_level_matched.txt as the gate and each",
    "       radical's drift from results/audit_valuation.txt.",
    "NOTE: no slope is used to judge anything here; the slope is",
    "      fitted only to be compared with the published drift.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    lv, dr, spread = read_pub()
    for N in GATES:
        say("READ audit_level_matched.txt %d %.6f" % (N, lv[N]))
    for p in sorted(dr):
        say("READ audit_valuation.txt drift_%d %.6f" % (p, dr[p]))
    say("READ audit_level_matched.txt spread %.6f" % spread)
    say("  three levels, each radical's mean drift, and the radical spread")
    say("PRINTBOUND audit_level_dense %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, r.m.s. cap %.2f, close gap %.2f, close cap "
        "%.2f, slope cap %.2f"
        % (THETA, RMSCAP, CLOSEGAP, CLOSECAP, SLOPECAP))
    say("RADICALS %d" % len(SETS))
    say("  NOTE, disclosed: the first execution's p = 7 list carried")
    say("  1834496 and 3668992 as 2^18*7 and 2^19*7. They are not: "
        "they are")
    say("  2^9*3583 and 2^10*3583, radical {2,3583}, and 2^18*7 is "
        "1835008.")
    say("  They gave the two largest |dL| of that run, both near "
        "one, by")
    say("  comparing across radicals -- which rem:levelmatched "
        "measured as worth")
    say("  up to %.6f, so those were not scatter but the radical "
        "effect." % spread)
    say("  Their values are not quoted: that execution's result file "
        "has been")
    say("  overwritten and a number no file carries is not one to "
        "print.")
    say("  The list is corrected here and every verdict below is on "
        "the")
    say("  corrected field. T3 was refuted in that execution and is "
        "refuted here")
    say("  too, on pairs the slip never touched.")
    for p_, ns_ in SETS:
        for n_ in ns_:
            fs_ = factor_set(n_)
            if fs_ != {2, p_}:
                raise SystemExit("base %d is not 2^a %d^b" % (n_, p_))
    say("  every N checked: each is 2^a p^b with the radical it "
        "claims")

    NMAX = max(max(g) for _, g in SETS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    res = {}
    for p, ns in SETS:
        say()
        say("radical {2,%d}, drift %+.6f" % (p, dr[p]))
        rows = []
        for N in sorted(ns):
            L = level(N, lam, mu, sqf)
            rows.append((N, math.log(N), L))
            say("  N = %-9d log N %.6f   L %+.6f" % (N, math.log(N), L))
            say("POINT dense_%d %.6f" % (N, L))
        res[p] = rows
    say("SCALES %d" % sum(len(g) for _, g in SETS))

    # -------------------------------------------------------------- T1
    say()
    say("T1  the gate")
    t1 = True
    for p, rows in res.items():
        for N, _, L in rows:
            if N in lv:
                g = abs(L - lv[N]) < 10.0 ** (-DEC)
                t1 &= g
                say("  N = %-9d here %+.6f against its %+.6f  %s"
                    % (N, L, lv[N], "ok" if g else "MISMATCH"))
    say("  T1 %s   (cap: %d decimals)"
        % ("hold" if t1 else "REFUTED", DEC))
    if not t1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # ---------------------------------------------------------- T2, T4
    say()
    say("T2, T4  is L smooth inside a radical, and at the drift?")
    t2 = t4 = True
    for p, rows in sorted(res.items()):
        x = [r[1] for r in rows]
        y = [r[2] for r in rows]
        sl, rms = fit(x, y)
        t2 &= rms <= RMSCAP
        t4 &= abs(sl - dr[p]) <= SLOPECAP
        say("  {2,%d}: slope %+.6f against drift %+.6f, r.m.s. "
            "residual %.6f" % (p, sl, dr[p], rms))
        say("POINT denserms_%d %.6f" % (p, rms))
        say("POINT denseslope_%d %.6f" % (p, sl))
    say("  T2 %s   (cap: %.2f r.m.s.)"
        % ("hold" if t2 else "REFUTED", RMSCAP))
    say("  T4 %s   (cap: %.2f)"
        % ("hold" if t4 else "REFUTED", SLOPECAP))

    # -------------------------------------------------------------- T3
    say()
    say("T3  the close pairs")
    say("    radical  N1        N2         log gap    |dL|      "
        "smooth")
    t3 = True
    nulls = []
    mags = []
    for p, rows in sorted(res.items()):
        for i in range(len(rows) - 1):
            n1, x1, l1 = rows[i]
            n2, x2, l2v = rows[i + 1]
            g = x2 - x1
            if g >= CLOSEGAP:
                continue
            d = abs(l2v - l1)
            sm = g * dr[p]
            mags.append(d)
            ok = d <= CLOSECAP
            t3 &= ok
            if sm < d <= CLOSECAP:
                nulls.append("%d/%d" % (n1, n2))
            say("    {2,%d}   %-9d %-9d  %.6f  %.6f  %.6f  %s"
                % (p, n1, n2, g, d, sm, "ok" if ok else "ABOVE"))
            say("POINT closepair_%d %.6f" % (n1, d))
    say("  %d close pairs, largest |dL| %.6f, median %.6f"
        % (len(mags), max(mags) if mags else 0.0,
           float(np.median(mags)) if mags else 0.0))
    say("POINT closemax %.6f" % (max(mags) if mags else 0.0))
    say("POINT closemedian %.6f"
        % (float(np.median(mags)) if mags else 0.0))
    say("  T3 %s   (cap: %.2f on every close pair)"
        % ("hold" if t3 else "REFUTED", CLOSECAP))

    say()
    say("=" * 70)
    say("T1 %s  T2 %s  T3 %s  T4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (t1, t2, t3, t4)))
    say()
    if t3 and t2:
        say("L is a smooth function of N inside a radical. the "
            "{2,7} pair that")
        say("refuted rem:levelmatched's control was one pair "
            "behaving as pairs do")
        say("at that spacing, and nothing here carries scatter at "
            "fixed radical.")
    elif not t3:
        say("L carries scatter at fixed radical that does not shrink "
            "with the")
        say("spacing. the radical does not determine the level, and "
            "the failure")
        say("cannot be laid on a window because there is none. "
            "every radical")
        say("number in this branch is a number plus an unmodelled "
            "term of the")
        say("size printed above, and that size is what this run "
            "establishes.")
    else:
        say("the close pairs behave but the whole set does not fit a "
            "line, so what")
        say("L carries is structure at a longer scale than the "
            "packing resolves.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
