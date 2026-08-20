# -*- coding: utf-8 -*-
r"""
The quantity itself at matched N, with no slope anywhere in it

WHAT IS AT STAKE

rem:valuation found that the statistic this branch has compared for
six ticks is not stable: refitting one family on its own first six and
last six points moves its drift by 0.03 to 0.08 -- base 25600 gives
+0.184332 and +0.114863 -- which is larger than the {2,3} against
{2,5} gap of 0.074374 the branch was chasing.  None of that is in the
printed standard errors, which measure scatter about a line and not
the movement of the line under its own window.  And since a different
base gives a different set of N, every cross-family comparison here
sat on a different window; ordered by how far their bases span in
size, the runs' disagreements go 2.837x -> chi-square 14.232684,
2.000x -> 0.035607, 1.600x -> 0.008540, 1.500x -> chi-square 2.688194.
**The disagreement tracked the size span, not the arithmetic.**

The fix is not a better slope.  It is to stop comparing slopes.
L(N) = log(|sum a| / l2) is measured at a single N with no fitting at
all, so comparing L at nearly equal N removes the window entirely.

Twelve N inside a band of factor 1.20, ten radicals, and **three
radicals appear twice inside the band** so the same-radical case is a
control rather than an assumption:

    {2}          2097152
    {2,3}        1889568   2125764
    {2,5}        2000000
    {2,7}        1882384   2151296
    {2,11}       1874048
    {2,13}       1827904
    {2,19}       2085136
    {2,3,5}      1800000   2025000
    {2,3,7}      2000376
    {2,3,5,7,11,13}  1921920

The band spans 0.182 in log N.  L rises with log N at something near
the drift this branch has measured, so N position can contribute at
most about 0.026 to any difference inside the band -- **no trend
correction is applied and none is needed**, and that bound is printed
rather than assumed.

BACKS: Remark {#rem:levelmatched} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  THE GATE.  |sum a| at N = 200000 reproduces the published
      87895.3236 to four decimals.
  S2  **The control.**  Each of the three radicals appearing twice
      gives two L within 0.05 of each other.
  S3  **The radicals separate.**  The spread of L across the twelve N
      exceeds 0.3, an order above the 0.026 the band's width can
      contribute.
  S4  And the ordering is the one six ticks of drift measurement
      implied: the primorial N has the largest L of the twelve.

REFUTATION RULE (fixed before the run)

  S1  REFUTED outside four decimals; nothing below is reported.
  S2  **REFUTED above 0.05 on any pair, and that is the outcome that
      costs the most.**  Two N of one radical, sixteen per cent apart
      in size, disagreeing in a quantity measured with no fitting
      would mean the radical does not determine the level either --
      and unlike every previous failure in this branch it could not
      be blamed on a window, because there is no window.  It would
      say the quantity depends on N in a way no radical statement can
      capture, and every remark from rem:whichfloor onward would be
      describing something that is not there.
  S3  REFUTED at or below 0.3.  Then the radicals do not separate in
      the level even though they separated in the drift, and the two
      measurements disagree about their own subject.
  S4  REFUTED if any N has a larger L than the primorial one.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  The band contributes
  up to about 0.026 to a pairwise difference through N position alone.
  **A difference below that number is not attributable to the radical
  whatever the verdict says**, and this run prints the actual log-N
  gap of every pair it judges so the bound is checkable rather than
  quoted.  S2's cap of 0.05 is deliberately set above it: a pair
  differing by less than 0.026 tests nothing, and if the control pairs
  come in under that, S2 holds for a reason that is not evidence and
  the remark must say so.

  WHAT THIS CANNOT DO.  Twelve N at one size.  This measures the
  level, not its trend, so nothing here speaks to how the level moves
  with N or to item 5's demand, which is a statement about exponents.
  Ten radicals do not determine a function of the radical, and no fit
  in the radical is attempted.
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
OUT = os.path.join(ROOT, "results", "audit_level_matched.txt")
SRC = os.path.join(ROOT, "results", "audit_deficit_direct.txt")

THETA = 0.56
NGATE = 200_000
BAND = [2_097_152, 1_889_568, 2_125_764, 2_000_000, 1_882_384,
        2_151_296, 1_874_048, 1_827_904, 2_085_136, 1_800_000,
        2_025_000, 2_000_376, 1_921_920]
DEC = 4
PAIRCAP = 0.05
SPREADCAP = 0.3
TRENDGUESS = 0.14


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
    """|sum a| and l2 at one N -- no fitting anywhere"""
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
    return sa, math.sqrt(l2sq)


def drift_range():
    """the drifts this branch has measured, from their own markers"""
    vals = []
    for fn in sorted(os.listdir(os.path.join(ROOT, "results"))):
        if not fn.endswith(".txt"):
            continue
        src = io.open(os.path.join(ROOT, "results", fn),
                      encoding="utf-8").read()
        for m in re.finditer(r"^POINT (?:raddrift|basedrift|"
                             r"blinddrift|valdrift)_\d+ "
                             r"([-+]?[\d.]+)\s*$", src, re.M):
            vals.append(float(m.group(1)))
    if not vals:
        raise SystemExit("no drift markers found")
    return min(vals), max(vals)


def read_pub():
    m = re.search(r"^POINT deficitdirect_%d ([\d.eE+-]+) " % NGATE,
                  io.open(SRC, encoding="utf-8").read(), re.M)
    if not m:
        raise SystemExit("no POINT marker for N = %d" % NGATE)
    return float(m.group(1))


HEAD = [
    "STATISTIC: L(N) = log(|sum a| / l2) measured at single N, with no",
    "           fitting of any kind, for twelve N inside a band of",
    "           factor 1.20 covering ten radicals, three of which",
    "           appear twice as the control.",
    "FIELD: N in %s; k over the squarefree k < N^%.2f coprime to N;"
    % ([min(BAND), max(BAND)], THETA),
    "       j over every index below N. |sum a| at N = %d is READ"
    % NGATE,
    "       from results/audit_deficit_direct.txt as the gate.",
    "NOTE: no slope is fitted here. rem:valuation measured that a",
    "      ten-point slope of this quantity moves 0.03 to 0.08 under",
    "      its own window; a level at one N has no window.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    say("READ audit_deficit_direct.txt %d %.4f" % (NGATE, pub))
    say("  |sum a| at the gate N")
    say("PRINTBOUND audit_level_matched %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, pair cap %.2f, spread cap %.2f"
        % (THETA, PAIRCAP, SPREADCAP))
    lo, hi = math.log(min(BAND)), math.log(max(BAND))
    say("  trend guess %.2f per log unit; band %.6f in log N, so N"
        % (TRENDGUESS, hi - lo))
    say("  position contributes at most %.6f on that guess"
        % (TRENDGUESS * (hi - lo)))
    say("POINT bandbound %.6f" % (TRENDGUESS * (hi - lo)))
    say("RADICALS %d"
        % len(set(tuple(sorted(factor_set(n))) for n in BAND)))

    NMAX = max(BAND + [NGATE])
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    # -------------------------------------------------------------- S1
    sag, _ = level(NGATE, lam, mu, sqf)
    say()
    say("S1  the gate")
    s1 = abs(round(sag, DEC) - round(pub, DEC)) < 10.0 ** (-DEC)
    say("  |sum a| at %d here %.4f against its %.4f  %s"
        % (NGATE, sag, pub, "ok" if s1 else "MISMATCH"))
    say("  S1 %s   (cap: %d decimals)"
        % ("hold" if s1 else "REFUTED", DEC))
    if not s1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rows = []
    for N in sorted(BAND):
        sa, l2 = level(N, lam, mu, sqf)
        rad = tuple(sorted(factor_set(N)))
        rows.append((N, rad, math.log(sa / l2)))
    say()
    say("    N          radical                  L = log(|sum a|/l2)")
    for N, rad, L in sorted(rows, key=lambda r: r[2]):
        say("  %-10d %-24s %+.6f" % (N, str(rad), L))
        say("POINT levelmatched_%d %.6f" % (N, L))
    say("SCALES %d" % len(rows))

    # -------------------------------------------------------------- S2
    say()
    say("S2  the control: do two N of one radical agree?")
    byrad = {}
    for N, rad, L in rows:
        byrad.setdefault(rad, []).append((N, L))
    s2 = True
    weak = []
    for rad, lst in sorted(byrad.items()):
        if len(lst) < 2:
            continue
        (n1, l1), (n2, l2v) = lst[0], lst[1]
        d = abs(l1 - l2v)
        g = abs(math.log(n1) - math.log(n2))
        s2 &= d <= PAIRCAP
        if d < TRENDGUESS * (hi - lo):
            weak.append(str(rad))
        say("  %-24s %d and %d: |dL| %.6f, log N gap %.6f"
            % (str(rad), n1, n2, d, g))
        say("POINT pairgap_%d %.6f" % (n1, d))
    say("  S2 %s   (cap: %.2f)"
        % ("hold" if s2 else "REFUTED", PAIRCAP))
    lodr, hidr = drift_range()
    say("  NOTE, disclosed: the band bound above uses one trend, "
        "%.2f. This" % TRENDGUESS)
    say("  branch's own drift measurements run %.6f to %.6f -- a "
        "factor of four" % (lodr, hidr))
    say("  -- so the bound is not one number either, and a pair's "
        "own bound is its")
    say("  log N gap times its own radical's drift. The rule is not "
        "rewritten; the")
    say("  registration used a single trend where the branch had "
        "measured a range,")
    say("  and that is recorded as the defect it is.")
    say("READ drift range %.6f %.6f" % (lodr, hidr))
    for rad, lst in sorted(byrad.items()):
        if len(lst) < 2:
            continue
        (n1, l1), (n2, l2v) = lst[0], lst[1]
        g = abs(math.log(n1) - math.log(n2))
        say("    %-24s |dL| %.6f against %.6f .. %.6f"
            % (str(rad), abs(l1 - l2v), g * lodr, g * hidr))
    if weak:
        say("  NOTE: %s came in under the band bound %.6f, so those "
            "pairs test" % (", ".join(weak), TRENDGUESS * (hi - lo)))
        say("  nothing and S2 holds for them by a margin that is not "
            "evidence,")
        say("  as the rule says")

    # -------------------------------------------------------------- S3
    say()
    say("S3  do the radicals separate in the level?")
    Ls = [r[2] for r in rows]
    spread = max(Ls) - min(Ls)
    s3 = spread > SPREADCAP
    say("  L runs %+.6f to %+.6f, spread %.6f"
        % (min(Ls), max(Ls), spread))
    say("POINT levelspread %.6f" % spread)
    say("  against the band's own contribution %.6f"
        % (TRENDGUESS * (hi - lo)))
    say("  S3 %s   (cap: above %.2f)"
        % ("hold" if s3 else "REFUTED", SPREADCAP))

    # -------------------------------------------------------------- S4
    say()
    say("S4  is the primorial N the largest?")
    prim = [r for r in rows if len(r[1]) == 6][0]
    top = max(rows, key=lambda r: r[2])
    s4 = top[0] == prim[0]
    say("  primorial %d has L %+.6f; the largest is %d at %+.6f"
        % (prim[0], prim[2], top[0], top[2]))
    say("  S4 %s" % ("hold" if s4 else "REFUTED"))

    say()
    say("=" * 70)
    say("S1 %s  S2 %s  S3 %s  S4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (s1, s2, s3, s4)))
    say()
    if s2 and s3:
        say("measured with no fitting anywhere, the level separates "
            "the radicals by")
        say("far more than the band can contribute, and two N of one "
            "radical agree.")
        say("that is the first statement in this branch about the "
            "radical that no")
        say("window can be blamed for, and it is about the level and "
            "not the drift.")
    elif not s2:
        say("two N of one radical, measured with no fitting, "
            "disagree. the radical")
        say("does not determine the level either, and this failure "
            "cannot be laid")
        say("on a window because there is no window. everything this "
            "branch has")
        say("said about radicals is describing something that is not "
            "there.")
    else:
        say("the radicals do not separate in the level though they "
            "separated in")
        say("the drift, so the two measurements disagree about their "
            "own subject")
        say("and neither is usable until that is resolved.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
