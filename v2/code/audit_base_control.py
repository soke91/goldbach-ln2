# -*- coding: utf-8 -*-
r"""
Is the drift a function of the radical, or of the base it was measured on

WHAT IS AT STAKE

rem:primecontrib ended by naming a confound it could not resolve: one
base per prime, so p is not separated from the base's own size or its
2-adic valuation.  Writing this run made the larger version of that
plain.  **Every measurement in this branch has used one base per
radical.**  rem:whichfloor compared 25000 against 30030, rem:radicallaw
compared six bases with six different radicals, rem:primecontrib nine
more -- and in none of them was the same radical measured twice.  So
"the drift depends on which primes divide N" has never been separated
from "the drift depends on the base", and the three remarks that read
a radical dependence read it from a design that could not tell the
difference.

This run is that control and nothing else.  Six bases with radical
{2,3}, spanning 2-adic valuation 2 to 12 and a factor 2.4 in size, and
three with radical {2,5}:

    {2,3}   2^6 3^5 = 15552    2^8 3^4 = 20736    2^2 3^8 = 26244
            2^10 3^3 = 27648   2^4 3^7 = 34992    2^12 3^2 = 36864
    {2,5}   2^7 5^3 = 16000    2^3 5^5 = 25000    2 5^6   = 31250

If the drift is a function of the radical, the six agree and the three
agree.  If it is not, three remarks lose their readings at once.

BACKS: Remark {#rem:basecontrol} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  P1  THE GATE.  Bases 20736 and 25000 reproduce rem:radicallaw's
      drifts +0.212384 and +0.138010 to six decimals -- same bases,
      same code path, one more doubling appended.
  P2  **THE CONTROL.**  The six {2,3} bases agree: their drifts have a
      range under 0.02, the cap this branch has used for agreement
      throughout.
  P3  The three {2,5} bases agree, by the same cap.
  P4  And the two groups stay apart: the {2,3} mean exceeds the {2,5}
      mean by more than 0.05.  rem:radicallaw put the single-base gap
      at +0.074374.

REFUTATION RULE (fixed before the run)

  P1  REFUTED outside six decimals; nothing below is reported.  Note
      that this run appends a tenth doubling to each family, so a
      drift computed here is over a longer window than
      rem:radicallaw's; **the gate is run on the nine N that remark
      used**, and the ten-N drifts are what P2 to P4 judge.
  P2  **REFUTED above a 0.02 range, and this is the most expensive
      outcome this branch can produce.**  It would mean the drift is
      not a function of the radical, and that rem:whichfloor,
      rem:radicallaw and rem:primecontrib each measured a quantity
      that moves with the base -- their numbers would stand and their
      readings would not.  If it fires, the spread is to be reported
      against 2-adic valuation and against base size, so the failure
      names what it is instead, and **no radical statement anywhere
      in this branch may be repeated without that qualification.**
  P3  REFUTED above 0.02.  Three bases is a weaker test than six and
      a failure here with P2 holding would say the {2,5} family in
      particular is base-sensitive.
  P4  REFUTED at or below 0.05.  Then even the largest radical
      difference this branch found does not survive averaging over
      bases.

  **THE UNRESOLVED CASE, NAMED.**  The single-family drift errors in
  rem:radicallaw ran 0.007 to 0.019, which is the size of P2's cap
  itself.  **A range under 0.02 is therefore not evidence that the
  bases agree if the errors are that large**: it is what six draws
  from one distribution with that error would give anyway.  This run
  prints the largest single drift error beside the range, and if the
  range is below it, **P2's verdict word stands and its reading is
  barred** -- the control would then be underpowered rather than
  passed, which is a different thing and must be said.  The same
  applies to P3.  This branch has now written a cap on the wrong
  quantity five times, so it is named before the run.

  WHAT THIS CANNOT DO.  Two radicals, nine bases.  The 2-adic
  valuation and the base size are varied but not independently, and a
  control that passes does not prove radical-only dependence for
  radicals not tested here.  Nothing in this run measures |sum a|
  against any floor or moves item 5's demand.
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
OUT = os.path.join(ROOT, "results", "audit_base_control.txt")
SRC = os.path.join(ROOT, "results", "audit_radical_law.txt")

THETA = 0.56
G23 = [15_552, 20_736, 26_244, 27_648, 34_992, 36_864]
G25 = [16_000, 25_000, 31_250]
NPER = 10
NGATE = 9
DEC = 6
RANGECAP = 0.02
GAPCAP = 0.05


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


def v2(n):
    e = 0
    while n % 2 == 0:
        n //= 2
        e += 1
    return e


def pair(N, lam, mu, sqf):
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


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se


def read_pub():
    src = io.open(SRC, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^POINT raddrift_(\d+) ([-+]?[\d.]+)\s*$",
                         src, re.M):
        out[int(m.group(1))] = float(m.group(2))
    return out


HEAD = [
    "STATISTIC: the drift of log(|sum a|/l2) against log N for nine",
    "           doubling families sharing two radicals -- six on",
    "           {2,3} and three on {2,5} -- to separate a dependence",
    "           on the radical from a dependence on the base.",
    "FIELD: base * 2^j for j < %d; bases %s on {2,3} and %s on {2,5},"
    % (NPER, G23, G25),
    "       spanning 2-adic valuation 2 to 12. k over the squarefree",
    "       k < N^%.2f coprime to N; j over every index below N."
    % THETA,
    "       The two drifts rem:radicallaw published for its shared",
    "       bases are READ from results/audit_radical_law.txt and",
    "       reproduced on the %d N that remark used." % NGATE,
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    for b in (20_736, 25_000):
        say("READ audit_radical_law.txt %d %.6f" % (b, pub[b]))
    say("  the two drifts this run shares with rem:radicallaw")
    say("PRINTBOUND audit_base_control %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, range cap %.2f, gap cap %.2f"
        % (THETA, RANGECAP, GAPCAP))
    say("RADICALS 2")

    NMAX = max(G23 + G25) * (1 << (NPER - 1))
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    res = {}
    for tag, group in (("{2,3}", G23), ("{2,5}", G25)):
        say()
        say("radical %s" % tag)
        rows = []
        for base in group:
            xs, ys = [], []
            for jj in range(NPER):
                N = base * (1 << jj)
                sa, l2 = pair(N, lam, mu, sqf)
                xs.append(math.log(N))
                ys.append(math.log(sa / l2))
            b10, se10 = fit(xs, ys)
            b9, _ = fit(xs[:NGATE], ys[:NGATE])
            rows.append((base, v2(base), b10, se10, b9))
            say("  base %-6d v2 %-3d drift %+.6f +- %.6f   (%d N: "
                "%+.6f)" % (base, v2(base), b10, se10, NGATE, b9))
            say("POINT basedrift_%d %.6f" % (base, b10))
            say("SPREAD basedrift_%d %.6f" % (base, se10))
        res[tag] = rows
    say("SCALES %d" % ((len(G23) + len(G25)) * NPER))

    # -------------------------------------------------------------- P1
    say()
    say("P1  the gate, on the %d N rem:radicallaw used" % NGATE)
    p1 = True
    for tag in res:
        for base, _, _, _, b9 in res[tag]:
            if base not in pub:
                continue
            g = abs(b9 - pub[base]) < 10.0 ** (-DEC)
            p1 &= g
            say("  base %-6d here %+.6f against its %+.6f  %s"
                % (base, b9, pub[base], "ok" if g else "MISMATCH"))
    say("  P1 %s   (cap: %d decimals)"
        % ("hold" if p1 else "REFUTED", DEC))
    if not p1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    worst = max(r[3] for tag in res for r in res[tag])
    say("  largest single drift error %.6f" % worst)

    # ---------------------------------------------------------- P2, P3
    verd = {}
    for name, tag in (("P2", "{2,3}"), ("P3", "{2,5}")):
        say()
        say("%s  do the %s bases agree?" % (name, tag))
        ds = [r[2] for r in res[tag]]
        rng = max(ds) - min(ds)
        ok = rng < RANGECAP
        verd[name] = ok
        say("  drifts %s" % " ".join("%+.6f" % d for d in ds))
        say("  range %.6f over %d bases" % (rng, len(ds)))
        say("POINT baserange_%s %.6f" % (name, rng))
        say("  %s %s   (cap: %.2f)"
            % (name, "hold" if ok else "REFUTED", RANGECAP))
        if rng < worst:
            say("  UNRESOLVED: the range is below the largest single "
                "drift error,")
            say("  so the control is underpowered rather than passed "
                "and the reading")
            say("  is barred, as the rule says")

    # -------------------------------------------------------------- P4
    say()
    say("P4  do the two radicals stay apart?")
    m23 = float(np.mean([r[2] for r in res["{2,3}"]]))
    m25 = float(np.mean([r[2] for r in res["{2,5}"]]))
    gap = m23 - m25
    p4 = gap > GAPCAP
    say("  {2,3} mean %+.6f, {2,5} mean %+.6f, gap %+.6f"
        % (m23, m25, gap))
    say("POINT basegap %.6f" % gap)
    say("  P4 %s   (cap: above %.2f)"
        % ("hold" if p4 else "REFUTED", GAPCAP))

    # what it is instead, if P2 fell
    if not verd["P2"]:
        say()
        say("  P2 fell, so the spread is reported against what it "
            "could be instead")
        rows = res["{2,3}"]
        vv = np.array([float(r[1]) for r in rows])
        ss = np.array([math.log(float(r[0])) for r in rows])
        dd = np.array([r[2] for r in rows])
        for nm, v in (("v2", vv), ("log base", ss)):
            b, se = fit(v, dd)
            say("    drift on %-9s slope %+.6f +- %.6f, t %+.2f"
                % (nm, b, se, b / se))
            say("TSTAT baseon%s %.2f" % (nm.replace(" ", ""), b / se))
            say("SPREAD baseon%s %.6f" % (nm.replace(" ", ""), se))

    # a diagnostic, after the verdicts and predicted by nothing.
    # P2's unresolved clause was written on "the largest single drift
    # error" without saying within which radical, and the largest is
    # in the other group entirely. the test it meant to make is
    # whether each group's own drifts are consistent with one
    # constant, using each family's own error.
    say()
    say("  a diagnostic, after the verdicts and predicted by nothing")
    say("  P2's clause compared the {2,3} range against the largest "
        "error over")
    say("  both radicals, which lies in the other group. the "
        "statistic it meant")
    say("  to make is each group against one constant, on its own "
        "errors.")
    for tag in ("{2,3}", "{2,5}"):
        rr = res[tag]
        vv = np.array([r[2] for r in rr])
        ee = np.array([r[3] for r in rr])
        w = 1.0 / ee ** 2
        wm = float((vv * w).sum() / w.sum())
        chi = float((((vv - wm) / ee) ** 2).sum())
        big = float(ee.max())
        rng = float(vv.max() - vv.min())
        say("    %-6s weighted mean %+.6f, chi-square %.3f on %d "
            "degrees" % (tag, wm, chi, len(rr) - 1))
        say("           range %.6f against this group's largest "
            "error %.6f" % (rng, big))
        say("POINT groupchi_%s %.6f"
            % (tag.replace("{", "").replace("}", "")
               .replace(",", "_"), chi))
    say("  read as a diagnostic and not as the registered control: "
        "the pre-")
    say("  registered P2 is barred by its own clause and a statistic "
        "chosen")
    say("  after the fact does not replace one chosen before it")

    say()
    say("=" * 70)
    say("P1 %s  P2 %s  P3 %s  P4 %s"
        % (("hold" if p1 else "REFUTED"),
           ("hold" if verd["P2"] else "REFUTED"),
           ("hold" if verd["P3"] else "REFUTED"),
           ("hold" if p4 else "REFUTED")))
    say()
    r23 = max(r[2] for r in res["{2,3}"]) - min(r[2] for r in
                                                res["{2,3}"])
    if verd["P2"] and r23 >= worst and p4:
        say("the drift is a function of the radical. six bases "
            "spanning 2-adic")
        say("valuation two to twelve agree inside a cap larger than "
            "their own")
        say("spread, and the two radicals stay apart. the readings of "
            "rem:whichfloor,")
        say("rem:radicallaw and rem:primecontrib survive the control "
            "they never had.")
    elif not verd["P2"]:
        say("the drift is not a function of the radical alone. six "
            "bases sharing")
        say("one radical disagree by more than the cap this branch "
            "has used for")
        say("agreement, so rem:whichfloor, rem:radicallaw and "
            "rem:primecontrib")
        say("each measured a quantity that moves with the base. their "
            "numbers")
        say("stand and their readings do not, and no radical "
            "statement in this")
        say("branch may be repeated without that qualification.")
    else:
        say("the control is underpowered. the bases neither agree nor "
            "disagree by")
        say("more than one family's own error, so this field does not "
            "establish")
        say("that the drift is a function of the radical -- and three "
            "remarks")
        say("continue to rest on a design that cannot tell.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
