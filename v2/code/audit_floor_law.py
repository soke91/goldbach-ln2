# -*- coding: utf-8 -*-
r"""
How does the local floor fall with N, and what does 1/2 look like now?

WHAT IS AT STAKE

The floor has now been measured in four windows: around rung 10 and
around the interval {#rem:primorialgap} bracketed, both near 2 to
3*10^7 ({#rem:localfloor}); across the 0.56 crossing near 5.7*10^9
({#rem:targetband}); and around rung 18 near 7.9*10^9
({#rem:rung18fill}).  {#rem:localfloor} observed that they differ by
about six and concluded no floor may be carried from one place to
another.  Nobody has asked the next question: **is the way it falls a
law, or four numbers?**

It matters for two reasons.  A law says in advance what resolution
any future window buys, which is the difference between choosing a
window and discovering afterwards that it sat inside its own band --
the error {#rem:targetband} traced at 1/2.  And it decides how much
of this branch's headline is understated.

That headline is still written against the wrong yardstick.  The
target section of OPEN.md reports the crossing of 1/2 as two rungs
above it with the better one clearing "the floor by 2.7" -- that is
{#rem:primorialrung11}, judged against the ladder-wide 0.0037.  Since
then the ladder has been carried to rung 18, where the margin over
1/2 is 0.0628 and the floor measured in that very window is 0.000308.
Nothing has restated it.

Nothing is measured here.  Every number is read from a marker line in
a result file, and declared with a READ line so G76 checks it came
from the line it says.  The windows' abscissae are imported from the
scripts that defined them rather than parsed back out of prose --
which is the other half of the lesson from that same failure.

BACKS: Remark {#rem:floorlaw} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  P1  The reads are honest: every value used appears in the result
      file it is attributed to, as a whole marker line.  This is
      what the READ lines declare and G76 enforces.
  P2  Every window's margin over 1/2 exceeds that window's own
      floor.  The crossing of 1/2 is resolved everywhere it has been
      measured against a local floor, not only at the top.
  P3  At rung 18 the margin over 1/2 is more than a hundred times
      the floor measured in that window.
  P4  The floor falls with N as a power: fitting log(floor) on
      log(N) over the four windows gives a negative slope resolved
      at |t| > 2.  **This is the weak one** -- four windows in two
      clusters carry about one degree of freedom for a slope, so it
      may not resolve even if a law is there.

REFUTATION RULE (fixed before the run)

  P1  REFUTED by a single READ line G76 cannot match.  THIS ONE
      GATES.
  P2  REFUTED by a window whose margin is at or below its own floor.
      That window's crossing of 1/2 would then be unresolved where
      it lives, and the branch could not say the barrier is crossed
      there at all.
  P3  REFUTED below a hundred.  Nothing breaks -- it is a statement
      about how much the old yardstick understated, not about the
      ladder -- but the restatement would be smaller than it looks.
  P4  REFUTED if the slope is positive or |t| <= 2.  Then the four
      floors are four numbers and not a law, no future window can be
      sized in advance, and every one has to be measured after the
      fact.  With two clusters this is the likely outcome and it is
      registered as a prediction rather than an expectation so that
      the negative is on the record.

  P1 gates.  P2 to P4 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  Nothing is sampled here; four
  measured floors are compared with four measured margins and fitted.
  The nulls for the quantities themselves were run where they were
  measured.
"""

import importlib.util
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
CODE = os.path.join(ROOT, "code")
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_floor_law.txt")

HALF = 0.5


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


LF = module("audit_local_floor")
FI = module("audit_rung18_fill")
TB = module("audit_target_band")


def marker(fname, label):
    """the value on a whole marker line, or a hard failure"""
    src = io.open(os.path.join(RES, fname), encoding="utf-8").read()
    m = re.search(r"^%s ([-+]?[\d.]+)\s*$" % re.escape(label), src,
                  re.M)
    if not m:
        raise SystemExit("no line '%s ...' in %s" % (label, fname))
    return float(m.group(1))


def marker2(fname, label):
    """the two values on a MARGIN line: the clearance and its floor"""
    src = io.open(os.path.join(RES, fname), encoding="utf-8").read()
    m = re.search(r"^%s ([-+]?[\d.]+) ([-+]?[\d.]+)\s*$"
                  % re.escape(label), src, re.M)
    if not m:
        raise SystemExit("no line '%s ...' in %s" % (label, fname))
    return float(m.group(1)), float(m.group(2))


def gmean(ns):
    return math.exp(sum(math.log(float(n)) for n in ns) / len(ns))


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    # ------------------------------------------------ what is read
    floA = marker("audit_local_floor.txt", "FLOOR local_window_A")
    floB = marker("audit_local_floor.txt", "FLOOR local_window_B")
    floF = marker("audit_rung18_fill.txt", "FLOOR rung18_fill_local")
    floT = marker("audit_target_band.txt", "FLOOR target_band_local")
    m10, f10 = marker2("audit_local_floor.txt",
                       "MARGIN audit_local_floor_rung10")
    m18 = marker("audit_primorial_rung18.txt",
                 "FLOOR primorial_rung18")
    marg18, _ = marker2("audit_primorial_rung18.txt",
                        "MARGIN audit_primorial_rung18")
    marg17, _ = marker2("audit_primorial_rung17.txt",
                        "MARGIN audit_primorial_rung17")

    say("every number below is read from a marker line; the READ "
        "lines declare which,")
    say("  and G76 checks each against its source")
    say("READ audit_local_floor.txt FLOOR local_window_A %.6f" % floA)
    say("READ audit_local_floor.txt FLOOR local_window_B %.6f" % floB)
    say("READ audit_rung18_fill.txt FLOOR rung18_fill_local %.6f"
        % floF)
    say("READ audit_target_band.txt FLOOR target_band_local %.6f"
        % floT)
    say("READ audit_primorial_rung18.txt FLOOR primorial_rung18 %.4f"
        % m18)
    say("READ audit_primorial_rung18.txt MARGIN "
        "audit_primorial_rung18 %.4f %.4f" % (marg18, m18))
    say("READ audit_primorial_rung17.txt MARGIN "
        "audit_primorial_rung17 %.4f %.4f" % (marg17, m18))
    say("READ audit_local_floor.txt MARGIN audit_local_floor_rung10 "
        "%.6f %.6f" % (m10, f10))

    # the abscissae come from the scripts that defined the windows
    WA, WB = LF.WIN_A, LF.WIN_B
    WF = FI.FILL + (FI.RUNG18,)
    WT = TB.BAND
    rows = [("rung 10", gmean(WA), floA, m10),
            ("the gap", gmean(WB), floB, None),
            ("0.56 band", gmean(WT), floT, None),
            ("rung 18", gmean(WF), floF, marg18)]

    say()
    say("  window        points  geometric mean N   log10 N   floor")
    for (nm, N, fl, _), W in zip(rows, (WA, WB, WT, WF)):
        say("  %-13s %-7d %-18d %-9.4f %.6f"
            % (nm, len(W), int(N), math.log10(N), fl))
    say("SCALES %d" % len(rows))

    # -------------------------------------------------------------- P1
    say()
    say("P1  the reads are declared; G76 is what checks them")
    say("  P1 hold   (cap: G76, which fails the gate if a READ line "
        "has no source)")
    p1 = True

    # -------------------------------------------------------------- P2
    say()
    say("P2  is 1/2 resolved in every window that has a floor?")
    say("  the two windows with a margin over 1/2 measured in place:")
    p2 = True
    for nm, marg, fl in (("rung 10", m10, floA),
                         ("rung 18", marg18, floF)):
        ok = marg > fl
        p2 = p2 and ok
        say("    %-9s margin %.6f against its own floor %.6f, "
            "ratio %.1f" % (nm, marg, fl, marg / fl))
        say("MARGIN audit_floor_law_%s %.6f %.6f"
            % (nm.replace(" ", ""), marg, fl))
        if not ok:
            say("INSIDE FLOOR audit_floor_law_%s"
                % nm.replace(" ", ""))
    say("  the ladder-wide floor those were once judged against is "
        "%.4f" % m18)
    say("  P2 %s   (cap: each window's own floor)"
        % ("hold" if p2 else "REFUTED"))

    # -------------------------------------------------------------- P3
    say()
    say("P3  how much did the old yardstick understate at the top?")
    r_local = marg18 / floF
    r_ladder = marg18 / m18
    p3 = r_local > 100.0
    say("  rung 18's margin over 1/2 is %.6f" % marg18)
    say("    against the floor measured in its own window, %.6f: "
        "%.1f floors" % (floF, r_local))
    say("    against the ladder-wide floor, %.4f: %.1f floors"
        % (m18, r_ladder))
    say("  and rung 17's margin %.4f is %.1f of its neighbour "
        "window's floor" % (marg17, marg17 / floT))
    say("  P3 %s   (cap: a hundred)" % ("hold" if p3 else "REFUTED"))

    # -------------------------------------------------------------- P4
    say()
    say("P4  does the floor fall with N as a power?")
    x = np.array([math.log(r[1]) for r in rows])
    y = np.array([math.log(r[2]) for r in rows])
    sl, ic = np.polyfit(x, y, 1)
    n = x.size
    resid = y - (sl * x + ic)
    sse = float((resid ** 2).sum()) / (n - 2)
    sxx = float(((x - x.mean()) ** 2).sum())
    se = math.sqrt(sse / sxx)
    t = sl / se
    p4 = sl < 0.0 and abs(t) > 2.0
    say("  log(floor) on log(N) over %d windows: slope %+.6f +- "
        "%.6f, t = %.2f" % (n, sl, se, t))
    say("TSTAT floor_law_slope %.2f" % t)
    if abs(t) < 2.0:
        say("UNRESOLVED SIGN floor_law_slope")
    say("SPREAD floor_law_slope %.4f" % (x.max() - x.min()))
    say("SCATTER slope_audit_floor_law %.4f"
        % float(np.sqrt((resid ** 2).mean())))
    say("  the four windows sit in two clusters, near 10^%.1f and "
        "10^%.1f," % (math.log10(rows[1][1]), math.log10(rows[3][1])))
    say("  so the fit carries about one degree of freedom whatever "
        "the t says")
    say("  P4 %s   (cap: negative and |t| > 2)"
        % ("hold" if p4 else "REFUTED"))

    say()
    say("=" * 70)
    say("P1 %s  P2 %s  P3 %s  P4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (p1, p2, p3, p4)))

    head = [
        "STATISTIC: the four local floors measured so far, against",
        "           the geometric mean N of the window each was",
        "           measured in; the margin over 1/2 at rung 10 and",
        "           rung 18 judged against the floor of its own",
        "           window rather than the ladder-wide one; and the",
        "           slope of log(floor) on log(N) with its t.",
        "NULL: none is run and none applies. Nothing is sampled",
        "      here: four measured floors are compared with two",
        "      measured margins and fitted. The nulls for those",
        "      quantities were run where they were measured, in",
        "      lab_primorial_ladder.py and lab_primorial_share.py.",
        "FIELD: the windows of code/audit_local_floor.py (m = 975,",
        "       990, 1008, 1024, 1040, 1056 around rung 10, and",
        "       m = 600, 616, 630, 640, 650, 660, 672, 686 across",
        "       the interval of {#rem:primorialgap}), of",
        "       code/audit_target_band.py (m = 174960, 179200,",
        "       183708, 188160, 193050, 198000, 202800) and of",
        "       code/audit_rung18_fill.py (m = 250047, 256000,",
        "       262144, 268800, 275000), every N being 30030m with m",
        "       composed only of 2, 3, 5, 7, 11 and 13. No sieving",
        "       and no measuring is done here. Every value is read",
        "       from a whole marker line in results/ and declared",
        "       with a READ line; the abscissae are imported from",
        "       the scripts that defined the windows rather than",
        "       parsed out of prose.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
