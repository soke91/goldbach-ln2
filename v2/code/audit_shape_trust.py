# -*- coding: utf-8 -*-
r"""
How far can a ladder be read before its shapes part?

WHAT IS AT STAKE

OPEN item 1 is that theta' = 0.56 is not located: on twelve rungs two
shapes survive at one standard error and put the crossing 2.78 decades
apart ({#rem:laddershape12}). Stated that way it sounds as though the
ladder says nothing beyond its top rung, which is not true -- the same
two shapes agree on where 1/2 is falls to 0.07 of a decade.

Between those two statements is a number nobody has computed: the N up
to which the surviving shapes are interchangeable, meaning they differ
by less than the ladder's own r.m.s. scatter. Below it the ladder can
be read; above it the forecast is the shape's and not the data's. That
number is what a reader needs, and it is arithmetic once the fits
exist.

The same question applies to the other shape adjudication this
repository runs, the flatness of {#rem:flatnessshape}, whose two
shapes are separated by 0.09 of the r.m.s.'s own standard error and
whose parameters its results file prints.

BACKS: Remark {#rem:shapetrust} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  The control: the twelve-rung fits reproduce
      {#rem:laddershape12}'s r.m.s. 0.00370 and 0.00433 and its
      crossings 11.0762 and 13.8607, to within 0.0001 and 0.001.
  C2  The ladder is readable well past where it was measured: the
      first log10 N at which the two survivors differ by more than
      the twelve-rung r.m.s. is at least two decades above the top
      rung.
  C3  But neither theta' forecast is inside that range: both 11.0762
      and 13.8607 lie above it.
  C4  And the flatness adjudication is the same shape of statement:
      its two shapes also part before the N at which either would
      reach its own bound.

REFUTATION RULE (fixed before the run)

  C1  REFUTED at either cap -- not the same fits, and nothing below
      may be compared with {#rem:laddershape12}.
  C2  REFUTED if the shapes part within two decades of the top rung.
      The ladder would then say nothing at all beyond what it
      measured, and Remark {#rem:primorialrung11}'s reading of the
      crossing of 1/2 would be the only thing it supports.
  C3  REFUTED if either forecast is inside the agreement range, which
      would mean the shape ambiguity does not reach the forecast and
      {#rem:laddershape12}'s 2.78 decades is not the operative
      uncertainty.
  C4  REFUTED if the flatness shapes agree out to the bound, which
      would make that adjudication a different kind of statement from
      this one.

  All four gate.

  NO NULL IS RUN and none applies. Fitted shapes are evaluated and
  compared with a measured scatter; there is no background to detect
  against. The scatter itself is the ladder's, measured in
  audit_primorial_rung11.py with the coin arms of
  lab_primorial_ladder.py behind it.
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
OUT = os.path.join(RES, "audit_shape_trust.txt")

THETA = 0.56
HALF = 0.5
UMAX = 400.0


def read_rungs():
    """the eleven published rungs and the twelfth"""
    src = io.open(os.path.join(RES, "audit_primorial_rung10.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log10 N   exponent   fitted     "
                  "residual")
    ns, ex = [], []
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        ns.append(int(f[0]))
        ex.append(float(f[2]))
    src2 = io.open(os.path.join(RES, "audit_primorial_rung11.txt"),
                   encoding="utf-8").read()
    j = src2.index("N            log10 N   exponent   "
                   "margin over 1/2  floor")
    f = src2[j:].splitlines()[1].split()
    ns.append(int(f[0]))
    ex.append(float(f[2]))
    return ns, ex


def read_shape12():
    """the published twelve-rung r.m.s. and crossings"""
    src = io.open(os.path.join(RES, "audit_ladder_shape12.txt"),
                  encoding="utf-8").read()
    i = src.index("shape                    r.m.s. 11  r.m.s. 12  "
                  "s.e. from best  log10 N at 0.56")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        m = re.match(r"^  (\S.*?)\s{2,}[\d.]+\s+([\d.]+)\s+[\d.]+"
                     r"\s+([\d.]+|none)\s*$", ln)
        if not m:
            break
        out[m.group(1).strip()] = (
            float(m.group(2)),
            None if m.group(3) == "none" else float(m.group(3)))
    return out


def read_flatness():
    """the flatness shapes' parameters and r.m.s."""
    src = io.open(os.path.join(RES, "audit_flatness_shape.txt"),
                  encoding="utf-8").read()
    m = re.search(r"F ~ N\^e\s+([\d.]+)\s+e = \+([\d.]+)", src)
    rp, e = float(m.group(1)), float(m.group(2))
    m2 = re.search(r"F = a \+ b/log N\s+([\d.]+)\s+a = ([\d.]+), "
                   r"b = ([-+][\d.]+)", src)
    rs, a, b = (float(m2.group(1)), float(m2.group(2)),
                float(m2.group(3)))
    m3 = re.search(r"at log10 N = ([\d.]+)", src)
    return rp, e, rs, a, b, float(m3.group(1))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ns, ex = read_rungs()
    pub = read_shape12()
    say("read %d rungs and %d published shape rows"
        % (len(ns), len(pub)))

    u = np.log(np.array(ns, dtype=np.float64))
    y = np.array(ex)

    def design(nm, uu):
        if nm == "line":
            return np.vstack([np.ones_like(uu), uu]).T
        return np.vstack([np.ones_like(uu), np.log(uu)]).T

    fits = {}
    for nm, lab in (("line", "a + b log N"),
                    ("loglog", "a + b log log N")):
        X = design(nm, u)
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ c
        fits[nm] = (c, float(np.sqrt((r ** 2).mean())), lab)

    def val(nm, uu):
        return float(design(nm, np.array([uu])) @ fits[nm][0])

    def crossing(nm, target):
        lo, hi = 2.0, UMAX
        for _ in range(300):
            mid = 0.5 * (lo + hi)
            if val(nm, mid) < target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi) / math.log(10.0)

    # ------------------------------------------------------------- C1
    say()
    say("C1  the control: the two surviving shapes on twelve rungs")
    say("  shape                    r.m.s.    published  0.56 here  "
        "published")
    c1 = True
    for nm in ("line", "loglog"):
        c, rms, lab = fits[nm]
        prms, pcr = pub[lab]
        cr = crossing(nm, THETA)
        if abs(rms - prms) >= 0.0001 or abs(cr - pcr) >= 0.001:
            c1 = False
        say("  %-24s %-9.5f %-10.5f %-10.4f %.4f"
            % (lab, rms, prms, cr, pcr))
    say("  C1 %s   (cap 0.0001 and 0.001)"
        % ("hold" if c1 else "REFUTED"))

    # ---------------------------------------------------------- C2/C3
    say()
    say("C2/C3  where the two surviving shapes part")
    scat = min(fits[nm][1] for nm in ("line", "loglog"))
    top = math.log10(max(ns))
    say("  the ladder's own r.m.s. on twelve rungs is %.5f, and the"
        % scat)
    say("  top rung is at log10 N = %.4f" % top)
    say("  log10 N     line      loglog    gap")
    part = None
    for d in [7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0]:
        uu = d * math.log(10.0)
        g = abs(val("line", uu) - val("loglog", uu))
        if part is None and d > top and g > scat:
            part = d
        say("  %-11.2f %-9.4f %-9.4f %.6f"
            % (d, val("line", uu), val("loglog", uu), g))
    lo, hi = top, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if abs(val("line", mid * math.log(10.0))
               - val("loglog", mid * math.log(10.0))) > scat:
            hi = mid
        else:
            lo = mid
    part = 0.5 * (lo + hi)
    say("  they first differ by more than %.5f at log10 N = %.4f,"
        % (scat, part))
    say("  which is %.4f decades above the top rung" % (part - top))
    c2 = (part - top) >= 2.0
    cr1, cr2 = crossing("line", THETA), crossing("loglog", THETA)
    c3 = (cr1 > part) and (cr2 > part)
    say("  C2 the agreement reaches two decades past the top   %s"
        % ("hold" if c2 else "REFUTED"))
    say("  the two theta' forecasts are at %.4f and %.4f" % (cr1, cr2))
    say("  C3 both lie above the parting point   %s"
        % ("hold" if c3 else "REFUTED"))
    say("TRUST ladder_theta %.4f %.4f" % (part, min(cr1, cr2)))
    if min(cr1, cr2) > part:
        say("FORECAST OUTSIDE ladder_theta")
    say("  both surviving shapes' answers, so the forecast is read")
    say("  as a range and not a point:")
    say("FORECAST BOTH ladder_theta %.4f %.4f %.4f"
        % (min(cr1, cr2), min(cr1, cr2), max(cr1, cr2)))
    say("  and where 1/2 is, which is behind the ladder:")
    say("    line %.4f, loglog %.4f, apart by %.4f"
        % (crossing("line", HALF), crossing("loglog", HALF),
           abs(crossing("line", HALF) - crossing("loglog", HALF))))

    # ------------------------------------------------------------- C4
    say()
    say("C4  the same for the flatness of {#rem:flatnessshape}")
    rp, e, rs, a, b, ucross = read_flatness()
    say("  read from results/audit_flatness_shape.txt: the power's")
    say("  r.m.s. %.6f with e = %.6f, the bounded shape's %.6f with"
        % (rp, e, rs))
    say("  a = %.6f, b = %+.6f, and the power reaching 1 at %.4f"
        % (a, b, ucross))
    f0 = a + b / (math.log(10.0) * 7.4)
    lnF0 = math.log(f0) - e * math.log(10.0) * 7.4
    sc2 = min(rp, rs)
    lo, hi = 7.4, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        uu = mid * math.log(10.0)
        if abs(math.exp(lnF0 + e * uu) - (a + b / uu)) > sc2:
            hi = mid
        else:
            lo = mid
    part2 = 0.5 * (lo + hi)
    c4 = part2 < ucross
    say("  matched at the top of its own range, the two part at")
    say("  log10 N = %.4f, against the bound reached at %.4f"
        % (part2, ucross))
    say("TRUST flatness %.4f %.4f" % (part2, ucross))
    if ucross > part2:
        say("FORECAST OUTSIDE flatness")
    say("  C4 they part before the bound   %s"
        % ("hold" if c4 else "REFUTED"))

    say()
    say("=" * 70)
    ok = c1 and c2 and c3 and c4
    say("the ladder is readable two decades past its top and not to "
        "theta'" if ok else "REFUTED")

    head = [
        "STATISTIC: the two shapes that survive at one standard error",
        "           on the primorial ladder's twelve rungs, their",
        "           r.m.s. residuals and their crossings of",
        "           theta' = " + str(THETA) + " and of 1/2; the log10 N",
        "           at which they first differ by more than that",
        "           r.m.s.; and the same parting point for the two",
        "           shapes of {#rem:flatnessshape}.",
        "NULL: none is run and none applies. Fitted shapes are",
        "      evaluated and compared with a measured scatter; there",
        "      is no background to detect against. The scatter is the",
        "      ladder's, measured in audit_primorial_rung11.py with",
        "      the coin arms of lab_primorial_ladder.py behind it.",
        "FIELD: the eleven rungs read from",
        "       results/audit_primorial_rung10.txt, the twelfth from",
        "       results/audit_primorial_rung11.txt, the published",
        "       shape table from results/audit_ladder_shape12.txt and",
        "       the flatness parameters from",
        "       results/audit_flatness_shape.txt. The ladder is",
        "       N = 30030*2^j, one odd radical throughout.",
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
