# -*- coding: utf-8 -*-
r"""
A twelfth rung: does it narrow the shape, or only move the fit?

WHAT IS AT STAKE

Remark {#rem:laddershape} is the reason theta' = 0.56 is undetermined.
Five shapes are fitted to the primorial ladder's eleven rungs; the
line is the best, but with nine degrees of freedom the r.m.s. carries
a 23.6 per cent standard error of its own, so three shapes survive at
one standard error and put the crossing of theta' at 11.2700, 14.6167
and 19.6207 -- 8.35 decades apart.

Remark {#rem:primorialrung11} has since measured a twelfth rung,
N = 61501440 at exponent 0.5099, and nothing has asked what it does to
that adjudication. It is the cheapest question left on the target: no
new arithmetic, one more point in a fit that already exists.

The shapes differ most where N is largest, so one more octave is
exactly where discrimination lives -- and it is also only one octave,
so the honest prior is that it moves the fits and not the ranking.

The implementation is independent of audit_ladder_shape.py's.

BACKS: Remark {#rem:laddershape12} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  J1  The control: on the eleven rungs the five r.m.s. reproduce the
      published 0.00372, 0.00412, 0.00432, 0.00473, 0.03764 to within
      0.00001 and the theta' crossings 11.2700, 14.6167, 19.6207,
      82.5771, 7.6204 to within 0.01 decades.
  J2  The twelfth rung does not change the ranking: the line is still
      the best of the five.
  J3  But it does narrow the field: at one standard error of the best
      r.m.s., now on ten degrees of freedom, fewer than three shapes
      survive.
  J4  And the extrapolation stays undetermined all the same: the
      surviving shapes' theta' crossings still span more than one
      decade.

REFUTATION RULE (fixed before the run)

  J1  REFUTED at either cap -- not the same statistic, and nothing
      below may be compared with {#rem:laddershape}.
  J2  REFUTED if any other shape fits better on twelve. That would be
      a real change: the line is what every published forecast on
      this ladder uses.
  J3  REFUTED if three or more still survive. One more octave would
      then have bought no discrimination at all, and the shape
      question would be out of reach of adding rungs one at a time --
      which is the practical question, since each rung costs twice
      the last.
  J4  REFUTED if the surviving crossings span a decade or less. The
      target theta' = 0.56 would then be located to within a factor
      ten, which no result in these papers currently claims.

  All four gate.

  NO NULL IS RUN and none applies. Five deterministic shapes are
  fitted to the same twelve measurements and compared by r.m.s.; there
  is no background to detect against. The measurements' own noise
  floor is the ladder scatter that audit_primorial_rung11.py declares,
  and the coin arms for the ladder were run in lab_primorial_ladder.py.
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
OUT = os.path.join(RES, "audit_ladder_shape12.txt")

THETA = 0.56
HALF = 0.5
UMAX = 400.0

NAMES = ["line", "saturating", "loglog", "heuristic2", "heuristic1"]
LABEL = {"line": "a + b log N",
         "saturating": "a + b / log N",
         "loglog": "a + b log log N",
         "heuristic2": "a + b log log N / log N",
         "heuristic1": "1 - c log log N / log N"}


def design(name, u):
    if name == "line":
        return np.vstack([np.ones_like(u), u]).T
    if name == "saturating":
        return np.vstack([np.ones_like(u), 1.0 / u]).T
    if name == "loglog":
        return np.vstack([np.ones_like(u), np.log(u)]).T
    if name == "heuristic2":
        return np.vstack([np.ones_like(u), np.log(u) / u]).T
    raise KeyError(name)


def fit_all(u, y):
    """the five shapes, each returning (coefficients, r.m.s.)"""
    out = {}
    for nm in NAMES[:-1]:
        X = design(nm, u)
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ c
        out[nm] = (c, float(np.sqrt((r ** 2).mean())))
    z = np.log(u) / u
    c1 = float(((1.0 - y) * z).sum() / (z * z).sum())
    r1 = y - (1.0 - c1 * z)
    out["heuristic1"] = (np.array([c1]),
                         float(np.sqrt((r1 ** 2).mean())))
    return out


def value(fits, nm, uu):
    c = fits[nm][0]
    if nm == "heuristic1":
        return 1.0 - c[0] * math.log(uu) / uu
    return float(design(nm, np.array([uu])) @ c)


def crossing(fits, nm, target):
    """where the fitted shape reaches the target, in log10 N"""
    if value(fits, nm, UMAX) < target:
        return None
    lo, hi = 2.0, UMAX
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if value(fits, nm, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi) / math.log(10.0)


def read_eleven():
    """the eleven published rungs"""
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
    return ns, ex


def read_twelfth():
    """the rung audit_primorial_rung11.py measured"""
    src = io.open(os.path.join(RES, "audit_primorial_rung11.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^  (\d+)\s+[\d.]+\s+([\d.]+)\s+[\d.]+\s+"
                  r"[\d.]+\s*$", src, re.M)
    i = src.index("N            log10 N   exponent   "
                  "margin over 1/2  floor")
    f = src[i:].splitlines()[1].split()
    return int(f[0]), float(f[2])


def read_published_shapes():
    """the five r.m.s. and theta' crossings on eleven rungs"""
    src = io.open(os.path.join(RES, "audit_ladder_shape.txt"),
                  encoding="utf-8").read()
    i = src.index("shape")
    rows = re.findall(r"^  (\S[^|]*?)\s{2,}([\d.]+)\s+([\d.]+)\s*$",
                      src[i:], re.M)
    return rows


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ns, ex = read_eleven()
    n12, e12 = read_twelfth()
    say("read %d rungs from results/audit_primorial_rung10.txt and "
        "the twelfth," % len(ns))
    say("  N = %d at exponent %.4f, from "
        "results/audit_primorial_rung11.txt" % (n12, e12))

    u11 = np.array([math.log(n) for n in ns])
    y11 = np.array(ex)
    u12 = np.append(u11, math.log(n12))
    y12 = np.append(y11, e12)

    f11 = fit_all(u11, y11)
    f12 = fit_all(u12, y12)

    # ------------------------------------------------------------- J1
    say()
    say("J1  the control: the five shapes on the eleven rungs")
    pubsrc = io.open(os.path.join(RES, "audit_ladder_shape.txt"),
                     encoding="utf-8").read()
    # the two published tables are read separately: H2 carries the
    # r.m.s. and H3/H4 the crossings, and the rows look alike
    h2 = pubsrc[pubsrc.index("shape                        r.m.s. "
                             "residual   vs the line"):]
    h3 = pubsrc[pubsrc.index("shape                        log10 N "
                             "at 0.50   at 0.56"):]
    say("  shape                    r.m.s.    published  log10 N at "
        "%.2f  published" % THETA)
    j1 = True
    for nm in NAMES:
        lab = LABEL[nm]
        pr = float(re.search(r"^\s*" + re.escape(lab) +
                             r"\s+([\d.]+)\s", h2, re.M).group(1))
        pc = float(re.search(r"^\s*" + re.escape(lab) +
                             r"\s+[\d.]+\s+([\d.]+)\s*$",
                             h3, re.M).group(1))
        rr = f11[nm][1]
        cc = crossing(f11, nm, THETA)
        d1 = abs(rr - pr)
        d2 = abs(cc - pc) if cc is not None else float("inf")
        if d1 >= 0.00001 or d2 >= 0.01:
            j1 = False
        say("  %-24s %-9.5f %-10.5f %-15.4f %.4f"
            % (lab, rr, pr, cc, pc))
    say("  J1 %s   (cap 0.00001 in r.m.s., cap 0.01 in decades)"
        % ("hold" if j1 else "REFUTED"))

    # ------------------------------------------------------------- J2
    say()
    say("J2/J3  the same five on twelve rungs")
    best11 = min(f11[nm][1] for nm in NAMES)
    best12 = min(f12[nm][1] for nm in NAMES)
    se11 = best11 / math.sqrt(2.0 * (len(ns) - 2))
    se12 = best12 / math.sqrt(2.0 * (len(ns) + 1 - 2))
    say("  the best r.m.s. carries a standard error of its own:")
    say("    eleven rungs  %.5f / sqrt(%d) = %.5f  (%.1f per cent)"
        % (best11, 2 * (len(ns) - 2), se11, 100 * se11 / best11))
    say("    twelve rungs  %.5f / sqrt(%d) = %.5f  (%.1f per cent)"
        % (best12, 2 * (len(ns) - 1), se12, 100 * se12 / best12))
    say("  shape                    r.m.s. 11  r.m.s. 12  "
        "s.e. from best  log10 N at %.2f" % THETA)
    surv12, cross12 = [], {}
    for nm in NAMES:
        rr = f12[nm][1]
        z = (rr - best12) / se12
        cc = crossing(f12, nm, THETA)
        cross12[nm] = cc
        if z <= 1.0:
            surv12.append(nm)
        say("  %-24s %-10.5f %-10.5f %-15.2f %s"
            % (LABEL[nm], f11[nm][1], rr, z,
               "none" if cc is None else "%.4f" % cc))
    bestnm12 = min(NAMES, key=lambda nm: f12[nm][1])
    j2 = bestnm12 == "line"
    say("  best of the five on twelve: %s" % LABEL[bestnm12])
    say("  J2 the line is still the best   %s"
        % ("hold" if j2 else "REFUTED"))
    surv11 = [nm for nm in NAMES
              if (f11[nm][1] - best11) / se11 <= 1.0]
    j3 = len(surv12) < 3
    say("  surviving at one standard error: %d on eleven, %d on twelve"
        % (len(surv11), len(surv12)))
    say("SHAPESURVIVE ladder_theta %d %d %.4f"
        % (len(ns), len(surv11),
           max(crossing(f11, nm, THETA) for nm in surv11)
           - min(crossing(f11, nm, THETA) for nm in surv11)))
    say("  J3 fewer than three survive on twelve   %s"
        % ("hold" if j3 else "REFUTED"))

    # ------------------------------------------------------------- J4
    say()
    say("J4  and what the survivors say about theta' = %.2f" % THETA)
    cs = [cross12[nm] for nm in surv12 if cross12[nm] is not None]
    spread = max(cs) - min(cs)
    j4 = spread > 1.0
    for nm in surv12:
        say("  %-24s %s"
            % (LABEL[nm], "none" if cross12[nm] is None
               else "%.4f" % cross12[nm]))
    say("  spread %.4f decades   (floor 1 decade)" % spread)
    two = sorted(f12[nm][1] for nm in NAMES)[:2]
    say("SHAPEGAP ladder_theta %.6f %.6f" % (two[1] - two[0], se12))
    if two[1] - two[0] <= se12:
        say("SHAPES TIED ladder_theta")
    say("SHAPESURVIVE ladder_theta %d %d %.4f"
        % (len(ns) + 1, len(surv12), spread))
    say("SHAPECURRENT ladder_theta %d" % (len(ns) + 1))
    say("  J4 %s" % ("hold" if j4 else "REFUTED"))

    say()
    say("  and where the same twelve put 1/2, which is behind the")
    say("  ladder rather than ahead of it and so is not an")
    say("  extrapolation:")
    for nm in surv12:
        c = crossing(f12, nm, HALF)
        say("  %-24s %s"
            % (LABEL[nm], "none" if c is None else "%.4f" % c))

    say()
    say("=" * 70)
    ok = j1 and j2 and j3 and j4
    say("one more rung narrows the field and not the forecast"
        if ok else "REFUTED")

    head = [
        "STATISTIC: five shapes -- a+b log N, a+b/log N,",
        "           a+b log log N, a+b log log N/log N and",
        "           1-c log log N/log N -- least-squares fitted to the",
        "           primorial ladder's exponent against log N, on the",
        "           eleven published rungs and on the twelve that",
        "           include audit_primorial_rung11.py's; each shape's",
        "           r.m.s. residual, its distance from the best in",
        "           standard errors of that r.m.s., and the log10 N at",
        "           which it reaches theta' = " + str(THETA) + ".",
        "NULL: none is run and none applies. Deterministic shapes are",
        "      fitted to the same measurements and compared by r.m.s.;",
        "      there is no background to detect against. The",
        "      measurements' noise floor is the ladder scatter",
        "      audit_primorial_rung11.py declares and the coin arms",
        "      were run in lab_primorial_ladder.py.",
        "FIELD: the eleven rungs read from",
        "       results/audit_primorial_rung10.txt and the twelfth",
        "       from results/audit_primorial_rung11.txt; the published",
        "       shape table from results/audit_ladder_shape.txt. The",
        "       ladder is N = 30030*2^j, one odd radical throughout,",
        "       so the threshold is constant and only the level moves.",
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
