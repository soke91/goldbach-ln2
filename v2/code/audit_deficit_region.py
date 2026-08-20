# -*- coding: utf-8 -*-
r"""
"Never closes" is one inequality, and it is excluded -- inside one family.

WHAT IS AT STAKE

{#rem:deficitshape} tested a single alternative -- a drift dying
exactly where the quadratic put the closure -- and excluded it at
3.28 standard errors.  One point of a region is a thin thing to
exclude, and its B5 had to be reported refuted because it asked
whether that alternative would be *detected* rather than whether the
measurement could be told from it.

The region has a boundary, and it is simpler than the single point.
Write the deficit as the derivative of the cubic,

    deficit(x) = c1 + c2 x + (1/2) c3 x^2,

a parabola in x = log N.  If c3 < 0 it opens downward, runs to minus
infinity, and the deficit reaches zero at some finite N whatever else
is true.  If c3 > 0 it opens upward with its vertex at x = -c2/c3;
with c2 > 0 that vertex is at negative x, below the field entirely,
so past the top of the field the deficit only rises and **never
reaches zero**.  So inside the cubic family the whole "never closes"
region is the single inequality c3 > 0, and excluding it is a
one-sided test on one coefficient.

That is a stronger statement than {#rem:deficitshape} made from the
same number, and it is bounded by the family.  A quartic term can
bend the parabola back up beyond the field and put the region back,
so whether one is resolved is part of the same question and is
tested here rather than assumed away -- which is what
{#rem:shapepower} exists to insist on.

Nothing is measured.  All 156 points are read from POINT markers.

BACKS: Remark {#rem:deficitregion} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  The control.  The cubic fitted here returns the coefficient
      and standard error results/audit_deficit_shape.txt declares,
      to the six decimals printed.
  C2  The premises of the reduction hold: the deficit at the top of
      the field is positive, c2 is positive, and a positive c3 would
      put the parabola's vertex below the field.
  C3  The region is excluded: the one-sided t on c3 is beyond -2, so
      c3 > 0 -- every shape in this family under which the deficit
      never closes -- is out.
  C4  The family is not obviously too small: a quartic term added to
      the fit is unresolved, |t| < 2.
  C5  And it does not overturn C3: with the quartic in the fit, the
      cubic coefficient stays negative and beyond -2.

REFUTATION RULE (fixed before the run)

  C1  REFUTED outside the printed decimals.  Same points, same
      model.  THIS ONE GATES.
  C2  REFUTED by any premise failing.  The reduction to one
      inequality would not hold and C3 would be about nothing.
      THIS ONE GATES.
  C3  REFUTED if the one-sided t does not pass -2, and that failure
      has two forms which mean different things.  **The sign could
      turn positive and resolve**, which would put the measurement
      inside the never-closes region and end the budget route.  **Or
      it could be unresolved** -- too noisy to tell -- in which case
      nothing is excluded and the honest report is that the region
      stands untouched, not that it is empty.
  C4  REFUTED if the quartic resolves.  Then the cubic family is too
      small to carry the reduction, the parabola argument does not
      describe the deficit past the field, and C3's exclusion is
      void whatever its t said.  This is the outcome that would cost
      the most, and it is not the one predicted.
  C5  REFUTED if the cubic coefficient loses its sign or its
      resolution once the quartic is in.  Unresolved here means the
      exclusion depended on the family being cut at three, which is
      an assumption and not a measurement.

  C1 and C2 gate.  C3 to C5 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  Nothing is sampled: one fit is
  compared with another on the same deterministic points, and the
  question is which coefficients the data resolves.  The nulls for
  the underlying quantity were run in lab_primorial_share.py and in
  {#rem:deficitshape}.
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
OUT = os.path.join(RES, "audit_deficit_region.txt")

DEC = 6


def read_points():
    out = {}
    for fn, tag in (("audit_deficit_direct.txt", "deficitdirect"),
                    ("audit_alpha_reach.txt", "alphareach")):
        src = io.open(os.path.join(RES, fn), encoding="utf-8").read()
        for m in re.finditer(r"^POINT %s_(\d+) ([\d.eE+-]+) "
                             r"([\d.eE+-]+)\s*$" % tag, src, re.M):
            out[int(m.group(1))] = (float(m.group(2)),
                                    float(m.group(3)))
    return out


def read_gamma():
    src = io.open(os.path.join(RES, "audit_deficit_shape.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^  cubic gamma ([-+][\d.]+) \+- ([\d.]+), "
                  r"t = ([-+]?[\d.]+)\s*$", src, re.M)
    if not m:
        raise SystemExit("no cubic gamma line in "
                         "results/audit_deficit_shape.txt")
    return float(m.group(1)), float(m.group(2))


def fit(x, y, deg):
    cols = [np.ones_like(x), x]
    fac = [1.0, 1.0]
    for d in range(2, deg + 1):
        cols.append(x ** d / math.factorial(d))
        fac.append(1.0)
    A = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    s2 = float((r ** 2).sum()) / (x.size - A.shape[1])
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, cov, float(np.sqrt((r ** 2).mean()))


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pts = read_points()
    pg, pse = read_gamma()
    rnd = 0.5 * 10.0 ** (-DEC)
    say("read %d POINT markers; nothing is measured here" % len(pts))
    say("  the cubic coefficient %+.8f +- %.8f comes from "
        "results/audit_deficit_shape.txt" % (pg, pse))
    say("PRINTBOUND audit_deficit_region %d %.8f" % (DEC, rnd))

    allN = sorted(pts)
    x = np.log(np.array(allN, dtype=np.float64))
    y = np.log(np.array([pts[N][0] / pts[N][1] for N in allN]))
    say("  the field is %d points over %.4f in log N, topping at "
        "%.4f" % (x.size, x.max() - x.min(), x.max()))
    say("SCALES 1")

    # -------------------------------------------------------------- C1
    say()
    say("C1  the control against {#rem:deficitshape}")
    c3v, cov3, rms3 = fit(x, y, 3)
    g, sg = float(c3v[3]), math.sqrt(float(cov3[3, 3]))
    c1ok = abs(g - pg) <= rnd and abs(sg - pse) <= rnd
    say("  cubic coefficient %+.8f +- %.8f against %+.8f +- %.8f"
        % (g, sg, pg, pse))
    say("  C1 %s   (cap: the printing bound)"
        % ("hold" if c1ok else "REFUTED"))

    # -------------------------------------------------------------- C2
    say()
    say("C2  do the premises of the reduction hold?")
    d1, d2 = float(c3v[1]), float(c3v[2])
    top = float(x.max())
    dtop = d1 + d2 * top + 0.5 * g * top * top
    vert = -d2 / abs(g) if g != 0 else float("nan")
    prem = (dtop > 0.0, d2 > 0.0, vert < x.min())
    c2ok = all(prem)
    say("  the deficit at the top is %+.6f  (needs > 0): %s"
        % (dtop, "yes" if prem[0] else "NO"))
    say("  c2 is %+.8f  (needs > 0): %s"
        % (d2, "yes" if prem[1] else "NO"))
    say("  a positive c3 would put the vertex at x = %.2f, below the "
        "field's %.4f: %s" % (vert, x.min(),
                              "yes" if prem[2] else "NO"))
    say("  C2 %s   (cap: all three)" % ("hold" if c2ok else "REFUTED"))
    if not (c1ok and c2ok):
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- C3
    say()
    say("C3  is the whole never-closes region excluded?")
    t3 = g / sg
    c3ok = t3 < -2.0
    say("  inside the cubic family the region is the single "
        "inequality c3 > 0,")
    say("  because c3 < 0 opens the parabola downward and the "
        "deficit must reach zero,")
    say("  while c3 > 0 puts its vertex below the field and the "
        "deficit only rises")
    say("  c3 = %+.8f +- %.8f, one-sided t = %.2f" % (g, sg, t3))
    say("TSTAT deficit_region_c3 %.2f" % t3)
    if abs(t3) < 2.0:
        say("UNRESOLVED SIGN deficit_region_c3")
    say("SPREAD deficit_region_c3 %.4f" % (x.max() - x.min()))
    say("SCATTER slope_audit_deficit_region %.6f" % rms3)
    say("  C3 %s   (cap: one-sided t = -2)"
        % ("hold" if c3ok else "REFUTED"))

    # --------------------------------------------------------- C4, C5
    say()
    say("C4/C5  is the cubic family big enough to carry the argument?")
    c4v, cov4, rms4 = fit(x, y, 4)
    q, sq = float(c4v[4]), math.sqrt(float(cov4[4, 4]))
    tq = q / sq
    g4, sg4 = float(c4v[3]), math.sqrt(float(cov4[3, 3]))
    t34 = g4 / sg4
    c4ok = abs(tq) < 2.0
    c5ok = g4 < 0.0 and t34 < -2.0
    say("  the quartic coefficient is %+.10f +- %.10f, t = %.2f"
        % (q, sq, tq))
    say("TSTAT deficit_region_quartic %.2f" % tq)
    if abs(tq) < 2.0:
        say("UNRESOLVED SIGN deficit_region_quartic")
    say("SPREAD deficit_region_quartic %.4f" % (x.max() - x.min()))
    say("  r.m.s. residual %.6f against the cubic's %.6f"
        % (rms4, rms3))
    say("  with it in the fit the cubic coefficient reads %+.8f +- "
        "%.8f, t = %.2f" % (g4, sg4, t34))
    say("  C4 %s   (cap: |t| = 2 on the quartic)"
        % ("hold" if c4ok else "REFUTED"))
    say("  C5 %s   (cap: the cubic stays negative beyond -2)"
        % ("hold" if c5ok else "REFUTED"))

    say()
    say("a diagnostic, run after C4 was refuted and registered as "
        "nothing")
    say("  C4 asked whether one more degree resolves. It does. The "
        "question that")
    say("  answer raises is whether the next one does too, and the "
        "one after that --")
    say("  because a family in which every added degree keeps "
        "resolving while the")
    say("  residual barely moves is not a family the data is choosing "
        "from. This is")
    say("  measured below and predicted by nothing; it is reported "
        "because C4's")
    say("  refutation is uninterpretable without it.")
    say("  degree   top coefficient t   r.m.s. residual   resolves?")
    for deg in range(2, 9):
        cc, vv, rr = fit(x, y, deg)
        tt = float(cc[deg]) / math.sqrt(float(vv[deg, deg]))
        say("  %-8d %+-19.2f %-17.6f %s"
            % (deg, tt, rr, "yes" if abs(tt) > 2.0 else "no"))
    say("DEGREES deficit_region 7")

    say()
    say("what this settles")
    if c3ok and c4ok and c5ok:
        say("  every shape in the cubic family under which the "
            "deficit never closes is out")
        say("  at %.2f standard errors, and no quartic bend is "
            "resolved that would put" % abs(t3))
        say("  the region back. The family is still an assumption: "
            "nothing here says the")
        say("  deficit is a polynomial in log N, only that within "
            "the polynomials the data")
        say("  admits, it closes.")
    elif not c4ok:
        say("  a quartic bend is resolved, so the parabola argument "
            "does not describe the")
        say("  deficit past the field and the exclusion is void "
            "whatever C3's t said")
    else:
        say("  the region is not excluded; whether that is because "
            "the measurement lies")
        say("  in it or because the field cannot tell is what C3's "
            "own sign says")

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  C3 %s  C4 %s  C5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (c1ok, c2ok, c3ok, c4ok, c5ok)))

    head = [
        "STATISTIC: the cubic and quartic coefficients of",
        "           log(|sum a|/l2) against log N over the whole",
        "           field, and the premises under which the region",
        "           of shapes where the deficit never reaches zero",
        "           reduces, inside the cubic family, to the single",
        "           inequality c3 > 0.",
        "NULL: none is run and none applies. Nothing is sampled here:",
        "      two fits are compared on the same deterministic",
        "      points and the question is which coefficients the",
        "      data resolves. The nulls for the underlying quantity",
        "      were run in lab_primorial_share.py and in",
        "      {#rem:deficitshape}.",
        "FIELD: N = 2^a 5^b with both a >= 1 and b >= 1 in",
        "       [2e5, 8e9], one coprimality class. Nothing is",
        "       measured here: every |sum a| and l2 is read from a",
        "       POINT marker, the 81 below 1.024e8 from",
        "       results/audit_deficit_direct.txt and the 75 above it",
        "       from results/audit_alpha_reach.txt; the cubic",
        "       coefficient C1 checks against comes from",
        "       results/audit_deficit_shape.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not (c1ok and c2ok):
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
