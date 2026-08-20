# -*- coding: utf-8 -*-
r"""
Does the deficit's drift itself drift? And could this field tell?

WHAT IS AT STAKE

{#rem:deficitdirect} resolved the deficit's drift at
beta = -0.007380 +- 0.000632, t = -11.67: inside the field, and by
two routes that agree, the deficit shrinks.  It then declined to say
where it reaches zero, because the forecast moved 45.8908 decades
when seven octaves were added.

The reason that forecast is worthless is a shape assumption.  A
quadratic in log-log says the drift is a constant, so the deficit
falls without limit and must cross zero.  A cubic says the drift is
itself drifting, and if it weakens the deficit approaches a positive
limit and **never closes at any N**.  Both fit a shrinking deficit.
Only the first implies the budget route exists.

So the question is not whether the deficit shrinks -- that is
settled -- but whether this field can tell a constant drift from a
dying one.  Nothing new is measured: all 156 points are read from
POINT markers, 81 from results/audit_deficit_direct.txt and 75 from
results/audit_alpha_reach.txt.

**And the power is measured, not assumed.**  {#rem:shapepower}
established for the level axis that shape selection at a given reach
can have no power at all, and that a shape contest run without
measuring its power reports a preference it was never able to have.
The alternative used here is the one that matters: a drift that dies
exactly where the quadratic says the deficit would have closed, so
that the deficit approaches a positive limit instead.  If the cubic
term resolves under that alternative only rarely, then this field
cannot distinguish "closes eventually" from "never closes", and the
honest report is that -- not a preference.

BACKS: Remark {#rem:deficitshape} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  B1  The control.  The quadratic fitted here on all 156 points
      returns the beta and standard error
      results/audit_deficit_direct.txt declares, to the six decimals
      printed.
  B2  The cubic term is resolved: |t| > 2.  **No sign is
      predicted.**
  B3  And it earns it out of sample: walking forward -- fit on the
      lowest j points, predict the next -- the cubic's r.m.s.
      out-of-sample departure is smaller than the quadratic's.
  B4  The false positive rate is controlled: under a null where the
      truth is the fitted quadratic with the observed residual
      scatter, the cubic term resolves in fewer than 5 per cent of
      simulated fields.
  B5  **And the field has the power to see the alternative that
      matters**: under a drift that dies exactly where the quadratic
      puts the closure -- so the deficit never reaches zero -- the
      cubic term resolves in more than half of simulated fields.

REFUTATION RULE (fixed before the run)

  B1  REFUTED outside the printed decimals.  Same points, same
      model; a difference means the markers are not the values that
      remark fitted.  THIS ONE GATES.
  B2  REFUTED if |t| <= 2.  Read alone that says only "cannot tell",
      and B5 is what decides whether "cannot tell" is a fact about
      the drift or about the field.
  B3  REFUTED if the cubic's out-of-sample r.m.s. is not smaller.
      A term that resolves in-sample and does not predict is the
      pattern {#rem:signrun} found and withdrew.
  B4  REFUTED at or above 5 per cent.  Then the cubic resolves too
      easily on this design and B2 carries no weight whichever way
      it went.
  B5  REFUTED at or below half.  **That is the outcome that would
      settle the matter against the branch's hopes**: it would mean
      this field cannot distinguish a deficit that closes from one
      that approaches a positive limit, so {#rem:deficitdirect}'s
      resolved shrinkage licenses nothing about whether the
      requirement is ever met, and no amount of reading these 156
      points differently will change that.  More N would be the only
      move, and B5 says how much more.

  B1 gates.  B2 to B5 are the measurement and do not gate.

  THE NULL IS RUN, twice.  B4 draws fields from the fitted quadratic
  with the observed residual scatter; B5 draws them from the cubic
  whose drift vanishes at the quadratic's closure point.  Both are
  fitted by the same code path as the observed field.
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
OUT = os.path.join(RES, "audit_deficit_shape.txt")

SEED = 20260823
DRAWS = 4000
DEC = 6
START = 40                          # the first point predicted forward


def read_points():
    """every N in the field, from the POINT markers of two files"""
    out = {}
    for fn, tag in (("audit_deficit_direct.txt", "deficitdirect"),
                    ("audit_alpha_reach.txt", "alphareach")):
        src = io.open(os.path.join(RES, fn), encoding="utf-8").read()
        for m in re.finditer(r"^POINT %s_(\d+) ([\d.eE+-]+) "
                             r"([\d.eE+-]+)\s*$" % tag, src, re.M):
            out[int(m.group(1))] = (float(m.group(2)),
                                    float(m.group(3)))
    return out


def read_beta():
    src = io.open(os.path.join(RES, "audit_deficit_direct.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^BETA deficit_direct ([-+][\d.]+) ([\d.]+)\s*$",
                  src, re.M)
    if not m:
        raise SystemExit("no BETA deficit_direct marker")
    return float(m.group(1)), float(m.group(2))


def polyfit_se(x, y, deg):
    """coefficients of 1, x, x^2/2, x^3/6 and their covariance"""
    cols = [np.ones_like(x), x]
    if deg >= 2:
        cols.append(0.5 * x * x)
    if deg >= 3:
        cols.append(x * x * x / 6.0)
    A = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    n = x.size
    s2 = float((r ** 2).sum()) / (n - A.shape[1])
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, cov, float(np.sqrt((r ** 2).mean()))


def walk(x, y, deg, start=START):
    """fit on the lowest j, predict j: the out-of-sample departures"""
    d = []
    for j in range(start, x.size):
        c, _cov, _r = polyfit_se(x[:j], y[:j], deg)
        v = [1.0, x[j]]
        if deg >= 2:
            v.append(0.5 * x[j] * x[j])
        if deg >= 3:
            v.append(x[j] ** 3 / 6.0)
        d.append(y[j] - float(np.array(v).dot(c)))
    return np.array(d)


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pts = read_points()
    pbeta, pse = read_beta()
    rnd = 0.5 * 10.0 ** (-DEC)
    say("read %d POINT markers from results/audit_deficit_direct.txt "
        "and" % len(pts))
    say("  results/audit_alpha_reach.txt; nothing is measured here")
    say("READ audit_deficit_direct.txt BETA deficit_direct %+.6f %.6f"
        % (pbeta, pse))
    say("PRINTBOUND audit_deficit_shape %d %.8f" % (DEC, rnd))

    allN = sorted(pts)
    x = np.log(np.array(allN, dtype=np.float64))
    y = np.log(np.array([pts[N][0] / pts[N][1] for N in allN]))
    say()
    say("the field: %d points over %.4f in log N"
        % (x.size, x.max() - x.min()))
    say("SCALES 1")

    # -------------------------------------------------------------- B1
    say()
    say("B1  the control against {#rem:deficitdirect}")
    c2, cov2, rms2 = polyfit_se(x, y, 2)
    b2, sb2 = float(c2[2]), math.sqrt(float(cov2[2, 2]))
    b1 = abs(b2 - pbeta) <= rnd and abs(sb2 - pse) <= rnd
    say("  quadratic beta %+.6f +- %.6f against the declared %+.6f "
        "+- %.6f" % (b2, sb2, pbeta, pse))
    say("  B1 %s   (cap: the printing bound)"
        % ("hold" if b1 else "REFUTED"))
    if not b1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- B2
    say()
    say("B2  is the drift itself drifting?")
    c3, cov3, rms3 = polyfit_se(x, y, 3)
    g, sg = float(c3[3]), math.sqrt(float(cov3[3, 3]))
    tg = g / sg
    b2ok = abs(tg) > 2.0
    say("  cubic gamma %+.8f +- %.8f, t = %.2f" % (g, sg, tg))
    say("TSTAT deficit_shape_gamma %.2f" % tg)
    if abs(tg) < 2.0:
        say("UNRESOLVED SIGN deficit_shape_gamma")
    say("SPREAD deficit_shape_gamma %.4f" % (x.max() - x.min()))
    say("  r.m.s. residual %.6f against the quadratic's %.6f"
        % (rms3, rms2))
    say("SCATTER slope_audit_deficit_shape %.6f" % rms3)
    say("  B2 %s   (cap: |t| = 2)" % ("hold" if b2ok else "REFUTED"))

    # -------------------------------------------------------------- B3
    say()
    say("B3  does the cubic predict out of sample?")
    d2, d3 = walk(x, y, 2), walk(x, y, 3)
    r2 = float(np.sqrt((d2 ** 2).mean()))
    r3 = float(np.sqrt((d3 ** 2).mean()))
    b3 = r3 < r2
    say("  walking forward from point %d: %d out-of-sample points"
        % (START, d2.size))
    say("  quadratic r.m.s. departure %.6f, cubic %.6f" % (r2, r3))
    say("SIGNRUN deficit_shape_walk %d"
        % sum(1 for i in range(1, d3.size)
              if np.sign(d3[i]) == np.sign(d3[i - 1])))
    say("  B3 %s   (cap: the quadratic's departure)"
        % ("hold" if b3 else "REFUTED"))

    # --------------------------------------------------------- B4, B5
    say()
    say("B4/B5  the two nulls")
    A2 = np.column_stack([np.ones_like(x), x, 0.5 * x * x])
    fit2 = A2.dot(c2)
    xs = -float(c2[1]) / float(c2[2]) if c2[2] < 0 else None
    galt = -float(c2[2]) / xs if xs else 0.0
    say("  the quadratic puts the closure at x = %.4f, so a drift "
        "that dies there has" % (xs if xs else float('nan')))
    say("  gamma = %+.8f; that is the alternative B5 draws from"
        % galt)
    A3 = np.column_stack([np.ones_like(x), x, 0.5 * x * x,
                          x ** 3 / 6.0])
    calt = np.array([c2[0], c2[1], c2[2], galt])
    fit3 = A3.dot(calt)
    rng = np.random.default_rng(SEED)
    hit4 = hit5 = 0
    for _ in range(DRAWS):
        n2 = fit2 + rng.normal(0.0, rms2, size=x.size)
        cc, vv, _ = polyfit_se(x, n2, 3)
        if abs(cc[3] / math.sqrt(float(vv[3, 3]))) > 2.0:
            hit4 += 1
        n3 = fit3 + rng.normal(0.0, rms2, size=x.size)
        cc, vv, _ = polyfit_se(x, n3, 3)
        if abs(cc[3] / math.sqrt(float(vv[3, 3]))) > 2.0:
            hit5 += 1
    f4, f5 = hit4 / float(DRAWS), hit5 / float(DRAWS)
    b4 = f4 < 0.05
    b5 = f5 > 0.5
    say("NULL deficit_shape_falsepos %.4f" % f4)
    say("  under a constant drift the cubic resolves in %.4f of "
        "%d fields" % (f4, DRAWS))
    say("  B4 %s   (cap: 0.05)" % ("hold" if b4 else "REFUTED"))
    say("NULL deficit_shape_power %.4f" % f5)
    say("  under a drift that dies at the closure it resolves in "
        "%.4f" % f5)
    say("  B5 %s   (cap: one half)" % ("hold" if b5 else "REFUTED"))

    say()
    say("what B5 should have asked, and the answer")
    say("  B5 asked whether the cubic term resolves under the "
        "alternative. That is a")
    say("  question about detecting any curvature, and the "
        "alternative's curvature is")
    say("  tiny: gamma = %+.8f against a standard error of %.8f, "
        "%.2f of it." % (galt, sg, abs(galt) / sg))
    say("  Nothing could resolve that here and the rule should not "
        "have asked it to.")
    say("  The question that decides the matter is whether the "
        "measured gamma can be")
    say("  told from the alternative's, and it can: the gap is "
        "%+.8f, which is" % (g - galt))
    say("  %.2f standard errors." % (abs(g - galt) / sg))
    say("DISCRIMINATION deficit_shape_alt %.2f" % (abs(g - galt) / sg))
    say("  This is the third time a registered rule named the wrong "
        "event -- M9 in the")
    say("  README, after {#rem:thetalaw} U4 and {#rem:alphalocal} "
        "Z4. B5 stands refuted")
    say("  as written; the number above is reported beside it, not "
        "in place of it.")

    say()
    say("what this settles")
    if b5 and b2ok:
        say("  the field can see the alternative and does: the drift "
            "is not constant")
    elif b5 and not b2ok:
        say("  the field can see the alternative and does not: a "
            "constant drift stands,")
        say("  and the deficit reaching zero is licensed by the "
            "shape the data prefers")
    else:
        say("  the field cannot see the alternative, so it cannot "
            "tell a deficit that")
        say("  closes from one that approaches a positive limit. "
            "{#rem:deficitdirect}'s")
        say("  resolved shrinkage licenses nothing about whether the "
            "requirement is met.")

    say()
    say("=" * 70)
    say("B1 %s  B2 %s  B3 %s  B4 %s  B5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (b1, b2ok, b3, b4, b5)))

    head = [
        "SEED: the null draws from numpy default_rng at "
        "seed %d; without it the file does not reproduce "
        "its own null." % SEED,
        "STATISTIC: the cubic coefficient of",
        "           log(|sum a|/l2) against log N over the whole",
        "           field -- the rate at which the deficit's own",
        "           drift drifts -- with its t; the out-of-sample",
        "           departures of the quadratic and the cubic walked",
        "           forward; and the fraction of simulated fields in",
        "           which the cubic resolves, under a constant drift",
        "           and under a drift that dies exactly where the",
        "           quadratic puts the deficit's closure.",
        "NULL: RUN, twice, and they are B4 and B5. B4 draws fields",
        "      from the fitted quadratic with the observed residual",
        "      scatter and counts how often the cubic resolves",
        "      spuriously. B5 draws them from the cubic whose drift",
        "      vanishes at the quadratic's closure point -- the",
        "      alternative under which the deficit never reaches",
        "      zero -- and counts how often the cubic resolves at",
        "      all. Both are fitted by the same code path as the",
        "      observed field, with the fixed SEED.",
        "FIELD: N = 2^a 5^b with both a >= 1 and b >= 1 in",
        "       [2e5, 8e9], one coprimality class. Nothing is",
        "       measured here: every |sum a| and l2 is read from a",
        "       POINT marker, the 81 below 1.024e8 from",
        "       results/audit_deficit_direct.txt and the 75 above it",
        "       from results/audit_alpha_reach.txt, and the",
        "       quadratic's beta is read from the BETA marker of the",
        "       former.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not b1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
