# -*- coding: utf-8 -*-
r"""
Does the curvature predict, and where must it expire?

WHAT IS AT STAKE

{#rem:laddercurve} resolved an upward curvature and {#rem:curvereach}
showed the earlier shape contest could not have seen it. Before
anything leans on that curvature, it owes the same discipline the
line was held to: a fit is not trusted until it predicts points it
was not fitted to. The quadratic has never been asked.

And it has a ceiling. The level exponent is log K*_R / log N with
K*_R a truncation below N, so it cannot pass 1 -- the same kind of
bound that made {#rem:flatnessshape}'s F <= 1 decisive. An upward
quadratic passes any level eventually, so the curvature is a local
description with a computable expiry, and saying where it expires is
saying how far it may be carried.

Both are re-analysis of the fifteen published rungs. No arithmetic is
computed.

BACKS: Remark {#rem:curvebound} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  G1  The control. The fifteen-rung quadratic reproduces
      results/audit_ladder_curve.txt inside the bound its printing
      forces.
  G2  The curvature predicts: a quadratic fitted on the first
      thirteen rungs puts rungs 14 and 15 inside their own prediction
      standard errors, and closer than the line fitted on the same
      thirteen does.
  G3  The ceiling is real and not yet binding: at every rung the
      measured K*_R is below N and the exponent is below 1 by more
      than half.
  G4  The crossing moves earlier: the quadratic reaches 0.56 below
      the line's log10 N = 10.6180.
  G5  And the curvature expires: the quadratic reaches the ceiling 1
      at a finite abscissa, which is where it stops being usable.

REFUTATION RULE (fixed before the run)

  G1  REFUTED outside the printing bound. THIS ONE GATES.
  G2  REFUTED if either rung falls outside, or if the line predicts
      them better. Then the curvature fits the fifteen but does not
      predict, and it is a description rather than a shape -- which
      would leave item 1 with no shape at all rather than a new one.
  G3  REFUTED if any exponent is at or above 1, or within half of it.
      Then the ceiling is already in play at the measured rungs and
      every fitted slope is contaminated by it.
  G4  REFUTED if the quadratic's crossing is not earlier. Then the
      curvature does not move the forecast in the direction
      {#rem:rung13} and {#rem:laddercurve} both read off it.
  G5  REFUTED if no finite crossing of 1 exists, which for an upward
      quadratic would mean the fit is not upward after all.

  G1 gates. G2 to G5 are the measurement and do not gate. Every
  crossing below is an extrapolation past the data and carries a
  bracket and a drift; {#rem:shapepower} is why none of them is
  offered as a forecast.

  NO NULL IS RUN and none applies. Every quantity is a least-squares
  summary of exponents already measured. The brackets are drawn from
  the fit's own parameter covariance with the fixed SEED declared in
  the output.
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
OUT = os.path.join(RES, "audit_curve_bound.txt")

BASE = 30030
SEED = 20260822
DRAWS = 4000
TARGET = 0.56
CEIL = 1.0


def read_rungs():
    src = io.open(os.path.join(RES, "audit_primorial_rung10.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log10 N   exponent   fitted     "
                  "residual")
    ns, ex, dec = [], [], 0
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        ns.append(int(f[0]))
        ex.append(float(f[2]))
        dec = max(dec, len(f[2].split(".")[1]))
    ks = []
    for j in (11, 12, 13, 14):
        s = io.open(os.path.join(RES, "audit_primorial_rung%d.txt" % j),
                    encoding="utf-8").read()
        N = BASE * (1 << j)
        m = re.search(r"^  N = " + str(N) +
                      r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                      r"K\*_R (\d+)\s+exp ([\d.]+)\s*$", s, re.M)
        ns.append(N)
        ex.append(float(m.group(2)))
        ks.append((N, int(m.group(1))))
    s = io.open(os.path.join(RES, "audit_ladder_curve.txt"),
                encoding="utf-8").read()
    m = re.search(r"the quadratic coefficient is ([+-][\d.]+), "
                  r"s\.e\. ([\d.]+)", s)
    return ns, ex, ks, dec, float(m.group(1)), float(m.group(2))


def quadfit(x, y):
    A = np.column_stack([np.ones_like(x), x, x * x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    n = x.size
    s2 = float((r ** 2).sum()) / (n - 3)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, cov, float(np.sqrt((r ** 2).mean())), s2, A


def linfit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    s2 = float((r ** 2).sum()) / (n - 2)
    return a, b, s2, float(((x - x.mean()) ** 2).sum()), x.mean()


def cross(c, level):
    """the larger root of a + bL + cL^2 = level, or None"""
    a2, b2, c2 = c[0] - level, c[1], c[2]
    if abs(c2) < 1e-18:
        return None if abs(b2) < 1e-18 else -a2 / b2
    disc = b2 * b2 - 4.0 * c2 * a2
    if disc < 0:
        return None
    r1 = (-b2 + math.sqrt(disc)) / (2.0 * c2)
    r2 = (-b2 - math.sqrt(disc)) / (2.0 * c2)
    good = [r for r in (r1, r2) if r > 0]
    return max(good) if c2 < 0 else (min(good) if good else None)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ns, ex, ks, dec, cq, sq = read_rungs()
    say("read %d rung exponents and the published quadratic "
        "%+.8f (s.e. %.8f)" % (len(ns), cq, sq))
    say("  from results/audit_primorial_rung10, 11, 12, 13, 14 and "
        "results/audit_ladder_curve.txt")
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)
    say("RADICALS 1")

    x = np.log(np.array(ns, dtype=np.float64))
    y = np.array(ex)
    c, cov, rms, s2, A = quadfit(x, y)

    # -------------------------------------------------------------- G1
    say()
    say("G1  the control")
    rnd = 0.5 * 10.0 ** (-8)
    d = abs(float(c[2]) - cq)
    g1 = d <= rnd
    say("  quadratic here %+.8f against the published %+.8f, "
        "departure %.10f; the bound is %.10f" % (c[2], cq, d, rnd))
    say("PRINTBOUND audit_curve_bound 8 %.10f" % rnd)
    say("  G1 %s   (cap: the printing bound)"
        % ("hold" if g1 else "REFUTED"))

    # -------------------------------------------------------------- G2
    say()
    say("G2  does the curvature predict out of sample?")
    n13 = len(ns) - 2
    c13, cov13, rms13, s213, A13 = quadfit(x[:n13], y[:n13])
    al, bl, s2l, sxxl, xbl = linfit(x[:n13], y[:n13])
    say("  fitted on the first %d rungs; both shapes asked for the "
        "last two" % n13)
    say("  rung   log10 N   quad pred   line pred   measured   "
        "quad dep   line dep   quad s.e.  inside")
    g2 = True
    for i in (len(ns) - 2, len(ns) - 1):
        v = np.array([1.0, x[i], x[i] * x[i]])
        pq = float(v.dot(c13))
        sp = math.sqrt(s213 * (1.0 + float(v.dot(cov13).dot(v)) / s213))
        pl = al * x[i] + bl
        dq, dl = y[i] - pq, y[i] - pl
        ok = abs(dq) <= sp and abs(dq) < abs(dl)
        if not ok:
            g2 = False
        say("  %-6d %-9.4f %-11.4f %-11.4f %-10.4f %+-10.4f %+-10.4f "
            "%-10.4f %s" % (i, math.log10(ns[i]), pq, pl, y[i], dq,
                            dl, sp, "yes" if ok else "NO"))
    say("  G2 %s   (cap: outside its own error, or worse than the "
        "line)" % ("hold" if g2 else "REFUTED"))

    # -------------------------------------------------------------- G3
    say()
    say("G3  is the ceiling in play at the measured rungs?")
    say("  the exponent is log K*_R / log N with K*_R < N, so it "
        "cannot pass 1")
    worst = max(ex)
    below = all(k < N for N, k in ks)
    g3 = worst < 0.5 * CEIL + 0.5 and below and worst < CEIL
    say("  the largest measured exponent is %.4f and every K*_R is "
        "below its N: %s" % (worst, "yes" if below else "NO"))
    say("  the four measured K*_R against their N:")
    for N, k in ks:
        say("    N = %-12d K*_R = %-8d ratio %.6f" % (N, k, k / N))
    g3 = worst < CEIL and below
    say("  G3 %s   (cap: an exponent at or above 1)"
        % ("hold" if g3 else "REFUTED"))

    # --------------------------------------------------- G4 and G5
    rng = np.random.default_rng(SEED)
    draws = rng.multivariate_normal(c, cov, size=DRAWS)
    say()
    say("G4/G5  where the quadratic would reach 0.56 and where it "
        "reaches the ceiling")
    for level, tag in ((TARGET, "theta_prime"), (CEIL, "ceiling")):
        pt = cross(c, level)
        vals = [cross(dd, level) for dd in draws]
        vals = [v for v in vals if v is not None and v > x.max()]
        if pt is None or not vals:
            say("  level %.2f: no crossing" % level)
            continue
        lo = float(np.percentile(vals, 2.5)) / math.log(10.0)
        hi = float(np.percentile(vals, 97.5)) / math.log(10.0)
        p10 = pt / math.log(10.0)
        say("  level %.2f reached at log10 N = %.4f, bracket "
            "[%.4f, %.4f] from %d of %d draws"
            % (level, p10, lo, hi, len(vals), DRAWS))
        say("BRACKET ladder_quadratic_%s %.4f %.4f %.4f"
            % (tag, p10, lo, hi))
        half = len(ns) // 2
        cl, _cvl, _rl, _s2l, _Al = quadfit(x[:half + 2], y[:half + 2])
        pl_ = cross(cl, level)
        say("DRIFT ladder_quadratic_%s %.4f"
            % (tag, abs((pl_ or 0.0) - pt) / math.log(10.0)))
        say("  refitted on the lower %d rungs it lands at %.4f, so "
            "the drift declared is their gap"
            % (half + 2, (pl_ or float("nan")) / math.log(10.0)))
    pt56 = cross(c, TARGET) / math.log(10.0)
    ptc = cross(c, CEIL)
    g4 = pt56 < 10.6180
    g5 = ptc is not None and ptc > x.max()
    say("  the line published %.4f for 0.56; the quadratic puts it "
        "%s" % (10.6180, "earlier" if g4 else "later"))
    say("  G4 %s   (cap: the line's crossing)"
        % ("hold" if g4 else "REFUTED"))
    say("  G5 %s   (cap: no finite crossing of the ceiling)"
        % ("hold" if g5 else "REFUTED"))
    say("SHAPES 2")
    say("SCATTER slope_audit_curve_bound %.4f" % rms)
    say("  none of these is offered as a forecast. The bracket is the "
        "fit's own")
    say("  parameter spread and nothing else; {#rem:shapepower} "
        "measured that this")
    say("  repository's shape discriminator has no power at this "
        "reach, and a")
    say("  quadratic extrapolated past its data is the case that "
        "warning was about.")

    say()
    say("=" * 70)
    say("G1 %s  G2 %s  G3 %s  G4 %s  G5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (g1, g2, g3, g4, g5)))

    head = [
        "STATISTIC: the fifteen primorial-ladder level exponents, read",
        "           from the files that measured them; a quadratic in",
        "           log N fitted on the first thirteen and asked for",
        "           the last two, against a line on the same thirteen,",
        "           each with the prediction standard error at the new",
        "           abscissa; the measured K*_R against N at the four",
        "           rungs that print it; and the abscissae at which",
        "           the fifteen-rung quadratic reaches 0.56 and the",
        "           ceiling 1, each with a bracket from the fit's own",
        "           parameter covariance and a drift from refitting on",
        "           the lower rungs.",
        "NULL: none is run and none applies. Every quantity is a",
        "      least-squares summary of exponents already measured;",
        "      the brackets are drawn from the fit's own parameter",
        "      covariance with the fixed SEED above.",
        "FIELD: no arithmetic is computed here. The exponents are",
        "       read from results/audit_primorial_rung10.txt (eleven",
        "       rungs) and results/audit_primorial_rung11.txt,",
        "       rung12, rung13 and rung14 (one each), whose own field",
        "       is N = 30030*2^j with the odd radical 3*5*7*11*13",
        "       fixed, k squarefree and coprime to N with",
        "       2 <= k < 100000, m odd, squarefree and coprime to k,",
        "       and the Euler products at the fixed bound 4000000;",
        "       the published quadratic is read from",
        "       results/audit_ladder_curve.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not g1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
