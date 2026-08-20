# -*- coding: utf-8 -*-
r"""
The deficit's drift, fitted on the points instead of on octave summaries.

WHAT IS AT STAKE

[rem:alphareach] carried the sign axis seven octaves further and left
the question where it found it: the deficit's drift is
-0.007519 +- 0.007837, unresolved, and dropping any one octave moves
it across zero.  It also named the reason the errors are what they
are, and it is not the field's size.

**The analysis summarises each octave into one number and then fits
sixteen numbers.**  Each octave's exponent is itself a two-parameter
fit on eleven to fourteen points, so a hundred and fifty-six
measurements become sixteen, and the drift is estimated from those
with thirteen degrees of freedom.  Nothing about the arithmetic
requires that.  The drift is the coefficient of the quadratic term in

    log(|sum a| / l2) = c + alpha*x + (1/2)*beta*x^2,   x = log N,

and fitting that directly uses every point, with a hundred and
fifty-three degrees of freedom.

Tried on the seventy-five points of {#rem:alphareach} alone, the
direct fit's standard error on beta was 5.5 times smaller than the
octave route's on the same data.  Over the whole field the reduction
should be larger, because the octave route's lever arm is the spread
of sixteen octave midpoints while the direct fit's is the spread of
the points themselves.

Nothing new is measured at large N.  The seventy-five points above
1.024e8 are read from the POINT markers {#rem:alphareach} printed for
exactly this purpose; only the eighty-one published points below it
are recomputed, which costs minutes.  Twelve of those eighty-one are
also POINT markers -- the control octave -- so the overlap is a
control that costs nothing.

**What is bought and what is risked.**  Bought: the drift may become
resolved for the first time on this axis.  Risked: a quadratic in
log-log assumes the drift is linear in log N, where the octave route
assumes nothing about its shape.  The two are reported side by side
and A5 is the check that they do not contradict; if they do, the
linear-drift model is wrong and that is the finding rather than the
tighter error.

No sign is predicted for beta.  The seventy-five-point trial gave
-0.002243, and the octave route gave +0.014678 on nine octaves and
-0.007519 on sixteen.  A quantity whose sign has moved three times
under three analyses is one to measure, not to guess.

BACKS: Remark {#rem:deficitdirect} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  A1  The overlap control.  The twelve points that are both
      recomputed here and printed as POINT markers by
      results/audit_alpha_reach.txt agree in the ten digits those
      markers carry.
  A2  The field control.  A straight line fitted to
      log(|sum a|/l2) on the eighty-one published points returns the
      global deficit results/audit_denominator.txt reports, to the
      precision that file prints it.
  A3  The direct fit is at least five times sharper: its standard
      error on beta is below a fifth of the 0.007837 the octave
      route gives over sixteen octaves.
  A4  And beta is resolved: |t| > 2.  **No sign is predicted.**
  A5  The two routes agree: the direct fit's beta lies within the
      octave route's own standard error of the octave route's slope.

REFUTATION RULE (fixed before the run)

  A1  REFUTED by a single point differing in the printed digits.
      Then the POINT markers are not the values this route produces
      and nothing below may use them.  THIS ONE GATES.
  A2  REFUTED outside the published precision.  Then this is not the
      quantity {#rem:alphalocal} and {#rem:alphareach} measured.
      THIS ONE GATES.
  A3  REFUTED at or above a fifth.  The octave summary would then be
      losing little, the reach really is the binding constraint on
      this axis, and the only remaining move is more N.
  A4  REFUTED if |t| <= 2, and the meaning depends on A3.  With A3
      holding it would say the drift is genuinely below what a
      hundred and fifty-six points can see, which is a stronger
      exclusion than {#rem:alphareach} could state.  With A3 also
      refuted it says only that nothing has changed.
  A5  REFUTED if the two disagree by more than the octave route's
      error.  **That is the outcome worth having and it is not the
      one predicted**: it would mean the deficit does not drift
      linearly in log N, so neither route measures a single drift,
      and the question has to be asked of a shape rather than a
      slope.  A tighter error on the wrong model is worth nothing.

  A1 and A2 gate.  A3 to A5 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  |sum a| and its l2 norm are
  deterministic once N is fixed; there is no sampling noise and no
  background to detect against.  The coin arms for the sign axis were
  run in lab_primorial_share.py.
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
OUT = os.path.join(RES, "audit_deficit_direct.txt")

SEED = 20260823
DRAWS = 4000


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


AR = module("audit_alpha_reach")
SPL = AR.SPL
LO, HI = AR.LO, AR.HI


def read_points():
    """the POINT markers {#rem:alphareach} left, N -> (|sum a|, l2)"""
    src = io.open(os.path.join(RES, "audit_alpha_reach.txt"),
                  encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^POINT alphareach_(\d+) ([\d.eE+-]+) "
                         r"([\d.eE+-]+)\s*$", src, re.M):
        out[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    return out


def octave_rows():
    """the sixteen octave summaries, from marker lines in both files"""
    out = []
    for fn, tag, js in (("audit_alpha_local.txt", "alphalocal",
                         range(9)),
                        ("audit_alpha_reach.txt", "alphareach",
                         range(9, 16))):
        src = io.open(os.path.join(RES, fn), encoding="utf-8").read()
        for j in js:
            m = re.search(r"^OCTAVE %s_%d ([\d.]+) ([-+][\d.]+) "
                          r"([\d.]+) ([-+][\d.]+) ([\d.]+)\s*$"
                          % (tag, j), src, re.M)
            if m:
                out.append((fn, "OCTAVE %s_%d" % (tag, j),
                            tuple(float(g) for g in m.groups())))
    return out


def read_published_deficit():
    """the global deficit and its printed precision"""
    src = io.open(os.path.join(RES, "audit_alpha_local.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^  and e\(l2\) here [-+][\d.]+, so the "
                  r"whole-field deficit is \+([\d.]+)\s*$",
                  src, re.M)
    p = re.search(r"^PRINTBOUND audit_alpha_local \d+ ([\d.]+)\s*$",
                  src, re.M)
    if not (m and p):
        raise SystemExit("no whole-field deficit line or PRINTBOUND "
                         "in results/audit_alpha_local.txt")
    return float(m.group(1)), float(p.group(1))


def quadfit(x, y):
    A = np.column_stack([np.ones_like(x), x, 0.5 * x * x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    n = x.size
    s2 = float((r ** 2).sum()) / (n - 3)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, cov, float(np.sqrt((r ** 2).mean()))


def linefit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), se, float(np.sqrt((r ** 2).mean()))


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pts = read_points()
    gdef, gprec = read_published_deficit()
    rows = octave_rows()
    say("read %d POINT markers from results/audit_alpha_reach.txt, "
        "and the %d" % (len(pts), len(rows)))
    say("  octave summaries the two files declare, each as a whole "
        "marker line:")
    for fn, lab, r in rows:
        say("READ %s %s %.4f %+.6f %.6f %+.6f %.6f" % ((fn, lab) + r))
    om = np.array([r[2][0] for r in rows])
    od = np.array([r[2][3] for r in rows])
    oct_slope, oct_se, _orms = linefit(om, od)
    say("  refitting those sixteen gives the octave route's drift "
        "%+.6f +- %.6f" % (oct_slope, oct_se))
    say("  and the global deficit %.6f to %.3f comes from "
        "results/audit_alpha_local.txt" % (gdef, gprec))

    old = AR.family(LO, HI)
    old81 = set(old)
    say()
    say("recomputing the %d published points below %d; the %d above "
        "are read" % (len(old), HI, len(pts) - 12))
    kind, mu = AR.kind_and_mu(HI)
    pv, plg = AR.power_table(HI)
    say("BYTES resident_arrays %d" % (kind.nbytes + mu.nbytes))

    S, L = {}, {}
    for N in old:
        _ks, a = AR.weighted(N, kind, mu, pv, plg)
        S[N] = abs(float(a.sum()))
        L[N] = float(np.sqrt((a ** 2).sum()))
        say("POINT deficitdirect_%d %.10e %.10e" % (N, S[N], L[N]))
    del kind, mu
    say("  the %d points below %d now carry markers too, so the whole"
        % (len(old), HI))
    say("  field is readable without measuring any of it again")

    # -------------------------------------------------------------- A1
    say()
    say("A1  the overlap: points computed here and printed there")
    over = sorted(set(old) & set(pts))
    worst, wn = 0.0, 0
    for N in over:
        for got, ref in ((S[N], pts[N][0]), (L[N], pts[N][1])):
            d = abs(got - ref) / max(abs(ref), 1e-300)
            if d > worst:
                worst, wn = d, N
    a1 = worst < 1e-9
    say("  %d points overlap; the worst relative difference is %.3e "
        "at N = %d" % (len(over), worst, wn))
    say("  A1 %s   (cap: the ten digits the markers carry)"
        % ("hold" if a1 else "REFUTED"))

    # -------------------------------------------------------------- A2
    say()
    say("A2  the field control on the published points")
    xo = np.log(np.array(old, dtype=np.float64))
    yo = np.log(np.array([S[N] / L[N] for N in old]))
    d0, se0, _r = linefit(xo, yo)
    a2 = abs(d0 - gdef) <= gprec
    say("  a line on the %d published points gives a deficit of "
        "%.6f (s.e. %.6f)" % (len(old), d0, se0))
    say("  against the published %.6f measured to %.3f"
        % (gdef, gprec))
    say("  A2 %s   (cap: the published precision)"
        % ("hold" if a2 else "REFUTED"))
    if not (a1 and a2):
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------- the whole field
    allN = sorted(set(old) | set(pts))
    for N in pts:
        if N not in S:
            S[N], L[N] = pts[N]
    x = np.log(np.array(allN, dtype=np.float64))
    y = np.log(np.array([S[N] / L[N] for N in allN]))
    say()
    say("the whole field: %d points over %.4f in log N"
        % (x.size, x.max() - x.min()))
    say("SCALES 1")

    # -------------------------------------------------------------- A3
    say()
    say("A3  is the direct fit sharper?")
    c, cov, rms = quadfit(x, y)
    beta, sb = float(c[2]), math.sqrt(float(cov[2, 2]))
    a3 = sb < oct_se / 5.0
    say("  direct quadratic: beta %+.6f +- %.6f, r.m.s. residual "
        "%.6f" % (beta, sb, rms))
    say("  the octave route over sixteen octaves gave %+.6f +- %.6f"
        % (oct_slope, oct_se))
    say("  the ratio of the errors is %.1f" % (oct_se / sb))
    say("  A3 %s   (cap: a fifth of the octave error)"
        % ("hold" if a3 else "REFUTED"))

    # -------------------------------------------------------------- A4
    say()
    say("A4  is the drift resolved?")
    t = beta / sb
    a4 = abs(t) > 2.0
    say("  beta = %+.6f, t = %.2f" % (beta, t))
    say("BETA deficit_direct %+.6f %.6f" % (beta, sb))
    say("TSTAT deficit_direct_beta %.2f" % t)
    if abs(t) < 2.0:
        say("UNRESOLVED SIGN deficit_direct_beta")
    say("SPREAD deficit_direct_beta %.4f" % (x.max() - x.min()))
    say("SCATTER slope_audit_deficit_direct %.6f" % rms)
    say("  A4 %s   (cap: |t| = 2)" % ("hold" if a4 else "REFUTED"))

    # -------------------------------------------------------------- A5
    say()
    say("A5  do the two routes agree?")
    gap = abs(beta - oct_slope)
    a5 = gap <= oct_se
    say("  direct %+.6f against octave %+.6f: a gap of %.6f against "
        "the octave error %.6f" % (beta, oct_slope, gap, oct_se))
    say("  A5 %s   (cap: the octave route's own error)"
        % ("hold" if a5 else "REFUTED"))

    # ------------------------------------------------- the consequence
    say()
    say("what this says about the budget route")
    alpha = float(c[1])
    say("  the fitted local deficit is %+.6f %+.6f x, so it reaches "
        "zero at" % (alpha, beta))
    if beta >= 0.0 or not a4:
        say("  no x, or none the sign supports; nothing is forecast")
        say("  the drift is not resolved negative, so no closure "
            "follows at any N")
    else:
        xs = -alpha / beta
        rng = np.random.default_rng(SEED)
        dr = rng.multivariate_normal([alpha, beta],
                                     cov[1:3, 1:3], size=DRAWS)
        vv = [-d[0] / d[1] for d in dr if d[1] < 0]
        vv = [v / math.log(10.0) for v in vv if v > x.max()]
        # the same fit on the published field alone, to measure how
        # far the forecast moved when the seven octaves were added
        m81 = np.array([n in old81 for n in allN])
        c8, cv8, _r8 = quadfit(x[m81], y[m81])
        x8 = (-float(c8[1]) / float(c8[2])
              if float(c8[2]) < 0 else float("nan"))
        say("  x* = %.4f in log N, which is 10^%.4f"
            % (xs, xs / math.log(10.0)))
        if vv:
            lo = float(np.percentile(vv, 2.5))
            hi = float(np.percentile(vv, 97.5))
            say("BRACKET deficit_closure %.4f %.4f %.4f"
                % (xs / math.log(10.0), lo, hi))
            say("  bracket [%.4f, %.4f] from %d of %d draws"
                % (lo, hi, len(vv), DRAWS))
        else:
            say("  no draw puts the closure above the field; no "
                "bracket is quoted")
        if x8 == x8:
            say("DRIFT deficit_closure %.4f"
                % abs(xs / math.log(10.0) - x8 / math.log(10.0)))
            say("  the published field alone puts it at 10^%.4f, so "
                "adding seven octaves moved it %.4f"
                % (x8 / math.log(10.0),
                   abs(xs / math.log(10.0) - x8 / math.log(10.0))))
        say("SHAPES 1")
        say("  the top of the field is 10^%.4f, so this is %.1f "
            "decades outside it"
            % (x.max() / math.log(10.0),
               (xs - x.max()) / math.log(10.0)))
        say("  no closure is published from this. {#rem:shapepower} "
            "forbids reading a level")
        say("  off an underived shape extrapolated past its data, "
            "and a quadratic in")
        say("  log-log is exactly that. What is measured is the "
            "drift, not the crossing.")
    say("  the sixteen-octave route could only exclude drifts above "
        "%.6f" % (2.0 * oct_se))

    say()
    say("=" * 70)
    say("A1 %s  A2 %s  A3 %s  A4 %s  A5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (a1, a2, a3, a4, a5)))

    head = [
        "SEED: the null draws from numpy default_rng at "
        "seed %d; without it the file does not reproduce "
        "its own null." % SEED,
        "STATISTIC: the coefficient beta of the quadratic term in",
        "           log(|sum a|/l2) = c + alpha x + (1/2) beta x^2",
        "           over every N in the field, x = log N -- the rate",
        "           at which the deficit drifts -- against the same",
        "           quantity estimated by fitting a line to sixteen",
        "           octave exponents in results/audit_alpha_reach.txt.",
        "NULL: none is run and none applies. |sum a| and its l2 norm",
        "      are deterministic once N is fixed; there is no",
        "      sampling noise and no background to detect against.",
        "      The coin arms for the sign axis were run in",
        "      lab_primorial_share.py.",
        "FIELD: N = 2^a 5^b with both a >= 1 and b >= 1 in",
        "       [2e5, 8e9], one coprimality class -- the class {2,5},",
        "       k coprime to 10 and N even; k squarefree with",
        "       2 <= k < N^0.56; m over 1 <= m < N/k with (m,k) = 1.",
        "       The points below 1.024e8 are recomputed here with the",
        "       kind-byte sieve and block sums of",
        "       code/audit_alpha_reach.py; those above it are read",
        "       from its POINT markers, and the twelve that are both",
        "       are A1's control. The global deficit A2 checks",
        "       against comes from results/audit_alpha_local.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not (a1 and a2):
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
