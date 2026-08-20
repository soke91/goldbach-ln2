# -*- coding: utf-8 -*-
r"""
The same degree question, asked of the ladder.

WHAT IS AT STAKE

{#rem:deficitregion} found that on the sign axis every polynomial
degree the data is offered keeps resolving -- t of -11.67, -2.96,
+2.61, -1.34, -1.16, +2.67, +7.41 for degrees two through eight --
while the residual falls only eight per cent across six added
parameters.  It read that as a family flexible enough to chase a
shape it does not contain, and concluded that no polynomial statement
about that quantity is stable past the field.

**The level axis rests on a quadratic.**  {#rem:laddercurve} resolved
its coefficient at t = 3.29 on fifteen rungs and {#rem:laddercap} at
t = 4.58 on seventeen, and {#rem:rung16}, {#rem:rung17} and
{#rem:rung18} each read a rung as the curvature predicting out of
sample.  If the ladder shows the same escalation, that quadratic is
not a shape the data chose; it is the first term anybody tested, and
the six out-of-sample hits are hits of a term rather than of a form.

The two axes are not the same measurement and the answer need not be
the same.  The ladder is nineteen points where the deficit had a
hundred and fifty-six, and fewer points resolve fewer degrees for
reasons that have nothing to do with the shape -- so a ladder that
stops escalating may be stopping for want of data rather than for
having found its form.  That is the reading D3 exists to keep
available.

Nothing is measured.  The nineteen uniform rungs are read from
results/audit_ladder_cap.txt for rungs 0 to 16 and from
results/audit_primorial_rung17.txt and rung18.txt for the top two.

BACKS: Remark {#rem:ladderdegree} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  D1  The control.  A quadratic on the seventeen rungs
      results/audit_ladder_cap.txt fitted returns the t that file
      declares, 4.58, to the two decimals printed.
  D2  On all nineteen rungs a cubic term is unresolved, |t| < 2, so
      the quadratic is not obviously the first of a sequence.
  D3  And the escalation does not run: across degrees two to six on
      the nineteen rungs, at most one coefficient beyond the second
      resolves.
  D4  The residual concedes what the coefficients claim: from the
      quadratic to degree six the r.m.s. residual falls by more than
      the eight per cent the sign axis conceded over the same span
      of added parameters.

REFUTATION RULE (fixed before the run)

  D1  REFUTED outside the printed decimals; then this is not the fit
      those remarks made.  THIS ONE GATES.
  D2  REFUTED if the cubic resolves, and that is the outcome that
      would cost the most: the ladder's quadratic would be in the
      position the sign axis's is, its six out-of-sample hits would
      be a term predicting rather than a form, and every reading
      that leans on the curvature -- {#rem:laddercurve},
      {#rem:curvebound}, the three rung remarks -- would have to say
      so.  If instead the cubic is **unresolved**, that is not proof
      the quadratic is right: nineteen points resolve fewer degrees
      than a hundred and fifty-six whatever the shape, and D3 and D4
      are what separate "the form was found" from "the data ran
      out".
  D3  REFUTED if two or more coefficients past the second resolve.
      The escalation would then be running on this axis too.
  D4  REFUTED if the residual falls by eight per cent or less.  Then
      the ladder buys as little per parameter as the deficit did,
      and D2's silence is the data running out rather than the shape
      being right -- **which is the reading that makes an unresolved
      cubic worth nothing**.

  D1 gates.  D2 to D4 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  Nothing is sampled: fits of
  different degree are compared on the same deterministic exponents,
  and the question is which coefficients the data resolves.  The coin
  arms for the ladder were run in lab_primorial_ladder.py and
  lab_primorial_share.py.
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
OUT = os.path.join(RES, "audit_ladder_degree.txt")

BASE = 30030
DEC = 2
SIGNAXIS_GAIN = 0.08                # what the deficit's residual bought


def read_ladder():
    """the nineteen uniform rungs at cap 10^6, from three files"""
    src = io.open(os.path.join(RES, "audit_ladder_cap.txt"),
                  encoding="utf-8").read()
    ns, ex = [], []
    for ln in src.splitlines():
        m = re.match(r"^  (\d+)\s+(\d+)\s+(.*)$", ln)
        if not m:
            continue
        f = [t for t in m.group(3).split() if t != "cap-invariant"]
        if len(f) < 8 or f[6] == "none":
            continue
        ns.append(int(m.group(2)))
        ex.append(float(f[7]))
    n17 = len(ns)
    for j, fn in ((17, "audit_primorial_rung17.txt"),
                  (18, "audit_primorial_rung18.txt")):
        s = io.open(os.path.join(RES, fn), encoding="utf-8").read()
        N = BASE * (1 << j)
        m = re.search(r"^  N = " + str(N) +
                      r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                      r"K\*_R \d+\s+exp ([\d.]+)\s*$", s, re.M)
        ns.append(N)
        ex.append(float(m.group(1)))
    return ns, ex, n17


def read_pub_t():
    src = io.open(os.path.join(RES, "audit_ladder_cap.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^TSTAT ladder_cap_quadratic ([-+]?[\d.]+)\s*$",
                  src, re.M)
    if not m:
        raise SystemExit("no TSTAT ladder_cap_quadratic marker")
    return float(m.group(1))


def fit(x, y, deg):
    A = np.column_stack([x ** d / math.factorial(d)
                         for d in range(deg + 1)])
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

    ns, ex, n17 = read_ladder()
    pubt = read_pub_t()
    say("read %d uniform rungs: %d from results/audit_ladder_cap.txt "
        "and the top two" % (len(ns), n17))
    say("  from results/audit_primorial_rung17.txt and rung18.txt; "
        "nothing is measured here")
    say("READ audit_ladder_cap.txt TSTAT ladder_cap_quadratic %.2f"
        % pubt)
    say("PRINTBOUND audit_ladder_degree %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    x = np.log(np.array(ns, dtype=np.float64))
    y = np.array(ex)
    say("  the ladder spans %.4f in log N" % (x.max() - x.min()))
    say("SCALES 1")

    # -------------------------------------------------------------- D1
    say()
    say("D1  the control on the seventeen {#rem:laddercap} fitted")
    c17, v17, _r17 = fit(x[:n17], y[:n17], 2)
    t17 = float(c17[2]) / math.sqrt(float(v17[2, 2]))
    d1 = abs(t17 - pubt) <= 0.5 * 10.0 ** (-DEC)
    say("  quadratic t here %.2f against the declared %.2f"
        % (t17, pubt))
    say("  D1 %s   (cap: the printing bound)"
        % ("hold" if d1 else "REFUTED"))
    if not d1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- D2
    say()
    say("D2  does a cubic term resolve on all nineteen?")
    c3, v3, r3 = fit(x, y, 3)
    g = float(c3[3])
    sg = math.sqrt(float(v3[3, 3]))
    tg = g / sg
    d2 = abs(tg) < 2.0
    say("  cubic coefficient %+.10f +- %.10f, t = %.2f" % (g, sg, tg))
    say("TSTAT ladder_degree_cubic %.2f" % tg)
    if abs(tg) < 2.0:
        say("UNRESOLVED SIGN ladder_degree_cubic")
    say("SPREAD ladder_degree_cubic %.4f" % (x.max() - x.min()))
    say("  D2 %s   (cap: |t| = 2)" % ("hold" if d2 else "REFUTED"))

    # -------------------------------------------------------------- D3
    say()
    say("D3  does the escalation run here?")
    say("  degree   top coefficient t   r.m.s. residual   resolves?")
    res = []
    for deg in range(2, 7):
        cc, vv, rr = fit(x, y, deg)
        tt = float(cc[deg]) / math.sqrt(float(vv[deg, deg]))
        res.append((deg, tt, rr))
        say("  %-8d %+-19.2f %-17.6f %s"
            % (deg, tt, rr, "yes" if abs(tt) > 2.0 else "no"))
    past = sum(1 for d, t, _ in res if d > 2 and abs(t) > 2.0)
    d3 = past <= 1
    say("DEGREES ladder_degree %d" % past)
    say("  %d coefficient(s) past the second resolve" % past)
    say("  D3 %s   (cap: at most one)"
        % ("hold" if d3 else "REFUTED"))

    # -------------------------------------------------------------- D4
    say()
    say("D4  does the residual concede what the coefficients claim?")
    r2 = res[0][2]
    r6 = res[-1][2]
    gain = (r2 - r6) / r2
    d4 = gain > SIGNAXIS_GAIN
    say("  the r.m.s. residual runs %.6f at degree two to %.6f at "
        "degree six," % (r2, r6))
    say("  a fall of %.4f against the %.2f the sign axis conceded "
        "over the same" % (gain, SIGNAXIS_GAIN))
    say("  span of added parameters")
    say("SCATTER slope_audit_ladder_degree %.6f" % r6)
    say("  D4 %s   (cap: the sign axis's gain)"
        % ("hold" if d4 else "REFUTED"))

    say()
    say("what this settles")
    if d2 and d3 and d4:
        say("  the ladder does not escalate and its residual pays "
            "for the parameters it")
        say("  takes, so the quadratic is a shape this axis chose "
            "rather than the first")
        say("  term anybody tried -- which is not the same as its "
            "being derived")
    elif not d2 or not d3:
        say("  the escalation runs here too: the ladder's quadratic "
            "is in the position")
        say("  the sign axis's was, and every reading that leans on "
            "the curvature has")
        say("  to say so")
    else:
        say("  the cubic is quiet but the residual concedes little, "
            "so the silence is the")
        say("  data running out rather than the shape being right")

    say()
    say("=" * 70)
    say("D1 %s  D2 %s  D3 %s  D4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (d1, d2, d3, d4)))

    head = [
        "STATISTIC: the highest coefficient of a polynomial of each",
        "           degree from two to six fitted to the nineteen",
        "           uniform rungs' level exponents against log N,",
        "           with its t and the r.m.s. residual, against the",
        "           same escalation measured on the sign axis in",
        "           results/audit_deficit_region.txt.",
        "NULL: none is run and none applies. Nothing is sampled:",
        "      fits of different degree are compared on the same",
        "      deterministic exponents and the question is which",
        "      coefficients the data resolves. The coin arms for the",
        "      ladder were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py.",
        "FIELD: N = 30030*2^j for j = 0..18, the odd radical",
        "       3*5*7*11*13 fixed so the threshold is constant; the",
        "       level exponent log K*_R / log N at the uniform cap",
        "       1000000. Nothing is measured here: rungs 0 to 16 are",
        "       read from results/audit_ladder_cap.txt and the top",
        "       two from results/audit_primorial_rung17.txt and",
        "       results/audit_primorial_rung18.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not d1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
