# -*- coding: utf-8 -*-
r"""
Offset or curve: what the fifteen rungs say about their own shape.

WHAT IS AT STAKE

{#rem:rung14} left one question open in the only place it can be
settled cheaply. Three rungs measured after the line was fitted all
sit above it, by 1.27, 2.02 and 1.29 prediction standard errors --
consistently above, not growing. Two readings fit that. The ladder may
be a straight line whose fit was pulled down by the small-N rungs, in
which case a line fitted on the upper rungs alone predicts the new
ones; or it may curve upward, in which case the local slope rises with
N and a quadratic term is resolved.

Nothing new has to be computed. Fifteen rung exponents are published
across four result files, and the two readings make opposite
predictions about them.

BACKS: Remark {#rem:laddercurve} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  The control. The fifteen rungs assembled here reproduce
      results/audit_primorial_rung14.txt's fifteen-rung slope and
      r.m.s. inside the bound its printing forces.
  E2  The ladder curves: the slope fitted on the upper rungs exceeds
      the slope on the lower rungs by more than two standard errors
      of the difference.
  E3  And a quadratic term in log N, fitted on all fifteen, is
      resolved positive.
  E4  So the offset reading fails: at least one of the three
      out-of-sample rungs departs from a line fitted on the upper
      rungs alone by more than its own prediction standard error.

REFUTATION RULE (fixed before the run)

  E1  REFUTED outside the printing bound. Then these are not the
      rungs those files measured. THIS ONE GATES.
  E2  REFUTED if the two slopes agree within two standard errors of
      their difference. Then the ladder is straight as far as its own
      rungs can tell and the three departures are a property of where
      the line was fitted, not of the ladder.
  E3  REFUTED if the quadratic term is not resolved positive. Same
      reading as E2 by a different route; both failing is the offset
      answer. Note that "not resolved" includes "too noisy to tell",
      which is not the same as "zero" (M9).
  E4  REFUTED if all three sit inside. Then a line fitted on the
      upper rungs predicts the out-of-sample points, the low rungs
      were what pulled the earlier fit down, and the published
      forecast has to be recomputed on the upper rungs rather than
      abandoned.

  E1 gates. E2 to E4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. Every quantity is a least-squares
  summary of exponents already measured, and the comparisons are
  between fits on subsets of the same points.
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
OUT = os.path.join(RES, "audit_ladder_curve.txt")

BASE = 30030


def read_rungs():
    """the fifteen rung exponents, from the four files that hold them"""
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
    for j, stem in ((11, "audit_primorial_rung11"),
                    (12, "audit_primorial_rung12"),
                    (13, "audit_primorial_rung13"),
                    (14, "audit_primorial_rung14")):
        s = io.open(os.path.join(RES, stem + ".txt"),
                    encoding="utf-8").read()
        N = BASE * (1 << j)
        m = re.search(r"^  N = " + str(N) +
                      r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                      r"K\*_R \d+\s+exp ([\d.]+)\s*$", s, re.M)
        ns.append(N)
        ex.append(float(m.group(1)))
    s = io.open(os.path.join(RES, "audit_primorial_rung14.txt"),
                encoding="utf-8").read()
    m = re.search(r"the fifteen rungs refitted.*?slope ([+-][\d.]+), "
                  r"r\.m\.s\. residual ([\d.]+)", s, re.S)
    return ns, ex, dec, float(m.group(1)), float(m.group(2))


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum())) \
        if n > 2 else float("inf")
    return float(a), float(b), float(np.sqrt((r ** 2).mean())), se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ns, ex, dec, slope15, rms15 = read_rungs()
    say("read %d rung exponents from results/audit_primorial_rung10, "
        "11, 12, 13 and 14" % len(ns))
    say("  with the published fifteen-rung slope %+.6f and r.m.s. "
        "residual %.4f" % (slope15, rms15))
    say("RADICALS 1")

    x = np.log(np.array(ns, dtype=np.float64))
    y = np.array(ex)
    say()
    say("  rung   N              log10 N   exponent")
    for i, (N, e) in enumerate(zip(ns, ex)):
        say("  %-6d %-14d %-9.4f %.4f" % (i, N, math.log10(N), e))

    # -------------------------------------------------------------- E1
    a, b, rms, se = fit(x, y)
    say()
    say("E1  the control")
    rnd = 0.5 * 10.0 ** (-dec)
    d = max(abs(a - slope15), abs(rms - rms15))
    e1 = d <= rnd
    say("  slope here %+.6f against the published %+.6f; r.m.s. %.4f "
        "against %.4f" % (a, slope15, rms, rms15))
    say("  worst departure %.6f; the tables print %d decimals, so the "
        "bound is %.8f" % (d, dec, rnd))
    say("PRINTBOUND audit_ladder_curve %d %.8f" % (dec, rnd))
    say("  E1 %s   (cap: the printing bound)"
        % ("hold" if e1 else "REFUTED"))

    # -------------------------------------------------------------- E2
    half = len(ns) // 2
    al, bl, rl, sl = fit(x[:half], y[:half])
    au, bu, ru, su = fit(x[half:], y[half:])
    dd = au - al
    sd = math.sqrt(sl * sl + su * su)
    e2 = dd > 2.0 * sd
    say()
    say("E2  do the two halves have the same slope?")
    say("  lower %d rungs: slope %+.6f, s.e. %.6f, r.m.s. %.4f"
        % (half, al, sl, rl))
    say("  upper %d rungs: slope %+.6f, s.e. %.6f, r.m.s. %.4f"
        % (len(ns) - half, au, su, ru))
    say("  the difference is %+.6f against %.6f, two standard errors "
        "of it being %.6f" % (dd, sd, 2.0 * sd))
    say("TSTAT slope_audit_ladder_curve %.2f" % (abs(dd) / sd))
    say("SPREAD slope_audit_ladder_curve %.4f"
        % float(x.max() - x.min()))
    if abs(dd) / sd < 2.0:
        say("UNRESOLVED SIGN slope_audit_ladder_curve")
    say("  E2 %s   (cap 2 standard errors of the difference)"
        % ("hold" if e2 else "REFUTED"))

    # -------------------------------------------------------------- E3
    say()
    say("E3  is a quadratic term resolved?")
    A = np.column_stack([np.ones_like(x), x, x * x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A.dot(c)
    n = x.size
    s2 = float((resid ** 2).sum()) / (n - 3)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    sq, ssq = float(c[2]), math.sqrt(float(cov[2, 2]))
    e3 = sq > 0.0 and abs(sq) / ssq >= 2.0
    say("  the quadratic coefficient is %+.8f, s.e. %.8f, t = %.2f"
        % (sq, ssq, abs(sq) / ssq))
    say("  the quadratic fit's r.m.s. is %.4f against the line's %.4f"
        % (float(np.sqrt((resid ** 2).mean())), rms))
    say("TSTAT slope_laddercurve_quad %.2f" % (abs(sq) / ssq))
    say("SPREAD slope_laddercurve_quad %.4f"
        % float(x.max() - x.min()))
    if abs(sq) / ssq < 2.0:
        say("UNRESOLVED SIGN slope_laddercurve_quad")
    say("  E3 %s   (cap 2 standard errors)"
        % ("hold" if e3 else "REFUTED"))

    # -------------------------------------------------------------- E4
    say()
    say("E4  does a line on the upper rungs predict the new ones?")
    say("  fitted on the rungs below the three out-of-sample points, "
        "and only those above the halfway mark:")
    keep = [i for i in range(len(ns)) if half <= i < len(ns) - 3]
    xf, yf = x[keep], y[keep]
    af, bf, rf, sf = fit(xf, yf)
    nf = xf.size
    rssf = float(((yf - (af * xf + bf)) ** 2).sum())
    s2f = rssf / (nf - 2)
    sxxf = float(((xf - xf.mean()) ** 2).sum())
    say("  %d rungs, slope %+.6f, s.e. %.6f, r.m.s. %.4f"
        % (nf, af, sf, rf))
    say("  rung   log10 N   predicted   measured   departure   "
        "pred s.e.   ratio")
    outside = 0
    for i in range(len(ns) - 3, len(ns)):
        pr = af * x[i] + bf
        rr = y[i] - pr
        sp = math.sqrt(s2f * (1.0 + 1.0 / nf
                              + (x[i] - xf.mean()) ** 2 / sxxf))
        if abs(rr) > sp:
            outside += 1
        say("  %-6d %-9.4f %-11.4f %-10.4f %+-11.4f %-11.4f %.2f"
            % (i, math.log10(ns[i]), pr, y[i], rr, sp, abs(rr) / sp))
    e4 = outside > 0
    say("  outside their own prediction standard error: %d of 3"
        % outside)
    say("  E4 %s   (cap: all three inside)"
        % ("hold" if e4 else "REFUTED"))

    say()
    say("=" * 70)
    say("E1 %s  E2 %s  E3 %s  E4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (e1, e2, e3, e4)))

    head = [
        "STATISTIC: the fifteen primorial-ladder level exponents",
        "           log K*_R / log N at N = 30030*2^j, read from the",
        "           result files that measured them; the",
        "           least-squares slope in log N on the whole set and",
        "           on each half; a quadratic term in log N fitted on",
        "           all fifteen with its own standard error; and the",
        "           three out-of-sample rungs against a line fitted",
        "           on the upper rungs below them, each with the",
        "           prediction standard error at its own abscissa.",
        "NULL: none is run and none applies. Every quantity is a",
        "      least-squares summary of exponents already measured,",
        "      and the comparisons are between fits on subsets of the",
        "      same points; there is no background to detect against.",
        "FIELD: no arithmetic is computed here. The exponents are",
        "       read from results/audit_primorial_rung10.txt (eleven",
        "       rungs) and results/audit_primorial_rung11.txt,",
        "       rung12, rung13 and rung14 (one each), whose own field",
        "       is N = 30030*2^j with the odd radical 3*5*7*11*13",
        "       fixed, k squarefree and coprime to N with",
        "       2 <= k < 100000, m odd, squarefree and coprime to k,",
        "       and the Euler products at the fixed bound 4000000.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not e1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
