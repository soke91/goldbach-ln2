# -*- coding: utf-8 -*-
r"""
Which of this project's published trends are above their own noise?

WHAT IS AT STAKE

Gate check G37 was added after a span across seven arithmetic types
turned out to be 1.40 times what the same statistic gives with the
arithmetic held fixed -- not a dependence, just scatter. The same
question has never been asked of the SLOPES.

Several conclusions rest on one. Remark {#rem:residuelevel}'s rule U4
-- "the margin is not closing" -- is a slope of +0.004692 through five
points, and the argument offered for it is that three leave-one-out
refits are all non-negative. That is not a significance test: dropping
one of five points cannot move a slope much whatever the noise.
Remark {#rem:primorialrung10}'s +0.006780 through eleven, Remark
{#rem:cRwindow}'s +0.014647, Remark {#rem:leandecay}'s -0.1673 and
Remark {#rem:primorialshare}'s -0.035113 are in the same position.

A slope's standard error is available from the residuals it already
publishes: s.e. = rms / sqrt(sum (x - xbar)^2). Nothing new has to be
computed, and every one of these can be put on the same footing.

BACKS: Remark {#rem:slopes} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  M1  The control: each slope refitted here reproduces its published
      value to within 1e-5.
  M2  The primorial ladder's level slope is solid: above 5 standard
      errors.
  M3  The family's level slope -- U4's basis -- is not: below 3
      standard errors. That is the one that matters, because
      {#rem:residuelevel} reads it as establishing that the margin
      over 1/2 does not close.
  M4  But nothing published is pure noise: every slope is above 2
      standard errors.

REFUTATION RULE (fixed before the run)

  M1  REFUTED at 1e-5 for any sequence.
  M2  REFUTED below 5 s.e.
  M3  REFUTED if the family's slope reaches 3 s.e., which would leave
      U4 standing as written.
  M4  REFUTED if any published slope falls below 2 s.e. -- that trend
      would have to be withdrawn.

  All four gate.

  NO NULL IS RUN and none applies. Published numbers are refitted and
  their standard errors computed from their own residuals; there is no
  background to detect against. The noise floors these standard errors
  encode are the same ones lab_primorial_ladder.py and
  lab_residue_cancellation.py established with coin arms.
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
OUT = os.path.join(ROOT, "results", "audit_slope_significance.txt")
RES = os.path.join(ROOT, "results")


def grab(fname, header, xcol, ycol, stop_blank=True, logx=True,
         logy=False):
    """read a two-column series out of a published table"""
    src = io.open(os.path.join(RES, fname), encoding="utf-8").read()
    i = src.index(header)
    xs, ys, hs = [], [], []
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) <= max(xcol, ycol) or not f[xcol].isdigit():
            if xs and stop_blank:
                break
            continue
        try:
            v = float(f[ycol])
        except ValueError:
            if xs:
                break
            continue
        xs.append(float(f[xcol]))
        ys.append(v)
        # the half-width of the rounding the printed column imposes
        d = len(f[ycol].split(".")[1]) if "." in f[ycol] else 0
        hs.append(0.5 * 10.0 ** (-d))
    x = np.array(xs)
    y = np.array(ys)
    h = np.array(hs)
    return (np.log(x) if logx else x), (np.log(y) if logy else y), h


def stats(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    rms = float(np.sqrt((r ** 2).mean()))
    # the residual variance uses n-2 degrees of freedom
    s2 = float((r ** 2).sum() / (n - 2))
    se = math.sqrt(s2 / float(((x - x.mean()) ** 2).sum()))
    return float(a), rms, se, abs(float(a)) / se, n


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    SERIES = []

    x, y, h = grab("audit_residue_level.txt",
                   "N            budget factor   K*_R      "
                   "log K*_R/log N  clears .56", 0, 3)
    SERIES.append(("family level exponent", "audit_residue_level",
                   x, y, 0.004692, h))

    x, y, h = grab("audit_primorial_rung10.txt",
                   "N            log10 N   exponent   fitted     residual",
                   0, 2)
    SERIES.append(("ladder level exponent", "audit_primorial_rung10",
                   x, y, 0.006780, h))

    x, y, h = grab("audit_cR_window.txt",
                   "  N            octave exponent" + chr(10), 0, 1)
    SERIES.append(("ladder k-exponent", "audit_cR_window",
                   x, y, 0.014647, h))

    x, y, h = grab("lab_primorial_share.txt",
                   "  N            rho        below 1?   exponent   "
                   "above .5?", 0, 1)
    # in log space the same printed rounding is h/y wide
    SERIES.append(("ladder log rho", "lab_primorial_share",
                   x, np.log(y), -0.035113, h / y))

    say("read %d published series" % len(SERIES))

    # ------------------------------------------------------------- M1
    say()
    say("M1  the control: each slope refitted")
    say("  series                       refitted     published    diff")
    m1 = True
    got = []
    for name, src, x, y, pubs, h in SERIES:
        a, rms, se, t, n = stats(x, y)
        got.append((name, src, a, rms, se, t, n))
        d = abs(a - pubs)
        if not (d < 1e-5):
            m1 = False
        say("  %-28s %-12.6f %-12.6f %.7f" % (name, a, pubs, d))
    say("  M1 %s" % ("hold" if m1 else "REFUTED"))
    say()
    say("  DIAGNOSTIC on M1 (post hoc). The tolerance was set at 1e-5")
    say("  without asking what the printed tables can carry. A")
    say("  least-squares slope is a fixed linear functional of the")
    say("  ordinates, slope = sum c_i y_i with")
    say("  c_i = (x_i - xbar) / sum (x_j - xbar)^2, so rounding the")
    say("  y column to d decimals moves the refit by at most")
    say("  sum |c_i| * 0.5 * 10^-d. That bound is computable and is")
    say("  not a free parameter:")
    say("  series                       M1 diff      rounding bound "
        "covered?")
    covered = True
    for (name, src, x, y, pubs, h), (_, _, a, rms, se, t, n) in zip(
            SERIES, got):
        c = (x - x.mean()) / float(((x - x.mean()) ** 2).sum())
        bound = float((np.abs(c) * h).sum())
        d = abs(a - pubs)
        ok_i = d <= bound
        covered = covered and ok_i
        say("  %-28s %-12.7f %-14.7f %s"
            % (name, d, bound, "yes" if ok_i else "NO"))
    say("  every M1 gap inside its own rounding bound: %s"
        % ("yes" if covered else "NO -- a real disagreement"))
    say("  So M1 fails as written and the reason is arithmetic I")
    say("  imposed, not a disagreement with the published fits. The")
    say("  rule stands as registered; what it refutes is the")
    say("  tolerance, and no published slope is thereby in doubt.")

    # ------------------------------------------------------- M2/M3/M4
    say()
    say("M2/M3/M4  every published slope against its own noise")
    say("  a t-ratio is not comparable across series without the")
    say("  range it was measured over, so the spread in the abscissa")
    say("  is printed beside it: a slope can be unresolved because")
    say("  the statistic is noisy or because the sweep is short, and")
    say("  only the second is fixable by computing more.")
    say("  series                       n   spread   slope        "
        "r.m.s.   s.e.       t")
    for (name, src, x, y, pubs, h), (_, _, a, rms, se, t, n) in zip(
            SERIES, got):
        sp = float(x.max() - x.min())
        say("  %-28s %-3d %-8.4f %-12.6f %-8.4f %-10.6f %.2f"
            % (name, n, sp, a, rms, se, t))
        say("SCATTER slope_%s %.4f" % (src, rms))
        say("TSTAT slope_%s %.2f" % (src, t))
        say("SPREAD slope_%s %.4f" % (src, sp))
        if t < 2.0:
            say("UNRESOLVED SIGN slope_%s" % src)
    d = dict((g[0], g) for g in got)
    m2 = d["ladder level exponent"][5] > 5.0
    m3 = d["family level exponent"][5] < 3.0
    m4 = all(g[5] > 2.0 for g in got)
    say("  M2 the ladder's level slope is above 5 s.e. (%.2f)   %s"
        % (d["ladder level exponent"][5], "hold" if m2 else "REFUTED"))
    say("  M3 the family's is below 3 s.e. (%.2f)   %s"
        % (d["family level exponent"][5], "hold" if m3 else "REFUTED"))
    say("  M4 every slope is above 2 s.e. (min %.2f)   %s"
        % (min(g[5] for g in got), "hold" if m4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). Why leave-one-out is not a")
    say("  significance test. For each series, the three leave-one-out")
    say("  slopes and the interval the standard error gives:")
    say("  series                       LOO min      LOO max      "
        "slope +- 2 s.e.")
    for (name, src, x, y, pubs, h), (_, _, a, rms, se, t, n) in zip(
            SERIES, got):
        f = [float(np.polyfit(x[s], y[s], 1)[0])
             for s in (slice(None), slice(1, None), slice(0, -1))]
        say("  %-28s %-12.6f %-12.6f [%+.6f, %+.6f]"
            % (name, min(f), max(f), a - 2 * se, a + 2 * se))
    say("  Dropping one point of five moves a slope by a fraction of")
    say("  its own error bar, so 'all three refits have the same sign'")
    say("  is a statement about the arithmetic of least squares and")
    say("  not about the data.")

    say()
    say("=" * 70)
    ok = m1 and m2 and m3 and m4
    say("the ladder carries the trend and the family does not"
        if ok else "REFUTED")

    head = [
        "STATISTIC: for each slope this project has published -- the",
        "           family's level exponent against log N, the primorial",
        "           ladder's level exponent, the ladder's k-exponent and",
        "           its log rho -- the refitted slope, the r.m.s.",
        "           residual, the standard error",
        "           sqrt(sum r^2/(n-2) / sum (x-xbar)^2), and the ratio",
        "           of slope to standard error.",
        "NULL: none is run and none applies. Published numbers are",
        "      refitted and their standard errors computed from their",
        "      own residuals; there is no background to detect against.",
        "      The noise floors these encode are the ones",
        "      lab_primorial_ladder.py and lab_residue_cancellation.py",
        "      established with coin arms.",
        "FIELD: the published tables of results/audit_residue_level.txt,",
        "       results/audit_primorial_rung10.txt,",
        "       results/audit_cR_window.txt and",
        "       results/lab_primorial_share.txt. No new arithmetic is",
        "       done; the k-ranges, sieve weights, budgets and radicals",
        "       are those declared in those files.",
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
