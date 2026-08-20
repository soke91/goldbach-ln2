# -*- coding: utf-8 -*-
r"""
What the level's residual is a function of, on a closed list of candidates

WHAT IS AT STAKE

rem:levelfine established, on thirty-nine adjacent pairs, that
|dL| does not follow the spacing: the scatter in L(N) =
log(|sum a| / l2) is not a scale effect.  Its size in that band is the
r.m.s. residual about the line, 0.007363.  What the scatter *is* was
left open.

Writing this run closed the list of things it could be.  Inside one
radical every construction-level quantity is fixed: the singular
series depends on N only through its radical, and so does which k the
range excludes, so neither varies here at all.  What varies from N to
N is only how log N is split among the three primes.  And that split
is not three free numbers either --

    v2 log 2 + v3 log 3 + v5 log 5 = log N

so once the log N trend is removed, **exactly two directions remain**.
Writing the shares s2 = v2 log 2 / log N and s5 = v5 log 5 / log N,
with s3 = 1 - s2 - s5 determined, a regression of the residual on
(s2, s5) tests the whole list.  **This is the first question in this
branch whose candidate set is complete rather than chosen.**

Nothing is measured here.  All forty L are read from the POINT markers
of results/audit_level_fine.txt.

BACKS: Remark {#rem:levelresidual} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  THE GATE.  A line on the forty read L returns rem:levelfine's
      published slope +0.251027 and r.m.s. residual 0.007363, to six
      decimals.
  V2  **The residual is the split.**  Regressing it on (s2, s5), at
      least one coefficient resolves at |t| above 3.
  V3  And the split explains a real share of it: the regression's
      R-squared exceeds 0.3.
  V4  And that is not what two regressors on forty points give by
      chance: the R-squared exceeds the ninety-fifth percentile of
      4000 permutations of the residuals against the same shares.

REFUTATION RULE (fixed before the run)

  V1  REFUTED outside six decimals on either; nothing below is
      reported, since a different line is a different residual.
  V2  **REFUTED at |t| of 3 or below on both.**  Then the residual is
      not a function of the exponent split, and since the split is
      the only thing that varies inside a radical at fixed N-trend,
      **the residual is not a function of N's arithmetic at all in
      any form this construction can express** -- it would be noise
      of the object rather than structure, and the 0.0074 would be
      irreducible for a reason, not merely by observation.  That is
      the outcome that costs the most and also the one that says the
      most; it must be stated in those words.
  V3  REFUTED at or below 0.3.
  V4  REFUTED at or below the permuted ninety-fifth percentile.  V3
      and V4 are separate on purpose: two regressors on forty points
      carry an expected R-squared near 0.05 under pure noise, so a
      value above 0.3 that fails V4 would mean the permutation null
      is wider than that expectation and V3's cap was the wrong
      threshold.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  The permutation
  null's own ninety-fifth percentile is printed, so V3's cap of 0.3
  can be checked against what chance actually gives here rather than
  against the textbook expectation.  **If the permuted percentile
  comes out above 0.3, V3 tests nothing** and its verdict stands
  without a reading, exactly as the underpowering clauses in
  rem:levelfine and rem:radicalblind did.  The seed is fixed and
  printed; without it this file does not reproduce its own null.

  WHAT THIS CANNOT DO.  One radical, one band, forty points.  The
  list of candidates is closed for *this construction* -- a quantity
  that depends on N through something the k-range and the singular
  series both ignore would not be on it, and no such quantity is
  known here.  Nothing in this run bounds |sum a| or moves item 5.
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
OUT = os.path.join(ROOT, "results", "audit_level_residual.txt")
SRC = os.path.join(ROOT, "results", "audit_level_fine.txt")

PUBSLOPE = 0.251027
PUBRMS = 0.007363
DEC = 6
TCAP = 3.0
R2CAP = 0.3
DRAWS = 4000
SEED = 20260921
PCTL = 95.0


def vp(n, p):
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def read_points():
    src = io.open(SRC, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^POINT fine_(\d+) ([-+]?[\d.]+)\s*$",
                         src, re.M):
        out[int(m.group(1))] = float(m.group(2))
    if not out:
        raise SystemExit("no fine markers in audit_level_fine.txt")
    return out


def lsq(cols, y):
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ coef
    n, p = A.shape
    dof = max(n - p, 1)
    s2 = float((r ** 2).sum()) / dof
    cov = s2 * np.linalg.pinv(A.T @ A)
    se = np.sqrt(np.abs(np.diag(cov)))
    ss = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((r ** 2).sum()) / ss if ss > 0 else 0.0
    return coef, se, r2, math.sqrt(float((r ** 2).mean()))


HEAD = [
    "STATISTIC: the residual of L(N) = log(|sum a| / l2) about a line",
    "           in log N, regressed on the shares of log N carried by",
    "           each prime of the radical, with a permutation null.",
    "FIELD: the forty N of radical {2,3,5} of rem:levelfine. Nothing",
    "       is measured here: every L is READ from the POINT markers",
    "       of results/audit_level_fine.txt.",
    "DERIVED: v2 log2 + v3 log3 + v5 log5 = log N, so after the log N",
    "         trend the split has exactly two free directions and",
    "         (s2, s5) exhausts the candidates this construction can",
    "         express.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pts = read_points()
    ns = sorted(pts)
    say("READ audit_level_fine.txt slope %.6f" % PUBSLOPE)
    say("READ audit_level_fine.txt rms %.6f" % PUBRMS)
    say("  read %d L markers" % len(ns))
    say("PRINTBOUND audit_level_residual %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("SEED: numpy default_rng at %d; the permutation null does not"
        % SEED)
    say("      reproduce without it")
    say("  |t| cap %.1f, R-squared cap %.2f, %d draws, percentile %.1f"
        % (TCAP, R2CAP, DRAWS, PCTL))
    say("RADICALS 1")

    x = np.array([math.log(n) for n in ns])
    y = np.array([pts[n] for n in ns])
    one = np.ones_like(x)

    # -------------------------------------------------------------- V1
    coef, _, _, rms = lsq([one, x], y)
    say()
    say("V1  the gate")
    a = abs(float(coef[1]) - PUBSLOPE) < 10.0 ** (-DEC)
    b = abs(rms - PUBRMS) < 10.0 ** (-DEC)
    v1 = a and b
    say("  slope here %+.6f against its %+.6f  %s"
        % (coef[1], PUBSLOPE, "ok" if a else "MISMATCH"))
    say("  r.m.s. here %.6f against its %.6f  %s"
        % (rms, PUBRMS, "ok" if b else "MISMATCH"))
    say("  V1 %s   (cap: %d decimals)"
        % ("hold" if v1 else "REFUTED", DEC))
    if not v1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    res = y - (coef[0] + coef[1] * x)
    s2 = np.array([vp(n, 2) * math.log(2) / math.log(n) for n in ns])
    s3 = np.array([vp(n, 3) * math.log(3) / math.log(n) for n in ns])
    s5 = np.array([vp(n, 5) * math.log(5) / math.log(n) for n in ns])
    say()
    say("  the shares close to one: max |s2+s3+s5-1| = %.2e"
        % float(np.abs(s2 + s3 + s5 - 1.0).max()))
    say("  s2 runs %.4f to %.4f, s5 runs %.4f to %.4f"
        % (s2.min(), s2.max(), s5.min(), s5.max()))
    say("    N          v2  v3  v5   s2      s5      residual")
    for n, a2, a5, r in zip(ns, s2, s5, res):
        say("  %-9d %-3d %-3d %-3d  %.4f  %.4f  %+.6f"
            % (n, vp(n, 2), vp(n, 3), vp(n, 5), a2, a5, r))
        say("POINT resid_%d %.6f" % (n, r))
    say("SCALES %d" % len(ns))

    # ---------------------------------------------------------- V2, V3
    coef2, se2, r2, _ = lsq([one, s2, s5], res)
    t2 = float(coef2[1]) / float(se2[1])
    t5 = float(coef2[2]) / float(se2[2])
    say()
    say("V2, V3  is the residual the split?")
    say("  on s2  %+.6f +- %.6f, t %+.2f" % (coef2[1], se2[1], t2))
    say("  on s5  %+.6f +- %.6f, t %+.2f" % (coef2[2], se2[2], t5))
    say("  R-squared %.6f" % r2)
    say("TSTAT residuals2 %.2f" % t2)
    say("SPREAD residuals2 %.6f" % float(se2[1]))
    say("TSTAT residuals5 %.2f" % t5)
    say("SPREAD residuals5 %.6f" % float(se2[2]))
    say("POINT residr2 %.6f" % r2)
    v2p = abs(t2) > TCAP or abs(t5) > TCAP
    v3 = r2 > R2CAP
    say("  V2 %s   (cap: |t| above %.1f on either)"
        % ("hold" if v2p else "REFUTED", TCAP))
    say("  V3 %s   (cap: above %.2f)"
        % ("hold" if v3 else "REFUTED", R2CAP))

    # -------------------------------------------------------------- V4
    say()
    say("V4  against a permutation null")
    rng = np.random.default_rng(SEED)
    nulls = np.empty(DRAWS)
    for i in range(DRAWS):
        _, _, rr, _ = lsq([one, s2, s5], rng.permutation(res))
        nulls[i] = rr
    p95 = float(np.percentile(nulls, PCTL))
    v4 = r2 > p95
    say("  %d permutations: median %.6f, %.0fth percentile %.6f"
        % (DRAWS, float(np.median(nulls)), PCTL, p95))
    say("  observed %.6f" % r2)
    say("POINT permp95 %.6f" % p95)
    say("  V4 %s   (cap: above the percentile)"
        % ("hold" if v4 else "REFUTED"))
    if p95 > R2CAP:
        say("  UNRESOLVED: the permuted percentile is above V3's cap "
            "of %.2f, so" % R2CAP)
        say("  V3 tested nothing and its verdict stands without a "
            "reading, as the")
        say("  rule says")

    say()
    say("=" * 70)
    say("V1 %s  V2 %s  V3 %s  V4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (v1, v2p, v3, v4)))
    say()
    if v2p and v4:
        say("the residual is a function of how log N is split among "
            "the primes of")
        say("the radical. that split is the whole of what varies "
            "inside a radical")
        say("at fixed N, so the scatter is arithmetic after all and "
            "this names it.")
    elif not v2p:
        say("the residual is not a function of the exponent split. "
            "the split is the")
        say("only thing that varies inside a radical once the N-trend "
            "is removed --")
        say("the singular series and the excluded k both depend on N "
            "through its")
        say("radical alone -- so the residual is not a function of "
            "N's arithmetic")
        say("in any form this construction can express. it is noise "
            "of the object,")
        say("and the %.6f rem:levelfine measured is irreducible for a "
            "reason" % PUBRMS)
        say("and not merely by observation.")
    else:
        say("a coefficient resolves but the permutation null covers "
            "the fit, so")
        say("what looks like structure is inside what chance gives on "
            "forty points")
        say("with two regressors, and nothing is claimed.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
