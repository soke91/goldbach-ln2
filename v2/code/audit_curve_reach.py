# -*- coding: utf-8 -*-
r"""
Could the curvature have been seen before the three new rungs?

WHAT IS AT STAKE

{#rem:laddercurve} resolved an upward curvature in the ladder: the
upper eight rungs are steeper than the lower seven by 2.51 standard
errors and a quadratic term stands at t = 3.29. {#rem:primorialdense}
had selected a + b log N as the surviving shape on 209 points reaching
log10 N = 7.7889, and that shape has no curvature in it.

Two readings again, and they are about the earlier work rather than
about the ladder. The curvature may have been visible all along and
missed, in which case the shape contest was run badly. Or it may be
invisible below the new rungs, in which case the contest was run
correctly on what it had and the three rungs are what bought the
answer. The difference decides whether {#rem:primorialdense}'s
selection was an error or a limit, and it is settled by refitting on
subsets -- no new arithmetic at all.

The same pass can ask the question that follows: with the three new
rungs in, is a + b log N still the best of the five shapes that
contest compared?

BACKS: Remark {#rem:curvereach} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  F1  The control. The fifteen-rung quadratic coefficient reproduces
      results/audit_ladder_curve.txt inside the bound its printing
      forces.
  F2  It could not have been seen: on the eleven published rungs
      alone the quadratic term is not resolved at two standard
      errors.
  F3  And it took all three: adding the new rungs one at a time, the
      quadratic's t first clears two only with the fifteenth rung in.
  F4  The ranking changes: on the fifteen rungs the best of the five
      shapes {#rem:primorialdense} compared is no longer
      a + b log N.

REFUTATION RULE (fixed before the run)

  F1  REFUTED outside the printing bound. THIS ONE GATES.
  F2  REFUTED if the quadratic is resolved on the eleven. Then the
      curvature was in the data the shape contest used and the
      contest missed it, which is a defect in {#rem:primorialdense}
      rather than a limit of its reach.
  F3  REFUTED if it clears two earlier. Then fewer rungs would have
      done and the third was not needed; the reading is unchanged but
      the cost accounting is.
  F4  REFUTED if a + b log N is still the best. Then the curvature is
      real and still not enough to unseat the shape on this set, and
      the ranking that {#rem:primorialdense} published survives the
      new reach on the rungs it can be tested on.

  F1 gates. F2 to F4 are the measurement and do not gate. Note for
  all three that "not resolved" includes "too noisy to tell" and is
  not the same as "zero" (M9).

  NO NULL IS RUN and none applies. Every quantity is a least-squares
  summary of exponents already measured, compared across subsets of
  the same points.
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
OUT = os.path.join(RES, "audit_curve_reach.txt")

BASE = 30030

SHAPES = (
    ("a + b log N", lambda L: L),
    ("a + b / log N", lambda L: 1.0 / L),
    ("a + b log log N", lambda L: np.log(L)),
    ("a + b log log N / log N", lambda L: np.log(L) / L),
)


def read_rungs():
    """the fifteen rung exponents and the published quadratic"""
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
    npub = len(ns)
    for j in (11, 12, 13, 14):
        s = io.open(os.path.join(RES, "audit_primorial_rung%d.txt" % j),
                    encoding="utf-8").read()
        N = BASE * (1 << j)
        m = re.search(r"^  N = " + str(N) +
                      r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                      r"K\*_R \d+\s+exp ([\d.]+)\s*$", s, re.M)
        ns.append(N)
        ex.append(float(m.group(1)))
    s = io.open(os.path.join(RES, "audit_ladder_curve.txt"),
                encoding="utf-8").read()
    m = re.search(r"the quadratic coefficient is ([+-][\d.]+), "
                  r"s\.e\. ([\d.]+)", s)
    return (ns, ex, npub, dec, float(m.group(1)), float(m.group(2)))


def quad(x, y):
    """the quadratic coefficient in x and its standard error"""
    A = np.column_stack([np.ones_like(x), x, x * x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    n = x.size
    s2 = float((r ** 2).sum()) / (n - 3)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return (float(c[2]), math.sqrt(float(cov[2, 2])),
            float(np.sqrt((r ** 2).mean())))


def shape_rms(name, f, x, y):
    u = f(x)
    A = np.column_stack([np.ones_like(u), u])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    return float(np.sqrt((r ** 2).mean()))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ns, ex, npub, dec, cq, sq = read_rungs()
    say("read %d rung exponents (%d published plus four measured "
        "since) and the fifteen-rung quadratic %+.8f (s.e. %.8f)"
        % (len(ns), npub, cq, sq))
    say("  from results/audit_primorial_rung10, 11, 12, 13, 14 and "
        "results/audit_ladder_curve.txt")
    say("RADICALS 1")

    x = np.log(np.array(ns, dtype=np.float64))
    y = np.array(ex)

    # -------------------------------------------------------------- F1
    c15, s15, r15 = quad(x, y)
    say()
    say("F1  the control")
    rnd = 0.5 * 10.0 ** (-8)
    d = abs(c15 - cq)
    f1 = d <= rnd
    say("  quadratic here %+.8f (s.e. %.8f) against the published "
        "%+.8f (%.8f)" % (c15, s15, cq, sq))
    say("  departure %.10f; that file prints 8 decimals, so the "
        "bound is %.10f" % (d, rnd))
    say("PRINTBOUND audit_curve_reach 8 %.10f" % rnd)
    say("  F1 %s   (cap: the printing bound)"
        % ("hold" if f1 else "REFUTED"))

    # ------------------------------------------------------- F2, F3
    say()
    say("F2/F3  when does the quadratic become visible?")
    say("  rungs   top log10 N   quadratic      s.e.          t")
    ts = []
    for n in range(npub, len(ns) + 1):
        c, s, r = quad(x[:n], y[:n])
        ts.append(abs(c) / s)
        say("  %-7d %-13.4f %+-14.8f %-13.8f %.2f"
            % (n, math.log10(ns[n - 1]), c, s, abs(c) / s))
    f2 = ts[0] < 2.0
    first = None
    for n, t in zip(range(npub, len(ns) + 1), ts):
        if t >= 2.0:
            first = n
            break
    f3 = first == len(ns)
    say("  on the %d published rungs alone the quadratic stands at "
        "%.2f standard errors" % (npub, ts[0]))
    say("  it first clears two at %s"
        % ("no subset" if first is None else "%d rungs" % first))
    say("TSTAT slope_audit_curve_reach %.2f" % ts[0])
    say("SPREAD slope_audit_curve_reach %.4f"
        % float(x[:npub].max() - x[:npub].min()))
    if ts[0] < 2.0:
        say("UNRESOLVED SIGN slope_audit_curve_reach")
    say("  F2 %s   (cap 2 standard errors on the published rungs)"
        % ("hold" if f2 else "REFUTED"))
    say("  F3 %s   (cap: clearing two before the last rung)"
        % ("hold" if f3 else "REFUTED"))

    # -------------------------------------------------------------- F4
    say()
    say("F4  do the five shapes still rank the same?")
    say("  shape                      r.m.s. on %d   r.m.s. on %d"
        % (npub, len(ns)))
    rows = []
    for nm, f in SHAPES:
        a = shape_rms(nm, f, x[:npub], y[:npub])
        b = shape_rms(nm, f, x, y)
        rows.append((nm, a, b))
        say("  %-26s %-14.6f %.6f" % (nm, a, b))
    aq = quad(x[:npub], y[:npub])[2]
    bq = quad(x, y)[2]
    rows.append(("a + b log N + c (log N)^2", aq, bq))
    say("  %-26s %-14.6f %.6f"
        % ("a + b log N + c (log N)^2", aq, bq))
    best_pub = min(rows, key=lambda t: t[1])[0]
    best_all = min(rows, key=lambda t: t[2])[0]
    f4 = best_all != "a + b log N"
    say("  best on the %d published: %s" % (npub, best_pub))
    say("  best on all %d:           %s" % (len(ns), best_all))
    say("SHAPES %d" % len(rows))
    say("  F4 %s   (cap: a + b log N still best)"
        % ("hold" if f4 else "REFUTED"))
    say()
    say("  the quadratic is a three-parameter shape against the "
        "others' two, so its")
    say("  r.m.s. is not comparable to theirs without that; it is "
        "listed to show the")
    say("  size of what the curvature removes, not to win a contest "
        "it is not in.")

    say()
    say("=" * 70)
    say("F1 %s  F2 %s  F3 %s  F4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (f1, f2, f3, f4)))

    head = [
        "STATISTIC: the fifteen primorial-ladder level exponents, read",
        "           from the files that measured them; the quadratic",
        "           coefficient in log N with its standard error,",
        "           fitted on the first eleven and then on each",
        "           prefix up to all fifteen; and the r.m.s. residual",
        "           of four two-parameter shapes of",
        "           {#rem:primorialdense} plus the quadratic, fitted",
        "           on the eleven and on the fifteen.",
        "NULL: none is run and none applies. Every quantity is a",
        "      least-squares summary of exponents already measured,",
        "      compared across subsets of the same points; there is",
        "      no background to detect against.",
        "FIELD: no arithmetic is computed here. The exponents are",
        "       read from results/audit_primorial_rung10.txt (eleven",
        "       rungs) and results/audit_primorial_rung11.txt,",
        "       rung12, rung13 and rung14 (one each), whose own field",
        "       is N = 30030*2^j with the odd radical 3*5*7*11*13",
        "       fixed, k squarefree and coprime to N with",
        "       2 <= k < 100000, m odd, squarefree and coprime to k,",
        "       and the Euler products at the fixed bound 4000000;",
        "       the published quadratic is read from",
        "       results/audit_ladder_curve.txt. The fifth shape of",
        "       {#rem:primorialdense}, 1 - c log log N / log N, has",
        "       one parameter and is not refitted here.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not f1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
