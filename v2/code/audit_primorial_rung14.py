# -*- coding: utf-8 -*-
r"""
The fifteenth rung: does the upward bend continue?

WHAT IS AT STAKE

{#rem:rung13} left the ladder bending upward. Its two out-of-sample
rungs both sit above the line the earlier rungs fit -- 1.27 prediction
standard errors at 30030*2^12 and 2.02 at 30030*2^13 -- the refit
scatter has begun to move (+0.0001, then +0.0005), and the second one
came out above BOTH candidate shapes rather than between them. Two
points are a direction, not a bend.

30030 * 2^14 = 492011520 at log10 N = 8.6919 is the third. It is the
test that separates the two readings: a noisy pair that comes back to
the line, or a curve the published shape does not have. Nothing else
in this repository can distinguish them, and the computation is
affordable where it was not before.

BACKS: Remark {#rem:rung14} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  The control. Rung 13 recomputed here reproduces
      results/audit_primorial_rung13.txt inside the bound its
      printing forces.
  C2  The margin keeps growing: the new exponent's margin over 1/2
      exceeds rung 13's 0.0283.
  C3  The bend continues: the new rung's departure from the
      fourteen-rung line exceeds the prediction standard error at its
      own abscissa, as rung 13's did.
  C4  And it grows: that departure, in units of its own prediction
      standard error, exceeds rung 13's 2.02.
  C5  The rung is above both candidate shapes again, as rung 13 was:
      the measured exponent exceeds both predictions there.

REFUTATION RULE (fixed before the run)

  C1  REFUTED outside the printing bound. THIS ONE GATES.
  C2  REFUTED if the margin does not grow. Four rungs have grown in a
      row; a fifth that does not would end the escalation.
  C3  REFUTED if the departure is inside the prediction standard
      error. That is the outcome that would settle it the other way:
      the third out-of-sample point returning to the line would make
      the two before it a noisy pair, and the published shape would
      stand. Note that "inside" here includes "too noisy to tell",
      and the ratio printed is what says which (M9).
  C4  REFUTED if the departure does not grow in those units. Then the
      bend is a fixed offset rather than a curve, which is a weaker
      and different statement from {#rem:rung13}'s.
  C5  REFUTED if the measurement falls below either shape. Then rung
      13's position above both was not a pattern, and the shapes are
      bracketing the ladder rather than trailing it.

  C1 gates. C2 to C5 are the measurement and do not gate.

  NO NULL IS RUN and none applies. A deterministic curve is located
  against a computed threshold. The coin arms for this statistic were
  run in lab_primorial_ladder.py and lab_primorial_share.py, and the
  scatter they left is the floor the margin is judged against.
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
OUT = os.path.join(RES, "audit_primorial_rung14.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R11 = module("audit_primorial_rung11")
primes_upto = R11.primes_upto
BASE = R11.BASE
CONTROL = BASE * (1 << 13)          # 246005760, the rung 13 point
NEW = BASE * (1 << 14)              # 492011520


def read_rung13():
    """rung 13's exponent, the fourteen-rung line and its ratio"""
    src = io.open(os.path.join(RES, "audit_primorial_rung13.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^  N = " + str(CONTROL) +
                  r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+K\*_R "
                  r"(\d+)\s+exp ([\d.]+)\s*$", src, re.M)
    dec = len(m.group(2).split(".")[1])
    m2 = re.search(r"the fourteen rungs refitted.*?slope "
                   r"([+-][\d.]+), r\.m\.s\. residual ([\d.]+)",
                   src, re.S)
    m3 = re.search(r"the new exponent is [\d.]+, margin ([\d.]+), "
                   r"against rung 12.s [\d.]+ and the scatter "
                   r"([\d.]+)", src)
    m4 = re.search(r"that is ([\d.]+) prediction standard errors",
                   src)
    return (int(m.group(1)), float(m.group(2)), dec,
            float(m2.group(1)), float(m2.group(2)),
            float(m3.group(1)), float(m3.group(2)),
            float(m4.group(1)))


def read_shapes():
    """the two best shapes, recovered from where each crosses"""
    src = io.open(os.path.join(RES, "audit_primorial_dense.txt"),
                  encoding="utf-8").read()
    rows = {}
    for m in re.finditer(r"^  (a \+ b log N|a \+ b log log N)\s+"
                         r"([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$",
                         src, re.M):
        rows[m.group(1)] = (float(m.group(2)), float(m.group(3)),
                            float(m.group(4)))
    m = re.search(r"top point ([\d.]+), parting at ([\d.]+)", src)
    return rows, float(m.group(1)), float(m.group(2))


def shape_at(name, half_at, target_at, x):
    """the two-parameter shape through (half_at, 0.5), (target_at, 0.56)

    x and the two crossings are log10 N. For a + b log N the inner log
    base is absorbed into b; for a + b log log N it is absorbed into a
    and b, so working in log10 throughout is the same family.
    """
    if name == "a + b log N":
        f = lambda t: t
    else:
        f = lambda t: math.log(t)
    u1, u2 = f(half_at), f(target_at)
    b = (0.56 - 0.5) / (u2 - u1)
    a = 0.5 - b * u1
    return a + b * f(x)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    (kpub, epub, dec, slope14, rms14pub, marg13, scatter,
     ratio13) = read_rung13()
    shapes, toppoint, parting = read_shapes()
    say("read rung 13 from results/audit_primorial_rung13.txt: "
        "K*_R %d, exponent %.4f, margin %.4f over one half"
        % (kpub, epub, marg13))
    say("  with the fourteen-rung line's slope %+.6f and r.m.s. "
        "residual %.4f, the ladder scatter %.4f, and rung 13's "
        "departure at %.2f prediction standard errors"
        % (slope14, rms14pub, scatter, ratio13))
    say("  and from results/audit_primorial_dense.txt the two best "
        "shapes, each recovered from where it crosses:")
    for nm, (r, h, t) in sorted(shapes.items()):
        say("    %-24s r.m.s. %.6f, 0.5 at %.4f, 0.56 at %.4f"
            % (nm, r, h, t))
    say("  the top measured point there is %.4f and the shapes part "
        "at %.4f" % (toppoint, parting))
    say("  the statistic, the sieve and the k-cap are imported from "
        "code/audit_primorial_rung11.py")

    ns, ex, sc0 = R11.read_ladder()
    say("  and %d published rungs from "
        "results/audit_primorial_rung10.txt" % len(ns))

    x0 = math.log10(NEW)
    say()
    say("the rungs: the control %d at log10 N = %.4f and the new %d "
        "at %.4f" % (CONTROL, math.log10(CONTROL), NEW, x0))
    say("  the new rung is %.4f decades past the parting point"
        % (x0 - parting))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in R11.factor_set(N) if q > 2))
                  for N in (CONTROL, NEW))))

    qs = [int(q) for q in primes_upto(R11.QSIEVE) if q > 2]
    say()
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NEW, ", ".join(map(str, qs))))
    lam, mu = R11.lambda_and_mu(NEW)
    sqf = mu != 0
    vmask = R11.residue_mask(NEW, qs)
    artin, twin = 1.0, 2.0
    assert CLIM == R11.CLIM
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    say("  the Euler products at the fixed bound %d: Artin %.9f, "
        "twin %.9f" % (CLIM, artin, twin))

    got = {}
    say()
    for N in (CONTROL, NEW):
        out = R11.measure(N, lam, mu, sqf, vmask, qs, artin, twin)
        if out is None:
            say("  N = %-12d no crossing below k = %d" % (N, R11.KCAP))
            continue
        kstar, e, bpn, beta, nk = out
        got[N] = (kstar, e)
        say("  N = %-12d thr %.6f  #k %-7d beta %.6f  K*_R %-8d "
            "exp %.4f" % (N, bpn, nk, beta, kstar, e))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))
    say("RADICALS 1")

    e13, e14 = got[CONTROL][1], got[NEW][1]

    # -------------------------------------------------------------- C1
    say()
    say("C1  the control at rung 13")
    rnd = 0.5 * 10.0 ** (-dec)
    d = abs(e13 - epub)
    c1 = d <= rnd
    say("  exponent here %.4f against the published %.4f, departure "
        "%.6f; the bound from %d decimals is %.8f"
        % (e13, epub, d, dec, rnd))
    say("  and K*_R here %d against %d" % (got[CONTROL][0], kpub))
    say("PRINTBOUND audit_primorial_rung14 %d %.8f" % (dec, rnd))
    say("  C1 %s   (cap: the printing bound)"
        % ("hold" if c1 else "REFUTED"))

    # -------------------------------------------------------------- C2
    say()
    say("C2  does the margin keep growing?")
    marg = e14 - 0.5
    c2 = marg > marg13
    say("  the new exponent is %.4f, margin %.4f, against rung 13's "
        "%.4f and the scatter %.4f" % (e14, marg, marg13, scatter))
    say("MARGIN audit_primorial_rung14 %.4f %.4f" % (marg, scatter))
    if marg <= scatter:
        say("INSIDE FLOOR audit_primorial_rung14")
    say("FLOOR primorial_rung14 %.4f" % scatter)
    say("  C2 %s   (cap: rung 13's margin)"
        % ("hold" if c2 else "REFUTED"))

    # -------------------------------------------------------------- C3
    say()
    say("C3  does the bend continue?")
    # the fourteen rungs are the eleven published and rungs 11-13
    x = np.log(np.array(ns + [BASE * (1 << 11), BASE * (1 << 12),
                              CONTROL], dtype=np.float64))
    src11 = io.open(os.path.join(RES, "audit_primorial_rung11.txt"),
                    encoding="utf-8").read()
    e11 = float(re.search(r"^  N = " + str(BASE * (1 << 11)) +
                          r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                          r"K\*_R \d+\s+exp ([\d.]+)\s*$", src11,
                          re.M).group(1))
    src12 = io.open(os.path.join(RES, "audit_primorial_rung12.txt"),
                    encoding="utf-8").read()
    e12 = float(re.search(r"^  N = " + str(BASE * (1 << 12)) +
                          r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                          r"K\*_R \d+\s+exp ([\d.]+)\s*$", src12,
                          re.M).group(1))
    y = np.array(ex + [e11, e12, epub])
    a, b = np.polyfit(x, y, 1)
    xn = math.log(NEW)
    pred = a * xn + b
    resid = e14 - pred
    n = x.size
    rss = float(((y - (a * x + b)) ** 2).sum())
    s2 = rss / (n - 2)
    sxx = float(((x - x.mean()) ** 2).sum())
    sepred = math.sqrt(s2 * (1.0 + 1.0 / n + (xn - x.mean()) ** 2
                             / sxx))
    c3 = abs(resid) > sepred
    ratio = abs(resid) / sepred
    say("  the fourteen-rung line refitted here: slope %+.6f against "
        "the published %+.6f" % (a, slope14))
    say("  it predicts %.4f; measured %.4f, departure %+.4f"
        % (pred, e14, resid))
    say("  in-sample r.m.s.        %.4f" % rms14pub)
    say("  prediction s.e. at x0   %.4f" % sepred)
    say("  that is %.2f prediction standard errors" % ratio)
    say("  C3 %s   (cap: the prediction standard error)"
        % ("hold" if c3 else "REFUTED"))
    say()
    say("C4  and does the departure grow?")
    c4 = ratio > ratio13
    say("  rung 13 stood at %.2f prediction standard errors, this "
        "rung at %.2f" % (ratio13, ratio))
    say("  C4 %s   (cap: rung 13's ratio)"
        % ("hold" if c4 else "REFUTED"))

    # --------------------------------------------------- C4 and C5
    say()
    say("C5  is the rung above both candidate shapes again?")
    nm1, nm2 = "a + b log N", "a + b log log N"
    r1, h1, t1 = shapes[nm1]
    r2, h2, t2 = shapes[nm2]
    p1 = shape_at(nm1, h1, t1, x0)
    p2 = shape_at(nm2, h2, t2, x0)
    say("  shape                    predicts here   measured minus it")
    say("  %-24s %.6f        %+.6f" % (nm1, p1, e14 - p1))
    say("  %-24s %.6f        %+.6f" % (nm2, p2, e14 - p2))
    say("  the two predictions differ by %.6f, against the surviving "
        "shape's r.m.s. %.6f" % (abs(p1 - p2), r1))
    c5 = e14 > p1 and e14 > p2
    say("  the measured %.4f is above both: %s"
        % (e14, "yes" if c5 else "NO"))
    say("SHAPES 2")
    say("  C5 %s   (cap: falling below either)"
        % ("hold" if c5 else "REFUTED"))

    say()
    say("  the fifteen rungs refitted, which is what any later "
        "extrapolation must use:")
    x14 = np.append(x, xn)
    y14 = np.append(y, e14)
    a14, b14 = np.polyfit(x14, y14, 1)
    r14 = y14 - (a14 * x14 + b14)
    rms14 = float(np.sqrt((r14 ** 2).mean()))
    se14 = math.sqrt(float((r14 ** 2).sum() / (x14.size - 2))
                     / float(((x14 - x14.mean()) ** 2).sum()))
    say("  slope %+.6f, r.m.s. residual %.4f, standard error %.6f, "
        "t = %.2f" % (a14, rms14, se14, abs(a14) / se14))
    say("SCATTER slope_audit_primorial_rung14 %.4f" % rms14)
    say("TSTAT slope_audit_primorial_rung14 %.2f" % (abs(a14) / se14))
    say("SPREAD slope_audit_primorial_rung14 %.4f"
        % float(x14.max() - x14.min()))
    if abs(a14) / se14 < 2.0:
        say("UNRESOLVED SIGN slope_audit_primorial_rung14")
    say("  and the r.m.s. against the fourteen-rung %.4f: %+.4f"
        % (rms14pub, rms14 - rms14pub))
    say("  no forecast is made from this.")

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  C3 %s  C4 %s  C5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (c1, c2, c3, c4, c5)))

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           at N = 30030*2^14 = 492011520 and, as a control,",
        "           at N = 30030*2^13; the margin over 1/2 against",
        "           rung 12's and against the ladder's scatter; the",
        "           departure from the thirteen-rung line against the",
        "           prediction standard error at the new abscissa;",
        "           the two best shapes of {#rem:primorialdense},",
        "           recovered from the crossings that file prints,",
        "           evaluated there and compared with the measured",
        "           exponent; and the fourteen-rung refit.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold; there is no",
        "      background to detect against. The coin arms for this",
        "      statistic were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py, and the scatter they left is",
        "      the floor the margin is judged against.",
        "FIELD: N = 246005760 and 492011520, the odd radical",
        "       3*5*7*11*13 fixed so the threshold is constant; k",
        "       squarefree and coprime to N with 2 <= k < 100000; m",
        "       odd, squarefree and coprime to k, m < N/k; the sieve",
        "       weight over the odd primes below 30; the Euler",
        "       products at the fixed bound 4000000. One odd radical,",
        "       as the RADICALS line declares. The statistic, the",
        "       sieve and the k-cap are imported from",
        "       code/audit_primorial_rung11.py; the published rungs",
        "       come from results/audit_primorial_rung10.txt, rung 11",
        "       from results/audit_primorial_rung11.txt, rung 12 and",
        "       the thirteen-rung line from",
        "       results/audit_primorial_rung12.txt, and the shapes",
        "       from results/audit_primorial_dense.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not c1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
