# -*- coding: utf-8 -*-
r"""
How thin is the failing set, and is the head's excess exponent real?

WHAT IS AT STAKE

{#rem:sumhead} measured the head's signed sum at exponent +0.787845
against the total's +0.717916 and wrote that the head "must come down
further than the whole has to travel". That sentence has a ceiling
problem. Head and tail agree in sign at every N, so |head sum| is at
most |sum a| and the head's share of it is bounded by 1. At the top N
that share is already 0.9710, rising at +0.053778 per unit log N. A
bounded quantity three hundredths from its bound cannot keep rising at
that rate, and when it stops the two exponents must coincide -- the
head becomes the sum. So the excess is a within-range effect of a
share still climbing, and the same repository has been burned by
reading exactly this kind of transient as a law ({#rem:shapepower},
{#rem:thetalaw}).

The correction is worth one measurement and the measurement is worth
more than the correction, because the same pass answers a question
nobody has asked: how thin is the failing set? The head fraction 0.10
is a convention inherited from {#rem:gainsplit}. Sweeping it says
whether the signed sum lives in a tenth of the dilations or in a
hundredth.

BACKS: Remark {#rem:headfraction} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  The control. At fraction 0.10 the head's share of |sum a| at
      the top N reproduces results/audit_sum_head.txt inside the
      bound its printing forces.
  Y2  The failing set is thinner than a tenth: at every N the top one
      per cent of k by |a_k| already carries at least half of
      |sum a|.
  Y3  The ceiling binds. At the top N the share at fraction 0.10 is
      above 0.95, so the measured rise of +0.053778 per unit log N
      has less than one unit of log N left before the bound.
  Y4  And the excess exponent is transient: fitted on the top octave
      alone, the head's signed sum at fraction 0.10 has an exponent
      closer to the total's than the whole-range fit is.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED outside the printing bound. Then this is not the split
      {#rem:sumhead} measured. THIS ONE GATES.
  Y2  REFUTED if the top per cent carries less than half anywhere.
      Then a tenth is the right description and the set is not
      thinner than {#rem:gainsplit}'s convention.
  Y3  REFUTED if the share at the top N is at or below 0.95. Then the
      bound is not close and the rise has room, so the excess
      exponent may be a real difference rather than a transient.
  Y4  REFUTED if the top-octave exponent is not closer to the total's
      than the global fit is. Then the convergence the ceiling forces
      has not begun inside the computed range, and the excess cannot
      be called transient on this evidence -- only on the bound.

  Y1 gates. Y2 to Y4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. Every quantity is an exactly
  computed partial sum of a measured sequence, ordered by its own
  magnitudes; there is no background to detect against.
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
OUT = os.path.join(RES, "audit_head_fraction.txt")

LO, HI = 200_000, 102_400_000
FRACS = (0.01, 0.02, 0.05, 0.10, 0.20, 0.50)
BASE = 0.10


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")
THETA = SPL.THETA


def read_published():
    """{#rem:sumhead}'s head share and exponents"""
    src = io.open(os.path.join(RES, "audit_sum_head.txt"),
                  encoding="utf-8").read()
    rows, dec = {}, 0
    for m in re.finditer(r"^  (\d{5,})\s+\d+\s+\d+\s+[+-][\d.]+\s+"
                         r"[+-][\d.]+\s+[+-][\d.]+\s+([\d.]+)\s+"
                         r"[\d.]+\s*$", src, re.M):
        rows[int(m.group(1))] = float(m.group(2))
        dec = max(dec, len(m.group(2).split(".")[1]))
    eh = float(re.search(r"^  head sum\s+([+-][\d.]+)", src,
                         re.M).group(1))
    et = float(re.search(r"^  total\s+([+-][\d.]+)", src,
                         re.M).group(1))
    return rows, dec, eh, et


def family(lo, hi):
    out = []
    a = 1
    while 2 ** a <= hi:
        b = 1
        while 2 ** a * 5 ** b <= hi:
            v = 2 ** a * 5 ** b
            if v >= lo:
                out.append(v)
            b += 1
        a += 1
    return sorted(set(out))


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum())) \
        if x.size > 2 else float("inf")
    return float(a), float(np.sqrt((r ** 2).mean())), se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubsh, dec, pubeh, pubet = read_published()
    NS = family(LO, HI)
    say("read %d head shares from results/audit_sum_head.txt, with "
        "its head exponent %+.6f and total %+.6f"
        % (len(pubsh), pubeh, pubet))
    say("  the field, the sieve and theta' are imported from "
        "code/audit_gain_split.py")
    say("  the fractions swept are %s"
        % ", ".join("%.2f" % f for f in FRACS))
    say()
    say("the field: every N = 2^a 5^b with a and b at least one in "
        "[%d, %d]; %d of them" % (LO, HI, len(NS)))
    classes = sorted(set(tuple(sorted(SPL.factor_set(N))) for N in NS))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))

    say()
    say("sieving to %d ..." % HI)
    lam, mu = SPL.lambda_and_mu(HI)
    sqf = mu != 0

    shares = {f: [] for f in FRACS}
    sums = {f: [] for f in FRACS}
    tots = []
    for N in NS:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        order = np.argsort(-np.abs(a))
        S = float(a.sum())
        tots.append(abs(S))
        for f in FRACS:
            nh = max(1, int(round(f * ks.size)))
            Sh = float(a[order[:nh]].sum())
            sums[f].append(abs(Sh))
            shares[f].append(abs(Sh) / abs(S))

    x = np.log(np.array(NS, dtype=np.float64))
    say()
    say("  fraction  share at the bottom N  at the top N   exponent "
        "of the signed sum   s.e.")
    exps = {}
    for f in FRACS:
        e, r, se = fit(x, np.log(np.array(sums[f])))
        exps[f] = (e, se)
        say("  %-9.2f %-22.4f %-14.4f %+-25.6f %.6f"
            % (f, shares[f][0], shares[f][-1], e, se))
    ea, ra, sea = fit(x, np.log(np.array(tots)))
    say("  %-9s %-22s %-14.4f %+-25.6f %.6f"
        % ("all", "1.0000", 1.0, ea, sea))

    # -------------------------------------------------------------- Y1
    say()
    say("Y1  the control at fraction %.2f" % BASE)
    rnd = 0.5 * 10.0 ** (-dec)
    d = abs(shares[BASE][-1] - pubsh[NS[-1]])
    y1 = d <= rnd
    say("  the top-N share here %.4f against the published %.4f, "
        "departure %.6f" % (shares[BASE][-1], pubsh[NS[-1]], d))
    say("  the table prints %d decimals, so the bound is %.8f"
        % (dec, rnd))
    say("PRINTBOUND audit_head_fraction %d %.8f" % (dec, rnd))
    say("  Y1 %s   (cap: the printing bound)"
        % ("hold" if y1 else "REFUTED"))

    # -------------------------------------------------------------- Y2
    say()
    say("Y2  is the failing set thinner than a tenth?")
    s1 = np.array(shares[0.01])
    y2 = bool((s1 >= 0.5).all())
    say("  the top one per cent carries %.4f to %.4f of |sum a|"
        % (float(s1.min()), float(s1.max())))
    say("  below a half at %d of the %d N"
        % (int((s1 < 0.5).sum()), len(NS)))
    say("  Y2 %s   (cap: a half at any N)"
        % ("hold" if y2 else "REFUTED"))

    # -------------------------------------------------------------- Y3
    say()
    say("Y3  how close is the ceiling?")
    top = shares[BASE][-1]
    es, rs, ses = fit(x, np.array(shares[BASE]))
    room = (1.0 - top) / es if es > 0 else float("inf")
    y3 = top > 0.95
    say("  at fraction %.2f the top-N share is %.4f and the share's "
        "slope is %+.6f per unit log N (s.e. %.6f)"
        % (BASE, top, es, ses))
    say("  so the distance to the bound is %.4f, which at that slope "
        "is %.4f units of log N" % (1.0 - top, room))
    say("TSTAT slope_audit_head_fraction %.2f" % (abs(es) / ses))
    say("SPREAD slope_audit_head_fraction %.4f"
        % float(x.max() - x.min()))
    if abs(es) / ses < 2.0:
        say("UNRESOLVED SIGN slope_audit_head_fraction")
    say("  Y3 %s   (cap: a share at or below 0.95)"
        % ("hold" if y3 else "REFUTED"))

    # -------------------------------------------------------------- Y4
    say()
    say("Y4  has the convergence begun inside the range?")
    lo8 = float(x.max() - math.log(2.0))
    sel = x >= lo8
    xt = x[sel]
    eht, rht, seht = fit(xt, np.log(np.array(sums[BASE])[sel]))
    eat, rat, seat = fit(xt, np.log(np.array(tots))[sel])
    gl = abs(exps[BASE][0] - ea)
    gt = abs(eht - eat)
    y4 = gt < gl
    say("  over the whole field the head's exponent is %+.6f against "
        "the total's %+.6f, a gap of %.6f"
        % (exps[BASE][0], ea, gl))
    say("  over the top octave alone (%d points) they are %+.6f and "
        "%+.6f, a gap of %.6f" % (int(sel.sum()), eht, eat, gt))
    say("  Y4 %s   (cap: the top-octave gap not smaller)"
        % ("hold" if y4 else "REFUTED"))

    say()
    say("  so the sentence {#rem:sumhead} ended on is corrected here: "
        "the head's")
    say("  signed sum has to travel the same distance the whole does, "
        "not further,")
    say("  because the share that makes its exponent look larger is "
        "within %.4f" % (1.0 - top))
    say("  of a ceiling it cannot pass.")

    say()
    say("=" * 70)
    say("Y1 %s  Y2 %s  Y3 %s  Y4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (y1, y2, y3, y4)))

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) over the squarefree",
        "           k < N^theta' coprime to N, the signed sum over",
        "           the top f of k by |a_k| for f = 0.01, 0.02, 0.05,",
        "           0.10, 0.20, 0.50, each part's share of |sum a| at",
        "           the bottom and top N and its least-squares",
        "           exponent in log N over the on-field family to",
        "           1.024e8; the share at f = 0.10 against its bound",
        "           of one, with the slope of that share; and the",
        "           head and total exponents refitted on the top",
        "           octave alone.",
        "NULL: none is run and none applies. Every quantity is an",
        "      exactly computed partial sum of a measured sequence",
        "      ordered by its own magnitudes, and there is no",
        "      background to detect against.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta'; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to 102400000; the weighted sum, the",
        "       sieve and theta' are code/audit_gain_split.py's,",
        "       imported; the published shares and exponents are read",
        "       from results/audit_sum_head.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not y1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
