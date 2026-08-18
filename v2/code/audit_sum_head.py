# -*- coding: utf-8 -*-
r"""
Which part of the k carries the signed sum, not the mass.

WHAT IS AT STAKE

Two lines of this investigation have not been joined. {#rem:splitreach}
found the deficit localised in the head -- the top tenth of k by
|a_k| -- and {#rem:denominator} found the demand is about |sum a|,
the signed total. Every head measurement so far has been about GAINS
and MASS: the head carries 0.3337 of sum|a| at the top N and its own
gain is 1.7957. Nobody has asked how much of the SIGNED sum the head
carries, and that is the quantity item 4(b) is about.

The published numbers already suggest an answer -- l1_head/G_head
against l1/G is about six sevenths -- but a ratio of two gains is not
a measurement of a signed sum, because the head's sum and the tail's
can point opposite ways and partly cancel. Measured directly, three
things follow or fail: whether the two parts agree in sign, whether
the head is most of the total, and whether the tail is of lower order.
If the tail is negligible, item 4(b) stops being a statement about
sum a and becomes one about the head's signed sum alone.

BACKS: Remark {#rem:sumhead} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  X1  The control. |sum a| reproduces results/audit_denominator.txt's
      exponent alpha inside the bound its printing forces.
  X2  The head's signed sum has the same sign as the total at every
      N, so the two parts do not oppose.
  X3  The head carries most of it: |sum over the head| is at least
      half of |sum a| at every N.
  X4  And the tail's signed sum is of lower order: its least-squares
      exponent in log N is below alpha at two standard errors.

REFUTATION RULE (fixed before the run)

  X1  REFUTED outside the printing bound. Then this is not the sum
      {#rem:denominator} measured. THIS ONE GATES.
  X2  REFUTED if the signs differ anywhere. Then the total is a
      difference of two opposing parts and the head's share of it is
      not a share at all.
  X3  REFUTED if the head carries less than half anywhere. Then the
      signed sum lives in the tail even though the mass and the gain
      deficit live in the head, and the two lines do not join.
  X4  REFUTED if the tail is not resolved below alpha. Then both
      parts grow at the same rate, item 4(b) stays a statement about
      the whole k-range, and localising it in the head buys nothing
      for the quantity that matters.

  X1 gates. X2 to X4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. The split is a deterministic
  partition of a measured sequence by its own magnitudes, and every
  part is an exactly computed sum; the arms that randomise signs on
  these magnitudes are audit_crossk_reference.py and
  lab_gain_opposition.py, and they are arms for the gain, not for
  this partition.
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
OUT = os.path.join(RES, "audit_sum_head.txt")

LO, HI = 200_000, 102_400_000


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")
THETA = SPL.THETA
HEAD = SPL.HEAD


def read_published():
    """alpha and e(l1) from {#rem:denominator}"""
    src = io.open(os.path.join(RES, "audit_denominator.txt"),
                  encoding="utf-8").read()
    out, dec = {}, 0
    for nm in ("|sum a|", "l1"):
        m = re.search(r"^  " + re.escape(nm) +
                      r"\s+([+-][\d.]+)\s+([\d.]+)\s+[\d.]+\s*$",
                      src, re.M)
        out[nm] = (float(m.group(1)), float(m.group(2)))
        dec = max(dec, len(m.group(1).split(".")[1]))
    return out, dec


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

    pub, dec = read_published()
    NS = family(LO, HI)
    say("read alpha = %+.6f (s.e. %.6f) and e(l1) = %+.6f (s.e. %.6f)"
        % (pub["|sum a|"][0], pub["|sum a|"][1], pub["l1"][0],
           pub["l1"][1]))
    say("  from results/audit_denominator.txt; the field, the sieve, "
        "theta' and the head fraction are imported from")
    say("  code/audit_gain_split.py")
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

    rows = []
    say()
    say("  N            #k      head   sum a           head sum     "
        "   tail sum        head share  head mass")
    for N in NS:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        w = np.abs(a)
        nh = max(1, int(round(HEAD * ks.size)))
        order = np.argsort(-w)
        hd, tl = order[:nh], order[nh:]
        S = float(a.sum())
        Sh = float(a[hd].sum())
        St = float(a[tl].sum())
        share = abs(Sh) / abs(S)
        mass = float(w[hd].sum() / w.sum())
        rows.append((N, int(ks.size), nh, S, Sh, St, share, mass))
        say("  %-12d %-7d %-6d %+-15.1f %+-15.1f %+-15.1f %-11.4f "
            "%.4f" % (N, ks.size, nh, S, Sh, St, share, mass))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))
    top = rows[-1]

    # -------------------------------------------------------------- X1
    say()
    say("X1  the control")
    ea, ra, sea = fit(x, np.log(np.array([abs(r[3]) for r in rows])))
    rnd = 0.5 * 10.0 ** (-dec)
    d = abs(ea - pub["|sum a|"][0])
    x1 = d <= rnd
    say("  alpha here %+.6f against the published %+.6f, departure "
        "%.8f" % (ea, pub["|sum a|"][0], d))
    say("  the table prints %d decimals, so the bound is %.8f"
        % (dec, rnd))
    say("PRINTBOUND audit_sum_head %d %.8f" % (dec, rnd))
    say("  X1 %s   (cap: the printing bound)"
        % ("hold" if x1 else "REFUTED"))

    # -------------------------------------------------------------- X2
    say()
    say("X2  do the two parts point the same way?")
    same = sum(1 for r in rows if (r[3] < 0) == (r[4] < 0))
    both = sum(1 for r in rows if (r[3] < 0) == (r[5] < 0))
    x2 = same == len(rows)
    say("  the head agrees in sign with the total at %d of %d N, and "
        "the tail at %d of %d" % (same, len(rows), both, len(rows)))
    say("  X2 %s" % ("hold" if x2 else "REFUTED"))

    # -------------------------------------------------------------- X3
    say()
    say("X3  how much of the signed sum is the head's?")
    sh = np.array([r[6] for r in rows])
    x3 = bool((sh >= 0.5).all())
    say("  the head's share of |sum a| runs %.4f to %.4f over the "
        "field" % (float(sh.min()), float(sh.max())))
    say("  at the top N it carries %+.1f of %+.1f while holding %.4f "
        "of the mass" % (top[4], top[3], top[7]))
    es, rs, ses = fit(x, sh)
    say("  its share has least-squares slope in log N = %+.6f, "
        "s.e. %.6f, t = %.2f" % (es, ses, abs(es) / ses))
    say("TSTAT slope_audit_sum_head %.2f" % (abs(es) / ses))
    say("SPREAD slope_audit_sum_head %.4f" % float(x.max() - x.min()))
    if abs(es) / ses < 2.0:
        say("UNRESOLVED SIGN slope_audit_sum_head")
    say("SITSIN sumhead_headshare %.4f %.4f"
        % (float(min(r[7] for r in rows)),
           float(max(r[7] for r in rows))))
    say("  X3 %s   (cap: a half at any N)"
        % ("hold" if x3 else "REFUTED"))

    # -------------------------------------------------------------- X4
    say()
    say("X4  and is the tail of lower order?")
    eh, rh, seh = fit(x, np.log(np.array([abs(r[4]) for r in rows])))
    et, rt, set_ = fit(x, np.log(np.array([abs(r[5]) for r in rows])))
    say("  part          exponent      s.e.        r.m.s.")
    for nm, e, se, rr in (("head sum", eh, seh, rh),
                          ("tail sum", et, set_, rt),
                          ("total", ea, sea, ra)):
        say("  %-13s %+-13.6f %-11.6f %.6f" % (nm, e, se, rr))
    dd = math.sqrt(set_ * set_ + sea * sea)
    x4 = (ea - et) > 2.0 * dd
    say("  the tail sits below the total by %+.6f, which is %.2f "
        "standard errors of the difference"
        % (et - ea, abs(et - ea) / dd))
    say("TSTAT slope_sumhead_tail %.2f" % (abs(et) / set_))
    say("SPREAD slope_sumhead_tail %.4f" % float(x.max() - x.min()))
    if abs(et) / set_ < 2.0:
        say("UNRESOLVED SIGN slope_sumhead_tail")
    say("  X4 %s   (cap 2 standard errors)"
        % ("hold" if x4 else "REFUTED"))

    say()
    say("  what that leaves item 4(b). The demand is |sum a| down to "
        "l2 order,")
    say("  which the identity e(G) = e(l1) - alpha puts at "
        "%+.6f against" % (pub["l1"][0] - 0.283586))
    say("  the measured %+.6f, a distance of %+.6f. On the split "
        "above that"
        % (ea, ea - (pub["l1"][0] - 0.283586)))
    say("  distance belongs to the head's signed sum, whose own "
        "exponent is %+.6f." % eh)

    say()
    say("=" * 70)
    say("X1 %s  X2 %s  X3 %s  X4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (x1, x2, x3, x4)))

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) over the squarefree",
        "           k < N^theta' coprime to N, the signed sum sum a",
        "           split into the top tenth of k by |a_k| and the",
        "           rest: each part's signed total, its sign against",
        "           the whole, the head's share of |sum a| and of",
        "           sum|a|, and each part's least-squares exponent in",
        "           log N over the on-field family to 1.024e8.",
        "NULL: none is run and none applies. The split is a",
        "      deterministic partition of a measured sequence by its",
        "      own magnitudes and every part is an exactly computed",
        "      sum. The arms that randomise signs on these magnitudes",
        "      are audit_crossk_reference.py and",
        "      lab_gain_opposition.py, and they are arms for the gain",
        "      rather than for this partition.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta'; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to 102400000; the weighted sum, the",
        "       sieve, theta' and the head fraction are",
        "       code/audit_gain_split.py's, imported; alpha and",
        "       e(l1) are read from results/audit_denominator.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not x1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
