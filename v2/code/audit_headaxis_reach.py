# -*- coding: utf-8 -*-
r"""
What selects the head, and what aligns it, on eighty-one points.

WHAT IS AT STAKE

Remark {#rem:splitreach} narrowed item 4(b) to the head: the top tenth
of k by |a_k| carries one sign four times in five, its alignment
decays too slowly to reach a coin anywhere computable, and its gain
falls further behind the square-root reference at every reach. A proof
would have to control that set, so what the set IS matters.

Two remarks answer that and both stand on eight doublings.
{#rem:headidentity} says the head is not the small-k end -- its
overlap with the smallest tenth of k is 0.2174 to 0.3263.
{#rem:headsign} says the alignment does not live on the k axis either:
cutting the same deciles on each factor of |a_k| = (log k) T_k |I_k|
separately, the top-minus-bottom spread of the negative share is
+0.3389 on |I|, +0.1537 on T and +0.0126 on k at the top of the eight.

Those are one N's deciles on eight N. This asks the same question on
the field to 1.024e8, and adds the one measurement neither remark
made: which of the three factors SELECTS the head. Alignment and
selection are different questions -- a factor can decide which k are
large without deciding which are negative -- and item 4(b) needs both,
because the set to be controlled is chosen by size and obstructs by
sign.

BACKS: Remark {#rem:headaxis} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  The control. At N = 25600000 the four negative-share spreads
      reproduce results/audit_head_sign.txt to 0.0001, and the head's
      one-sign fraction reproduces results/audit_split_reach.txt.
  Y2  The alignment ordering survives the reach: at the top N the
      negative-share spread is larger on |I| than on T, and larger on
      T than on k.
  Y3  The alignment does not wash out: the |I| spread's least-squares
      slope in log N is not resolved negative.
  Y4  Selection follows alignment. In the variance decomposition of
      log|a_k| = log log k + log T_k + log|I_k|, the log|I| term
      carries the largest share of the variance at every N.
  Y5  And k stays out of it: at the top N the k-decile spread of the
      negative share is inside two binomial deviations of zero.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED above 0.0001 on any of the five. Then this is not the
      statistic those remarks measured. THIS ONE GATES.
  Y2  REFUTED if either inequality fails at the top N. Then the axis
      reading of {#rem:headsign} was a property of the doublings and
      the head's alignment has to be re-attributed.
  Y3  REFUTED if the slope is resolved negative. Then the imbalance
      axis is losing its grip on the sign as N grows, which would put
      a decay mechanism where {#rem:splitreach} found only a slow one
      -- and would be the first thing in this programme pointing at
      the head dissolving.
  Y4  REFUTED if T or log log k carries more variance anywhere. Then
      the head is selected by mass and aligned by imbalance, two
      different sets overlapping, and no single axis describes it.
  Y5  REFUTED if the k spread clears two binomial deviations. Then k
      does carry some of the sign after all, at a reach the doublings
      could not see, and the arithmetic of the dilation is back in
      play.

  Y1 gates. Y2 to Y5 are the measurement and do not gate.

  NO NULL IS RUN for the decile cuts, which are deterministic
  partitions of a measured sequence, and Y5's cap IS the null for the
  sign: under independent signs every decile's negative share is
  binomial about one half, and that binomial deviation is the cap. The
  coin arms for the gain are audit_crossk_reference.py and
  lab_gain_opposition.py.
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
OUT = os.path.join(RES, "audit_headaxis_reach.txt")

LO, HI = 200_000, 102_400_000
BLOCKS = 10


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HS = module("audit_head_sign")
SPL = HS.SPL
THETA = SPL.THETA
HEAD = SPL.HEAD


def read_published():
    """the four spreads {#rem:headsign} printed, and the head fraction"""
    src = io.open(os.path.join(RES, "audit_head_sign.txt"),
                  encoding="utf-8").read()
    i = src.index("spread from half, top decile minus bottom decile:")
    sp = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) != 2:
            break
        sp[f[0]] = float(f[1])
    src2 = io.open(os.path.join(RES, "audit_split_reach.txt"),
                   encoding="utf-8").read()
    ag = {}
    for m in re.finditer(r"^  (\d{5,})\s+\d+\s+\d+\s+[\d.]+\s+[\d.]+"
                         r"\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s+[\d.]+\s*$",
                         src2, re.M):
        ag[int(m.group(1))] = float(m.group(2))
    return sp, ag


def family(lo, hi):
    """the field: N = 2^a 5^b with a and b at least one"""
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


def decile_spread(key, neg):
    """top-minus-bottom spread of the negative share on one axis

    The cut is audit_head_sign.py's exactly: sort DESCENDING on the
    axis, edges at int(round(d*n/BLOCKS)), decile 0 the largest, and
    the spread decile 0 minus decile 9. A floor instead of the round
    moves single k across an edge and the control will not reproduce.
    """
    order = np.argsort(-key)
    n = order.size
    edges = [int(round(d * n / BLOCKS)) for d in range(BLOCKS + 1)]
    sh = [float(neg[order[edges[d]:edges[d + 1]]].mean())
          for d in range(BLOCKS)]
    return sh[0] - sh[-1], sh


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubsp, pubag = read_published()
    NS = family(LO, HI)
    CTRL = 25_600_000

    say("read the four spreads from results/audit_head_sign.txt: %s"
        % ", ".join("%s %+.4f" % (k, v) for k, v in pubsp.items()))
    say("  and %d head one-sign fractions from "
        "results/audit_split_reach.txt" % len(pubag))
    say("  the field, the sieve, theta' and the head fraction are "
        "imported from code/audit_gain_split.py, and the signed split")
    say("  of the inner sum from code/audit_head_sign.py")
    say()
    say("the field: every N = 2^a 5^b with a and b at least one in "
        "[%d, %d]" % (LO, HI))
    classes = sorted(set(tuple(sorted(SPL.factor_set(N))) for N in NS))
    say("  %d of them; the control N is %d" % (len(NS), CTRL))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))

    say()
    say("sieving to %d ..." % HI)
    lam, mu = SPL.lambda_and_mu(HI)
    sqf = mu != 0

    rows = []
    tops = {}
    say()
    say("  N            #k      onesign  spread |I|   spread T     "
        "spread k     var log|I|  var T   var loglogk")
    for N in NS:
        ks, P, M, C = HS.signed_parts(N, lam, mu, sqf)
        H = P - M
        T = P + M
        a = np.log(ks.astype(np.float64)) * H
        w = np.abs(a)
        order = np.argsort(-w)
        nh = max(1, int(round(HEAD * ks.size)))
        hd = order[:nh]
        sh = np.sign(a[hd])
        agree = max(float((sh > 0).mean()), float((sh < 0).mean()))
        ok = T > 0
        I = np.zeros_like(a)
        I[ok] = H[ok] / T[ok]
        neg = (I < 0).astype(np.float64)
        sI, dI = decile_spread(np.abs(I), neg)
        sT, dT = decile_spread(T, neg)
        sK, dK = decile_spread(ks.astype(np.float64), neg)

        # the variance decomposition of log|a|, on the k where it exists
        live = (w > 0) & ok
        u1 = np.log(np.log(ks[live].astype(np.float64)))
        u2 = np.log(T[live])
        u3 = np.log(np.abs(I[live]))
        s = u1 + u2 + u3
        vs = float(np.var(s))
        c1 = float(np.cov(u1, s, bias=True)[0, 1]) / vs
        c2 = float(np.cov(u2, s, bias=True)[0, 1]) / vs
        c3 = float(np.cov(u3, s, bias=True)[0, 1]) / vs
        rows.append((N, int(ks.size), nh, agree, sI, sT, sK,
                     c3, c2, c1, int(live.sum())))
        if N == CTRL:
            tops[N] = (dI, dT, dK)
        say("  %-12d %-7d %-8.4f %+-12.4f %+-12.4f %+-12.4f "
            "%-11.4f %-7.4f %.4f"
            % (N, ks.size, agree, sI, sT, sK, c3, c2, c1))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))
    top = rows[-1]

    # -------------------------------------------------------------- Y1
    say()
    say("Y1  the control at N = %d" % CTRL)
    ctrl = [r for r in rows if r[0] == CTRL][0]
    say("  axis     here        published")
    worst = 0.0
    for nm, here in (("|I|", ctrl[4]), ("T", ctrl[5]), ("k", ctrl[6])):
        pv = pubsp.get(nm, float("nan"))
        worst = max(worst, abs(here - pv))
        say("  %-8s %+-11.4f %+.4f" % (nm, here, pv))
    da = abs(ctrl[3] - pubag.get(CTRL, float("nan")))
    worst = max(worst, da)
    say("  the head's one-sign fraction %.4f against the published "
        "%.4f" % (ctrl[3], pubag.get(CTRL, float("nan"))))
    y1 = worst <= 0.0001
    say("  worst departure %.6f" % worst)
    say("  Y1 %s   (tol 0.0001)" % ("hold" if y1 else "REFUTED"))

    # -------------------------------------------------------------- Y2
    say()
    say("Y2  does the axis ordering survive at the new top?")
    say("  at N = %d the spreads are |I| %+.4f, T %+.4f, k %+.4f"
        % (top[0], top[4], top[5], top[6]))
    y2 = top[4] > top[5] > top[6]
    say("  Y2 %s   (cap: either inequality failing)"
        % ("hold" if y2 else "REFUTED"))

    # -------------------------------------------------------------- Y3
    say()
    say("Y3  is the imbalance axis losing its grip?")
    eI, rI, seI = fit(x, np.array([r[4] for r in rows]))
    say("  the |I| spread runs %+.4f to %+.4f over the field"
        % (min(r[4] for r in rows), max(r[4] for r in rows)))
    say("  its least-squares slope in log N = %+.6f, s.e. %.6f, "
        "t = %.2f, r.m.s. %.6f" % (eI, seI, abs(eI) / seI, rI))
    y3 = not (eI < 0.0 and abs(eI) / seI >= 2.0)
    say("TSTAT slope_audit_headaxis_reach %.2f" % (abs(eI) / seI))
    say("SPREAD slope_audit_headaxis_reach %.4f"
        % float(x.max() - x.min()))
    if abs(eI) / seI < 2.0:
        say("UNRESOLVED SIGN slope_audit_headaxis_reach")
    eT, rT, seT = fit(x, np.array([r[5] for r in rows]))
    say("  for comparison the T spread's slope is %+.6f, s.e. %.6f, "
        "t = %.2f" % (eT, seT, abs(eT) / seT))
    say("TSTAT slope_headaxis_T %.2f" % (abs(eT) / seT))
    say("SPREAD slope_headaxis_T %.4f" % float(x.max() - x.min()))
    if abs(eT) / seT < 2.0:
        say("UNRESOLVED SIGN slope_headaxis_T")
    say("  Y3 %s   (cap: resolved negative at two standard errors)"
        % ("hold" if y3 else "REFUTED"))

    # -------------------------------------------------------------- Y4
    say()
    say("Y4  which factor SELECTS the head?")
    say("  the variance of log|a| = log log k + log T + log|I| "
        "decomposed")
    say("  by covariance, so the three shares sum to one:")
    say("  N            log|I|      log T       log log k   sum")
    for r in (rows[0], rows[len(rows) // 2], top):
        say("  %-12d %-11.4f %-11.4f %-11.4f %.4f"
            % (r[0], r[7], r[8], r[9], r[7] + r[8] + r[9]))
    y4 = all(r[7] > r[8] and r[7] > r[9] for r in rows)
    worstI = min(r[7] for r in rows)
    worstT = max(r[8] for r in rows)
    say("  over the field the log|I| share runs %.4f to %.4f and the "
        "log T share %.4f to %.4f"
        % (worstI, max(r[7] for r in rows),
           min(r[8] for r in rows), worstT))
    eV, rV, seV = fit(x, np.array([r[7] for r in rows]))
    say("  the log|I| share's slope in log N = %+.6f, s.e. %.6f, "
        "t = %.2f" % (eV, seV, abs(eV) / seV))
    say("TSTAT slope_headaxis_varI %.2f" % (abs(eV) / seV))
    say("SPREAD slope_headaxis_varI %.4f" % float(x.max() - x.min()))
    if abs(eV) / seV < 2.0:
        say("UNRESOLVED SIGN slope_headaxis_varI")
    say("  Y4 %s   (cap: T or log log k carrying more at any N)"
        % ("hold" if y4 else "REFUTED"))

    # -------------------------------------------------------------- Y5
    say()
    say("Y5  does k carry any of the sign at the new reach?")
    nd = top[1] // BLOCKS
    sd = math.sqrt(0.25 / nd) * math.sqrt(2.0)
    say("  at the top N each decile holds about %d of the k, so the "
        "binomial deviation of a top-minus-bottom difference is %.6f"
        % (nd, sd))
    say("  the k spread there is %+.4f, which is %.2f of those"
        % (top[6], abs(top[6]) / sd))
    y5 = abs(top[6]) <= 2.0 * sd
    say("  and on |I| the same ratio is %.2f" % (abs(top[4]) / sd))
    say("  Y5 %s   (cap: two binomial deviations)"
        % ("hold" if y5 else "REFUTED"))

    say()
    say("  the deciles at the control N, cut on each axis, as "
        "{#rem:headsign} printed them:")
    dI, dT, dK = tops[CTRL]
    say("  decile   by |I|      by T        by k")
    for d in range(BLOCKS):
        say("  %-8d %-11.4f %-11.4f %.4f" % (d + 1, dI[d], dT[d],
                                             dK[d]))
    say("PERN headaxis_spread_I %d %.4f %.4f"
        % (len(rows), min(r[4] for r in rows),
           max(r[4] for r in rows)))
    say("PERN headaxis_varshare_I %d %.4f %.4f"
        % (len(rows), worstI, max(r[7] for r in rows)))
    rr = [r[4] / r[7] for r in rows]
    say("RATIO headaxis_spread_I headaxis_varshare_I %.4f %.4f"
        % (min(rr), max(rr)))

    say()
    say("=" * 70)
    say("Y1 %s  Y2 %s  Y3 %s  Y4 %s  Y5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (y1, y2, y3, y4, y5)))

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) with the inner sum split",
        "           by the sign of mu(m) into P_k and M_k, so",
        "           H = P - M, T = P + M and I = H/T, the fraction of",
        "           k with a_k < 0 in each decile of the k cut",
        "           separately on |I|, on T and on k, and the",
        "           top-minus-bottom spread of that fraction on each",
        "           axis; the head's one-sign fraction; and the",
        "           covariance decomposition of the variance of",
        "           log|a| = log log k + log T + log|I| into the three",
        "           shares it sums to; all at every on-field N to",
        "           1.024e8, with each spread's least-squares slope in",
        "           log N.",
        "NULL: none is run for the decile cuts, which are",
        "      deterministic partitions of a measured sequence. For",
        "      the sign the null IS the cap Y5 uses: under independent",
        "      signs every decile's negative share is binomial about",
        "      one half. The coin arms for the gain are",
        "      audit_crossk_reference.py and lab_gain_opposition.py.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta'; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to 102400000; the signed split of the",
        "       inner sum is code/audit_head_sign.py's signed_parts,",
        "       imported; the field, the sieve, theta' and the head",
        "       fraction come from code/audit_gain_split.py; the",
        "       published spreads are read from",
        "       results/audit_head_sign.txt and the published one-sign",
        "       fractions from results/audit_split_reach.txt.",
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
