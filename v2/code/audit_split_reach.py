# -*- coding: utf-8 -*-
r"""
The head/tail split of the gain, on the field and at the new reach.

WHAT IS AT STAKE

Remark {#rem:gainsplit} located item 4(b)'s deficit by mass rank: on
the eight doublings the top tenth of k by |a_k| has gain exponent
+0.077963 and the bottom nine tenths +0.340006, so the small terms
cancel BETTER than square root and the large ones hardly cancel at
all. Its sharpest sentence is about the head's signs -- the fraction
of the head carrying one sign runs 1.0000, 0.9783, 0.9118, 0.8800,
0.8716, 0.8721, 0.8545, 0.8274 -- "a positive proportion of the
dilations, a tenth of them, growing like N^theta', carries one sign".

Every one of those numbers is eight points of a doubling family, and
that fraction is falling monotonically across them. Whether it is
falling to a positive limit or to one half is the difference between
an obstruction that persists and one that is a small-N effect, and
eight points spanning a factor of 128 cannot tell. Remark
{#rem:fieldreach} has just shown the field runs to 1.024e8 for seconds
per N, so the same split can be put on 81 points over a spread of
6.2383 in log N.

This is the axis OPEN item 4(b) is actually on: e(G) has to close
+0.134019 onto the measured e(l1/l2) = +0.283586, and {#rem:gainsplit}
is the only remark that says WHERE in the k that closing would have to
happen.

BACKS: Remark {#rem:splitreach} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  The control. At the eight doublings the whole-range gain, the
      head's gain, the tail's gain, the head's mass share and the
      head's one-sign fraction reproduce results/audit_gain_split.txt
      to the precision it prints.
  Z2  The ordering survives the field and the reach: on the 81
      on-field N, e(head) < e(whole) < e(tail).
  Z3  The tail still beats square root: e(tail) exceeds theta'/2 by
      more than two of its own standard errors.
  Z4  The head's one-sign fraction keeps falling, and does not reach
      a coin. Its least-squares slope in log N is negative at two
      standard errors, AND at the top N it still exceeds one half by
      more than two binomial standard deviations of the head's size.
  Z5  The head's mass share is flat, as {#rem:gainprofile} measured
      the top block's share to be: its exponent in log N is within
      two standard errors of zero.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED on any departure above 0.0001. Then this is not the
      statistic {#rem:gainsplit} measured and nothing below may be
      compared with it. THIS ONE GATES.
  Z2  REFUTED if the ordering breaks anywhere. Then the deficit is
      not localised by size on the field, and {#rem:gainsplit}'s
      reading was a property of the doublings.
  Z3  REFUTED if the tail's excess falls to two standard errors or
      less. Then the small terms do not beat square root after all
      and the split says only that the head is worse than the tail,
      not that the tail is already good enough.
  Z4  REFUTED either way, and they mean opposite things. If the slope
      is not resolved negative, the fraction is not falling and the
      one-sign head is a fixed feature. If the top N is within two
      binomial standard deviations of one half, the head's signs are
      a coin at the reach we can compute and {#rem:gainsplit}'s
      "carries one sign" is a small-N statement.
  Z5  REFUTED if the share's exponent is resolved away from zero.
      Then the head is not a fixed proportion of the mass and the
      exponents on it are not comparable across N.

  Z1 gates. Z2 to Z5 are the measurement and do not gate.

  NO NULL IS RUN for the split, which is a deterministic partition of
  a measured sequence, and none is needed for the one-sign fraction:
  its null is the binomial one half, which Z4 uses as its own cap.
  The coin arms for the gain itself were run in
  audit_crossk_reference.py and lab_gain_opposition.py.
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
OUT = os.path.join(RES, "audit_split_reach.txt")

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
    """{#rem:gainsplit}'s eight rows and its three exponents"""
    src = io.open(os.path.join(RES, "audit_gain_split.txt"),
                  encoding="utf-8").read()
    rows = {}
    for m in re.finditer(r"^  N = (\d+)\s+#k = (\d+)\s+head (\d+)\s+"
                         r"G ([\d.]+)\s+head ([\d.]+)\s+tail "
                         r"([\d.]+)\s+mass ([\d.]+)\s*$", src, re.M):
        rows[int(m.group(1))] = (int(m.group(2)), int(m.group(3)),
                                 float(m.group(4)), float(m.group(5)),
                                 float(m.group(6)), float(m.group(7)))
    ag = {}
    i = src.index("N            top-decile share  same sign in the "
                  "head")
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) != 3 or not f[0].isdigit():
            break
        ag[int(f[0])] = float(f[2])
    m = re.search(r"^GAINSPLIT crossk ([+-][\d.]+) ([+-][\d.]+) "
                  r"([+-][\d.]+)\s*$", src, re.M)
    return (rows, ag, float(m.group(1)), float(m.group(2)),
            float(m.group(3)))


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


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub, pubag, pubh, pubt, pubw = read_published()
    ceil_ = THETA / 2.0
    NS = family(LO, HI)
    RUNGS = [N for N in SPL.NS if N in NS]

    say("read %d published rows from results/audit_gain_split.txt"
        % len(pub))
    say("  with its head, tail and whole exponents %+.6f %+.6f %+.6f"
        % (pubh, pubt, pubw))
    say("  the field, the sieve, theta' and the head fraction are "
        "imported from code/audit_gain_split.py")
    say("  the arithmetic reference for a gain exponent is "
        "theta'/2 = %.4f" % ceil_)
    say()
    say("the field: every N = 2^a 5^b with a and b at least one in "
        "[%d, %d]" % (LO, HI))
    classes = sorted(set(tuple(sorted(SPL.factor_set(N))) for N in NS))
    say("  %d of them, and %d of the doublings sit inside for the "
        "control" % (len(NS), len(RUNGS)))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))

    say()
    say("sieving to %d ..." % HI)
    lam, mu = SPL.lambda_and_mu(HI)
    sqf = mu != 0

    rows = []
    say()
    say("  N            #k      head   G        head     tail     "
        "mass    onesign  overlap")
    for N in NS:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        w = np.abs(a)
        order = np.argsort(-w)
        nh = max(1, int(round(HEAD * ks.size)))
        hd, tl = order[:nh], order[nh:]

        def gain(idx):
            s = abs(float(a[idx].sum()))
            return float(np.abs(a[idx]).sum()) / s if s > 0 else \
                float("inf")

        sh = np.sign(a[hd])
        agree = max(float((sh > 0).mean()), float((sh < 0).mean()))
        small = set(int(v) for v in np.argsort(ks)[:nh])
        ovl = len(small & set(int(v) for v in hd)) / float(nh)
        rows.append((N, int(ks.size), nh, gain(np.arange(ks.size)),
                     gain(hd), gain(tl), float(w[hd].sum() / w.sum()),
                     agree, ovl))
        say("  %-12d %-7d %-6d %-8.4f %-8.4f %-8.4f %-7.4f %-8.4f "
            "%.4f" % rows[-1])

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))

    # -------------------------------------------------------------- Z1
    say()
    say("Z1  the control at the doublings")
    say("  N            G here / pub          head           tail     "
        "      mass         onesign")
    worst = 0.0
    for r in rows:
        N = r[0]
        if N not in pub:
            continue
        _pk, _pn, pg, ph, pt, pm = pub[N]
        pa = pubag.get(N, float("nan"))
        worst = max(worst, abs(r[3] - pg), abs(r[4] - ph),
                    abs(r[5] - pt), abs(r[6] - pm), abs(r[7] - pa))
        say("  %-12d %.4f/%.4f  %.4f/%.4f  %.4f/%.4f  %.4f/%.4f  "
            "%.4f/%.4f"
            % (N, r[3], pg, r[4], ph, r[5], pt, r[6], pm, r[7], pa))
    z1 = worst <= 0.0001
    say("  worst departure %.6f" % worst)
    say("  Z1 %s   (tol 0.0001)" % ("hold" if z1 else "REFUTED"))

    # -------------------------------------------------------- Z2 / Z3
    ew, rw, sw = fit(x, np.log(np.array([r[3] for r in rows])))
    eh, rh, sh_ = fit(x, np.log(np.array([r[4] for r in rows])))
    et, rt, st = fit(x, np.log(np.array([r[5] for r in rows])))
    say()
    say("Z2/Z3  the exponent on each part, on %d points over %.4f in "
        "log N" % (len(rows), float(x.max() - x.min())))
    say("  part           exponent     s.e.       r.m.s.     "
        "published (8 doublings)")
    for nm, e, se, rms, pv in (("whole", ew, sw, rw, pubw),
                               ("head tenth", eh, sh_, rh, pubh),
                               ("tail", et, st, rt, pubt)):
        say("  %-14s %+-12.6f %-10.6f %-10.6f %+.6f"
            % (nm, e, se, rms, pv))
        lab = "slope_splitreach_" + nm.replace(" ", "_")
        say("SCATTER %s %.4f" % (lab, rms))
        say("TSTAT %s %.2f" % (lab, abs(e) / se))
        say("SPREAD %s %.4f" % (lab, float(x.max() - x.min())))
        if abs(e) / se < 2.0:
            say("UNRESOLVED SIGN %s" % lab)
    z2 = eh < ew < et
    z3 = (et - ceil_) > 2.0 * st
    say("  the tail's excess over theta'/2 is %+.6f = %.2f standard "
        "errors" % (et - ceil_, (et - ceil_) / st))
    say("  Z2 the ordering head < whole < tail   %s"
        % ("hold" if z2 else "REFUTED"))
    say("  Z3 the tail beats square root   %s   (cap 2 standard "
        "errors)" % ("hold" if z3 else "REFUTED"))
    say("GAINSPLIT crossk_reach %+.6f %+.6f %+.6f" % (eh, et, ew))
    say("SPLITOVERLAP crossk_reach %.4f %.4f"
        % (min(r[8] for r in rows), max(r[8] for r in rows)))

    # -------------------------------------------------------------- Z4
    say()
    say("Z4  is the head's one sign falling, and to what?")
    ag = np.array([r[7] for r in rows])
    ea, ra, sa = fit(x, ag)
    top = rows[-1]
    sd = math.sqrt(0.25 / top[2])
    say("  the fraction runs %.4f to %.4f over the field"
        % (float(ag.min()), float(ag.max())))
    say("  its least-squares slope in log N = %+.6f, s.e. %.6f, "
        "t = %.2f, r.m.s. %.6f" % (ea, sa, abs(ea) / sa, ra))
    say("  at the top N the head has %d members, so a coin's standard "
        "deviation is %.6f" % (top[2], sd))
    say("  the top fraction is %.4f, which is %.2f of those above one "
        "half" % (top[7], (top[7] - 0.5) / sd))
    z4a = ea < 0.0 and abs(ea) / sa >= 2.0
    z4b = (top[7] - 0.5) > 2.0 * sd
    z4 = z4a and z4b
    say("TSTAT slope_splitreach_onesign %.2f" % (abs(ea) / sa))
    say("SPREAD slope_splitreach_onesign %.4f"
        % float(x.max() - x.min()))
    say("TSTAT slope_audit_split_reach %.2f" % (abs(ea) / sa))
    say("SPREAD slope_audit_split_reach %.4f"
        % float(x.max() - x.min()))
    if abs(ea) / sa < 2.0:
        say("UNRESOLVED SIGN slope_splitreach_onesign")
    say("  the falling part %s and the above-a-coin part %s"
        % ("holds" if z4a else "is REFUTED",
           "holds" if z4b else "is REFUTED"))
    say("  Z4 %s   (cap 2 standard errors and 2 binomial deviations)"
        % ("hold" if z4 else "REFUTED"))

    # -------------------------------------------------------------- Z5
    say()
    say("Z5  is the head a fixed share of the mass?")
    ms = np.array([r[6] for r in rows])
    em, rm, sm = fit(x, np.log(ms))
    z5 = abs(em) <= 2.0 * sm
    say("  the share runs %.4f to %.4f, exponent %+.6f, s.e. %.6f, "
        "t = %.2f" % (float(ms.min()), float(ms.max()), em, sm,
                      abs(em) / sm))
    say("TSTAT slope_splitreach_mass %.2f" % (abs(em) / sm))
    say("SPREAD slope_splitreach_mass %.4f" % float(x.max() - x.min()))
    if abs(em) / sm < 2.0:
        say("UNRESOLVED SIGN slope_splitreach_mass")
    say("  Z5 %s   (cap 2 standard errors)"
        % ("hold" if z5 else "REFUTED"))
    say("PERN splitreach_head_mass %d %.4f %.4f"
        % (len(rows), float(ms.min()), float(ms.max())))
    say("PERN splitreach_head_onesign %d %.4f %.4f"
        % (len(rows), float(ag.min()), float(ag.max())))
    rr = [r[6] / r[7] for r in rows]
    say("RATIO splitreach_head_mass splitreach_head_onesign %.4f %.4f"
        % (min(rr), max(rr)))

    # ------------------------------------------- not pre-registered
    say()
    say("X1  how far the head is from cancelling at all")
    say("  (written after Z2; not pre-registered)")
    say("  N            head G    sqrt(#head)  ratio")
    for r in (rows[0], rows[len(rows) // 2], rows[-1]):
        say("  %-12d %-9.4f %-12.4f %.4f"
            % (r[0], r[4], math.sqrt(r[2]), r[4] / math.sqrt(r[2])))
    hg = np.array([r[4] for r in rows])
    hs = np.array([math.sqrt(r[2]) for r in rows])
    er, rr2, sr = fit(x, np.log(hg / hs))
    say("  the ratio head G / sqrt(#head) has exponent %+.6f, "
        "s.e. %.6f, t = %.2f" % (er, sr, abs(er) / sr))
    say("TSTAT slope_splitreach_headref %.2f" % (abs(er) / sr))
    say("SPREAD slope_splitreach_headref %.4f"
        % float(x.max() - x.min()))
    if abs(er) / sr < 2.0:
        say("UNRESOLVED SIGN slope_splitreach_headref")
    say("REFERENCE audit_split_reach %d %.4f %.4f"
        % (len(rows), float((hg / hs).min()), float((hg / hs).max())))

    say()
    say("=" * 70)
    say("Z1 %s  Z2 %s  Z3 %s  Z4 %s  Z5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (z1, z2, z3, z4, z5)))

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) over the squarefree",
        "           k < N^theta' coprime to N, the gain",
        "           sum|a| / |sum a| on the whole range, on the top",
        "           tenth of k by |a_k| and on the bottom nine",
        "           tenths; that head's share of sum|a|, the fraction",
        "           of the head carrying one sign, and its overlap",
        "           with the smallest tenth of k by size; each of the",
        "           three gains' least-squares exponent in log N over",
        "           the on-field family to 1.024e8, the one-sign",
        "           fraction's slope and the mass share's exponent.",
        "NULL: none is run for the split, a deterministic partition",
        "      of a measured sequence, and none is needed for the",
        "      one-sign fraction: its null is the binomial one half",
        "      and Z4 uses it as the cap. The coin arms for the gain",
        "      itself are audit_crossk_reference.py and",
        "      lab_gain_opposition.py.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta'; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to 102400000; the head fraction, the",
        "       field, the sieve and theta' are imported from",
        "       code/audit_gain_split.py; the published rows and",
        "       exponents are read from results/audit_gain_split.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not z1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
