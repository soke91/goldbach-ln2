# -*- coding: utf-8 -*-
r"""
The seventy points are three fields, not one.

WHAT IS AT STAKE

audit_flatness_fill.py filled the family with every N = 2^a 5^b in
[2e5, 2.56e7] and its docstring says "Every N = 2^a 5^b has odd radical
5, so the admissible k-set and the threshold are fixed exactly as along
the doublings". That sentence is false at the two edges of the family.
The loop starts at a = 0 and at b = 0, so the seventy points contain

  * N = 2^a with b = 0, whose odd radical is 1, not 5. There k runs
    over the squarefree k coprime to 2 alone, and every k divisible by
    5 -- which the doublings exclude -- is admitted.
  * N = 5^b with a = 0, which is ODD. There k runs over the squarefree
    k coprime to 5 alone, and every EVEN k is admitted.

Three coprimality classes were averaged as one. The file's own RADICALS
line says 2 and the gate accepted it, because G34 counts odd radicals
and 2^a contributes the empty one, so nothing rang.

This matters because the fits are what the last cycle read. The gain
column of that file runs 1.72 to 3.76 on the mixed N and 3.79 to 12.70
on the ten that are not mixed, and the r.m.s. of the seventy-point fits
is large enough that the excess over the ceiling theta'/2, which stood
at 3.15 standard errors on eight doublings, was reported as 0.02 and
declared gone. The conclusion drawn from that -- that item 4(b)'s target
is now an arithmetic constant and F has stopped rising -- rests on a
mixture of three fields.

WHAT WAS ALREADY KNOWN WHEN THIS WAS WRITTEN

An unregistered scratch refit of that file's own printed table, split by
class, was run before this script existed. It is the reason this script
exists, and it fixed e(l1/l2), e(G) and e(#k) on each class. So V3 and
V4 below are CONFIRMATORY, not blind, and are marked as such. V1, V2,
V5 and V6 were not computed before this ran.

BACKS: Remark {#rem:fillfield} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The control, and the mechanism. Recomputing #k here from the
      squarefree k coprime to N reproduces every one of the seventy
      counts in that file exactly; and the ten off-field N are exactly
      the ones whose admissible k-set is not the coprime-to-10 one.
  V2  The ten are not a rescaling of the sixty. Every off-field N's
      gain exceeds the largest gain of any on-field N within a factor
      of two of it, both ways.
  V3  CONFIRMATORY. On the sixty on-field N, e(l1/l2) exceeds the
      arithmetic ceiling theta'/2 by more than two of its own standard
      errors -- the eight-point excess returns rather than dissolving.
  V4  CONFIRMATORY. And it exceeds the MEASURED ceiling e(#k)/2 by more
      than two standard errors of the difference, so F is resolved
      rising on the field, against that file's U4.
  V5  The two shapes of {#rem:flatnessshape} separate on the sixty:
      the r.m.s. gap between the power law and a + b/log N exceeds the
      r.m.s.'s own standard error.
  V6  F is smooth at sub-doubling spacing on the field: every octave's
      r.m.s. residual of log F over the on-field points of that octave
      is below the published eight-point r.m.s.

REFUTATION RULE (fixed before the run)

  V1  REFUTED on any departure in any count, or if the ten off-field N
      are not exactly the ones with a different coprimality class.
      Then this is not the partition being claimed and nothing below
      may be read. THIS ONE GATES.
  V2  REFUTED if any off-field gain sits inside the range of on-field
      gains in its own window. Then the ten are a loud sample of the
      same field rather than a different one, the partition is
      cosmetic, and the seventy-point fits stand as published.
  V3  REFUTED if the excess falls within two standard errors of zero
      or below. Then the dissolution reported on seventy points was not
      an artefact of the mixture and audit_flatness_fill's U3 stands.
  V4  REFUTED if the difference against the measured ceiling is within
      two standard errors. Then F's rise is unresolved on the field
      too, and that file's U4 stands for the right reason.
  V5  REFUTED if the gap is at or below the r.m.s.'s own standard
      error. Then removing the mixture does not decide the shape
      either, and {#rem:flatnessshape}'s tie survives on a clean field.
  V6  REFUTED if any octave's r.m.s. reaches the published eight-point
      value. Then F's scatter on the clean field is still its own
      noise, and no denser sampling of this range settles the shape.

  V1 gates. V2 to V6 are the measurement and do not gate.

  NO NULL IS RUN and none applies. The partition is arithmetic, not
  statistical: which primes divide N is not a measured quantity and has
  no background to detect against. The comparison in V2 is between
  measured gains at neighbouring N, and the arm that randomises signs
  cannot move any of these -- |a_k| and #k are properties of the
  magnitudes and the k-set alone, as audit_flatness_fill's own NULL
  says.
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
OUT = os.path.join(RES, "audit_fill_field.txt")

LO, HI = 200_000, 25_600_000
ROW = re.compile(r"^  (\d{5,})\s+([\d.]+)\s+(\d+)\s+([\d.]+)\s+"
                 r"([\d.]+)\s+([\d.]+)\s*$", re.M)


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")
THETA = SPL.THETA
RUNGS = SPL.NS


def read_filled():
    """the seventy rows and the published eight-point r.m.s."""
    src = io.open(os.path.join(RES, "audit_flatness_fill.txt"),
                  encoding="utf-8").read()
    rows = []
    for m in ROW.finditer(src):
        rows.append((int(m.group(1)), int(m.group(3)),
                     float(m.group(4)), float(m.group(5)),
                     float(m.group(6))))
    src2 = io.open(os.path.join(RES, "audit_flatness_shape.txt"),
                   encoding="utf-8").read()
    m = re.search(r"best r\.m\.s\. ([\d.]+), standard error ([\d.]+)",
                  src2)
    return rows, float(m.group(1)), float(m.group(2))


def squarefree_upto(n):
    sf = np.ones(n + 2, dtype=bool)
    sf[0] = False
    d = 2
    while d * d <= n:
        sf[d * d::d * d] = False
        d += 1
    return sf


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum())) \
        if x.size > 2 else float("inf")
    cor = abs(float(np.corrcoef(x, y)[0, 1])) if x.size > 2 else 0.0
    return float(a), float(np.sqrt((r ** 2).mean())), se, cor


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    rows, prms, pse = read_filled()
    ceil_ = THETA / 2.0
    say("read from results/audit_flatness_fill.txt: %d rows" % len(rows))
    say("  and the published eight-point r.m.s. %.6f with its own "
        "standard error %.6f from results/audit_flatness_shape.txt"
        % (prms, pse))
    say("  theta' and the k-set rule are imported from "
        "code/audit_gain_split.py")
    say("  the arithmetic ceiling on e(l1/l2) is theta'/2 = %.4f"
        % ceil_)

    NS = [r[0] for r in rows]
    nk = {r[0]: r[1] for r in rows}
    L12 = {r[0]: r[2] for r in rows}
    F = {r[0]: r[3] for r in rows}
    G = {r[0]: r[4] for r in rows}

    KMAX = int(HI ** THETA) + 2
    sf = squarefree_upto(KMAX)

    def kset(N, primes):
        K = int(N ** THETA)
        c = 0
        for k in range(2, K):
            if not sf[k]:
                continue
            if any(k % q == 0 for q in primes):
                continue
            c += 1
        return c

    cls = {}
    for N in NS:
        cls[N] = tuple(sorted(SPL.factor_set(N)))
    classes = sorted(set(cls.values()))
    onfield = [N for N in NS if cls[N] == (2, 5)]
    off = [N for N in NS if cls[N] != (2, 5)]

    say()
    say("the classes inside the family")
    say("  primes of N   how many   what k is coprime to   N even?")
    for c in classes:
        got = [N for N in NS if cls[N] == c]
        say("  %-13s %-10d %-22s %s"
            % ("{" + ",".join(str(q) for q in c) + "}", len(got),
               "coprime to " + "*".join(str(q) for q in c),
               "yes" if 2 in c else "NO"))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))
    say("COPRIME FOR audit_flatness_fill %d" % len(classes))
    if len(classes) > 1:
        say("FIELD SPLIT audit_fill_field")
        say("FIELD SPLIT audit_flatness_fill")

    # -------------------------------------------------------------- V1
    say()
    say("V1  the control: is #k what the coprimality class says?")
    say("  N            #k read   #k here   coprime to")
    worst = 0
    bad_class = []
    for N in NS:
        c = kset(N, cls[N])
        worst = max(worst, abs(c - nk[N]))
        if N in off or N in onfield[:0]:
            say("  %-12d %-9d %-9d %s"
                % (N, nk[N], c, "*".join(str(q) for q in cls[N])))
    for N in RUNGS:
        if N in nk:
            c = kset(N, cls[N])
            worst = max(worst, abs(c - nk[N]))
    ten = sorted(off)
    v1a = worst == 0
    v1b = all(cls[N] != (2, 5) for N in ten) and len(ten) == len(off)
    v1 = v1a and v1b
    say("  worst departure over all %d rows: %d" % (len(NS), worst))
    say("  the off-field N are %d of them and none has the "
        "coprime-to-10 k-set: %s" % (len(ten), "yes" if v1b else "no"))
    say("  V1 %s   (cap 0 on any count)" % ("hold" if v1 else "REFUTED"))

    # -------------------------------------------------------------- V2
    say()
    say("V2  are the ten a rescaling of the sixty, or another field?")
    say("  N            class    G         on-field G in [N/2, 2N]  "
        "inside?")
    v2 = True
    for N in ten:
        win = [G[M] for M in onfield if N / 2.0 <= M <= 2.0 * N]
        if not win:
            continue
        inside = min(win) <= G[N] <= max(win)
        if inside:
            v2 = False
        say("  %-12d %-8s %-9.4f %.4f to %-16.4f %s"
            % (N, "*".join(str(q) for q in cls[N]), G[N],
               min(win), max(win), "INSIDE" if inside else "outside"))
    say("  on-field gain runs %.4f to %.4f over the whole family"
        % (min(G[N] for N in onfield), max(G[N] for N in onfield)))
    say("  off-field gain runs %.4f to %.4f"
        % (min(G[N] for N in off), max(G[N] for N in off)))
    say("REFERENCE audit_fill_field %d %.4f %.4f"
        % (len(onfield), min(F[N] for N in onfield),
           max(F[N] for N in onfield)))
    say("  V2 %s   (cap: any one inside its own window)"
        % ("hold" if v2 else "REFUTED"))

    # ------------------------------------------------- the three fits
    xall = np.log(np.array(NS, dtype=np.float64))
    xon = np.log(np.array(onfield, dtype=np.float64))
    say()
    say("  the same three exponents, per class and on the mixture")
    say("  set             n    e(#k)        e(l1/l2)     e(G)        "
        " r.m.s. of log G")
    fits = {}
    for c in classes + [None]:
        got = NS if c is None else [N for N in NS if cls[N] == c]
        if len(got) < 3:
            say("  %-15s %-4d too few points to fit"
                % ("{" + ",".join(str(q) for q in c) + "}", len(got)))
            continue
        xx = np.log(np.array(got, dtype=np.float64))
        ek, _rk, sk, _ck = fit(xx, np.log(np.array(
            [float(nk[N]) for N in got])))
        ec, _rc, sc, _cc = fit(xx, np.log(np.array(
            [L12[N] for N in got])))
        eg, rg, sg, _cg = fit(xx, np.log(np.array([G[N] for N in got])))
        fits[c] = (len(got), ek, sk, ec, sc, eg, sg, rg)
        say("  %-15s %-4d %+-12.6f %+-12.6f %+-12.6f %.6f"
            % ("the mixture" if c is None
               else "{" + ",".join(str(q) for q in c) + "}",
               len(got), ek, ec, eg, rg))

    # -------------------------------------------------------------- V3
    n_on, ek, sk, ec, sc, eg, sg, rg = fits[(2, 5)]
    say()
    say("V3  does the excess over the ceiling return on the field?")
    say("  on the %d on-field N, e(l1/l2) = %+.6f, s.e. %.6f"
        % (n_on, ec, sc))
    say("  the excess over theta'/2 is %+.6f, which is %.2f standard "
        "errors" % (ec - ceil_, (ec - ceil_) / sc))
    _emix, _r, smix, _c = fit(xall, np.log(np.array(
        [L12[N] for N in NS])))
    emix = fits[None][3]
    say("  on the mixture it was %+.6f, s.e. %.6f, an excess of %+.6f "
        "= %.2f standard errors"
        % (emix, smix, emix - ceil_, (emix - ceil_) / smix))
    v3 = (ec - ceil_) > 2.0 * sc
    say("TSTAT slope_audit_fill_field %.2f" % (abs(ec) / sc))
    say("SPREAD slope_audit_fill_field %.4f"
        % float(xon.max() - xon.min()))
    if abs(ec) / sc < 2.0:
        say("UNRESOLVED SIGN slope_audit_fill_field")
    say("  V3 %s   (cap 2 standard errors)"
        % ("hold" if v3 else "REFUTED"))

    # -------------------------------------------------------------- V4
    say()
    say("V4  and against the measured ceiling e(#k)/2?")
    ef, rf, sfe, _cf = fit(xon, np.log(np.array([F[N] for N in onfield])))
    say("  e(#k) on the field %+.6f, s.e. %.6f, so the measured "
        "ceiling is %+.6f" % (ek, sk, ek / 2.0))
    say("  e(l1/l2) - e(#k)/2 = %+.6f" % (ec - ek / 2.0))
    say("  and the same difference fitted directly as e(F): %+.6f, "
        "s.e. %.6f, t = %.2f" % (ef, sfe, abs(ef) / sfe))
    v4 = ef > 0.0 and abs(ef) / sfe >= 2.0
    say("TSTAT slope_fillfield_flat %.2f" % (abs(ef) / sfe))
    say("SPREAD slope_fillfield_flat %.4f"
        % float(xon.max() - xon.min()))
    if abs(ef) / sfe < 2.0:
        say("UNRESOLVED SIGN slope_fillfield_flat")
    emixF, _rmF, smixF, _cmF = fit(xall, np.log(np.array(
        [F[N] for N in NS])))
    say("  on the mixture the same slope was %+.6f, s.e. %.6f, t = "
        "%.2f, which is what that file read as unresolved"
        % (emixF, smixF, abs(emixF) / smixF))
    say("TSTAT slope_fillfield_mixture %.2f" % (abs(emixF) / smixF))
    say("SPREAD slope_fillfield_mixture %.4f"
        % float(xall.max() - xall.min()))
    if abs(emixF) / smixF < 2.0:
        say("UNRESOLVED SIGN slope_fillfield_mixture")
    say("  V4 %s   (cap 2 standard errors)"
        % ("hold" if v4 else "REFUTED"))

    # -------------------------------------------------------------- V5
    say()
    say("V5  do the two shapes separate once the mixture is out?")
    yF = np.array([F[N] for N in onfield])
    A1 = np.column_stack([np.ones_like(xon), xon])
    c1, *_ = np.linalg.lstsq(A1, np.log(yF), rcond=None)
    r1abs = float(np.sqrt(((yF - np.exp(A1.dot(c1))) ** 2).mean()))
    A2 = np.column_stack([np.ones_like(xon), 1.0 / xon])
    c2, *_ = np.linalg.lstsq(A2, yF, rcond=None)
    r2 = float(np.sqrt(((yF - A2.dot(c2)) ** 2).mean()))
    se_r = min(r1abs, r2) / math.sqrt(2.0 * (len(onfield) - 2))
    say("  shape                 r.m.s.")
    say("  F ~ N^e               %.6f   (e = %+.6f)" % (r1abs, c1[1]))
    say("  F = a + b/log N       %.6f   (a = %.6f, b = %.6f)"
        % (r2, c2[0], c2[1]))
    say("  the gap %.6f against the r.m.s.'s own standard error %.6f "
        "on %d points" % (abs(r1abs - r2), se_r, len(onfield)))
    v5 = abs(r1abs - r2) > se_r
    say("  V5 %s   (cap: the r.m.s.'s own standard error)"
        % ("hold" if v5 else "REFUTED"))
    if not v5:
        say("  so a clean field does not decide the shape either -- the "
            "tie is not the mixture's doing")

    # -------------------------------------------------------------- V6
    nocts = 0
    while LO * 2 ** (nocts + 1) <= HI:
        nocts += 1
    octs = []
    for j in range(nocts):
        lo, hi = LO * 2 ** j, LO * 2 ** (j + 1)
        got = [N for N in onfield if lo <= N < hi]
        if j == nocts - 1:
            got = [N for N in onfield if lo <= N <= hi]
        if len(got) >= 3:
            octs.append((j, got))
    say()
    say("V6  is F smooth at sub-doubling spacing on the field?")
    say("  octave   points   slope        s.e.        r.m.s.     |corr|")
    v6 = True
    locs, cors = [], []
    for jj, got in octs:
        xo = np.log(np.array(got, dtype=np.float64))
        yo = np.log(np.array([F[N] for N in got]))
        a_, r_, s_, c_ = fit(xo, yo)
        locs.append(a_)
        cors.append(c_)
        if r_ >= prms:
            v6 = False
        say("  %-8d %-8d %+-12.6f %-11.6f %-10.6f %.5f"
            % (jj, xo.size, a_, s_, r_, c_))
    say("  against the published eight-point r.m.s. %.6f" % prms)
    say("SWEPT fillfield_localslope octave-range %.6f"
        % (max(locs) - min(locs)))
    say("POP fillfield_localslope %d" % min(len(g) for _j, g in octs))
    say("CORR fillfield_localslope %.5f" % min(cors))
    say("  V6 %s   (cap: the published eight-point r.m.s.)"
        % ("hold" if v6 else "REFUTED"))

    # ------------------------------------------------- what is left
    say()
    say("  and the deficit item 4(b) has to close, on the field alone")
    say("  quantity                       value        s.e.")
    say("  e(l1/l2)                       %+-12.6f %.6f" % (ec, sc))
    say("  e(G)                           %+-12.6f %.6f" % (eg, sg))
    say("  the deficit e(l1/l2) - e(G)    %+-12.6f" % (ec - eg))
    say("  the distance to the arithmetic ceiling, theta'/2 - e(G)")
    say("                                 %+-12.6f" % (ceil_ - eg))
    say("TSTAT slope_fillfield_gain %.2f" % (abs(eg) / sg))
    say("SPREAD slope_fillfield_gain %.4f"
        % float(xon.max() - xon.min()))
    if abs(eg) / sg < 2.0:
        say("UNRESOLVED SIGN slope_fillfield_gain")
    say("  so the target is the measured e(l1/l2) and not the ceiling: "
        "the")
    say("  concentration is still above theta'/2 on the field, and "
        "what")
    say("  closes that gap is F reaching its bound, not a fit.")

    say()
    say("=" * 70)
    say("V1 %s  V2 %s  V3 %s  V4 %s  V5 %s  V6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (v1, v2, v3, v4, v5, v6)))

    head = [
        "STATISTIC: the seventy rows of results/audit_flatness_fill.txt",
        "           partitioned by the set of primes dividing N, which",
        "           fixes the admissible k-set; the count #k",
        "           recomputed here from the squarefree k coprime to N",
        "           against the count that file printed; the gain of",
        "           each off-field N against the range of on-field",
        "           gains within a factor of two of it; and the",
        "           least-squares exponents of #k, l1/l2, F and G in",
        "           log N on each class separately and on the mixture,",
        "           with the two shapes of {#rem:flatnessshape} and the",
        "           octave-wise r.m.s. of log F refitted on the",
        "           on-field points alone.",
        "NULL: none is run and none applies. The partition is",
        "      arithmetic -- which primes divide N is not a measured",
        "      quantity and has no background to detect against -- and",
        "      the quantities compared are the ones",
        "      audit_flatness_fill.py already declared no sign arm can",
        "      move, being properties of |a_k| and of the k-set alone.",
        "FIELD: the N are read from results/audit_flatness_fill.txt,",
        "       that is every N = 2^a 5^b in [2e5, 2.56e7] with a and b",
        "       from zero, so the family carries three coprimality",
        "       classes and not one, as COPRIME and FIELD SPLIT say:",
        "       {2,5} with k coprime to 10 and N even, which is the",
        "       field of the doublings; {2} with k coprime to 2 alone;",
        "       and {5} with k coprime to 5 alone and N odd. k is",
        "       squarefree with 2 <= k < N^theta' throughout; theta',",
        "       factor_set and the k-set rule are imported from",
        "       code/audit_gain_split.py; the published eight-point",
        "       r.m.s. is read from results/audit_flatness_shape.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not v1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
