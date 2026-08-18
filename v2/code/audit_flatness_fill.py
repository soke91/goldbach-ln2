# -*- coding: utf-8 -*-
r"""
Seventy points on the flatness, to find where F saturates.

WHAT IS AT STAKE

Remark {#rem:levelmagnitude} left item 4(b) with one computational axis.
The demand is e(G) -> e(l1/l2), and the concentration exponent is
+0.287798 with a standard error of 0.002472 -- 3.15 standard errors
ABOVE the ceiling theta'/2 = 0.28 that #k ~ N^theta' imposes. Remark
{#rem:flatnessshape} explained why that is possible without being a
contradiction: F = (l1/l2)/sqrt(#k) is 0.6622 to 0.6986 and still
rising, and Cauchy-Schwarz gives F <= 1, so the rise has to stop and
the exponent has to come down to the ceiling. Where it stops is not
known. Eight points cannot say: a power law beats the bounded shape
a + b/log N by 0.000163 in r.m.s. against an r.m.s. whose own standard
error is 0.001870, which is 0.09 of it, and the gate marks the two
shapes tied.

The eight points also wobble. F runs 0.6760, 0.6622, 0.6854, 0.6764,
0.6802, 0.6872, 0.6986, 0.6909 -- not monotone -- and the fits sit
0.006479 r.m.s. away from them. **Whether that 0.006479 is the
statistic's own noise or a shape sampled once per doubling is the thing
that decides whether the saturation question is answerable by
computing.** Remark {#rem:primorialdense} asked exactly this of the
primorial ladder's level exponent and got noise, which closed that
route. F is a different kind of object: a ratio of two norms over the
same k, with no location and nothing to quantise, so if any statistic
in this programme is smooth at sub-doubling spacing it is this one.

The family can be filled. Every N = 2^a 5^b has odd radical 5, so the
admissible k-set and the threshold are fixed exactly as along the
doublings, and there are 70 such N between 2e5 and 2.56e7 against the
8 that have been used.

BACKS: Remark {#rem:flatnessfill} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  The controls. At the eight doubling N, #k, l1/l2 and F reproduce
      results/audit_flatness_shape.txt to 0.001 and the gain
      reproduces results/audit_gain_split.txt to 0.01.
  U2  F is smooth at sub-doubling spacing. Fitted inside each octave
      on that octave's own points, the r.m.s. residual is below the
      published eight-point 0.006479 in every octave.
  U3  The excess over the ceiling survives the denser set. e(l1/l2) on
      the 70 points stays above theta'/2 by more than two of its own
      standard errors.
  U4  And F is still rising at the top. The local slope of log F over
      the top octave is positive at two standard errors.

REFUTATION RULE (fixed before the run)

  U1  REFUTED at 0.001 on any of the three or 0.01 on a gain. Any would
      mean this is not the field {#rem:flatnessshape} measured.
  U2  REFUTED if any octave's r.m.s. reaches the published 0.006479.
      Then F's scatter is its own noise, the shape question cannot be
      settled by more points in this range, and
      {#rem:flatnessshape}'s tie stands at this reach whatever else is
      computed.
  U3  REFUTED if the excess falls within two standard errors of zero
      or below. That is the outcome worth having: the 3.15 standard
      errors would be a small-sample artefact, the concentration
      exponent would sit at the ceiling where the arithmetic says it
      must, and {#rem:flatnessshape}'s puzzle would dissolve rather
      than being solved.
  U4  REFUTED if the top octave's local slope is not resolved
      positive. F may then already have saturated, and the ceiling
      would be being approached from below no longer.

  U1 gates: without it this is not that field.
  U2, U3 and U4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. F is a deterministic ratio of two
  norms of a measured vector and the comparisons are between its own
  values at different N; there is no background to detect against. The
  coin arm for this vector was run in audit_crossk_reference.py and the
  sign-vector floor in audit_lean_floor.py, and a coin's F is what the
  Cauchy-Schwarz bound is about: random signs do not change |a_k| at
  all, so F is a property of the magnitudes alone and no sign arm can
  move it.
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
OUT = os.path.join(RES, "audit_flatness_fill.txt")

LO, HI = 200_000, 25_600_000


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")
THETA = SPL.THETA
RUNGS = SPL.NS


def family(lo, hi):
    out = []
    a = 0
    while 2 ** a <= hi:
        b = 0
        while 2 ** a * 5 ** b <= hi:
            v = 2 ** a * 5 ** b
            if v >= lo:
                out.append(v)
            b += 1
        a += 1
    return sorted(set(out))


def read_published():
    src = io.open(os.path.join(RES, "audit_flatness_shape.txt"),
                  encoding="utf-8").read()
    rows = {}
    for m in re.finditer(r"^  N = (\d+)\s+#k = (\d+)\s+l1/l2 = "
                         r"([\d.]+)\s+F = ([\d.]+)\s*$", src, re.M):
        rows[int(m.group(1))] = (int(m.group(2)), float(m.group(3)),
                                 float(m.group(4)))
    m = re.search(r"best r\.m\.s\. ([\d.]+), standard error ([\d.]+)",
                  src)
    src2 = io.open(os.path.join(RES, "audit_gain_split.txt"),
                   encoding="utf-8").read()
    g = {}
    for mm in re.finditer(r"^  N = (\d+)\s+#k = \d+\s+head \d+\s+"
                          r"G ([\d.]+)\s", src2, re.M):
        g[int(mm.group(1))] = float(mm.group(2))
    m3 = re.search(r"^  whole\s+([+-][\d.]+)\s", src2, re.M)
    src3 = io.open(os.path.join(RES, "audit_lean_identity.txt"),
                   encoding="utf-8").read()
    m4 = re.search(r"l1/l2\s+([+-][\d.]+)", src3)
    return (rows, float(m.group(1)), float(m.group(2)), g,
            float(m3.group(1)), float(m4.group(1)))


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

    pub, prms, pse, pubg, pubeg, pubec = read_published()
    NS = family(LO, HI)
    say("read from results/audit_flatness_shape.txt: %d published rows, "
        "best" % len(pub))
    say("  r.m.s. %.6f with its own standard error %.6f"
        % (prms, pse))
    say("  and from results/audit_gain_split.txt %d published gains"
        % len(pubg))
    say("  the field and the sieve are imported from "
        "code/audit_gain_split.py")
    ceil_ = THETA / 2.0
    say("  the ceiling #k ~ N^theta' puts on e(l1/l2) is theta'/2 = "
        "%.4f" % ceil_)
    say()
    say("the family: every N = 2^a 5^b in [%d, %d], odd radical 5 "
        "throughout" % (LO, HI))
    say("  %d of them against the %d doublings that have been used"
        % (len(NS), len(RUNGS)))

    NMAX = max(NS)
    say()
    say("sieving to %d ..." % NMAX)
    lam, mu = SPL.lambda_and_mu(NMAX)
    sqf = mu != 0
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in SPL.factor_set(N) if q > 2))
                  for N in NS)))

    nk, L12, F, G = {}, {}, {}, {}
    say()
    say("  N            log10 N   #k      l1/l2      F         G")
    for N in NS:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        w = np.abs(a)
        l1 = float(w.sum())
        l2 = float(np.sqrt((a ** 2).sum()))
        s = abs(float(a.sum()))
        nk[N] = int(ks.size)
        L12[N] = l1 / l2
        F[N] = (l1 / l2) / math.sqrt(ks.size)
        G[N] = l1 / s if s > 0 else float("inf")
        say("  %-12d %-9.4f %-7d %-10.4f %-9.4f %.4f"
            % (N, math.log10(N), nk[N], L12[N], F[N], G[N]))
    say("REFERENCE audit_flatness_fill %d %.4f %.4f"
        % (len(NS), min(F.values()), max(F.values())))

    x = np.log(np.array(NS, dtype=np.float64))

    # ------------------------------------------------------------- U1
    say()
    say("U1  the controls at the eight doublings")
    say("  N            #k here  #k pub  l1/l2 here  pub        F here  "
        "  pub       G here   pub")
    worst, wg = 0.0, 0.0
    for N in RUNGS:
        if N not in pub:
            continue
        pk, pl, pf = pub[N]
        d = max(abs(nk[N] - pk), abs(L12[N] - pl), abs(F[N] - pf))
        worst = max(worst, d)
        if N in pubg:
            wg = max(wg, abs(G[N] - pubg[N]))
        say("  %-12d %-8d %-7d %-11.4f %-10.4f %-9.4f %-9.4f %-8.4f %.4f"
            % (N, nk[N], pk, L12[N], pl, F[N], pf, G[N],
               pubg.get(N, float("nan"))))
    u1 = worst < 0.001 and wg < 0.01
    say("  worst departure %.6f on the three, %.6f on a gain"
        % (worst, wg))
    say("  U1 %s   (cap 0.001 and 0.01)" % ("hold" if u1 else "REFUTED"))

    # ------------------------------------------------------------- U2
    nocts = 0
    while LO * 2 ** (nocts + 1) <= HI:
        nocts += 1
    octs = []
    for j in range(nocts):
        lo, hi = LO * 2 ** j, LO * 2 ** (j + 1)
        got = [N for N in NS if lo <= N < hi]
        if j == nocts - 1:
            got = [N for N in NS if lo <= N <= hi]
        if len(got) >= 3:
            octs.append((j, got))
    say()
    say("U2  is F smooth at sub-doubling spacing?")
    say("  octave   points   slope        s.e.        r.m.s.     |corr|")
    u2 = True
    locs, rmss, cors = [], [], []
    for jj, got in octs:
        xo = np.log(np.array(got, dtype=np.float64))
        yo = np.log(np.array([F[N] for N in got]))
        a_, r_, s_, c_ = fit(xo, yo)
        locs.append(a_)
        rmss.append(r_)
        cors.append(c_)
        if r_ >= prms:
            u2 = False
        say("  %-8d %-8d %+-12.6f %-11.6f %-10.6f %.5f"
            % (jj, xo.size, a_, s_, r_, c_))
    say("  the published eight-point r.m.s. is %.6f" % prms)
    say("  U2 %s" % ("hold" if u2 else "REFUTED"))
    say("SWEPT flatness_fill_localslope octave-range %.6f"
        % (max(locs) - min(locs)))
    say("POP flatness_fill_localslope %d"
        % min(len(g) for _j, g in octs))
    say("CORR flatness_fill_localslope %.5f" % min(cors))

    # ------------------------------------------------------------- U3
    say()
    say("U3  does the excess over the ceiling survive?")
    ec, rc, sc, _cc = fit(x, np.log(np.array([L12[N] for N in NS])))
    u3 = (ec - ceil_) > 2.0 * sc
    say("  e(l1/l2) on %d points %+.6f, s.e. %.6f, r.m.s. %.6f"
        % (len(NS), ec, sc, rc))
    say("  the published eight-point value was +0.287798 and the "
        "ceiling is %.4f" % ceil_)
    say("  the excess is %+.6f, which is %.2f standard errors"
        % (ec - ceil_, (ec - ceil_) / sc))
    say("TSTAT slope_flatnessfill_conc %.2f" % (abs(ec) / sc))
    say("SPREAD slope_flatnessfill_conc %.4f" % float(x.max() - x.min()))
    if abs(ec) / sc < 2.0:
        say("UNRESOLVED SIGN slope_flatnessfill_conc")
    say("  U3 %s   (cap 2 standard errors)"
        % ("hold" if u3 else "REFUTED"))

    # ------------------------------------------------------------- U4
    say()
    say("U4  is F still rising at the top?")
    jt, top = octs[-1]
    xt = np.log(np.array(top, dtype=np.float64))
    yt = np.log(np.array([F[N] for N in top]))
    at, rt, st, ct = fit(xt, yt)
    u4 = at > 0.0 and abs(at) / st >= 2.0
    say("  the top octave has %d points over %.4f in log N"
        % (xt.size, float(xt.max() - xt.min())))
    say("  its slope %+.6f, s.e. %.6f, t = %.2f" % (at, st, abs(at) / st))
    say("TSTAT slope_flatnessfill_top %.2f" % (abs(at) / st))
    say("SPREAD slope_flatnessfill_top %.4f"
        % float(xt.max() - xt.min()))
    if abs(at) / st < 2.0:
        say("UNRESOLVED SIGN slope_flatnessfill_top")
    ea, ra, sa, _ca = fit(x, np.log(np.array([F[N] for N in NS])))
    say("  over the whole sweep F rises at %+.6f, s.e. %.6f, t = %.2f"
        % (ea, sa, abs(ea) / sa))
    say("TSTAT slope_flatnessfill_all %.2f" % (abs(ea) / sa))
    say("SPREAD slope_flatnessfill_all %.4f" % float(x.max() - x.min()))
    if abs(ea) / sa < 2.0:
        say("UNRESOLVED SIGN slope_flatnessfill_all")
    say("  U4 %s" % ("hold" if u4 else "REFUTED"))

    # the two shapes, on 70 points
    say()
    say("  and the two shapes {#rem:flatnessshape} could not separate, "
        "refitted:")
    yF = np.array([F[N] for N in NS])
    A1 = np.column_stack([np.ones_like(x), x])
    c1, *_ = np.linalg.lstsq(A1, np.log(yF), rcond=None)
    r1 = float(np.sqrt(((np.log(yF) - A1.dot(c1)) ** 2).mean()))
    A2 = np.column_stack([np.ones_like(x), 1.0 / x])
    c2, *_ = np.linalg.lstsq(A2, yF, rcond=None)
    r2 = float(np.sqrt(((yF - A2.dot(c2)) ** 2).mean()))
    r1abs = float(np.sqrt(((yF - np.exp(A1.dot(c1))) ** 2).mean()))
    se_r = min(r1abs, r2) / math.sqrt(2.0 * (len(NS) - 2))
    say("  shape                 r.m.s.")
    say("  F ~ N^e               %.6f   (e = %+.6f)" % (r1abs, c1[1]))
    say("  F = a + b/log N       %.6f   (a = %.6f, b = %.6f)"
        % (r2, c2[0], c2[1]))
    say("  the gap %.6f against the r.m.s.'s own standard error %.6f "
        "on %d points" % (abs(r1abs - r2), se_r, len(NS)))
    say("  on eight points the same gap stood against a published "
        "standard")
    say("  error of %.6f (published), which is what left the two tied"
        % pse)
    if abs(r1abs - r2) > se_r:
        say("  so the denser set SEPARATES them, which eight points "
            "could not")
    else:
        say("  so the denser set does not separate them either: more "
            "points in")
        say("  the same range do not decide a shape, only a longer "
            "lever would")

    # ------------------------------------------- not pre-registered
    say()
    say("X1  and the other side of the deficit, on the same 70 points")
    say("  (written after U3 fell; not pre-registered). The deficit")
    say("  {#rem:leanidentity} asks to close is e(l1/l2) - e(G), and "
        "only")
    say("  one side of it was refitted above. Both sides have to be on "
        "the")
    say("  same points or the difference is not a difference.")
    eg, rg, sg, _cg = fit(x, np.log(np.array([G[N] for N in NS])))
    say("  quantity     exponent     s.e.       r.m.s.     published")
    say("  e(G)         %+-12.6f %-10.6f %-10.6f %+.6f"
        % (eg, sg, rg, pubeg))
    say("  e(l1/l2)     %+-12.6f %-10.6f %-10.6f %+.6f"
        % (ec, sc, rc, pubec))
    say("  the deficit on 70 points is %+.6f against the eight-point "
        "%+.6f" % (ec - eg, pubec - pubeg))
    say("TSTAT slope_flatnessfill_gain %.2f" % (abs(eg) / sg))
    say("SPREAD slope_flatnessfill_gain %.4f" % float(x.max() - x.min()))
    if abs(eg) / sg < 2.0:
        say("UNRESOLVED SIGN slope_flatnessfill_gain")
    say("  and whether the gain is smooth at sub-doubling spacing, "
        "which is")
    say("  the same question U2 asked of F:")
    say("  octave   points   r.m.s. of log G")
    for jj, got in octs:
        xo = np.log(np.array(got, dtype=np.float64))
        yo = np.log(np.array([G[N] for N in got]))
        _a2, r2_, _s2, _c2 = fit(xo, yo)
        say("  %-8d %-8d %.6f" % (jj, xo.size, r2_))
    say("  so the target of item 4(b) is now an arithmetic constant "
        "and not")
    say("  a fitted one: e(G) has to reach theta'/2 = %.4f, and the "
        "distance" % ceil_)
    say("  it has to travel is %+.6f." % (ceil_ - eg))

    say()
    say("=" * 70)
    say("U1 %s  U2 %s  U3 %s  U4 %s"
        % tuple("hold" if v else "REFUTED" for v in (u1, u2, u3, u4)))

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) over the squarefree",
        "           k < N^theta' coprime to N, the counts #k, the",
        "           concentration l1/l2, the flatness",
        "           F = (l1/l2)/sqrt(#k) and the gain l1/|sum a|, at",
        "           every N = 2^a 5^b in [2e5, 2.56e7]; the",
        "           least-squares slope of log F inside each octave of N",
        "           with its r.m.s. residual against the published",
        "           eight-point r.m.s.; the exponent of l1/l2 on the",
        "           whole set against the ceiling theta'/2; the top",
        "           octave's local slope; and the two shapes of",
        "           {#rem:flatnessshape} refitted on all of them.",
        "NULL: none is run and none applies. F is a deterministic ratio",
        "      of two norms of a measured vector and the comparisons are",
        "      between its own values at different N; there is no",
        "      background to detect against. Random signs do not change",
        "      |a_k|, so F is a property of the magnitudes alone and no",
        "      sign arm can move it -- the coin arm of",
        "      audit_crossk_reference.py and the sign-vector floor of",
        "      audit_lean_floor.py are the arms for the gain, not for",
        "      this.",
        "FIELD: N = 2^a 5^b in [2e5, 2.56e7], every one with odd radical",
        "       5, as RADICALS says, so the admissible k-set and the",
        "       threshold are fixed exactly as along the doublings; k",
        "       squarefree and coprime to N with 2 <= k < N^theta';",
        "       m over 1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + "; the field, the sieve",
        "       and theta' are imported from code/audit_gain_split.py;",
        "       the published rows are read from",
        "       results/audit_flatness_shape.txt and the published gains",
        "       from results/audit_gain_split.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not u1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
