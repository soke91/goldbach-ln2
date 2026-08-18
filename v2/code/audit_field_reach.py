# -*- coding: utf-8 -*-
r"""
Two more octaves of the field, which is what V5 said it would take.

WHAT IS AT STAKE

Remark {#rem:fillfield} cleaned the family: of the 70 points
{#rem:flatnessfill} fitted, ten were not the field, and on the 60 that
were, the concentration exponent sits 4.98 standard errors above the
ceiling theta'/2 and F is resolved rising at t = 4.76. Two of its
pre-registered predictions were refuted, and both said the same thing.
V5: the power law and the bounded shape a + b/log N are still tied on a
clean field, the gap being 0.000041 against the r.m.s.'s own standard
error 0.000602. V6: F's octave-wise scatter is at or above the
published eight-point r.m.s. in four octaves of seven. Both readings
end in the same sentence -- more points inside [2e5, 2.56e7] cannot
decide where F saturates, only a longer lever can.

The lever is affordable. The whole field to 2.56e7 costs seconds per
N, and the sieve to 1.024e8 fits in memory, so the family can be run
two octaves further: [2e5, 1.024e8] carries 80 on-field N against the
60 already measured, and the spread in log N goes from 4.8520 to
6.2384. This is the same move {#rem:slopereach} made when a family's
level slope was unresolved at five points -- not a new argument, a
longer sweep.

Nothing here is known in advance. No point above 2.56e7 has been
computed in this repository, on this statistic or any other.

BACKS: Remark {#rem:fieldreach} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The control. At every on-field N that
      results/audit_fill_field.txt carries, #k reproduces exactly and
      l1/l2, F and G reproduce to the precision those tables print.
  W2  The excess over the ceiling survives the longer lever:
      e(l1/l2) - theta'/2 stays above two of its own standard errors.
  W3  F is still resolved rising: e(F) positive at two standard errors
      over the extended field.
  W4  The two shapes separate. On the extended field the r.m.s. gap
      between F ~ N^e and F = a + b/log N exceeds the r.m.s.'s own
      standard error -- the thing V5 said only a longer lever could do.
  W5  The deficit is unchanged: e(l1/l2) - e(G) stays within two
      standard errors of the +0.137891 measured on the 60.
  W6  And the rise has not stopped at the new top: the local slope of
      log F over the topmost octave is positive at two standard errors.

REFUTATION RULE (fixed before the run)

  W1  REFUTED on any departure in a count, or above 0.0001 on a
      printed ratio. Then this is not the field those tables measured
      and nothing below may be read. THIS ONE GATES.
  W2  REFUTED if the excess falls to two standard errors or less. Then
      the concentration exponent is at its ceiling after all and
      {#rem:fillfield}'s 4.98 was the short lever, exactly as
      {#rem:fillfield} said the mixture's 0.02 was the wrong field.
  W3  REFUTED if e(F) is not resolved positive. Then F has stopped
      rising somewhere inside the new range and the saturation
      {#rem:flatnessshape} predicted has been reached, not forecast.
  W4  REFUTED if the gap is at or below the r.m.s.'s own standard
      error. Then two octaves are not enough lever either, and the
      shape question is not answerable at any reach this programme can
      compute -- which retires it rather than leaving it open.
  W5  REFUTED if the deficit moves by more than two standard errors.
      Then the quantity item 4(b) has to close is not stable under the
      reach it is measured at, and no statement about closing it is
      safe.
  W6  REFUTED if the top octave's local slope is not resolved
      positive. Unlike W3 this is local: F may be rising over the
      whole sweep while already flat at the top, and that is where
      saturation would first show.

  W1 gates. W2 to W6 are the measurement and do not gate.

  NO NULL IS RUN and none applies, for the reason
  audit_flatness_fill.py gave: F and l1/l2 are deterministic ratios of
  norms of the measured vector, random signs do not move |a_k|, and
  the comparisons here are between the same statistic at different N.
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
OUT = os.path.join(RES, "audit_field_reach.txt")

LO, HI = 200_000, 102_400_000
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


def read_published():
    """the on-field rows and the deficit {#rem:fillfield} measured"""
    src = io.open(os.path.join(RES, "audit_flatness_fill.txt"),
                  encoding="utf-8").read()
    rows, dec = {}, 0
    for m in ROW.finditer(src):
        rows[int(m.group(1))] = (int(m.group(3)), float(m.group(4)),
                                 float(m.group(5)), float(m.group(6)))
        for g in (4, 5, 6):
            dec = max(dec, len(m.group(g).split(".")[1]))
    src2 = io.open(os.path.join(RES, "audit_fill_field.txt"),
                   encoding="utf-8").read()
    m = re.search(r"the deficit e\(l1/l2\) - e\(G\)\s+([+-][\d.]+)",
                  src2)
    m2 = re.search(r"e\(l1/l2\)\s+([+-][\d.]+)\s+([\d.]+)\s*$", src2,
                   re.M)
    m3 = re.search(r"e\(G\)\s+([+-][\d.]+)\s+([\d.]+)\s*$", src2, re.M)
    src3 = io.open(os.path.join(RES, "audit_flatness_shape.txt"),
                   encoding="utf-8").read()
    m4 = re.search(r"best r\.m\.s\. ([\d.]+), standard error ([\d.]+)",
                   src3)
    return (rows, float(m.group(1)), float(m2.group(1)),
            float(m2.group(2)), float(m3.group(1)), float(m3.group(2)),
            float(m4.group(1)), dec)


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
    cor = abs(float(np.corrcoef(x, y)[0, 1])) if x.size > 2 else 0.0
    return float(a), float(np.sqrt((r ** 2).mean())), se, cor


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    (pub, pubdef, pubec, pubsc, pubeg, pubsg,
     prms, dec) = read_published()
    ceil_ = THETA / 2.0
    NS = family(LO, HI)
    old = [N for N in NS if N <= 25_600_000]
    new = [N for N in NS if N > 25_600_000]

    say("read from results/audit_flatness_fill.txt: %d rows to "
        "compare against" % len(pub))
    say("  and from results/audit_fill_field.txt the on-field "
        "deficit %+.6f" % pubdef)
    say("  with e(l1/l2) %+.6f (s.e. %.6f) and e(G) %+.6f (s.e. %.6f)"
        % (pubec, pubsc, pubeg, pubsg))
    say("  and the published eight-point r.m.s. %.6f from "
        "results/audit_flatness_shape.txt" % prms)
    say("  theta' and the field are imported from "
        "code/audit_gain_split.py")
    say("  the arithmetic ceiling on e(l1/l2) is theta'/2 = %.4f"
        % ceil_)
    say()
    say("the field, extended: every N = 2^a 5^b with a and b at least "
        "one,")
    say("  in [%d, %d]; %d of them, %d of which are new above "
        "the old top" % (LO, HI, len(NS), len(new)))
    classes = sorted(set(tuple(sorted(SPL.factor_set(N))) for N in NS))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))
    say("  every N carries the prime set %s, so the k-set rule is one"
        % ("{" + ",".join(str(q) for q in classes[0]) + "}"))

    say()
    say("sieving to %d ..." % HI)
    lam, mu = SPL.lambda_and_mu(HI)
    sqf = mu != 0

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
    say("REFERENCE audit_field_reach %d %.4f %.4f"
        % (len(NS), min(F.values()), max(F.values())))

    x = np.log(np.array(NS, dtype=np.float64))
    xold = np.log(np.array(old, dtype=np.float64))

    # -------------------------------------------------------------- W1
    say()
    say("W1  the control against the rows already published")
    worstk, worstr = 0, 0.0
    seen = 0
    for N in old:
        if N not in pub:
            continue
        seen += 1
        pk, pl, pf, pg = pub[N]
        worstk = max(worstk, abs(nk[N] - pk))
        worstr = max(worstr, abs(L12[N] - pl), abs(F[N] - pf),
                     abs(G[N] - pg))
    w1 = worstk == 0 and worstr <= 0.0001
    say("  %d rows compared; worst departure %d on a count and %.6f "
        "on a printed ratio" % (seen, worstk, worstr))
    say("  W1 %s   (cap 0 and tol 0.0001)" % ("hold" if w1 else
                                              "REFUTED"))
    say("  the table it is judged against prints %d decimals, so a "
        "single value is within %.8f of what produced it"
        % (dec, 0.5 * 10.0 ** (-dec)))
    say("PRINTBOUND audit_field_reach %d %.8f"
        % (dec, 0.5 * 10.0 ** (-dec)))

    # -------------------------------------------------------------- W2
    say()
    say("W2  does the excess over the ceiling survive the lever?")
    ec, rc, sc, _cc = fit(x, np.log(np.array([L12[N] for N in NS])))
    eco, _ro, sco, _co = fit(xold, np.log(np.array(
        [L12[N] for N in old])))
    say("  on %d points over a spread of %.4f in log N, e(l1/l2) = "
        "%+.6f, s.e. %.6f" % (len(NS), float(x.max() - x.min()), ec,
                              sc))
    say("  the excess over theta'/2 is %+.6f = %.2f standard errors"
        % (ec - ceil_, (ec - ceil_) / sc))
    say("  on the old reach alone, %d points over %.4f: %+.6f, s.e. "
        "%.6f, excess %.2f standard errors"
        % (len(old), float(xold.max() - xold.min()), eco, sco,
           (eco - ceil_) / sco))
    w2 = (ec - ceil_) > 2.0 * sc
    say("TSTAT slope_audit_field_reach %.2f" % (abs(ec) / sc))
    say("SPREAD slope_audit_field_reach %.4f"
        % float(x.max() - x.min()))
    if abs(ec) / sc < 2.0:
        say("UNRESOLVED SIGN slope_audit_field_reach")
    say("  W2 %s   (cap 2 standard errors)" % ("hold" if w2 else
                                               "REFUTED"))

    # -------------------------------------------------------------- W3
    say()
    say("W3  is F still resolved rising?")
    ef, rf, sf_, _cf = fit(x, np.log(np.array([F[N] for N in NS])))
    ek, _rk, sk, _ck = fit(x, np.log(np.array(
        [float(nk[N]) for N in NS])))
    say("  e(#k) %+.6f, s.e. %.6f, so the measured ceiling is %+.6f"
        % (ek, sk, ek / 2.0))
    say("  e(F) %+.6f, s.e. %.6f, t = %.2f" % (ef, sf_, abs(ef) / sf_))
    say("  F runs %.4f to %.4f over the extended field"
        % (min(F.values()), max(F.values())))
    w3 = ef > 0.0 and abs(ef) / sf_ >= 2.0
    say("TSTAT slope_fieldreach_flat %.2f" % (abs(ef) / sf_))
    say("SPREAD slope_fieldreach_flat %.4f" % float(x.max() - x.min()))
    if abs(ef) / sf_ < 2.0:
        say("UNRESOLVED SIGN slope_fieldreach_flat")
    say("  W3 %s   (cap 2 standard errors)" % ("hold" if w3 else
                                               "REFUTED"))

    # -------------------------------------------------------------- W4
    say()
    say("W4  do the two shapes separate with the longer lever?")
    yF = np.array([F[N] for N in NS])
    A1 = np.column_stack([np.ones_like(x), x])
    c1, *_ = np.linalg.lstsq(A1, np.log(yF), rcond=None)
    r1abs = float(np.sqrt(((yF - np.exp(A1.dot(c1))) ** 2).mean()))
    A2 = np.column_stack([np.ones_like(x), 1.0 / x])
    c2, *_ = np.linalg.lstsq(A2, yF, rcond=None)
    r2 = float(np.sqrt(((yF - A2.dot(c2)) ** 2).mean()))
    se_r = min(r1abs, r2) / math.sqrt(2.0 * (len(NS) - 2))
    say("  shape                 r.m.s.")
    say("  F ~ N^e               %.6f   (e = %+.6f)" % (r1abs, c1[1]))
    say("  F = a + b/log N       %.6f   (a = %.6f, b = %.6f)"
        % (r2, c2[0], c2[1]))
    say("  the gap %.6f against the r.m.s.'s own standard error %.6f "
        "on %d points" % (abs(r1abs - r2), se_r, len(NS)))
    w4 = abs(r1abs - r2) > se_r
    say("SHAPES 2")
    say("  W4 %s   (cap: the r.m.s.'s own standard error)"
        % ("hold" if w4 else "REFUTED"))

    # -------------------------------------------------------------- W5
    say()
    say("W5  is the deficit stable under the reach?")
    eg, rg, sg, _cg = fit(x, np.log(np.array([G[N] for N in NS])))
    def_ = ec - eg
    sd = math.sqrt(sc * sc + sg * sg)
    say("  e(G) %+.6f, s.e. %.6f, against the published %+.6f (s.e. "
        "%.6f)" % (eg, sg, pubeg, pubsg))
    say("  the deficit e(l1/l2) - e(G) = %+.6f, against the published "
        "%+.6f" % (def_, pubdef))
    say("  the move is %+.6f against a combined standard error of "
        "%.6f, which is %.2f of it"
        % (def_ - pubdef, sd, abs(def_ - pubdef) / sd))
    w5 = abs(def_ - pubdef) <= 2.0 * sd
    say("TSTAT slope_fieldreach_gain %.2f" % (abs(eg) / sg))
    say("SPREAD slope_fieldreach_gain %.4f" % float(x.max() - x.min()))
    if abs(eg) / sg < 2.0:
        say("UNRESOLVED SIGN slope_fieldreach_gain")
    say("  W5 %s   (cap 2 standard errors of the move)"
        % ("hold" if w5 else "REFUTED"))

    # -------------------------------------------------------------- W6
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
    say("W6  and is the rise still going at the new top?")
    say("  octave   points   slope        s.e.        r.m.s.     |corr|")
    locs, cors = [], []
    for jj, got in octs:
        xo = np.log(np.array(got, dtype=np.float64))
        yo = np.log(np.array([F[N] for N in got]))
        a_, r_, s_, c_ = fit(xo, yo)
        locs.append(a_)
        cors.append(c_)
        say("  %-8d %-8d %+-12.6f %-11.6f %-10.6f %.5f"
            % (jj, xo.size, a_, s_, r_, c_))
    jt, top = octs[-1]
    xt = np.log(np.array(top, dtype=np.float64))
    yt = np.log(np.array([F[N] for N in top]))
    at, rt, st, _ct = fit(xt, yt)
    w6 = at > 0.0 and abs(at) / st >= 2.0
    say("  the top octave has %d points over %.4f in log N, slope "
        "%+.6f, s.e. %.6f, t = %.2f"
        % (xt.size, float(xt.max() - xt.min()), at, st, abs(at) / st))
    say("  and the octave r.m.s. against the published eight-point "
        "%.6f" % prms)
    say("TSTAT slope_fieldreach_top %.2f" % (abs(at) / st))
    say("SPREAD slope_fieldreach_top %.4f" % float(xt.max() - xt.min()))
    if abs(at) / st < 2.0:
        say("UNRESOLVED SIGN slope_fieldreach_top")
    say("SWEPT fieldreach_localslope octave-range %.6f"
        % (max(locs) - min(locs)))
    say("POP fieldreach_localslope %d" % min(len(g) for _j, g in octs))
    say("CORR fieldreach_localslope %.5f" % min(cors))
    say("  W6 %s   (cap 2 standard errors)" % ("hold" if w6 else
                                               "REFUTED"))

    say()
    say("=" * 70)
    say("W1 %s  W2 %s  W3 %s  W4 %s  W5 %s  W6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (w1, w2, w3, w4, w5, w6)))

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) over the squarefree",
        "           k < N^theta' coprime to N, the counts #k, the",
        "           concentration l1/l2, the flatness",
        "           F = (l1/l2)/sqrt(#k) and the gain l1/|sum a|, at",
        "           every on-field N = 2^a 5^b with a, b >= 1 in",
        "           [2e5, 1.024e8]; the least-squares exponents of",
        "           l1/l2, #k, F and G in log N over that whole set",
        "           and over the old reach alone; the two shapes of",
        "           {#rem:flatnessshape} refitted on it; and the",
        "           octave-wise local slope of log F including the two",
        "           octaves above the old top.",
        "NULL: none is run and none applies, for the reason",
        "      audit_flatness_fill.py gave: F and l1/l2 are",
        "      deterministic ratios of norms of the measured vector,",
        "      random signs do not move |a_k|, and every comparison",
        "      here is between the same statistic at different N.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], which is one coprimality class as",
        "       COPRIME says -- the class {2,5} of the doublings, k",
        "       coprime to 10 and N even; the degenerate members",
        "       {#rem:fillfield} removed, N = 2^a and N = 5^b, are not",
        "       generated here at all; k squarefree with",
        "       2 <= k < N^theta'; m over 1 <= m < N/k with (m,k) = 1;",
        "       Lambda and mu from an integer sieve to 102400000; the",
        "       field, the sieve and theta' are imported from",
        "       code/audit_gain_split.py; the rows compared against are",
        "       read from results/audit_flatness_fill.txt and the",
        "       published exponents from results/audit_fill_field.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not w1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
