# -*- coding: utf-8 -*-
r"""
Is alpha a constant, or is the deficit closing on its own?

WHAT IS AT STAKE

Six cycles have quoted alpha = +0.717916 as if it were a number the
field has. {#rem:headfraction} refitted the top octave alone and got
+0.664374 for the same quantity, and declined to read it because a
twelve-point fit came with no error printed. That is the right refusal
and the wrong place to stop: if alpha drifts downward, the distance
item 4(b) has to travel is shrinking on its own, and every statement
made about "+0.134019" is a statement about one window.

The quantity that matters is the deficit alpha - e(l2), not alpha:
{#rem:denominator}'s identity makes the demand |sum a| ~ l2, so a
drift in alpha buys nothing if e(l2) drifts with it. Both are
measurable octave by octave, with the errors the earlier refusal
lacked.

Three outcomes are worth separating. The local exponents may be flat,
in which case the top-octave value was noise and the global fit
stands. They may drift with the deficit closing, which would be the
first quantity in this programme moving towards the demand rather
than away. Or alpha may drift with e(l2) drifting equally, in which
case nothing has changed and the drift is a property of the fit rather
than of the problem.

BACKS: Remark {#rem:alphalocal} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  The control. The whole-field alpha reproduces
      results/audit_denominator.txt inside the bound its printing
      forces.
  Z2  The top octave is genuinely lower: its local exponent for
      |sum a| sits below the whole-field alpha by more than two of
      its own standard errors.
  Z3  And it is a drift, not one odd octave: the local exponents'
      least-squares slope against log N is resolved negative.
  Z4  But the deficit does not close: alpha_local - e(l2)_local stays
      positive by more than two standard errors in every octave.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED outside the printing bound. THIS ONE GATES.
  Z2  REFUTED if the top octave is within two standard errors of the
      global alpha. Then {#rem:headfraction} was right to refuse and
      the top-octave numbers it printed are noise.
  Z3  REFUTED if the slope is not resolved negative. Then whatever
      the top octave is doing is local to it, and alpha may be
      treated as a constant of the range after all.
  Z4  REFUTED if the deficit closes in any octave. That is the
      outcome worth having and the one this programme has not seen
      once: the demand met, at some scale, by the field itself.

  Z1 gates. Z2 to Z4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. Every quantity is an exactly
  computed norm of a measured vector and the comparisons are between
  its own local fits; there is no background to detect against.
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
OUT = os.path.join(RES, "audit_alpha_local.txt")

LO, HI = 200_000, 102_400_000


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")
THETA = SPL.THETA


def read_published():
    src = io.open(os.path.join(RES, "audit_denominator.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^  \|sum a\|\s+([+-][\d.]+)\s+([\d.]+)\s+[\d.]+"
                  r"\s*$", src, re.M)
    dec = len(m.group(1).split(".")[1])
    return float(m.group(1)), float(m.group(2)), dec


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
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum())) \
        if n > 2 else float("inf")
    cor = abs(float(np.corrcoef(x, y)[0, 1])) if n > 2 else 0.0
    return float(a), float(np.sqrt((r ** 2).mean())), se, cor


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    puba, pubse, dec = read_published()
    NS = family(LO, HI)
    say("read alpha = %+.6f (s.e. %.6f) from "
        "results/audit_denominator.txt" % (puba, pubse))
    say("  the field, the sieve and theta' are imported from "
        "code/audit_gain_split.py")
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

    S, L2 = [], []
    for N in NS:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        S.append(abs(float(a.sum())))
        L2.append(float(np.sqrt((a ** 2).sum())))
    x = np.log(np.array(NS, dtype=np.float64))
    ys = np.log(np.array(S))
    yl = np.log(np.array(L2))

    # -------------------------------------------------------------- Z1
    ea, ra, sea, _c = fit(x, ys)
    el, rl, sel, _c2 = fit(x, yl)
    say()
    say("Z1  the control")
    rnd = 0.5 * 10.0 ** (-dec)
    d = abs(ea - puba)
    z1 = d <= rnd
    say("  alpha here %+.6f against the published %+.6f, departure "
        "%.8f; the bound from %d decimals is %.8f"
        % (ea, puba, d, dec, rnd))
    say("  and e(l2) here %+.6f, so the whole-field deficit is %+.6f"
        % (el, ea - el))
    say("PRINTBOUND audit_alpha_local %d %.8f" % (dec, rnd))
    say("  Z1 %s   (cap: the printing bound)"
        % ("hold" if z1 else "REFUTED"))

    # ------------------------------------------------------- octaves
    nocts = 0
    while LO * 2 ** (nocts + 1) <= HI:
        nocts += 1
    octs = []
    for j in range(nocts):
        lo, hi = LO * 2 ** j, LO * 2 ** (j + 1)
        sel_ = [i for i, N in enumerate(NS)
                if (lo <= N < hi or (j == nocts - 1 and N == hi))]
        if len(sel_) >= 3:
            octs.append((j, sel_))
    say()
    say("  octave   points   mid log N   alpha_local   s.e.        "
        "e(l2)_local   s.e.        deficit     s.e.")
    locs, defs, cors, mids = [], [], [], []
    for j, idx in octs:
        xo = x[idx]
        a1, r1, s1, c1 = fit(xo, ys[idx])
        a2, r2, s2, c2 = fit(xo, yl[idx])
        dd = a1 - a2
        sd = math.sqrt(s1 * s1 + s2 * s2)
        locs.append((a1, s1))
        defs.append((dd, sd))
        cors.append(c1)
        mids.append(float(xo.mean()))
        say("  %-8d %-8d %-11.4f %+-13.6f %-11.6f %+-13.6f %-11.6f "
            "%+-11.6f %.6f"
            % (j, len(idx), float(xo.mean()), a1, s1, a2, s2, dd, sd))
        say("OCTAVE alphalocal_%d %.4f %+.6f %.6f %+.6f %.6f"
            % (j, float(xo.mean()), a1, s1, dd, sd))
    say("SWEPT alphalocal_octave octave-range %.6f"
        % (max(a for a, _ in locs) - min(a for a, _ in locs)))
    say("POP alphalocal_octave %d" % min(len(i) for _j, i in octs))
    say("CORR alphalocal_octave %.5f" % min(cors))

    # -------------------------------------------------------------- Z2
    say()
    say("Z2  is the top octave genuinely below the global alpha?")
    at, st = locs[-1]
    z2 = (ea - at) > 2.0 * st
    say("  the top octave gives %+.6f (s.e. %.6f) against the global "
        "%+.6f, a departure of %+.6f = %.2f of its own standard errors"
        % (at, st, ea, at - ea, abs(at - ea) / st))
    say("  Z2 %s   (cap 2 standard errors)"
        % ("hold" if z2 else "REFUTED"))

    # -------------------------------------------------------------- Z3
    say()
    say("Z3  is it a drift?")
    mm = np.array(mids)
    aa = np.array([a for a, _ in locs])
    sl, rr, sse, _cc = fit(mm, aa)
    z3 = sl < 0.0 and abs(sl) / sse >= 2.0
    say("  the local exponents run %+.6f to %+.6f across %d octaves"
        % (float(aa.min()), float(aa.max()), len(locs)))
    say("  their least-squares slope against mid log N is %+.6f, "
        "s.e. %.6f, t = %.2f" % (sl, sse, abs(sl) / sse))
    say("TSTAT slope_audit_alpha_local %.2f" % (abs(sl) / sse))
    say("SPREAD slope_audit_alpha_local %.4f"
        % float(mm.max() - mm.min()))
    if abs(sl) / sse < 2.0:
        say("UNRESOLVED SIGN slope_audit_alpha_local")
    say("  Z3 %s   (cap 2 standard errors)"
        % ("hold" if z3 else "REFUTED"))

    # -------------------------------------------------------------- Z4
    say()
    say("Z4  does the deficit close anywhere?")
    z4 = True
    closed = []
    for (j, _idx), (dd, sd) in zip(octs, defs):
        if not (dd > 2.0 * sd):
            z4 = False
            closed.append(j)
    say("  the deficits run %+.6f to %+.6f"
        % (min(d for d, _ in defs), max(d for d, _ in defs)))
    say("  octaves where it is not resolved positive: %s"
        % (", ".join(str(j) for j in closed) if closed else "none"))
    dl, rdl, sdl, _cdl = fit(mm, np.array([d for d, _ in defs]))
    say("  and the deficit's own slope against mid log N is %+.6f, "
        "s.e. %.6f, t = %.2f" % (dl, sdl, abs(dl) / sdl))
    say("TSTAT slope_alphalocal_deficit %.2f" % (abs(dl) / sdl))
    say("SPREAD slope_alphalocal_deficit %.4f"
        % float(mm.max() - mm.min()))
    if abs(dl) / sdl < 2.0:
        say("UNRESOLVED SIGN slope_alphalocal_deficit")
    say("  Z4 %s   (cap 2 standard errors in every octave)"
        % ("hold" if z4 else "REFUTED"))

    say()
    say("=" * 70)
    say("Z1 %s  Z2 %s  Z3 %s  Z4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (z1, z2, z3, z4)))

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) over the squarefree",
        "           k < N^theta' coprime to N, the exponents of",
        "           |sum a| and of the l2 norm in log N, fitted over",
        "           the whole on-field family to 1.024e8 and again",
        "           inside each octave of N on that octave's own",
        "           points; the top octave's exponent against the",
        "           whole-field one; the local exponents' slope",
        "           against mid log N; and the local deficit",
        "           alpha - e(l2) in each octave with its own",
        "           standard error and its slope.",
        "NULL: none is run and none applies. Every quantity is an",
        "      exactly computed norm of a measured vector and the",
        "      comparisons are between its own local fits; there is no",
        "      background to detect against.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta'; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to 102400000; the weighted sum, the",
        "       sieve and theta' are code/audit_gain_split.py's,",
        "       imported; alpha is read from",
        "       results/audit_denominator.txt.",
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
