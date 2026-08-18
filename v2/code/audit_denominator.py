# -*- coding: utf-8 -*-
r"""
What the gain's denominator is, and which restriction makes it that.

WHAT IS AT STAKE

Item 4(b) is stated as an exponent: e(G) = +0.149567 must travel
+0.134019 to reach the measured e(l1/l2) = +0.283586, where
G = l1/|sum a| for a_k = (log k)H(N;k). Everything this repository has
done to it treats |sum a| as a cancelling sum whose smallness is the
thing to be improved. It has never been asked what that sum IS.

It is a truncated Chebyshev-Goldbach correlation. Because mu * log =
Lambda as Dirichlet convolution, an UNRESTRICTED double sum collapses:

    sum_k log k sum_m mu(m) Lambda(N - mk)
        = sum_j Lambda(N - j) sum_{mk = j} mu(m) log k
        = sum_j Lambda(N - j) Lambda(j),

which for even N is the Goldbach-Chebyshev correlation, positive and
of order N times a singular series. The repository's sum is that with
four restrictions on the index: k < N^theta', k squarefree, k coprime
to N, and m coprime to k. Under them the sum comes out NEGATIVE and a
fifth to a quarter of the size, so at least one restriction does not
merely shrink the main term but removes it and overshoots.

Which one, and what the survivor's size is, restates item 4(b) in
terms of an object with arithmetic content instead of a fitted
exponent. If |sum a| is a power N^alpha, the demand e(G) -> theta'/2
is exactly the demand alpha -> e(l1) - theta'/2, and that is a
statement about a correlation sum rather than about a ratio.

BACKS: Remark {#rem:denominator} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  T1  The control, and the framing. With every restriction removed --
      all k >= 1, all m >= 1 with mk < N -- the double sum equals
      sum_j Lambda(j) Lambda(N - j) to 1e-12 relative, at three small
      N. This is the identity the whole reading rests on.
  T2  With the repository's restrictions the sum is negative at every
      N of the field.
  T3  And the main term is gone, not merely reduced: |sum a| is a
      resolved power of N with exponent alpha below 1 at two standard
      errors, while sum_j Lambda(j)Lambda(N-j) divided by N is flat.
  T4  The truncation is the restriction that does it: imposing
      k < N^theta' alone turns the sum negative, while each of the
      other three alone leaves it positive.

REFUTATION RULE (fixed before the run)

  T1  REFUTED above 1e-12 relative at any of the three. Then the
      convolution identity is not what this sum is and every reading
      below is void. THIS ONE GATES.
  T2  REFUTED if the sum is positive anywhere on the field. Then the
      sign is not a property of the restricted sum and the "overshoot"
      reading is wrong.
  T3  REFUTED if alpha is not resolved below 1. Then |sum a| is still
      of main-term order, e(G) is pinned by that order, and the demand
      is a demand about a main term -- which would be a stronger and
      worse result than the one predicted.
  T4  REFUTED if any other single restriction flips the sign, or if
      the truncation alone does not. Then the sign is not the
      truncation's doing and the four restrictions are not separable
      this way.

  T1 gates. T2 to T4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. T1 is an identity, T2 and T4 are
  signs of exactly computed sums, and T3 is a fit to them; there is no
  background to detect against and no threshold that a null would
  calibrate. The sign arms for the gain itself are
  audit_crossk_reference.py and lab_gain_opposition.py.
"""

import importlib.util
import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CODE = os.path.join(ROOT, "code")
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_denominator.txt")

LO, HI = 200_000, 102_400_000
SMALL = (20_000, 50_000, 100_000)


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")
THETA = SPL.THETA


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


def chebyshev(N, lam):
    """sum_j Lambda(j) Lambda(N - j), the untruncated correlation"""
    j = np.arange(1, N, dtype=np.int64)
    return float((lam[j] * lam[N - j]).sum())


def restricted(N, lam, mu, sqf, trunc, sqfree, coprimeN, coprimek):
    """the double sum under any subset of the four restrictions"""
    PN = SPL.factor_set(N)
    K = int(N ** THETA) if trunc else N
    tot = 0.0
    for k in range(2, K):
        if sqfree and not sqf[k]:
            continue
        if coprimeN and any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 1:
            continue
        ms = np.arange(1, M + 1, dtype=np.int64)
        if coprimek:
            for q in SPL.factor_set(k):
                ms = ms[ms % q != 0]
        tot += math.log(k) * float(
            (lam[N - ms * k] * mu[ms].astype(np.float64)).sum())
    return tot


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

    say("the field, the sieve and theta' are imported from "
        "code/audit_gain_split.py")
    say("  theta' = %.2f, so the truncation is k < N^theta'" % THETA)

    # -------------------------------------------------------------- T1
    say()
    say("T1  the identity, with every restriction removed")
    say("  N          unrestricted double sum   sum Lam(j)Lam(N-j)"
        "   relative gap")
    lam, mu = SPL.lambda_and_mu(max(SMALL))
    sqf = mu != 0
    worst = 0.0
    for N in SMALL:
        u = restricted(N, lam, mu, sqf, False, False, False, False)
        c = chebyshev(N, lam)
        rel = abs(u - c) / max(abs(c), 1e-12)
        worst = max(worst, rel)
        say("  %-10d %+-25.4f %+-20.4f %.3e" % (N, u, c, rel))
    t1 = worst <= 1e-12
    say("  worst relative gap %.3e" % worst)
    say("  T1 %s   (tol 1e-12 relative)"
        % ("hold" if t1 else "REFUTED"))
    say("  so the repository's sum a is that correlation with four "
        "restrictions on the index, and nothing else.")

    # -------------------------------------------------------------- T4
    say()
    say("T4  which restriction changes the sign, one at a time")
    say("  N          none          trunc only    squarefree    "
        "coprime to N  m coprime k")
    t4 = True
    for N in SMALL:
        row = [restricted(N, lam, mu, sqf, False, False, False, False),
               restricted(N, lam, mu, sqf, True, False, False, False),
               restricted(N, lam, mu, sqf, False, True, False, False),
               restricted(N, lam, mu, sqf, False, False, True, False),
               restricted(N, lam, mu, sqf, False, False, False, True)]
        if not (row[1] < 0 and all(v > 0 for v in row[2:])):
            t4 = False
        say("  %-10d %+-13.1f %+-13.1f %+-13.1f %+-13.1f %+.1f"
            % (N, row[0], row[1], row[2], row[3], row[4]))
    say("  T4 %s   (cap: the truncation alone negative, each other "
        "alone positive)" % ("hold" if t4 else "REFUTED"))

    # --------------------------------------------------- T2 and T3
    say()
    say("sieving to %d ..." % HI)
    lam, mu = SPL.lambda_and_mu(HI)
    sqf = mu != 0
    NS = family(LO, HI)
    classes = sorted(set(tuple(sorted(SPL.factor_set(N))) for N in NS))
    say("the field: every N = 2^a 5^b with a and b at least one in "
        "[%d, %d]; %d of them" % (LO, HI, len(NS)))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))

    rows, kcounts, flats = [], [], []
    say()
    say("  N            sum a           l1              G        "
        "|sum a|/N   Lam*Lam/N")
    for N in NS:
        ks, a = SPL.weighted(N, lam, mu, sqf)
        S = float(a.sum())
        l1 = float(np.abs(a).sum())
        l2 = float(np.sqrt((a ** 2).sum()))
        flats.append((l1 / l2) / math.sqrt(ks.size))
        c = chebyshev(N, lam)
        rows.append((N, S, l1, l1 / abs(S), abs(S) / N, c / N))
        kcounts.append(float(ks.size))
        say("  %-12d %+-15.1f %-15.1f %-8.4f %-11.6f %.6f"
            % (N, S, l1, l1 / abs(S), abs(S) / N, c / N))

    say("  the count reference sqrt(#k) is used in X1 below, so the "
        "magnitude one is declared here: l1/l2 over sqrt(#k) runs")
    say("  %.4f to %.4f over the field" % (min(flats), max(flats)))
    say("REFERENCE audit_denominator %d %.4f %.4f"
        % (len(rows), min(flats), max(flats)))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))

    say()
    say("T2  is the restricted sum negative on the whole field?")
    neg = sum(1 for r in rows if r[1] < 0)
    t2 = neg == len(rows)
    say("  negative at %d of the %d N" % (neg, len(rows)))
    say("  T2 %s" % ("hold" if t2 else "REFUTED"))

    say()
    say("T3  and how big is what survives?")
    ea, ra, sea = fit(x, np.log(np.array([abs(r[1]) for r in rows])))
    el, rl, sel = fit(x, np.log(np.array([r[2] for r in rows])))
    eg, rg, seg = fit(x, np.log(np.array([r[3] for r in rows])))
    ec, rc, sec = fit(x, np.log(np.array([r[5] for r in rows])))
    say("  quantity        exponent in log N   s.e.        t")
    for nm, e, se in (("|sum a|", ea, sea), ("l1", el, sel),
                      ("G", eg, seg), ("Lam*Lam / N", ec, sec)):
        say("  %-15s %+-19.6f %-11.6f %.2f" % (nm, e, se, abs(e) / se))
    t3 = (1.0 - ea) > 2.0 * sea
    say("  the untruncated correlation over N is flat to %.2f "
        "standard errors, as it must be" % (abs(ec) / sec))
    say("  |sum a| falls below main-term order by %+.6f = %.2f "
        "standard errors" % (ea - 1.0, (1.0 - ea) / sea))
    say("TSTAT slope_audit_denominator %.2f" % (abs(ea) / sea))
    say("SPREAD slope_audit_denominator %.4f"
        % float(x.max() - x.min()))
    if abs(ea) / sea < 2.0:
        say("UNRESOLVED SIGN slope_audit_denominator")
    say("  T3 %s   (cap 2 standard errors below 1)"
        % ("hold" if t3 else "REFUTED"))

    say()
    say("  item 4(b) restated in these terms. G = l1/|sum a| gives")
    say("  e(G) = e(l1) - alpha identically, and the three fitted "
        "here are")
    say("  %+.6f = %+.6f - %+.6f, closing to %.1e"
        % (eg, el, ea, abs(eg - (el - ea))))
    need = el - THETA / 2.0
    say("  so the demand e(G) -> theta'/2 = %.4f is the demand"
        % (THETA / 2.0))
    say("  alpha -> e(l1) - theta'/2 = %+.6f, against the measured "
        "%+.6f" % (need, ea))
    say("  that is: |sum a| must be smaller than it is by a factor "
        "growing like")
    say("  N^%+.6f, on a sum that is a truncated Lambda*Lambda "
        "correlation." % (ea - need))

    say()
    say("  X1  what the shortfall from main-term order equals")
    say("  (written after T3; not pre-registered). The correlation is")
    say("  of order N and the restricted sum is of order N^alpha, so "
        "the")
    say("  truncation costs N^(1-alpha). Measured against the "
        "arithmetic")
    say("  ceiling theta'/2 that #k ~ N^theta' puts on a square root:")
    ek, rk, sek = fit(x, np.log(np.array(kcounts)))
    gap = (1.0 - ea) - THETA / 2.0
    say("  1 - alpha = %+.6f against theta'/2 = %.4f, a gap of "
        "%+.6f" % (1.0 - ea, THETA / 2.0, gap))
    say("  which is %.2f standard errors of alpha" % (abs(gap) / sea))
    say("  and against half the measured #k exponent, %+.6f, the gap "
        "is %+.6f" % (ek / 2.0, (1.0 - ea) - ek / 2.0))
    say("  so the truncation leaves the correlation divided by the "
        "square root of")
    say("  the index set it kept: |sum a| ~ (Lambda*Lambda) / "
        "sqrt(#k).")
    say()
    say("  That turns the demand into one line. e(G) = e(l1) - alpha "
        "and")
    say("  alpha = 1 - e(#k)/2 to the precision just measured, so")
    say("  e(G) -> e(#k)/2 is e(l1) -> 1: the demand of item 4(b) is "
        "that")
    say("  the l1 norm of the weighted dilation sums reach main-term "
        "order.")
    say("  It is measured at %+.6f, short by %+.6f." % (el, 1.0 - el))

    say()
    say("=" * 70)
    say("T1 %s  T2 %s  T3 %s  T4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (t1, t2, t3, t4)))

    head = [
        "STATISTIC: the double sum sum_k log k sum_m mu(m)",
        "           Lambda(N - mk) under each subset of the four",
        "           restrictions (k < N^theta', k squarefree, k",
        "           coprime to N, m coprime to k), against the",
        "           untruncated Chebyshev correlation",
        "           sum_j Lambda(j) Lambda(N - j); and, on the",
        "           on-field family to 1.024e8, the restricted sum a,",
        "           its l1 norm, the gain l1/|sum a| and the",
        "           correlation over N, with each one's least-squares",
        "           exponent in log N.",
        "NULL: none is run and none applies. T1 is an identity, T2 and",
        "      T4 are signs of exactly computed sums, T3 is a fit to",
        "      them; there is no background to detect against and no",
        "      threshold a null would calibrate. The sign arms for the",
        "      gain are audit_crossk_reference.py and",
        "      lab_gain_opposition.py.",
        "FIELD: for T1 and T4, N = 20000, 50000, 100000, small enough",
        "       that the unrestricted double sum over every k < N is",
        "       affordable; for T2 and T3, N = 2^a 5^b with BOTH",
        "       a >= 1 and b >= 1 in [2e5, 1.024e8], one coprimality",
        "       class as COPRIME says. Lambda and mu from an integer",
        "       sieve; the restricted sum, the field and theta' are",
        "       code/audit_gain_split.py's, imported.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not t1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
