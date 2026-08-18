# -*- coding: utf-8 -*-
r"""
How the three exponents move with theta', and where the demand meets.

WHAT IS AT STAKE

{#rem:denominator} left item 4(b) as one line: e(G) = e(l1) - alpha
identically, alpha = 1 - e(#k)/2 to 0.88 standard errors, so the
demand e(G) -> theta'/2 is the demand e(l1) -> 1. Put beside
{#rem:headsign}, which measured the inner sums cancelling at exactly
the independent-sign rate (|I| times sqrt(n_k) is 1.0310 to 1.0995
with an unresolved drift), that suggests a law rather than a fit:

    |H(N;k)| ~ sqrt(N/k) up to logs, so
    l1 = sum_k (log k)|H| ~ N^{(1+theta')/2},
    alpha = 1 - theta'/2,   e(G) = e(l1) - alpha = theta' - 1/2.

The demand is e(G) = theta'/2. Those two lines in theta' meet only at
theta' = 1 -- the whole range, where the reduction is vacuous. If that
is right, item 4(b) is not merely hard at theta' = 0.56; it is
unreachable at every theta' the Huang-Li reduction can use, and the
reason is arithmetic rather than a shortfall to be closed.

It is testable without any new theory: sweep theta', refit the three
exponents at each, and see whether they move on the predicted lines.
The log factors are large over the computable range, so the LEVELS
will not match the model; the SLOPES in theta' are the test, since
logs do not depend on theta' at leading order.

BACKS: Remark {#rem:thetalaw} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  The control. At theta' = 0.56 the three exponents reproduce
      results/audit_denominator.txt inside the bound its printing
      forces.
  U2  e(l1) rises with theta', resolved at two standard errors.
  U3  alpha falls with theta', resolved at two standard errors.
  U4  And the demand recedes: the gap e(G) - theta'/2 is increasing
      in theta', resolved at two standard errors, so raising the
      level makes the demand harder to meet and not easier.
  U5  The crossing where e(G) meets theta'/2 lies above theta' = 0.9,
      with a bracket from the fits' own errors.
  U6  The slopes are the model's: d e(l1)/d theta' = 1/2,
      d alpha/d theta' = -1/2 and d e(G)/d theta' = 1, each within
      two standard errors of the fitted slope.

REFUTATION RULE (fixed before the run)

  U1  REFUTED outside the printing bound. Then this is not the sweep
      {#rem:denominator} measured. THIS ONE GATES.
  U2  REFUTED if e(l1) does not rise. Then l1 is not built from more
      k as the level grows and the whole picture is wrong.
  U3  REFUTED if alpha does not fall. Then the truncation's
      square-root leaving is not what {#rem:denominator} measured it
      to be, off theta' = 0.56.
  U4  REFUTED if the gap is flat or shrinking. That is the outcome
      worth having: the demand would then be met by raising the
      level, and the level is the one axis the reduction lets us
      choose.
  U5  REFUTED if the crossing is at or below 0.9. Then some usable
      level meets the demand and item 4(b) has an address rather than
      only a size.
  U6  REFUTED on any of the three. Then the square-root law is not
      what governs the theta'-dependence, the levels' disagreement is
      not only logs, and U5's extrapolation rests on a line that is
      not the model's.

  U1 gates. U2 to U6 are the measurement and do not gate.

  NO NULL IS RUN and none applies. Every quantity is an exactly
  computed sum and the comparisons are between its exponents at
  different theta'; there is no background to detect against. The
  sign arms for the gain are audit_crossk_reference.py and
  lab_gain_opposition.py, and no sign is randomised here.
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
OUT = os.path.join(RES, "audit_theta_law.txt")

LO, HI = 200_000, 102_400_000
THETAS = (0.40, 0.46, 0.52, 0.56, 0.60, 0.64)
BASE = 0.56


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")


def read_published():
    """{#rem:denominator}'s three exponents at theta' = 0.56"""
    src = io.open(os.path.join(RES, "audit_denominator.txt"),
                  encoding="utf-8").read()
    out, dec = {}, 0
    for nm in ("|sum a|", "l1", "G"):
        m = re.search(r"^  " + re.escape(nm) +
                      r"\s+([+-][\d.]+)\s+([\d.]+)\s+[\d.]+\s*$",
                      src, re.M)
        out[nm] = float(m.group(1))
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
    return float(a), float(b), float(np.sqrt((r ** 2).mean())), se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub, dec = read_published()
    NS = family(LO, HI)
    say("read the three exponents at theta' = %.2f from "
        "results/audit_denominator.txt: %s"
        % (BASE, ", ".join("%s %+.6f" % (k, v)
                           for k, v in pub.items())))
    say("  the field, the sieve and the weighted sum are imported "
        "from code/audit_gain_split.py; only its THETA is varied")
    say()
    say("the field: every N = 2^a 5^b with a and b at least one in "
        "[%d, %d]; %d of them" % (LO, HI, len(NS)))
    classes = sorted(set(tuple(sorted(SPL.factor_set(N))) for N in NS))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))
    say("  the levels swept are theta' = %s"
        % ", ".join("%.2f" % t for t in THETAS))

    say()
    say("sieving to %d ..." % HI)
    lam, mu = SPL.lambda_and_mu(HI)
    sqf = mu != 0

    x = np.log(np.array(NS, dtype=np.float64))
    table = {}
    say()
    say("  theta'   e(#k)        e(l1)        alpha        e(G)      "
        "   s.e.(e(G))  #k at the top N")
    for th in THETAS:
        SPL.THETA = th
        l1s, sas, nks = [], [], []
        for N in NS:
            ks, a = SPL.weighted(N, lam, mu, sqf)
            l1s.append(float(np.abs(a).sum()))
            sas.append(abs(float(a.sum())))
            nks.append(float(ks.size))
        ek = fit(x, np.log(np.array(nks)))
        el = fit(x, np.log(np.array(l1s)))
        ea = fit(x, np.log(np.array(sas)))
        eg = fit(x, np.log(np.array(l1s) / np.array(sas)))
        table[th] = (ek, el, ea, eg, nks[-1])
        say("  %-8.2f %+-12.6f %+-12.6f %+-12.6f %+-12.6f %-11.6f %d"
            % (th, ek[0], el[0], ea[0], eg[0], eg[3], int(nks[-1])))
    SPL.THETA = BASE

    # -------------------------------------------------------------- U1
    say()
    say("U1  the control at theta' = %.2f" % BASE)
    rnd = 0.5 * 10.0 ** (-dec)
    got = {"|sum a|": table[BASE][2][0], "l1": table[BASE][1][0],
           "G": table[BASE][3][0]}
    worst = max(abs(got[k] - pub[k]) for k in pub)
    u1 = worst <= rnd
    for k in ("|sum a|", "l1", "G"):
        say("  %-9s here %+.6f against the published %+.6f"
            % (k, got[k], pub[k]))
    say("  worst departure %.8f; the table prints %d decimals, so "
        "the bound is %.8f" % (worst, dec, rnd))
    say("PRINTBOUND audit_theta_law %d %.8f" % (dec, rnd))
    say("  U1 %s   (cap: the printing bound)"
        % ("hold" if u1 else "REFUTED"))

    # ------------------------------------------------- U2, U3, U6
    t = np.array(THETAS)
    say()
    say("U2/U3/U6  how each exponent moves with theta'")
    say("  quantity   slope in theta'   s.e.        model   "
        "departure in s.e.")
    slopes = {}
    for nm, j, model in (("e(l1)", 1, 0.5), ("alpha", 2, -0.5),
                         ("e(G)", 3, 1.0)):
        y = np.array([table[th][j][0] for th in THETAS])
        s, b, r, se = fit(t, y)
        slopes[nm] = (s, se, b)
        say("  %-10s %+-17.6f %-11.6f %+-7.2f %.2f"
            % (nm, s, se, model, abs(s - model) / se))
        say("TSTAT slope_thetalaw_%s %.2f"
            % (nm.replace("(", "").replace(")", "").replace("1", "one"),
               abs(s) / se))
        say("SPREAD slope_thetalaw_%s %.4f"
            % (nm.replace("(", "").replace(")", "").replace("1", "one"),
               float(t.max() - t.min())))
        if abs(s) / se < 2.0:
            say("UNRESOLVED SIGN slope_thetalaw_%s"
                % nm.replace("(", "").replace(")", "").replace("1",
                                                               "one"))
    u2 = slopes["e(l1)"][0] > 0 and \
        abs(slopes["e(l1)"][0]) / slopes["e(l1)"][1] >= 2.0
    u3 = slopes["alpha"][0] < 0 and \
        abs(slopes["alpha"][0]) / slopes["alpha"][1] >= 2.0
    u6 = all(abs(slopes[nm][0] - mo) <= 2.0 * slopes[nm][1]
             for nm, mo in (("e(l1)", 0.5), ("alpha", -0.5),
                            ("e(G)", 1.0)))
    say("  U2 %s   U3 %s   U6 %s   (caps 2 standard errors)"
        % ("hold" if u2 else "REFUTED", "hold" if u3 else "REFUTED",
           "hold" if u6 else "REFUTED"))

    # -------------------------------------------------------------- U4
    say()
    say("U4  does the demand recede as the level rises?")
    gap = np.array([table[th][3][0] - th / 2.0 for th in THETAS])
    sg, bg, rg, seg = fit(t, gap)
    u4 = sg > 0 and abs(sg) / seg >= 2.0
    say("  theta'   e(G)         theta'/2     gap")
    for i, th in enumerate(THETAS):
        say("  %-8.2f %+-12.6f %-12.6f %+.6f"
            % (th, table[th][3][0], th / 2.0, gap[i]))
    say("  the gap's slope in theta' is %+.6f, s.e. %.6f, t = %.2f"
        % (sg, seg, abs(sg) / seg))
    say("TSTAT slope_audit_theta_law %.2f" % (abs(sg) / seg))
    say("SPREAD slope_audit_theta_law %.4f" % float(t.max() - t.min()))
    if abs(sg) / seg < 2.0:
        say("UNRESOLVED SIGN slope_audit_theta_law")
    say("  U4 %s   (cap 2 standard errors)"
        % ("hold" if u4 else "REFUTED"))

    # -------------------------------------------------------------- U5
    say()
    say("U5  where do the two lines meet?")
    sG, seG, bG = slopes["e(G)"]
    cross = (0.0 - bG) / (sG - 0.5) if abs(sG - 0.5) > 1e-12 \
        else float("inf")
    lo = (0.0 - (bG + 2.0 * seG * t.mean())) / \
        ((sG + 2.0 * seG) - 0.5)
    hi = (0.0 - (bG - 2.0 * seG * t.mean())) / \
        ((sG - 2.0 * seG) - 0.5)
    a_, b_ = min(lo, hi), max(lo, hi)
    u5 = cross > 0.9
    say("  e(G) fitted as %+.6f theta' %+.6f, against the demand "
        "0.5 theta'" % (sG, bG))
    say("  they meet at theta' = %.4f" % cross)
    say("BRACKET theta_where_gain_meets_ceiling %.4f %.4f %.4f"
        % (cross, a_, b_))
    lohalf = fit(t[:3], np.array([table[th][3][0]
                                  for th in THETAS[:3]]))
    hihalf = fit(t[3:], np.array([table[th][3][0]
                                  for th in THETAS[3:]]))
    say("DRIFT theta_where_gain_meets_ceiling %.6f"
        % abs(lohalf[0] - hihalf[0]))
    say("  the slope refitted on the lower three levels is %+.6f and "
        "on the upper three %+.6f, so the drift declared above is "
        "their difference" % (lohalf[0], hihalf[0]))
    say("  U5 %s   (cap: a crossing at or below 0.9)"
        % ("hold" if u5 else "REFUTED"))

    say()
    say("=" * 70)
    say("U1 %s  U2 %s  U3 %s  U4 %s  U5 %s  U6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (u1, u2, u3, u4, u5, u6)))

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) over the squarefree",
        "           k < N^theta' coprime to N, the counts #k, the l1",
        "           norm, |sum a| and the gain l1/|sum a|, each",
        "           fitted as a power of N over the on-field family",
        "           to 1.024e8, at every theta' in 0.40, 0.46, 0.52,",
        "           0.56, 0.60, 0.64; then each exponent's slope in",
        "           theta' against the square-root law's 1/2, -1/2",
        "           and 1; the gap e(G) - theta'/2 and its slope; and",
        "           the theta' at which the fitted e(G) meets the",
        "           demand theta'/2, with a bracket from the fit's own",
        "           errors and a drift from refitting on each half of",
        "           the swept levels.",
        "NULL: none is run and none applies. Every quantity is an",
        "      exactly computed sum and the comparisons are between",
        "      its exponents at different theta'; there is no",
        "      background to detect against and no sign is randomised.",
        "      The sign arms for the gain are",
        "      audit_crossk_reference.py and lab_gain_opposition.py.",
        "FIELD: N = 2^a 5^b with BOTH a >= 1 and b >= 1 in",
        "       [2e5, 1.024e8], one coprimality class as COPRIME",
        "       says -- the class {2,5}, k coprime to 10 and N even;",
        "       k squarefree with 2 <= k < N^theta' for each swept",
        "       theta'; m over 1 <= m < N/k with (m,k) = 1; Lambda",
        "       and mu from an integer sieve to 102400000; the",
        "       weighted sum, the sieve and the field are",
        "       code/audit_gain_split.py's, imported, with its THETA",
        "       set to each swept level in turn and restored to 0.56",
        "       afterwards; the published exponents are read from",
        "       results/audit_denominator.txt.",
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
