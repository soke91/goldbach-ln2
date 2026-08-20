# -*- coding: utf-8 -*-
r"""
The radical control, registered in advance and run on radicals never used

WHAT IS AT STAKE

rem:basecontrol found that every measurement in this branch had used
one base per radical, ran the missing control, and then could not
claim it: P2's unresolved clause fired on an error taken from the
other radical group, so the reading stayed barred, and the statistic
that answered the question -- each group's drifts against one constant
on their own errors -- was chosen after the numbers were on the page.
That remark named what would settle it: **register the chi-square in
advance and run it on a radical none of these used**, so the test is
both correct and blind.

Used already: {2}, {2,3}, {2,5}, {2,7}, {2,11}, {2,13}, {2,17},
{2,23}, {2,47}, {2,101}, {2,3,5}, {2,3,5,7,11,13}.  Fresh here:

    {2,19}    2 19^3 = 13718     2^10 19 = 19456    2^6 19^2 = 23104
              2^2 19^3 = 27436   2^11 19 = 38912
    {2,3,7}   2^2 3^4 7^2 = 15876   2^8 3^2 7 = 16128
              2^4 3^3 7^2 = 21168   2^10 3 7 = 21504
              2 3^5 7^2 = 23814

Five bases each, 2-adic valuation spanning 1 to 11 and 1 to 10.  Base
20736 from rem:basecontrol rides along as the gate on the code path.

**And the power is registered too.**  A chi-square that comes in under
its critical value has failed to reject a constant; it has not proved
one.  What size of base-to-base scatter this design could have caught
is arithmetic on the errors, and it is printed rather than left
implicit: an extra common scatter s inflates the expected chi-square
by s^2 * sum(1/sigma_i^2), so the smallest detectable s is
sqrt((chi2crit - dof) / sum(1/sigma_i^2)).  Anything below that is not
excluded by a passing test and this run says so in its own output.

BACKS: Remark {#rem:radicalblind} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  THE GATE.  Base 20736 reproduces rem:basecontrol's ten-N drift
      +0.213528 to six decimals -- same base, same N, same code path.
  Q2  **The {2,19} bases are one constant.**  Their chi-square against
      a weighted mean, on their own errors, is at or below 9.49, the
      five per cent point for four degrees of freedom.
  Q3  The {2,3,7} bases likewise, same statistic and same cap.
  Q4  And the two fresh radicals are apart: their weighted means
      differ at |t| above 3, on the errors of those means.  {2,3,7}
      contains the 3 that rem:primecontrib found carrying almost all
      of the prime dependence, and {2,19} does not.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED outside six decimals; nothing below is reported.
  Q2  **REFUTED above 9.49, and this is the most expensive outcome
      this branch can produce.**  Five bases sharing one radical and
      disagreeing beyond their own errors would mean the drift is not
      a function of the radical, and rem:whichfloor, rem:radicallaw
      and rem:primecontrib would each have measured a quantity that
      moves with the base -- numbers standing, readings withdrawn.
      If it fires, the drifts are to be reported against 2-adic
      valuation and base size so the failure names what it is instead.
  Q3  REFUTED above 9.49, same consequence for the three-prime case.
  Q4  REFUTED at |t| of 3 or below.  Then two radicals that differ in
      whether 3 divides N are not separated, which would contradict
      rem:radicallaw's largest measured effect and put the
      disagreement with that run rather than with the control.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  A chi-square below
  its critical value is a failure to reject, not a proof.  **This run
  prints the smallest base-to-base scatter it could have detected**,
  and a passing Q2 or Q3 licenses only "no base dependence larger
  than that number", never "the drift is a function of the radical"
  full stop.  The remark must quote the number beside the verdict.
  Five bases give four degrees of freedom and that is a weak test; it
  is the test this branch registered and its weakness is part of the
  result.  This is the seventh time a cap here has been written on a
  quantity that was not the one at issue, and the correction this
  time is to state the power rather than to widen the cap.

  WHAT THIS CANNOT DO.  Two radicals, five bases each.  The 2-adic
  valuation and the base size are varied together and not separated
  from one another.  Nothing here measures |sum a| against any floor
  or moves item 5's demand; rem:shapepower and rem:deficitlog stand.
"""

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
OUT = os.path.join(ROOT, "results", "audit_radical_blind.txt")
SRC = os.path.join(ROOT, "results", "audit_base_control.txt")

THETA = 0.56
GROUPS = [("{2,19}", [13_718, 19_456, 23_104, 27_436, 38_912]),
          ("{2,3,7}", [15_876, 16_128, 21_168, 21_504, 23_814])]
GATEBASE = 20_736
NPER = 10
DEC = 6
CHICRIT = 9.49
TCAP = 3.0


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int64)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def factor_set(n):
    v, out, d = n, set(), 2
    while d * d <= v:
        if v % d == 0:
            out.add(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        out.add(v)
    return out


def v2(n):
    e = 0
    while n % 2 == 0:
        n //= 2
        e += 1
    return e


def pair(N, lam, mu, sqf):
    PN = factor_set(N)
    K = int(N ** THETA)
    lk = np.zeros(N, dtype=np.float64)
    l2sq = 0.0
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        lg = math.log(k)
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        mm = mu[ms].astype(np.float64)
        l2sq += (lg * float((lam[N - ms * k] * mm).sum())) ** 2
        lk[ms * k] += lg * mm
        del ms, mm
    j = np.arange(1, N, dtype=np.int64)
    sa = abs(float((lam[N - j] * lk[1:]).sum()))
    del j, lk
    return sa, math.sqrt(l2sq)


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se


def drift(base, lam, mu, sqf):
    xs, ys = [], []
    for jj in range(NPER):
        N = base * (1 << jj)
        sa, l2 = pair(N, lam, mu, sqf)
        xs.append(math.log(N))
        ys.append(math.log(sa / l2))
    return fit(xs, ys)


def read_pub():
    m = re.search(r"^POINT basedrift_%d ([-+]?[\d.]+)\s*$" % GATEBASE,
                  io.open(SRC, encoding="utf-8").read(), re.M)
    if not m:
        raise SystemExit("no basedrift marker for %d" % GATEBASE)
    return float(m.group(1))


HEAD = [
    "STATISTIC: the drift of log(|sum a|/l2) against log N for five",
    "           bases on each of two radicals never used in this",
    "           branch, their chi-square against one constant per",
    "           radical, and the smallest base scatter that test",
    "           could have detected.",
    "FIELD: base * 2^j for j < %d; bases %s on %s and %s on %s,"
    % (NPER, GROUPS[0][1], GROUPS[0][0], GROUPS[1][1], GROUPS[1][0]),
    "       spanning 2-adic valuation 1 to 11. Base %d rides along as"
    % GATEBASE,
    "       the gate. k over the squarefree k < N^%.2f coprime to N;"
    % THETA,
    "       j over every index below N. The gate drift is READ from",
    "       results/audit_base_control.txt.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    say("READ audit_base_control.txt %d %.6f" % (GATEBASE, pub))
    say("  the ten-N drift rem:basecontrol published for the gate base")
    say("PRINTBOUND audit_radical_blind %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, chi-square cap %.2f on 4 degrees, |t| cap %.1f"
        % (THETA, CHICRIT, TCAP))
    say("RADICALS %d" % len(GROUPS))

    NMAX = max([GATEBASE] + [b for _, g in GROUPS for b in g]) \
        * (1 << (NPER - 1))
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    # -------------------------------------------------------------- Q1
    gb, gse = drift(GATEBASE, lam, mu, sqf)
    say()
    say("Q1  the gate")
    q1 = abs(gb - pub) < 10.0 ** (-DEC)
    say("  base %-6d here %+.6f against its %+.6f  %s"
        % (GATEBASE, gb, pub, "ok" if q1 else "MISMATCH"))
    say("  Q1 %s   (cap: %d decimals)"
        % ("hold" if q1 else "REFUTED", DEC))
    if not q1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    res = {}
    for tag, bases in GROUPS:
        say()
        say("radical %s" % tag)
        rows = []
        for base in bases:
            b, se = drift(base, lam, mu, sqf)
            rows.append((base, v2(base), b, se))
            say("  base %-6d v2 %-3d drift %+.6f +- %.6f"
                % (base, v2(base), b, se))
            say("POINT blinddrift_%d %.6f" % (base, b))
            say("SPREAD blinddrift_%d %.6f" % (base, se))
        res[tag] = rows
    say("SCALES %d" % (sum(len(g) for _, g in GROUPS) * NPER))

    # ---------------------------------------------------------- Q2, Q3
    stats = {}
    for name, (tag, _) in zip(("Q2", "Q3"), GROUPS):
        rows = res[tag]
        v = np.array([r[2] for r in rows])
        e = np.array([r[3] for r in rows])
        w = 1.0 / e ** 2
        wm = float((v * w).sum() / w.sum())
        wse = float(1.0 / math.sqrt(w.sum()))
        chi = float((((v - wm) / e) ** 2).sum())
        dof = len(rows) - 1
        detect = math.sqrt(max(CHICRIT - dof, 0.0) / float(w.sum()))
        ok = chi <= CHICRIT
        stats[tag] = (wm, wse, chi, detect, ok)
        say()
        say("%s  are the %s bases one constant?" % (name, tag))
        say("  weighted mean %+.6f +- %.6f" % (wm, wse))
        say("  chi-square %.3f on %d degrees, cap %.2f"
            % (chi, dof, CHICRIT))
        say("POINT blindchi_%s %.6f"
            % (tag.replace("{", "").replace("}", "").replace(",", "_"),
               chi))
        say("  smallest base scatter this test could detect: %.6f"
            % detect)
        say("POINT blinddetect_%s %.6f"
            % (tag.replace("{", "").replace("}", "").replace(",", "_"),
               detect))
        say("  %s %s   (cap: %.2f)"
            % (name, "hold" if ok else "REFUTED", CHICRIT))
        if ok:
            say("  a failure to reject, not a proof: base dependence "
                "below %.6f" % detect)
            say("  is not excluded, as the rule says")
        else:
            vv = np.array([float(r[1]) for r in rows])
            ss = np.array([math.log(float(r[0])) for r in rows])
            for nm, x in (("v2", vv), ("log base", ss)):
                b, se = fit(x, v)
                say("    drift on %-9s slope %+.6f +- %.6f, t %+.2f"
                    % (nm, b, se, b / se))
                say("TSTAT blindon%s_%s %.2f"
                    % (nm.replace(" ", ""),
                       tag.replace("{", "").replace("}", "")
                       .replace(",", "_"), b / se))
                say("SPREAD blindon%s_%s %.6f"
                    % (nm.replace(" ", ""),
                       tag.replace("{", "").replace("}", "")
                       .replace(",", "_"), se))

    # -------------------------------------------------------------- Q4
    say()
    say("Q4  are the two fresh radicals apart?")
    (m1, s1, _, _, _) = stats[GROUPS[0][0]]
    (m2, s2, _, _, _) = stats[GROUPS[1][0]]
    d = m2 - m1
    sd = math.sqrt(s1 ** 2 + s2 ** 2)
    q4 = abs(d) / sd > TCAP
    say("  %s %+.6f, %s %+.6f" % (GROUPS[0][0], m1, GROUPS[1][0], m2))
    say("  difference %+.6f +- %.6f, t %+.2f" % (d, sd, d / sd))
    say("TSTAT blindgap %.2f" % (d / sd))
    say("SPREAD blindgap %.6f" % sd)
    say("  Q4 %s   (cap: |t| above %.1f)"
        % ("hold" if q4 else "REFUTED", TCAP))

    # a diagnostic, after the verdicts and predicted by nothing:
    # in bases of a fixed size the 2-adic valuation cannot rise
    # without the odd primes' valuations falling, so a v2 slope is
    # not separable from a v_p slope in this design.
    say()
    say("  a diagnostic, after the verdicts and predicted by nothing")
    say("  a base of fixed size cannot raise v2 without lowering the "
        "odd primes'")
    say("  valuations, so any v2 slope here is not separable from a "
        "v_p one.")
    for tag, bases in GROUPS:
        say("    %s" % tag)
        vv, pp = [], []
        for base in bases:
            n = base
            e2 = v2(n)
            odd = []
            for q in sorted(factor_set(base)):
                if q == 2:
                    continue
                e = 0
                m = base
                while m % q == 0:
                    m //= q
                    e += 1
                odd.append((q, e))
            vv.append(float(e2))
            pp.append(float(sum(e for _, e in odd)))
            say("      base %-6d v2 %-3d odd %s"
                % (base, e2, ", ".join("%d^%d" % t for t in odd)))
        c = float(np.corrcoef(np.array(vv), np.array(pp))[0, 1])
        say("      corr(v2, total odd valuation) %+.6f" % c)
        say("CORR blindvalu_%s_regressors %.6f"
            % (tag.replace("{", "").replace("}", "")
               .replace(",", "_"), abs(c)))
        if abs(c) >= 0.99:
            say("COEFF NOT SEPARABLE blindvalu_%s"
                % tag.replace("{", "").replace("}", "")
                .replace(",", "_"))
    say("  so the v2 slope printed for %s names a direction and not "
        "a cause;" % GROUPS[0][0])
    say("  what Q2 establishes is that the bases disagree, not what "
        "they disagree")
    say("  along")

    say()
    say("=" * 70)
    q2 = stats[GROUPS[0][0]][4]
    q3 = stats[GROUPS[1][0]][4]
    say("Q1 %s  Q2 %s  Q3 %s  Q4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (q1, q2, q3, q4)))
    say()
    if q2 and q3 and q4:
        say("the control passes, blind and registered in advance, at "
            "two radicals")
        say("this branch had never touched. what it licenses is what "
            "its power")
        say("allows: no base dependence larger than the detection "
            "sizes printed")
        say("above, and two radicals separated at the t printed for "
            "Q4. that is")
        say("the first controlled statement this branch has about "
            "the radical, and")
        say("it is a bound on base dependence rather than a proof of "
            "its absence.")
    elif not (q2 and q3):
        say("bases sharing one radical disagree beyond their own "
            "errors. the drift")
        say("is not a function of the radical, and rem:whichfloor, "
            "rem:radicallaw")
        say("and rem:primecontrib each measured a quantity that moves "
            "with the")
        say("base -- their numbers stand and their readings are "
            "withdrawn. what it")
        say("moves with is reported above.")
    else:
        say("each radical is internally consistent and the two are "
            "not separated.")
        say("that contradicts rem:radicallaw's largest measured "
            "effect, so the")
        say("disagreement is with that run and not with the control, "
            "and this")
        say("branch has two measurements of one quantity that do not "
            "agree.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
