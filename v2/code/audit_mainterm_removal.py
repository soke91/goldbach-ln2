# -*- coding: utf-8 -*-
r"""
The demand restated: how exactly a main term has to be removed

WHAT IS AT STAKE

rem:bilinearcancel found that sum a's smallness is cancellation over
the modulus d, that the cancellation runs 2.96 to 4.30 times better
than square-root and improving, and that the single largest piece is
d = 1 at 104 times l2.  That last number locates the whole problem in
one place, so split the sum there.  Writing

    A = I(1),    B = sum_{d >= 2} mu(d) I(d),    S = -(A + B),

A is a Goldbach-type count over m >= K -- a main term, of order N --
and B has to cancel it.  The demand |S| <~ l2 then says

    |A + B| / |A|  <~  N^(-1/2)

while rem:denominator's alpha puts the achieved ratio at
N^(alpha - 1) = N^(-0.282084).  **So item 5 is the statement that a
main term of size N is removed to a relative accuracy of N^(-1/2),
and what is achieved is N^(-0.282).**  That is the sharpest form the
requirement has taken in this repository, and it is measurable.

Eight N over 2.1 decades, all of it computed from the Type II pieces
rem:bilinear verified.

BACKS: Remark {#rem:maintermremoval} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  THE GATE.  -(A + B) at N = 200000 reproduces rem:bilinear's
      covered part, read from its full-precision POINT marker, to a
      relative 1e-12.
  Z2  **B fights A.**  A and B have opposite signs at every N.
  Z3  **The residue's exponent is alpha - 1.**  Fitting
      log(|A + B| / |A|) on log N gives a slope within 0.05 of
      -0.282084.
  Z4  And A is main-term order: the exponent of |A| is within 0.05
      of 1.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED outside 1e-12; nothing below is reported.  The marker
      is read at full double precision because ten digits was not
      enough last time -- see rem:bilinearcancel's two TOL BELOW
      PRINT failures.
  Z2  REFUTED by any N where the signs agree.  Then B is not
      cancelling A and the split is not the structure this run
      assumes; the reading below would not be available.
  Z3  **REFUTED outside 0.05.**  Then the residue does not shrink at
      the rate alpha says, and either alpha or this split is
      measuring something else -- the two would have to be
      reconciled before either is used.
  Z4  REFUTED outside 0.05 of exponent 1.  Then A is not a main term
      and calling the requirement "a main term removed to relative
      N^(-1/2)" is the wrong description, however the rest comes out.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  Eight points give
  six degrees of freedom on a slope and the standard errors are
  printed.  **A slope inside 0.05 of the target whose own error
  exceeds 0.05 has not confirmed anything** -- it has failed to
  distinguish the target from a range that includes it.  Z3 and Z4
  are to be read against their printed errors, and if an error
  exceeds the cap the verdict word stands without a reading, as in
  rem:levelfine and rem:radicalblind.

  WHAT THIS CANNOT DO.  The exponent -1/2 for the demand comes from
  l2 being of order sqrt(N) times logs, which is an asymptotic
  reading and not something measured here; the fitted e(l2) on this
  repository's field is +0.583897, and over eight N a log power is
  not separable from a small shift in the exponent.  **The two
  statements of the demand -- relative N^(-1/2) and relative
  N^(e(l2)-1) -- differ by 0.084 and this run does not choose between
  them.**  Nothing here bounds anything.
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
OUT = os.path.join(ROOT, "results", "audit_mainterm_removal.txt")
SRCB = os.path.join(ROOT, "results", "audit_bilinear.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000]
NGATE = 200_000
RELID = 1e-12
ALPHAM1 = -0.282084
SLOPECAP = 0.05
MAINEXP = 1.0


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


def split(N, lam, mu, sqf):
    """A = I(1), B = sum over d >= 2, and l2 over the k-range"""
    PN = factor_set(N)
    K = int(N ** THETA)
    A = 0.0
    B = 0.0
    for d in range(1, (N - 1) // K + 1):
        md = int(mu[d])
        if md == 0 or any(d % q == 0 for q in PN):
            continue
        ms = np.arange(K, (N - 1) // d + 1, dtype=np.int64)
        if ms.size == 0:
            continue
        keep = sqf[ms]
        for q in factor_set(d) | PN:
            keep &= (ms % int(q)) != 0
        if d == 1:
            keep &= lam[ms] == 0.0
        ms = ms[keep]
        if ms.size == 0:
            continue
        c = float((lam[N - d * ms]
                   * np.log(ms.astype(np.float64))).sum())
        if d == 1:
            A = md * c
        else:
            B += md * c
        del ms, keep
    l2sq = 0.0
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        l2sq += (math.log(k)
                 * float((lam[N - ms * k]
                          * mu[ms].astype(np.float64)).sum())) ** 2
        del ms
    return A, B, math.sqrt(l2sq), K


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se


def read_pub():
    m = re.search(r"^POINT covered_%d ([-+]?[\d.eE+-]+)\s*$" % NGATE,
                  io.open(SRCB, encoding="utf-8").read(), re.M)
    if not m:
        raise SystemExit("no covered marker for N = %d" % NGATE)
    return float(m.group(1))


HEAD = [
    "STATISTIC: the Type II sum split at the modulus, A = I(1) against",
    "           B = sum over d >= 2, their signs, the relative residue",
    "           |A+B|/|A| and its exponent, and the exponent of A.",
    "FIELD: N = %s; d over the squarefree d <= N/K coprime to N and m"
    % NS,
    "       over the squarefree m in [K, N/d) coprime to dN, with",
    "       K = floor(N^%.2f), the ranges of rem:bilinear. The covered"
    % THETA,
    "       part at N = %d is READ from results/audit_bilinear.txt at"
    % NGATE,
    "       full double precision.",
    "DERIVED: the demand |S| <~ l2 with l2 of order sqrt(N) times logs",
    "         is |A+B|/|A| <~ N^(-1/2); alpha puts the achieved ratio",
    "         at N^(alpha-1).",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    say("READ audit_bilinear.txt %d %.17e" % (NGATE, pub))
    say("  the covered part, at full double precision")
    say("PRINTBOUND audit_mainterm_removal %d %.20f" % (17, 5e-18))
    say("  theta %.2f, tolerance %.0e, alpha-1 %.6f, slope cap %.2f"
        % (THETA, RELID, ALPHAM1, SLOPECAP))
    say("  the demand's exponent for the relative residue: %.2f"
        % (-0.5))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    rows = []
    for N in NS:
        A, B, l2, K = split(N, lam, mu, sqf)
        S = -(A + B)
        rows.append((N, A, B, S, l2, K))
        say("  N = %-9d A %+14.2f  B %+14.2f  A+B %+12.2f  |S|/l2 "
            "%8.4f" % (N, A, B, A + B, abs(S) / l2))
        say("POINT mainA_%d %.6e" % (N, A))
        say("POINT mainres_%d %.6e" % (N, abs(A + B) / abs(A)))
    say("SCALES %d" % len(rows))

    # -------------------------------------------------------------- Z1
    say()
    say("Z1  the gate")
    g = [r for r in rows if r[0] == NGATE][0]
    rel = abs(g[3] - pub) / max(abs(pub), 1.0)
    z1 = rel <= RELID
    say("  -(A+B) here %.17e" % g[3])
    say("  against its       %.17e   relative %.2e" % (pub, rel))
    say("  Z1 %s   (cap: %.0e relative)"
        % ("hold" if z1 else "REFUTED", RELID))
    if not z1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- Z2
    say()
    say("Z2  do A and B fight?")
    bad = [N for N, A, B, _, _, _ in rows if A * B >= 0]
    z2 = not bad
    for N, A, B, _, _, _ in rows:
        say("  N = %-9d sign A %+d, sign B %+d, |B|/|A| %.6f"
            % (N, 1 if A > 0 else -1, 1 if B > 0 else -1,
               abs(B) / abs(A)))
    say("  Z2 %s   (cap: opposite at every N)"
        % ("hold" if z2 else "REFUTED"))

    x = np.array([math.log(r[0]) for r in rows])

    # -------------------------------------------------------------- Z3
    say()
    say("Z3  the relative residue's exponent")
    y = np.array([math.log(abs(r[1] + r[2]) / abs(r[1]))
                  for r in rows])
    b3, se3 = fit(x, y)
    z3 = abs(b3 - ALPHAM1) <= SLOPECAP
    say("  slope %+.6f +- %.6f against alpha-1 %.6f, gap %+.6f"
        % (b3, se3, ALPHAM1, b3 - ALPHAM1))
    say("  the demand wants %.2f, so the gap to the demand is %+.6f"
        % (-0.5, b3 - (-0.5)))
    say("TSTAT mainres_slope %.2f" % (b3 / se3))
    say("SPREAD mainres_slope %.6f" % se3)
    say("POINT mainres_exp %.6f" % b3)
    say("  Z3 %s   (cap: %.2f of alpha-1)"
        % ("hold" if z3 else "REFUTED", SLOPECAP))
    if se3 > SLOPECAP:
        say("  UNRESOLVED: the slope's own error exceeds the cap, so "
            "the verdict")
        say("  stands without a reading, as the rule says")

    # -------------------------------------------------------------- Z4
    say()
    say("Z4  is A a main term?")
    y4 = np.array([math.log(abs(r[1])) for r in rows])
    b4, se4 = fit(x, y4)
    z4 = abs(b4 - MAINEXP) <= SLOPECAP
    say("  exponent of |A| %+.6f +- %.6f against %.2f, gap %+.6f"
        % (b4, se4, MAINEXP, b4 - MAINEXP))
    say("TSTAT mainA_slope %.2f" % ((b4 - MAINEXP) / se4))
    say("SPREAD mainA_slope %.6f" % se4)
    say("POINT mainA_exp %.6f" % b4)
    say("  Z4 %s   (cap: %.2f of %.1f)"
        % ("hold" if z4 else "REFUTED", SLOPECAP, MAINEXP))
    if se4 > SLOPECAP:
        say("  UNRESOLVED: the slope's own error exceeds the cap, so "
            "the verdict")
        say("  stands without a reading, as the rule says")

    # a diagnostic, after the verdicts and predicted by nothing:
    # I(1) carries a log m weight, which this run's own formula shows
    # and its Z4 cap forgot.
    say()
    say("  a diagnostic, after the verdicts and predicted by nothing")
    say("  I(1) = sum_m mu^2(m) Lambda(N-dm) log m carries a log m "
        "weight, so a")
    say("  main term for it is N log N and not N. Z4's cap was "
        "written on N.")
    yl = np.array([math.log(n * math.log(n)) for n in NS])
    bl, sel = fit(x, yl)
    say("    N log N over these same %d N fits exponent %+.6f"
        % (len(NS), bl))
    say("    |A| fits %+.6f +- %.6f, so it stands %+.6f from N log N"
        % (b4, se4, b4 - bl))
    say("POINT nlogn_exp %.6f" % bl)
    say("POINT mainA_vs_nlogn %.6f" % (b4 - bl))
    say("  and the model-free statement this run makes needs no "
        "exponent at all:")
    say("    |S|/l2 runs %.4f to %.4f across the field, and the "
        "demand needs it"
        % (abs(rows[0][3]) / rows[0][4],
           abs(rows[-1][3]) / rows[-1][4]))
    say("    bounded")
    say("POINT sl2_first %.6f" % (abs(rows[0][3]) / rows[0][4]))
    say("POINT sl2_last %.6f" % (abs(rows[-1][3]) / rows[-1][4]))

    say()
    say("=" * 70)
    say("Z1 %s  Z2 %s  Z3 %s  Z4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (z1, z2, z3, z4)))
    say()
    if z2 and z4:
        say("item 5 is the removal of a main term to a relative "
            "accuracy. A is of")
        say("main-term order and B fights it at every N; what the "
            "demand asks is")
        say("that the leftover be a relative N^-0.5 of A, and what is "
            "achieved is")
        say("the exponent printed for Z3. that is the sharpest form "
            "this")
        say("requirement has taken here, and it is a restatement, not "
            "a bound.")
    elif not z2:
        say("A and B do not oppose at every N, so the split is not "
            "the structure")
        say("this run assumed and the restatement is not available.")
    else:
        say("A is not of main-term order, so calling the requirement "
            "a main-term")
        say("removal is the wrong description whatever the residue "
            "does.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
