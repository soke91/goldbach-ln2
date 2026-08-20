# -*- coding: utf-8 -*-
r"""
The target of the sign axis is derived, not fitted

WHAT IS AT STAKE

rem:deficitregion closed the computational branch of the sign axis
without an answer, and named what would open it: "a derivation that
says which family the deficit belongs to, and no measurement in this
repository has ever supplied one."  The deficit is a difference of two
exponents,

    deficit = e(l1/l2) - e(G),   currently +0.283586 - 0.149567,

and both have been treated as fitted quantities for six cycles.  One
of them need not be.

rem:leanidentity's W4 already wrote the ingredient down and drew the
wrong conclusion from it.  It observed that l1 <= sqrt(#k) * l2 --
Cauchy-Schwarz, an identity of norms that cannot be violated -- and
then reported the measured exponent +0.287798 as "3.15 standard errors
above" its ceiling of theta/2 = 0.28, calling W4 refuted.  A ratio
bounded by one cannot exceed its ceiling asymptotically.  The measured
ratio (l1/l2)/sqrt(#k) runs 0.6760 to 0.6909: **below one at every
point, and rising.**  An exponent above the ceiling on a bounded
rising ratio is a transient, and rem:headfraction made exactly this
argument about the head's share of |sum a| -- it was never applied
here.

And #k is not an empirical quantity either.  The k-range is the
squarefree k < N^theta coprime to N, so

    #k = (6/pi^2) * prod_{p | N} p/(p+1) * N^theta * (1 + o(1)),

leading order derived, no fit.  **theta = 0.56 is a parameter of the
construction -- the k-range -- and not the unknown level exponent**,
so the ceiling theta/2 = 0.28 is an exactly known number and not a
number that moves with what the program is trying to prove.

If this holds, the demand of rem:leanidentity stops being "e(G) must
reach a fitted +0.283586" and becomes **"e(G) must reach 0.28"**, the
deficit's target side is pinned by derivation, and what happens past
the field is a question about G alone rather than about a difference
of two free exponents.  That is one of the two terms rem:deficitregion
asked for.

This run measures the doubling family to the published field top,
1.024e8, one radical throughout, and asks whether the derivation
describes what is there.

BACKS: Remark {#rem:targetderived} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  THE GATE.  On the eight N of audit_lean_identity.py this
      reproduces its #k exactly and its l1/l2 to four decimals, and
      the ratio r = (l1/l2)/sqrt(#k) reproduces 0.6760 .. 0.6909.
  E2  THE DERIVATION CHECK.  #k agrees with its derived leading order
      (6/pi^2) prod_{p|N} p/(p+1) N^theta within 3 per cent at every
      N.  This tests the derivation, not the data: a failure means the
      count is not what the algebra says and the rest does not follow.
  E3  r <= 1 at every N.  Cauchy-Schwarz forces it; a violation is an
      implementation error in this script and nothing below is
      reported.
  E4  r rises with decaying increments: its total rise over the top
      half of the field is smaller than over the bottom half.  That is
      what a bounded rising ratio approaching a limit does.
  E5  **THE ONE THAT MATTERS.**  The local exponent of l1/l2, fitted
      on a sliding window, falls toward the ceiling: the top window's
      value is below the bottom window's, and closer to theta/2.
  E6  The budget is finite and small: the total remaining rise
      available to log(l1/l2) above (theta/2) log N is -log r at the
      top N, and it is below 0.5.

REFUTATION RULE (fixed before the run)

  E1  REFUTED outside four decimals or on any #k; nothing is reported.
  E2  REFUTED outside 3 per cent.  Then the derived count is wrong and
      the ceiling argument loses the half that is not Cauchy-Schwarz.
      The Cauchy-Schwarz half would survive -- r <= 1 is independent
      of any asymptotic for #k -- and the remark must say so rather
      than claim the whole.
  E3  REFUTED by any r > 1.  Implementation error, reported as such.
  E4  REFUTED if the top half's rise is the larger.  Then r is not
      visibly settling and the transient is not seen to decay inside
      the field, though it would still be bounded.
  E5  **REFUTED if the top window's exponent is not below the bottom
      window's.**  Then the approach to the ceiling is not visible in
      this field, and the derivation, while still true, would have no
      measured content here: the target would be pinned in the limit
      and unpinned everywhere the program can compute.  That is the
      outcome that costs the most and it must be stated plainly.
      **The unresolved case is named**: a sliding-window exponent on
      ten points carries an error near the 0.034 that rem:alphalocal
      measured for octave fits, so a fall smaller than that is not
      distinguishable from no fall, and E5 is to be read as refuted
      unless the fall exceeds the printed window error.
  E6  REFUTED at or above 0.5.  Then the transient's budget is not
      small next to the deficit's +0.134019 and pinning the target
      buys less than this run claims.

  WHAT THIS CANNOT DO.  Ten N on one radical.  The derived count's
  arithmetic factor is constant along a doubling family, so this run
  cannot see it vary and does not test that half of E2 -- it tests
  only the power.  Nothing here measures G, and nothing here says
  where e(G) goes.
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
OUT = os.path.join(ROOT, "results", "audit_target_derived.txt")
SRC = os.path.join(ROOT, "results", "audit_lean_identity.txt")

THETA = 0.56
NS = [200_000 * (1 << j) for j in range(10)]
WIN = 5
DEC = 4
TOLPC = 3.0
BUDGETCAP = 0.5
SIXPI2 = 6.0 / math.pi ** 2


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


def krange(N, sqf):
    """the squarefree k < N^theta coprime to N, as lean_identity has"""
    PN = factor_set(N)
    K = int(N ** THETA)
    keep = sqf[2:K].copy()
    idx = np.arange(2, K, dtype=np.int64)
    for q in PN:
        keep &= (idx % int(q)) != 0
    return idx[keep]


def norms(N, lam, mu, ks):
    """l1, l2 and |sum a| for a_k = (log k) H(N;k)"""
    s1 = s2 = ss = 0.0
    for k in ks:
        k = int(k)
        M = (N - 1) // k
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % int(q) != 0]
        h = float((lam[N - ms * k]
                   * mu[ms].astype(np.float64)).sum())
        del ms
        a = math.log(k) * h
        s1 += abs(a)
        s2 += a * a
        ss += a
    return s1, math.sqrt(s2), abs(ss)


def fit(x, y):
    n = len(x)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)[0], np.polyfit(x, y, 1)[1]
    r = y - (b * x + a0)
    if n > 2:
        se = math.sqrt(float((r ** 2).sum() / (n - 2))
                       / float(((x - x.mean()) ** 2).sum()))
    else:
        se = float("nan")
    return float(b), se


def read_pub():
    src = io.open(SRC, encoding="utf-8").read()
    out = {}
    for m in re.finditer(
            r"N = (\d+)\s+#k = (\d+)\s+G = ([\d.]+)\s+"
            r"l1/l2 = ([\d.]+)", src):
        out[int(m.group(1))] = (int(m.group(2)), float(m.group(4)))
    return out


HEAD = [
    "STATISTIC: the ratio r = (l1/l2)/sqrt(#k) for the sign axis'",
    "           a_k = (log k)H(N;k), its rise, the local exponent of",
    "           l1/l2 on sliding windows against the derived ceiling",
    "           theta/2, and the remaining budget -log r.",
    "FIELD: N = 200000*2^j for j < %d, one radical throughout; k over"
    % len(NS),
    "       the squarefree k < N^%.2f coprime to N, the k-range of"
    % THETA,
    "       code/audit_lean_identity.py. #k, l1/l2 for the eight N",
    "       that script published are READ from",
    "       results/audit_lean_identity.txt as the gate.",
    "DERIVED: #k = (6/pi^2) prod_{p|N} p/(p+1) N^theta (1+o(1)); the",
    "         ceiling theta/2 = %.2f is a fact about the k-range and"
    % (THETA / 2),
    "         not about the level exponent the program is proving.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    for N in sorted(pub):
        say("READ audit_lean_identity.txt %d %.4f" % (N, pub[N][1]))
    say("  l1/l2 at the N that script published, its gate here")
    say("PRINTBOUND audit_target_derived %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  window %d, tolerance %.1f per cent, budget cap %.1f"
        % (WIN, TOLPC, BUDGETCAP))
    say("  theta %.2f, ceiling %.2f, 6/pi^2 %.6f"
        % (THETA, THETA / 2, SIXPI2))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    rows = []
    for N in NS:
        ks = krange(N, sqf)
        l1, l2, sa = norms(N, lam, mu, ks)
        c = SIXPI2
        for p in sorted(factor_set(N)):
            c *= p / (p + 1.0)
        pred = c * N ** THETA
        r = (l1 / l2) / math.sqrt(ks.size)
        rows.append((N, ks.size, l1, l2, sa, pred, r))
        say("  N = %-11d #k = %-7d l1/l2 = %10.4f  r = %.4f"
            % (N, ks.size, l1 / l2, r))
        say("POINT targetratio_%d %.6f" % (N, r))
    say("SCALES %d" % len(rows))

    x = np.array([math.log(t[0]) for t in rows])

    # -------------------------------------------------------------- E1
    say()
    say("E1  the gate: does this reproduce the published rows?")
    e1 = True
    for N, nk, l1, l2, sa, pred, r in rows:
        if N not in pub:
            continue
        pk, pl = pub[N]
        g = (nk == pk
             and abs(round(l1 / l2, DEC) - round(pl, DEC))
             < 10.0 ** (-DEC) / 2)
        e1 &= g
        say("  N = %-11d #k %d/%d  l1/l2 %.4f/%.4f  %s"
            % (N, nk, pk, l1 / l2, pl, "ok" if g else "MISMATCH"))
    say("  E1 %s   (cap: #k exact, l1/l2 to %d decimals)"
        % ("hold" if e1 else "REFUTED", DEC))
    if not e1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- E2
    say()
    say("E2  does #k match its derived leading order?")
    say("    N            #k        derived      per cent")
    worst = 0.0
    for N, nk, l1, l2, sa, pred, r in rows:
        d = 100.0 * (nk - pred) / pred
        worst = max(worst, abs(d))
        say("  %-12d %-9d %12.1f %+11.3f" % (N, nk, pred, d))
    e2 = worst <= TOLPC
    say("POINT derivedcount_worst %.4f" % worst)
    say("  worst %.3f per cent" % worst)
    say("  E2 %s   (cap: %.1f per cent)"
        % ("hold" if e2 else "REFUTED", TOLPC))

    # -------------------------------------------------------------- E3
    say()
    say("E3  is r <= 1 everywhere, as Cauchy-Schwarz forces?")
    rs = [t[6] for t in rows]
    e3 = max(rs) <= 1.0
    say("  r runs %.4f to %.4f" % (min(rs), max(rs)))
    say("  E3 %s   (cap: 1)" % ("hold" if e3 else "REFUTED"))
    if not e3:
        say("  that is an implementation error in this script and "
            "nothing")
        say("  below is reported, as the rule says")
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- E4
    say()
    say("E4  does r rise with decaying increments?")
    h = len(rows) // 2
    lo = rs[h] - rs[0]
    hi = rs[-1] - rs[h]
    e4 = hi < lo
    say("  bottom half rise %+.5f, top half rise %+.5f" % (lo, hi))
    say("POINT riselow %.5f" % lo)
    say("POINT risehigh %.5f" % hi)
    say("  E4 %s   (cap: top half the smaller)"
        % ("hold" if e4 else "REFUTED"))

    # -------------------------------------------------------------- E5
    say()
    say("E5  does the local exponent of l1/l2 fall toward %.2f?"
        % (THETA / 2))
    y = np.array([math.log(t[2] / t[3]) for t in rows])
    say("    window            slope        s.e.   above ceiling")
    wins = []
    for i in range(0, len(rows) - WIN + 1):
        b, se = fit(x[i:i + WIN], y[i:i + WIN])
        wins.append((b, se))
        say("  %5.2f - %5.2f   %+.6f   %.6f   %+.6f"
            % (math.log10(rows[i][0]),
               math.log10(rows[i + WIN - 1][0]),
               b, se, b - THETA / 2))
    b0, se0 = wins[0]
    b1, se1 = wins[-1]
    drop = b0 - b1
    e5 = (b1 < b0) and (drop > max(se0, se1))
    say("TSTAT targetdrop %.2f" % (drop / max(se0, se1)))
    say("SPREAD targetdrop %.6f" % max(se0, se1))
    say("  fall %+.6f against the larger window error %.6f"
        % (drop, max(se0, se1)))
    say("  E5 %s   (cap: a fall larger than the window error)"
        % ("hold" if e5 else "REFUTED"))

    # -------------------------------------------------------------- E6
    say()
    say("E6  how much rise is left to log(l1/l2) above the ceiling?")
    budget = -math.log(rs[-1])
    e6 = budget < BUDGETCAP
    say("  r at the top N is %.4f, so the budget is %.4f"
        % (rs[-1], budget))
    say("POINT budget %.5f" % budget)
    say("  E6 %s   (cap: %.1f)"
        % ("hold" if e6 else "REFUTED", BUDGETCAP))

    say()
    say("=" * 70)
    say("E1 %s  E2 %s  E3 %s  E4 %s  E5 %s  E6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (e1, e2, e3, e4, e5, e6)))
    say()
    if e2 and e5:
        say("the target side of the deficit is derived. e(l1/l2) is "
            "not a free")
        say("exponent: it is %.2f plus a transient of a ratio bounded "
            "by one," % (THETA / 2))
        say("the ceiling is a fact about the k-range rather than "
            "about the level")
        say("the program is proving, and the fall toward it is "
            "visible here.")
        say("rem:leanidentity's demand becomes e(G) -> %.2f and what "
            "happens past" % (THETA / 2))
        say("the field is a question about G alone.")
    elif e2:
        say("the derivation stands and its content is not visible in "
            "this field.")
        say("r is bounded by one and #k is the derived count, so "
            "e(l1/l2) -> %.2f" % (THETA / 2))
        say("in the limit; but the fall is not resolved on these ten "
            "points, so the")
        say("target is pinned in the limit and unpinned everywhere "
            "this program can")
        say("compute. That is the outcome E5's rule said would cost "
            "the most.")
    else:
        say("the derived count does not describe #k here, so the half "
            "of the")
        say("ceiling argument that is not Cauchy-Schwarz does not "
            "stand. r <= 1")
        say("survives -- it needs no asymptotic -- and nothing "
            "stronger is claimed.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
