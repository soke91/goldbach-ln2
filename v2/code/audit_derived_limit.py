# -*- coding: utf-8 -*-
r"""
A derived limit for the corrections, with no free parameter

WHAT IS AT STAKE

rem:driftsigns established the one window-free fact this branch has
about the weight corrections: c(3, N) falls at all ten doublings, a
two-sided coin tail of 0.001953 against four primes named in advance.
It could say the direction and not the destination, and asked where
c(3) goes.

Writing this run answered that from the construction instead of from
a fit.  Count the densities in I(p)/I(1):

  * m runs over [K, N/p) instead of [K, N), a factor 1/p;
  * m must be coprime to p, and among squarefree m that costs a
    factor p/(p+1);
  * N - pm is confined to the class N mod p, where primes sit at
    p/(p-1) times their unrestricted density.

So I(p)/I(1) -> (1/p)(p/(p+1))(p/(p-1)) = p/(p^2-1), and dividing by
the model's own 1/(p-1),

    c(p)  ->  p / (p + 1).

**No parameter is fitted anywhere in that.**  It is the first derived
target this branch has had for the corrections, and it is what a
description needs if it is to be more than a curve.

The prediction is checkable on the eleven N already used: writing
g(p, N) = c(p, N) (p+1)/p, the claim is g -> 1.

BACKS: Remark {#rem:derivedlimit} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  JJ1 THE GATE.  A at N = 200000 reproduces rem:maintermremoval's
      POINT mainA marker to a relative 1e-12, and the eleven-N drift
      of c(3) reproduces rem:driftpower's -0.009622 to six decimals.
  JJ2 **THE DERIVED TARGET.**  At the largest N, g(p) is within 0.05
      of 1 for p = 7, 11 and 13.
  JJ3 **And c(3) is above its target and coming down to it**: g(3)
      exceeds 1 at every one of the eleven N, and is smaller at the
      last than at the first.
  JJ4 The approach is general, not three lucky primes: over every
      prime that contributes at all eleven N, the median |g - 1| at
      the largest N is under 0.10.

REFUTATION RULE (fixed before the run)

  JJ1 REFUTED outside either tolerance; nothing below is reported.
  JJ2 **REFUTED outside 0.05 on any of the three.**  Then p/(p+1) is
      not the limit and the density count above is wrong somewhere --
      most likely in the class-density factor, which is the only one
      of the three that is not elementary counting.  The remark must
      then say the derivation failed and not adjust it: a factor
      fitted after the fact is the thing this whole branch keeps
      refusing to do.
  JJ3 REFUTED if g(3) dips below 1 at any N, or if it is not smaller
      at the last than at the first.  A g(3) below 1 would mean c(3)
      has passed its derived limit and is still falling, which no
      limit argument survives.
  JJ4 REFUTED above 0.10.  Then the target holds at the three primes
      chosen and not as a rule, which is a weaker claim and must be
      reported as one.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  rem:leveldense
  measured a scatter of about 0.03 in a related weight quantity and
  rem:driftsigns' per-step differences run to 0.016, so **an
  agreement inside 0.05 is agreement inside roughly one scatter and
  is not a sharp confirmation**.  This run prints |g - 1| for every
  prime so the distribution is visible rather than summarised, and
  JJ2's cap is to be read as "not excluded at the resolution this
  field has" rather than as "confirmed".  A derivation that survives
  is not thereby proved.

  WHAT THIS CANNOT DO.  One radical family, eleven N, and a heuristic
  density count -- the class-density factor p/(p-1) is a statement
  about primes in progressions on average and is not proved for the
  particular ranges here.  Nothing in this run bounds |sum a| or
  moves item 5's demand; it describes a correction to a model
  rem:residuemodel already refuted as a whole.
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
OUT = os.path.join(ROOT, "results", "audit_derived_limit.txt")
SRCM = os.path.join(ROOT, "results", "audit_mainterm_removal.txt")
SRCP = os.path.join(ROOT, "results", "audit_drift_power.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000, 6_400_000, 12_800_000, 25_600_000]
NGATE = 200_000
RELID = 1e-12
DEC = 6
TARGET = (7, 11, 13)
NEAR = 0.05
MEDCAP = 0.10
NOISEFRAC = 0.001


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


def phi(n):
    r = n
    for p in factor_set(n):
        r = r // p * (p - 1)
    return r


def pieces(N, lam, mu, sqf):
    PN = factor_set(N)
    K = int(N ** THETA)
    D = (N - 1) // K
    out = {}
    for d in range(1, D + 1):
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
        out[d] = (md, float((lam[N - d * ms]
                             * np.log(ms.astype(np.float64))).sum()))
        del ms, keep
    return out, D


def model(d, D):
    return (1.0 / phi(d)) * (1.0 - d / D) / (1.0 - 1.0 / D)


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se


def read_pub():
    m = re.search(r"^POINT mainA_%d (\S+)\s*$" % NGATE,
                  io.open(SRCM, encoding="utf-8").read(), re.M)
    d = re.search(r"^POINT drift11_3 (\S+)\s*$",
                  io.open(SRCP, encoding="utf-8").read(), re.M)
    if not m or not d:
        raise SystemExit("a published value is missing")
    return m.group(1), d.group(1)


HEAD = [
    "STATISTIC: g(p, N) = c(p, N) (p+1)/p, the weight correction",
    "           against the limit p/(p+1) that the density count",
    "           derives, at eleven N and over every prime the field",
    "           supports.",
    "FIELD: N = %s; d over the squarefree d <= D coprime to N, with"
    % NS,
    "       D = floor((N-1)/K), K = floor(N^%.2f). A at N = %d and"
    % (THETA, NGATE),
    "       the eleven-N drift of c(3) are READ from",
    "       results/audit_mainterm_removal.txt and",
    "       results/audit_drift_power.txt.",
    "DERIVED: I(p)/I(1) -> (1/p)(p/(p+1))(p/(p-1)), from the m-range,",
    "         the coprimality of m to p among squarefree m, and the",
    "         class density of primes in N mod p; divided by the",
    "         model's 1/(p-1) this gives c(p) -> p/(p+1). No",
    "         parameter is fitted.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    rawA, rawd3 = read_pub()
    pubA, pubd3 = float(rawA), float(rawd3)
    say("READ audit_mainterm_removal.txt POINT mainA_%d %s"
        % (NGATE, rawA))
    say("READ audit_drift_power.txt POINT drift11_3 %s" % rawd3)
    say("PRINTBOUND audit_derived_limit 17 5e-18")
    say("PRINTBOUND FOR audit_derived_limit_table %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  the finest printing any tolerance here judges is the")
    say("  seventeen-decimal A marker; the six-decimal drift marker")
    say("  is judged at its own last place.")
    say("  theta %.2f, target primes %s, near %.2f, median cap %.2f,"
        % (THETA, TARGET, NEAR, MEDCAP))
    say("  noise floor %.3f of |w(3)|" % NOISEFRAC)
    say("  the derived limits: %s"
        % ", ".join("c(%d) -> %.6f" % (p, p / (p + 1.0))
                    for p in (3,) + TARGET))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    cs, A200 = {}, None
    for N in NS:
        ps, D = pieces(N, lam, mu, sqf)
        A = ps[1][1]
        if N == NGATE:
            A200 = A
        ws = {d: v / A for d, (md, v) in ps.items() if d != 1}
        cs[N] = {d: (ws[d] / model(d, D), abs(ws[d]))
                 for d in ws if d != D}
        say("  N = %-10d D = %-5d %d contributing d" % (N, D, len(ws)))
    say("SCALES %d" % len(NS))

    # ------------------------------------------------------------- JJ1
    x = np.array([math.log(N) for N in NS])
    b3, _ = fit(x, [cs[N][3][0] for N in NS])
    say()
    say("JJ1  the gate")
    ra = abs(A200 - pubA) / max(abs(pubA), 1.0)
    ok3 = abs(round(b3, DEC) - round(pubd3, DEC)) < 10.0 ** (-DEC)
    jj1 = ra <= RELID and ok3
    say("  A relative %.2e against %.0e   %s"
        % (ra, RELID, "ok" if ra <= RELID else "MISMATCH"))
    say("  c(3) eleven-N drift %+.6f against its %+.6f  %s"
        % (b3, pubd3, "ok" if ok3 else "MISMATCH"))
    say("  JJ1 %s" % ("hold" if jj1 else "REFUTED"))
    if not jj1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    NTOP = NS[-1]
    say()
    say("      p   target      c(p) at the top N     g = c(p+1)/p"
        "     |g-1|")
    for p in (3,) + TARGET:
        c = cs[NTOP][p][0]
        g = c * (p + 1.0) / p
        say("  %5d   %.6f    %.6f            %.6f      %.6f"
            % (p, p / (p + 1.0), c, g, abs(g - 1.0)))
        say("POINT gtop_%d %.6f" % (p, g))

    # ------------------------------------------------------------- JJ2
    say()
    say("JJ2  is the derived target hit at 7, 11 and 13?")
    jj2 = True
    for p in TARGET:
        g = cs[NTOP][p][0] * (p + 1.0) / p
        ok = abs(g - 1.0) <= NEAR
        jj2 &= ok
        say("  p = %-3d |g - 1| = %.6f   %s"
            % (p, abs(g - 1.0), "ok" if ok else "OUT"))
    say("  JJ2 %s   (cap: %.2f)"
        % ("hold" if jj2 else "REFUTED", NEAR))

    # ------------------------------------------------------------- JJ3
    say()
    say("JJ3  is c(3) above its target and coming down?")
    g3 = [cs[N][3][0] * 4.0 / 3.0 for N in NS]
    above = all(v > 1.0 for v in g3)
    closer = abs(g3[-1] - 1.0) < abs(g3[0] - 1.0)
    jj3 = above and closer
    say("  g(3) %s" % " ".join("%.4f" % v for v in g3))
    say("  above one at every N: %s; closer at the last: %s"
        % (above, closer))
    say("POINT g3first %.6f" % g3[0])
    say("POINT g3last %.6f" % g3[-1])
    say("  JJ3 %s" % ("hold" if jj3 else "REFUTED"))

    # ------------------------------------------------------------- JJ4
    say()
    say("JJ4  does it hold as a rule, not at three primes?")
    common = None
    for N in NS:
        here = {p for p, cw in cs[N].items()
                if len(factor_set(p)) == 1
                and cw[1] >= NOISEFRAC * abs(cs[N][3][1])}
        common = here if common is None else (common & here)
    gs = [(p, cs[NTOP][p][0] * (p + 1.0) / p) for p in sorted(common)]
    devs = sorted(abs(g - 1.0) for _, g in gs)
    med = float(np.median(devs))
    jj4 = med <= MEDCAP
    say("  %d primes contribute at all eleven N" % len(gs))
    say("      p    g       |g-1|")
    for p, g in gs:
        say("  %5d  %.4f  %.4f" % (p, g, abs(g - 1.0)))
    say("  median |g-1| %.6f, worst %.6f" % (med, devs[-1]))
    say("POINT gmedian %.6f" % med)
    say("  JJ4 %s   (cap: %.2f)"
        % ("hold" if jj4 else "REFUTED", MEDCAP))

    say()
    say("=" * 70)
    say("JJ1 %s  JJ2 %s  JJ3 %s  JJ4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (jj1, jj2, jj3, jj4)))
    say()
    if jj2 and jj3:
        say("the corrections have a derived limit and the field does "
            "not exclude it.")
        say("c(p) -> p/(p+1) comes from counting three densities and "
            "fits no")
        say("parameter, c(3) sits above it and descends toward it, "
            "and the others")
        say("sit within one scatter of it. that is the first derived "
            "target this")
        say("branch has had for these weights -- and a derivation "
            "that survives is")
        say("not thereby proved.")
    elif not jj2:
        say("p/(p+1) is not the limit. the density count is wrong "
            "somewhere, most")
        say("likely in the class-density factor, and no factor is "
            "fitted after the")
        say("fact to rescue it.")
    else:
        say("c(3) has passed its derived limit or is not descending, "
            "which no limit")
        say("argument survives, and the target is withdrawn.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
