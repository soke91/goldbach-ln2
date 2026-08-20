# -*- coding: utf-8 -*-
r"""
Which way the residual in p resolves, decided by a sign

WHAT IS AT STAKE

rem:derivedlimit put a limit under the corrections that nothing fitted
-- c(p) -> p/(p+1) -- and JJ2 and JJ3 held at it while JJ4 broke on the
shape of its own failure: over twenty primes g = c(p+1)/p descends
monotonically from 1.154331 to 0.8609, crossing 1 between 11 and 13.
That remark named two ways out and could not choose between them.

  (a) The approach is in N/p, not N.  Then the large primes are simply
      far from their own asymptotics and g(13) must climb to 1 as N
      grows.
  (b) p/(p+1) is the limit for small p and not a law.  Then g(13) keeps
      falling.

**(a) is falsifiable by a sign alone and needs no rate.**  If
g(p, N) = G(N/p) for a single function G, then G is fixed by the
p-direction: g falls as p rises at fixed N, so G rises in its
argument, so **g must rise as N rises at fixed p.**  rem:driftsigns
measured c(3) falling at all ten doublings and rem:derivedlimit
measured g(3) falling from 1.243394 to 1.154331.  A collapse in N/p
predicts the opposite sign of the one fact this branch has that no
window can bend.

This run states that as a rule before looking, counts the same signs
for the primes below 1 where the question actually lives, and puts the
two directions on the same scale as finite differences -- no line is
fitted anywhere, because rem:driftpower showed slopes here are the
window's.

BACKS: Remark {#rem:limitdirection} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  KK1 THE GATE.  g(3) at the eleven N reproduces rem:derivedlimit's
      POINT g3first and POINT g3last to six decimals.
  KK2 **THE PRIMES BELOW ONE MOVE AWAY.**  For p = 13, 17 and 19,
      g(p, N) falls at 8 or more of the ten doublings.  These sit
      below 1, so falling is moving away from the derived limit and
      is (b).
  KK3 **THE COLLAPSE IS REFUTED BY SIGN.**  At the top N, the change
      of g per doubling of N and the change of g per halving of p
      have **opposite** signs, which no single G(N/p) permits.
  KK4 And the mismatch is not marginal: the p-direction moves at
      least 5 times as fast per octave as the N-direction.
  KK5 THE d=1 ASYMMETRY IS A p-INDEPENDENT SHIFT.  I(1) alone drops
      prime-power m while I(p) does not, a subset of relative size
      about 1/log N.  Recomputing every g with that exclusion lifted,
      the shift is nearly the same for all primes: the spread of the
      twenty shifts is under a fifth of their mean magnitude.

REFUTATION RULE (fixed before the run)

  KK1 REFUTED outside six decimals; nothing below is reported.
  KK2 REFUTED if any of the three falls at 5 or fewer of the ten
      doublings.  That prime is then climbing toward 1 and (a) is
      alive for it.  Three primes are named in advance, so the
      adjusted two-sided bound is 0.05/3 = 0.016667 and an 8/10 --
      tail 0.109375 -- **does not clear it**; only a 9/10 or 10/10
      does.  A KK2 that holds at exactly 8/10 is a lean and must be
      written as one.
  KK3 REFUTED if the two signs agree.  Then a collapse in N/p is not
      excluded and (a) survives; the run then says so and the choice
      stays open.
  KK4 REFUTED below 5.  Then the two directions are comparable and
      the sign disagreement could be a small quantity crossing zero
      rather than a structural mismatch.
  KK5 REFUTED if the spread of the shifts exceeds a fifth of their
      mean magnitude.  Then the d=1 asymmetry carries p-dependence
      of its own and is a candidate for the residual rather than a
      common offset.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  KK2 and KK3 can
  both hold and still leave (b) unproved: rem:leveldense's scatter is
  0.03 and the whole p-descent of g spans 0.293431 from 1.154331 to
  0.860900, so **only about ten scatters separate the two ends of the
  effect being explained** -- and KK2's 8/10 case has a tail of
  0.109375 against an adjusted bound of 0.016667, which is not
  evidence.  If KK2 lands at 8/10 for any of the three, that prime is
  undecided and the remark must say so rather than pool it with the
  others.

  WHAT THIS CANNOT DO.  Refuting the N/p collapse says the residual
  is not a finite-N artefact of that one form; it does not say what
  the residual is.  Nothing here bounds |sum a|, and KK5 tests the
  model's normalisation, not the object -- lifting the exclusion
  changes what I(1) counts and is a diagnostic, never a redefinition.
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
OUT = os.path.join(ROOT, "results", "audit_limit_direction.txt")
SRCD = os.path.join(ROOT, "results", "audit_derived_limit.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000, 6_400_000, 12_800_000, 25_600_000]
DEC = 6
WATCH = (13, 17, 19)
MOVEMIN = 8
LEAN = 9
RATIOMIN = 5.0
SPREADFRAC = 0.2
PLO, PHI = 3, 7


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
    """Both normalisations at once: A drops prime-power m as the
    object does, Araw keeps them.  Only d = 1 differs."""
    PN = factor_set(N)
    K = int(N ** THETA)
    D = (N - 1) // K
    out, araw = {}, None
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
            lg = np.log(ms[keep].astype(np.float64))
            araw = float((lam[N - ms[keep]] * lg).sum())
            keep &= lam[ms] == 0.0
        ms = ms[keep]
        if ms.size == 0:
            continue
        out[d] = (md, float((lam[N - d * ms]
                             * np.log(ms.astype(np.float64))).sum()))
        del ms, keep
    return out, D, araw


def model(d, D):
    return (1.0 / phi(d)) * (1.0 - d / D) / (1.0 - 1.0 / D)


def read_pub():
    src = io.open(SRCD, encoding="utf-8").read()
    a = re.search(r"^POINT g3first (\S+)\s*$", src, re.M)
    b = re.search(r"^POINT g3last (\S+)\s*$", src, re.M)
    if not a or not b:
        raise SystemExit("a published value is missing")
    return a.group(1), b.group(1)


def tail(k, tot):
    hi = max(k, tot - k)
    return 2.0 * sum(math.comb(tot, j)
                     for j in range(hi, tot + 1)) / 2.0 ** tot


HEAD = [
    "STATISTIC: the sign of the motion of g(p, N) = c(p, N)(p+1)/p,",
    "           in N at fixed p and in p at fixed N, tested against",
    "           the one thing a collapse g = G(N/p) would force.",
    "FIELD: N = %s;" % NS,
    "       d over the squarefree d <= D coprime to N, with",
    "       D = floor((N-1)/K), K = floor(N^%.2f). g(3) at the first"
    % THETA,
    "       and last N are READ from results/audit_derived_limit.txt.",
    "DERIVED: if g(p, N) = G(N/p) then G rises, because g falls in p;",
    "         so g must rise in N. It does not. No rate is needed and",
    "         no line is fitted -- every quantity below is a finite",
    "         difference.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    rawf, rawl = read_pub()
    say("READ audit_derived_limit.txt POINT g3first %s" % rawf)
    say("READ audit_derived_limit.txt POINT g3last %s" % rawl)
    say("PRINTBOUND audit_limit_direction %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, watched primes %s, moves at least %d of ten,"
        % (THETA, WATCH, MOVEMIN))
    say("  a lean below %d, speed ratio at least %.1f, spread under"
        % (LEAN, RATIOMIN))
    say("  %.1f of the mean; the two primes bracketing the"
        % SPREADFRAC)
    say("  p-direction are %d and %d" % (PLO, PHI))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    G, GR = {}, {}
    for N in NS:
        ps, D, araw = pieces(N, lam, mu, sqf)
        A = ps[1][1]
        G[N] = {d: (ps[d][1] / A) / model(d, D) * (d + 1.0) / d
                for d in ps if d != 1 and d != D}
        GR[N] = {d: (ps[d][1] / araw) / model(d, D) * (d + 1.0) / d
                 for d in ps if d != 1 and d != D}
        G[N]["_raw"] = araw / A
        say("  N = %-10d D = %-5d A/Araw = %.6f"
            % (N, D, 1.0 / (araw / A)))
    say("SCALES %d" % len(NS))

    # ------------------------------------------------------------- KK1
    say()
    say("KK1  the gate")
    g3 = [G[N][3] for N in NS]
    okf = abs(round(g3[0], DEC) - float(rawf)) < 10.0 ** (-DEC)
    okl = abs(round(g3[-1], DEC) - float(rawl)) < 10.0 ** (-DEC)
    kk1 = okf and okl
    say("  g(3) first %.6f against its %s   %s"
        % (g3[0], rawf, "ok" if okf else "MISMATCH"))
    say("  g(3) last  %.6f against its %s   %s"
        % (g3[-1], rawl, "ok" if okl else "MISMATCH"))
    say("  KK1 %s" % ("hold" if kk1 else "REFUTED"))
    if not kk1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------------------- KK2
    say()
    say("KK2  do the primes below one move away from it?")
    say("      p   g first   g last    falls   tail    verdict")
    kk2, leaned = True, []
    for p in WATCH:
        vs = [G[N][p] for N in NS]
        dn = sum(1 for i in range(len(vs) - 1) if vs[i + 1] < vs[i])
        tot = len(vs) - 1
        ok = dn >= MOVEMIN
        kk2 &= ok
        if ok and dn < LEAN:
            leaned.append(p)
        say("  %5d   %.6f  %.6f  %2d/%d  %.6f  %s"
            % (p, vs[0], vs[-1], dn, tot, tail(dn, tot),
               "ok" if ok else "OUT"))
        say("POINT falls_%d %d" % (p, dn))
        say("POINT glast_%d %.6f" % (p, vs[-1]))
    say("  adjusted two-sided bound %.6f over %d named primes"
        % (0.05 / len(WATCH), len(WATCH)))
    say("  KK2 %s%s"
        % ("hold" if kk2 else "REFUTED",
           "" if not leaned else
           "  -- but %s land at %d of ten, a lean and not evidence"
           % (leaned, MOVEMIN)))

    # ------------------------------------------------------------- KK3
    say()
    say("KK3  can one G(N/p) carry both directions?")
    NTOP = NS[-1]
    dN = (g3[-1] - g3[0]) / (len(NS) - 1)
    dp = ((G[NTOP][PLO] - G[NTOP][PHI])
          / (math.log(PHI / float(PLO)) / math.log(2.0)))
    say("  per doubling of N at p = %d      %+.6f" % (PLO, dN))
    say("  per halving  of p at the top N   %+.6f" % dp)
    say("POINT perN %+.6f" % dN)
    say("POINT perP %+.6f" % dp)
    kk3 = (dN < 0.0) != (dp < 0.0)
    say("  a collapse needs the same sign in both;  they %s"
        % ("differ" if kk3 else "AGREE"))
    say("  KK3 %s" % ("hold" if kk3 else "REFUTED"))

    # ------------------------------------------------------------- KK4
    say()
    say("KK4  is the mismatch structural or a crossing of zero?")
    ratio = abs(dp) / abs(dN) if dN else float("inf")
    kk4 = ratio >= RATIOMIN
    say("  |per p| / |per N| = %.6f" % ratio)
    say("POINT speedratio %.6f" % ratio)
    say("  KK4 %s   (floor: %.1f)"
        % ("hold" if kk4 else "REFUTED", RATIOMIN))

    # ------------------------------------------------------------- KK5
    say()
    say("KK5  is the d=1 asymmetry a common offset?")
    common = None
    for N in NS:
        here = {p for p in G[N]
                if isinstance(p, int) and len(factor_set(p)) == 1}
        common = here if common is None else (common & here)
    common = sorted(p for p in common
                    if all(abs(G[N][p]) > 0 for N in NS))
    sh = [(p, GR[NTOP][p] - G[NTOP][p]) for p in common]
    vals = np.array([v for _, v in sh])
    mean, sd = float(abs(vals.mean())), float(vals.std(ddof=1))
    kk5 = sd <= SPREADFRAC * mean
    say("  A/Araw at the top N = %.6f, 1/log N = %.6f"
        % (1.0 / G[NTOP]["_raw"], 1.0 / math.log(NTOP)))
    say("  %d primes; shift mean %+.6f, spread %.6f, spread/mean %.6f"
        % (len(sh), vals.mean(), sd, sd / mean))
    say("POINT shiftmean %+.6f" % vals.mean())
    say("POINT shiftspread %.6f" % sd)
    say("  KK5 %s   (cap: %.1f of the mean)"
        % ("hold" if kk5 else "REFUTED", SPREADFRAC))
    say("      p    g       g without the exclusion   shift")
    for p, v in sh:
        say("  %5d  %.6f  %.6f  %+.6f" % (p, G[NTOP][p], GR[NTOP][p], v))

    say()
    say("=" * 70)
    say("KK1 %s  KK2 %s  KK3 %s  KK4 %s  KK5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (kk1, kk2, kk3, kk4, kk5)))
    say()
    if kk3 and kk4:
        say("the residual in p is not a finite-N artefact of the N/p "
            "form. g falls")
        say("in p and falls in N, and one function of N/p cannot do "
            "both, by a")
        say("factor of many rather than by a hair. what supplies the "
            "p-shape is")
        say("still unnamed -- refuting one form is not producing "
            "another.")
    else:
        say("a collapse in N/p is not excluded and the choice stays "
            "open, which is")
        say("what the rule said this outcome would mean.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
