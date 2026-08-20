# -*- coding: utf-8 -*-
r"""
The split that was always there: how many hits, and how heavy their log

WHAT IS AT STAKE

rem:conditionedlimit measured, above a floor that eq:condfloor demands
rather than chooses, that g descends in p with zero inversions in nine
steps and that its centre climbs by 0.006052 across three doublings
with 0.039116 still to go.  Two motions, both tiny, both unexplained.

They need not be explained by a fit.  I(p) is a sum of
Lambda(N - pm) log m over m in [K, N/p), so it factors **exactly**:

    W(p) = sum_m Lambda(N - pm)          the prime-hit mass
    L(p) = I(p) / W(p)                   the weighted mean of log m
    I(p) = W(p) L(p)

and therefore g = gcount * hlog with

    gcount(p) = (W(p)/W(1)) / w_model(p) * (p+1)/p
    hlog(p)   = L(p) / L(1) .

**rem:derivedlimit's density count was a prediction about W alone.**
It counted the m-range, the coprimality of m to p and the class
density of primes in N mod p, and every one of those is a statement
about how many hits there are.  The second factor was never counted,
and its size is forced by the construction: m runs to N/p against N,
so L(p) is about log(N/p) - 1 against log N - 1, which **falls in p
and rises to 1 in N** -- the two motions, in the right directions,
with no parameter.

And the level is forced too.  rem:limitdirection found I(1) alone
drops prime-power m, that A/A_raw climbs from 0.763851 to 0.861216,
and that lifting the exclusion moves every g by roughly a common
amount.  If W carries only the densities, gcount should sit at the
reciprocal of that ratio and be flat in p.

BACKS: Remark {#rem:logweightsplit} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  NN1 THE GATE.  At N = 25600000, g reproduces rem:derivedlimit's
      POINT gtop_3 and gtop_11 and rem:limitdirection's POINT
      glast_17 and POINT glast_19 to six decimals; the rectangle
      eq:condfloor leaves reproduces rem:conditionedlimit's POINT
      rectprimes and POINT rectrungs exactly; and gcount * hlog
      reproduces g to a relative 1e-12, because the split is an
      identity and not a model.
  NN2 **hlog CARRIES THE DESCENT.**  Over the rectangle's primes at
      the top rung, the span of gcount is less than half the span of
      g.
  NN3 **AND hlog IS THE DERIVED ONE.**  For every rectangle prime at
      the top rung, |hlog(p) - (log(N/p)-1)/(log N -1)| <= 0.02.
  NN4 THE LEVEL OF gcount IS THE d=1 EXCLUSION.  The median of gcount
      over the rectangle's primes is within 0.05 of A_raw/A at the
      top rung.
  NN5 AND gcount IS FLAT.  Its span over the rectangle's primes is
      under 0.10.

REFUTATION RULE (fixed before the run)

  NN1 REFUTED on any of the three legs; nothing below is reported.
      The 1e-12 leg cannot fail unless the code is wrong, and that
      is what it is for.
  NN2 REFUTED at half or more.  **Then the p-descent is not the log
      weight** and rem:derivedlimit's count is missing something
      else; the run says so and names no replacement, because a
      factor invented after the fact is what this branch keeps
      refusing.
  NN3 REFUTED outside 0.02 on any prime.  Then L(p) is not the mean
      of log m over its range -- most likely because Lambda(N-pm)
      anti-correlates with log m, small N-pm carrying small log q --
      and the derived form is wrong even if NN2 holds.  **NN2 and
      NN3 can split**, and that outcome is the informative one: the
      descent would be the log weight, measured, without the closed
      form being right.
  NN4 REFUTED outside 0.05.  Then W does not carry only the
      densities and the offset is something other than the d=1
      exclusion.
  NN5 REFUTED at 0.10 or more.  Then gcount keeps p-structure of its
      own and NN2 holding would only mean hlog carries the larger
      half.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  This run tests a
  decomposition at one rung, so it cannot see N-dependence at all:
  hlog -> 1 and A/A_raw -> 1 are both statements about the limit and
  **both are invisible here**.  rem:conditionedlimit's centre moved
  0.006052 in three doublings against 0.039116 remaining, so the
  motion this split would explain is a hundredth while the distance
  is four hundredths.  A hold means the shape at the top rung is
  accounted for, not that the limit is reached; that is a different
  run and eq:condfloor says it needs an N this machine does not have.

  WHAT THIS CANNOT DO.  The split is exact and its two halves are
  derivable, but c(p, N) is a correction to rem:residuemodel's
  elementary Mobius sum, which was refuted as a whole.  Explaining
  the correction's shape does not bound |sum a| and does not move
  item 5's demand of +0.134019.
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
OUT = os.path.join(ROOT, "results", "audit_logweight_split.txt")
SRCD = os.path.join(ROOT, "results", "audit_derived_limit.txt")
SRCL = os.path.join(ROOT, "results", "audit_limit_direction.txt")
SRCC = os.path.join(ROOT, "results", "audit_conditioned_limit.txt")

THETA = 0.56
NS = [6_400_000, 12_800_000, 25_600_000, 51_200_000,
      102_400_000, 204_800_000]
CAND = (3, 7, 11, 13, 17, 19, 23, 29, 31, 37)
NGATE = 25_600_000
RESOL = 0.01
DEC = 6
RELID = 1e-12
HALF = 0.5
TOLD = 0.02
TOLLVL = 0.05
FLATCAP = 0.10


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
    del pr, lgp
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
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


def split(N, lam, sqf):
    """W, L and g for CAND, plus the term counts and A_raw/A."""
    PN = factor_set(N)
    K = int(N ** THETA)
    D = (N - 1) // K

    def survivors(lo, hi, drop, nonprime=False):
        keep = sqf[lo:hi].copy()
        if nonprime:
            keep &= lam[lo:hi] == 0.0
        for q in drop:
            keep[((-lo) % int(q))::int(q)] = False
        return np.flatnonzero(keep).astype(np.int64) + lo

    m1 = survivors(K, N, PN, nonprime=True)
    lm1 = np.log(m1.astype(np.float64))
    hit1 = lam[N - m1]
    w1 = float(hit1.sum())
    i1 = float((hit1 * lm1).sum())
    del hit1, lm1
    mr = survivors(K, N, PN)
    araw = float((lam[N - mr] * np.log(mr.astype(np.float64))).sum())
    del mr, m1

    out = {}
    for p in CAND:
        if p > D or p in PN:
            out[p] = None
            continue
        mp = survivors(K, (N - 1) // p + 1, {p} | PN)
        lmp = np.log(mp.astype(np.float64))
        hit = lam[N - p * mp]
        wp = float(hit.sum())
        ip = float((hit * lmp).sum())
        del hit, lmp
        wm = (1.0 / phi(p)) * (1.0 - p / D) / (1.0 - 1.0 / D)
        out[p] = dict(n=int(mp.size), W=wp, L=ip / wp,
                      g=(ip / i1) / wm * (p + 1.0) / p,
                      gcount=(wp / w1) / wm * (p + 1.0) / p,
                      hlog=(ip / wp) / (i1 / w1))
        del mp
    return out, i1 / araw, K, D


def read_pub():
    out = {}
    for path, keys in ((SRCD, ("gtop_3", "gtop_11")),
                       (SRCL, ("glast_17", "glast_19")),
                       (SRCC, ("rectprimes", "rectrungs"))):
        src = io.open(path, encoding="utf-8").read()
        for k in keys:
            m = re.search(r"^POINT %s (\S+)\s*$" % k, src, re.M)
            if not m:
                raise SystemExit("missing published value: " + k)
            out[k] = m.group(1)
    return out


HEAD = [
    "STATISTIC: the exact factorisation of I(p) into a prime-hit mass",
    "           W(p) and a weighted mean of log m, and which of the",
    "           two carries the descent of g in p at the top rung.",
    "FIELD: N = %s;" % NS,
    "       p over CAND = %s," % (CAND,),
    "       K = floor(N^%.2f), D = floor((N-1)/K); the rectangle is"
    % THETA,
    "       the one eq:condfloor leaves at resolution %.2f. Six" % RESOL,
    "       values are READ from results/audit_derived_limit.txt,",
    "       results/audit_limit_direction.txt and",
    "       results/audit_conditioned_limit.txt.",
    "DERIVED: I(p) = W(p) L(p) is an identity. W is what the density",
    "         count of rem:derivedlimit predicted; L is about",
    "         log(N/p) - 1, never counted, falling in p and rising to",
    "         1 in N. Neither half is fitted.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    for f, k in ((SRCD, "gtop_3"), (SRCD, "gtop_11"),
                 (SRCL, "glast_17"), (SRCL, "glast_19"),
                 (SRCC, "rectprimes"), (SRCC, "rectrungs")):
        say("READ %s POINT %s %s"
            % (os.path.basename(f), k, pub[k]))
    say("PRINTBOUND audit_logweight_split %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, resolution %.2f, identity to %.0e,"
        % (THETA, RESOL, RELID))
    say("  span share under %.1f, derived form within %.2f,"
        % (HALF, TOLD))
    say("  level within %.2f, flatness under %.2f"
        % (TOLLVL, FLATCAP))
    say("  candidates %s" % (CAND,))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    del mu

    S, RAT = {}, {}
    for N in NS:
        s, rat, K, D = split(N, lam, sqf)
        S[N], RAT[N] = s, rat
        say("  N = %-10d K = %-6d D = %-5d A/Araw %.6f"
            % (N, K, D, rat))
    say("SCALES %d" % len(NS))

    floor = {N: math.log(N) / RESOL ** 2 for N in NS}
    best = None
    for i in range(len(NS)):
        rungs = NS[i:]
        for j in range(1, len(CAND) + 1):
            ps = CAND[:j]
            if all(S[N][p] and S[N][p]["n"] >= floor[N]
                   for N in rungs for p in ps):
                c = (len(ps) * len(rungs), len(rungs), ps, rungs)
                if best is None or c[:2] > best[:2]:
                    best = c
    _, _, PS, RUNGS = best
    NT = RUNGS[-1]

    # ------------------------------------------------------------- NN1
    say()
    say("NN1  the gate")
    nn1 = True
    for lab, p, key in (("g(3)", 3, "gtop_3"), ("g(11)", 11, "gtop_11"),
                        ("g(17)", 17, "glast_17"),
                        ("g(19)", 19, "glast_19")):
        got = S[NGATE][p]["g"]
        ok = abs(round(got, DEC) - float(pub[key])) < 10.0 ** (-DEC)
        nn1 &= ok
        say("  %-6s %.6f against its %s   %s"
            % (lab, got, pub[key], "ok" if ok else "MISMATCH"))
    okr = (len(PS) == int(pub["rectprimes"])
           and len(RUNGS) == int(pub["rectrungs"]))
    nn1 &= okr
    say("  rectangle %d by %d against its %s by %s   %s"
        % (len(PS), len(RUNGS), pub["rectprimes"], pub["rectrungs"],
           "ok" if okr else "MISMATCH"))
    worst = max(abs(S[NT][p]["gcount"] * S[NT][p]["hlog"]
                    - S[NT][p]["g"]) / abs(S[NT][p]["g"]) for p in PS)
    oki = worst <= RELID
    nn1 &= oki
    say("  gcount * hlog against g, worst relative %.2e   %s"
        % (worst, "ok" if oki else "MISMATCH"))
    say("  NN1 %s" % ("hold" if nn1 else "REFUTED"))
    if not nn1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    lg = math.log(NT)
    say()
    say("      p      g         gcount    hlog      derived   diff")
    for p in PS:
        d = (lg - math.log(p) - 1.0) / (lg - 1.0)
        r = S[NT][p]
        say("  %5d  %.6f  %.6f  %.6f  %.6f  %+.6f"
            % (p, r["g"], r["gcount"], r["hlog"], d, r["hlog"] - d))
        say("POINT gcount_%d %.6f" % (p, r["gcount"]))
        say("POINT hlog_%d %.6f" % (p, r["hlog"]))

    def span(key):
        vs = [S[NT][p][key] for p in PS]
        return max(vs) - min(vs)

    # ------------------------------------------------------------- NN2
    say()
    say("NN2  which half carries the descent?")
    sg, sc = span("g"), span("gcount")
    share = sc / sg
    nn2 = share < HALF
    say("  span of g      %.6f" % sg)
    say("  span of gcount %.6f" % sc)
    say("  span of hlog   %.6f" % span("hlog"))
    say("  share %.6f" % share)
    say("POINT spang %.6f" % sg)
    say("POINT spangcount %.6f" % sc)
    say("POINT share %.6f" % share)
    say("  NN2 %s   (cap: %.1f)"
        % ("hold" if nn2 else "REFUTED", HALF))

    # ------------------------------------------------------------- NN3
    say()
    say("NN3  is hlog the derived one?")
    devs = [abs(S[NT][p]["hlog"]
                - (lg - math.log(p) - 1.0) / (lg - 1.0)) for p in PS]
    nn3 = max(devs) <= TOLD
    say("  worst |hlog - derived| %.6f over %d primes"
        % (max(devs), len(PS)))
    say("POINT hlogworst %.6f" % max(devs))
    say("  NN3 %s   (cap: %.2f)"
        % ("hold" if nn3 else "REFUTED", TOLD))

    # ------------------------------------------------------------- NN4
    say()
    say("NN4  is gcount's level the d=1 exclusion?")
    medc = float(np.median([S[NT][p]["gcount"] for p in PS]))
    want = 1.0 / RAT[NT]
    nn4 = abs(medc - want) <= TOLLVL
    say("  median gcount %.6f, Araw/A %.6f, difference %+.6f"
        % (medc, want, medc - want))
    say("POINT medgcount %.6f" % medc)
    say("POINT arawovera %.6f" % want)
    say("  NN4 %s   (cap: %.2f)"
        % ("hold" if nn4 else "REFUTED", TOLLVL))

    # ------------------------------------------------------------- NN5
    say()
    say("NN5  is gcount flat?")
    nn5 = sc < FLATCAP
    say("  span of gcount %.6f" % sc)
    say("  NN5 %s   (cap: %.2f)"
        % ("hold" if nn5 else "REFUTED", FLATCAP))

    say()
    say("=" * 70)
    say("NN1 %s  NN2 %s  NN3 %s  NN4 %s  NN5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (nn1, nn2, nn3, nn4, nn5)))
    say()
    if nn2 and nn3:
        say("the descent of g in p is the log weight, and the log "
            "weight is the one")
        say("the construction forces. nothing is fitted on either "
            "side of an exact")
        say("identity. this accounts for the shape at one rung and "
            "says nothing")
        say("about the limit, which is a different run and a larger "
            "N.")
    elif nn2 and not nn3:
        say("the descent is the log weight and the closed form for it "
            "is wrong.")
        say("that split is the informative outcome the rule named, "
            "and no factor")
        say("is invented after the fact to close it.")
    else:
        say("the p-descent is not the log weight. the density count "
            "is missing")
        say("something else and this run does not name it.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
