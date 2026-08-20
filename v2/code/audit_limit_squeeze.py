# -*- coding: utf-8 -*-
r"""
Does the curve close on one, and is the closing its own?

WHAT IS AT STAKE

rem:limitdirection killed both ways out of rem:derivedlimit's refuted
JJ4 and left a third picture in their place: g(3) and g(7) come down
toward 1, g(11) has stopped at 0.991622, and g(17) and g(19) come up
from 0.916563 and 0.926432 toward it.  **The curve looks like it is
closing on the derived value from both sides.**  That was a reading of
five primes chosen for other reasons, and a reading is not a count.

Closing is one number per N and needs no window and no fit: over a
prime set fixed before the run, the span max(g) - min(g).  If it falls
at every doubling the closing is measured.  But a span that shrinks
says nothing about where -- a curve can close on 0.86 as easily as on
1 -- so the centre must be counted beside it, and the bulk must be
counted beside that, because at the top N only 3 and 7 sit above 1 and
eighteen of twenty sit below.

And the closing might not be its own.  rem:limitdirection found that
I(1) alone drops prime-power m, that the resulting A/A_raw is climbing
across the eleven N from 0.763851 to 0.861216, and that lifting the
exclusion shifts g by an amount that is itself monotone in p, from
-0.160203 to -0.118401.  **A shift that varies with p by 0.041802 is a
compression of exactly the kind being measured**, and it is growing.
So the same span is computed under both normalisations.

BACKS: Remark {#rem:limitsqueeze} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  LL1 THE GATE.  g(3) at the first and last N reproduce
      rem:derivedlimit's POINT g3first and POINT g3last to six
      decimals, and g(17) and g(19) at the last N reproduce
      rem:limitdirection's POINT glast_17 and POINT glast_19.
  LL2 **THE SPAN CLOSES.**  Over the twenty primes named in PSET, the
      span max(g) - min(g) falls at 9 or more of the ten doublings.
  LL3 **AND IT CLOSES ON ONE.**  The median of g over PSET is nearer
      to 1 at the last N than at the first, and at the last N is
      within 0.10 of 1.
  LL4 THE CLOSING IS NOT THE d=1 SHADOW.  Recomputed with the
      prime-power exclusion lifted, the span still falls at 9 or more
      of the ten doublings.
  LL5 **THE BULK CLIMBS.**  At least 15 of the 20 primes have g
      higher at the last N than at the first.  Eighteen of twenty sit
      below 1, so a closing on 1 requires most of them to rise.

REFUTATION RULE (fixed before the run)

  LL1 REFUTED outside six decimals on any of the four; nothing below
      is reported.
  LL2 REFUTED at 6 or fewer.  Then the curve does not close and
      rem:limitdirection's five-prime reading was a selection.  A
      span at exactly 7 or 8 of ten is a lean, not evidence: the
      two-sided coin tails there are 0.343750 and 0.109375 and the
      set was named in advance but the twenty primes are not
      independent, so no multiplicity correction is honest here and
      only a 9 or 10 is reported as a count.
  LL3 REFUTED if the median is not nearer at the last N, **or** if it
      is outside 0.10 there.  Then the span may close but not on the
      derived value, and eq:derivedlimit is a coincidence of the
      first few primes.  The top-N median is already visible in
      rem:derivedlimit's twenty-prime table as sitting near 0.90, so
      **this cap of 0.10 is a near thing by construction and a hold
      is not a confirmation** -- it is the statement that the centre
      has not yet excluded 1.
  LL4 REFUTED at 6 or fewer under the lifted exclusion.  Then the
      closing lives in the d=1 normalisation rather than in the
      weights, and the whole reading transfers to a question about
      A/A_raw.
  LL5 REFUTED below 15.  Then the closing is the head coming down to
      meet a tail that stays put, the centre cannot reach 1, and LL3
      passing would have to be an accident of the median's position.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  LL2, LL3, LL4 and
  LL5 can all hold and still not give the limit.  The span at the top
  N is about 0.30 and the centre sits about 0.10 below 1; at the rate
  rem:derivedlimit measured for g(3) -- 0.089063 over ten doublings --
  closing that 0.10 takes another ten doublings, which is N beyond
  2.6e10 and beyond this machine.  **So a hold here says the field is
  consistent with eq:derivedlimit at the resolution it has, and the
  eleven N cannot say more.**  If LL3 lands between 0.09 and 0.10 the
  cap was decided by rounding and must be reported that way.

  WHAT THIS CANNOT DO.  Nothing here bounds |sum a| or moves item 5's
  demand.  A closing span is a statement about a model correction,
  and the model it corrects -- rem:residuemodel's elementary Mobius
  sum -- was refuted as a whole; what is being learned is the shape
  of the error, not a route to the exponent.
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
OUT = os.path.join(ROOT, "results", "audit_limit_squeeze.txt")
SRCD = os.path.join(ROOT, "results", "audit_derived_limit.txt")
SRCL = os.path.join(ROOT, "results", "audit_limit_direction.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000, 6_400_000, 12_800_000, 25_600_000]
PSET = (3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
        61, 67, 71, 73, 79)
DEC = 6
FALLMIN = 9
NEAR = 0.10
RISEMIN = 15


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


def gvals(N, lam, mu, sqf):
    """g(p, N) for the fixed prime set, under both normalisations."""
    PN = factor_set(N)
    K = int(N ** THETA)
    D = (N - 1) // K
    ms = np.arange(K, N, dtype=np.int64)
    keep = sqf[ms]
    for q in PN:
        keep &= (ms % int(q)) != 0
    m1 = ms[keep]
    lg1 = np.log(m1.astype(np.float64))
    araw = float((lam[N - m1] * lg1).sum())
    good = lam[m1] == 0.0
    a = float((lam[N - m1[good]] * lg1[good]).sum())
    del ms, keep, m1, lg1, good

    g, gr = {}, {}
    for p in PSET:
        if p > D or p in PN:
            raise SystemExit("prime %d is not in the field at N = %d"
                             % (p, N))
        mp = np.arange(K, (N - 1) // p + 1, dtype=np.int64)
        kp = sqf[mp]
        for q in {p} | PN:
            kp &= (mp % int(q)) != 0
        mp = mp[kp]
        v = float((lam[N - p * mp]
                   * np.log(mp.astype(np.float64))).sum())
        w = (1.0 / phi(p)) * (1.0 - p / D) / (1.0 - 1.0 / D)
        g[p] = (v / a) / w * (p + 1.0) / p
        gr[p] = (v / araw) / w * (p + 1.0) / p
        del mp, kp
    return g, gr, a / araw, D


def read_pub():
    d = io.open(SRCD, encoding="utf-8").read()
    l = io.open(SRCL, encoding="utf-8").read()
    out = {}
    for src, key in ((d, "g3first"), (d, "g3last"),
                     (l, "glast_17"), (l, "glast_19")):
        m = re.search(r"^POINT %s (\S+)\s*$" % key, src, re.M)
        if not m:
            raise SystemExit("a published value is missing: " + key)
        out[key] = m.group(1)
    return out


def tail(k, tot):
    hi = max(k, tot - k)
    return 2.0 * sum(math.comb(tot, j)
                     for j in range(hi, tot + 1)) / 2.0 ** tot


HEAD = [
    "STATISTIC: the span max(g) - min(g) and the median of g over a",
    "           prime set fixed before the run, at eleven N, under",
    "           both normalisations of I(1); and how many of the",
    "           twenty primes are higher at the last N.",
    "FIELD: N = %s;" % NS,
    "       p over PSET = %s," % (PSET,),
    "       every one of them squarefree, at most D = floor((N-1)/K)",
    "       and coprime to every N, with K = floor(N^%.2f)." % THETA,
    "       g(3) at the ends and g(17), g(19) at the last N are READ",
    "       from results/audit_derived_limit.txt and",
    "       results/audit_limit_direction.txt.",
    "DERIVED: a span is one number per N and needs no window; a span",
    "         that closes says nothing about where, so the median and",
    "         the count of risers are computed beside it.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    say("READ audit_derived_limit.txt POINT g3first %s" % pub["g3first"])
    say("READ audit_derived_limit.txt POINT g3last %s" % pub["g3last"])
    say("READ audit_limit_direction.txt POINT glast_17 %s"
        % pub["glast_17"])
    say("READ audit_limit_direction.txt POINT glast_19 %s"
        % pub["glast_19"])
    say("PRINTBOUND audit_limit_squeeze %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, %d primes, falls at least %d of ten,"
        % (THETA, len(PSET), FALLMIN))
    say("  centre within %.2f of one, risers at least %d"
        % (NEAR, RISEMIN))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    G, GR, RAT = {}, {}, {}
    for N in NS:
        g, gr, rat, D = gvals(N, lam, mu, sqf)
        G[N], GR[N], RAT[N] = g, gr, rat
        vs = np.array([g[p] for p in PSET])
        say("  N = %-10d D = %-5d span %.6f  median %.6f  A/Araw %.6f"
            % (N, D, vs.max() - vs.min(), float(np.median(vs)), rat))
    say("SCALES %d" % len(NS))

    # ------------------------------------------------------------- LL1
    say()
    say("LL1  the gate")
    checks = [("g(3) first", G[NS[0]][3], pub["g3first"]),
              ("g(3) last", G[NS[-1]][3], pub["g3last"]),
              ("g(17) last", G[NS[-1]][17], pub["glast_17"]),
              ("g(19) last", G[NS[-1]][19], pub["glast_19"])]
    ll1 = True
    for lab, got, want in checks:
        ok = abs(round(got, DEC) - float(want)) < 10.0 ** (-DEC)
        ll1 &= ok
        say("  %-11s %.6f against its %s   %s"
            % (lab, got, want, "ok" if ok else "MISMATCH"))
    say("  LL1 %s" % ("hold" if ll1 else "REFUTED"))
    if not ll1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    def spans(H):
        return [float(np.max([H[N][p] for p in PSET])
                      - np.min([H[N][p] for p in PSET])) for N in NS]

    # ------------------------------------------------------------- LL2
    say()
    say("LL2  does the span close?")
    sp = spans(G)
    dn = sum(1 for i in range(len(sp) - 1) if sp[i + 1] < sp[i])
    tot = len(sp) - 1
    ll2 = dn >= FALLMIN
    say("  span  %s" % " ".join("%.4f" % v for v in sp))
    say("  falls %d/%d, two-sided coin tail %.6f"
        % (dn, tot, tail(dn, tot)))
    say("POINT spanfirst %.6f" % sp[0])
    say("POINT spanlast %.6f" % sp[-1])
    say("POINT spanfalls %d" % dn)
    say("  LL2 %s   (floor: %d of ten)"
        % ("hold" if ll2 else "REFUTED", FALLMIN))

    # ------------------------------------------------------------- LL3
    say()
    say("LL3  does it close on one?")
    med = [float(np.median([G[N][p] for p in PSET])) for N in NS]
    nearer = abs(med[-1] - 1.0) < abs(med[0] - 1.0)
    inside = abs(med[-1] - 1.0) <= NEAR
    ll3 = nearer and inside
    say("  median %s" % " ".join("%.4f" % v for v in med))
    say("  |median - 1| first %.6f, last %.6f"
        % (abs(med[0] - 1.0), abs(med[-1] - 1.0)))
    say("POINT medfirst %.6f" % med[0])
    say("POINT medlast %.6f" % med[-1])
    say("  nearer at the last: %s; inside the cap: %s"
        % (nearer, inside))
    say("  LL3 %s   (cap: %.2f)"
        % ("hold" if ll3 else "REFUTED", NEAR))

    # ------------------------------------------------------------- LL4
    say()
    say("LL4  is the closing the d=1 shadow?")
    spr = spans(GR)
    dnr = sum(1 for i in range(len(spr) - 1) if spr[i + 1] < spr[i])
    ll4 = dnr >= FALLMIN
    say("  span without the exclusion  %s"
        % " ".join("%.4f" % v for v in spr))
    say("  falls %d/%d" % (dnr, tot))
    say("  A/Araw %s" % " ".join("%.4f" % RAT[N] for N in NS))
    say("POINT spanrawfalls %d" % dnr)
    say("POINT spanrawlast %.6f" % spr[-1])
    say("  LL4 %s   (floor: %d of ten)"
        % ("hold" if ll4 else "REFUTED", FALLMIN))

    # ------------------------------------------------------------- LL5
    say()
    say("LL5  does the bulk climb?")
    up = [p for p in PSET if G[NS[-1]][p] > G[NS[0]][p]]
    ll5 = len(up) >= RISEMIN
    say("      p     g first   g last    change")
    for p in PSET:
        say("  %5d   %.6f  %.6f  %+.6f"
            % (p, G[NS[0]][p], G[NS[-1]][p],
               G[NS[-1]][p] - G[NS[0]][p]))
    say("  %d of %d rise: %s" % (len(up), len(PSET), up))
    say("POINT risers %d" % len(up))
    say("  LL5 %s   (floor: %d)"
        % ("hold" if ll5 else "REFUTED", RISEMIN))

    say()
    say("=" * 70)
    say("LL1 %s  LL2 %s  LL3 %s  LL4 %s  LL5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (ll1, ll2, ll3, ll4, ll5)))
    say()
    if ll2 and ll3 and ll5:
        say("the curve closes, it closes on one, and the bulk is what "
            "moves. the")
        say("derived value survives a count and not a reading -- at "
            "the resolution")
        say("eleven N carry, and the rule said in advance that this "
            "is consistency")
        say("and not the limit.")
    elif ll2 and not ll3:
        say("the curve closes and not on one. eq:derivedlimit is then "
            "a coincidence")
        say("of the first few primes and the centre says so.")
    else:
        say("the curve does not close by a count, and "
            "rem:limitdirection's reading")
        say("of five primes was a selection.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
