# -*- coding: utf-8 -*-
r"""
Does the closed form transport, and is its level the Goldbach series?

WHAT IS AT STAKE

rem:logweightsplit factored I(p) = W(p) L(p) exactly and found both
halves derivable at the top rung: hlog sits within 0.001766 of
(log(N/p)-1)/(log N -1) and gcount's median sits within 0.001152 of
A_raw/A.  That gives

    c(p, N) ~ (p/(p+1)) (A_raw/A) (log(N/p)-1)/(log N -1)

with nothing fitted.  Two things are untested in it.

**One rung is not a law.**  Everything above was read at
N = 204800000.  A form that holds at one N and drifts at the others
is a coincidence, and the rungs below are on disk already.

**And the level has a law of its own that can be derived rather than
measured.**  A_raw - A is exactly the part of the d = 1 sum carried by
prime m, and for prime m one has log m = Lambda(m), so

    A_raw - A  =  sum_{K <= m < N, (m,N)=1, m prime} Lambda(N-m) log m

which is the Goldbach sum, and its size is SS(N) N with the singular
series SS(N) = 2 C_2 prod_{q | N, q > 2} (q-1)/(q-2).  Since
A_raw is about its own density times N (log N - 1), the prediction is
that **(1 - A/A_raw)(log N - 1) is constant across rungs** -- and
that constant, not a fitted one, is what makes A/A_raw climb to 1 and
hence what makes eq:derivedlimit's limit slow.

rem:logweightsplit also localised the one leftover: gcount is flat at
about 1.1417 from p = 7 up and p = 3 alone sits at 1.217140.  Whether
that excess belongs to p or to N is a prediction here, not an
exclusion: p = 3 is tested by its own rule and is not dropped from the
field.

BACKS: Remark {#rem:closedform} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  OO1 THE GATE.  POINT medgcount, POINT arawovera, POINT hlogworst
      and POINT share reproduce rem:logweightsplit to six decimals.
  OO2 **THE FORM TRANSPORTS.**  Over every cell (p, N) that clears
      eq:condfloor with p >= 7, |g - (A_raw/A)(log(N/p)-1)/(log N-1)|
      <= 0.01.
  OO3 **p = 3's EXCESS BELONGS TO p, NOT N.**  At every rung its
      deviation from the same form is at least 0.05, and the spread
      of those deviations across the rungs is at most 0.02.
  OO4 **THE LEVEL OBEYS A LAW.**  (1 - A/A_raw)(log N - 1) has a
      spread of at most 0.02 across the six rungs.
  OO5 **AND THE LAW IS THE GOLDBACH SERIES.**  At every rung,
      |(A_raw - A)/(SS(N) N) - 1| <= 0.05.

REFUTATION RULE (fixed before the run)

  OO1 REFUTED outside six decimals on any of the four; nothing below
      is reported.
  OO2 REFUTED above 0.01 on any clearing cell.  Then eq:cclosed is a
      statement about one rung and must be written as one; the
      remark then withdraws the word "form" and keeps only the top
      rung's numbers.
  OO3 REFUTED if any rung's deviation is under 0.05, or if the
      spread exceeds 0.02.  A deviation that shrinks with N would
      make p = 3 a finite-N effect and not a property of 3, which is
      the opposite reading of [rem:driftsigns]'s 10/10 and would
      have to be said plainly.
  OO4 REFUTED above 0.02.  Then A/A_raw does not climb at the rate
      the Goldbach sum sets and eq:cclosed's second factor has no
      derived law, only a measured drift.
  OO5 REFUTED outside 0.05 at any rung.  Then the identification of
      A_raw - A with the Goldbach sum is wrong somewhere -- most
      likely in the density that divides it -- and no factor is
      fitted afterwards to close the gap.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  OO4 and OO5 are the
  same quantity judged two ways and are not independent evidence:
  OO4 asks whether it is flat and OO5 whether it is right, and a
  single systematic error in the density would move both together.
  **Six rungs spanning five doublings cannot separate a constant from
  a term of size 1/log N**, which changes by only about a fifth
  across that range, so a flat OO4 is consistent with a form this run
  cannot see.  If OO4's spread lands between 0.015 and 0.02 the cap
  decided it and that must be reported.

  WHAT THIS CANNOT DO.  Every rule here is about c(p, N), a
  correction to rem:residuemodel's elementary Mobius sum, which was
  refuted as a whole.  A derived closed form for a refuted model's
  error is a description, not a route; it does not bound |sum a| and
  does not move item 5's demand of +0.134019.
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
OUT = os.path.join(ROOT, "results", "audit_closed_form.txt")
SRCS = os.path.join(ROOT, "results", "audit_logweight_split.txt")

THETA = 0.56
NS = [6_400_000, 12_800_000, 25_600_000, 51_200_000,
      102_400_000, 204_800_000]
CAND = (3, 7, 11, 13, 17, 19, 23, 29, 31, 37)
RESOL = 0.01
DEC = 6
TWINC2 = 0.66016181584686957392
TOLFORM = 0.01
MINEX = 0.05
EXSPREAD = 0.02
LVLSPREAD = 0.02
TOLSS = 0.05


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


def singular(N):
    """SS(N) = 2 C_2 prod_{q | N, q > 2} (q-1)/(q-2), with C_2 the
    published twin-prime constant -- no Euler product is built here."""
    s = 2.0 * TWINC2
    for q in factor_set(N):
        if q > 2:
            s *= (q - 1.0) / (q - 2.0)
    return s


def rung(N, lam, sqf):
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
    hit1 = lam[N - m1]
    w1 = float(hit1.sum())
    i1 = float((hit1 * np.log(m1.astype(np.float64))).sum())
    del hit1, m1
    mr = survivors(K, N, PN)
    araw = float((lam[N - mr] * np.log(mr.astype(np.float64))).sum())
    del mr

    out = {}
    for p in CAND:
        if p > D or p in PN:
            out[p] = None
            continue
        mp = survivors(K, (N - 1) // p + 1, {p} | PN)
        hit = lam[N - p * mp]
        wp = float(hit.sum())
        ip = float((hit * np.log(mp.astype(np.float64))).sum())
        del hit
        wm = (1.0 / phi(p)) * (1.0 - p / D) / (1.0 - 1.0 / D)
        out[p] = dict(n=int(mp.size),
                      g=(ip / i1) / wm * (p + 1.0) / p,
                      gcount=(wp / w1) / wm * (p + 1.0) / p,
                      hlog=(ip / wp) / (i1 / w1))
        del mp
    return out, i1, araw, K, D


def read_pub():
    src = io.open(SRCS, encoding="utf-8").read()
    out = {}
    for k in ("medgcount", "arawovera", "hlogworst", "share"):
        m = re.search(r"^POINT %s (\S+)\s*$" % k, src, re.M)
        if not m:
            raise SystemExit("missing published value: " + k)
        out[k] = m.group(1)
    return out


HEAD = [
    "STATISTIC: the departure of g from the closed form",
    "           (A_raw/A)(log(N/p)-1)/(log N -1) over every cell that",
    "           clears eq:condfloor, and the constancy and value of",
    "           (1 - A/A_raw)(log N - 1) across six rungs.",
    "FIELD: N = %s;" % NS,
    "       p over CAND = %s," % (CAND,),
    "       K = floor(N^%.2f), D = floor((N-1)/K); a cell is kept"
    % THETA,
    "       when n >= log N / %.2f^2. Four values are READ from"
    % RESOL,
    "       results/audit_logweight_split.txt.",
    "DERIVED: A_raw - A is the part of the d=1 sum carried by prime m,",
    "         where log m = Lambda(m), so it is the Goldbach sum and",
    "         is SS(N) N with SS(N) = 2 C_2 prod_{q|N, q>2} (q-1)/(q-2);",
    "         C_2 = %.17f is the published twin-prime constant and no"
    % TWINC2,
    "         Euler product is built at the measurement range.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    for k in ("medgcount", "arawovera", "hlogworst", "share"):
        say("READ audit_logweight_split.txt POINT %s %s" % (k, pub[k]))
    say("PRINTBOUND audit_closed_form %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, resolution %.2f, form within %.2f,"
        % (THETA, RESOL, TOLFORM))
    say("  p = 3 excess at least %.2f with spread at most %.2f,"
        % (MINEX, EXSPREAD))
    say("  level spread at most %.2f, series within %.2f"
        % (LVLSPREAD, TOLSS))
    say("  candidates %s" % (CAND,))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    del mu

    R, A1, ARAW = {}, {}, {}
    for N in NS:
        r, i1, araw, K, D = rung(N, lam, sqf)
        R[N], A1[N], ARAW[N] = r, i1, araw
        say("  N = %-10d K = %-6d D = %-5d A/Araw %.6f  SS(N) %.6f"
            % (N, K, D, i1 / araw, singular(N)))
    say("SCALES %d" % len(NS))

    def form(p, N):
        lg = math.log(N)
        return ((ARAW[N] / A1[N])
                * (lg - math.log(p) - 1.0) / (lg - 1.0))

    floor = {N: math.log(N) / RESOL ** 2 for N in NS}

    # ------------------------------------------------------------- OO1
    say()
    say("OO1  the gate")
    NT = NS[-1]
    lgT = math.log(NT)
    gc = [R[NT][p]["gcount"] for p in CAND]
    oo1 = True
    for k, v in (("medgcount", float(np.median(gc))),
                 ("arawovera", ARAW[NT] / A1[NT])):
        ok = abs(round(v, DEC) - float(pub[k])) < 10.0 ** (-DEC)
        oo1 &= ok
        say("  %-10s %.6f against its %s   %s"
            % (k, v, pub[k], "ok" if ok else "MISMATCH"))
    sg = (max(R[NT][p]["g"] for p in CAND)
          - min(R[NT][p]["g"] for p in CAND))
    share = (max(gc) - min(gc)) / sg
    oksh = abs(round(share, DEC) - float(pub["share"])) < 10.0 ** (-DEC)
    oo1 &= oksh
    say("  share      %.6f against its %s   %s"
        % (share, pub["share"], "ok" if oksh else "MISMATCH"))
    hw = max(abs(R[NT][p]["hlog"]
                 - (lgT - math.log(p) - 1.0) / (lgT - 1.0))
             for p in CAND)
    okh = abs(round(hw, DEC) - float(pub["hlogworst"])) < 10.0 ** (-DEC)
    oo1 &= okh
    say("  hlogworst  %.6f against its %s   %s"
        % (hw, pub["hlogworst"], "ok" if okh else "MISMATCH"))
    say("  OO1 %s" % ("hold" if oo1 else "REFUTED"))
    if not oo1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------------------- OO2
    say()
    say("OO2  does the form transport across the rungs?")
    say("      p " + "".join("%12d" % N for N in NS))
    worst, ncell = 0.0, 0
    for p in CAND:
        row = ""
        for N in NS:
            r = R[N][p]
            if not r or r["n"] < floor[N]:
                row += "           ."
                continue
            d = r["g"] - form(p, N)
            row += "  %+.6f" % d
            if p >= 7:
                worst = max(worst, abs(d))
                ncell += 1
        say("  %5d %s" % (p, row))
    oo2 = worst <= TOLFORM
    say("  %d clearing cells with p >= 7, worst |departure| %.6f"
        % (ncell, worst))
    say("POINT cells %d" % ncell)
    say("POINT formworst %.6f" % worst)
    say("  OO2 %s   (cap: %.2f)"
        % ("hold" if oo2 else "REFUTED", TOLFORM))

    # ------------------------------------------------------------- OO3
    say()
    say("OO3  does p = 3's excess belong to p or to N?")
    e3 = [R[N][3]["g"] - form(3, N) for N in NS]
    sp3 = max(e3) - min(e3)
    oo3 = min(abs(v) for v in e3) >= MINEX and sp3 <= EXSPREAD
    say("  departures %s" % " ".join("%+.6f" % v for v in e3))
    say("  smallest %.6f, spread %.6f"
        % (min(abs(v) for v in e3), sp3))
    say("POINT ex3min %.6f" % min(abs(v) for v in e3))
    say("POINT ex3spread %.6f" % sp3)
    say("  OO3 %s   (floor %.2f, cap %.2f)"
        % ("hold" if oo3 else "REFUTED", MINEX, EXSPREAD))

    # ------------------------------------------------------------- OO4
    say()
    say("OO4  is the level constant against log N - 1?")
    lev = [(1.0 - A1[N] / ARAW[N]) * (math.log(N) - 1.0) for N in NS]
    spl = max(lev) - min(lev)
    oo4 = spl <= LVLSPREAD
    say("  (1 - A/Araw)(log N - 1)  %s"
        % " ".join("%.6f" % v for v in lev))
    say("  spread %.6f" % spl)
    say("POINT levfirst %.6f" % lev[0])
    say("POINT levlast %.6f" % lev[-1])
    say("POINT levspread %.6f" % spl)
    say("  OO4 %s   (cap: %.2f)"
        % ("hold" if oo4 else "REFUTED", LVLSPREAD))

    # ------------------------------------------------------------- OO5
    say()
    say("OO5  is that level the Goldbach series?")
    rs = [(ARAW[N] - A1[N]) / (singular(N) * N) for N in NS]
    oo5 = max(abs(v - 1.0) for v in rs) <= TOLSS
    say("  (Araw - A) / (SS(N) N)  %s"
        % " ".join("%.6f" % v for v in rs))
    say("  worst |ratio - 1| %.6f" % max(abs(v - 1.0) for v in rs))
    say("POINT ssfirst %.6f" % rs[0])
    say("POINT sslast %.6f" % rs[-1])
    say("  OO5 %s   (cap: %.2f)"
        % ("hold" if oo5 else "REFUTED", TOLSS))

    say()
    say("=" * 70)
    say("OO1 %s  OO2 %s  OO3 %s  OO4 %s  OO5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (oo1, oo2, oo3, oo4, oo5)))
    say()
    if oo2 and oo4:
        say("the closed form is a form and not one rung, and its "
            "level moves at a")
        say("rate the goldbach series sets. that is why "
            "eq:derivedlimit's limit is")
        say("slow: the approach is a reciprocal logarithm and not a "
            "power. nothing")
        say("here bounds the sum this repository is trying to bound.")
    else:
        say("the closed form does not survive its own rungs, and what "
            "was read at")
        say("the top is a description of one N.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
