# -*- coding: utf-8 -*-
r"""
The seed rem:identitydeep said was owed

WHAT IS AT STAKE

rem:identitydeep found the one thing in this repository that a
resolved random multiplicative ensemble fails to reach: T/(S(N)N) is
1.0039, 0.9865, 0.9893 for mu and no draw of 512 comes below 1.9149,
2.0030, 2.1084.  It wrote its own caveat -- one seed per ensemble, and
another owed before the result is leaned on -- and this run pays it.

Two things are being checked and they are different.  Whether the
refutation reproduces at all is one: a second and third seed either
also fail to reach mu or they do not.  Whether the *bound* is stable
is the other, and it is the weaker of the two: 1.9149 is a minimum
over 512 draws, a single order statistic, and a second seed could
easily put it somewhat lower without any draw coming near mu.

Only the multiplicative arm is run.  The iid arm's own minimum was
4.8387, four and a half times mu, and nothing about its
reproducibility is in question; running it again would cost as much as
this whole script and answer nothing that was asked.

BACKS: Remark {#rem:identityseeds} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  M1  THE GATE.  mu reproduces rem:identitydeep's three ratios to
      four decimals.  mu's arm is deterministic, so a mismatch means
      the construction here is not that one.
  M2  **The refutation reproduces.**  At each of the two new seeds and
      each of the three N, none of 512 multiplicative draws reaches
      mu.
  M3  The bound is stable in size if not in value: each new seed's
      smallest ratio is within the ensemble's own standard deviation
      of the first seed's, which rem:identitydeep put at 2.0400 to
      2.3376.
  M4  And the bulk agrees more tightly than the tail: each new seed's
      median is within a tenth of the first seed's, which was 4.2771,
      4.5269, 4.7875.

REFUTATION RULE (fixed before the run)

  M1  REFUTED outside four decimals; nothing below is reported.
  M2  **REFUTED if any draw at any seed reaches mu.**  Then the first
      seed was lucky and the identity separation is not what
      rem:identitydeep recorded -- that remark would need its claim
      cut to the seeds that show it, and the count of reaching draws
      is the measure.
  M3  REFUTED if a new minimum sits more than one ensemble standard
      deviation from the first.  **The unresolved case is named**:
      these are minima over 512 draws, and the spread of a minimum is
      not the spread of the ensemble, so this rule is generous by
      construction and failing it would mean the tail is much less
      stable than the bulk -- not that the separation is in doubt.
      M3 gates nothing.
  M4  REFUTED outside a tenth.  A median that moves by more than that
      across seeds would mean 512 draws do not pin the bulk either,
      which would put M3 and the whole design in question rather than
      just the tail.

  WHAT THIS CANNOT DO.  Three of the five N that block publishes, and
  three seeds.  A separation that reproduces here is not shown to
  reproduce at the two larger N, which remain untested.
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
OUT = os.path.join(ROOT, "results",
                   "audit_directidentity_seeds.txt")
SRC = os.path.join(ROOT, "results", "audit_directidentity_deep.txt")

NS = [200_000, 400_000, 800_000]
NMAX = max(NS)
CLIM = 4_000_000
PASSES = 8
PERPASS = 64
DRAWS = PASSES * PERPASS
SEEDS = [20260907, 20260908]
DEC = 4
MEDTOL = 0.1


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n):
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
    return pr, lam, mu


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


def read_prev():
    """mu's ratios and the first seed's min, sd and median"""
    src = io.open(SRC, encoding="utf-8").read()
    mus, mins, meds, sds = [], [], [], []
    for n in NS:
        m = re.search(r"^\s*multiplicative\s+%d\s+([\d.]+)\s+"
                      r"([\d.]+)\s+([\d.]+)\s" % n, src, re.M)
        if not m:
            raise SystemExit("no multiplicative row for N = %d" % n)
        mus.append(float(m.group(1)))
        mins.append(float(m.group(2)))
        meds.append(float(m.group(3)))
        s = re.search(r"^  their spreads are ([\d.]+) and", src, re.M)
        sds.append(s)
    sds = [float(x) for x in re.findall(
        r"^  their spreads are ([\d.]+) and", src, re.M)]
    if len(sds) != len(NS):
        raise SystemExit("expected %d spread lines, found %d"
                         % (len(NS), len(sds)))
    return mus, mins, meds, sds


HEAD = [
    "STATISTIC: |T|/(S(N)N) for mu and for %d random multiplicative"
    % DRAWS,
    "           sign patterns at each of two fresh seeds, counting",
    "           draws that reach mu and comparing each seed's",
    "           minimum and median with the first seed's.",
    "FIELD: N = %s; k over every squarefree 2 <= k < N; S(N) from"
    % (NS,),
    "       the twin product at the fixed bound %d; field, k-range"
    % CLIM,
    "       and S(N) taken unchanged from",
    "       code/audit_directidentity_null.py. The first seed's",
    "       numbers are READ from",
    "       results/audit_directidentity_deep.txt. Only the",
    "       multiplicative arm is run; the iid arm's minimum was four",
    "       and a half times mu and is not in question.",
    "SEEDS: numpy default_rng at %s; without them the file does not"
    % (SEEDS,),
    "       reproduce its own control.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    mus, pmin, pmed, psd = read_prev()
    for n, a, b, c in zip(NS, mus, pmin, pmed):
        say("READ audit_directidentity_deep.txt %d %.4f" % (n, a))
    say("  the first seed gave minima %s and medians %s"
        % (", ".join("%.4f" % x for x in pmin),
           ", ".join("%.4f" % x for x in pmed)))
    say("  with ensemble spreads %s"
        % ", ".join("%.4f" % x for x in psd))
    say("PRINTBOUND audit_directidentity_seeds %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d draws at each of %d fresh seeds" % (DRAWS, len(SEEDS)))

    say("sieving to %d" % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    muf = mu.astype(np.float64)
    del mu
    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    ks, lgs, idxs, Ss = {}, {}, {}, {}
    for N in NS:
        kk = np.flatnonzero(sqf[2:N]).astype(np.int64) + 2
        ks[N] = kk
        lgs[N] = np.log(kk.astype(np.float64))
        idxs[N] = np.arange(1, N, dtype=np.int64)
        S = twin
        for q in sorted(factor_set(N)):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        Ss[N] = S

    def ratio_of(sg):
        out = []
        for N in NS:
            f0 = np.zeros(N, dtype=np.float64)
            f0[1:] = lam[1:N] * sg[N - idxs[N]]
            kk = ks[N]
            lg = lgs[N]
            s = 0.0
            for i, k in enumerate(kk):
                k = int(k)
                r = N % k
                a = f0[r::k].sum() if r else f0[k::k].sum()
                s += lg[i] * sg[k] * a
            del f0
            out.append(abs(s) / (Ss[N] * N))
        return out

    real = ratio_of(muf)

    # -------------------------------------------------------------- M1
    say()
    say("M1  does this construction reproduce mu's ratios?")
    ok = True
    for N, r, v in zip(NS, real, mus):
        g = abs(round(r, DEC) - round(v, DEC)) < 10.0 ** (-DEC) / 2
        ok &= g
        say("  N = %-9d here %.4f against its %.4f  %s"
            % (N, r, v, "ok" if g else "MISMATCH"))
    say("  M1 %s   (cap: %d decimals at all three)"
        % ("hold" if ok else "REFUTED", DEC))
    if not ok:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    pools = {}
    for sd_ in SEEDS:
        rng = np.random.default_rng(sd_)
        pool = np.zeros((DRAWS, len(NS)))
        row = 0
        for k in range(PASSES):
            t = np.zeros(NMAX + 1, dtype=np.uint64)
            hi = rng.integers(0, 1 << 32, size=len(pr),
                              dtype=np.uint64)
            lo = rng.integers(0, 1 << 32, size=len(pr),
                              dtype=np.uint64)
            pm = (hi << np.uint64(32)) | lo
            for i, p in enumerate(pr):
                t[int(p)::int(p)] ^= pm[i]
            del hi, lo, pm
            for d in range(PERPASS):
                bit = (t >> np.uint64(d)) & np.uint64(1)
                sg = np.where(sqf,
                              1.0 - 2.0 * bit.astype(np.float64), 0.0)
                pool[row] = ratio_of(sg)
                row += 1
                del sg, bit
            del t
            say("  seed %d pass %d of %d, %d draws"
                % (sd_, k + 1, PASSES, row))
        pools[sd_] = pool

    say()
    say("      seed          N        mu    draw min   draw median"
        "    draws reaching mu")
    reach, mins, meds = {}, {}, {}
    for sd_ in SEEDS:
        pool = pools[sd_]
        rr, mn, md = [], [], []
        for i, N in enumerate(NS):
            c = int((pool[:, i] <= real[i]).sum())
            rr.append(c)
            mn.append(float(pool[:, i].min()))
            md.append(float(np.median(pool[:, i])))
            say("      %-12d %8d  %.4f  %9.4f  %11.4f  %14d"
                % (sd_, N, real[i], mn[-1], md[-1], c))
            say("POINT identseed_%d_%d %.5f" % (sd_, N, mn[-1]))
        reach[sd_], mins[sd_], meds[sd_] = rr, mn, md
    say("SCALES %d" % len(NS))

    # -------------------------------------------------------------- M2
    say()
    say("M2  does the refutation reproduce at both seeds?")
    m2 = all(c == 0 for sd_ in SEEDS for c in reach[sd_])
    for sd_ in SEEDS:
        say("  seed %d: %s of %d at each N"
            % (sd_, ", ".join(str(c) for c in reach[sd_]), DRAWS))
    say("COUNT identseed_reach %d"
        % sum(c for sd_ in SEEDS for c in reach[sd_]))
    say("  M2 %s   (cap: none at any seed or N)"
        % ("hold" if m2 else "REFUTED"))

    # -------------------------------------------------------------- M3
    say()
    say("M3  is the bound stable?")
    m3 = True
    say("      seed          N    first min   this min      gap   "
        "ensemble sd")
    for sd_ in SEEDS:
        for i, N in enumerate(NS):
            g = abs(mins[sd_][i] - pmin[i])
            m3 &= g <= psd[i]
            say("      %-12d %8d  %9.4f  %9.4f  %7.4f  %11.4f"
                % (sd_, N, pmin[i], mins[sd_][i], g, psd[i]))
    say("  M3 %s   (cap: within the ensemble's own spread)"
        % ("hold" if m3 else "REFUTED"))

    # -------------------------------------------------------------- M4
    say()
    say("M4  does the bulk agree more tightly?")
    m4 = True
    worst = 0.0
    for sd_ in SEEDS:
        for i, N in enumerate(NS):
            g = abs(meds[sd_][i] - pmed[i])
            worst = max(worst, g)
            m4 &= g <= MEDTOL
    say("  the largest median gap across seeds and N is %.4f" % worst)
    say("SPREAD identseed_median %.4f" % worst)
    say("  M4 %s   (cap: %.1f)" % ("hold" if m4 else "REFUTED",
                                   MEDTOL))

    say()
    say("=" * 70)
    say("M1 %s  M2 %s  M3 %s  M4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (ok, m2, m3, m4)))
    say()
    if m2:
        say("the seed rem:identitydeep said was owed is paid, twice. "
            "No draw of")
        say("512 at either fresh seed reaches mu at any of the three "
            "N, so that")
        say("remark's refutation is not one seed's luck. What it is "
            "still bounded")
        say("by is the smallest ratio seen, now over three seeds, and "
            "the two")
        say("larger N of that block remain untested.")
    else:
        say("a fresh seed reaches mu where the first did not, so "
            "rem:identitydeep's")
        say("refutation was that seed's. The counts above say how "
            "often, and that")
        say("remark's claim has to be cut to what reproduces.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
