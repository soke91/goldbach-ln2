# -*- coding: utf-8 -*-
r"""
The last two N of rem:identitynull, and what depth buys at the low end

WHAT IS AT STAKE

rem:identitydeep tested three of the five N that rem:identitynull
publishes and found no multiplicative draw reaching mu;
rem:identityseeds paid the seed that remark owed, twice, with the same
answer.  Both said the same thing about what was left: **the two
larger N are untested.**  This is item 8's last piece.

The reduction that makes it affordable rests on a measurement rather
than on hope.  The statistic runs over every squarefree k below N, so
a draw at these two N costs about ten seconds against under two at the
three already done, and 512 draws would be an hour and a half.  This
run uses 256.  rem:identityseeds found the minimum to be the *stable*
statistic of this ensemble -- across three seeds the minima moved by
at most 0.1702 while the medians moved by 1.1706, a factor of seven --
and the minimum is exactly what the question needs.  **Halving the
depth still costs something at the low end, and rather than assume it
is small this run measures it**: the minimum over the first 64, 128
and 256 draws is printed, so the drift with depth is visible.

BACKS: Remark {#rem:identitybig} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  N1  THE GATE.  mu reproduces rem:identitynull's published ratios at
      the two N, 1.0017 and 0.9968, to four decimals.
  N2  **The separation holds where it has never been looked**: none of
      256 multiplicative draws reaches mu at either N.
  N3  And it widens with N, as the three smaller N did -- their minima
      were 1.9149, 2.0030, 2.1084 at the first seed.  The minimum at
      1.6e6 exceeds 2.1084 and the minimum at 3.2e6 exceeds the one at
      1.6e6.
  N4  Depth buys little at the low end: the minimum over 256 draws is
      within 0.2 of the minimum over the first 64.

REFUTATION RULE (fixed before the run)

  N1  REFUTED outside four decimals; nothing below is reported.
  N2  **REFUTED if any draw reaches mu at either N.**  Then the
      identity separation is a property of small N and
      rem:identitydeep's claim has to be cut to the range that shows
      it -- which would matter, because that remark is the only
      surviving escape in this paper.
  N3  REFUTED if either minimum fails to exceed its predecessor.
      **This comparison is biased in N3's favour and the bias is
      named**: 256 draws give a higher minimum than 512 do, so a rise
      here is partly the depth and not only the N. N4 measures how
      much, and if N4 fails then N3 holds for a reason it did not
      intend and must be read as unmeasured rather than as shown.
  N4  REFUTED outside 0.2.  Then the low end is not as depth-stable
      as rem:identityseeds found across seeds, the reduction to 256
      was not justified by that finding, and N3 is what suffers.

  WHAT THIS CANNOT DO.  One seed at these two N, and 256 draws.  A
  count of zero bounds a tail no better than 256 draws can, and the
  smallest ratio seen is printed as the bound this run gives.
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
OUT = os.path.join(ROOT, "results", "audit_directidentity_bigN.txt")
SRC = os.path.join(ROOT, "results", "audit_directidentity_null.txt")
SRCD = os.path.join(ROOT, "results",
                    "audit_directidentity_deep.txt")

NS = [1_600_000, 3_200_000]
NMAX = max(NS)
CLIM = 4_000_000
PASSES = 4
PERPASS = 64
DRAWS = PASSES * PERPASS
CHECKS = (64, 128, 256)
SEED = 20260909
DEC = 4
DEPTHTOL = 0.2


assert DRAWS == max(CHECKS)


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


def read_smallmins():
    """the three smaller N's multiplicative minima, from that file"""
    src = io.open(SRCD, encoding="utf-8").read()
    out = [float(x) for x in re.findall(
        r"^\s*multiplicative\s+\d+\s+[\d.]+\s+([\d.]+)\s", src,
        re.M)]
    if len(out) != 3:
        raise SystemExit("expected 3 multiplicative rows, got %d"
                         % len(out))
    return out


def read_pub():
    src = io.open(SRC, encoding="utf-8").read()
    blk = src[src.index("V2  T / (S(N)N)"):]
    out = []
    for n in NS:
        m = re.search(r"^\s*%d\s+([\d.]+)\s" % n, blk, re.M)
        if not m:
            raise SystemExit("no V2 row for N = %d" % n)
        out.append(float(m.group(1)))
    return out


HEAD = [
    "STATISTIC: |T|/(S(N)N) at the two largest N of",
    "           rem:identitynull, for mu and for %d random" % DRAWS,
    "           multiplicative sign patterns, counting draws that",
    "           reach mu and printing the minimum over the first 64,",
    "           128 and 256 draws so the drift with depth is visible.",
    "FIELD: N = %s; k over every squarefree 2 <= k < N; S(N) from"
    % (NS,),
    "       the twin product at the fixed bound %d; field, k-range"
    % CLIM,
    "       and S(N) taken unchanged from",
    "       code/audit_directidentity_null.py. mu's published ratios",
    "       are READ from results/audit_directidentity_null.txt and",
    "       re-measured here as the gate. Only the multiplicative arm",
    "       is run and only 256 draws, against 512 at the three",
    "       smaller N, because a draw here costs about ten seconds.",
    "SEED: the draws come from numpy default_rng at seed %d; without"
    % SEED,
    "      it the file does not reproduce its own control.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    small = read_smallmins()
    prevmin = small[-1]
    for n, v in zip(NS, pub):
        say("READ audit_directidentity_null.txt %d %.4f" % (n, v))
    for v in small:
        say("READ audit_directidentity_deep.txt multiplicative %.4f"
            % v)
    say("  mu's T/(S(N)N) at the two largest N, and the three smaller")
    say("  N's multiplicative minima, both read from those files")
    say("  the one to beat is %.4f" % prevmin)
    say("PRINTBOUND audit_directidentity_bigN %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d draws in %d passes of %d" % (DRAWS, PASSES, PERPASS))

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
        say("  N = %-9d #k = %-9d S = %.6f" % (N, kk.size, S))

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

    # -------------------------------------------------------------- N1
    say()
    say("N1  does this construction reproduce mu's ratios?")
    ok = True
    for N, r, v in zip(NS, real, pub):
        g = abs(round(r, DEC) - round(v, DEC)) < 10.0 ** (-DEC) / 2
        ok &= g
        say("  N = %-9d here %.4f against its %.4f  %s"
            % (N, r, v, "ok" if g else "MISMATCH"))
    say("  N1 %s   (cap: %d decimals at both)"
        % ("hold" if ok else "REFUTED", DEC))
    if not ok:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    pool = np.zeros((DRAWS, len(NS)))
    row = 0
    for k in range(PASSES):
        t = np.zeros(NMAX + 1, dtype=np.uint64)
        hi = rng.integers(0, 1 << 32, size=len(pr), dtype=np.uint64)
        lo = rng.integers(0, 1 << 32, size=len(pr), dtype=np.uint64)
        pm = (hi << np.uint64(32)) | lo
        for i, p in enumerate(pr):
            t[int(p)::int(p)] ^= pm[i]
        del hi, lo, pm
        for d in range(PERPASS):
            bit = (t >> np.uint64(d)) & np.uint64(1)
            sg = np.where(sqf, 1.0 - 2.0 * bit.astype(np.float64), 0.0)
            pool[row] = ratio_of(sg)
            row += 1
            del sg, bit
        del t
        say("  pass %d of %d, %d draws" % (k + 1, PASSES, row))

    say()
    say("           N        mu    draw min   draw median"
        "    draws reaching mu")
    reach, mins = [], []
    for i, N in enumerate(NS):
        c = int((pool[:, i] <= real[i]).sum())
        reach.append(c)
        mins.append(float(pool[:, i].min()))
        say("      %8d  %.4f  %9.4f  %11.4f  %14d"
            % (N, real[i], mins[i], float(np.median(pool[:, i])), c))
        say("POINT identbig_%d %.5f" % (N, mins[i]))
    say("SCALES %d" % len(NS))

    # -------------------------------------------------------------- N2
    say()
    say("N2  does the separation hold where it has never been looked?")
    n2 = all(c == 0 for c in reach)
    say("  draws reaching mu: %s of %d at each N"
        % (", ".join(str(c) for c in reach), DRAWS))
    say("COUNT identbig_reach %d" % sum(reach))
    say("  N2 %s   (cap: none at either N)"
        % ("hold" if n2 else "REFUTED"))

    # -------------------------------------------------------------- N3
    say()
    say("N3  does the separation widen with N?")
    say("  the three smaller N gave %s"
        % ", ".join("%.4f" % v for v in small))
    say("  here: %s" % ", ".join("%.4f" % m for m in mins))
    n3 = mins[0] > prevmin and mins[1] > mins[0]
    say("  N3 %s   (cap: each above its predecessor)"
        % ("hold" if n3 else "REFUTED"))

    # -------------------------------------------------------------- N4
    say()
    say("N4  what does depth buy at the low end?")
    say("      draws        %s" % "        ".join("N=%d" % N
                                                  for N in NS))
    drift = 0.0
    first = None
    for c in CHECKS:
        row_ = [float(pool[:c, i].min()) for i in range(len(NS))]
        if first is None:
            first = row_
        drift = max(drift, max(abs(a - b) for a, b in zip(row_, first)))
        say("      %5d    %s" % (c, "    ".join("%.4f" % x
                                                for x in row_)))
    say("SCALES %d" % len(CHECKS))
    say("  the minimum drifts by at most %.4f between %d and %d draws"
        % (drift, min(CHECKS), max(CHECKS)))
    say("SPREAD identbig_depth %.4f" % drift)
    n4 = drift <= DEPTHTOL
    say("  N4 %s   (cap: %.1f)" % ("hold" if n4 else "REFUTED",
                                   DEPTHTOL))
    if not n4:
        say("  N4 failed, so N3's rise is partly the depth and is "
            "recorded as")
        say("  unmeasured rather than as shown")

    say()
    say("=" * 70)
    say("N1 %s  N2 %s  N3 %s  N4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (ok, n2, n3, n4)))
    say()
    if n2:
        say("the last two N go the way the three smaller ones did, so "
            "the identity")
        say("separation holds across the whole range rem:identitynull "
            "publishes and")
        say("item 8's work list is finished. The bound is the "
            "smallest ratio seen")
        say("at each N and the depth here is half what the smaller N "
            "had.")
    else:
        say("the separation fails at the larger N, so it is a "
            "property of small N")
        say("and rem:identitydeep's claim has to be cut to the range "
            "that shows it.")
        say("The counts above are where and how often.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
