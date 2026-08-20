# -*- coding: utf-8 -*-
r"""
The last untested block: rem:identitynull against a multiplicative null

WHAT IS AT STAKE

rem:coinsurface left two blocks outside the C(N) branch standing on
eight draws.  rem:decaydeep and rem:decaysweep took one of them and
both of its halves went the same way: the separation survives an iid
coin at full depth and fails against a random multiplicative ensemble.
**rem:identitynull is the last one that has not been asked.**

Its separation is the largest anywhere in this paper.  T/(S(N)N) is
1.0039, 0.9865, 0.9893 for mu at the three smallest N and 4.9835 to
5.6877 across the eight published draws -- a factor of five, not a
tail.

And unlike the decay block there is a mechanism, which is why this run
is worth making rather than assuming.  Writing v = N - i,

    T = sum_v Lambda(N-v) s(v) sum_{k | v, k >= 2, squarefree}
        log(k) s(k),

and for s = mu the inner sum is the classical -Lambda(v), which
collapses T to R exactly -- the published |T-R|/R is 3.294e-16.  No
other sign function has that identity.  **But a random multiplicative
f is not structureless either**: sum_{d | n} f(d) = prod_{p | n}
(1 + f(p)), which vanishes whenever any f(p) = -1, so a multiplicative
draw has a cancellation an iid draw has none of.  Whether that
cancellation is enough to bring T down by a factor of five is the
question, and the algebra does not answer it.

Three N rather than five, because the k-range here is every squarefree
k below N and one draw over the three smallest costs 2.6 seconds
against roughly nine over all five.  512 draws of each null over three
N is under an hour; over five it is hours.  The two largest N are
therefore untested and that is stated rather than hidden.

Field, k-range, S(N) and the twin constant are taken unchanged from
audit_directidentity_null.py.

BACKS: Remark {#rem:identitydeep} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  L1  THE GATE.  mu reproduces both published quantities: |T-R|/R
      below 1e-12 at all three N, and T/(S(N)N) equal to 1.0039,
      0.9865, 0.9893 to four decimals.
  L2  The published claim survives depth: none of 512 iid draws
      reaches mu, that is has |T|/(S(N)N) at or below mu's value, at
      any of the three N.
  L3  **And the multiplicative null closes it**: at each of the three
      N, at least one of 512 multiplicative draws reaches mu.  This
      is the deflationary prediction and the one this branch's
      history points at.
  L4  The multiplicative ensemble is the lower of the two: its
      smallest |T|/(S(N)N) is below the iid ensemble's at every N.

REFUTATION RULE (fixed before the run)

  L1  REFUTED outside those tolerances; nothing below is reported.
  L2  REFUTED if any iid draw reaches mu.  Then eight draws were
      hiding a tail and the block's own control was too small,
      independently of the null type -- a plainer fault than L3 finds
      and the one to report first.
  L3  **REFUTED if no multiplicative draw reaches mu at some N.**
      Then the identity separation survives the right null, and it is
      the only thing in this repository that does -- larger than
      anything the C(N) branch produced, since everything there
      closed, and larger than the decay block, which did not.  It
      would need another seed before being leaned on and the remark
      must say so, together with the smallest ratio any draw achieved
      as the bound this run actually gives.
  L4  REFUTED if the iid ensemble goes as low or lower at any N.
      **The unresolved case is named**: these are two minima over 512
      draws each, single order statistics, so a gap smaller than the
      ensembles' own spread is not a gap -- if they are within each
      other's spread the reading is "not resolved", never "the same".
      L4 is a direction and it gates nothing.

  WHAT THIS CANNOT DO.  Three of five N, one seed per ensemble, 512
  draws of a heavy-tailed object.  If L3 is refuted the bound is the
  smallest ratio seen and not zero, and it is printed.
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
OUT = os.path.join(ROOT, "results", "audit_directidentity_deep.txt")
SRC = os.path.join(ROOT, "results", "audit_directidentity_null.txt")

NS = [200_000, 400_000, 800_000]
NMAX = max(NS)
CLIM = 4_000_000
PASSES = 8
PERPASS = 64
DRAWS = PASSES * PERPASS
SEED = 20260906
DEC = 4
TOLID = 1e-12


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


def read_pub():
    """mu's T/(S N) at the three N, from that remark's own file"""
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
    "STATISTIC: the untruncated sum",
    "           T = sum_k (log k) s(k) A_s(N;k) over every squarefree",
    "           k with 2 <= k < N, as the ratio |T|/(S(N)N), for mu",
    "           and for two ensembles of %d sign patterns -- iid on" % DRAWS,
    "           supp(mu^2), and random multiplicative with f(m) the",
    "           product of f(p) over p | m on squarefree m; and mu's",
    "           |T-R|/R, which the Moebius identity makes zero.",
    "FIELD: N = %s; k over every squarefree 2 <= k < N; S(N) from"
    % (NS,),
    "       the twin product with the odd prime factors of N, at the",
    "       fixed bound %d. Field, k-range, S(N) and the twin" % CLIM,
    "       constant are taken unchanged from",
    "       code/audit_directidentity_null.py so the numbers are",
    "       comparable. mu's published ratios are READ from",
    "       results/audit_directidentity_null.txt and re-measured here",
    "       as the gate. The two largest N of that file are not run.",
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
    for n, v in zip(NS, pub):
        say("READ audit_directidentity_null.txt %d %.4f" % (n, v))
    say("  mu's T/(S(N)N) at the three N, read from that file")
    say("PRINTBOUND audit_directidentity_deep %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d draws of each null, in %d passes of %d"
        % (DRAWS, PASSES, PERPASS))

    say("sieving to %d" % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    muf = mu.astype(np.float64)
    isp = np.zeros(NMAX + 1, dtype=np.float64)
    isp[primes_upto(NMAX)] = 1.0
    del mu
    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    say("  twin constant %.6f, %d primes to %d" % (twin, len(pr), NMAX))

    ks, lgs, idxs, Ss, Rs = {}, {}, {}, {}, {}
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
        Rs[N] = float((lam[1:N] * lam[N - 1:0:-1]
                       * isp[N - 1:0:-1]).sum())
        say("  N = %-9d #k = %-9d S = %.6f  R/N = %.4f"
            % (N, kk.size, S, Rs[N] / N))

    def tot_of(sg):
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
            out.append(s)
        return out

    real = tot_of(muf)
    ratio = [abs(t) / (Ss[N] * N) for t, N in zip(real, NS)]
    ident = [abs(t - Rs[N]) / Rs[N] for t, N in zip(real, NS)]

    # -------------------------------------------------------------- L1
    say()
    say("L1  does this construction reproduce both published "
        "quantities?")
    ok = True
    for N, r, d, v in zip(NS, ratio, ident, pub):
        g1 = d < TOLID
        g2 = abs(round(r, DEC) - round(v, DEC)) < 10.0 ** (-DEC) / 2
        ok &= g1 and g2
        say("  N = %-9d |T-R|/R %.3e   T/(S N) here %.4f against its "
            "%.4f  %s" % (N, d, r, v, "ok" if g1 and g2 else "MISMATCH"))
    say("  L1 %s   (cap: %.0e on the identity, %d decimals on the "
        "ratio)" % ("hold" if ok else "REFUTED", TOLID, DEC))
    if not ok:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    pools = {}
    for kind in ("iid", "multiplicative"):
        pool = np.zeros((DRAWS, len(NS)))
        row = 0
        for k in range(PASSES):
            if kind == "iid":
                hi = rng.integers(0, 1 << 32, size=NMAX + 1,
                                  dtype=np.uint64)
                lo = rng.integers(0, 1 << 32, size=NMAX + 1,
                                  dtype=np.uint64)
                t = (hi << np.uint64(32)) | lo
                del hi, lo
            else:
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
                tt = tot_of(sg)
                pool[row] = [abs(x) / (Ss[N] * N)
                             for x, N in zip(tt, NS)]
                row += 1
                del sg, bit
            del t
            say("  %s pass %d of %d, %d draws"
                % (kind, k + 1, PASSES, row))
        pools[kind] = pool

    say()
    say("      null              N        mu    draw min   draw median"
        "    draws reaching mu")
    reach = {}
    for kind in ("iid", "multiplicative"):
        pool = pools[kind]
        cnt = []
        for i, N in enumerate(NS):
            c = int((pool[:, i] <= ratio[i]).sum())
            cnt.append(c)
            say("      %-16s %8d  %.4f  %9.4f  %11.4f  %14d"
                % (kind, N, ratio[i], float(pool[:, i].min()),
                   float(np.median(pool[:, i])), c))
            say("POINT identdeep_%s_%d %.5f"
                % (kind[:4], N, float(pool[:, i].min())))
        reach[kind] = cnt
    say("SCALES %d" % len(NS))

    # -------------------------------------------------------------- L2
    say()
    say("L2  does the published claim survive depth?")
    l2 = all(c == 0 for c in reach["iid"])
    say("  iid draws reaching mu: %s of %d at each N"
        % (", ".join(str(c) for c in reach["iid"]), DRAWS))
    say("COUNT identdeep_iid_reach %d" % sum(reach["iid"]))
    say("  L2 %s   (cap: none at any N)"
        % ("hold" if l2 else "REFUTED"))

    # -------------------------------------------------------------- L3
    say()
    say("L3  does the multiplicative null close it?")
    l3 = all(c > 0 for c in reach["multiplicative"])
    say("  multiplicative draws reaching mu: %s of %d at each N"
        % (", ".join(str(c) for c in reach["multiplicative"]), DRAWS))
    say("COUNT identdeep_mult_reach %d"
        % sum(reach["multiplicative"]))
    say("  the smallest ratio any multiplicative draw reached: %s"
        % ", ".join("%.4f" % float(pools["multiplicative"][:, i].min())
                    for i in range(len(NS))))
    say("  L3 %s   (cap: at least one at every N)"
        % ("hold" if l3 else "REFUTED"))

    # -------------------------------------------------------------- L4
    say()
    say("L4  which ensemble goes lower?")
    l4 = True
    for i, N in enumerate(NS):
        a = float(pools["multiplicative"][:, i].min())
        b = float(pools["iid"][:, i].min())
        sa = float(pools["multiplicative"][:, i].std(ddof=1))
        sb = float(pools["iid"][:, i].std(ddof=1))
        l4 &= a < b
        say("  N = %-9d multiplicative min %.4f, iid min %.4f, gap "
            "%+.4f" % (N, a, b, b - a))
        say("  their spreads are %.4f and %.4f, so a gap below about "
            "%.4f is" % (sa, sb, max(sa, sb)))
        say("  not resolved by this design")
    say("  L4 %s   (cap: multiplicative lower at every N)"
        % ("hold" if l4 else "REFUTED"))

    say()
    say("=" * 70)
    say("L1 %s  L2 %s  L3 %s  L4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (ok, l2, l3, l4)))
    say()
    if l2 and l3:
        say("the last block goes the way the others did. The identity "
            "separation")
        say("survives an iid coin at full depth and a random "
            "multiplicative")
        say("ensemble reaches it, so every place rem:coinsurface "
            "named now rests")
        say("on a control of the wrong kind, and none of them on one "
            "too small.")
    elif l2 and not l3:
        say("the identity separation survives the right null. It is "
            "the only")
        say("thing in this repository that a random multiplicative "
            "ensemble has")
        say("failed to reach, and the mechanism is named: no sign "
            "function but mu")
        say("has sum_{k|v} s(k) log k = -Lambda(v). One seed and 512 "
            "draws; the")
        say("smallest ratio any draw achieved is printed above as the "
            "bound this")
        say("run gives, and another seed is owed before it is leaned "
            "on.")
    else:
        say("iid draws reach mu once there are enough of them, so "
            "that block's")
        say("own control was too small and the null type is a second "
            "question")
        say("rather than the first. The counts above are how badly.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
