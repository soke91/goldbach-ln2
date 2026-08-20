# -*- coding: utf-8 -*-
r"""
rem:decaynull's separation, at depth and against the right null

WHAT IS AT STAKE

rem:coinsurface counted what this paper stands on having escaped a
coin: eleven blocks, eight of them in the C(N) branch that found the
problem, one a false hit, and two outside -- rem:identitynull and
rem:decaynull.  **Both rest on ensembles of eight draws.**
rem:cnmultdeep measured a multiplicative ensemble's own spread more
than doubling between 32 draws and 512, and rem:cncoindeep found an
iid coin narrower than a multiplicative one at full resolution.  Eight
draws is below both of those numbers, so neither block has been asked
the question this branch learned to ask.

This asks it of rem:decaynull, whose claim is a separation in size:
mu's |1/2 - f(N)| reads 0.2727 and 0.2772 at the two smallest N
against a coin maximum of 0.0826 over the whole sweep.  Two things are
changed and nothing else: the draw count goes from 8 to 512, and a
random multiplicative ensemble is run beside the iid one.

**The band is cut to the two smallest N and that is a real
restriction.**  The statistic costs a strided sum for every squarefree
k < N^0.56 at every N and every draw, so 512 draws over the published
six-point sweep is hours; over N = 2e5 and 4e5 it is minutes.  What
this run can therefore say is whether the separation survives depth
and the null type *at those two N*, which is where mu's two largest
values sit.  It says nothing about the four larger N and nothing about
the alpha sweep, which is a different claim of that remark.

The field, weight, k-range and theta' are taken from
audit_decayfamily_null.py unchanged, so that the numbers are
comparable rather than merely similar.

BACKS: Remark {#rem:decaydeep} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  J1  THE GATE.  mu's |1/2 - f| reproduces rem:decaynull's published
      0.2727 and 0.2772 at N = 2e5 and 4e5, to four decimals.  If the
      construction here is not that one, nothing below compares.
  J2  **The published separation survives depth.**  At 512 iid draws,
      none reaches mu at either N.
  J3  And it survives the null type: at 512 random multiplicative
      draws, none reaches mu at either N.
  J4  The multiplicative ensemble is the wider of the two here, as
      rem:cncoindeep found for C(N): its largest |1/2 - f| over draws
      exceeds the iid ensemble's at both N.

REFUTATION RULE (fixed before the run)

  J1  REFUTED outside four decimals; nothing below is reported.
  J2  **REFUTED if any of the 512 iid draws reaches mu.**  Then eight
      draws were hiding a tail and rem:decaynull's separation is not
      what it was recorded as -- the block would need its numbers
      re-stated, and the count of draws that reach is the measure of
      how badly.
  J3  REFUTED on the same terms for the multiplicative ensemble.  J2
      and J3 can part: iid holding while multiplicative fails is
      exactly the shape rem:cncoindeep found, and would say the null
      type matters here too and that the block's control was the
      wrong one rather than too small.
  J4  REFUTED if the iid ensemble is as wide or wider at either N.
      **The unresolved case is named**: with 512 draws each, the two
      maxima are single order statistics and a difference smaller
      than their own scatter is not a difference -- if the two are
      within each other's spread the reading is "not resolved", never
      "the same width".  J4 is a direction, not a measurement, and it
      gates nothing.

  WHAT THIS CANNOT DO.  Two of six N, one statistic of the two that
  remark makes, and 512 draws of a heavy-tailed object.  A separation
  that survives here is not shown to survive at the four larger N,
  and nothing here touches rem:identitynull, whose separation is a
  factor of five and rests on an identity rather than on a tail.
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
OUT = os.path.join(ROOT, "results", "audit_decayfamily_deep.txt")
SRC = os.path.join(ROOT, "results", "audit_decayfamily_null.txt")

NS = [200_000, 400_000]
THETA = 0.56
NMAX = max(NS)
PASSES = 8
PERPASS = 64
DRAWS = PASSES * PERPASS
SEED = 20260904
DEC = 4


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


def read_pub():
    """mu's |1/2-f| at the two smallest N, from that remark's file"""
    src = io.open(SRC, encoding="utf-8").read()
    out = []
    for n in NS:
        m = re.search(r"^\s*%d\s+([\d.]+)\s" % n, src, re.M)
        if not m:
            raise SystemExit("no row for N = %d in %s" % (n, SRC))
        out.append(float(m.group(1)))
    return out


HEAD = [
    "STATISTIC: |1/2 - f(N)| for the decay family at the two smallest",
    "           N of rem:decaynull's sweep, for mu and for two",
    "           ensembles of %d sign patterns -- iid on supp(mu^2)," % DRAWS,
    "           and random multiplicative, f(p) = +-1 iid with",
    "           f(m) the product over p | m on squarefree m.",
    "FIELD: N in %s with theta' = %.2f, k over the squarefree" % (NS, THETA),
    "       k < N^theta' coprime to N, the field, weight and k-range",
    "       taken unchanged from code/audit_decayfamily_null.py so the",
    "       numbers are comparable. mu's two values are READ from",
    "       results/audit_decayfamily_null.txt and re-measured here as",
    "       the gate.",
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
        say("READ audit_decayfamily_null.txt %d %.4f" % (n, v))
    say("  mu's |1/2 - f| at the two smallest N, read from that file")
    say("PRINTBOUND audit_decayfamily_deep %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d draws of each null, in %d passes of %d"
        % (DRAWS, PASSES, PERPASS))

    say("sieving to %d" % NMAX)
    lam, mu = sieves(NMAX)
    sqf = mu != 0
    muf = mu.astype(np.float64)
    del mu
    pr = primes_upto(NMAX)
    say("  %d squarefree m carry the signs, %d primes"
        % (int(sqf.sum()), len(pr)))

    ks, lgs, idxs = {}, {}, {}
    for N in NS:
        PN = factor_set(N)
        K = int(N ** THETA)
        kk = np.array([k for k in range(2, K)
                       if sqf[k] and all(k % q for q in PN)])
        ks[N] = kk
        lgs[N] = np.log(kk.astype(float))
        idxs[N] = np.arange(1, N, dtype=np.int64)
        say("  N = %-8d  K = %-8d  #k = %d" % (N, K, kk.size))

    def dev_of(sg):
        out = []
        for N in NS:
            f0 = np.zeros(N, dtype=np.float64)
            f0[1:] = lam[1:N] * sg[N - idxs[N]]
            kk = ks[N]
            A = np.empty(kk.size)
            for i, k in enumerate(kk):
                k = int(k)
                r = N % k
                A[i] = f0[r::k].sum() if r else f0[k::k].sum()
            del f0
            H = sg[kk] * A
            w = lgs[N] * np.abs(H)
            fr = float(w[H > 0].sum() / w.sum())
            out.append(abs(0.5 - fr))
        return out

    real = dev_of(muf)

    # -------------------------------------------------------------- J1
    say()
    say("J1  does this construction reproduce the published values?")
    ok = True
    for n, v, r in zip(NS, pub, real):
        good = abs(round(r, DEC) - round(v, DEC)) < 10.0 ** (-DEC) / 2
        ok &= good
        say("  N = %-8d here %.4f against its %.4f  %s"
            % (n, r, v, "ok" if good else "MISMATCH"))
    say("  J1 %s   (cap: %d decimals at both)"
        % ("hold" if ok else "REFUTED", DEC))
    if not ok:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    pools = {}
    for kind in ("iid", "multiplicative"):
        pool = np.zeros((DRAWS, len(NS)))
        row = 0
        for _ in range(PASSES):
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
                bit = ((t >> np.uint64(d)) & np.uint64(1))
                sg = np.where(sqf, 1.0 - 2.0 * bit.astype(np.float64),
                              0.0)
                pool[row] = dev_of(sg)
                row += 1
                del sg, bit
            del t
        pools[kind] = pool
        say("  %s ensemble done, %d draws" % (kind, row))

    say()
    say("      null              N      mu      draw max   draw median"
        "    draws reaching mu")
    reach = {}
    for kind in ("iid", "multiplicative"):
        pool = pools[kind]
        cnt = []
        for i, N in enumerate(NS):
            c = int((pool[:, i] >= real[i]).sum())
            cnt.append(c)
            say("      %-16s %8d  %.4f  %9.4f  %11.4f  %14d"
                % (kind, N, real[i], float(pool[:, i].max()),
                   float(np.median(pool[:, i])), c))
            say("POINT deepdecay_%s_%d %.5f"
                % (kind[:4], N, float(pool[:, i].max())))
        reach[kind] = cnt
    say("SCALES %d" % len(NS))

    # -------------------------------------------------------------- J2
    say()
    say("J2  does the separation survive depth, against iid?")
    j2 = all(c == 0 for c in reach["iid"])
    say("  draws reaching mu: %s of %d at each N"
        % (", ".join(str(c) for c in reach["iid"]), DRAWS))
    say("COUNT deepdecay_iid_reach %d" % sum(reach["iid"]))
    say("  J2 %s   (cap: none at either N)"
        % ("hold" if j2 else "REFUTED"))

    # -------------------------------------------------------------- J3
    say()
    say("J3  does it survive the multiplicative null?")
    j3 = all(c == 0 for c in reach["multiplicative"])
    say("  draws reaching mu: %s of %d at each N"
        % (", ".join(str(c) for c in reach["multiplicative"]), DRAWS))
    say("COUNT deepdecay_mult_reach %d" % sum(reach["multiplicative"]))
    say("  J3 %s   (cap: none at either N)"
        % ("hold" if j3 else "REFUTED"))

    # -------------------------------------------------------------- J4
    say()
    say("J4  which null is wider here?")
    j4 = True
    for i, N in enumerate(NS):
        a = float(pools["multiplicative"][:, i].max())
        b = float(pools["iid"][:, i].max())
        sa = float(pools["multiplicative"][:, i].std(ddof=1))
        sb = float(pools["iid"][:, i].std(ddof=1))
        j4 &= a > b
        say("  N = %-8d multiplicative max %.4f, iid max %.4f, "
            "gap %+.4f" % (N, a, b, a - b))
        say("  their spreads are %.4f and %.4f, so a gap below about "
            "%.4f is" % (sa, sb, max(sa, sb)))
        say("  not resolved by this design")
    say("  J4 %s   (cap: multiplicative wider at both N)"
        % ("hold" if j4 else "REFUTED"))

    say()
    say("=" * 70)
    say("J1 %s  J2 %s  J3 %s  J4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (ok, j2, j3, j4)))
    say()
    if j2 and j3:
        say("rem:decaynull's separation is not an eight-draw "
            "artefact. At the two")
        say("N this run reaches, neither 512 iid draws nor 512 random "
            "multiplicative")
        say("ones come to mu, so the block survives both the depth "
            "and the null")
        say("type that rem:coinsurface put in question. The four "
            "larger N are")
        say("untouched and so is that remark's alpha sweep.")
    elif j2 and not j3:
        say("the null type matters here as it did for C(N). Depth "
            "alone leaves the")
        say("separation standing and the multiplicative ensemble does "
            "not, so")
        say("rem:decaynull's control was the wrong one rather than "
            "too small.")
    else:
        say("eight draws were hiding a tail. The separation "
            "rem:decaynull recorded")
        say("is reached by draws of its own null once there are "
            "enough of them,")
        say("and the counts above are how badly.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
