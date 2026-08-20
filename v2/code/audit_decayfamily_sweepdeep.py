# -*- coding: utf-8 -*-
r"""
rem:decaynull's other half: the alpha sweep, at depth and both nulls

WHAT IS AT STAKE

rem:decaydeep tested one of rem:decaynull's two claims and flipped it.
The size separation -- mu's |1/2-f| against the draws' -- survives 512
iid draws untouched and fails against 512 random multiplicative ones,
67 and 54 of them reaching mu at the two smallest N.  It said in so
many words what it had not done: only two of six N, and nothing at all
about that remark's second and stronger claim, the alpha sweep.

That claim is the sharper of the two.  Fitting
log|1/2-f| = a - c (log N)^alpha over alpha in 0.05 to 1.50, mu gives
an interior minimiser at 1.45 with residual 0.011556, while **not one**
of the eight published draws has an interior minimiser at all, every
one pinning at a grid end with residuals 2.991540 to 10.078929 -- two
to three orders in fit quality.  If that survives a multiplicative
ensemble it is the strongest thing outside the C(N) branch; if it does
not, both halves of the block go the same way and the pattern
rem:cncoindeep found holds wherever it has been looked for.

The cost objection rem:decaydeep raised was wrong and is corrected
here.  One draw over the whole six-point sweep takes under half a
second -- the k-range is squarefree and coprime to N, which leaves
2186 of them at the top N rather than the 6570 that N^0.56 suggests --
so 512 draws of each null over all six N is minutes, and this run does
the size comparison at all six N as well.

The field, weight, k-range, theta', alpha grid and the sweep's own
one-per-cent band rule are taken from audit_decayfamily_null.py
unchanged.

BACKS: Remark {#rem:decaysweep} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  K1  THE GATE.  mu reproduces that remark's published sweep exactly:
      an interior minimiser at alpha* = 1.45 with residual 0.011556 to
      six decimals, and |1/2-f| = 0.2727 at the smallest N to four.
  K2  **The size separation fails against the multiplicative null at
      every N**, not only the two rem:decaydeep reached: at each of
      the six, at least one of 512 multiplicative draws reaches mu.
  K3  **And the sweep separation fails too**: at least one
      multiplicative draw has an interior minimiser with residual at
      or below mu's 0.011556.  This is the deflationary prediction and
      the one this branch's history points at.
  K4  While the published iid claim survives depth: none of 512 iid
      draws has an interior minimiser with residual at or below mu's.

REFUTATION RULE (fixed before the run)

  K1  REFUTED outside those tolerances; nothing below is reported,
      because the sweep is a fit and a fit that does not reproduce is
      not the same fit.
  K2  REFUTED if any N has no multiplicative draw reaching mu.  Then
      the size separation is real at some N and not at others, and
      which N is the finding -- rem:decaydeep only reached the two
      smallest and mu's values fall from 0.2727 to 0.1624 across the
      sweep, so the larger N are where it could survive.
  K3  **REFUTED if no multiplicative draw manages both.**  Then the
      sweep separation is real against the right null and it is the
      strongest surviving claim outside the C(N) branch -- larger than
      anything that branch produced, since everything there closed.
      It would need repeating at another seed before being leaned on,
      and the remark must say so.
  K4  REFUTED if any iid draw manages both.  Then rem:decaynull's
      published sweep claim was an eight-draw artefact independently
      of the null type, which would be a plainer fault than K3 finds
      and would have to be reported first.

  K3 and K4 can both fail, and that combination means the sweep claim
  was thin in both directions at once; the counts are printed either
  way so the two can be told apart.

  WHAT THIS CANNOT DO.  One seed per ensemble and 512 draws of a
  heavy-tailed object.  A count of zero here bounds a tail no better
  than 512 draws can, and the largest residual seen is printed so the
  bound is visible rather than implied.
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
                   "audit_decayfamily_sweepdeep.txt")
SRC = os.path.join(ROOT, "results", "audit_decayfamily_null.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
ALPHAS = [round(0.05 * i, 2) for i in range(1, 31)]
NMAX = max(NS)
PASSES = 8
PERPASS = 64
DRAWS = PASSES * PERPASS
SEED = 20260905
DEC4 = 4
DEC6 = 6


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


def sweep(dev):
    """RSS of log dev = a - c (log N)^alpha, as that script does it"""
    y = np.log(np.array(dev))
    L = np.log(np.array(NS, dtype=float))
    rss = []
    for a in ALPHAS:
        x = L ** a
        c = np.polyfit(x, y, 1)
        rss.append(float(((y - (c[0] * x + c[1])) ** 2).sum()))
    rss = np.array(rss)
    j = int(np.argmin(rss))
    interior = 0 < j < len(ALPHAS) - 1
    keep = rss <= rss[j] * 1.01
    band = (ALPHAS[int(np.flatnonzero(keep)[-1])]
            - ALPHAS[int(np.flatnonzero(keep)[0])])
    return ALPHAS[j], band, interior, float(rss[j])


def read_pub():
    src = io.open(SRC, encoding="utf-8").read()
    devs = []
    for n in NS:
        m = re.search(r"^\s*%d\s+([\d.]+)\s" % n, src, re.M)
        if not m:
            raise SystemExit("no row for N = %d in %s" % (n, SRC))
        devs.append(float(m.group(1)))
    m = re.search(r"^\s*mu\s+([\d.]+)\s+([\d.]+)\s+True\s+([\d.]+)\s*$",
                  src, re.M)
    if not m:
        raise SystemExit("no mu sweep row in %s" % SRC)
    return devs, float(m.group(1)), float(m.group(3))


HEAD = [
    "STATISTIC: |1/2 - f(N)| at six N and the alpha sweep of",
    "           log|1/2-f| = a - c (log N)^alpha over that series, for",
    "           mu and for two ensembles of %d sign patterns -- iid" % DRAWS,
    "           on supp(mu^2), and random multiplicative with",
    "           f(m) the product of f(p) over p | m on squarefree m --",
    "           counting draws that reach mu in size at each N and",
    "           draws whose sweep has an interior minimiser at or",
    "           below mu's residual.",
    "FIELD: N = 2e5 through 6.4e6 by doubling, theta' = %.2f, k over"
    % THETA,
    "       the squarefree k < N^theta' coprime to N, alpha over 0.05",
    "       to 1.50 in steps of 0.05 with the one-per-cent band rule;",
    "       field, weight, k-range and sweep taken unchanged from",
    "       code/audit_decayfamily_null.py. mu's published series and",
    "       sweep are READ from results/audit_decayfamily_null.txt and",
    "       re-measured here as the gate.",
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

    pubdev, puba, pubrss = read_pub()
    for n, v in zip(NS, pubdev):
        say("READ audit_decayfamily_null.txt %d %.4f" % (n, v))
    say("READ audit_decayfamily_null.txt mu %.2f" % puba)
    say("  mu's published series and its sweep minimiser, read from "
        "that file")
    say("PRINTBOUND audit_decayfamily_sweepdeep %d %.8f"
        % (DEC6, 0.5 * 10.0 ** (-DEC6)))
    say("  %d draws of each null, in %d passes of %d"
        % (DRAWS, PASSES, PERPASS))

    say("sieving to %d" % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    muf = mu.astype(np.float64)
    del mu
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
    say("  #k per N: %s" % ", ".join(str(ks[N].size) for N in NS))

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
            out.append(max(abs(0.5 - fr), 1e-12))
        return out

    real = dev_of(muf)
    ra, rb, rint, rrss = sweep(real)

    # -------------------------------------------------------------- K1
    say()
    say("K1  does this construction reproduce the published sweep?")
    okdev = abs(round(real[0], DEC4) - round(pubdev[0], DEC4)) \
        < 10.0 ** (-DEC4) / 2
    oka = rint and abs(ra - puba) < 1e-9
    okr = abs(round(rrss, DEC6) - round(pubrss, DEC6)) \
        < 10.0 ** (-DEC6) / 2
    say("  |1/2-f| at the smallest N: here %.4f against its %.4f"
        % (real[0], pubdev[0]))
    say("  sweep: alpha* %.2f interior %s residual %.6f against its "
        "%.2f and %.6f" % (ra, rint, rrss, puba, pubrss))
    k1 = okdev and oka and okr
    say("  K1 %s   (cap: %d decimals on the series, %d on the "
        "residual)" % ("hold" if k1 else "REFUTED", DEC4, DEC6))
    if not k1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    res = {}
    for kind in ("iid", "multiplicative"):
        devs = np.zeros((DRAWS, len(NS)))
        swp = np.zeros((DRAWS, 2))
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
                bit = (t >> np.uint64(d)) & np.uint64(1)
                sg = np.where(sqf,
                              1.0 - 2.0 * bit.astype(np.float64), 0.0)
                dv = dev_of(sg)
                devs[row] = dv
                a_, b_, in_, r_ = sweep(dv)
                swp[row] = (1.0 if in_ else 0.0, r_)
                row += 1
                del sg, bit
            del t
        res[kind] = (devs, swp)
        say("  %s ensemble done, %d draws" % (kind, row))

    # -------------------------------------------------------------- K2
    say()
    say("K2  does the size separation fail at every N?")
    say("      null              N       mu     draw max"
        "   draws reaching mu")
    k2 = True
    for kind in ("iid", "multiplicative"):
        devs = res[kind][0]
        for i, N in enumerate(NS):
            c = int((devs[:, i] >= real[i]).sum())
            if kind == "multiplicative":
                k2 &= c > 0
            say("      %-16s %8d  %.4f   %.4f  %14d"
                % (kind, N, real[i], float(devs[:, i].max()), c))
            say("POINT sweepdeep_%s_%d %.5f"
                % (kind[:4], N, float(devs[:, i].max())))
    say("SCALES %d" % len(NS))
    say("  K2 %s   (cap: a multiplicative draw at every N)"
        % ("hold" if k2 else "REFUTED"))

    # -------------------------------------------------------------- K3
    say()
    say("K3  does the sweep separation fail against the "
        "multiplicative null?")
    say("      null              interior   residual <= mu's"
        "   both     best residual")
    counts = {}
    for kind in ("iid", "multiplicative"):
        swp = res[kind][1]
        ni = int(swp[:, 0].sum())
        nr = int((swp[:, 1] <= rrss).sum())
        nb = int(((swp[:, 0] > 0) & (swp[:, 1] <= rrss)).sum())
        counts[kind] = (ni, nr, nb, float(swp[:, 1].min()))
        say("      %-16s %8d   %16d   %4d   %13.6f"
            % (kind, ni, nr, nb, float(swp[:, 1].min())))
        say("POINT sweepdeep_best_%s %.6f" % (kind[:4],
                                              float(swp[:, 1].min())))
    k3 = counts["multiplicative"][2] > 0
    say("COUNT sweepdeep_mult_both %d" % counts["multiplicative"][2])
    say("  mu's own residual is %.6f at an interior alpha* of %.2f"
        % (rrss, ra))
    say("  K3 %s   (cap: at least one multiplicative draw with both)"
        % ("hold" if k3 else "REFUTED"))

    # -------------------------------------------------------------- K4
    say()
    say("K4  does the published iid claim survive depth?")
    k4 = counts["iid"][2] == 0
    say("  iid draws with an interior minimiser at or below mu's "
        "residual: %d of %d" % (counts["iid"][2], DRAWS))
    say("COUNT sweepdeep_iid_both %d" % counts["iid"][2])
    say("  K4 %s   (cap: none)" % ("hold" if k4 else "REFUTED"))

    say()
    say("=" * 70)
    say("K1 %s  K2 %s  K3 %s  K4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (k1, k2, k3, k4)))
    say()
    if k3 and k2:
        say("both halves of rem:decaynull go the same way. Size and "
            "sweep alike")
        say("survive an iid coin at full depth and fail against a "
            "random")
        say("multiplicative ensemble, so that block's control was the "
            "wrong kind")
        say("throughout and the pattern rem:cncoindeep found holds "
            "wherever it")
        say("has been looked for.")
    elif not k3:
        say("the sweep separation survives the right null. It is the "
            "strongest")
        say("claim outside the C(N) branch and the only thing in this "
            "repository")
        say("that a multiplicative ensemble has failed to reach. One "
            "seed and 512")
        say("draws; it needs another seed before it is leaned on, and "
            "the best")
        say("residual any draw achieved is printed above as the bound "
            "this run")
        say("actually gives.")
    else:
        say("the size separation holds at some N and not others, "
            "which neither")
        say("rem:decaydeep nor this run predicted. The per-N counts "
            "above are the")
        say("finding and the larger N are where to look.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
