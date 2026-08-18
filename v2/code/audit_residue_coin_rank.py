# -*- coding: utf-8 -*-
r"""
Is mu just one more coin against delta, or the most extreme of them?

WHAT IS AT STAKE

Remark {#rem:residuecancel} is the statement the route turns on: the
residue buys exactly a coin's worth of cancellation against
delta = Lambda(N-mk) - beta w(m,k), so "nothing in R is doing better
than random signs, so no argument that mu is special against delta can
help". The evidence is a band: |R|/coin runs over [0.78, 1.33] across
thirty-four octaves and "sits on 1".

A band is not a test. The published ratio divides mu's octave mean by
the mean of sixteen coins, and the mean of sixteen draws has its own
scatter, which was never computed. Nor was the sign counted: reading
down that column, most entries are below 1 and at the largest N all
six are, which is what a small systematic advantage would look like
and also what noise would look like. The two are separated by treating
mu as a seventeenth draw.

That is an exact test and needs no distributional assumption. The
sixteen sign vectors are global -- one vector per draw, held across
every k, exactly as mu is -- so under the hypothesis that mu is just
another sign pattern, the seventeen are exchangeable and mu's rank
among them is uniform.

The implementation is the bitmask one of audit_level_slope_reach.py,
independent of lab_residue_cancellation.py's residue-class masking,
so V1 is a cross-check of the published table.

BACKS: Remark {#rem:coinrank} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The control: the l1 exponents reproduce the published
      0.9745, 0.9919, 0.9962, 0.9963, 0.9979 to within 0.001, and
      every published |R|/l2 octave mean to within 0.01.
  V2  The published band holds: every |R|/coin octave ratio is inside
      [0.78, 1.33].
  V3  mu is not distinguishable from a coin: its pooled statistic --
      the mean over octaves of log of its octave mean against the
      mean of the OTHER draws -- lies inside the range the sixteen
      coins give for themselves, at least one above and one below.
  V4  And the sign is not systematic: over the octaves, the number at
      which mu's octave mean falls below the coins' median is within
      two standard deviations of half.

REFUTATION RULE (fixed before the run)

  V1  REFUTED at either cap -- not the same statistic, and nothing
      below may be compared with {#rem:residuecancel}.
  V2  REFUTED if any ratio leaves the band, which would mean the
      published band is not what this measurement gives.
  V3  REFUTED if mu's pooled statistic is the most extreme of the
      seventeen. Under exchangeability that is a 2-in-17 event, and
      it would say mu is NOT just another sign pattern against delta
      -- the sentence "no argument that mu is special against delta
      can help" would have to be qualified, and in the direction that
      matters, since a systematic advantage is what such an argument
      would need.
  V4  REFUTED if the count is outside two standard deviations of
      half, which is the same conclusion by a route that uses no
      magnitudes at all.

  All four gate.

  THE NULL IS THE POINT HERE and it is run: sixteen global sign
  vectors on the odd squarefree m, each held across all k as mu is,
  summed against the IDENTICAL delta. Same convention and same seed as
  lab_residue_cancellation.py, so the two are comparable.
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
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_residue_coin_rank.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
OCT = [2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288,
       2097152]
MINPTS = 10
COINS = 16
SEED = 20260808
CLIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
    """von Mangoldt and Moebius, the cofactor kept in int32"""
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
    cof = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        cof[p::p] //= p
        pk = p * p
        while pk <= n:
            cof[pk::pk] //= p
            if pk > n // p:
                break
            pk *= p
    big = cof > 1
    del cof
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def residue_mask(n, qs):
    """bit i of mask[v] is set exactly when qs[i] divides v"""
    m = np.zeros(n + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        m[0::q] |= np.uint16(1 << i)
    return m


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


def read_published():
    """the l1 exponents and the |R|/l2 octave means"""
    src = io.open(os.path.join(RES, "lab_residue_cancellation.txt"),
                  encoding="utf-8").read()
    i = src.index("N            l1 exponent   correlation   "
                  "thinnest bin")
    l1 = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            if f and f[0].startswith("POP"):
                continue
            break
        l1[int(f[0])] = float(f[1])
    j = src.index("V3  |R| against the l2 norm of its own summands")
    row = re.compile(r"^\s*\[\s*(\d+)\s*,\s*(\d+)\s*\)\s+\d+"
                     r"\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$")
    tab, cur = {}, None
    for ln in src[j:].splitlines():
        f = ln.split()
        if len(f) == 3 and f[0] == "N" and f[1] == "=":
            cur = int(f[2])
            continue
        m = row.match(ln)
        if cur is not None and m:
            tab[(cur, int(m.group(1)))] = (float(m.group(3)),
                                           float(m.group(5)))
        elif cur is not None and ln.strip().startswith("the scale"):
            break
    return l1, tab


def one_N(N, lam, mu, sqf, vmask, qs, coin, artin, twin):
    """per admissible k: |R|, l1, l2 and the sixteen coin sums"""
    PN = factor_set(N)
    ks, Rv, L1, L2, CS, Hs, Ps = [], [], [], [], [], [], []
    for k in range(2, KCAP):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[sqf[ms]]
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        for q in factor_set(k):
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        g = mu[ms].astype(np.float64)
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        ks.append(k)
        Hs.append(float((lam[vals] * g).sum()))
        Ps.append(ck * float(g[keep].sum()))
        CS.append((ms, vals, g, keep, ck))
    ks = np.array(ks, dtype=np.int64)
    H = np.array(Hs)
    P = np.array(Ps)
    beta = float((H * P).sum() / (P * P).sum())

    Rv, L1, L2, Cv = [], [], [], []
    for ms, vals, g, keep, ck in CS:
        w = np.where(keep, ck, 0.0)
        d = lam[vals] - beta * w
        Rv.append(abs(float((g * d).sum())))
        L1.append(float(np.abs(d).sum()))
        L2.append(float(np.sqrt((d * d).sum())))
        Cv.append(np.abs(coin[:, ms].astype(np.float64) @ d))
    return (ks, np.array(Rv), np.array(L1), np.array(L2),
            np.array(Cv).T, beta)


def octaves(ks, N, L2):
    """the closed octaves of N/k; a k whose delta vanishes carries
    no ratio and is dropped, as the published script does"""
    inner = (N // ks).astype(np.float64)
    out = []
    for a, b in zip(OCT[:-1], OCT[1:]):
        sel = (inner >= a) & (inner < b) & (L2 > 0)
        if int(sel.sum()) >= MINPTS:
            out.append((a, b, sel))
    return out


def octfit(ks, N, v):
    """the published octave fit: mean of N/k against mean of v"""
    inner = (N // ks).astype(np.float64)
    cx, cy = [], []
    for a, b in zip(OCT[:-1], OCT[1:]):
        sel = (inner >= a) & (inner < b)
        if int(sel.sum()) < MINPTS:
            continue
        m = float(v[sel].mean())
        if not (m > 0):
            continue
        cx.append(math.log(float(inner[sel].mean())))
        cy.append(math.log(m))
    x = np.array(cx)
    y = np.array(cy)
    return float(np.polyfit(x, y, 1)[0]), x.size


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    publ1, pubtab = read_published()
    say("read %d published l1 exponents and %d octave rows from "
        "results/lab_residue_cancellation.txt"
        % (len(publ1), len(pubtab)))

    NMAX = max(NS)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    vmask = residue_mask(NMAX, qs)

    rng = np.random.default_rng(SEED)
    coin = np.zeros((COINS, NMAX + 1), dtype=np.int8)
    live = np.flatnonzero(sqf)
    live = live[live % 2 == 1]
    coin[:, live] = (rng.integers(0, 2, size=(COINS, live.size))
                     .astype(np.int8) * 2 - 1)
    say("%d global sign vectors on the %d odd squarefree m, seed %d"
        % (COINS, live.size, SEED))

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    data = []
    for N in NS:
        ks, Rv, L1, L2, Cv, beta = one_N(
            N, lam, mu, sqf, vmask, qs, coin, artin, twin)
        data.append((N, ks, Rv, L1, L2, Cv, beta))
        say("  N = %-10d #k = %-7d beta = %.6f" % (N, ks.size, beta))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(d[0]) if q > 2))
                  for d in data)))

    # ------------------------------------------------------------- V1
    say()
    say("V1  the control: the l1 exponent and the octave means")
    say("  N            l1 exponent  published    worst |R|/l2 diff")
    v1 = True
    for N, ks, Rv, L1, L2, Cv, beta in data:
        a, nb = octfit(ks, N, L1)
        if not (abs(a - publ1[N]) < 0.001):
            v1 = False
        worst = 0.0
        for lo, hi, sel in octaves(ks, N, L2):
            m = float((Rv[sel] / L2[sel]).mean())
            key = (N, lo)
            if key in pubtab:
                worst = max(worst, abs(m - pubtab[key][0]))
        if worst >= 0.01:
            v1 = False
        say("  %-12d %-12.4f %-12.4f %.5f" % (N, a, publ1[N], worst))
    say("  V1 %s   (cap 0.001 and cap 0.01)"
        % ("hold" if v1 else "REFUTED"))

    # ---------------------------------------------------------- V2/V3
    say()
    say("V2/V3  mu as a seventeenth draw")
    say("  for each octave: mu's mean of |R|/l2, the mean of the")
    say("  sixteen coins, and mu's z against the coin-to-coin scatter")
    say("  N          octave       #k     mu       coins    ratio   z")
    ratios, zs, below = [], [], 0
    # per-draw octave means, mu first then the sixteen coins
    pooled = [[] for _ in range(COINS + 1)]
    means = []
    v2 = True
    for N, ks, Rv, L1, L2, Cv, beta in data:
        for lo, hi, sel in octaves(ks, N, L2):
            a = float((Rv[sel] / L2[sel]).mean())
            cm = (Cv[:, sel] / L2[sel]).mean(axis=1)
            b = float(cm.mean())
            s = float(cm.std(ddof=1))
            r = a / b
            ratios.append(r)
            zs.append((a - b) / s)
            if a < float(np.median(cm)):
                below += 1
            if not (0.78 <= r <= 1.33):
                v2 = False
            allm = np.concatenate(([a], cm))
            means.append(allm)
            for j in range(COINS + 1):
                others = np.delete(allm, j)
                pooled[j].append(math.log(allm[j] / float(others.mean())))
            say("  %-10d [%-6d,%-7d) %-6d %-8.4f %-8.4f %-7.4f %+.2f"
                % (N, lo, hi, int(sel.sum()), a, b, r, (a - b) / s))
    say("  V2 every ratio inside the band   %s"
        % ("hold" if v2 else "REFUTED"))
    say("  the band this measurement gives: [%.4f, %.4f]"
        % (min(ratios), max(ratios)))
    say("PERN residue_coin_ratio %d %.4f %.4f"
        % (len(ratios), min(ratios), max(ratios)))

    T = [float(np.mean(p)) for p in pooled]
    Tmu, Tc = T[0], T[1:]
    nabove = sum(1 for t in Tc if t > Tmu)
    nbelow = sum(1 for t in Tc if t < Tmu)
    v3 = nabove > 0 and nbelow > 0
    say()
    say("  the pooled statistic, each draw against the mean of the")
    say("  other sixteen, averaged over the %d octaves:" % len(ratios))
    say("  mu      %+.6f" % Tmu)
    say("  coins   " + ", ".join("%+.4f" % t for t in sorted(Tc)))
    say("  %d coins above mu, %d below   rank %d of %d"
        % (nabove, nbelow, nbelow + 1, COINS + 1))
    say("  V3 mu is inside the coins' own range   %s"
        % ("hold" if v3 else "REFUTED"))

    # ------------------------------------------------------------- V4
    say()
    say("V4  the sign, counted, with no magnitudes used")
    n = len(ratios)
    exp = n / 2.0
    sd = math.sqrt(n) / 2.0
    v4 = abs(below - exp) <= 2.0 * sd
    say("  mu below the coins' median at %d of %d octaves" % (below, n))
    say("  expected %.1f, two standard deviations %.1f" % (exp, 2 * sd))
    say("  V4 %s" % ("hold" if v4 else "REFUTED"))
    say("MARGIN audit_residue_coin_rank %.4f %.4f"
        % (abs(below - exp), 2.0 * sd))
    say("SCALES audit_residue_coin_rank 1")
    say("ONE SCALE audit_residue_coin_rank")

    say()
    say("  DIAGNOSTIC (post hoc). V4's null assumes the %d octaves"
        % n)
    say("  are independent, and they are not: the same sixteen sign")
    say("  vectors are used at every octave and at every N, and the")
    say("  octaves at different N run over overlapping m. The")
    say("  exchangeable version costs nothing -- count the same thing")
    say("  for each coin, each against the median of the other")
    say("  sixteen draws:")
    M = np.array(means)
    counts = []
    for j in range(COINS + 1):
        others = np.delete(M, j, axis=1)
        counts.append(int((M[:, j] < np.median(others, axis=1)).sum()))
    say("  mu      %d of %d" % (counts[0], n))
    say("  coins   " + ", ".join(str(c) for c in sorted(counts[1:])))
    ge = sum(1 for c in counts[1:] if c >= counts[0])
    say("  %d of the %d coins lean at least as much as mu, so the"
        % (ge, COINS))
    say("  count is a %d-in-%d event among exchangeable draws, not the"
        % (ge + 1, COINS + 1))
    say("  %.1f-sigma the binomial gives. The mean z over the octaves"
        % (abs(below - exp) / sd))
    say("  is %+.4f; the pooled test above is that comparison done"
        % float(np.mean(zs)))
    say("  exactly, and it puts mu at rank %d of %d."
        % (nbelow + 1, COINS + 1))
    say("EXCHANGE audit_residue_coin_rank %d %d %d"
        % (counts[0], ge + 1, COINS + 1))
    say("  So V4 fails as registered and the failure is the null, not")
    say("  the data: with the dependence taken out, mu's lean is")
    say("  %s the coins' own." % ("not separated from"
                                  if ge > 0 else "outside"))

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("mu is one more sign pattern against delta"
        if ok else "REFUTED")

    head = [
        "STATISTIC: for each admissible k, R = sum_m mu(m) delta(m,k)",
        "           with delta = Lambda(N-mk) - beta w(m,k), the l1 and",
        "           l2 norms of the same summands, and the sums of the",
        "           same delta against sixteen global sign vectors;",
        "           then, per closed octave of N/k, mu's mean of",
        "           |R|/l2 against the sixteen coins' means, mu's z",
        "           against the coin-to-coin scatter, a pooled",
        "           leave-one-out statistic ranking mu among the",
        "           seventeen draws, and the count of octaves at which",
        "           mu falls below the coins' median.",
        "NULL: sixteen global sign vectors on the odd squarefree m,",
        "      each held across all k as mu is, summed against the",
        "      IDENTICAL delta; same convention and seed as",
        "      lab_residue_cancellation.py. Under the hypothesis that",
        "      mu is one more sign pattern the seventeen are",
        "      exchangeable, so the rank test is exact.",
        "FIELD: N = 2e5 to 3.2e6 by doubling; k squarefree and coprime",
        "       to N with 2 <= k < " + str(KCAP) + "; m odd, squarefree",
        "       and coprime to k, m <= (N-1)/k; the sieve weight over",
        "       the odd primes below " + str(QSIEVE) + "; beta refitted",
        "       as sum(H P)/sum(P^2) on the same k-range; octaves",
        "       closed at both ends and used only when they hold at",
        "       least " + str(MINPTS) + " k; numpy default_rng seed "
        + str(SEED) + ".",
        "       Every N is 2^a 5^b, one odd radical, as RADICALS says.",
        "       The published exponents and octave means are read from",
        "       results/lab_residue_cancellation.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not ok:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
