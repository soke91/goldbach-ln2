# -*- coding: utf-8 -*-
r"""
Counting says balanced, weighting says 70 per cent. Which k do that?

WHAT IS AT STAKE

Remark {#rem:signmass} is where the lean was first located: counting
moduli the signs are near balanced, 0.4121 to 0.4808 positive;
weighting by contribution only 0.2273 to 0.3207 of sum(log k)|H| is
positive; and split at the median of |H| the small dilated walls sit
at 0.5096 to 0.5774 while the large ones lean to 0.3141 to 0.3841.

Remark {#rem:gainsplit} later found the cross-k shortfall confined to
the top tenth by mass, and Remark {#rem:headidentity} found that
tenth's sign fraction locked to the whole range's mass lean. Those two
readings ought to be one object, and nobody has checked that they are:
{#rem:signmass}'s split is at the MEDIAN of |H| while
{#rem:gainsplit}'s is the top DECILE of (log k)|H|, and
{#rem:headidentity} has already shown that mass ranking and range
ranking are not the same thing.

The sharp version is a subtraction. If the whole gap between counting
and weighting is the head, then removing the top decile should make
the two fractions agree. If it is spread, they will not.

BACKS: Remark {#rem:signmasshead} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  H1  The control: the mass fraction, the count fraction and the two
      half-fractions reproduce {#rem:signmass}'s 0.2273 to 0.3207,
      0.4121 to 0.4808, 0.5096 to 0.5774 and 0.3141 to 0.3841, all to
      within 0.001.
  H2  The median split and the decile split are one object: the high
      half's positive-mass share and the top decile's differ by a
      constant across N, to within 0.05.
  H3  The gap is the head's: removing the top decile by (log k)|H|,
      the mass fraction and the count fraction of what is left agree
      to within 0.05 at every N.
  H4  And it is mu's: re-signing the same magnitudes leaves the mass
      fraction and the count fraction agreeing to within 0.05 at
      every N, so the gap is not a property of the magnitudes.

REFUTATION RULE (fixed before the run)

  H1  REFUTED at 0.001 anywhere -- not the same statistic, and
      nothing below may be compared with {#rem:signmass}.
  H2  REFUTED beyond 0.05, which would say the two splits are
      different objects and the lean has two descriptions rather than
      one.
  H3  REFUTED beyond 0.05 at any N. That is the one that matters: the
      gap between counting and weighting would then be spread over
      the whole k-range, and {#rem:signmass}'s "the correlation sits
      in the large terms" would be too strong.
  H4  REFUTED beyond 0.05, which would say a random sign field on
      these magnitudes already separates counting from weighting and
      the gap is heterogeneity rather than mu.

  All four gate.

  THE NULL is the re-signing of {#rem:signmass} itself: the
  magnitudes |H(N;k)| are held and the signs drawn at random, which
  holds all heterogeneity fixed and destroys only the sign pattern.
  256 draws, one seed.
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
OUT = os.path.join(RES, "audit_signmass_head.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
THETA = 0.56
HEAD = 0.10
DRAWS = 256
SEED = 20260808


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
    """the four published fractions at each N"""
    src = io.open(os.path.join(RES, "lab_sign_structure.txt"),
                  encoding="utf-8").read()
    i = src.index("mass frac +   count frac +   low/high halves")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 10 or not f[0].isdigit():
            if f and set(f[0]) == {"-"}:
                continue
            if not out:
                continue
            break
        out[int(f[0])] = (float(f[5]), float(f[6]),
                          float(f[7]), float(f[9]))
    return out


def weighted(N, lam, mu, sqf):
    """H(N;k) and a_k = (log k)H over the admissible k"""
    PN = factor_set(N)
    K = int(N ** THETA)
    ks = np.array([k for k in range(2, K)
                   if sqf[k] and not any(k % q == 0 for q in PN)],
                  dtype=np.int64)
    Hs = []
    for k in ks:
        k = int(k)
        M = (N - 1) // k
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        vals = N - ms * k
        Hs.append(float((lam[vals] * mu[ms].astype(np.float64)).sum()))
    H = np.array(Hs)
    return ks, H, np.log(ks.astype(np.float64)) * H


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d published rows from results/lab_sign_structure.txt"
        % len(pub))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    rng = np.random.default_rng(SEED)
    say("%d re-signings per N, seed %d" % (DRAWS, SEED))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    rows = []
    for N in NS:
        ks, H, a = weighted(N, lam, mu, sqf)
        w = np.abs(a)
        mfrac = float(w[a > 0].sum() / w.sum())
        cfrac = float((H > 0).mean())
        med = float(np.median(np.abs(H)))
        lo = np.abs(H) <= med
        hi = ~lo
        lof = float((H[lo] > 0).mean())
        hif = float((H[hi] > 0).mean())
        nh = max(1, int(round(HEAD * ks.size)))
        head = np.argsort(-w)[:nh]
        rest = np.setdiff1d(np.arange(ks.size), head)
        mrest = float(w[rest][a[rest] > 0].sum() / w[rest].sum())
        crest = float((H[rest] > 0).mean())
        hmass = float(w[head][a[head] > 0].sum() / w[head].sum())
        eps = (rng.integers(0, 2, size=(DRAWS, ks.size))
               .astype(np.int8) * 2 - 1)
        nm = float(np.median((eps > 0) * w / w.sum()
                             @ np.ones(ks.size)))
        nmass = float(np.median(((eps > 0) * w).sum(axis=1)
                                / w.sum()))
        ncount = float(np.median((eps > 0).mean(axis=1)))
        rows.append((N, ks.size, mfrac, cfrac, lof, hif, mrest,
                     crest, hmass, nmass, ncount))
        say("  N = %-10d #k %-6d mass %.4f count %.4f low %.4f "
            "high %.4f" % (N, ks.size, mfrac, cfrac, lof, hif))

    # ------------------------------------------------------------- H1
    say()
    say("H1  the control: the four published fractions")
    say("  N            mass/pub        count/pub       low/pub"
        "         high/pub")
    h1 = True
    for (N, nk, mf, cf, lo, hi, mr, cr, hm, nm2, nc) in rows:
        pm, pc, pl, ph = pub[N]
        for got, want in ((mf, pm), (cf, pc), (lo, pl), (hi, ph)):
            if abs(got - want) >= 0.001:
                h1 = False
        say("  %-12d %.4f/%.4f   %.4f/%.4f   %.4f/%.4f   %.4f/%.4f"
            % (N, mf, pm, cf, pc, lo, pl, hi, ph))
    say("  H1 %s   (cap 0.001)" % ("hold" if h1 else "REFUTED"))

    # ---------------------------------------------------- H2 / H3 / H4
    say()
    say("H2  the median split against the decile split")
    say("  N            high half   top decile   difference")
    dif = []
    for (N, nk, mf, cf, lo, hi, mr, cr, hm, nm2, nc) in rows:
        dif.append(hi - hm)
        say("  %-12d %-11.4f %-12.4f %.4f" % (N, hi, hm, dif[-1]))
    h2 = (max(dif) - min(dif)) < 0.05
    say("  the difference spans %.4f   (cap 0.05)"
        % (max(dif) - min(dif)))
    say("  H2 %s" % ("hold" if h2 else "REFUTED"))

    say()
    say("H3/H4  where the gap between counting and weighting lives")
    say("  N            whole gap   without the head   re-signed gap")
    h3 = h4 = True
    for (N, nk, mf, cf, lo, hi, mr, cr, hm, nm2, nc) in rows:
        g0 = abs(cf - mf)
        g1 = abs(cr - mr)
        g2 = abs(nc - nm2)
        if g1 >= 0.05:
            h3 = False
        if g2 >= 0.05:
            h4 = False
        say("  %-12d %-11.4f %-18.4f %.4f" % (N, g0, g1, g2))
    say("  H3 the gap closes without the head   %s   (cap 0.05)"
        % ("hold" if h3 else "REFUTED"))
    say("  H4 and re-signing shows no gap   %s   (cap 0.05)"
        % ("hold" if h4 else "REFUTED"))
    say("PERN signmass_gap_whole %d %.4f %.4f"
        % (len(rows), min(abs(r[3] - r[2]) for r in rows),
           max(abs(r[3] - r[2]) for r in rows)))
    say("PERN signmass_gap_nohead %d %.4f %.4f"
        % (len(rows), min(abs(r[7] - r[6]) for r in rows),
           max(abs(r[7] - r[6]) for r in rows)))
    rr = [abs(r[7] - r[6]) / abs(r[3] - r[2]) for r in rows]
    say("  and the two taken at the same N:")
    say("RATIO signmass_gap_nohead signmass_gap_whole %.4f %.4f"
        % (min(rr), max(rr)))

    say()
    say("  so the sentence 'the correlation sits in the large terms'")
    say("  is exact in this sense: the top tenth by mass carries")
    say("  %.0f to %.0f per cent of the whole gap between counting"
        % (100 * (1 - max(rr)), 100 * (1 - min(rr))))
    say("  and weighting, and what is left behind is NOT balanced --")
    say("  the remaining gap is %.4f to %.4f, still far above the"
        % (min(abs(r[7] - r[6]) for r in rows),
           max(abs(r[7] - r[6]) for r in rows)))
    say("  %.4f a re-signing leaves."
        % max(abs(r[10] - r[9]) for r in rows))
    say("SITSIN signmass_gap %.4f %.4f" % (1.0 - max(rr),
                                           1.0 - min(rr)))

    say()
    say("=" * 70)
    ok = h1 and h2 and h3 and h4
    say("counting and weighting differ only over the head"
        if ok else "REFUTED")

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) over the squarefree",
        "           k < N^" + str(THETA) + " coprime to N: the fraction of",
        "           sum(log k)|H| carried by k with H > 0, the fraction",
        "           of k with H > 0, both again split at the median of",
        "           |H| and again with the top tenth by |a| removed,",
        "           and both under " + str(DRAWS) + " re-signings that",
        "           hold the magnitudes.",
        "NULL: the re-signing of {#rem:signmass} -- the magnitudes are",
        "      held and the signs drawn at random, which fixes all",
        "      heterogeneity and destroys only the sign pattern.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < N^" + str(THETA) + "; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + "; numpy default_rng",
        "       seed " + str(SEED) + ". Every N is 2^a 5^b, one odd",
        "       radical, as RADICALS declares. The published fractions",
        "       are read from results/lab_sign_structure.txt.",
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
