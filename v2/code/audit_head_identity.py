# -*- coding: utf-8 -*-
r"""
Is the head a new object, or the small-k end of the old one?

WHAT IS AT STAKE

Remark {#rem:gainsplit} localised the cross-k shortfall: the bottom
nine tenths of the k cancel better than square root, the top tenth
does not cancel at all, and the fraction of that top tenth carrying a
single sign falls 1.0000 to 0.8274 over a factor 128. OPEN item 5 now
turns on whether that fraction goes to 1/2 -- which would release
e(G) -- or stops above it, and eight points cannot separate the two.

Before fitting a shape to it, it is worth asking whether the head is
a new object at all. |a_k| = (log k)|H(N;k)| is largest where the
inner sum is longest, so the top decile by mass should be very nearly
the smallest decile of k, and the sign it carries should be the lean
Remark {#rem:signmass} measured on the whole range. If both hold, the
head's fraction is not a free quantity: it is the lean read on a
sub-range, and Remarks {#rem:sievedepth} and {#rem:levelweighted}
have already established what carries the lean -- primality, at a
level no bounded modulus reaches.

The published tables already show the head's fraction and the whole
range's mass lean moving in step, so Z4 is registered in that
direction rather than against it. What they do not show, and what is
tested here, is whether the head is the small-k end and whether its
sign is the lean's.

BACKS: Remark {#rem:headidentity} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  The control: the head's gain reproduces {#rem:gainsplit}'s
      1.0000 to 1.4639 and its single-sign fraction its 1.0000 to
      0.8274, both to within 0.001.
  Z2  The head is the small-k end: the overlap between the top decile
      by |a| and the smallest decile of k is above 0.8 at every N.
  Z3  And it carries the lean's sign: the head's majority sign is
      negative at every N.
  Z4  And its fraction is the lean, not a new quantity: the head's
      single-sign fraction minus the whole range's mass-weighted
      negative share is constant across N to within 0.05.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED at 0.001 -- not the same statistic, and nothing below
      may be compared with {#rem:gainsplit}.
  Z2  REFUTED below 0.8 at any N. The head would then be selected by
      something other than the length of the inner sum, and would be
      a genuinely new set to explain.
  Z3  REFUTED if the majority sign is ever positive, which would sever
      the head from the lean {#rem:signmass} measured.
  Z4  REFUTED if the difference varies by more than 0.05. That is the
      one that matters: the head's fraction would then be a free
      quantity with its own limit, and OPEN item 5 would need a shape
      fitted to it. If it holds, item 5 is the lean question again
      and has no independent degree of freedom.

  All four gate.

  NO NULL IS RUN and none applies to Z1-Z4: a deterministic set is
  compared with another deterministic set and a sign is counted. The
  coin arms for the lean were run in audit_lean_floor.py and for the
  gain in audit_crossk_reference.py.
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
OUT = os.path.join(RES, "audit_head_identity.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000,
      12_800_000, 25_600_000]
THETA = 0.56
HEAD = 0.10


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
    """the head's gain and single-sign fraction"""
    src = io.open(os.path.join(RES, "audit_gain_split.txt"),
                  encoding="utf-8").read()
    i = src.index("  head tenth   ")
    gains = [float(v) for v in src[i:].splitlines()[0].split()[2:]]
    j = src.index("N            top-decile share  same sign in the "
                  "head")
    agr = {}
    for ln in src[j:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        agr[int(f[0])] = float(f[2])
    return gains, agr


def weighted(N, lam, mu, sqf):
    """(log k)H(N;k) over the squarefree k < N^theta coprime to N"""
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
    return ks, np.log(ks.astype(np.float64)) * np.array(Hs)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubg, pubagr = read_published()
    say("read %d head gains and %d single-sign fractions from "
        "results/audit_gain_split.txt" % (len(pubg), len(pubagr)))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    rows = []
    for N in NS:
        ks, a = weighted(N, lam, mu, sqf)
        w = np.abs(a)
        nh = max(1, int(round(HEAD * ks.size)))
        hd = np.argsort(-w)[:nh]
        small = np.arange(nh)          # ks is increasing in k
        ov = float(len(set(hd.tolist()) & set(small.tolist())) / nh)
        sh = np.sign(a[hd])
        neg = float((sh < 0).mean())
        frac = max(neg, 1.0 - neg)
        s = abs(float(a[hd].sum()))
        g = float(w[hd].sum()) / s if s > 0 else float("inf")
        negmass = float(w[a < 0].sum() / w.sum())
        hH = np.argsort(-np.abs(a / np.log(ks.astype(np.float64))))[:nh]
        ovH = float(len(set(hd.tolist()) & set(hH.tolist())) / nh)
        medrank = float(np.median(hd) / ks.size)
        rows.append((N, ks.size, nh, g, frac, ov, neg > 0.5, negmass,
                     ovH, medrank))
        say("  N = %-10d #k %-6d head %-5d gain %-8.4f frac %-8.4f "
            "overlap %-7.4f neg mass %.4f"
            % (N, ks.size, nh, g, frac, ov, negmass))

    # ------------------------------------------------------------- Z1
    say()
    say("Z1  the control: the head's gain and single-sign fraction")
    say("  N            gain here  published  frac here  published")
    z1 = True
    for (N, nk, nh, g, fr, ov, isneg, nm, oh, mr), pg in zip(
            rows, pubg):
        pa = pubagr[N]
        if abs(g - pg) >= 0.001 or abs(fr - pa) >= 0.001:
            z1 = False
        say("  %-12d %-10.4f %-10.4f %-10.4f %.4f" % (N, g, pg, fr, pa))
    say("  Z1 %s   (cap 0.001)" % ("hold" if z1 else "REFUTED"))

    # -------------------------------------------------- Z2 / Z3 / Z4
    say()
    say("Z2/Z3  is the head the small-k end, and is its sign the "
        "lean's?")
    say("  N            overlap with the smallest decile  majority")
    z2 = z3 = True
    for N, nk, nh, g, fr, ov, isneg, nm, oh, mr in rows:
        if ov <= 0.8:
            z2 = False
        if not isneg:
            z3 = False
        say("  %-12d %-33.4f %s"
            % (N, ov, "negative" if isneg else "POSITIVE"))
    say("  Z2 the overlap is above 0.8 at every N   %s   (floor 0.8)"
        % ("hold" if z2 else "REFUTED"))
    say("  Z3 the majority sign is negative at every N   %s"
        % ("hold" if z3 else "REFUTED"))
    say("PERN headidentity_overlap %d %.4f %.4f"
        % (len(rows), min(r[5] for r in rows),
           max(r[5] for r in rows)))
    say()
    say("  DIAGNOSTIC on Z2 (post hoc). The head is not the small-k")
    say("  end, so what is it? |a_k| = (log k)|H(N;k)| and the log k")
    say("  factor grows with k while |H| falls, so the two pull")
    say("  against each other. Dropping the weight:")
    say("  N            overlap with the top decile by |H|  "
        "median rank of the head")
    for N, nk, nh, g, fr, ov, isneg, nm, oh, mr in rows:
        say("  %-12d %-35.4f %.4f" % (N, oh, mr))
    say("SPLITOVERLAP crossk %.4f %.4f"
        % (min(r[5] for r in rows), max(r[5] for r in rows)))
    say("  the head sits at the %.2f to %.2f point of the k-order, so"
        % (min(r[9] for r in rows), max(r[9] for r in rows)))
    say("  it is spread across the range and not a restriction of it.")

    say()
    say("Z4  is the head's fraction the lean read on a sub-range?")
    say("  N            head fraction  negative mass  difference")
    dif = []
    for N, nk, nh, g, fr, ov, isneg, nm, oh, mr in rows:
        dif.append(fr - nm)
        say("  %-12d %-14.4f %-14.4f %.4f" % (N, fr, nm, dif[-1]))
    sp = max(dif) - min(dif)
    z4 = sp < 0.05
    say("  the difference runs %.4f to %.4f, a spread of %.4f"
        % (min(dif), max(dif), sp))
    say("  (cap 0.05)")
    say("  and the head's fraction is a proportion over %d to %d"
        % (min(r[2] for r in rows), max(r[2] for r in rows)))
    say("  terms, so it carries a binomial error of its own:")
    say("  N            head size  s.e. of the fraction")
    ses = []
    for N, nk, nh, g, fr, ov, isneg, nm, oh, mr in rows:
        se = math.sqrt(max(fr * (1.0 - fr), 0.0) / nh)
        ses.append(se)
        say("  %-12d %-10d %.6f" % (N, nh, se))
    rng = np.random.default_rng(20260808)
    z = rng.normal(0.0, 1.0, size=(20000, len(rows)))
    span = float((z.max(axis=1) - z.min(axis=1)).mean())
    samp = span * float(np.mean(ses))
    say("  the expected span of %d such draws is %.4f times the"
        % (len(rows), span))
    say("  standard error, i.e. %.4f against the observed %.4f"
        % (samp, sp))
    say("CONSTSPREAD headidentity_diff %.6f %.6f" % (sp, samp))
    if sp > 2.0 * samp:
        say("CONST DRIFTS headidentity_diff")
    say("PERN headidentity_headfrac %d %.4f %.4f"
        % (len(rows), min(r[4] for r in rows),
           max(r[4] for r in rows)))
    say("PERN headidentity_negmass %d %.4f %.4f"
        % (len(rows), min(r[7] for r in rows),
           max(r[7] for r in rows)))
    rr = [r[4] / r[7] for r in rows]
    say("RATIO headidentity_headfrac headidentity_negmass %.4f %.4f"
        % (min(rr), max(rr)))
    say("  Z4 %s" % ("hold" if z4 else "REFUTED"))

    say()
    say("  what that settles. If the head is the small-k end and its")
    say("  sign is the lean's, then its single-sign fraction is not a")
    say("  free quantity with a limit of its own -- it is the lean of")
    say("  {#rem:signmass} read over the longest inner sums, and")
    say("  {#rem:sievedepth} has already found what carries that: the")
    say("  primality of N - mk, at a level no bounded modulus reaches.")

    say()
    say("=" * 70)
    ok = z1 and z2 and z3 and z4
    say("the head is the lean, read on the longest inner sums"
        if ok else "REFUTED")

    head = [
        "STATISTIC: for a_k = (log k)H(N;k) over the squarefree",
        "           k < N^" + str(THETA) + " coprime to N: the top tenth",
        "           by |a|, its cross-k gain, the fraction of it",
        "           carrying a single sign, its overlap with the",
        "           smallest tenth of k, and the whole range's",
        "           mass-weighted share of negative terms.",
        "NULL: none applies -- a deterministic set is compared with",
        "      another and a sign is counted. The coin arms were run",
        "      in audit_lean_floor.py for the lean and in",
        "      audit_crossk_reference.py for the gain.",
        "FIELD: N = 2e5 through 2.56e7 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < N^" + str(THETA) + "; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + ". Every N is 2^a 5^b,",
        "       one odd radical, as RADICALS declares. The published",
        "       head gains and fractions are read from",
        "       results/audit_gain_split.txt.",
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
