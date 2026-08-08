# -*- coding: utf-8 -*-
r"""
Does the sieve-weighted predictor reach the lean's own k?

WHAT IS AT STAKE

Remark {#rem:oddmertensrange} withdrew Modd's claim to explain the
lean: on the k-range the lean is measured on, Modd predicts a lean
nearly twice too deep and predicts no decay at all. What that leaves
standing is the other predictor. Remark {#rem:survivors} reports that
the sieve-weighted P(N;k) = sum_m mu(m) w(m,k) agrees with sign H at
0.9274 down to 0.9080 -- "nine tenths of the sign of the dilated wall
is a sieve-weighted Moebius sum over its own inner range".

It is measured on the same window as Modd's: 2 <= N/k <= 1000, short
inner sums. Remark {#rem:predictable} says so in as many words. The
lean f runs over k < N^0.56, whose inner lengths reach 2.1e6. So the
same question has to be put to P, and it is now the only elementary
mechanism left for the lean where the lean is measured.

The control here uses every admissible k in the published window
rather than the 60000-per-N subsample, so the cap on S1 is 0.005
rather than the 0.001 used elsewhere: a subsample of about 15000 k
carries a binomial standard error near 0.0023 on its own.

The implementation is independent of lab_survivor_selection.py's.

BACKS: Remark {#rem:survivorrange} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  The control: on 2 <= N/k <= 1000 the agreement of sign H with
      sign P reproduces the published 0.9274, 0.9224, 0.9204, 0.9154,
      0.9080 and with sign Modd the published 0.8377 ... 0.8238, both
      to within 0.005.
  S2  P survives the longer window where Modd did not: on k < N^0.56
      the agreement stays above 0.85 at every N.
  S3  And it carries the lean: replacing each sign H by sign P on
      mu's own magnitudes gives a predicted |0.5 - f| within a factor
      1.5 of mu's at every N.
  S4  And its trend: the predicted lean's slope against log N is
      within two standard errors of mu's.

REFUTATION RULE (fixed before the run)

  S1  REFUTED at 0.005 anywhere -- not the same statistic, and
      nothing below may be compared with {#rem:survivors}.
  S2  REFUTED below 0.85 at any N. P would then be a statement about
      short inner sums only, like Modd, and "nine tenths of the sign
      of the dilated wall" would need the window attached.
  S3  REFUTED outside a factor 1.5 at any N. That is the one that
      matters: with Modd already withdrawn there, it would leave NO
      elementary predictor carrying the lean where the lean is
      measured, and {#rem:survivors}'s closing claim -- that the lean
      is a sieve-weighted Moebius sum over a short range -- would be
      a statement about short ranges and not about the lean.
  S4  REFUTED beyond two standard errors, the same conclusion for the
      trend rather than the level.

  All four gate.

  THE NULL for S1 is the published one: sixteen permutations of the
  predictor's signs across the sampled k. For S3 and S4 the floor is
  the 256 global sign vectors of audit_lean_floor.py on the identical
  magnitudes, so the predicted lean is read on the same scale as mu's
  and as Modd's in {#rem:oddmertensrange}.
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
OUT = os.path.join(RES, "audit_survivor_range.txt")

CTRL = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
XHI = 1000
QSIEVE = 30
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
    """the published P and Modd agreements"""
    src = io.open(os.path.join(RES, "lab_survivor_selection.txt"),
                  encoding="utf-8").read()
    i = src.index("N            sign P    sign Modd   draws max   #k")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            if not out:
                continue
            break
        out[int(f[0])] = (float(f[1]), float(f[2]))
    return out


def scan(N, klo, khi, lam, mu, sqf, oddsqf, vmask, qs):
    """H and the sieve-weighted P over a k-window, in the convention
    of lab_survivor_selection.py: m odd, squarefree, coprime to k and
    at most N//k, with H taken on that same set"""
    PN = factor_set(N)
    ks, Hs, Ps = [], [], []
    for k in range(max(2, klo), khi):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = N // k
        if M < 1:
            continue
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[oddsqf[ms]]
        for q in factor_set(k):
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        Hs.append(float((lam[vals] * mu[ms].astype(np.float64)).sum()))
        kb = 0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        # the sieve weight is C_k on the survivors, so its sign is
        # the sign of the unweighted Moebius sum over them
        Ps.append(float(mu[ms[keep]].sum(dtype=np.int64)))
        ks.append(k)
    return (np.array(ks, dtype=np.int64), np.array(Hs),
            np.array(Ps))


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se, abs(a) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    wsrc = io.open(os.path.join(RES, "audit_oddmertens_range.txt"),
                   encoding="utf-8").read()
    say("read %d published agreements from "
        "results/lab_survivor_selection.txt" % len(pub))

    NMAX = max(NS)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    oddsqf = sqf.copy()
    oddsqf[::2] = False
    vmask = residue_mask(NMAX, qs)
    om = np.zeros(NMAX + 1, dtype=np.int64)
    om[1::2] = mu[1::2]
    modd = np.cumsum(om)
    del om

    rng = np.random.default_rng(SEED)
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    # ------------------------------------------------------------- S1
    say()
    say("S1  the control: agreement over 2 <= N/k <= %d, every k"
        % XHI)
    say("  N            #k       sign P   published  sign Modd  "
        "published")
    s1 = True
    for N in CTRL:
        ks, H, P = scan(N, N // XHI, N // 2, lam, mu, sqf, oddsqf,
                        vmask, qs)
        inner = N // ks
        # the published filter: a k is compared only where all three
        # signs are defined
        sel = ((inner >= 2) & (inner <= XHI) & (H != 0) & (P != 0)
               & (modd[inner] != 0))
        sh = np.sign(H[sel])
        sp = np.sign(P[sel])
        sm = np.sign(modd[inner[sel]])
        ap = float((sh == sp).mean())
        am = float((sh == sm).mean())
        pp, pm = pub[N]
        if abs(ap - pp) >= 0.005 or abs(am - pm) >= 0.005:
            s1 = False
        say("  %-12d %-8d %-8.4f %-10.4f %-10.4f %.4f"
            % (N, int(sel.sum()), ap, pp, am, pm))
    say("  S1 %s   (cap 0.005; the published figures are a "
        "60000-k subsample of this set)" % ("hold" if s1 else "REFUTED"))

    # -------------------------------------------------- S2 / S3 / S4
    say()
    say("S2/S3/S4  on the k-range the lean is measured on, "
        "k < N^%.2f" % THETA)
    say("  N            #k     agreement  mu lean   P lean    ratio")
    s2 = s3 = True
    xs, mus, prd, flo, agrs = [], [], [], [], []
    for N in NS:
        ks, H, P = scan(N, 2, int(N ** THETA), lam, mu, sqf, oddsqf,
                        vmask, qs)
        a = np.log(ks.astype(np.float64)) * H
        l1 = float(np.abs(a).sum())
        sh, sp = np.sign(H), np.sign(P)
        ok = (sh != 0) & (sp != 0)
        # where P vanishes it predicts nothing, so it is scored as a
        # coin for the lean rather than dropped: the lean needs a
        # sign at every k or the two f are not comparable
        sp = np.where(sp == 0, 1.0, sp)
        agr = float((sh[ok] == sp[ok]).mean())
        w = np.abs(a)
        f_mu = float(a[a > 0].sum() / l1)
        f_pr = float(w[sp > 0].sum() / l1)
        lm, lp = abs(0.5 - f_mu), abs(0.5 - f_pr)
        eps = (rng.integers(0, 2, size=(DRAWS, ks.size))
               .astype(np.int8) * 2 - 1)
        fl = float(np.median(np.abs((eps @ w) / (2.0 * l1))))
        agrs.append(agr)
        xs.append(math.log(N))
        mus.append(lm)
        prd.append(lp)
        flo.append(fl)
        if agr <= 0.85:
            s2 = False
        r = lp / lm
        if not (1.0 / 1.5 <= r <= 1.5):
            s3 = False
        say("  %-12d %-6d %-10.4f %-9.4f %-9.4f %.4f"
            % (N, ks.size, agr, lm, lp, r))
    say("  S2 the agreement stays above 0.85 there   %s"
        % ("hold" if s2 else "REFUTED"))
    say("  S3 the predicted lean is within a factor 1.5   %s"
        % ("hold" if s3 else "REFUTED"))

    x = np.array(xs)
    bm, rm, sem, tm = fit(x, np.log(np.array(mus)))
    bp, rp, sep, tp = fit(x, np.log(np.array(prd)))
    sd = math.sqrt(sem ** 2 + sep ** 2)
    s4 = abs(bm - bp) <= 2.0 * sd
    say()
    say("S4  and the trend")
    say("  mu        slope %+.6f, standard error %.6f" % (bm, sem))
    say("  P         slope %+.6f, standard error %.6f" % (bp, sep))
    say("  difference %+.6f against %.6f = 2 s.e., i.e. %.2f s.e."
        % (bm - bp, 2.0 * sd, abs(bm - bp) / sd))
    say("SCATTER slope_audit_survivor_range %.4f" % rp)
    say("TSTAT slope_audit_survivor_range %.2f" % tp)
    say("SPREAD slope_audit_survivor_range %.4f"
        % float(x.max() - x.min()))
    if tp < 2.0:
        say("UNRESOLVED SIGN slope_audit_survivor_range")
    say("  S4 %s" % ("hold" if s4 else "REFUTED"))

    say()
    say("  against the floor, and against what Modd gave on the same")
    say("  window in {#rem:oddmertensrange}:")
    say("  N            mu/floor   P/floor")
    for N, lm, lp, fl in zip(NS, mus, prd, flo):
        say("  %-12d %-10.4f %.4f" % (N, lm / fl, lp / fl))
    b1, _r1, s1e, _t1 = fit(x, np.log(np.array(mus) / np.array(flo)))
    b2, _r2, s2e, _t2 = fit(x, np.log(np.array(prd) / np.array(flo)))
    say("  mu %+.6f (s.e. %.6f);  P %+.6f (s.e. %.6f); "
        "difference %.2f s.e."
        % (b1, s1e, b2, s2e,
           abs(b1 - b2) / math.sqrt(s1e ** 2 + s2e ** 2)))

    say()
    say("  DIAGNOSTIC (post hoc). S3 and S4 fail narrowly and in one")
    say("  direction, so how far P is from carrying the lean is worth")
    say("  putting beside Modd, which failed the same tests widely.")
    say("  Modd's ratios are read from")
    say("  results/audit_oddmertens_range.txt:")
    orows = re.findall(r"^  (\d+)\s+\d+\s+\d+\s+[\d.]+\s+[\d.]+"
                       r"\s+[\d.]+\s+([\d.]+)\s*$", wsrc, re.M)
    omr = [float(v) for _n, v in orows]
    prat = [lp / lm for lm, lp in zip(mus, prd)]
    oagr = [float(v) for v in re.findall(
        r"^  \d+\s+\d+\s+\d+\s+([\d.]+)\s+[\d.]+\s+[\d.]+"
        r"\s+[\d.]+\s*$", wsrc, re.M)]
    om_slope = float(re.search(r"Modd      slope ([-+][\d.]+)",
                               wsrc).group(1))
    say("  predictor   ratio to mu's lean   slope       agreement")
    say("  Modd        %.4f to %.4f       %+.6f   %.4f to %.4f"
        % (min(omr), max(omr), om_slope, min(oagr), max(oagr)))
    say("  P           %.4f to %.4f       %+.6f   %.4f to %.4f"
        % (min(prat), max(prat), bp, min(agrs), max(agrs)))
    say("  mu          %.4f by definition   %+.6f"
        % (float(np.mean(np.array(mus) / np.array(mus))), bm))
    say("CARRIES lab_survivor_selection lab_lean_decay %.4f %.4f"
        % (min(prat), max(prat)))
    say("CARRIES lab_lean_oddmertens lab_lean_decay %.4f %.4f"
        % (min(omr), max(omr)))
    say("  so P misses by %.0f to %.0f per cent and decays %.2f"
        % (100 * (1 - max(prat)), 100 * (1 - min(prat)),
           abs(bm - bp) / sd))
    say("  standard errors too fast, where Modd overshoots by %.0f to"
        % (100 * (min(omr) - 1)))
    say("  %.0f per cent and does not decay at all. Neither carries"
        % (100 * (max(omr) - 1)))
    say("  the lean on its own window; P is much the closer.")

    say()
    say("  the two windows of inner length N/k, side by side:")
    say("    the agreement was demonstrated on   [2, %d]" % XHI)
    wm = re.search(r"^WINDOW lab_lean_decay (\d+) (\d+)\s*$",
                   wsrc, re.M)
    wlo, whi = int(wm.group(1)), int(wm.group(2))
    say("    the lean is measured on             [%d, %d]," % (wlo, whi))
    say("    read from results/audit_oddmertens_range.txt")
    say("WINDOW lab_survivor_selection 2 %d" % XHI)
    say("WINDOW lab_lean_decay %d %d" % (wlo, whi))
    if 2 <= wlo and whi <= XHI:
        say("EXPLAINS lab_survivor_selection lab_lean_decay")
    else:
        say("WINDOWS DISJOINT lab_survivor_selection lab_lean_decay")

    say()
    say("=" * 70)
    ok = s1 and s2 and s3 and s4
    say("the sieve-weighted predictor carries the lean where the "
        "lean is measured" if ok else "REFUTED")

    head = [
        "STATISTIC: the agreement of sign H(N;k) with the sign of the",
        "           sieve-weighted P(N;k) = sum_m mu(m) w(m,k), first",
        "           on the published window 2 <= N/k <= " + str(XHI),
        "           over every admissible k rather than a subsample,",
        "           then on k < N^" + str(THETA) + ", the range the",
        "           mass-weighted lean f is measured on; the lean that",
        "           results from replacing each sign H by sign P on",
        "           mu's own magnitudes; and both against the median",
        "           lean of " + str(DRAWS) + " sign vectors on those",
        "           magnitudes.",
        "NULL: for S1 the published permutation control stands; for S3",
        "      and S4 the floor is " + str(DRAWS) + " sign vectors on",
        "      the identical magnitudes, the convention of",
        "      audit_lean_floor.py, so the predicted lean is read on",
        "      the same scale as mu's and as Modd's in",
        "      audit_oddmertens_range.py.",
        "FIELD: N = 2e5 through 6.4e6 by doubling; k squarefree and",
        "       coprime to N; m over 1 <= m < N/k with (m,k) = 1 for",
        "       H and over the odd squarefree such m for P, as",
        "       published; the sieve weight over the odd primes below "
        + str(QSIEVE) + ";",
        "       Lambda, mu and Modd from an integer sieve to "
        + str(NMAX) + ";",
        "       numpy default_rng seed " + str(SEED) + ". Every N is",
        "       2^a 5^b, one odd radical, as RADICALS declares. The",
        "       published agreements are read from",
        "       results/lab_survivor_selection.txt and the lean's",
        "       window from results/audit_oddmertens_range.txt.",
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
