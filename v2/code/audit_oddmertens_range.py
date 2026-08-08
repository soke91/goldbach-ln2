# -*- coding: utf-8 -*-
r"""
Does the Modd mechanism reach the k-range the lean is measured on?

WHAT IS AT STAKE

Remark {#rem:leanodd} names the sign lean: sign H(N;k) agrees with
sign Modd(floor(N/k)) at 0.7657, against 0.5738 for the best of
sixteen permutations and 0.5563 for the full Mertens function. It then
credits that mechanism with explaining Remark {#rem:leandecay}'s
decay -- "a fixed k-range sees a shallower lean as N grows and N/k
with it".

The two statements are measured on different k. The agreement runs
over 2 <= N/k <= 1000 with k up to N/2, so it is a statement about
SHORT inner sums. The lean f whose decay is being explained runs over
k < N^0.56, so its inner sums are N/k > N^0.44 -- at N = 3.2e6 that is
N/k > 428, and the agreement's window stops at 1000. The mechanism
was demonstrated where the statistic mostly is not.

That matters more now than when it was written. Remarks
{#rem:leanfloor} and {#rem:leanextended} have shown that the decay
being explained is the floor moving, not the lean going away, so what
Modd has to account for is no longer what it was credited with.

The implementation is independent of lab_lean_oddmertens.py's.

BACKS: Remark {#rem:oddmertensrange} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  P1  The control: over 2 <= N/k <= 1000 with k < N/2 the agreement
      reproduces the published 0.7704, 0.7692, 0.7693, 0.7669,
      0.7657 and the full-Mertens 0.5472 ... 0.5563, both to within
      0.001.
  P2  On the k-range the lean is actually measured on, k < N^0.56,
      the agreement is materially weaker: below 0.70 at every N,
      because those inner sums are long and Modd is shallow there.
  P3  The mechanism nonetheless carries the lean: replacing each
      sign H by sign Modd(floor(N/k)) on mu's own magnitudes gives a
      predicted |0.5 - f| within a factor two of mu's at every N.
  P4  And its trend: the predicted lean's slope against log N is
      within two standard errors of mu's.

REFUTATION RULE (fixed before the run)

  P1  REFUTED at 0.001 anywhere -- not the same statistic, and
      nothing below may be compared with {#rem:leanodd}.
  P2  REFUTED if the agreement stays at or above 0.70 on that range.
      The mechanism would then transfer to the lean's own k without
      qualification and nothing in {#rem:leanodd} needs changing.
  P3  REFUTED outside a factor two at any N. That is the one that
      matters: it would say Modd does not carry the lean where the
      lean is measured, and {#rem:leanodd}'s "it also explains the
      decay without appeal" would have to be withdrawn.
  P4  REFUTED beyond two standard errors, which is the same
      conclusion for the trend rather than the level.

  All four gate.

  THE NULL for P1 is the published one: sixteen permutations of the
  predictor's signs among the distinct values of floor(N/k),
  preserving both marginal sign distributions. For P3 and P4 the
  floor is the 256 global sign vectors of audit_lean_floor.py on the
  identical magnitudes, so the predicted lean is read on the same
  scale as mu's.
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
OUT = os.path.join(RES, "audit_oddmertens_range.txt")

CTRL = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
XHI = 1000
PERMS = 16
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
    """the published agreements, read from the results file"""
    src = io.open(os.path.join(RES, "lab_lean_oddmertens.txt"),
                  encoding="utf-8").read()
    i = src.index("N            #k        odd Mertens   draws max   "
                  "full Mertens")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            if not out:
                continue
            break
        out[int(f[0])] = (float(f[2]), float(f[4]))
    return out


def hvec(N, kmax, lam, mu, sqf):
    """H(N;k) over the squarefree k < kmax coprime to N"""
    PN = factor_set(N)
    ks, Hs = [], []
    for k in range(2, kmax):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 1:
            continue
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        vals = N - ms * k
        Hs.append(float((lam[vals] * mu[ms].astype(np.float64)).sum()))
        ks.append(k)
    return (np.array(ks, dtype=np.int64), np.array(Hs))


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
    say("read %d published agreements from "
        "results/lab_lean_oddmertens.txt" % len(pub))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    # Modd(x) = sum_{m<=x, m odd} mu(m), and the full Mertens
    om = np.zeros(NMAX + 1, dtype=np.int64)
    om[1::2] = mu[1::2]
    modd = np.cumsum(om)
    mert = np.cumsum(mu.astype(np.int64))
    del om
    say("Modd and the full Mertens built to %d" % NMAX)

    rng = np.random.default_rng(SEED)

    # ------------------------------------------------------------- P1
    say()
    say("P1  the control: agreement over 2 <= N/k <= %d, k < N/2"
        % XHI)
    say("  N            #k        odd Mertens  published   "
        "full Mertens  published")
    p1 = True
    for N in CTRL:
        ks, H = hvec(N, N // 2, lam, mu, sqf)
        inner = N // ks
        sel = (inner >= 2) & (inner <= XHI) & (H != 0)
        s1 = np.sign(H[sel])
        pm = np.sign(modd[inner[sel]])
        pf = np.sign(mert[inner[sel]])
        # the selection is the predictor's own: Modd = 0 says
        # nothing and is dropped, while a full-Mertens zero is a
        # failure to predict on a k the other side does predict, so
        # it counts as a miss -- the published convention
        ok = pm != 0
        a1 = float((s1[ok] == pm[ok]).mean())
        a2 = float((s1[ok] == pf[ok]).mean())
        po, pfu = pub[N]
        if abs(a1 - po) >= 0.001 or abs(a2 - pfu) >= 0.001:
            p1 = False
        say("  %-12d %-9d %-12.4f %-11.4f %-13.4f %.4f"
            % (N, int(ok.sum()), a1, po, a2, pfu))
    say("  P1 %s   (cap 0.001)" % ("hold" if p1 else "REFUTED"))

    # -------------------------------------------------- P2 / P3 / P4
    say()
    say("P2/P3/P4  on the k-range the lean is measured on, "
        "k < N^%.2f" % THETA)
    say("  N            #k     N/k from   agreement  mu lean   "
        "Modd lean  ratio")
    p2 = p3 = True
    xs, mus, prd, flo = [], [], [], []
    inlo, inhi = [], []
    for N in NS:
        ks, H = hvec(N, int(N ** THETA), lam, mu, sqf)
        a = np.log(ks.astype(np.float64)) * H
        l1 = float(np.abs(a).sum())
        inner = N // ks
        pm = np.sign(modd[inner])
        s1 = np.sign(H)
        ok = (pm != 0) & (s1 != 0)
        agr = float((s1[ok] == pm[ok]).mean())
        f_mu = float(a[a > 0].sum() / l1)
        w = np.abs(a)
        f_pr = float(w[pm > 0].sum() / l1)
        lm, lp = abs(0.5 - f_mu), abs(0.5 - f_pr)
        eps = (rng.integers(0, 2, size=(DRAWS, ks.size))
               .astype(np.int8) * 2 - 1)
        fl = float(np.median(np.abs((eps @ w) / (2.0 * l1))))
        inlo.append(int(inner.min()))
        inhi.append(int(inner.max()))
        xs.append(math.log(N))
        mus.append(lm)
        prd.append(lp)
        flo.append(fl)
        if agr >= 0.70:
            p2 = False
        r = lp / lm
        if not (0.5 <= r <= 2.0):
            p3 = False
        say("  %-12d %-6d %-10d %-10.4f %-9.4f %-10.4f %.4f"
            % (N, ks.size, int(inner.min()), agr, lm, lp, r))
    say("  P2 the agreement is below 0.70 there   %s"
        % ("hold" if p2 else "REFUTED"))
    say("  P3 the predicted lean is within a factor two   %s"
        % ("hold" if p3 else "REFUTED"))
    say()
    say("  the two windows of inner length N/k, side by side:")
    say("    the agreement was demonstrated on   [%d, %d]" % (2, XHI))
    say("    the lean is measured on             [%d, %d]"
        % (min(inlo), max(inhi)))
    say("WINDOW lab_lean_oddmertens %d %d" % (2, XHI))
    say("WINDOW lab_lean_decay %d %d" % (min(inlo), max(inhi)))
    if 2 <= min(inlo) and max(inhi) <= XHI:
        say("  the second sits inside the first:")
        say("EXPLAINS lab_lean_oddmertens lab_lean_decay")
    else:
        say("  the second does not sit inside the first -- it runs a")
        say("  factor %.0f past its upper end, so the mechanism was"
            % (max(inhi) / float(XHI)))
        say("  demonstrated where the statistic mostly is not:")
        say("WINDOWS DISJOINT lab_lean_oddmertens lab_lean_decay")

    x = np.array(xs)
    bm, rm, sem, tm = fit(x, np.log(np.array(mus)))
    bp, rp, sep, tp = fit(x, np.log(np.array(prd)))
    sd = math.sqrt(sem ** 2 + sep ** 2)
    p4 = abs(bm - bp) <= 2.0 * sd
    say()
    say("P4  and the trend")
    say("  mu        slope %+.6f, standard error %.6f" % (bm, sem))
    say("  Modd      slope %+.6f, standard error %.6f" % (bp, sep))
    say("  difference %+.6f against %.6f = 2 s.e., i.e. %.2f s.e."
        % (bm - bp, 2.0 * sd, abs(bm - bp) / sd))
    say("SCATTER slope_audit_oddmertens_range %.4f" % rp)
    say("TSTAT slope_audit_oddmertens_range %.2f" % tp)
    say("SPREAD slope_audit_oddmertens_range %.4f"
        % float(x.max() - x.min()))
    if tp < 2.0:
        say("UNRESOLVED SIGN slope_audit_oddmertens_range")
    say("  P4 %s" % ("hold" if p4 else "REFUTED"))

    say()
    say("  and the same two against the floor, which is what")
    say("  {#rem:leanfloor} showed the decay to be:")
    say("  N            mu/floor   Modd/floor")
    for N, lm, lp, fl in zip(NS, mus, prd, flo):
        say("  %-12d %-10.4f %.4f" % (N, lm / fl, lp / fl))
    b1, _r1, s1e, _t1 = fit(x, np.log(np.array(mus) / np.array(flo)))
    b2, _r2, s2e, _t2 = fit(x, np.log(np.array(prd) / np.array(flo)))
    say("  mu   %+.6f (s.e. %.6f);  Modd %+.6f (s.e. %.6f)"
        % (b1, s1e, b2, s2e))
    say("  difference %+.6f, i.e. %.2f s.e."
        % (b1 - b2, abs(b1 - b2) / math.sqrt(s1e ** 2 + s2e ** 2)))

    say()
    say("=" * 70)
    ok = p1 and p2 and p3 and p4
    say("the mechanism is weaker on the lean's own k and still "
        "carries it" if ok else "REFUTED")

    head = [
        "STATISTIC: the fraction of k at which sign H(N;k) equals",
        "           sign Modd(floor(N/k)), first on the published",
        "           window 2 <= N/k <= " + str(XHI) + " with k < N/2,",
        "           then on k < N^" + str(THETA) + ", the range the",
        "           mass-weighted lean f is measured on; the lean that",
        "           results from replacing each sign H by",
        "           sign Modd(floor(N/k)) on mu's own magnitudes; and",
        "           both against the median lean of " + str(DRAWS),
        "           sign vectors on those magnitudes.",
        "NULL: for P1 the published permutation control stands; for",
        "      P3 and P4 the floor is " + str(DRAWS) + " sign vectors",
        "      on the identical magnitudes, the convention of",
        "      audit_lean_floor.py, so the predicted lean is read on",
        "      the same scale as mu's.",
        "FIELD: N = 2e5 through 6.4e6 by doubling; k squarefree and",
        "       coprime to N, below N/2 for P1 and below N^"
        + str(THETA),
        "       for the rest; m over 1 <= m < N/k with (m,k) = 1;",
        "       Lambda, mu, Modd and the full Mertens from an integer",
        "       sieve to " + str(NMAX) + "; k with Modd = 0 or H = 0",
        "       excluded from the agreement; numpy default_rng seed "
        + str(SEED) + ".",
        "       The published agreements are read from",
        "       results/lab_lean_oddmertens.txt.",
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
