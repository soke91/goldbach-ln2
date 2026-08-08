# -*- coding: utf-8 -*-
r"""
Where does the cross-k cancellation fall short -- head or tail?

WHAT IS AT STAKE

Remark {#rem:flatnessshape} reduced OPEN item 5 to one number. The
sign lean grows against its floor at theta'/2 - e(G), and asks whether
e(G) has any reason to rise from its measured 0.153911 to the 0.28
that would stop it. G = l1/|sum a| is the gain cancellation across
dilations buys, and Remark {#rem:nocrossk} read its shortfall as the
dilated walls moving together.

If they move together, the shortfall should be localised. Remark
{#rem:nocrossk}'s own rule T4 measured the top decile of k carrying
0.3486 to 0.3587 of sum(log k)|H| -- concentrated but not dominant --
and read that as ruling out a heavy tail. It did not ask the
complementary question: whether the many small terms cancel at the
square-root rate while a few large ones do not. Splitting the k by
mass rank answers it, and the split is by a FRACTION of the k so that
each part still has #S of order N^theta' and the same square-root
reference theta'/2.

If the tail's exponent sits at 0.28 and the head's does not, the
obstruction is a bounded set of dilations and the asymptotics are the
head's. If both fall short, the shortfall is spread over the whole
range and no rearrangement of the k reaches square root.

BACKS: Remark {#rem:gainsplit} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  The control: G over the whole range reproduces
      {#rem:leandecay}'s 1.834, 1.804, 2.207, 2.588, 2.789, 3.079 to
      within 0.01, and its exponent reproduces
      {#rem:leanidentity}'s 0.153911 to within 0.001.
  Y2  The head carries the shortfall: on the top tenth of k by |a|,
      G's exponent is below the whole range's.
  Y3  And the tail does better: on the bottom nine tenths it is above
      the whole range's.
  Y4  But not square root: the tail's exponent is still below
      theta'/2 = 0.28 by more than two standard errors.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED at either cap -- not the same statistic, and nothing
      below may be compared with {#rem:nocrossk} or
      {#rem:leanidentity}.
  Y2  REFUTED if the head's exponent is at or above the whole
      range's, which would say the largest terms are not where the
      cancellation fails.
  Y3  REFUTED if the tail's is at or below the whole range's. Then
      the shortfall would not be localised by size at all and
      {#rem:nocrossk}'s "the dilated walls move together" would be a
      statement about every k equally.
  Y4  REFUTED if the tail reaches theta'/2. That is the one that
      matters: the small terms would then cancel at square root and
      the whole obstruction would be a bounded head, which is a much
      weaker statement than the one the papers make.

  All four gate.

  NO NULL IS RUN for the split, which is a deterministic partition of
  a measured sequence. The coin arms for G were run in
  audit_crossk_reference.py, where random signs on mu's own
  magnitudes gave 9.94 to 12.98 times mu's gain.
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
OUT = os.path.join(RES, "audit_gain_split.txt")

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
    """the published gain G and the exponent e(G)"""
    src = io.open(os.path.join(RES, "lab_lean_decay.txt"),
                  encoding="utf-8").read()
    i = src.index("N          #k     mass frac +   |0.5 - f|   "
                  "G = 1/|2f-1|")
    g = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            if f and set(f[0]) == {"-"}:
                continue
            if not g:
                continue
            break
        g[int(f[0])] = float(f[4])
    src2 = io.open(os.path.join(RES, "audit_lean_identity.txt"),
                   encoding="utf-8").read()
    e = float(re.search(r"^  G\s+([-+][\d.]+)", src2, re.M).group(1))
    return g, e


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

    pubg, pube = read_published()
    say("read %d published gains from results/lab_lean_decay.txt and "
        "the exponent" % len(pubg))
    say("  %+.6f from results/audit_lean_identity.txt" % pube)

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
        order = np.argsort(-w)
        nh = max(1, int(round(HEAD * ks.size)))
        hd, tl = order[:nh], order[nh:]

        def gain(idx):
            s = abs(float(a[idx].sum()))
            return float(np.abs(a[idx]).sum()) / s if s > 0 else \
                float("inf")

        sh = np.sign(a[hd])
        agree = max(float((sh > 0).mean()), float((sh < 0).mean()))
        rows.append((N, ks.size, nh, gain(np.arange(ks.size)),
                     gain(hd), gain(tl),
                     float(w[hd].sum() / w.sum()), agree))
        say("  N = %-10d #k = %-6d head %-5d G %-8.4f head %-8.4f "
            "tail %-8.4f mass %.4f"
            % (N, ks.size, nh, rows[-1][3], rows[-1][4], rows[-1][5],
               rows[-1][6]))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))

    # ------------------------------------------------------------- Y1
    say()
    say("Y1  the control: the whole-range gain")
    say("  N            G here     published  diff")
    y1 = True
    for N, nk, nh, gf, gh, gt, ms, ag in rows:
        if N not in pubg:
            continue
        d = abs(gf - pubg[N])
        if not (d < 0.01):
            y1 = False
        say("  %-12d %-10.4f %-10.4f %.5f" % (N, gf, pubg[N], d))
    ef, rf, sef, tf = fit(x, np.log(np.array([r[3] for r in rows])))
    if not (abs(ef - pube) < 0.001):
        y1 = False
    say("  exponent %+.6f against the published %+.6f, diff %.6f"
        % (ef, pube, abs(ef - pube)))
    say("  Y1 %s   (cap 0.01 and cap 0.001)"
        % ("hold" if y1 else "REFUTED"))

    # --------------------------------------------------- Y2 / Y3 / Y4
    say()
    say("Y2/Y3/Y4  the gain on each part, and its exponent")
    say("  part           " + "".join("%-11d" % r[0] for r in rows))
    for nm, j in (("whole", 3), ("head tenth", 4), ("tail", 5)):
        say("  %-14s %s" % (nm, "".join("%-11.4f" % r[j]
                                        for r in rows)))
    eh, rh, seh, th = fit(x, np.log(np.array([r[4] for r in rows])))
    et, rt, set_, tt = fit(x, np.log(np.array([r[5] for r in rows])))
    say("  part           exponent     s.e.       t")
    for nm, e, se, t, rms in (("whole", ef, sef, tf, rf),
                              ("head tenth", eh, seh, th, rh),
                              ("tail", et, set_, tt, rt)):
        say("  %-14s %+-12.6f %-10.6f %.2f" % (nm, e, se, t))
        say("SCATTER slope_gainsplit_%s %.4f"
            % (nm.replace(" ", "_"), rms))
        say("TSTAT slope_gainsplit_%s %.2f" % (nm.replace(" ", "_"), t))
        say("SPREAD slope_gainsplit_%s %.4f"
            % (nm.replace(" ", "_"), float(x.max() - x.min())))
        if t < 2.0:
            say("UNRESOLVED SIGN slope_gainsplit_%s"
                % nm.replace(" ", "_"))
    y2 = eh < ef
    y3 = et > ef
    half = THETA / 2.0
    y4 = (half - et) > 2.0 * set_
    say("  Y2 the head is below the whole range   %s"
        % ("hold" if y2 else "REFUTED"))
    say("  Y3 the tail is above it   %s" % ("hold" if y3 else "REFUTED"))
    say("  Y4 and still below theta'/2 = %.2f by more than two s.e. "
        "(%.2f)   %s" % (half, (half - et) / set_,
                         "hold" if y4 else "REFUTED"))

    say()
    say("  and the head's share of the mass, which is what")
    say("  {#rem:nocrossk}'s rule T4 measured:")
    say("  N            top-decile share  same sign in the head")
    for N, nk, nh, gf, gh, gt, ms, ag in rows:
        say("  %-12d %-17.4f %.4f" % (N, ms, ag))
    say("  the head's gain is 1 exactly where that fraction is 1:")
    say("  a positive proportion of the dilations -- a tenth of them,")
    say("  growing like N^%.2f -- carries one sign." % THETA)
    say("GAINSPLIT crossk %+.6f %+.6f %+.6f" % (eh, et, ef))
    say("PERN gainsplit_head_mass %d %.4f %.4f"
        % (len(rows), min(r[6] for r in rows),
           max(r[6] for r in rows)))
    say("PERN gainsplit_head_agree %d %.4f %.4f"
        % (len(rows), min(r[7] for r in rows),
           max(r[7] for r in rows)))
    rr = [r[6] / r[7] for r in rows]
    say("RATIO gainsplit_head_mass gainsplit_head_agree %.4f %.4f"
        % (min(rr), max(rr)))

    say()
    say("=" * 70)
    ok = y1 and y2 and y3 and y4
    say("the shortfall is spread by size, not localised in a head"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the cross-k gain G = sum|a| / |sum a| for",
        "           a_k = (log k)H(N;k) over the squarefree",
        "           k < N^" + str(THETA) + " coprime to N, computed on",
        "           the whole range, on the top tenth of k by |a| and",
        "           on the remaining nine tenths; each part's",
        "           least-squares exponent against log N with its",
        "           standard error; and the head's share of",
        "           sum(log k)|H|.",
        "NULL: none for the split, which is a deterministic partition",
        "      of a measured sequence. The coin arms for G were run in",
        "      audit_crossk_reference.py, where random signs on mu's",
        "      own magnitudes gave 9.94 to 12.98 times mu's gain.",
        "FIELD: N = 2e5 through 2.56e7 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < N^" + str(THETA) + "; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + ". The split is by a",
        "       FIXED FRACTION of the k, so each part has #S of order",
        "       N^theta' and the same square-root reference theta'/2.",
        "       Every N is 2^a 5^b, one odd radical, as RADICALS says.",
        "       The published gains are read from",
        "       results/lab_lean_decay.txt and the exponent from",
        "       results/audit_lean_identity.txt.",
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
