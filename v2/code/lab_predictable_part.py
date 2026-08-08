# -*- coding: utf-8 -*-
r"""
How much of the demand the elementary predictor actually accounts for,
over the k-range the reduction uses.

WHAT IS AT STAKE

Remark {#rem:survivors} found that the sign of H(N;k) is given, at
about 91 percent, by the sieve-weighted Mobius sum
P(N;k) = sum_m mu(m) w(m,k) over its own inner range, and called the
remaining 9 percent the place any saving would have to live.  That
measurement was taken over 2 <= N/k <= 1000 -- short inner sums, which
is large k.

The demand the reduction has to bound is
B_H(N;K) = sum_{k<K}(log k)|H(N;k)| with K near N^{0.7}, and there the
inner sums are long: at N = 3.2e6 the crossing sits at k of order
2.3e4, where N/k is of order 1e2, while the sum runs down to k = 2
where N/k is of order 1e6.  Over a long range sum_m mu(m) cancels and
P should be small, while |H| keeps its square-root size.  So the
predictor may be explaining the sign exactly where the mass is not.

That is the question: where does the mass of B_H sit, and does the
elementary part survive there?

BACKS: Remark {#rem:predictable} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  The predictor is a short-range phenomenon: over k with
      N/k <= 1000 the correlation of H with P exceeds 0.5, and over k
      with N/k >= 1e4 it is under 0.2, at every N.
  S2  The mass is at the long end: the share of sum(log k)|H| over
      2 <= k < 30000 that comes from k with N/k >= 1e4 exceeds 0.5 at
      every N.
  S3  So subtracting the elementary part barely helps: with beta the
      least-squares scale through the origin,
      sum(log k)|H - beta P| / sum(log k)|H| stays above 0.9 at
      every N.
  S4  And the predictor is void where the mass is: the sign agreement
      of H with P restricted to N/k >= 1e4 lies in [0.45, 0.55].

REFUTATION RULE (fixed before the run)

  S1  REFUTED if either bound is crossed at any N.
  S2  REFUTED if the share is at or under 0.5 at any N, in which case
      the mass is short-range after all and [rem:survivors] speaks
      about the demand rather than beside it.
  S3  REFUTED if the ratio drops to 0.9 or below at any N. A
      refutation would be the good outcome: it would mean an
      elementary sum carries a definite share of the demand.
  S4  REFUTED if the agreement leaves that band, which would say some
      structure survives into the long range.

  All four gate.

  NO SIGN CONTROL IS RUN and none is needed. S1 and S4 are already
  comparisons against the no-information value for a sign predictor,
  and S2 and S3 are decompositions of one measured sum, where a
  randomisation would move both parts together. The permutation nulls
  for this predictor were run in lab_survivor_selection.py.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_predictable_part.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
SHORT, LONG = 1_000, 10_000


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


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("  sieve weight over the odd primes %s"
        % ", ".join(map(str, QS)))

    res = []
    for N in NS:
        PN = factor_set(N)
        ks, Hs, Ps = [], [], []
        for k in range(2, KCAP):
            if not sqf[k] or any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
            if M < 2:
                continue
            ms = np.arange(1, M + 1, 2, dtype=np.int64)
            ms = ms[sqf[ms]]
            for q in factor_set(k):
                if q > 2:
                    ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            vals = N - ms * k
            g = mu[ms].astype(np.float64)
            H = float((lam[vals] * g).sum())
            w = np.ones(ms.size, dtype=np.float64)
            for q in QS:
                if k % q == 0:
                    continue
                w *= np.where(vals % q == 0, 0.0, q / (q - 1.0))
            P = float((g * w).sum())
            ks.append(k)
            Hs.append(H)
            Ps.append(P)
        ks = np.array(ks, dtype=np.int64)
        Hs = np.array(Hs)
        Ps = np.array(Ps)
        Minner = N // ks
        beta = float((Hs * Ps).sum() / (Ps * Ps).sum())
        res.append((N, ks, Hs, Ps, Minner, beta))
        say("  N = %-10d  #k = %-7d beta = %.4f" % (N, ks.size, beta))

    say()
    say("S1  correlation of H with P, by the length of the inner sum")
    say("  N            N/k <= %-6d N/k >= %-6d" % (SHORT, LONG))
    s1 = True
    for N, ks, Hs, Ps, Minner, beta in res:
        out = []
        for sel in ((Minner <= SHORT), (Minner >= LONG)):
            if sel.sum() > 2 and Ps[sel].std() > 0:
                out.append(float(np.corrcoef(Hs[sel], Ps[sel])[0, 1]))
            else:
                out.append(float("nan"))
        if not (out[0] > 0.5 and out[1] < 0.2):
            s1 = False
        say("  %-12d %-14.4f %.4f" % (N, out[0], out[1]))
    say("  S1 %s" % ("hold" if s1 else "REFUTED"))

    say()
    say("S2  where the mass of sum (log k)|H| sits")
    say("  N            share at N/k >= %-6d  share at N/k <= %d"
        % (LONG, SHORT))
    s2 = True
    for N, ks, Hs, Ps, Minner, beta in res:
        w = np.log(ks.astype(np.float64)) * np.abs(Hs)
        tot = float(w.sum())
        a = float(w[Minner >= LONG].sum()) / tot
        b = float(w[Minner <= SHORT].sum()) / tot
        if a <= 0.5:
            s2 = False
        say("  %-12d %-22.4f %.4f" % (N, a, b))
    say("  S2 %s" % ("hold" if s2 else "REFUTED"))

    say()
    say("S3  what subtracting the elementary part leaves")
    say("  N            beta      residual share   |beta P| share")
    s3 = True
    for N, ks, Hs, Ps, Minner, beta in res:
        lw = np.log(ks.astype(np.float64))
        tot = float((lw * np.abs(Hs)).sum())
        r = float((lw * np.abs(Hs - beta * Ps)).sum()) / tot
        q = float((lw * np.abs(beta * Ps)).sum()) / tot
        if r <= 0.9:
            s3 = False
        say("  %-12d %-9.4f %-16.4f %.4f" % (N, beta, r, q))
    say("  S3 %s" % ("hold" if s3 else "REFUTED"))

    say()
    say("S4  is the predictor void where the mass is?")
    say("  N            sign agreement at N/k >= %d" % LONG)
    s4 = True
    for N, ks, Hs, Ps, Minner, beta in res:
        sel = (Minner >= LONG) & (Hs != 0) & (Ps != 0)
        a = float((np.sign(Hs[sel]) == np.sign(Ps[sel])).mean())
        if not (0.45 <= a <= 0.55):
            s4 = False
        say("  %-12d %.4f" % (N, a))
    say("  S4 %s" % ("hold" if s4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). The two regimes side by side, so the")
    say("  split can be read directly:")
    say("  N            |H| mean short   |H| mean long   |P| mean short"
        "   |P| mean long")
    for N, ks, Hs, Ps, Minner, beta in res:
        sh = Minner <= SHORT
        lo = Minner >= LONG
        say("  %-12d %-16.2f %-15.2f %-16.4f %.4f"
            % (N, float(np.abs(Hs[sh]).mean()),
               float(np.abs(Hs[lo]).mean()),
               float(np.abs(Ps[sh]).mean()),
               float(np.abs(Ps[lo]).mean())))
    say("  |H| grows towards the long end because the inner sum is")
    say("  longer there, while |P| collapses because sum_m mu(m)")
    say("  cancels over a long range. The elementary predictor is")
    say("  strong exactly where the terms are small.")

    say()
    say("  Cross-check lines. lab_residue_size.py fits the same beta")
    say("  on the same k-range before splitting off the residue.")
    for N, ks, Hs, Ps, Minner, beta in res:
        say("AGREE beta_HP N=%d %.6f 0.01" % (N, beta))

    say()
    say("=" * 70)
    ok = s1 and s2 and s3 and s4
    say("the elementary predictor explains the sign where the mass is "
        "not" if ok else "REFUTED")

    head = [
        "STATISTIC: the correlation of H(N;k) with the sieve-weighted",
        "           P(N;k) = sum_m mu(m)w(m,k), split by the length N/k of",
        "           the inner sum; the share of sum(log k)|H| carried by",
        "           each side of that split; the residual share after",
        "           subtracting beta P with beta the least-squares scale",
        "           through the origin; and the sign agreement where the",
        "           inner sums are long.",
        "NULL: none is run and none is needed here. S1 and S4 are already",
        "      stated against the no-information value for a sign",
        "      predictor, and S2 and S3 decompose one measured sum, where",
        "      a randomisation would move both parts together. The",
        "      permutation nulls for this predictor were run in",
        "      lab_survivor_selection.py.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree, coprime",
        "       to N, 2 <= k < 30000, which covers the crossing K* at",
        "       every N; m odd squarefree, coprime to k, m <= (N-1)/k; the",
        "       sieve weight uses the odd primes up to 30.",
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
