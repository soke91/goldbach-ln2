# -*- coding: utf-8 -*-
r"""
Does the affordable level cut the demand, or only the sign?

WHAT IS AT STAKE

Remark {#rem:levelthreshold} located the level at which the wall's
sign stops slipping: between N^0.2 and N^0.3, well below the
theta' = 0.56 the reduction consumes. What it also found is that
holding the sign is not the same as capturing it -- at alpha = 0.3 the
predicted lean is 0.87 to 0.91 of mu's, flat but short. OPEN's
question is whether that shortfall is usable or whether it is the
residue.

The currency is the demand, not the lean. Remark {#rem:predictable}
measures what subtracting an elementary predictor takes off
B_H = sum(log k)|H|: with beta the least-squares scale through the
origin, sum(log k)|H - beta P| is 0.6310 down to 0.5307 of it, at the
level 29. The same statistic at level N^alpha says what a higher level
would buy, and its far end is fixed by arithmetic rather than by
fitting: at alpha = 1/2 the survivors are exactly the primes, so
H and P_{1/2} run over the same terms and differ only by the weights
log p. Whatever residual share survives there is the spread of log p,
not a Moebius-prime correlation at all.

So the sweep brackets the question. Its far end is the floor an
unbounded level could reach; alpha = 0.3 is what the route can afford.

BACKS: Remark {#rem:leveldemand} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  The control: at Q = 29 the beta and the residual share
      reproduce {#rem:predictable}'s 2.6992, 2.8588, 2.9952, 3.1437,
      3.0473 and 0.6310, 0.5979, 0.5554, 0.5421, 0.5307 to within
      0.0001 and 0.001.
  U2  Depth cuts the demand: the residual share falls with alpha at
      every N.
  U3  And the far end is far: at alpha = 0.5, where the terms are the
      same and only the weights differ, the residual share is below
      0.35 at every N.
  U4  The affordable level gets most of the way: at alpha = 0.3 the
      residual share is below 0.45 at every N.

REFUTATION RULE (fixed before the run)

  U1  REFUTED at either cap -- not the same statistic, and nothing
      below may be compared with {#rem:predictable}.
  U2  REFUTED if the share rises anywhere. Level would then not be
      the axis for the demand as it is for the sign, and the two
      statistics would be measuring different things.
  U3  REFUTED at or above 0.35. The spread of log p alone would then
      account for more than a third of the demand, and no predictor
      of any level could cut it further -- the residue would be
      mostly the weights and not the correlation.
  U4  REFUTED at or above 0.45 at any N. That is the one that
      matters: the level the route can afford would then buy little
      of the demand even though it holds the sign, and the shortfall
      of {#rem:levelthreshold} would be the residue rather than
      something a bounded effort reaches.

  All four gate.

  THE NULL is lab_predictable_part.py's, audited and run in
  audit_predictable_null.py: breaking either half of the predictor
  alone reaches at most 15 per cent of the cut, so the cut is a fact
  about the per-k pairing. It is not repeated at every level here.
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
OUT = os.path.join(RES, "audit_level_demand.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
ALPHAS = [0.1, 0.2, 0.3, 0.4, 0.5]
QFIXED = 29


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


def survivors(n, q):
    """not divisible by any odd prime at or below q"""
    s = np.ones(n + 1, dtype=bool)
    s[0] = False
    for p in primes_upto(max(q, 3)):
        p = int(p)
        if p == 2 or p > q:
            continue
        s[p::p] = False
    return s


def ckconst(k, q):
    """the sieve weight's constant, prod q/(q-1) over q <= Q, q not
    dividing k -- a positive scale, so it cannot change a sign"""
    c = 1.0
    for p in primes_upto(max(q, 3)):
        p = int(p)
        if p == 2 or p > q or k % p == 0:
            continue
        c *= p / (p - 1.0)
    return c


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
    """beta and the residual share of {#rem:predictable}"""
    src = io.open(os.path.join(RES, "lab_predictable_part.txt"),
                  encoding="utf-8").read()
    i = src.index("N            beta      residual share   "
                  "|beta P| share")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            if not out:
                continue
            break
        out[int(f[0])] = (float(f[1]), float(f[2]))
    return out


def share(H, P, w):
    """refit beta through the origin and return the residual share"""
    b = float((H * P).sum() / (P * P).sum())
    return b, float((w * np.abs(H - b * P)).sum()
                    / (w * np.abs(H)).sum())


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d published beta and residual shares from "
        "results/lab_predictable_part.txt" % len(pub))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    oddsqf = sqf.copy()
    oddsqf[::2] = False
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    lvl = {}
    for N in NS:
        lvl[(N, 0.0)] = QFIXED
        for al in ALPHAS:
            lvl[(N, al)] = max(3, int(round(N ** al)))
    say("the levels swept, Q = N^alpha:")
    say("  N            " + "".join("a=%-8.1f" % a
                                    for a in [0.0] + ALPHAS))
    for N in NS:
        say("  %-12d %s" % (N, "".join(
            "%-10d" % lvl[(N, a)] for a in [0.0] + ALPHAS)))
    surv = {}
    for q in sorted(set(lvl.values())):
        surv[q] = survivors(NMAX, q)

    betas, shares = {}, {}
    for N in NS:
        PN = factor_set(N)
        here = sorted(set(lvl[(N, a)] for a in [0.0] + ALPHAS))
        ks, Hs = [], []
        Ps = dict((q, []) for q in here)
        for k in range(2, KCAP):
            if not sqf[k]:
                continue
            if any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
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
            g = mu[ms].astype(np.float64)
            ks.append(k)
            Hs.append(float((lam[vals] * g).sum()))
            for q in here:
                Ps[q].append(ckconst(k, q)
                             * float(g[surv[q][vals]].sum()))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        w = np.log(ks.astype(np.float64))
        for al in [0.0] + ALPHAS:
            b, sh = share(H, np.array(Ps[lvl[(N, al)]]), w)
            betas[(N, al)] = b
            shares[(N, al)] = sh
        say("  N = %-10d #k = %-6d done" % (N, ks.size))

    # ------------------------------------------------------------- U1
    say()
    say("U1  the control: beta and the residual share at Q = %d"
        % QFIXED)
    say("  N            beta here   beta pub   share here  share pub")
    u1 = True
    for N in NS:
        pb, ps = pub[N]
        if (abs(betas[(N, 0.0)] - pb) >= 0.0001
                or abs(shares[(N, 0.0)] - ps) >= 0.001):
            u1 = False
        say("  %-12d %-11.4f %-10.4f %-11.4f %.4f"
            % (N, betas[(N, 0.0)], pb, shares[(N, 0.0)], ps))
    say("  U1 %s   (cap 0.0001 in beta, cap 0.001 in the share)"
        % ("hold" if u1 else "REFUTED"))

    # -------------------------------------------------- U2 / U3 / U4
    say()
    say("U2/U3/U4  the residual share as the level rises")
    say("  N            " + "".join("a=%-8.1f" % a
                                    for a in [0.0] + ALPHAS))
    u2 = u3 = u4 = True
    for N in NS:
        row = [shares[(N, a)] for a in [0.0] + ALPHAS]
        if any(row[i + 1] > row[i] for i in range(len(row) - 1)):
            u2 = False
        if row[-1] >= 0.35:
            u3 = False
        if shares[(N, 0.3)] >= 0.45:
            u4 = False
        say("  %-12d %s" % (N, "".join("%-10.4f" % v for v in row)))
    say("  U2 the share falls with alpha at every N   %s"
        % ("hold" if u2 else "REFUTED"))
    if not u2:
        mono = all(shares[(N, ALPHAS[i + 1])] <= shares[(N, ALPHAS[i])]
                   for N in NS for i in range(len(ALPHAS) - 1))
        say("  DIAGNOSTIC on U2 (post hoc). The ladder it is checked")
        say("  on is not ordered by its own index: the fixed Q = %d"
            % QFIXED)
        say("  is a deeper sieve than N^0.1 and N^0.2 at every N here")
        say("  -- %s against %s and %s -- so the first steps go"
            % (QFIXED,
               ", ".join(str(lvl[(N, 0.1)]) for N in NS),
               ", ".join(str(lvl[(N, 0.2)]) for N in NS)))
        say("  backwards in depth, not forwards. On the alpha ladder")
        say("  proper the share is monotone at every N: %s."
            % ("yes" if mono else "no"))
        say("LADDER leveldemand_with_fixed %d unordered"
            % (len(ALPHAS) + 1))
        say("LADDER leveldemand_alpha %d %s"
            % (len(ALPHAS), "ordered" if mono else "unordered"))
    say("  U3 alpha = 0.5 is below 0.35   %s   (cap 0.35)"
        % ("hold" if u3 else "REFUTED"))
    say("  U4 alpha = 0.3 is below 0.45   %s   (cap 0.45)"
        % ("hold" if u4 else "REFUTED"))
    say("PERN leveldemand_a03 %d %.4f %.4f"
        % (len(NS), min(shares[(N, 0.3)] for N in NS),
           max(shares[(N, 0.3)] for N in NS)))
    say("PERN leveldemand_a05 %d %.4f %.4f"
        % (len(NS), min(shares[(N, 0.5)] for N in NS),
           max(shares[(N, 0.5)] for N in NS)))
    say("  and the two taken at the same N, which is the factor by")
    say("  which the affordable level still leaves more than the")
    say("  unbounded one:")
    rr = [shares[(N, 0.3)] / shares[(N, 0.5)] for N in NS]
    for N, v in zip(NS, rr):
        say("  %-12d %.4f" % (N, v))
    say("RATIO leveldemand_a03 leveldemand_a05 %.4f %.4f"
        % (min(rr), max(rr)))

    say()
    say("  and the scale each level fits, which says what the")
    say("  predictor is being asked to be:")
    say("  N            " + "".join("a=%-8.1f" % a
                                    for a in [0.0] + ALPHAS))
    for N in NS:
        say("  %-12d %s" % (N, "".join(
            "%-10.4f" % betas[(N, a)] for a in [0.0] + ALPHAS)))
    say("  at alpha = 0.5 the terms are the same and only the")
    say("  weights differ, so beta there is a mean of log p; log N")
    say("  runs %.4f to %.4f over this sweep."
        % (math.log(min(NS)), math.log(max(NS))))

    say()
    say("  how far the affordable level goes, as a fraction of the")
    say("  whole distance from Q = %d to alpha = 0.5:" % QFIXED)
    say("  N            share at 29  at 0.3     at 0.5     fraction")
    fr = []
    for N in NS:
        a0, a3, a5 = (shares[(N, 0.0)], shares[(N, 0.3)],
                      shares[(N, 0.5)])
        v = (a0 - a3) / (a0 - a5)
        fr.append(v)
        say("  %-12d %-12.4f %-10.4f %-10.4f %.4f"
            % (N, a0, a3, a5, v))
    say("  the affordable level buys %.0f to %.0f per cent of what an"
        % (100 * min(fr), 100 * max(fr)))
    say("  unbounded one would.")

    say()
    say("  two thresholds on one axis, each taken by its own rule.")
    say("  The sign's is read from")
    say("  results/audit_level_threshold.txt; the demand's is the")
    say("  smallest alpha whose share is below half the level-%d"
        % QFIXED)
    say("  share at every N:")
    tsrc = io.open(os.path.join(RES, "audit_level_threshold.txt"),
                   encoding="utf-8").read()
    sthr = re.search(r"^LEVEL sieve_alpha_threshold (\S+)\s*$",
                     tsrc, re.M).group(1)
    dthr = "none"
    for al in ALPHAS:
        if all(shares[(N, al)] < 0.5 * shares[(N, 0.0)] for N in NS):
            dthr = "%.1f" % al
            break
    say("    the sign's lean stops slipping at alpha = %s" % sthr)
    say("    the demand halves at alpha = %s" % dthr)
    say("AXIS sieve_alpha leanratio %s" % sthr)
    say("AXIS sieve_alpha demandshare %s" % dthr)
    if sthr != dthr:
        say("THRESHOLDS DIFFER sieve_alpha")
        say("  so the two are not the same level, and the cheaper one")
        say("  may not be carried over to the other statistic.")

    say()
    say("=" * 70)
    ok = u1 and u2 and u3 and u4
    say("the affordable level cuts the demand, not only the sign"
        if ok else "REFUTED")

    head = [
        "STATISTIC: with beta the least-squares scale through the",
        "           origin, the residual share",
        "           sum(log k)|H - beta P_Q| / sum(log k)|H| over",
        "           2 <= k < " + str(KCAP) + ", where P_Q is the",
        "           sieve-weighted Moebius sum whose survivors are the",
        "           m with N - mk free of odd prime factors at or",
        "           below Q, at the fixed Q = " + str(QFIXED)
        + " and at Q = N^alpha",
        "           for alpha = " + ", ".join("%.1f" % a
                                              for a in ALPHAS) + ".",
        "NULL: lab_predictable_part.py's, audited and run in",
        "      audit_predictable_null.py -- breaking either half of",
        "      the predictor alone reaches at most 15 per cent of the",
        "      cut. It is not repeated at every level here.",
        "FIELD: N = 2e5 through 3.2e6 by doubling, the field of",
        "       lab_predictable_part.py; k squarefree and coprime to N",
        "       with 2 <= k < " + str(KCAP) + "; m odd, squarefree and",
        "       coprime to k, m <= (N-1)/k; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + ". The sieve weight's",
        "       constant is kept so that beta is on the published",
        "       scale; the skip of primes dividing k is vacuous for",
        "       the survivor set, as audit_sieve_depth.py shows. Every",
        "       N is 2^a 5^b, one odd radical, as RADICALS declares.",
        "       The published beta and shares are read from",
        "       results/lab_predictable_part.txt.",
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
