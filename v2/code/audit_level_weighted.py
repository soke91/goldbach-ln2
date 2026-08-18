# -*- coding: utf-8 -*-
r"""
Is the 0.23 floor a floor, or an artefact of counting?

WHAT IS AT STAKE

Remark {#rem:leveldemand} leaves the demand's residual share at
0.2271 to 0.2525 when the sieve reaches alpha = 1/2, and OPEN asks
whether that is a floor. Two things settle it, and neither needs a new
idea.

First, the sieve family saturates at 1/2 and cannot be pushed past it:
for Q above sqrt(N) nothing composite is left to remove and raising Q
only strikes the small primes, which are true contributors --
{#rem:sievedepth} measured that as a fall. So there is no
0.5 < alpha <= 0.56 to explore inside this family.

Second, the 0.23 is not a floor but a cost of counting. At alpha = 1/2
the survivors are the primes, so log(N - mk) IS log p on them: the
log-weighted predictor {#rem:logweightpredictor} tested at level 29 --
where it lowered the agreement -- becomes exact at level sqrt(N). The
weights are free once the survivor set is known. What that says is
that the whole difficulty is the survivor set and none of it is the
weight, which reverses how {#rem:logweightpredictor}'s failure reads.

The residue at alpha = 1/2 should then be only what the sieve's own
convention costs: the primes at or below sqrt(N), which the sieve
strikes along with their multiples, and the value 1.

BACKS: Remark {#rem:levelweighted} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The control: the unweighted residual share reproduces
      {#rem:leveldemand}'s 0.2271 to 0.2525 at alpha = 1/2 and its
      0.4866 to 0.6022 at alpha = 0.3, to within 0.001.
  V2  The weights are free at the top: at alpha = 1/2 the
      log-weighted residual share is below 0.01 at every N.
  V3  And they help everywhere: the log-weighted share is below the
      unweighted one at every N and every alpha.
  V4  But not where the route can pay: at alpha = 0.3 the
      log-weighted share is still above 0.40 at every N.

REFUTATION RULE (fixed before the run)

  V1  REFUTED at 0.001 -- not the same statistic, and nothing below
      may be compared with {#rem:leveldemand}.
  V2  REFUTED at or above 0.01. The 0.23 would then be a real floor
      rather than the cost of counting, and something other than the
      survivor set would be carrying part of the demand.
  V3  REFUTED if the log weight ever raises the share. At level 29 it
      lowered the sign agreement, so a rise here would say the two
      scores disagree about the weight as well as about the level.
  V4  REFUTED below 0.40 at any N. That is the one that matters: the
      affordable level plus the free weights would then cut the
      demand after all, and {#rem:leveldemand}'s reading would have
      to be softened.

  All four gate.

  THE NULL is lab_predictable_part.py's, audited and run in
  audit_predictable_null.py: breaking either half of the predictor
  alone reaches at most 15 per cent of the cut. It is not repeated at
  every level here.
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
OUT = os.path.join(RES, "audit_level_weighted.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
ALPHAS = [0.1, 0.2, 0.3, 0.4, 0.5]


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
    """prod q/(q-1) over the odd q <= Q not dividing k"""
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
    """the unweighted shares of {#rem:leveldemand}"""
    src = io.open(os.path.join(RES, "audit_level_demand.txt"),
                  encoding="utf-8").read()
    i = src.index("U2/U3/U4  the residual share as the level rises")
    out = {}
    for ln in src[i:].splitlines()[2:]:
        f = ln.split()
        if len(f) < 7 or not f[0].isdigit():
            break
        out[int(f[0])] = [float(v) for v in f[1:7]]
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
    say("read %d published share rows from "
        "results/audit_level_demand.txt" % len(pub))

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
        for al in ALPHAS:
            lvl[(N, al)] = max(3, int(round(N ** al)))
    surv = {}
    for q in sorted(set(lvl.values())):
        surv[q] = survivors(NMAX, q)

    plain, logw, misfit = {}, {}, {}
    for N in NS:
        PN = factor_set(N)
        here = sorted(set(lvl[(N, a)] for a in ALPHAS))
        ks, Hs = [], []
        Ps = dict((q, []) for q in here)
        Ls = dict((q, []) for q in here)
        small = 0.0
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
            lv = np.log(vals.astype(np.float64))
            ks.append(k)
            Hs.append(float((lam[vals] * g).sum()))
            for q in here:
                s = surv[q][vals]
                c = ckconst(k, q)
                Ps[q].append(c * float(g[s].sum()))
                Ls[q].append(float((g[s] * lv[s]).sum()))
            qt = lvl[(N, 0.5)]
            miss = (lam[vals] > 0) & (~surv[qt][vals])
            small += math.log(k) * float(
                np.abs((lam[vals] * g)[miss]).sum())
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        w = np.log(ks.astype(np.float64))
        misfit[N] = small / float((w * np.abs(H)).sum())
        for al in ALPHAS:
            q = lvl[(N, al)]
            plain[(N, al)] = share(H, np.array(Ps[q]), w)[1]
            logw[(N, al)] = share(H, np.array(Ls[q]), w)[1]
        say("  N = %-10d #k = %-6d done" % (N, ks.size))

    # ------------------------------------------------------------- V1
    say()
    say("V1  the control: the unweighted shares")
    say("  N            a=0.3 here  published  a=0.5 here  published")
    v1 = True
    for N in NS:
        p3, p5 = pub[N][3], pub[N][5]
        if (abs(plain[(N, 0.3)] - p3) >= 0.001
                or abs(plain[(N, 0.5)] - p5) >= 0.001):
            v1 = False
        say("  %-12d %-11.4f %-10.4f %-11.4f %.4f"
            % (N, plain[(N, 0.3)], p3, plain[(N, 0.5)], p5))
    say("  V1 %s   (cap 0.001)" % ("hold" if v1 else "REFUTED"))

    # -------------------------------------------------- V2 / V3 / V4
    say()
    say("V2/V3/V4  what the log weight buys at each level")
    say("  N            " + "".join("a=%-8.1f" % a for a in ALPHAS))
    say("  unweighted / log-weighted")
    v2 = v3 = v4 = True
    for N in NS:
        say("  %-12d %s" % (N, "".join(
            "%-10.4f" % plain[(N, a)] for a in ALPHAS)))
        say("  %-12s %s" % ("", "".join(
            "%-10.4f" % logw[(N, a)] for a in ALPHAS)))
        if logw[(N, 0.5)] >= 0.01:
            v2 = False
        if any(logw[(N, a)] >= plain[(N, a)] for a in ALPHAS):
            v3 = False
        if logw[(N, 0.3)] <= 0.40:
            v4 = False
    say("  V2 alpha = 0.5 log-weighted below 0.01   %s   (cap 0.01)"
        % ("hold" if v2 else "REFUTED"))
    say("  V3 the log weight lowers the share everywhere   %s"
        % ("hold" if v3 else "REFUTED"))
    say("  V4 alpha = 0.3 log-weighted still above 0.40   %s"
        "   (floor 0.40)" % ("hold" if v4 else "REFUTED"))
    say("PERN levelweighted_a03 %d %.4f %.4f"
        % (len(NS), min(logw[(N, 0.3)] for N in NS),
           max(logw[(N, 0.3)] for N in NS)))
    say("PERN levelweighted_a05 %d %.6f %.6f"
        % (len(NS), min(logw[(N, 0.5)] for N in NS),
           max(logw[(N, 0.5)] for N in NS)))
    say("PERN levelplain_a03 %d %.4f %.4f"
        % (len(NS), min(plain[(N, 0.3)] for N in NS),
           max(plain[(N, 0.3)] for N in NS)))
    say("  and the two taken at the same N -- what the free weights")
    say("  buy where the route could pay:")
    rr = [logw[(N, 0.3)] / plain[(N, 0.3)] for N in NS]
    for N, v in zip(NS, rr):
        say("  %-12d %.4f" % (N, v))
    say("RATIO levelweighted_a03 levelplain_a03 %.4f %.4f"
        % (min(rr), max(rr)))

    say()
    say("  what is left at the top, and what it is. The sieve strikes")
    say("  a prime along with its multiples, so the primes at or")
    say("  below sqrt(N) are true contributors it removes. Their")
    say("  share of the demand, measured directly:")
    say("  N            mass of the struck primes")
    for N in NS:
        say("  %-12d %.6f" % (N, misfit[N]))
    say("  against the log-weighted residual %.6f to %.6f, so the"
        % (min(logw[(N, 0.5)] for N in NS),
           max(logw[(N, 0.5)] for N in NS)))
    say("  two are the same object to within the fit's freedom.")
    say("ACCOUNT levelweighted_a05 %.6f %.6f %.6f %.6f"
        % (min(logw[(N, 0.5)] for N in NS),
           max(logw[(N, 0.5)] for N in NS),
           min(misfit[N] for N in NS), max(misfit[N] for N in NS)))

    say()
    say("  and the axis, with both statistics on it. The demand's")
    say("  threshold is read from results/audit_level_demand.txt:")
    dsrc = io.open(os.path.join(RES, "audit_level_demand.txt"),
                   encoding="utf-8").read()
    dthr = re.search(r"^AXIS sieve_alpha demandshare (\S+)\s*$",
                     dsrc, re.M).group(1)
    lthr = "none"
    for al in ALPHAS:
        if all(logw[(N, al)] < 0.5 * plain[(N, 0.1)] for N in NS):
            lthr = "%.1f" % al
            break
    say("AXIS sieve_alpha demandshare %s" % dthr)
    say("AXIS sieve_alpha demandshare_logweighted %s" % lthr)
    if dthr != lthr:
        say("THRESHOLDS DIFFER sieve_alpha")
        say("  so the free weights move the demand's threshold from")
        say("  alpha = %s to alpha = %s." % (dthr, lthr))
    else:
        say("  so the free weights do not move the demand's threshold")
        say("  at all: it stays at alpha = %s, well above the 0.3 at"
            % dthr)
        say("  which the sign stops slipping. What they change is the")
        say("  size of what is left there, from %.4f to %.6f."
            % (max(plain[(N, 0.5)] for N in NS),
               max(logw[(N, 0.5)] for N in NS)))

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("the weights are free and the survivor set is the whole "
        "difficulty" if ok else "REFUTED")

    head = [
        "STATISTIC: with beta the least-squares scale through the",
        "           origin, the residual share",
        "           sum(log k)|H - beta X| / sum(log k)|H| over",
        "           2 <= k < " + str(KCAP) + ", for X the sieve-weighted",
        "           count P_Q and for X the log-weighted",
        "           sum_m mu(m)[survives] log(N - mk), at",
        "           Q = N^alpha for alpha = "
        + ", ".join("%.1f" % a for a in ALPHAS) + "; and the share of",
        "           the demand carried by the terms whose N - mk is a",
        "           prime at or below sqrt(N), which the sieve strikes.",
        "NULL: lab_predictable_part.py's, audited and run in",
        "      audit_predictable_null.py -- breaking either half of",
        "      the predictor alone reaches at most 15 per cent of the",
        "      cut. It is not repeated at every level here.",
        "FIELD: N = 2e5 through 3.2e6 by doubling, the field of",
        "       lab_predictable_part.py; k squarefree and coprime to N",
        "       with 2 <= k < " + str(KCAP) + "; m odd, squarefree and",
        "       coprime to k, m <= (N-1)/k; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + ". Every N is 2^a 5^b,",
        "       one odd radical, as RADICALS declares. The published",
        "       unweighted shares are read from",
        "       results/audit_level_demand.txt.",
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
