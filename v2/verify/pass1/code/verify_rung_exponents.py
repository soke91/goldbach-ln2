# -*- coding: utf-8 -*-
r"""
Pass 1, blind: the primorial ladder's level exponent, reimplemented.

WHAT IS BEING VERIFIED

Item 1 of OPEN now rests on the rung exponents: {#rem:laddercurve}
reads curvature out of them, {#rem:curvebound} and {#rem:rung15} put
theta' a decade earlier than the published shape on the strength of
them, and each rung above the eleventh has exactly one witness -- one
script, run once. That is the protocol's first priority: numbers with
a single witness that a step rests on.

This script does not import or re-run any of that code. It builds its
own sieves and its own loops from the definitions the paper states,
and asks whether the same numbers come out. Where the paper leaves a
convention implicit, the choice made here is written down so the
disagreement, if any, can be attributed.

THE DEFINITIONS, AS READ FROM paper/theorem_A.md

  k    squarefree, coprime to N, 2 <= k < 100000.
  m    odd, squarefree, coprime to k, m <= (N-1)/k.
  H    sum_m Lambda(N - mk) mu(m).
  S    the m above with N - mk not divisible by any odd prime q <= 29
       that does not divide k -- {#rem:provablehalf}'s
       "m not congruent to N k^{-1} mod q for all q <= 29".
  P    C_k sum_{m in S} mu(m), with C_k the product of q/(q-1) over
       those same q.
  beta the least-squares coefficient sum(HP)/sum(P^2).
  R    H - beta P.
  K*_R the first k at which sum_{k'<=k} (log k')|R(N;k')| reaches
       S(N)(1 - A(N))N.
  S(N) 2 prod_{p>2} (1 - 1/(p-1)^2) times prod_{q|N, q>2} (1 + 1/(q-2)).
  A(N) prod_p (1 - 1/(p(p-1))) divided by prod_{q|N} the same factor.
  The Euler products are cut at the fixed bound 4000000.

The exponent is log K*_R / log N.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The three rungs reproduce the published 0.5023, 0.5099 and
      0.5178 to the bound their printing forces, 0.00005.
  V2  And so does K*_R itself, exactly: 5773, 9367, 15461.
  V3  beta reproduces to 0.000001 at each rung.

REFUTATION RULE (fixed before the run)

  V1  REFUTED outside 0.00005 at any rung. Then the exponent this
      repository publishes is not the one its own definitions give,
      and every reading built on the ladder is suspect until the
      difference is found.
  V2  REFUTED on any mismatch. K*_R is an integer and a reimplementation
      that agrees on the exponent but not on it has a different
      threshold or a different cumulative sum.
  V3  REFUTED outside 0.000001. beta is a ratio of two sums over the
      same index set; a disagreement there localises the difference to
      P, since H is the simpler object.

  All three gate: this is a verification pass and a failure is its
  finding.

  NO NULL IS RUN and none applies. Two implementations of one
  definition are compared.
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

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "verify_rung_exponents.txt")

BASE = 30030
RUNGS = (10, 11, 12)
KCAP = 100_000
QMAX = 29
CLIM = 4_000_000
PUB = {10: (5773, 0.5023), 11: (9367, 0.5099), 12: (15461, 0.5178)}


def smallest_prime_factor(n):
    """spf[i] = the least prime dividing i, by a plain sieve"""
    spf = np.zeros(n + 1, dtype=np.int32)
    spf[2::2] = 2
    i = 3
    while i * i <= n:
        if spf[i] == 0:
            spf[i * i::2 * i] = np.where(spf[i * i::2 * i] == 0, i,
                                         spf[i * i::2 * i])
        i += 2
    odd = np.arange(3, n + 1, 2)
    z = odd[spf[odd] == 0]
    spf[z] = z
    if n >= 1:
        spf[1] = 1
    return spf


def lam_from_spf(n, spf):
    """Lambda by walking each prime's powers -- independent of mu"""
    lam = np.zeros(n + 1, dtype=np.float64)
    for p in range(2, n + 1):
        if int(spf[p]) != p:
            continue
        lp = math.log(p)
        q = p
        while q <= n:
            lam[q] = lp
            if q > n // p:
                break
            q *= p
    return lam


def mu_from_spf(n, spf):
    mu = np.zeros(n + 1, dtype=np.int8)
    mu[1] = 1
    for x in range(2, n + 1):
        p = int(spf[x])
        y = x // p
        mu[x] = 0 if (y % p == 0) else -mu[y]
    return mu


def odd_primes_upto(m):
    out = []
    for p in range(3, m + 1, 2):
        ok = True
        for d in range(3, int(p ** 0.5) + 1, 2):
            if p % d == 0:
                ok = False
                break
        if ok:
            out.append(p)
    return out


def factorset(x):
    out, d = set(), 2
    while d * d <= x:
        if x % d == 0:
            out.add(d)
            while x % d == 0:
                x //= d
        d += 1
    if x > 1:
        out.add(x)
    return out


def euler(N):
    """S(N) and A(N), cut at the fixed bound"""
    spf = smallest_prime_factor(CLIM)
    twin, artin = 2.0, 1.0
    for p in range(2, CLIM + 1):
        if int(spf[p]) != p:
            continue
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    for q in factorset(N):
        artin /= 1.0 - 1.0 / (q * (q - 1.0))
        if q > 2:
            twin *= 1.0 + 1.0 / (q - 2.0)
    return twin, artin


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("verify pass 1 -- the ladder's level exponent, reimplemented")
    say("  sealed targets (sha256 at the top of FINDINGS.md)")
    say("  no code is imported from v2/code; the sieves, the sums and")
    say("  the Euler products are written here from the definitions")
    say("  quoted in the docstring.")

    NMAX = BASE * (1 << max(RUNGS))
    say()
    say("building the smallest-prime-factor table to %d ..." % NMAX)
    spf = smallest_prime_factor(NMAX)
    lam = lam_from_spf(NMAX, spf)
    mu = mu_from_spf(NMAX, spf)
    sqf = mu != 0
    qs = odd_primes_upto(QMAX)
    say("  odd primes to %d: %s" % (QMAX, ", ".join(map(str, qs))))

    say()
    say("  rung  N              beta        K*_R      exponent   "
        "published K*_R / exponent")
    rows = []
    worst_e, worst_b, kmiss = 0.0, 0.0, 0
    for j in RUNGS:
        N = BASE * (1 << j)
        PN = factorset(N)
        twin, artin = euler(N)
        thr = twin * (1.0 - artin) * N
        ks, Hs, Ps = [], [], []
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
            for q in factorset(k):
                if q > 2:
                    ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            v = N - ms * k
            g = mu[ms].astype(np.float64)
            ck = 1.0
            keep = np.ones(ms.size, dtype=bool)
            for q in qs:
                if k % q == 0:
                    continue
                ck *= q / (q - 1.0)
                keep &= (v % q) != 0
            ks.append(k)
            Hs.append(float((lam[v] * g).sum()))
            Ps.append(ck * float(g[keep].sum()))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        P = np.array(Ps)
        beta = float((H * P).sum() / (P * P).sum())
        R = H - beta * P
        cum = np.cumsum(np.log(ks.astype(np.float64)) * np.abs(R))
        idx = int(np.searchsorted(cum, thr))
        kstar = int(ks[idx])
        e = math.log(kstar) / math.log(N)
        pk, pe = PUB[j]
        worst_e = max(worst_e, abs(e - pe))
        if kstar != pk:
            kmiss += 1
        rows.append((j, N, beta, kstar, e, pk, pe))
        say("  %-5d %-14d %-11.6f %-9d %-10.4f %d / %.4f"
            % (j, N, beta, kstar, e, pk, pe))

    say()
    say("V1  the exponents")
    v1 = worst_e <= 0.00005
    say("  worst departure %.6f against the cap 0.00005" % worst_e)
    say("  V1 %s" % ("hold" if v1 else "REFUTED"))

    say()
    say("V2  K*_R itself")
    v2 = kmiss == 0
    say("  mismatches: %d of %d" % (kmiss, len(RUNGS)))
    say("  V2 %s" % ("hold" if v2 else "REFUTED"))

    say()
    say("V3  beta")
    say("  beta is not published to compare against; the three values")
    say("  are printed above and this prediction is withdrawn as")
    say("  unfalsifiable against the paper as sealed. The published")
    say("  betas live in results/, which a blind pass does not read.")
    v3 = None
    say("  V3 WITHDRAWN")

    say()
    say("=" * 70)
    say("V1 %s  V2 %s  V3 WITHDRAWN"
        % ("hold" if v1 else "REFUTED", "hold" if v2 else "REFUTED"))

    head = [
        "STATISTIC: the primorial ladder's level exponent",
        "           log K*_R / log N at N = 30030*2^j for j = 10, 11,",
        "           12, computed from an independent implementation of",
        "           the definitions quoted in the script's docstring:",
        "           an own smallest-prime-factor sieve, own Lambda and",
        "           mu built from it, own sifted set and own Euler",
        "           products at the fixed bound 4000000; against the",
        "           K*_R and exponents the paper publishes.",
        "NULL: none is run and none applies. Two implementations of",
        "      one definition are compared.",
        "FIELD: N = 30030*2^j, j = 10, 11, 12; k squarefree and",
        "       coprime to N with 2 <= k < 100000; m odd, squarefree",
        "       and coprime to k with m <= (N-1)/k; the sifted set S",
        "       removes m with N - mk divisible by an odd prime",
        "       q <= 29 not dividing k, and C_k is the product of",
        "       q/(q-1) over those q; beta = sum(HP)/sum(P^2);",
        "       R = H - beta P; the threshold is S(N)(1-A(N))N with",
        "       both Euler products cut at 4000000.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % os.path.abspath(OUT))
    if not (v1 and v2):
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
