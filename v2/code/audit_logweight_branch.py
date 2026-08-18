# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, "Numerical verification", bullet 4 -- the log-weight
branch, and an independent re-test of the sign fixed in Remark
{#rem:sign}.

WHAT IS UNDER TEST

The paper prints, at theta' = 0.56 and t = N-1:

    N            5e4    1e5    2e5    4e5    8e5
    R_log/N      2.28   2.26   2.21   2.13   2.07
    predicted    2.33   2.26   2.19   2.13   2.08

calls this "an object of size ~ N, consistent with Theorem [thm:C]",
and records S(N) = 1.7604 for these N.  Neither script it names for
this bullet exists in this repository.

Definitions used here, taken from the paper:

    R_log(t) = sum_{N-t<=u<N} Lambda(N-u) mu(u)
                 sum_{k|u, k>=K, (k,N)=1} mu(k) log k
             = sum_{N-t<=u<N} Lambda(N-u) mu(u)
                 [ CL(u) - sigmaL_K(u) ],

    CL(u)  = sum_{k|u,(k,N)=1} mu(k) log k = -Lambda(u'),   Lemma
             [lem:completelog], u' the part of u prime to N,
    MT_log = sum_{m<M,(m,N)=1} mu(m) c(m) I(m),   c(m) = A(N)lambda(m)/m,
    I(m)   = int_{mK}^{N} log(v/m) dv,                       Section {#sec:C}(i)
    K      = floor(N^theta'),  M = N/K.

WHY THE SIGN MATTERS HERE

Version 3 of the note evaluated I(m)'s coefficient as
A(N) Gtilde(1) = +S(N), giving a residual main term -S(N)N.  Remark
{#rem:sign} corrects that to A(N) Gtilde(1) = -S(N) and +S(N)N.  The
printed row above is positive and near 2, i.e. near S(N) = 1.7604 from
above, so the printed measurement already contradicts the version-3
sign.  That is checked here as a separate prediction rather than
inferred, because it is a second, independent witness to the same
correction.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  R_log/N reproduces 2.28, 2.26, 2.21, 2.13, 2.07.
  V2  MT_log/N reproduces 2.33, 2.26, 2.19, 2.13, 2.08.
  V3  The log branch does not decay: R_log/N stays above 1.7 at every
      N, whereas the w=1 branch's R/N falls from 0.1140 to 0.0483 over
      the same range.  Concretely R_log/N > 20 * (R/N) at N = 8e5.
  V4  MT_log/N is monotone decreasing across the five N and stays above
      S(N) = 1.7604 at each, i.e. it approaches the singular series from
      above.
  V5  MT_log/N > 0 at every N.  The version-3 sign predicts it negative,
      so this discriminates the two.
  V6  S(N) = 1.7604 at all five N, since every one of them has
      P(N) = {2,5}.

REFUTATION RULE (fixed before the run)

  V1, V2  REFUTED if any entry differs from the printed value by more
          than 0.005, i.e. by more than the printed two decimals allow.
  V3      REFUTED if R_log/N <= 1.7 at any N, or if the ratio at 8e5 is
          at most 20.
  V4      REFUTED if MT_log/N is not monotone decreasing, or dips to or
          below 1.7604 at any of the five N.
  V5      REFUTED if MT_log/N <= 0 at any N -- which would reinstate the
          version-3 sign and refute Remark {#rem:sign}.
  V6      REFUTED if S(N) differs from 1.7604 by more than 0.00005 at
          any N.

  All six gate.
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
OUT = os.path.join(ROOT, "results", "audit_logweight_branch.txt")

THETA = 0.56
NS = [50_000, 100_000, 200_000, 400_000, 800_000]
PUB_R = [2.28, 2.26, 2.21, 2.13, 2.07]
PUB_MT = [2.33, 2.26, 2.19, 2.13, 2.08]
PUB_W1 = [0.1140, 0.0965, 0.0785, 0.0618, 0.0483]     # the w=1 branch
PLIM = 4_000_000


def sieves(n):
    spf = np.zeros(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if spf[p] == 0:
            blk = spf[p::p]
            spf[p::p] = np.where(blk == 0, p, blk)
    mu = np.ones(n + 1, dtype=np.int64)
    mu[0] = 0
    for v in range(2, n + 1):
        p = int(spf[v])
        w = v // p
        mu[v] = 0 if w % p == 0 else -mu[w]
    lam = np.zeros(n + 1, dtype=np.float64)
    for p in range(2, n + 1):
        if int(spf[p]) != p:
            continue
        q, lg = p, math.log(p)
        while q <= n:
            lam[q] = lg
            if q > n // p:
                break
            q *= p
    return spf, mu, lam


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s)


def prime_set(N, spf):
    s = set()
    while N > 1:
        p = int(spf[N])
        s.add(p)
        while N % p == 0:
            N //= p
    return s


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ps = primes_upto(PLIM)
    twin = 2.0
    for p in ps[1:]:
        twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    spf, mu, lam = sieves(NMAX)

    say()
    say("  N        K     M       S(N)     R_log/N   pub     MT_log/N  pub  "
        "   |R-MT|/R   R_{w=1}/N")
    say("  " + "-" * 92)
    Rn, MTn, Ss = [], [], []
    for N in NS:
        PN = prime_set(N, spf)
        A = 1.0
        for p in ps:
            p = int(p)
            if p == 2 or p in PN:
                continue
            A *= 1.0 - 1.0 / (p * (p - 1.0))
        S = twin
        for p in sorted(PN):
            if p > 2:
                S *= 1.0 + 1.0 / (p - 2.0)
        Ss.append(S)

        K = int(N ** THETA)
        M = N / K
        t = N - 1

        # sigmaL_K(u) = sum_{k|u, k<K, (k,N)=1} mu(k) log k
        sigL = np.zeros(N, dtype=np.float64)
        sig = np.zeros(N, dtype=np.int64)
        for k in range(2, K):
            if mu[k] == 0:
                continue
            if any(k % p == 0 for p in PN):
                continue
            sigL[k::k] += int(mu[k]) * math.log(k)
        for k in range(1, K):
            if mu[k] == 0:
                continue
            if any(k % p == 0 for p in PN):
                continue
            sig[k::k] += int(mu[k])

        # CL(u) = -Lambda(u'), u' the part of u prime to N; and 1_{rad(u)|N}
        CL = np.zeros(N, dtype=np.float64)
        comp = np.zeros(N, dtype=np.int64)
        for u in range(1, N):
            if mu[u] == 0:
                continue
            v, outs = u, []
            while v > 1:
                p = int(spf[v])
                if p not in PN:
                    outs.append(p)
                while v % p == 0:
                    v //= p
            comp[u] = 1 if not outs else 0
            if len(outs) == 1:
                CL[u] = -math.log(outs[0])

        n = np.arange(1, t + 1, dtype=np.int64)
        f = lam[1:t + 1] * mu[N - n]
        Rlog = float((f * (CL[N - n] - sigL[N - n])).sum())
        Rw1 = float((f * (comp[N - n] - sig[N - n])).sum())

        # MT_log = A(N) sum_{m<M,(m,N)=1} mu(m) lambda(m)/m * I(m)
        tot = 0.0
        m = 1
        while m < M:
            if mu[m] != 0:
                ok, lm, v = True, 1.0, m
                while v > 1:
                    p = int(spf[v])
                    if p in PN:
                        ok = False
                    lm /= (1.0 - 1.0 / (p * (p - 1.0)))
                    while v % p == 0:
                        v //= p
                if ok:
                    I = (N * math.log(N / m) - N
                         - m * K * math.log(K) + m * K)
                    tot += int(mu[m]) * lm / m * I
            m += 1
        MTlog = A * tot

        i = NS.index(N)
        Rn.append(Rlog / N)
        MTn.append(MTlog / N)
        say("  %-8d %-5d %-7.1f %-8.4f %-9.4f %-7.2f %-9.4f %-7.2f %-10.4f "
            "%.4f" % (N, K, M, S, Rlog / N, PUB_R[i], MTlog / N, PUB_MT[i],
                      abs(Rlog - MTlog) / abs(Rlog), Rw1 / N))

    say()
    e1 = max(abs(a - b) for a, b in zip(Rn, PUB_R))
    e2 = max(abs(a - b) for a, b in zip(MTn, PUB_MT))
    v1, v2 = e1 <= 0.005, e2 <= 0.005
    say("V1  max |R_log/N recomputed - printed|   = %.5f  (tol 0.005)  %s"
        % (e1, "hold" if v1 else "REFUTED"))
    say("V2  max |MT_log/N recomputed - printed|  = %.5f  (tol 0.005)  %s"
        % (e2, "hold" if v2 else "REFUTED"))

    ratio = Rn[-1] / PUB_W1[-1]
    v3 = all(r > 1.7 for r in Rn) and ratio > 20
    say("V3  min R_log/N = %.4f (floor 1.7); R_log/N over R_{w=1}/N at "
        "8e5 = %.1f (floor 20)   %s"
        % (min(Rn), ratio, "hold" if v3 else "REFUTED"))

    mono = all(MTn[i] > MTn[i + 1] for i in range(len(MTn) - 1))
    v4 = mono and all(m > 1.7604 for m in MTn)
    say("V4  MT_log/N monotone decreasing: %s; min = %.4f above S = %.4f"
        "   %s" % (mono, min(MTn), Ss[0], "hold" if v4 else "REFUTED"))

    v5 = all(m > 0 for m in MTn)
    say("V5  MT_log/N > 0 at every N: %s -- the version-3 sign predicts "
        "negative   %s" % (v5, "hold" if v5 else "REFUTED"))

    v6 = all(abs(s - 1.7604) <= 5e-5 for s in Ss)
    say("V6  S(N) = %s   against printed 1.7604   %s"
        % (", ".join("%.5f" % s for s in Ss), "hold" if v6 else "REFUTED"))

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4 and v5 and v6
    say("V1 %s  V2 %s  V3 %s  V4 %s  V5 %s  V6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (v1, v2, v3, v4, v5, v6)))
    say("the log-weight bullet reproduces, and Remark {#rem:sign} gets "
        "a second witness" if ok else "REFUTED")

    head = [
        "STATISTIC: R_log(t)/N with R_log the k>=K residual of the log-",
        "           weighted switch at t=N-1, evaluated as",
        "           sum Lambda(N-u) mu(u) [CL(u) - sigmaL_K(u)] with",
        "           CL(u) = -Lambda(u') of Lemma {#lem:completelog};",
        "           MT_log/N with MT_log = A(N) sum mu(m)lambda(m)/m I(m),",
        "           I(m) = int_{mK}^{N} log(v/m) dv; |R_log-MT_log|/R_log;",
        "           the singular series S(N); and R/N for the w=1 branch",
        "           computed on the same field for contrast.",
        "FIELD: N = 5e4, 1e5, 2e5, 4e5, 8e5 with theta' = 0.56, so",
        "       K = floor(N^0.56) and M = N/K; t = N-1 throughout;",
        "       Lambda and mu from an integer sieve to 8e5; A(N) and the",
        "       twin constant as Euler products over p < 4e6.",
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
