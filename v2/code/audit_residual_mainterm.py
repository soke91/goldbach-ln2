# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, "Numerical verification", bullets 2 and 3 --
the residual R(t) against its predicted main term MT(t), and R's scale.

WHAT IS UNDER TEST

The paper prints, at theta' = 0.56 and t = N-1:

    N        5e4      1e5      2e5      4e5      8e5
    R/N      0.1140   0.0965   0.0785   0.0618   0.0483
    MT/N     0.1190   0.0964   0.0770   0.0616   0.0491

calls this "1-4% agreement", says the main term decays through the
cancellation of Lemma [lem:mu], and adds that the observed decay is
R ~ N^{1-theta'/2}, the ratio R/N^{1-theta'/2} staying in [2.17, 2.46]
over the whole range, rather than R ~ N.  Neither script it names for
these two bullets exists in this repository.

Definitions used here, taken from the paper and not re-derived:

    R(t) = sum_{N-t<=u<N} Lambda(N-u) mu(u)
             sum_{k|u, k>=K, (k,N)=1} mu(k)                    [eq:PR]
    MT(t) = sum_{m<M, (m,N)=1} mu(m) T_m(t) c(m),              [eq:R3]
    c(m)  = A(N) lambda(m) / m,                          Lemma [lem:density]
    T_m(t) = min(t, N - mK),   M = N/K,   K = floor(N^theta').    [eq:Tm]

The (m,N)=1 restriction in MT is the one Remark [rem:trap] says is
load-bearing, so it is tested rather than assumed.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  R/N reproduces 0.1140, 0.0965, 0.0785, 0.0618, 0.0483.
  U2  MT/N reproduces 0.1190, 0.0964, 0.0770, 0.0616, 0.0491.
  U3  |R - MT| / R <= 0.04 at every N -- the paper's "1-4%".
  U4  R / N^{1-theta'/2} lies in [2.17, 2.46] at every N.
  U5  The (m,N)=1 restriction is load-bearing: dropping it moves MT/N
      outside the 4% band of U3 at every N.
  U6  R is not of size N: R/N falls monotonically across the five N,
      and the fitted exponent b in R ~ N^b is closer to
      1 - theta'/2 = 0.72 than to 1.

REFUTATION RULE (fixed before the run)

  U1, U2  REFUTED if any entry differs from the printed value by more
          than 0.0001, i.e. by more than the printed precision allows.
  U3      REFUTED if |R - MT|/R > 0.04 at any N.
  U4      REFUTED if the ratio leaves [2.17, 2.46] at any N.
  U5      REFUTED if the unrestricted MT stays within 4% of R at any N.
  U6      REFUTED if R/N is not monotone, or if |b - 1| < |b - 0.72|.

  All six gate: the script exits non-zero if any fails.  U1 and U2 are
  the transcription check; U3-U6 are the claims the bullets actually
  make.

CITED BY: {#rem:band} in paper/.
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
OUT = os.path.join(ROOT, "results", "audit_residual_mainterm.txt")

THETA = 0.56
NS = [50_000, 100_000, 200_000, 400_000, 800_000]
PUB_R = [0.1140, 0.0965, 0.0785, 0.0618, 0.0483]
PUB_MT = [0.1190, 0.0964, 0.0770, 0.0616, 0.0491]
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


def A_of_N(PN):
    # the bound is the module constant, not the caller's list:
    # audit_constants.py shows the truncation reaches the sixth
    # printed decimal, so it must not depend on the call site.
    v = 1.0
    for p in primes_upto(PLIM):
        p = int(p)
        if p == 2 or p in PN:
            continue
        v *= 1.0 - 1.0 / (p * (p - 1.0))
    return v


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ps = primes_upto(PLIM)
    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    spf, mu, lam = sieves(NMAX)

    say()
    say("  N        K     M       R/N       pub       MT/N      pub      "
        " |R-MT|/R   MT_all/N")
    say("  " + "-" * 86)
    Rn, MTn, MTall, Rabs = [], [], [], []
    for N in NS:
        PN = prime_set(N, spf)
        A = A_of_N(PN)
        K = int(N ** THETA)
        M = N / K
        t = N - 1

        # ---- R(t) by [eq:PR], on the u side
        sig = np.zeros(N, dtype=np.int64)          # sigma_K(u), k < K
        for k in range(1, K):
            if mu[k] == 0:
                continue
            if any(k % p == 0 for p in PN):
                continue
            sig[k::k] += int(mu[k])
        comp = np.zeros(N, dtype=np.int64)         # 1_{rad(u)|N}
        for u in range(1, N):
            if mu[u] == 0:
                continue
            v, ok = u, True
            while v > 1:
                p = int(spf[v])
                if p not in PN:
                    ok = False
                    break
                while v % p == 0:
                    v //= p
            comp[u] = 1 if ok else 0
        n = np.arange(1, t + 1, dtype=np.int64)
        f = lam[1:t + 1] * mu[N - n]
        R = float((f * (comp[N - n] - sig[N - n])).sum())

        # ---- MT(t) by [eq:R3] with c(m) = A(N) lambda(m)/m
        def mt(restrict):
            tot = 0.0
            m = 1
            while m < M:
                if mu[m] != 0:
                    ok = True
                    lm = 1.0
                    v = m
                    while v > 1:
                        p = int(spf[v])
                        if p in PN:
                            ok = False
                        lm /= (1.0 - 1.0 / (p * (p - 1.0)))
                        while v % p == 0:
                            v //= p
                    if ok or not restrict:
                        Tm = min(t, N - m * K)
                        tot += int(mu[m]) * Tm * A * lm / m
                m += 1
            return tot

        MT = mt(True)
        MT0 = mt(False)
        Rn.append(R / N)
        MTn.append(MT / N)
        MTall.append(MT0 / N)
        Rabs.append(R)
        i = NS.index(N)
        say("  %-8d %-5d %-7.1f %-9.4f %-9.4f %-9.4f %-9.4f %-10.4f %.4f"
            % (N, K, M, R / N, PUB_R[i], MT / N, PUB_MT[i],
               abs(R - MT) / abs(R), MT0 / N))

    say()
    e1 = max(abs(a - b) for a, b in zip(Rn, PUB_R))
    e2 = max(abs(a - b) for a, b in zip(MTn, PUB_MT))
    u1, u2 = e1 <= 1e-4, e2 <= 1e-4
    say("U1  max |R/N recomputed - printed|   = %.6f   (tol 0.0001)  %s"
        % (e1, "hold" if u1 else "REFUTED"))
    say("U2  max |MT/N recomputed - printed|  = %.6f   (tol 0.0001)  %s"
        % (e2, "hold" if u2 else "REFUTED"))

    rel = [abs(r - m) / abs(r) for r, m in zip(Rn, MTn)]
    u3 = max(rel) <= 0.04
    say("U3  max |R-MT|/R = %.4f   (cap 0.04)                      %s"
        % (max(rel), "hold" if u3 else "REFUTED"))

    say()
    exp = 1.0 - THETA / 2.0
    sc = [r / N ** exp for r, N in zip(Rabs, NS)]
    u4 = all(2.17 <= s <= 2.46 for s in sc)
    say("U4  R / N^%.2f = %s" % (exp, ", ".join("%.3f" % s for s in sc)))
    say("    band [2.17, 2.46]                                     %s"
        % ("hold" if u4 else "REFUTED"))

    relall = [abs(r - m) / abs(r) for r, m in zip(Rn, MTall)]
    u5 = all(x > 0.04 for x in relall)
    say("U5  |R-MT_all|/R = %s"
        % ", ".join("%.4f" % x for x in relall))
    say("    all above the 4%% band                                 %s"
        % ("hold" if u5 else "REFUTED"))

    mono = all(Rn[i] > Rn[i + 1] for i in range(len(Rn) - 1))
    x = np.log(np.array(NS, dtype=float))
    y = np.log(np.array(Rabs, dtype=float))
    b = float(np.polyfit(x, y, 1)[0])
    u6 = mono and abs(b - exp) < abs(b - 1.0)
    say("U6  R/N monotone decreasing: %s;  fitted R ~ N^%.4f"
        % (mono, b))
    say("    against 1-theta'/2 = %.2f and 1.00                     %s"
        % (exp, "hold" if u6 else "REFUTED"))

    say()
    say("=" * 70)
    ok = u1 and u2 and u3 and u4 and u5 and u6
    say("U1 %s  U2 %s  U3 %s  U4 %s  U5 %s  U6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (u1, u2, u3, u4, u5, u6)))
    say("the residual-vs-main-term bullets reproduce" if ok else "REFUTED")

    head = [
        "STATISTIC: R(t)/N with R the k>=K residual of [eq:PR] at t=N-1;",
        "           MT(t)/N with MT the main term of [eq:R3] and",
        "           c(m)=A(N)lambda(m)/m of Lemma {#lem:density};",
        "           |R-MT|/R; the same with the (m,N)=1 restriction on MT",
        "           dropped; R/N^{1-theta'/2}; and the exponent b in a",
        "           least-squares fit of log R against log N.",
        "FIELD: N = 5e4, 1e5, 2e5, 4e5, 8e5 with theta' = 0.56, so",
        "       K = floor(N^0.56) and M = N/K; t = N-1 throughout;",
        "       Lambda and mu from an integer sieve to 8e5; A(N) as an",
        "       Euler product over p < 4e6 with p=2 and p|N omitted.",
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
