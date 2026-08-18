# -*- coding: utf-8 -*-
r"""
paper/wall_v3.md, Proposition {#prop:W} and the absolute budget
[eq:budget] under it.

WHAT IS UNDER TEST

With c(h) = sum_{p'-p=h, p,p'<=X} (log p)(log p') and V as in
Proposition {#prop:V}, the paper defines

    Gamma(N) = ( sum_{h != 0} c(h) ) / V(N)

and asserts Gamma(N) ~ N / (A(N) log N), so that a hypothesis
|S(h)| <= eps yields nothing better than |rho - 1| <= eps Gamma(N).
It prints

    Gamma = 1.5489e3, 1.8517e4, 3.5798e5   at N = 1e4, 1.6e5, 4e6,
    Gamma log N / N = 1.4266, 1.3868, 1.3605   against 1/A -> 1.270.

It then defines the absolute budget

    B(X) = (1/V(X)) sum_{0<|h|<X} c(h) |S(h)|,
    S(h) = (1/(X-|h|)) sum_{|h|<u<=X} mu(u) mu(u-|h|),

and prints B(X) = 13.3, 17.3, 23.1, 30.5 at X = 2e4, 4e4, 8e4, 1.6e5,
adding that normalising S by X rather than by the number X-|h| of terms
gives 7.9, 10.3, 13.8, 18.2.  No script for any of this exists here.

Note sum_{h != 0} c(h) = theta(X)^2 - sum_{p<=X} (log p)^2 exactly, so
the numerator of Gamma needs no correlation at all; the correlation is
needed only for B.  Both are computed here, the first in closed form
and the second by FFT autocorrelation, which makes them independent.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  Gamma reproduces 1.5489e3, 1.8517e4, 3.5798e5.
  Y2  Gamma log N / N reproduces 1.4266, 1.3868, 1.3605, and the row is
      decreasing.
  Y3  1/A(N) = 1.27020 at all three N, since P(N) = {2,5} for each.
  Y4  B(X) reproduces 13.3, 17.3, 23.1, 30.5; every entry exceeds 1 and
      the row is increasing.  This is what closes the absolute-value
      route, so it is the substantive one.
  Y5  The X-normalised variant reproduces 7.9, 10.3, 13.8, 18.2.
  Y6  Gamma A(N) log N / N -> 1: it is within 15% of 1 at N = 4e6.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED if any entry differs by more than 0.05% relative.
  Y2  REFUTED if any entry differs by more than 0.00005, or if the row
      is not decreasing.
  Y3  REFUTED if 1/A(N) differs from 1.27020 by more than 0.000005 at
      any of the three N.
  Y4  REFUTED if any entry differs by more than 0.05, or if any entry is
      at most 1, or if the row is not increasing.
  Y5  REFUTED if any entry differs by more than 0.05.
  Y6  REFUTED if the ratio is outside [0.85, 1.15] at N = 4e6.

  All six gate.  Y2 and Y6 are two forms of the same asymptotic and both
  are kept: Y2 is the printed row, Y6 is the claim the row is offered
  as evidence for, and they can fail separately.

CITED BY: {#rem:cdef} in paper/.
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
OUT = os.path.join(ROOT, "results", "audit_amplification.txt")

GAMMA_NS = [10_000, 160_000, 4_000_000]
PUB_GAMMA = [1.5489e3, 1.8517e4, 3.5798e5]
PUB_GLN = [1.4266, 1.3868, 1.3605]
BUDGET_XS = [20_000, 40_000, 80_000, 160_000]
PUB_B = [13.3, 17.3, 23.1, 30.5]
PUB_B2 = [7.9, 10.3, 13.8, 18.2]
ARTIN_LIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def mobius_upto(n):
    mu = np.ones(n + 1, dtype=np.int64)
    prim = np.ones(n + 1, dtype=bool)
    prim[:2] = False
    for p in range(2, n + 1):
        if prim[p]:
            prim[p * p::p] = False
            mu[p::p] *= -1
            mu[p * p::p * p] = 0
    return mu


def autocorr(a, n):
    """A(h) = sum_u a(u) a(u-h), h >= 0, by zero-padded rfft."""
    L = 1
    while L < 2 * n + 2:
        L <<= 1
    b = np.zeros(L, dtype=np.float64)
    b[:len(a)] = a
    F = np.fft.rfft(b)
    return np.fft.irfft(F * np.conj(F), L)[:n + 1]


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(max(GAMMA_NS), max(BUDGET_XS))
    say("sieving to %d ..." % NMAX)
    pr = primes_upto(NMAX)
    lgp = np.log(pr.astype(np.float64))
    mu = mobius_upto(NMAX)

    artin = 1.0
    for p in primes_upto(ARTIN_LIM):
        artin *= 1.0 - 1.0 / (int(p) * (int(p) - 1.0))

    def A_of(N):
        v, fac, d = N, set(), 2
        while d * d <= v:
            if v % d == 0:
                fac.add(d)
                while v % d == 0:
                    v //= d
            d += 1
        if v > 1:
            fac.add(v)
        a = artin
        for q in fac:
            a /= (1.0 - 1.0 / (q * (q - 1.0)))
        return a, sorted(fac)

    def V_of(N):
        """sum_{v<N} mu^2(v) Lambda(N-v)^2, by prime powers w = N-v."""
        tot = 0.0
        for i, p in enumerate(pr):
            p = int(p)
            if p >= N:
                break
            q, lg2 = p, lgp[i] ** 2
            while q < N:
                if mu[N - q] != 0:
                    tot += lg2
                if q > (N - 1) // p:
                    break
                q *= p
        return tot

    say()
    say("Y1/Y2/Y3/Y6   the amplification")
    say("=" * 74)
    say("  N          P(N)      1/A(N)     V(N)          Gamma        pub"
        "          GlogN/N   pub")
    say("  " + "-" * 88)
    gs, glns, invAs = [], [], []
    for i, N in enumerate(GAMMA_NS):
        A, fac = A_of(N)
        m = pr < N
        th = float(lgp[m].sum())
        s2 = float((lgp[m] ** 2).sum())
        num = th * th - s2                       # sum_{h != 0} c(h), exact
        V = V_of(N)
        G = num / V
        gln = G * math.log(N) / N
        gs.append(G)
        glns.append(gln)
        invAs.append(1.0 / A)
        say("  %-10d %-9s %-10.6f %-13.4f %-12.4e %-12.4e %-9.4f %.4f"
            % (N, fac, 1.0 / A, V, G, PUB_GAMMA[i], gln, PUB_GLN[i]))

    e1 = max(abs(a - b) / b for a, b in zip(gs, PUB_GAMMA))
    y1 = e1 <= 5e-4
    say("  Y1  max relative deviation = %.5f  (tol 0.0005)   %s"
        % (e1, "hold" if y1 else "REFUTED"))
    e2 = max(abs(a - b) for a, b in zip(glns, PUB_GLN))
    dec = all(glns[i] > glns[i + 1] for i in range(len(glns) - 1))
    y2 = e2 <= 5e-5 and dec
    say("  Y2  max |deviation| = %.6f  (tol 0.00005), decreasing %s   %s"
        % (e2, dec, "hold" if y2 else "REFUTED"))
    y3 = all(abs(a - 1.27020) <= 5e-6 for a in invAs)
    say("  Y3  1/A(N) = %s   %s"
        % (", ".join("%.6f" % a for a in invAs),
           "hold" if y3 else "REFUTED"))
    r6 = glns[-1] / invAs[-1]
    y6 = 0.85 <= r6 <= 1.15
    say("  Y6  Gamma A log N / N at N = 4e6 = %.6f   %s"
        % (r6, "hold" if y6 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc, not a pre-registered test). The gap in Y1")
    say("  is -2.33%, -0.57%, -0.11%, i.e. it shrinks like N^{-1/2}. That")
    say("  is the signature of prime powers: psi(x) - theta(x) ~ sqrt x,")
    say("  so replacing the text's c(h) = sum_{p'-p=h}(log p)(log p') by")
    say("  the Lambda-weighted sum_{n} Lambda(n)Lambda(n+h) multiplies the")
    say("  numerator by (psi/theta)^2 = 1 + 2/sqrt(N) + ...  Both are")
    say("  computed below against the printed row.")
    say("  N          primes (text)   Lambda (prime powers)   published")
    for i, N in enumerate(GAMMA_NS):
        m = pr < N
        th = float(lgp[m].sum())
        s2 = float((lgp[m] ** 2).sum())
        psi, t2 = 0.0, 0.0
        for j, p in enumerate(pr):
            p = int(p)
            if p >= N:
                break
            q = p
            while q < N:
                psi += lgp[j]
                t2 += lgp[j] ** 2
                if q > (N - 1) // p:
                    break
                q *= p
        V = V_of(N)
        say("  %-10d %-15.4e %-23.4e %.4e"
            % (N, (th * th - s2) / V, (psi * psi - t2) / V, PUB_GAMMA[i]))

    say()
    say("Y4/Y5   the absolute budget")
    say("=" * 74)
    say("  X          V(X)         B(X)      pub     B_X(X)    pub"
        "     max|S(h)|")
    say("  " + "-" * 74)
    bs, b2s = [], []
    for i, Xv in enumerate(BUDGET_XS):
        a = np.zeros(Xv + 1, dtype=np.float64)
        m = pr <= Xv
        a[pr[m]] = lgp[m]
        ch = autocorr(a, Xv)                     # c(h), h = 0..X
        am = np.zeros(Xv + 1, dtype=np.float64)
        am[1:] = mu[1:Xv + 1]
        Ah = autocorr(am, Xv)                    # sum mu(u)mu(u-h)
        h = np.arange(1, Xv, dtype=np.float64)
        S = Ah[1:Xv] / (Xv - h)
        S2 = Ah[1:Xv] / Xv
        V = V_of(Xv)
        B = 2.0 * float((ch[1:Xv] * np.abs(S)).sum()) / V
        B2 = 2.0 * float((ch[1:Xv] * np.abs(S2)).sum()) / V
        bs.append(B)
        b2s.append(B2)
        say("  %-10d %-12.2f %-9.4f %-7.1f %-9.4f %-7.1f %.6f"
            % (Xv, V, B, PUB_B[i], B2, PUB_B2[i], float(np.abs(S).max())))

    e4 = max(abs(a - b) for a, b in zip(bs, PUB_B))
    inc = all(bs[i] < bs[i + 1] for i in range(len(bs) - 1))
    y4 = e4 <= 0.05 and all(b > 1 for b in bs) and inc
    say("  Y4  max |deviation| = %.4f  (tol 0.05), all > 1 %s, "
        "increasing %s   %s"
        % (e4, all(b > 1 for b in bs), inc, "hold" if y4 else "REFUTED"))
    e5 = max(abs(a - b) for a, b in zip(b2s, PUB_B2))
    y5 = e5 <= 0.05
    say("  Y5  max |deviation| = %.4f  (tol 0.05)   %s"
        % (e5, "hold" if y5 else "REFUTED"))

    say()
    say("=" * 74)
    ok = y1 and y2 and y3 and y4 and y5 and y6
    say("Y1 %s  Y2 %s  Y3 %s  Y4 %s  Y5 %s  Y6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (y1, y2, y3, y4, y5, y6)))
    say("Proposition {#prop:W} and [eq:budget] reproduce" if ok
        else "REFUTED")

    head = [
        "STATISTIC: Gamma(N) = (sum_{h!=0} c(h)) / V(N) with the numerator",
        "           in the closed form theta(N)^2 - sum_{p<N}(log p)^2;",
        "           Gamma log N / N against 1/A(N); and the absolute budget",
        "           B(X) = (1/V(X)) sum_{0<|h|<X} c(h)|S(h)| with",
        "           S(h) = (1/(X-|h|)) sum mu(u)mu(u-|h|), together with",
        "           the variant normalising S by X.  c(h) and the mu",
        "           autocorrelation are computed by zero-padded rfft, so",
        "           the numerator of Gamma and the budget use independent",
        "           routes.",
        "FIELD: N = 1e4, 1.6e5, 4e6 for Gamma; X = 2e4, 4e4, 8e4, 1.6e5",
        "       for the budget; primes and mu from an integer sieve to",
        "       4e6; Artin's constant as an Euler product over p < 4e6;",
        "       V by direct summation over prime powers w < N.",
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
