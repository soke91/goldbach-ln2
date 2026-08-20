# -*- coding: utf-8 -*-
r"""
paper/wall_v3.md, Proposition {#prop:V} and the paragraph
"The local factor is A, not S" under it.

WHAT IS UNDER TEST

    V(N) = sum_{v<N} mu^2(v) Lambda(N-v)^2,   W(N) = sum_{w<N} Lambda(w)^2,
    A(N) = prod_{q not| N} (1 - 1/(q(q-1))),

and the proposition asserts V(N) = W(N) A(N) (1+o(1)).  The paragraph
under it prints five figures and names no script:

    (a) the residual standard deviation over even N in [1e5, 1.6e7],
        after rescaling each candidate to the measured mean, is
        0.000323 for A against 0.245235 for S -- a factor of 759;
    (b) V(N)/W(N) against A(N) has mean 1.000000 and standard
        deviation 0.000166 in the top octave;
    (c) it agrees to six decimals in each of the six radical cells from
        2|N through 2*3*5*7*11*13 | N;
    (d) at N = 4e6 the ratio W/V is 1.27080 against 1/A(N) = 1.27020;
    (e) taking the lower cutoff below 1e5 degrades the figure for A to
        0.000582.

THE STATISTIC IS NOT FIXED BY THE TEXT, SO IT IS FIXED HERE

"Rescaling each candidate to the measured mean, so that only its shape
in N is judged" admits at least two readings.  Both are computed:

    reading R1 (primary): z_c(N) = (V/W)(N) / c(N); the reported figure
        is sd(z_c / mean(z_c)).
    reading R2:           choose alpha with mean(alpha c) = mean(V/W);
        the reported figure is sd(V/W - alpha c) / mean(V/W).

R1 is primary because (b) reports "mean 1.000000, sd 0.000166", which
is R1's shape with c = A, where the mean is already 1.

Cell convention, likewise not fixed by the text and fixed here: cell j
is the set of even N divisible by the product of 2 and the first j-1 of
3,5,7,11,13, for j = 1..6.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  X1  1/A(4e6) = 1.27020.
  X2  W(4e6)/V(4e6) = 1.27080.
  X3  Over the top octave (8e6, 1.6e7], V/W divided by A has mean
      1.000000 and standard deviation 0.000166.
  X4  Under reading R1 over even N in [1e5, 1.6e7], the figure is
      0.000323 for A and 0.245235 for S, ratio 759.
  X5  Each of the six radical cells has |mean(V/(W A)) - 1| < 5e-7,
      i.e. agrees to six decimals.
  X6  Some lower cutoff in {2, 1e3, 1e4} reproduces 0.000582 for A.
  X7  V(N)/(A(N) N log N) -> 1: it lies within 5% of 1 at N = 1.6e7.

REFUTATION RULE (fixed before the run)

  X1, X2  REFUTED if the recomputed value differs by more than 0.000005.
  X3      REFUTED if |mean - 1| > 5e-7 or |sd - 0.000166| > 5e-6.
  X4      REFUTED if either figure differs by more than 0.00002 under
          reading R1, or if the ratio differs from 759 by more than 5.
  X5      REFUTED if any of the six cells has |mean - 1| >= 5e-7.
  X6      REFUTED if no cutoff in {2, 1e3, 1e4} lands within 0.00001 of
          0.000582.
  X7      REFUTED if the ratio is outside [0.95, 1.05] at N = 1.6e7.

  X1, X2, X3, X7 gate: they are the proposition and its two sharp
  numbers.  X4, X5, X6 are reported but do not gate, because a
  disagreement there is evidence that the statistic is under-specified
  in the source text rather than that the proposition is wrong -- and
  reading R2 is computed alongside so the disagreement is diagnosable.

CITED BY: {#rem:secondorder} in paper/.
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
OUT = os.path.join(ROOT, "results", "lab_second_moment.txt")

X = 16_000_000
LOW = 100_000
NCHK = 4_000_000
ARTIN_LIM = 10_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("sieving to %d ..." % X)
    pr = primes_upto(X)
    lgp = np.log(pr.astype(np.float64))

    # mu^2 : squarefree indicator
    sqf = np.ones(X, dtype=np.float64)
    sqf[0] = 0.0
    for i, p in enumerate(pr):
        p = int(p)
        if p * p >= X:
            break
        sqf[p * p::p * p] = 0.0

    # Lambda^2
    lam2 = np.zeros(X, dtype=np.float64)
    lam2[pr[pr < X]] = lgp[pr < X] ** 2
    for i, p in enumerate(pr):
        p = int(p)
        if p * p >= X:
            break
        q = p * p
        while q < X:
            lam2[q] = lgp[i] ** 2
            if q > (X - 1) // p:
                break
            q *= p

    W = np.cumsum(lam2)                    # W(N) = sum_{w<N} = W[N-1]

    say("convolving (length 2^25 rfft) ...")
    n = 1 << 25
    a = np.zeros(n, dtype=np.float64)
    a[:X] = lam2
    del lam2
    F = np.fft.rfft(a)
    del a
    b = np.zeros(n, dtype=np.float64)
    b[:X] = sqf
    del sqf
    F *= np.fft.rfft(b)
    del b
    V = np.fft.irfft(F, n)[:X + 1]
    del F

    # A(N) = Artin * prod_{q|N} (1 - 1/(q(q-1)))^{-1}
    say("building A(N) ...")
    artin = 1.0
    for p in primes_upto(ARTIN_LIM):
        artin *= 1.0 - 1.0 / (int(p) * (int(p) - 1.0))
    AN = np.full(X, artin, dtype=np.float64)
    twin = 2.0
    for p in pr[1:]:
        twin *= 1.0 - 1.0 / (int(p) - 1.0) ** 2
    SN = np.full(X, twin, dtype=np.float64)
    for p in pr:
        p = int(p)
        AN[p::p] /= (1.0 - 1.0 / (p * (p - 1.0)))
        if p > 2:
            SN[p::p] *= (1.0 + 1.0 / (p - 2.0))

    ev = np.arange(2, X, 2, dtype=np.int64)
    y = V[ev] / W[ev - 1]                  # V(N)/W(N)

    def r1(c):
        z = y / c[ev]
        return float(np.std(z / z.mean())), float(z.mean())

    def r2(c):
        alpha = y.mean() / c[ev].mean()
        return float(np.std(y - alpha * c[ev]) / y.mean())

    say()
    say("X1/X2  the single point N = %d" % NCHK)
    say("=" * 70)
    invA = 1.0 / AN[NCHK]
    WV = W[NCHK - 1] / V[NCHK]
    say("  1/A(N)  = %.6f   published 1.27020" % invA)
    say("  W/V     = %.6f   published 1.27080" % WV)
    x1 = abs(invA - 1.27020) <= 5e-6
    x2 = abs(WV - 1.27080) <= 5e-6
    say("  X1 %s   X2 %s" % ("hold" if x1 else "REFUTED",
                             "hold" if x2 else "REFUTED"))

    say()
    say("X3  the top octave (8e6, 1.6e7]")
    say("=" * 70)
    top = ev[(ev > X // 2)]
    zt = V[top] / W[top - 1] / AN[top]
    mt, st = float(zt.mean()), float(zt.std())
    say("  mean = %.7f   published 1.000000" % mt)
    say("  sd   = %.7f   published 0.000166" % st)
    x3 = abs(mt - 1.0) <= 5e-7 and abs(st - 0.000166) <= 5e-6
    say("  X3 %s" % ("hold" if x3 else "REFUTED"))

    say()
    say("X4  shape residual over even N in [1e5, 1.6e7], both readings")
    say("=" * 70)
    m = ev >= LOW
    evf, ysave = ev, y
    ev, y = ev[m], y[m]
    sA, meanA = r1(AN)
    sS, meanS = r1(SN)
    say("  reading R1 (primary):  A %.6f   S %.6f   ratio %.1f"
        % (sA, sS, sS / sA if sA else float("inf")))
    say("                          published: A 0.000323   S 0.245235"
        "   ratio 759")
    say("  reading R2:            A %.6f   S %.6f   ratio %.1f"
        % (r2(AN), r2(SN), r2(SN) / r2(AN)))
    say("  mean of V/(W A) on this range = %.7f;  of V/(W S) = %.7f"
        % (meanA, meanS))
    x4 = (abs(sA - 0.000323) <= 2e-5 and abs(sS - 0.245235) <= 2e-5
          and abs(sS / sA - 759) <= 5)
    say("  X4 %s" % ("hold" if x4 else "REFUTED"))

    say()
    say("X5  the six radical cells")
    say("=" * 70)
    say("  cell   modulus   count      mean V/(W A)      |mean - 1|")
    x5 = True
    mod = 2
    for j, q in enumerate([1, 3, 5, 7, 11, 13]):
        mod *= q
        sel = ev[ev % mod == 0]
        if sel.size == 0:
            say("  %-6d %-9d (empty)" % (j + 1, mod))
            x5 = False
            continue
        z = V[sel] / W[sel - 1] / AN[sel]
        d = abs(float(z.mean()) - 1.0)
        say("  %-6d %-9d %-10d %.9f      %.2e"
            % (j + 1, mod, sel.size, float(z.mean()), d))
        if d >= 5e-7:
            x5 = False
    say("  X5 %s" % ("hold" if x5 else "REFUTED"))

    say()
    say("X6  what lower cutoff reproduces the published 0.000582 for A?")
    say("=" * 70)
    x6 = False
    for lo in (2, 1_000, 10_000, 50_000):
        mm = evf >= lo
        ev, y = evf[mm], ysave[mm]
        s, _ = r1(AN)
        hit = abs(s - 0.000582) <= 1e-5
        say("  cutoff %-8d -> %.6f %s" % (lo, s, "<- matches" if hit else ""))
        if hit and lo in (2, 1_000, 10_000):
            x6 = True
    say("  X6 %s" % ("hold" if x6 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc, not a pre-registered test): the sweep is")
    say("  monotone in the cutoff, so the cutoff the source used can be")
    say("  read off.")
    for lo in (20, 50, 100, 200, 400, 600, 800):
        mm = evf >= lo
        ev, y = evf[mm], ysave[mm]
        s, _ = r1(AN)
        say("    cutoff %-8d -> %.6f %s"
            % (lo, s, "<- 0.000582" if abs(s - 0.000582) <= 5e-6 else ""))

    say()
    say("X7  V(N) / (A(N) N log N) at the top")
    say("=" * 70)
    for N in (1_000_000, 4_000_000, 16_000_000 - 2):
        r = V[N] / (AN[N] * N * math.log(N))
        say("  N = %-10d  ratio = %.6f" % (N, r))
    Nt = X - 2
    rt = V[Nt] / (AN[Nt] * Nt * math.log(Nt))
    x7 = 0.95 <= rt <= 1.05
    say("  X7 %s" % ("hold" if x7 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc, not a pre-registered test): partial")
    say("  summation from theta(x) ~ x gives sum_{p<=x}(log p)^2 ~ x log x")
    say("  - x, so the second-order form is A(N)(N log N - N), not")
    say("  A(N) N log N.  The band above ignored the -N term, which is a")
    say("  factor 1 - 1/log N = %.6f at the top." % (1 - 1 / math.log(Nt)))
    for N in (1_000_000, 4_000_000, 16_000_000 - 2):
        r = V[N] / (AN[N] * (N * math.log(N) - N))
        say("    N = %-10d  V/(A(N)(N log N - N)) = %.6f" % (N, r))

    say()
    say("=" * 70)
    say("X1 %s  X2 %s  X3 %s  X4 %s  X5 %s  X6 %s  X7 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (x1, x2, x3, x4, x5, x6, x7)))
    gating = x1 and x2 and x3 and x7
    say("Proposition {#prop:V} stands" if gating else "REFUTED")

    head = [
        "STATISTIC: V(N) = sum_{v<N} mu^2(v) Lambda(N-v)^2 by exact FFT",
        "           convolution; W(N) = sum_{w<N} Lambda(w)^2; the ratio",
        "           V/W divided by A(N) and by S(N); the standard",
        "           deviation of that ratio normalised to mean one",
        "           (reading R1) and the residual after a least-mean",
        "           rescaling (reading R2); the same ratio's mean inside",
        "           each of six radical cells; and V/(A N log N).",
        "FIELD: even N; the octave (8e6, 1.6e7] for the top-octave",
        "       figures, [1e5, 1.6e7] for the shape residual, and cutoffs",
        "       2, 1e3, 1e4, 5e4 for the cutoff sweep; cell j is the even",
        "       N divisible by 2 times the first j-1 of 3,5,7,11,13;",
        "       Lambda and the squarefree indicator from a sieve to 1.6e7;",
        "       Artin's constant as an Euler product over p < 1e7.",
        'NULL: none applies. V(N) depends on mu only through mu^2, so a coin',
        '      on the same support gives V byte-identical -- Lemma {#lem:coin}',
        '      itself. Every statistic here is a support statistic by',
        '      construction and is claimed as one.',
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not gating:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
