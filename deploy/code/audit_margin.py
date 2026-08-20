# -*- coding: utf-8 -*-
r"""
paper/wall_v3.md, Section "The margin, and where the difficulty is".

WHAT IS UNDER TEST

The section quotes the requirement C(N) = o(N) at the extreme rather
than at a typical N, and prints:

  (a) max |C|/N falls from 0.114 in the octave at 3e4 to 0.0101 in the
      octave at 1.6e7, so the margin at the top of the computed range is
      a factor 99;
  (b) the fitted decay is max|C|/N ~ N^{-0.43}, which extrapolates to
      99 * 6.25^{0.43}, a factor near 220, at N = 10^8;
  (c) with max|C| ~ a_n sqrt(V(N)) and a_n the Gumbel location for the
      number of even N below the point,
        N / max_{N<=X} |C(N)| ~ sqrt(N) / (a_n sqrt(A log N)),
      which is 10^{4.5} at N = 10^{12} and 10^{22.9} at N = 10^{50};
  (d) Cauchy-Schwarz gives N sqrt(6(log N - 1)/pi^2), a factor
      (log N)^{A+1/2} above the target N (log N)^{-A}.

No script for any of it exists here.

CONVENTIONS FIXED HERE

"The octave at X" is read as (X/2, X], which is the reading
Section {#sec:floor} uses for its top octave (8e6, 1.6e7].  The Gumbel
location is read as a_n = sqrt(2 log n) with n the number of even N
below the point, and A(N) at the value 0.787275 that P(N) = {2,5}
gives -- both are stated because neither is fixed by the text.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  M1  max|C|/N = 0.114 on the octave (1.5e4, 3e4] and 0.0101 on the
      octave (8e6, 1.6e7]; the ratio is 99.
  M2  Fitting max|C|/N ~ N^{-b} across the octaves from (1.5e4, 3e4] to
      (8e6, 1.6e7] gives b = 0.43.
  M3  99 * 6.25^{0.43} lies in [215, 225].
  M4  The Gumbel formula gives 10^{4.5} at N = 10^{12} and 10^{22.9} at
      N = 10^{50}.
  M5  The Cauchy-Schwarz bound equals N sqrt(6(log N - 1)/pi^2): at
      N = 2^20 that is 2.7963 N, matching the 2.795 that
      Proposition {#prop:E}'s table reports.

REFUTATION RULE (fixed before the run)

  M1  REFUTED if the first differs by more than 0.0005, the second by
      more than 0.00005, or the ratio by more than 0.5.
  M2  REFUTED if b differs from 0.43 by more than 0.005.
  M3  REFUTED if the product falls outside [215, 225].
  M4  REFUTED if either exponent differs by more than 0.1.
  M5  REFUTED if the closed form differs from the measured
      ||S_Lambda||_2 ||S_mu||_2 / N at N = 2^20 by more than 0.005.

  All five gate.

CITED BY: {#rem:grid} in paper/.
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
OUT = os.path.join(ROOT, "results", "audit_margin.txt")

X = 16_000_000
A_TYP = 0.787275


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def pow2(n):
    L = 1
    while L < n:
        L <<= 1
    return L


def loo(x, y, name, say):
    """Refit dropping each end in turn, and report the spread.

    Every exponent this repository quotes is a slope over four to eight
    values of N, and audit_truncation_exponent.py showed what such a
    slope is worth when nobody varies the free parameter that defines
    it. For a direct fit the free parameter is the N-range, so the
    cheapest honest check is to refit without the smallest N and
    without the largest and print how far the answer moves.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    sp = max(f) - min(f)
    say("  leave-one-out on %s: full %.4f, without the smallest N "
        "%.4f," % (name, f[0], f[1]))
    say("  without the largest %.4f -- spread %.4f" % (f[2], sp))
    say("SWEPT %s N-range %.4f" % (name, sp))
    return sp


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("sieving to %d ..." % X)
    pr = primes_upto(X)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(X + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > X:
            break
        q = p * p
        while q <= X:
            lam[q] = lgp[i]
            if q > X // p:
                break
            q *= p

    mu = np.ones(X + 1, dtype=np.int8)
    rem = np.arange(X + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(X))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= X:
            mu[p * p::p * p] = 0
        q = p
        while q <= X:
            rem[q::q] //= p
            if q > X // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0

    say("convolving C = mu * Lambda ...")
    n = pow2(2 * (X + 1))
    a = np.zeros(n, dtype=np.float64)
    a[:X + 1] = lam
    F = np.fft.rfft(a)
    a[:] = 0.0
    a[:X + 1] = mu
    C = np.fft.irfft(F * np.fft.rfft(a), n)[:X + 1]
    del a, F

    say()
    say("M1/M2   max |C|/N by octave, read as (X/2, X]")
    say("  octave top     count      max |C|/N     argmax N")
    say("  " + "-" * 56)
    tops, vals = [], []
    t = 30_000
    while t <= X:
        lo = t // 2
        Ns = np.arange(lo + 2 - (lo % 2), t + 1, 2, dtype=np.int64)
        r = np.abs(C[Ns]) / Ns
        j = int(np.argmax(r))
        tops.append(t)
        vals.append(float(r[j]))
        say("  %-14d %-10d %-13.6f %d" % (t, Ns.size, r[j], int(Ns[j])))
        t *= 2

    first, last = vals[0], vals[-1]
    ratio = 1.0 / last
    m1 = (abs(first - 0.114) <= 5e-4 and abs(last - 0.0101) <= 5e-5
          and abs(ratio - 99.0) <= 0.5)
    say("  M1  first %.6f (pub 0.114), last %.6f (pub 0.0101), "
        "1/last = %.3f (pub 99)   %s"
        % (first, last, ratio, "hold" if m1 else "REFUTED"))

    xs = np.log(np.array(tops, dtype=float))
    ys = np.log(np.array(vals, dtype=float))
    b = -float(np.polyfit(xs, ys, 1)[0])
    m2 = abs(b - 0.43) <= 5e-3
    say("  M2  fitted max|C|/N ~ N^{-%.6f}   published 0.43   %s"
        % (b, "hold" if m2 else "REFUTED"))
    loo(xs, ys, "wall_max_decay", say)

    say("  DIAGNOSTIC (post hoc). The top octave matches but the bottom")
    say("  does not, so the octave grid's anchor is the free choice. The")
    say("  grid above is anchored at 3e4 and doubles; \\S{#sec:floor} uses")
    say("  'eight octaves from 6.25e4 to 1.6e7', i.e. a grid anchored at")
    say("  1.6e7 and halving. On that grid:")
    say("  octave top     count      max |C|/N     argmax N")
    tops2, vals2 = [], []
    t = X
    while t >= 30_000:
        lo = t // 2
        Ns = np.arange(lo + 2 - (lo % 2), t + 1, 2, dtype=np.int64)
        r = np.abs(C[Ns]) / Ns
        j = int(np.argmax(r))
        tops2.append(t)
        vals2.append(float(r[j]))
        say("  %-14d %-10d %-13.6f %d" % (t, Ns.size, r[j], int(Ns[j])))
        t //= 2
    tops2, vals2 = tops2[::-1], vals2[::-1]
    b2 = -float(np.polyfit(np.log(np.array(tops2, dtype=float)),
                           np.log(np.array(vals2, dtype=float)), 1)[0])
    say("  bottom %.6f (pub 0.114), top %.6f (pub 0.0101), "
        "1/top = %.3f (pub 99), fitted b = %.6f (pub 0.43)"
        % (vals2[0], vals2[-1], 1.0 / vals2[-1], b2))

    ext = ratio * 6.25 ** 0.43
    m3 = 215.0 <= ext <= 225.0
    say()
    say("M3   extrapolation to N = 10^8: %.3f * 6.25^0.43 = %.2f   "
        "published 'near 220'   %s"
        % (ratio, ext, "hold" if m3 else "REFUTED"))

    say()
    say("M4   the Gumbel form  sqrt(N) / (a_n sqrt(A log N))")
    say("  N          a_n        sqrt(A log N)   value        log10   pub")
    m4 = True
    for Nv, pub in ((1e12, 4.5), (1e50, 22.9)):
        nn = Nv / 2.0
        an = math.sqrt(2.0 * math.log(nn))
        s = math.sqrt(A_TYP * math.log(Nv))
        v = math.sqrt(Nv) / (an * s)
        lg = math.log10(v)
        ok = abs(lg - pub) <= 0.1
        m4 = m4 and ok
        say("  10^%-8d %-10.4f %-15.4f %-12.4e %-7.3f %.1f  %s"
            % (int(round(math.log10(Nv))), an, s, v, lg, pub,
               "ok" if ok else "MISMATCH"))
    say("  M4 %s" % ("hold" if m4 else "REFUTED"))

    say()
    say("M5   the Cauchy-Schwarz closed form at N = 2^20")
    Nv = 1 << 20
    closed = math.sqrt(6.0 * (math.log(Nv) - 1.0) / math.pi ** 2)
    meas = (math.sqrt(float((lam[1:Nv + 1] ** 2).sum()))
            * math.sqrt(float((mu[1:Nv + 1].astype(np.float64) ** 2).sum()))
            / Nv)
    m5 = abs(closed - meas) <= 5e-3
    say("  closed form N sqrt(6(log N - 1)/pi^2) : %.6f N" % closed)
    say("  measured   ||S_Lambda||_2 ||S_mu||_2  : %.6f N   "
        "(prop:E table prints 2.795)" % meas)
    say("  M5 %s" % ("hold" if m5 else "REFUTED"))

    say()
    say("=" * 70)
    ok = m1 and m2 and m3 and m4 and m5
    say("M1 %s  M2 %s  M3 %s  M4 %s  M5 %s"
        % tuple("hold" if v else "REFUTED" for v in (m1, m2, m3, m4, m5)))
    say("the margin section reproduces" if ok else "REFUTED")

    head = [
        "STATISTIC: max |C(N)|/N over even N in each octave (T/2, T], with",
        "           C = mu * Lambda by exact FFT convolution; the exponent",
        "           b in a least-squares fit of log max|C|/N against",
        "           log T; the extrapolation 99 * 6.25^0.43; the Gumbel",
        "           form sqrt(N)/(a_n sqrt(A log N)) with",
        "           a_n = sqrt(2 log(N/2)) and A = 0.787275; and the",
        "           Cauchy-Schwarz closed form sqrt(6(log N - 1)/pi^2)",
        "           against ||S_Lambda||_2 ||S_mu||_2 / N.",
        "FIELD: even N; octaves with tops 3e4, 6e4, ... up to 1.6e7;",
        "       Lambda and mu from an integer sieve to 1.6e7; the Gumbel",
        "       evaluations at N = 10^12 and 10^50; the Cauchy-Schwarz",
        "       check at N = 2^20.",
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
