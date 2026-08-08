# -*- coding: utf-8 -*-
r"""
Whether the quantitative match in Remark {#rem:directlevel} is what it
looks like.

WHAT IS AT STAKE

That remark reports that the measured crossing K*_H is not merely of
the right exponent but of the right SIZE: solving
K = S^2 N / (4 log^2 K) -- the crossing that |H| ~ (N/k)^{1/2} would
give -- reproduces the measurement to within 23 percent at the bottom
and 4 percent at the top.

Two later measurements say that heuristic omitted two things, each
large.  audit_forecast_null.py measured the mean of |A|/sqrt(N/k),
which is |H|/sqrt(N/k), at 3.5647 to 3.7497 -- so the constant in the
square-root law is not 1 but about sqrt(log N).  And the k in the sum
are squarefree and coprime to N, a density near 0.34, where the
heuristic integrated over every k.  One factor is about 3.7 up, the
other about 0.34 down.

If those two nearly cancel then the reported agreement is the product
of two omissions and not evidence about the square-root law, and
putting them back should not improve the prediction.  That is what is
tested here, by rebuilding both predictions and comparing them with
the measurement.

BACKS: Remark {#rem:heuristic} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  X1  The naive heuristic is correctly reconstructed: solving
      K = S^2 N/(4 log^2 K) reproduces the miss ratios printed in
      results/lab_direct_level.txt to better than 0.01 at every N.
  X2  The two omissions nearly cancel: the product of the measured
      constant c(N) and the measured density of admissible k lies in
      [0.9, 1.5] at every N.
  X3  Putting them back does not help: the corrected prediction's
      miss is not closer to 1 than the naive one at a majority of N.
  X4  And the naive agreement is therefore not evidence: its miss
      drifts monotonically across the sweep rather than settling.

REFUTATION RULE (fixed before the run)

  X1  REFUTED at 0.01 at any N -- it would mean this script is not
      computing what the remark computed.
  X2  REFUTED if the product leaves [0.9, 1.5] at any N.
  X3  REFUTED if the corrected prediction is closer at three or more
      of the five N. That is the good outcome: it would mean the
      corrected heuristic is the right one and the remark's agreement
      is a weaker version of a real match.
  X4  REFUTED if the naive miss is not monotone.

  All four gate.

  NO NULL IS RUN and none applies. Nothing here is a detection against
  a background: two deterministic predictions are compared with one
  measured crossing. The sign controls for K*_H were run in
  lab_direct_level.py, whose mu-squared reference established that the
  level is bought by cancellation at all.
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
OUT = os.path.join(ROOT, "results", "audit_directlevel_heuristic.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 60_000
CLIM = 4_000_000


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


def read_published():
    """the miss ratios lab_direct_level.py prints, read not copied"""
    p = os.path.join(ROOT, "results", "lab_direct_level.txt")
    if not os.path.exists(p):
        return {}
    src = io.open(p, encoding="utf-8").read()
    i = src.find("S^2 N/(4 log^2 K*)")
    if i < 0:
        return {}
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            if out:
                break
            continue
        out[int(f[0])] = (int(f[1]), float(f[3]))
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

    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    pub = read_published()
    say("  read %d published rows from results/lab_direct_level.txt"
        % len(pub))

    rows = []
    for N in NS:
        PN = factor_set(N)
        S = twin
        for q in sorted(PN):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        ks = np.array([k for k in range(2, KCAP)
                       if sqf[k] and all(k % q for q in PN)],
                      dtype=np.int64)
        lw = np.log(ks.astype(np.float64))
        f0 = np.zeros(N, dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        f0[1:] = lam[1:N] * mu[N - idx].astype(np.float64)
        A = np.empty(ks.size)
        for i, k in enumerate(ks):
            r = N % int(k)
            A[i] = f0[r::int(k)].sum() if r else f0[int(k)::int(k)].sum()
        del f0
        aH = np.abs(A)
        cum = np.cumsum(lw * aH)
        thr = S * N
        j = int(np.searchsorted(cum, thr))
        kstar = int(ks[min(j, ks.size - 1)])
        sel = ks <= kstar
        c = float((aH[sel] / np.sqrt(N / ks[sel].astype(float))).mean())
        dens = float(sel.sum()) / kstar
        rows.append((N, S, ks, lw, kstar, c, dens))
        say("  N = %-10d  K* = %-8d c(N) = %.4f  density = %.4f"
            % (N, kstar, c, dens))

    say()
    say("X1  the naive heuristic, K = S^2 N / (4 log^2 K)")
    say("  N            K*        naive K   ratio     published")
    x1 = True
    naive = []
    for N, S, ks, lw, kstar, c, dens in rows:
        pk = S * S * N / (4.0 * math.log(kstar) ** 2)
        r = kstar / pk
        naive.append(r)
        p = pub.get(N, (0, float("nan")))[1]
        if not (abs(r - p) < 0.01):
            x1 = False
        say("  %-12d %-9d %-9.0f %-9.4f %.4f" % (N, kstar, pk, r, p))
    say("  X1 %s" % ("hold" if x1 else "REFUTED"))

    say()
    say("X2  the two omitted factors")
    say("  N            c(N)      density   product")
    x2 = True
    for N, S, ks, lw, kstar, c, dens in rows:
        pr2 = c * dens
        if not (0.9 <= pr2 <= 1.5):
            x2 = False
        say("  %-12d %-9.4f %-9.4f %.4f" % (N, c, dens, pr2))
    say("  X2 %s" % ("hold" if x2 else "REFUTED"))

    say()
    say("X3  putting them back: solve sum_{k<K} (log k) c sqrt(N/k) = S N")
    say("  over the ACTUAL admissible k")
    say("  N            corrected K   ratio     naive ratio   closer?")
    better = 0
    for i, (N, S, ks, lw, kstar, c, dens) in enumerate(rows):
        model = np.cumsum(lw * c * np.sqrt(N / ks.astype(float)))
        j = int(np.searchsorted(model, S * N))
        ck = int(ks[min(j, ks.size - 1)])
        r = kstar / ck
        good = abs(r - 1.0) < abs(naive[i] - 1.0)
        better += int(good)
        say("  %-12d %-13d %-9.4f %-13.4f %s"
            % (N, ck, r, naive[i], "yes" if good else "no"))
    x3 = better < 3
    say("  X3 corrected is closer at %d of %d   (cap 2)   %s"
        % (better, len(rows), "hold" if x3 else "REFUTED"))

    say()
    x4 = (all(naive[i] > naive[i + 1] for i in range(len(naive) - 1))
          or all(naive[i] < naive[i + 1] for i in range(len(naive) - 1)))
    say("X4  the naive miss drifts: %s   %s"
        % (", ".join("%.4f" % v for v in naive),
           "hold" if x4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). What the naive form was really doing.")
    say("  It replaced |H| by sqrt(N/k) with constant 1 and summed over")
    say("  every k; the truth is a constant near sqrt(log N) over a")
    say("  k-set of density about a third. Their product against 1:")
    say("  N            c(N)/sqrt(log N)   c(N)*density   1/(c*density)")
    for N, S, ks, lw, kstar, c, dens in rows:
        say("  %-12d %-18.4f %-14.4f %.4f"
            % (N, c / math.sqrt(math.log(N)), c * dens, 1.0 / (c * dens)))

    say()
    say("=" * 70)
    ok = x1 and x2 and x3 and x4
    say("the quantitative agreement is the product of two omissions "
        "that nearly cancel" if ok else "REFUTED")

    head = [
        "STATISTIC: the measured crossing K*_H; the naive prediction",
        "           K = S^2 N/(4 log^2 K) and its miss ratio, against the",
        "           ratios printed by lab_direct_level.py; the measured",
        "           constant c(N) = mean |H|/sqrt(N/k) below K* and the",
        "           density of admissible k; and the prediction from",
        "           sum_{k<K}(log k) c sqrt(N/k) = S(N)N over the actual",
        "           admissible k.",
        "NULL: none is run and none applies. Two deterministic",
        "      predictions are compared with one measured crossing;",
        "      there is no background to detect against. The sign",
        "      controls for K*_H were run in lab_direct_level.py, whose",
        "      mu-squared reference established that the level is bought",
        "      by cancellation at all.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree, coprime",
        "       to N, 2 <= k < 60000; S(N) from an Euler product at the",
        "       fixed bound 4e6; the published ratios are read from",
        "       results/lab_direct_level.txt.",
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
