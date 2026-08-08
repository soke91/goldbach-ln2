# -*- coding: utf-8 -*-
r"""
The control lab_decay_family.py points at, against what it actually
does.

WHAT IS AT STAKE

lab_decay_family.py declines a control by pointing elsewhere: "the
field and its coin reference were measured in lab_sign_structure.py
and lab_lean_decay.py, and this script only fits a functional form to
numbers already controlled there."

Fitting a functional form is not nothing.  The script sweeps alpha in
log|1/2 - f| = a - c (log N)^alpha, reports the minimiser, the span of
alpha within one percent of it, and the range of N over which three
canonical laws put |1/2 - f| = 0.01 -- nine orders.  None of that is
controlled by having measured f itself: the question is whether the
SWEEP can discriminate, and eight noisy points will produce a
minimiser, a band and a wild extrapolation whatever they contain.

For a coin f is 1/2 in expectation, so |1/2 - f| is a near-zero
fluctuation and, by the criterion of {#rem:weightgapnull}, its
MAGNITUDE is ill conditioned under the control.  But the claims here
are not about magnitude: they are about whether a fit exists, how
tight its band is, and how far it extrapolates.  Those a coin can
answer, and this script asks it.

BACKS: Remark {#rem:decaynull} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  The lean itself is mu's: |1/2 - f| for mu exceeds every one of
      8 draws at every N.
  Y2  And its decay is: mu's |1/2 - f| falls across the sweep, and at
      most 2 of the 8 draws fall monotonically.
  Y3  The sweep does not discriminate: at least 6 of the 8 draws also
      produce an interior minimiser of the residual sum of squares in
      alpha, so "there is a minimiser" is not evidence.
  Y4  Nor does the band: the median draw's one-percent band in alpha
      is at least half mu's, so the band's width is not evidence
      either -- which is what lab_decay_family.py's own diagnostic
      says about its threshold, now measured against a control.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED if any draw reaches mu at any N. This would mean the
      coin construction is wrong, since a coin has no lean.
  Y2  REFUTED if 3 or more draws fall monotonically.
  Y3  REFUTED if fewer than 6 draws have an interior minimiser -- the
      good outcome, since it would mean the sweep does discriminate.
  Y4  REFUTED if the median band is under half mu's, likewise good.

  All four gate.

  NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere, with
  the field, the weight log k, the k-range and theta' identical.
  Eight draws.
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
OUT = os.path.join(ROOT, "results", "audit_decayfamily_null.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
ALPHAS = [round(0.05 * i, 2) for i in range(1, 31)]
COINS = 8
SEED = 20260808


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


def sweep(dev):
    """RSS of log dev = a - c (log N)^alpha over the alpha grid."""
    y = np.log(np.array(dev))
    L = np.log(np.array(NS, dtype=float))
    rss = []
    for a in ALPHAS:
        x = L ** a
        c = np.polyfit(x, y, 1)
        rss.append(float(((y - (c[0] * x + c[1])) ** 2).sum()))
    rss = np.array(rss)
    j = int(np.argmin(rss))
    interior = 0 < j < len(ALPHAS) - 1
    keep = rss <= rss[j] * 1.01
    band = (ALPHAS[int(np.flatnonzero(keep)[-1])]
            - ALPHAS[int(np.flatnonzero(keep)[0])])
    return ALPHAS[j], band, interior, float(rss[j])


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0

    rng = np.random.default_rng(SEED)
    signs = [mu.astype(np.float64)]
    for _ in range(COINS):
        c = np.zeros(NMAX + 1, dtype=np.float64)
        c[sqf] = rng.choice([-1.0, 1.0], size=int(sqf.sum()))
        signs.append(c)
    say("  1 mu and %d coins on supp(mu^2), field and weights identical"
        % COINS)

    dev = [[] for _ in signs]
    for N in NS:
        PN = factor_set(N)
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if sqf[k] and all(k % q for q in PN)])
        lg = np.log(ks.astype(float))
        idx = np.arange(1, N, dtype=np.int64)
        for j, sg in enumerate(signs):
            f0 = np.zeros(N, dtype=np.float64)
            f0[1:] = lam[1:N] * sg[N - idx]
            A = np.empty(ks.size)
            for i, k in enumerate(ks):
                r = N % int(k)
                A[i] = (f0[r::int(k)].sum() if r
                        else f0[int(k)::int(k)].sum())
            del f0
            H = sg[ks] * A
            w = lg * np.abs(H)
            fr = float(w[H > 0].sum() / w.sum())
            dev[j].append(max(abs(0.5 - fr), 1e-12))
        say("  N = %-10d  #k = %d" % (N, ks.size))

    say()
    say("Y1/Y2  |1/2 - f| for mu and the draws")
    say("  N            mu        coin min   coin median   coin max")
    y1 = True
    for i, N in enumerate(NS):
        c = [dev[j][i] for j in range(1, len(signs))]
        if max(c) >= dev[0][i]:
            y1 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, dev[0][i], min(c), float(np.median(c)), max(c)))
    say("  Y1 mu exceeds every draw everywhere   %s"
        % ("hold" if y1 else "REFUTED"))

    def falling(v):
        return all(v[i] > v[i + 1] for i in range(len(v) - 1))

    nfall = sum(falling(dev[j]) for j in range(1, len(signs)))
    y2 = falling(dev[0]) and nfall <= 2
    say("  Y2 mu falls monotonically: %s; draws that do: %d of %d   %s"
        % (falling(dev[0]), nfall, COINS, "hold" if y2 else "REFUTED"))

    say()
    say("Y3/Y4  what the alpha sweep does on each series")
    say("  who       alpha*    1%% band   interior   RSS at the minimum")
    am, bm, im, rm = sweep(dev[0])
    say("  %-9s %-9.2f %-9.2f %-10s %.6f" % ("mu", am, bm, im, rm))
    bands, nint = [], 0
    for j in range(1, len(signs)):
        a, b, it, r = sweep(dev[j])
        bands.append(b)
        nint += int(it)
        say("  coin %-4d %-9.2f %-9.2f %-10s %.6f"
            % (j - 1, a, b, it, r))
    y3 = nint >= 6
    say("  Y3 draws with an interior minimiser: %d of %d   %s"
        % (nint, COINS, "hold" if y3 else "REFUTED"))
    med = float(np.median(bands))
    y4 = med >= 0.5 * bm
    say("  Y4 median draw band %.2f against mu's %.2f, ratio %.2f   %s"
        % (med, bm, med / bm if bm else float("nan"),
           "hold" if y4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). What the pointer covers and what it")
    say("  does not. The control it names measured f itself, and Y1")
    say("  confirms that measurement: mu's lean is real and no draw")
    say("  comes near it. What it does not cover is the sweep, and Y3")
    say("  and Y4 say the sweep is uninformative on its own -- a series")
    say("  with no lean at all still yields a minimiser and a band of")
    say("  comparable width. So lab_decay_family.py's own remark, that")
    say("  its one-percent band is a threshold and not a null, is")
    say("  correct and is now measured rather than conceded.")

    say()
    say("=" * 70)
    ok = y1 and y2 and y3 and y4
    say("the lean is mu's and the alpha sweep is not evidence about it"
        if ok else "REFUTED")

    head = [
        "STATISTIC: |1/2 - f(N)| with f the mass-weighted fraction of k",
        "           with H(N;k) > 0, for mu and 8 coins; whether each",
        "           series falls monotonically; and for each series the",
        "           minimiser alpha* of the residual sum of squares of",
        "           log|1/2-f| = a - c (log N)^alpha, the span of alpha",
        "           within one percent of that minimum, and whether the",
        "           minimiser is interior to the grid.",
        "NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere,",
        "      field, weight log k, k-range and theta' identical. This is",
        "      the control lab_decay_family.py declined by pointing at",
        "      lab_sign_structure.py and lab_lean_decay.py, which measure",
        "      f but not the alpha sweep that this script's claims are",
        "      about.",
        "FIELD: N = 2e5 through 6.4e6 by doubling, theta' = 0.56; k over",
        "       the squarefree k < N^0.56 coprime to N; alpha swept over",
        "       0.05 to 1.50 in steps of 0.05; seed 20260808.",
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
