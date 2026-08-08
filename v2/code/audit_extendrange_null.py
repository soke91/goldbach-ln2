# -*- coding: utf-8 -*-
r"""
The control lab_extend_range.py declined, run.

WHAT IS AT STAKE

lab_extend_range.py backs Remark {#rem:extendrange}, whose content is
an extrapolation: B(N)/N is fitted against N and the fit is followed
out to where it crosses the Goldbach threshold S(N)(1-A(N)), giving
the bracket this program quotes.  Its result file declines a control:
"none is run here and none is needed for what is claimed ... No new
detection is claimed, only a longer lever."

Remark {#rem:weightgapnull} showed that a declined control is worth
running when the statistic is well conditioned under it, and worth
NOT running when it is not.  B(N)/N is a sum of absolute values, so a
coin cannot send it near zero: it is exactly the well-conditioned
case.  And the claim is not "B is small" but "B decays at such a rate
that it crosses at 10^{8.4}", which is a statement about a fitted
decay -- precisely what a control can speak to.

BACKS: Remark {#rem:extendnull} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  Direction, as {#rem:whycoinwins} requires: the coin's B/N is
      below mu's at every N and every draw, because the coin gets
      square-root cancellation inside the progression sum and mu does
      not.
  W2  A good fit is not evidence: the coin's own log-law fit of B/N
      against log N has correlation past 0.98 in magnitude at every
      draw, just as mu's does. So the quality of the fit says nothing
      about mu.
  W3  The rate is mu's: mu's fitted log-law exponent lies outside the
      range spanned by the 8 draws.
  W4  And therefore so is the crossing: mu's crossing point, in
      log10 N, lies outside the range spanned by the 8 draws.

REFUTATION RULE (fixed before the run)

  W1  REFUTED if any draw reaches mu's B/N at any N. A failure would
      contradict [rem:whycoinwins] and mean the coin construction is
      wrong.
  W2  REFUTED if any draw's correlation is under 0.98 in magnitude --
      which would be the good outcome, since it would mean fit quality
      IS diagnostic here.
  W3  REFUTED if mu's exponent falls inside the draws' range. That is
      the one that decides whether the extrapolated rate is mu's.
  W4  REFUTED if mu's crossing falls inside the draws' range.

  All four gate.

  NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere, with
  the field, the weight log k, the k-range and theta' identical, so
  the sign pattern on the long variable is the only difference. Eight
  draws.
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
OUT = os.path.join(ROOT, "results", "audit_extendrange_null.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
CLIM = 4_000_000
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


def phi_of(k):
    v, phi, d = k, 1, 2
    while d * d <= v:
        if v % d == 0:
            phi *= (d - 1)
            v //= d
            while v % d == 0:
                phi *= d
                v //= d
        d += 1
    if v > 1:
        phi *= (v - 1)
    return phi


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


def crossing(c, thr):
    """where the log-law fit exp(c1 + c0 log log N) meets thr"""
    lo, hi = math.log(1e5), 900.0
    for _ in range(400):
        mid = 0.5 * (lo + hi)
        if math.exp(c[1] + c[0] * math.log(mid)) > thr:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi) / math.log(10)


def loo(x, y, name, say):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    f = [float(np.polyfit(x[t], y[t], 1)[0])
         for t in (slice(None), slice(1, None), slice(0, -1))]
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

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    rng = np.random.default_rng(SEED)
    sup = mu != 0
    signs = [mu.astype(np.float64)]
    for _ in range(COINS):
        c = np.zeros(NMAX + 1, dtype=np.float64)
        c[sup] = rng.choice([-1.0, 1.0], size=int(sup.sum()))
        signs.append(c)
    say("  1 mu and %d coins on supp(mu^2), field and weights identical"
        % COINS)

    B = [[] for _ in signs]
    thrs = []
    for N in NS:
        PN = factor_set(N)
        A_, S_ = artin, twin
        for q in sorted(PN):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))
        thrs.append(S_ * (1.0 - A_))
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if mu[k] != 0 and all(k % q for q in PN)])
        lg = np.log(ks.astype(float))
        iph = np.array([phi_of(int(k)) for k in ks], dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        for j, sg in enumerate(signs):
            f0 = np.zeros(N, dtype=np.float64)
            f0[1:] = lam[1:N] * sg[N - idx]
            C = float(f0.sum())
            A = np.empty(ks.size)
            for i, k in enumerate(ks):
                r = N % int(k)
                A[i] = (f0[r::int(k)].sum() if r
                        else f0[int(k)::int(k)].sum())
            del f0
            B[j].append(float((lg * np.abs(A - C / iph)).sum()) / N)
        say("  N = %-10d  K = %-7d #k = %d  threshold %.6f"
            % (N, K, ks.size, thrs[-1]))

    say()
    say("W1  B(N)/N, mu against the draws")
    say("  N            mu        coin min   coin median   coin max")
    w1 = True
    for i, N in enumerate(NS):
        c = [B[j][i] for j in range(1, len(signs))]
        if max(c) >= B[0][i]:
            w1 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, B[0][i], min(c), float(np.median(c)), max(c)))
    say("  W1 %s" % ("hold" if w1 else "REFUTED"))

    say()
    LL = np.log(np.log(np.array(NS, dtype=float)))
    fits, cors, cross = [], [], []
    for j in range(len(signs)):
        y = np.log(np.array(B[j]))
        c = np.polyfit(LL, y, 1)
        fits.append(float(c[0]))
        cors.append(float(np.corrcoef(LL, y)[0, 1]))
        cross.append(crossing(c, thrs[-1]))
    say("W2  the log-law fit, B/N ~ (log N)^e")
    say("  who       exponent   correlation   crossing at 10^x")
    say("  %-9s %-10.4f %-13.5f %.2f" % ("mu", fits[0], cors[0], cross[0]))
    w2 = True
    for j in range(1, len(signs)):
        if abs(cors[j]) <= 0.98:
            w2 = False
        say("  coin %-4d %-10.4f %-13.5f %.2f"
            % (j - 1, fits[j], cors[j], cross[j]))
    say("  W2 every draw fits past 0.98 in magnitude   %s"
        % ("hold" if w2 else "REFUTED"))

    say()
    ce = fits[1:]
    w3 = not (min(ce) <= fits[0] <= max(ce))
    say("W3  mu's exponent %.4f against the draws' range [%.4f, %.4f]"
        % (fits[0], min(ce), max(ce)))
    say("    %s" % ("hold" if w3 else "REFUTED"))
    cc = cross[1:]
    w4 = not (min(cc) <= cross[0] <= max(cc))
    say("W4  mu's crossing 10^%.2f against the draws' range "
        "[10^%.2f, 10^%.2f]" % (cross[0], min(cc), max(cc)))
    say("    %s" % ("hold" if w4 else "REFUTED"))

    say()
    loo(LL, np.log(np.array(B[0])), "extendrange_Bdecay", say)

    say()
    say("  DIAGNOSTIC (post hoc). Why this control is usable where the")
    say("  one in [rem:weightgapnull] was not: B/N is a sum of absolute")
    say("  values, so a coin cannot send it near zero. The draws'")
    say("  spread as a fraction of their median:")
    say("  N            median    spread/median")
    for i, N in enumerate(NS):
        c = [B[j][i] for j in range(1, len(signs))]
        m = float(np.median(c))
        say("  %-12d %-9.4f %.4f" % (N, m, (max(c) - min(c)) / m))
    say("  A well-conditioned statistic under the control is the")
    say("  condition for the control to mean anything, and it is what")
    say("  separates this audit from the ratio in [rem:weightgapnull].")

    say()
    say("=" * 70)
    ok = w1 and w2 and w3 and w4
    say("the extrapolated rate and crossing are mu's, and the quality "
        "of the fit is not evidence" if ok else "REFUTED")

    head = [
        "STATISTIC: B(N)/N = sum_{k<K}(log k)|E(N;k)|/N for mu and for 8",
        "           coins; the log-law fit of each against log N with its",
        "           correlation; and the N at which each fit crosses the",
        "           Goldbach threshold S(N)(1-A(N)).",
        "NULL: the coin, eps(v) = +-1 on supp(mu^2) and zero elsewhere,",
        "      with the field, the weight log k, the k-range and theta'",
        "      identical. This is the control lab_extend_range.py",
        "      declined. It is usable here, unlike the ratio of",
        "      [rem:weightgapnull], because B/N is a sum of absolute",
        "      values and a coin cannot drive it near zero.",
        "FIELD: N = 2e5 through 6.4e6 by doubling, theta' = 0.56; k over",
        "       the squarefree k < N^0.56 coprime to N; S and A from",
        "       Euler products at the fixed bound 4e6; seed 20260808.",
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
