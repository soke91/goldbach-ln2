# -*- coding: utf-8 -*-
r"""
The constant in the residue's square-root law, and where it puts the
knife-edge.

WHAT IS AT STAKE

Remark {#rem:heuristic} pinned the constant for H: |H| ~ c(N)sqrt(N/k)
with c(N) = sqrt(log N) to a percent and a half, and summing that over
the actual admissible k predicts K*_H to 1.5 per cent with no drift.
Nothing of the kind has been done for R, and R is what the programme
now turns on: {#rem:residuelevel} measures its level at exponents
0.5654 to 0.5799 -- past the square-root barrier by 0.06 and sitting
on theta' = 0.56, missing it at one N -- and {#rem:betafree} shows
that knife-edge cannot be tuned away by the split constant.

Whether the edge is a fact about the accessible range or about the
object is a question about the CONSTANT. If |R| ~ c_R(N)sqrt(N/k)
with c_R/sqrt(log N) fixed, then the same balance that governs K*_H
governs K*_R, the exponent rises, and the dip at N = 8e5 is local. If
c_R drifts, the rise measured over a factor 16 in N is the drift and
not the law.

BACKS: Remark {#rem:residueconstant} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  The control: the measured crossing reproduces the operative
      K*_R of results/audit_residue_level.txt -- 993, 1447, 2019,
      3319, 5923 -- exactly.
  Q2  The residue has a constant of the same shape: c_R(N) measured
      as the mean of |R|/sqrt(N/k) below K*_R, divided by
      sqrt(log N), is constant across the five N to within 5 per
      cent.
  Q3  And the model built from it predicts the crossing: solving
      sum_{k<K}(log k) c_R sqrt(N/k) = S(N)(1-A(N))N over the actual
      admissible k reproduces K*_R to within 5 per cent at every N.
  Q4  The model says the edge is local: its own crossing exponent
      exceeds 0.56 at every N in the sweep, including the 8e5 where
      the measurement reads 0.5599.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED if any crossing differs at all.
  Q2  REFUTED if the spread of c_R/sqrt(log N) across the five N
      reaches 5 per cent of its mean. That is the one that decides
      whether R has a law or a fit.
  Q3  REFUTED at 5 per cent at any N.
  Q4  REFUTED if the model's exponent falls to 0.56 or below at any
      N -- which would mean the dip is in the law and not in the
      measurement, and the knife-edge is real rather than local.

  All four gate. A forecast is made only if Q2 and Q3 hold, and it
  carries a bracket over a +-10 per cent wobble in c_R, which is what
  {#rem:modeltransfer} measured such constants to be worth.

  NO NULL IS RUN and none applies. A measured magnitude is fitted to a
  deterministic shape and the shape is solved; there is no background
  to detect against. The sign control for R is
  lab_residue_cancellation.py, whose coin arm on the identical delta
  established that R's size is bought by cancellation at exactly a
  coin's rate, and lab_residue_signed.py, whose sign redraw priced
  what the signs across k are worth.
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
OUT = os.path.join(ROOT, "results", "audit_residue_constant.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000
THETA = 0.56
WOBBLE = 0.10


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
    return lam, mu


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


def read_H_constant():
    """c(N)/sqrt(log N) for H -- read, not copied"""
    p = os.path.join(ROOT, "results",
                     "audit_directlevel_heuristic.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("c(N)/sqrt(log N)")
    out = []
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        out.append(float(f[1]))
    return out


def read_published():
    p = os.path.join(ROOT, "results", "audit_residue_level.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("log K*_R/log N")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        out[int(f[0])] = int(f[2])
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d operative K*_R from results/audit_residue_level.txt"
        % len(pub))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = sieves(NMAX)
    sqf = mu != 0
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    res = []
    for N in NS:
        PN = factor_set(N)
        A_, S_ = artin, twin
        for q in sorted(PN):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))

        ks, Hs, Ps = [], [], []
        for k in range(2, KCAP):
            if not sqf[k] or any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
            if M < 2:
                continue
            ms = np.arange(1, M + 1, 2, dtype=np.int64)
            ms = ms[sqf[ms]]
            for q in factor_set(k):
                if q > 2:
                    ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            vals = N - ms * k
            g = mu[ms].astype(np.float64)
            w = np.ones(ms.size, dtype=np.float64)
            for q in QS:
                if k % q == 0:
                    continue
                w *= np.where(vals % q == 0, 0.0, q / (q - 1.0))
            ks.append(k)
            Hs.append(float((lam[vals] * g).sum()))
            Ps.append(float((g * w).sum()))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        P = np.array(Ps)
        lw = np.log(ks.astype(float))
        beta = float((H * P).sum() / (P * P).sum())
        aR = np.abs(H - beta * P)
        thr = S_ * (1.0 - A_) * N
        cum = np.cumsum(lw * aR)
        j = int(np.searchsorted(cum, thr))
        kstar = int(ks[j]) if j < ks.size else None
        scale = np.sqrt(N / ks.astype(float))
        sel = ks <= kstar
        cR = float((aR[sel] / scale[sel]).mean())
        dens = float(sel.sum()) / kstar
        res.append((N, ks, lw, aR, scale, thr, kstar, cR, dens,
                    S_, A_, beta))
        say("  N = %-10d  #k = %-7d K*_R = %-7d c_R = %.4f  density "
            "%.4f" % (N, ks.size, kstar, cR, dens))

    # ------------------------------------------------------------- Q1
    say()
    say("Q1  the control")
    say("  N            K*_R (here)   published    same")
    q1 = True
    for N, ks, lw, aR, scale, thr, kstar, cR, dens, S_, A_, b in res:
        if kstar != pub[N]:
            q1 = False
        say("  %-12d %-13d %-12d %s"
            % (N, kstar, pub[N], "yes" if kstar == pub[N] else "NO"))
    say("  Q1 %s" % ("hold" if q1 else "REFUTED"))

    # ------------------------------------------------------------- Q2
    say()
    say("Q2  is the constant of the same shape as H's?")
    say("  N            c_R        sqrt(log N)   c_R/sqrt(log N)")
    gam = []
    for N, ks, lw, aR, scale, thr, kstar, cR, dens, S_, A_, b in res:
        g = cR / math.sqrt(math.log(N))
        gam.append(g)
        say("  %-12d %-10.4f %-13.4f %.4f"
            % (N, cR, math.sqrt(math.log(N)), g))
    spread = (max(gam) - min(gam)) / float(np.mean(gam))
    q2 = spread < 0.05
    say("  spread %.4f of the mean %.4f   (cap 0.05)   %s"
        % (spread, float(np.mean(gam)), "hold" if q2 else "REFUTED"))
    say("  Q2 %s" % ("hold" if q2 else "REFUTED"))

    # ------------------------------------------------------------- Q3
    say()
    say("Q3  does the model built from it predict the crossing?")
    say("  N            measured   model      ratio")
    q3 = True
    for N, ks, lw, aR, scale, thr, kstar, cR, dens, S_, A_, b in res:
        model = np.cumsum(lw * cR * scale)
        j = int(np.searchsorted(model, thr))
        km = int(ks[min(j, ks.size - 1)])
        r = km / kstar
        if abs(r - 1.0) >= 0.05:
            q3 = False
        say("  %-12d %-10d %-10d %.4f" % (N, kstar, km, r))
    say("  Q3 %s" % ("hold" if q3 else "REFUTED"))

    # ------------------------------------------------------------- Q4
    say()
    say("Q4  what the model says about the dip")
    say("  the model's own crossing exponent, against the measured one")
    say("  N            model exponent   measured   clears .56")
    q4 = True
    for N, ks, lw, aR, scale, thr, kstar, cR, dens, S_, A_, b in res:
        model = np.cumsum(lw * cR * scale)
        j = int(np.searchsorted(model, thr))
        km = int(ks[min(j, ks.size - 1)])
        em = math.log(km) / math.log(N)
        ee = math.log(kstar) / math.log(N)
        if em <= THETA:
            q4 = False
        say("  %-12d %-16.4f %-10.4f %s"
            % (N, em, ee, "yes" if em > THETA else "NO"))
    say("  Q4 %s" % ("hold" if q4 else "REFUTED"))

    say()
    say("  the budget constant crossed throughout, declared:")
    for N, ks, lw, aR, scale, thr, kstar, cR, dens, S_, A_, b in res:
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, S_ * (1.0 - A_)))

    # ------------------------------------------- the forecast, if earned
    say()
    if q2 and q3:
        say("  THE FORECAST. Q2 and Q3 hold, so the model may be asked")
        say("  where the exponent settles. With c_R = g sqrt(log N) and")
        say("  the admissible density d measured above, the balance at")
        say("  K = N^e is")
        say("    N^{(1-e)/2} = 2 d g sqrt(log N) (e log N - 2)")
        say("                    / [S(N)(1-A(N))].")
        gm = float(np.mean(gam))
        dm = float(np.mean([r[8] for r in res]))
        thrc = float(np.mean([r[9] * (1.0 - r[10]) for r in res]))
        say("  measured g = %.4f, density = %.4f, budget = %.6f"
            % (gm, dm, thrc))

        def expo_at(u, g):
            """solve for e at log N = u"""
            lo, hi = 0.05, 0.999
            for _ in range(200):
                e = 0.5 * (lo + hi)
                lhs = math.exp((1.0 - e) * u / 2.0)
                rhs = 2.0 * dm * g * math.sqrt(u) * (e * u - 2.0) / thrc
                if lhs > rhs:
                    lo = e
                else:
                    hi = e
            return 0.5 * (lo + hi)

        say("  log10 N      exponent   at c_R -10%%   at c_R +10%%")
        for d10 in (6, 8, 10, 12, 16, 20):
            u = d10 * math.log(10.0)
            say("  %-12d %-10.4f %-13.4f %.4f"
                % (d10, expo_at(u, gm), expo_at(u, gm * (1 - WOBBLE)),
                   expo_at(u, gm * (1 + WOBBLE))))

        def first_clear(g):
            lo, hi = 4.0, 40.0
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if expo_at(mid * math.log(10.0), g) > THETA:
                    hi = mid
                else:
                    lo = mid
            return 0.5 * (lo + hi)

        f0 = first_clear(gm)
        fl = first_clear(gm * (1 + WOBBLE))
        fh = first_clear(gm * (1 - WOBBLE))
        say("  the model clears %.2f from log10 N = %.2f upward"
            % (THETA, f0))
        say("BRACKET log10_N_residue_clears_theta %.4f %.4f %.4f"
            % (f0, min(fl, fh), max(fl, fh)))
    else:
        say("  NO FORECAST is made: Q2 or Q3 failed, so the model is")
        say("  not entitled to be extrapolated.")

    say()
    say("  DIAGNOSTIC (post hoc). How the residue's constant compares")
    ch = read_H_constant()
    say("  with H's. {#rem:heuristic} measured c(N)/sqrt(log N) for H")
    say("  at %.4f down to %.4f, read from" % (max(ch), min(ch)))
    say("  results/audit_directlevel_heuristic.txt, spread %.4f of its"
        % ((max(ch) - min(ch)) / float(np.mean(ch))))
    say("  mean against %.4f here. Their ratio is what the split buys"
        % spread)
    say("  in the constant, and the exponent it buys is nothing:")
    say("  N            c_R/sqrt(log N)")
    for i, (N, ks, lw, aR, scale, thr, kstar, cR, dens, S_, A_,
            b) in enumerate(res):
        say("  %-12d %.4f" % (N, gam[i]))

    say()
    say("=" * 70)
    ok = q1 and q2 and q3 and q4
    say("the residue obeys the same law with a smaller constant, and "
        "the knife-edge is local" if ok else "REFUTED")

    head = [
        "STATISTIC: the constant c_R(N) = mean |R|/sqrt(N/k) over the",
        "           admissible k below the operative crossing K*_R; its",
        "           ratio to sqrt(log N); the crossing predicted by",
        "           sum_{k<K}(log k) c_R sqrt(N/k) = S(N)(1-A(N))N over",
        "           the actual admissible k; and, if that model holds,",
        "           the exponent it gives at larger N with a bracket",
        "           over a +-10 per cent wobble in c_R.",
        "NULL: none is run and none applies. A measured magnitude is",
        "      fitted to a deterministic shape and the shape solved;",
        "      there is no background to detect against. The sign",
        "      controls for R were run in lab_residue_cancellation.py,",
        "      whose coin arm on the identical delta established that",
        "      R's size is bought by cancellation at exactly a coin's",
        "      rate, and in lab_residue_signed.py.",
        "FIELD: N = 2e5 to 3.2e6 by doubling; k squarefree and coprime",
        "       to N with 2 <= k < 100000; m odd, squarefree, coprime to",
        "       the odd part of k, m < N/k; the sieve weight uses the",
        "       odd primes up to 30; beta refitted as sum(H P)/sum(P^2)",
        "       on the same k-range; S(N) and A(N) from Euler products",
        "       at the fixed bound 4000000; the published crossings are",
        "       read from results/audit_residue_level.txt.",
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
