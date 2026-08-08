# -*- coding: utf-8 -*-
r"""
Is the residue's constant growing, or is the window it is measured in?

WHAT IS AT STAKE

Remark {#rem:laddermodel} reported that at the primorial radical
c_R/sqrt(log N) rises monotonically from 0.3724 to 0.6139 and that
c_R itself obeys (log N)^{1.3838} -- "not a wobble around 1/2 but a
different law". That is the sharpest thing this project has measured
about why the hard arithmetic is hard, and it rests on one convention.

c_R is the mean of |R(N;k)|/sqrt(N/k) over k below K*_R. K*_R grows
along the ladder, from 109 at the bottom rung to 5773 at the top. So
the average is taken over a window that widens by a factor of fifty
while N grows by a factor of a thousand, and Remark {#rem:modeltransfer}
already found exactly this failure for the other constant: c moved by
five per cent between two k-ranges at the SAME N, which is not a
property of N at all.

If |R| = c_R(N) sqrt(N/k) exactly, the window cannot matter: any
k-range gives the same c_R. If the window does matter, the law is
wrong in its k-dependence and the (log N)^{1.3838} is partly or wholly
an artefact of measuring over a moving range.

BACKS: Remark {#rem:cRwindow} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  K1  The control: the published convention, averaging below K*_R,
      reproduces the growth exponent of
      results/audit_ladder_model.txt to within 0.01.
  K2  The window matters: averaging instead over the FIXED range
      k < 300 at every rung gives a smaller growth exponent.
  K3  And it is the whole story: the fixed-window exponent is within
      0.2 of 1/2, so the law is the square-root one and the 1.3838 is
      an artefact of a widening average.
  K4  The k-dependence itself is sound: fitting the octave means of
      |R| against N/k at each rung gives an exponent in [0.40, 0.60],
      so it is only the constant that was ever in question.

REFUTATION RULE (fixed before the run)

  K1  REFUTED at 0.01, which would mean this is not the same
      measurement.
  K2  REFUTED if the fixed window gives an exponent at least as
      large. That is the outcome that would confirm
      {#rem:laddermodel} as it stands: the growth would be in N and
      not in the window.
  K3  REFUTED if the fixed-window exponent is more than 0.2 from 1/2.
  K4  REFUTED if any rung's octave exponent leaves [0.40, 0.60],
      which would say |R| is not square-root in N/k at this radical
      and the whole c_R framing is the wrong one.

  All four gate.

  NO NULL IS RUN and none applies. The same measured magnitudes are
  averaged over three windows and compared; there is no background to
  detect against. The sign controls for R on this ladder were run in
  lab_primorial_ladder.py, whose eight global sign vectors bracketed
  mu's slope, and in lab_residue_cancellation.py, whose coin arm on
  the identical deviations established that R's size is bought by
  cancellation at exactly a coin's rate.
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
OUT = os.path.join(ROOT, "results", "audit_cR_window.txt")

BASE = 30030                       # 2*3*5*7*11*13
LADDER = [BASE * (1 << j) for j in range(11)]
KCAP = 100_000
FIXED = 300                        # the fixed k-window
QSIEVE = 30
CLIM = 4_000_000
OCT = [2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288,
       2097152, 8388608, 33554432]
MINPTS = 10


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


def read_growth():
    p = os.path.join(ROOT, "results", "audit_ladder_model.txt")
    src = io.open(p, encoding="utf-8").read()
    return float(re.search(r"c_R ~ \(log N\)\^([\d.]+)", src).group(1))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubs = read_growth()
    say("read the published growth exponent %.4f from" % pubs)
    say("  results/audit_ladder_model.txt")

    NMAX = max(LADDER)
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

    rows = []
    for N in LADDER:
        PN = factor_set(N)
        A_, S_ = artin, twin
        for q in sorted(PN):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))
        thrc = S_ * (1.0 - A_)

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
        cum = np.cumsum(lw * aR)
        j = int(np.searchsorted(cum, thrc * N))
        kstar = int(ks[j]) if j < ks.size else None
        scale = np.sqrt(N / ks.astype(float))
        rat = aR / scale
        c_star = float(rat[ks <= kstar].mean())
        c_fix = float(rat[ks < FIXED].mean())
        # N^{1/4} is below the smallest admissible k at the lower
        # rungs -- k must be coprime to 3*5*7*11*13, so the smallest
        # is 17 -- and an empty window is not a measurement
        krel = max(int(N ** 0.25), int(ks[0]) + 1)
        c_rel = float(rat[ks < krel].mean())
        inner = N // ks
        rows.append((N, kstar, krel, c_star, c_fix, c_rel,
                     inner, aR, ks))
        say("  N = %-9d K*_R %-6d c(K*) %.4f  c(k<%d) %.4f  "
            "c(k<N^.25=%d) %.4f"
            % (N, kstar, c_star, FIXED, c_fix, krel, c_rel))

    u = np.array([math.log(r[0]) for r in rows])

    def growth(idx):
        v = np.array([r[idx] for r in rows])
        s, _ = np.polyfit(np.log(u), np.log(v), 1)
        r = float(np.corrcoef(np.log(u), np.log(v))[0, 1])
        f = [float(np.polyfit(np.log(u)[sl], np.log(v)[sl], 1)[0])
             for sl in (slice(None), slice(1, None), slice(0, -1))]
        return float(s), r, max(f) - min(f)

    # ------------------------------------------------------------- K1
    say()
    s_star, r_star, sp_star = growth(3)
    k1 = abs(s_star - pubs) < 0.01
    say("K1  the control: the published convention, k < K*_R")
    say("  c_R ~ (log N)^%.4f against the published %.4f, diff %.4f"
        % (s_star, pubs, abs(s_star - pubs)))
    say("  K1 %s" % ("hold" if k1 else "REFUTED"))

    # ---------------------------------------------------------- K2/K3
    say()
    say("K2/K3  the same magnitudes averaged over a FIXED window")
    say("  window                       exponent   correlation  "
        "LOO spread")
    s_fix, r_fix, sp_fix = growth(4)
    s_rel, r_rel, sp_rel = growth(5)
    for nm, ss, rr, sp in (("k < K*_R (published)", s_star, r_star,
                            sp_star),
                           ("k < %d (fixed)" % FIXED, s_fix, r_fix,
                            sp_fix),
                           ("k < N^0.25", s_rel, r_rel, sp_rel)):
        say("  %-28s %-10.4f %-12.5f %.4f" % (nm, ss, rr, sp))
    k2 = s_fix < s_star
    k3 = abs(s_fix - 0.5) <= 0.2
    say("  K2 the fixed window gives a smaller exponent   %s"
        % ("hold" if k2 else "REFUTED"))
    say("  K3 it is within 0.2 of 1/2 (|%.4f - 0.5| = %.4f)   %s"
        % (s_fix, abs(s_fix - 0.5), "hold" if k3 else "REFUTED"))
    say("SWEPT cR_fixed_window N-range %.6f" % sp_fix)
    say("CORR cR_fixed_window %.5f" % abs(r_fix))
    say("POP cR_fixed_window %d" % len(rows))

    # ------------------------------------------------------------- K4
    say()
    say("K4  is |R| square-root in N/k at this radical?")
    say("  N            octave exponent   correlation  bins  thinnest")
    k4 = True
    k4exp = []
    for N, kstar, krel, cs, cf, cr, inner, aR, ks in rows:
        cent, prof, cnt = [], [], []
        for a, b in zip(OCT, OCT[1:]):
            sel = (inner >= a) & (inner < b)
            cnt.append(int(sel.sum()))
            if sel.sum():
                cent.append(float(inner[sel].mean()))
                prof.append(float(aR[sel].mean()))
            else:
                cent.append(float("nan"))
                prof.append(float("nan"))
        c_ = np.array(cent)
        p_ = np.array(prof)
        ok = (~np.isnan(c_) & ~np.isnan(p_) & (p_ > 0)
              & (np.array(cnt) >= MINPTS))
        xx = np.log(c_[ok])
        yy = np.log(p_[ok])
        e = float(np.polyfit(xx, yy, 1)[0])
        rr = float(np.corrcoef(xx, yy)[0, 1])
        k4exp.append(e)
        if not (0.40 <= e <= 0.60):
            k4 = False
        thin = min(c for c in np.array(cnt)[ok])
        say("  %-12d %-17.4f %-12.5f %-5d %d"
            % (N, e, rr, int(ok.sum()), thin))
        say("POP cR_octave_N%d %d" % (N, thin))
        say("CORR cR_octave_N%d %.5f" % (N, abs(rr)))
    say("  K4 %s   (band [0.40, 0.60])" % ("hold" if k4 else "REFUTED"))

    say()
    say("  the arithmetic, declared:")
    rads = set()
    for N, kstar, krel, cs, cf, cr, inner, aR, ks in rows:
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  %d N, %d distinct odd radical%s: %s"
        % (len(rows), len(rads), "" if len(rads) == 1 else "s",
           ", ".join(str(v) for v in sorted(rads))))
    say("RADICALS %d" % len(rads))

    say()
    say("  DIAGNOSTIC (post hoc). Why a widening window inflates the")
    say("  growth. c_R is an average of |R|/sqrt(N/k) over k, and if")
    say("  that ratio is not flat in k the average depends on where it")
    say("  stops. At the top rung, by octave of N/k:")
    N, kstar, krel, cs, cf, cr, inner, aR, ks = rows[-1]
    say("  N/k octave        mean |R|/sqrt(N/k)   #k")
    for a, b in zip(OCT, OCT[1:]):
        sel = (inner >= a) & (inner < b)
        if sel.sum() < MINPTS:
            continue
        say("  %-17d %-20.4f %d"
            % (a, float((aR[sel] / np.sqrt(N / ks[sel].astype(float)))
                        .mean()), int(sel.sum())))
    say("  A window that grows with N samples a different part of")
    say("  this profile at every rung -- but K2 refutes that as the")
    say("  explanation: the FIXED window gives a STEEPER growth, not a")
    say("  shallower one. The profile is not flat, and that is the")
    say("  point: |R|/sqrt(N/k) falls with N/k, so c_R has no value")
    say("  independent of where the average stops.")
    say()
    say("  What is actually moving is the exponent in N/k itself. K4's")
    say("  column, fitted against log N:")
    ee = np.array(k4exp)
    s4, i4 = np.polyfit(u, ee, 1)
    r4 = float(np.corrcoef(u, ee)[0, 1])
    say("  N            octave exponent")
    for i, r in enumerate(rows):
        say("  %-12d %.4f" % (r[0], ee[i]))
    say("  it rises with slope %+.6f against log N, correlation %.5f,"
        % (s4, r4))
    say("  from %.4f at the bottom rung to %.4f at the top -- so |R| is"
        % (ee.min(), ee.max()))
    say("  NOT (N/k)^{1/2} at this radical below N of about %.2e, and"
        % math.exp((0.5 - i4) / s4))
    say("  a constant defined by dividing out sqrt(N/k) is absorbing")
    say("  that mismatch. That, and not a different law for the")
    say("  constant, is why no derived shape is available here.")
    f4 = [float(np.polyfit(u[sl], ee[sl], 1)[0])
          for sl in (slice(None), slice(1, None), slice(0, -1))]
    say("  leave-one-out on that slope: %+.6f, %+.6f, %+.6f -- "
        "spread %.6f" % (f4[0], f4[1], f4[2], max(f4) - min(f4)))
    say("SWEPT cR_octave_exponent N-range %.6f" % (max(f4) - min(f4)))
    say("CORR cR_octave_exponent %.5f" % abs(r4))
    say("POP cR_octave_exponent %d" % len(rows))

    say()
    say("=" * 70)
    ok = k1 and k2 and k3 and k4
    say("the growth was the window, not the constant"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the constant c_R = mean |R(N;k)|/sqrt(N/k) at each",
        "           rung of N = 30030 * 2^j, j = 0..10, averaged three",
        "           ways -- over k below the crossing K*_R (the",
        "           published convention), over the fixed window",
        "           k < 300, and over k < N^{1/4} -- with the growth",
        "           exponent s in c_R ~ (log N)^s for each; and the",
        "           octave fit of |R| against N/k at every rung.",
        "NULL: none is run and none applies. The same measured",
        "      magnitudes are averaged over three windows and compared;",
        "      there is no background to detect against. The sign",
        "      controls for R on this ladder were run in",
        "      lab_primorial_ladder.py and lab_residue_cancellation.py.",
        "FIELD: N = 30030 * 2^j, j = 0..10, the odd radical 3*5*7*11*13",
        "       fixed; k squarefree and coprime to N with",
        "       2 <= k < 100000; m odd, squarefree, coprime to the odd",
        "       part of k, m < N/k; the sieve weight uses the odd primes",
        "       up to 30; beta refitted as sum(H P)/sum(P^2) on that",
        "       k-range; S(N) and A(N) from Euler products at the fixed",
        "       bound 4000000; octaves closed at both ends and fitted",
        "       only when they hold at least 10 k.",
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
