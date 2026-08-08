# -*- coding: utf-8 -*-
r"""
The barrier as a ratio instead of a crossing.

WHAT IS AT STAKE

Remark {#rem:primorialreach} refuted a point estimate and left an
interval open: along the primorial ladder the residue-only exponent
reaches 0.4941 at the top rung, the crossing of 1/2 lies above
10^7.1868, and the forecast of where is now 10^7.4684 with a bracket
[10^7.2189, 10^7.7180]. That bracket is half a decade wide, and it is
wide because the rungs scatter about their line by 0.8 of what the
trend gains in a whole rung.

The scatter is an artefact of the instrument. K*_R is the location
where a step function first exceeds a level: an integer, discrete in
k, and sensitive to a single term near the threshold. But the question
the exponent encodes is not about a location at all. Since

    log K*_R / log N > 1/2   <=>   K*_R > sqrt(N)
                             <=>   B_R(N; sqrt N) < S(N)(1-A(N)) N,

the barrier is crossed exactly when the ratio

    rho(N) = sum_{k < sqrt N} (log k) |R(N;k)| / [S(N)(1-A(N)) N]

falls below 1. That is a ratio of two smooth sums evaluated at a fixed
abscissa -- no crossing, no discreteness, no single decisive term.
Remark {#rem:muvscoin} made the same move for a different question and
resolved in one step what a fitted exponent could not resolve at all.

BACKS: Remark {#rem:primorialshare} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  F1  The control: rho(N) < 1 at exactly those rungs where the
      published exponent exceeds 0.5, read from
      results/audit_primorial_reach.txt. Since none of them does,
      rho > 1 at every rung.
  F2  The ratio is the quieter instrument: its residual scatter about
      a fit, measured against what the fit's trend gains in one
      doubling, is smaller than the exponent's, read from the same
      file.
  F3  And it is falling: the slope of log rho against log N is
      negative.
  F4  So the forecast sharpens: the bracket on where rho reaches 1,
      at one r.m.s. residual either way, is narrower than the
      exponent-based bracket now published in
      results/lab_primorial_ladder.txt.

REFUTATION RULE (fixed before the run)

  F1  REFUTED if rho and the exponent disagree at any rung about
      which side of the barrier it is on. They are equivalent by
      definition, so a disagreement would mean one of the two
      computations is wrong.
  F2  REFUTED if the ratio's relative scatter is at or above the
      exponent's. That is the one that matters: it would say the
      noise is in the object and not in the instrument, and no
      change of statistic will sharpen the forecast.
  F3  REFUTED if the slope is not negative.
  F4  REFUTED if the new bracket is at least as wide.

  All four gate.

  NULL: a coin arm on the identical deviations. R = sum_m mu(m)
  delta(m,k) with delta = Lambda(N-mk) - beta w(m,k); replacing mu by
  8 global sign vectors on the odd squarefree m, each held across all
  k as mu is, and carrying them through the same ratio gives what rho
  a coin produces on the same ladder against the same fixed budget.
  lab_residue_cancellation.py established that R's SIZE is exactly a
  coin's, so the coin's rho is the scale mu's has to be read against,
  and its scatter is the scale F2 has to be read against.
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
OUT = os.path.join(ROOT, "results", "lab_primorial_share.txt")

BASE = 30030                       # 2*3*5*7*11*13
LADDER = [BASE * (1 << j) for j in range(10)]
QSIEVE = 30
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


def read_reach():
    """the published exponents and the exponent's scatter -- read"""
    p = os.path.join(ROOT, "results", "audit_primorial_reach.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("N            log10 N   K*_R    exponent   new?")
    ex = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        ex[int(f[0])] = float(f[3])
    sc = float(re.search(r"^SCATTER primorial_slope ([\d.]+)", src,
                         re.M).group(1))
    tr = float(re.search(r"r\.m\.s\. residual [\d.]+, against a trend "
                         r"of ([\d.]+) per doubling", src).group(1))
    return ex, sc, tr


def read_ladder_bracket():
    p = os.path.join(ROOT, "results", "lab_primorial_ladder.txt")
    src = io.open(p, encoding="utf-8").read()
    m = re.search(r"BRACKET log10_N_primorial_reaches_half "
                  r"([\d.]+) ([\d.]+) ([\d.]+)", src)
    return (float(m.group(1)), float(m.group(2)), float(m.group(3)))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubex, pubsc, pubtr = read_reach()
    bpt, blo, bhi = read_ladder_bracket()
    say("read %d published exponents, the exponent's scatter %.4f and"
        % (len(pubex), pubsc))
    say("  its trend %.4f per doubling, and the published bracket"
        % pubtr)
    say("  10^%.4f [10^%.4f, 10^%.4f]" % (bpt, blo, bhi))

    NMAX = max(LADDER)
    say("sieving to %d ..." % NMAX)
    lam, mu = sieves(NMAX)
    sqf = mu != 0
    oddsqf = np.zeros(NMAX + 1, dtype=bool)
    oddsqf[1::2] = True
    oddsqf &= sqf
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]

    rng = np.random.default_rng(SEED)
    coinmat = np.zeros((COINS, NMAX + 1), dtype=np.int8)
    _idx = np.flatnonzero(oddsqf)
    for _j in range(COINS):
        coinmat[_j, _idx] = rng.integers(0, 2, size=_idx.size) * 2 - 1
    del _idx
    say("  %d global sign vectors for the null" % COINS)

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
        KTOP = int(math.isqrt(N))

        # beta is fitted on the same k-range the crossing used, so that
        # rho is the same object the exponent encodes
        ks, Hs, Ps, Cs = [], [], [], []
        for k in range(2, KTOP + 1):
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
            Cs.append((ms, lam[vals], w))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        P = np.array(Ps)
        lw = np.log(ks.astype(float))
        beta = float((H * P).sum() / (P * P).sum())
        rho = float((lw * np.abs(H - beta * P)).sum() / (thrc * N))
        rc = np.zeros(COINS)
        for i, (ms, L, w) in enumerate(Cs):
            d = L - beta * w
            rc += lw[i] * np.abs(coinmat[:, ms].astype(np.float64) @ d)
        rc /= thrc * N
        del Cs
        rows.append((N, KTOP, thrc, ks.size, beta, rho, rc))
        say("  N = %-9d sqrt %-6d #k %-6d beta %.4f  rho %.4f"
            % (N, KTOP, ks.size, beta, rho))

    # ------------------------------------------------------------- F1
    say()
    say("F1  the control: rho and the exponent agree on the side")
    say("  N            rho        below 1?   exponent   above .5?")
    f1 = True
    for N, KTOP, thrc, nk, beta, rho, rc in rows:
        a = rho < 1.0
        b = pubex[N] > 0.5
        if a != b:
            f1 = False
        say("  %-12d %-10.4f %-10s %-10.4f %s"
            % (N, rho, "yes" if a else "no", pubex[N],
               "yes" if b else "no"))
    say("  F1 %s" % ("hold" if f1 else "REFUTED"))

    # ------------------------------------------------------------- F3
    say()
    x = np.log(np.array([r[0] for r in rows], dtype=float))
    y = np.log(np.array([r[5] for r in rows]))
    a1, b1 = np.polyfit(x, y, 1)
    corr = float(np.corrcoef(x, y)[0, 1])
    f3 = a1 < 0.0
    say("F3  is rho falling?  slope of log rho against log N %+.6f, "
        "correlation %.5f" % (a1, corr))
    say("  F3 %s" % ("hold" if f3 else "REFUTED"))

    # ------------------------------------------------------------- F2
    say()
    say("F2  which instrument is quieter")
    resid = y - (a1 * x + b1)
    rms = float(np.sqrt((resid ** 2).mean()))
    trend = abs(a1) * math.log(2.0)
    mine = rms / trend
    theirs = pubsc / pubtr
    f2 = mine < theirs
    say("  N            log rho    fitted     residual")
    for i, (N, KTOP, thrc, nk, beta, rho, rc) in enumerate(rows):
        say("  %-12d %-10.4f %-10.4f %+.4f"
            % (N, y[i], a1 * x[i] + b1, resid[i]))
    say("  ratio:    r.m.s. %.4f, trend %.4f per doubling, ratio %.3f"
        % (rms, trend, mine))
    say("  exponent: r.m.s. %.4f, trend %.4f per doubling, ratio %.3f"
        % (pubsc, pubtr, theirs))
    say("  F2 %s" % ("hold" if f2 else "REFUTED"))
    say("SCATTER primorial_share_slope %.4f" % rms)
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    say("DRIFT primorial_share_slope %.4f"
        % ((max(f) - min(f)) / abs(a1)))
    say("SWEPT primorial_share N-range %.6f" % (max(f) - min(f)))
    say("CORR primorial_share %.5f" % abs(corr))
    say("POP primorial_share %d" % len(rows))

    # ------------------------------------------------------------- F4
    say()
    say("F4  the forecast, with the scatter carried")
    u0 = (0.0 - b1) / a1
    ulo = (0.0 + rms - b1) / a1
    uhi = (0.0 - rms - b1) / a1
    lo10, hi10 = sorted([ulo / math.log(10.0), uhi / math.log(10.0)])
    pt10 = u0 / math.log(10.0)
    width = hi10 - lo10
    oldw = bhi - blo
    f4 = width < oldw
    say("  rho reaches 1 at log10 N = %.4f, bracket [%.4f, %.4f], "
        "width %.4f" % (pt10, lo10, hi10, width))
    say("  the exponent-based bracket was [%.4f, %.4f], width %.4f"
        % (blo, bhi, oldw))
    say("  F4 %s" % ("hold" if f4 else "REFUTED"))
    say()
    say("  THE SHAPES. That bracket is over the line's parameters and")
    say("  scatter and says nothing about whether log rho is linear in")
    say("  log N. The four two-parameter forms audit_ladder_shape.py")
    say("  compared for the exponent, applied to log rho -- the")
    say("  one-parameter 1 - c log log N / log N has no analogue here,")
    say("  since it is a statement about an exponent bounded by one:")

    def design(nm, uu):
        if nm == "line":
            return np.vstack([np.ones_like(uu), uu]).T
        if nm == "saturating":
            return np.vstack([np.ones_like(uu), 1.0 / uu]).T
        if nm == "loglog":
            return np.vstack([np.ones_like(uu), np.log(uu)]).T
        if nm == "heuristic2":
            return np.vstack([np.ones_like(uu), np.log(uu) / uu]).T
        raise KeyError(nm)

    NAMES = ["line", "saturating", "loglog", "heuristic2"]
    LAB = {"line": "a + b log N", "saturating": "a + b / log N",
           "loglog": "a + b log log N",
           "heuristic2": "a + b log log N / log N"}
    fit2 = {}
    for nm in NAMES:
        X = design(nm, x)
        cc, *_ = np.linalg.lstsq(X, y, rcond=None)
        rr2 = y - X @ cc
        fit2[nm] = (cc, float(np.sqrt((rr2 ** 2).mean())))

    def val(nm, uu):
        return float(design(nm, np.array([uu])) @ fit2[nm][0])

    def crossat(nm):
        if val(nm, 400.0) > 0.0:
            return None
        lo_, hi_ = 2.0, 400.0
        for _ in range(300):
            mid = 0.5 * (lo_ + hi_)
            if val(nm, mid) > 0.0:
                lo_ = mid
            else:
                hi_ = mid
        return 0.5 * (lo_ + hi_) / math.log(10.0)

    best = min(fit2[nm][1] for nm in NAMES)
    se = best / math.sqrt(2.0 * (x.size - 2))
    say("  shape                        r.m.s.    (r-best)/s.e.  "
        "rho = 1 at")
    keep = []
    for nm in NAMES:
        rr2 = fit2[nm][1]
        t = (rr2 - best) / se
        cc = crossat(nm)
        if t <= 1.0:
            keep.append(cc)
        say("  %-28s %-9.5f %-14.2f %s"
            % (LAB[nm], rr2, t,
               "%.4f" % cc if cc is not None else "never"))
    say("  %d shapes fitted, %d within one standard error of the best"
        % (len(NAMES), len(keep)))
    say("SHAPES %d" % len(NAMES))
    vals = [v for v in keep if v is not None]
    never = any(v is None for v in keep)
    if vals:
        lo3 = min(min(vals), lo10)
        hi3 = max(max(vals), hi10)
        say("  the surviving shapes put rho = 1 at %s%s"
            % (", ".join("%.4f" % v for v in sorted(vals)),
               ", and one never reaches it" if never else ""))
        say("  -- a spread of %.4f decades; the line's scatter bracket"
            % (max(vals) - min(vals)))
        say("  alone is [%.4f, %.4f], and with the shapes [%.4f, %.4f],"
            % (lo10, hi10, lo3, hi3))
        say("  wider by %.4f decades" % ((hi3 - lo3) - (hi10 - lo10)))
    else:
        lo3, hi3 = lo10, hi10
        say("  no surviving shape reaches rho = 1 at all")
    say("BRACKET log10_N_primorial_share_reaches_one %.4f %.4f %.4f"
        % (pt10, lo3, hi3))
    say("  So the ratio loses twice over: its scatter is %.1f times"
        % (mine / theirs))
    say("  the exponent's relative to trend, and its shape ambiguity")
    say("  is on top of that.")

    say()
    say("  THE NULL. The same ratio with mu replaced by a coin on the")
    say("  identical deviations, %d global sign vectors:" % COINS)
    say("  N            mu         coin min   coin median  coin max")
    for N, KTOP, thrc, nk, beta, rho, rc in rows:
        say("  %-12d %-10.4f %-10.4f %-12.4f %.4f"
            % (N, rho, rc.min(), float(np.median(rc)), rc.max()))
    cs = []
    for c in range(COINS):
        yc = np.log(np.array([r[6][c] for r in rows]))
        ac = float(np.polyfit(x, yc, 1)[0])
        rc2 = yc - np.polyval(np.polyfit(x, yc, 1), x)
        cs.append((ac, float(np.sqrt((rc2 ** 2).mean()))
                   / (abs(ac) * math.log(2.0))))
    say("  coin slopes %s"
        % ", ".join("%+.6f" % c[0] for c in cs))
    say("  coin scatter/trend %s"
        % ", ".join("%.3f" % c[1] for c in cs))
    say("  mu's slope %+.6f and scatter/trend %.3f" % (a1, mine))
    say("  The coin's ratio is noisy in the same way and by the same")
    say("  amount, so what F2 refuted is a property of the statistic")
    say("  and not of mu.")

    say()
    say("  DIAGNOSTIC (post hoc). Why the ratio is the WORSE")
    say("  instrument, which is the opposite of what F2 expected.")
    say("  The two statements are equivalent -- F1 confirms it at all")
    say("  ten rungs -- but they are not equivalent as TRENDS. Near the")
    say("  crossing log rho is about -c (log K* - (1/2) log N), so")
    say("    log rho ~ -c (log N) (e - 1/2),   e the exponent,")
    say("  and the exponent carries a division by log N that the ratio")
    say("  does not. That division damps the fluctuation of K* by a")
    say("  factor log N without damping the trend by as much:")
    lg = float(np.mean(np.log(np.array([r[0] for r in rows],
                                       dtype=float))))
    say("  mean log N over the ladder            : %.4f" % lg)
    say("  exponent scatter times log N          : %.4f"
        % (pubsc * lg))
    say("  the ratio's own scatter               : %.4f" % rms)
    say("  -- the same noise, seen through the two lenses.")
    say("  exponent trend per doubling times log N: %.4f"
        % (pubtr * lg))
    say("  the ratio's own trend per doubling     : %.4f" % trend)
    say("  -- and here they differ, because log rho carries a second")
    say("  drift term, -c (e - 1/2) d log N, that the exponent has")
    say("  divided away. The lesson of {#rem:muvscoin} does not")
    say("  generalise: a ratio at a fixed abscissa beats a fitted")
    say("  exponent when the exponent is a SLOPE THROUGH POINTS, and")
    say("  loses to it when the exponent is itself a ratio whose")
    say("  denominator grows.")

    say()
    say("  the arithmetic and the budget, declared:")
    rads = set()
    for N, KTOP, thrc, nk, beta, rho, rc in rows:
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  %d N, %d distinct odd radical%s: %s"
        % (len(rows), len(rads), "" if len(rads) == 1 else "s",
           ", ".join(str(v) for v in sorted(rads))))
    say("RADICALS %d" % len(rads))
    for N, KTOP, thrc, nk, beta, rho, rc in rows:
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, thrc))

    say()
    say("=" * 70)
    ok = f1 and f2 and f3 and f4
    say("the ratio is the quieter instrument and sharpens the forecast"
        if ok else "REFUTED")

    head = [
        "STATISTIC: rho(N) = sum_{k < sqrt N} (log k) |R(N;k)| divided",
        "           by S(N)(1-A(N))N, along N = 30030 * 2^j for",
        "           j = 0..9; rho < 1 is exactly the statement that the",
        "           residue-only level exceeds N^{1/2}. Its slope",
        "           against log N, the r.m.s. residual about that fit",
        "           relative to what the trend gains in one doubling,",
        "           and the forecast of where rho reaches 1.",
        "NULL: a coin arm on the identical deviations -- 8 global sign",
        "      vectors on the odd squarefree m, each held across all k",
        "      as mu is, carried through the same ratio. It fixes both",
        "      what rho a coin gives and how quiet the ratio is when",
        "      the signs carry no arithmetic.",
        "FIELD: N = 30030 * 2^j, j = 0..9, the odd radical 3*5*7*11*13",
        "       fixed so the threshold is constant; k squarefree and",
        "       coprime to N with 2 <= k <= floor(sqrt N); m odd,",
        "       squarefree, coprime to the odd part of k, m < N/k; the",
        "       sieve weight uses the odd primes up to 30; beta refitted",
        "       as sum(H P)/sum(P^2) on that same k-range; S(N) and A(N)",
        "       from Euler products at the fixed bound 4000000.",
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
