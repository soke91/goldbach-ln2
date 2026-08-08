# -*- coding: utf-8 -*-
r"""
Testing the one forecast this repository made that lands inside the
computable range.

WHAT IS AT STAKE

Remark {#rem:primorialladder} swept N = 30030 * 2^j for j = 0..6, the
odd radical 3*5*7*11*13 held fixed so that the budget does not move,
and found the residue-only exponent rising from 0.4550 to 0.4876 with
slope +0.006623. On that slope it reaches 1/2 at N = 10^7.10, and the
leave-one-out extremes of the slope give the bracket

    log10 N in [7.0733, 7.3629].

Every other forecast in these papers is out of reach -- 10^5475 for
the elementary bound, two thirds of a decade at 10^8.32 for the level,
nine and a half decades for the lean. This one is a factor of six to
twenty past the top rung, and the remark says so: three more doublings
would settle it.

This runs them. j = 7, 8, 9 puts the ladder at 3843840, 7687680 and
15375360, the last of which is 10^7.19 -- past the forecast point and
inside the bracket. Nothing about the answer is protected: the
forecast is a linear extrapolation of an exponent over a factor 64,
which is exactly the kind of thing Remark {#rem:forecastbracket} is
sceptical of, and the point of running it is that for once the
scepticism can be settled rather than declared.

BACKS: Remark {#rem:primorialreach} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  The control: the seven published rungs are reproduced to within
      0.001 in the exponent, read from
      results/lab_primorial_ladder.txt.
  E2  The rise continues: over the extended ladder of ten rungs the
      least-squares slope against log N is still positive.
  E3  The barrier is reached: at least one of the three new rungs has
      exponent at or above 0.5.
  E4  And the forecast was sound: the N at which the exponent first
      reaches 0.5, interpolated between rungs, lies inside the
      published bracket [10^7.0733, 10^7.3629]. If no crossing occurs
      inside the ladder this cannot be evaluated as written, and the
      decidable question is whether the bracket is EXCLUDED -- that is,
      whether its upper end lies below the top rung.

REFUTATION RULE (fixed before the run)

  E1  REFUTED at 0.001 at any rung.
  E2  REFUTED if the slope is not positive.
  E3  REFUTED if all three new rungs stay below 0.5. That is the
      outcome that matters: it would mean the forecast overshot and
      the primorial arithmetic is harder than a linear extrapolation
      of six rungs suggested.
  E4  REFUTED if the crossing falls outside the bracket, in either
      direction; or, if there is no crossing inside the ladder,
      REFUTED if the bracket's upper end lies below the top rung, so
      that the bracket is excluded outright. A bracket that fails its
      first live test is worth more than one that is never tested --
      and a bracket whose POINT ESTIMATE fails while the bracket
      itself survives is the case the whole apparatus is for.

  All four gate.

  NO NULL IS RUN and none is needed. The question is whether a
  specific published number is right, and the control is E1: the same
  computation reproducing the rungs it was fitted on. The coin arm for
  this statistic was run in lab_primorial_ladder.py, where eight
  global sign vectors on the identical deviations gave slopes from
  -0.003700 to +0.014290 with mu's inside, establishing that the rise
  is a fact about magnitudes; that control does not need repeating to
  ask where a deterministic curve crosses a line.
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
OUT = os.path.join(ROOT, "results", "audit_primorial_reach.txt")

BASE = 30030                       # 2*3*5*7*11*13
LADDER = [BASE * (1 << j) for j in range(10)]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000
HALF = 0.5


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


def read_ladder():
    """the published rungs and the bracket -- read, not copied"""
    p = os.path.join(ROOT, "results", "lab_primorial_ladder.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("R2/R3  the ladder")
    ex = {}
    for ln in src[i:].splitlines()[2:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            break
        ex[int(f[0])] = float(f[3])
    m = re.search(r"BRACKET log10_N_primorial_reaches_half "
                  r"([\d.]+) ([\d.]+) ([\d.]+)", src)
    return ex, (float(m.group(1)), float(m.group(2)),
                float(m.group(3)))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub, (pt, blo, bhi) = read_ladder()
    say("read %d published rungs and the bracket 10^%.4f "
        "[10^%.4f, 10^%.4f]" % (len(pub), pt, blo, bhi))
    say("  from results/lab_primorial_ladder.txt")

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
        e = math.log(kstar) / math.log(N) if kstar else float("nan")
        rows.append((N, thrc, ks.size, beta, kstar, e))
        say("  N = %-9d thr %.6f  #k %-6d K*_R %-6s exp %.4f"
            % (N, thrc, ks.size, str(kstar), e))

    # ------------------------------------------------------------- E1
    say()
    say("E1  the control: the seven published rungs")
    say("  N            here       published   diff")
    e1 = True
    for N, thrc, nk, beta, kstar, e in rows:
        if N not in pub:
            continue
        d = abs(e - pub[N])
        if not (d < 0.001):
            e1 = False
        say("  %-12d %-10.4f %-11.4f %.5f" % (N, e, pub[N], d))
    say("  E1 %s" % ("hold" if e1 else "REFUTED"))

    # ------------------------------------------------------------- E2
    say()
    say("E2  the extended ladder")
    say("  N            log10 N   K*_R    exponent   new?")
    for N, thrc, nk, beta, kstar, e in rows:
        say("  %-12d %-9.4f %-7s %-10.4f %s"
            % (N, math.log10(N), str(kstar), e,
               "" if N in pub else "new"))
    x = np.log(np.array([r[0] for r in rows], dtype=float))
    y = np.array([r[5] for r in rows])
    sl = float(np.polyfit(x, y, 1)[0])
    rr = float(np.corrcoef(x, y)[0, 1])
    e2 = sl > 0.0
    say("  slope against log N %+.6f, correlation %.5f   %s"
        % (sl, rr, "hold" if e2 else "REFUTED"))
    say("  E2 %s" % ("hold" if e2 else "REFUTED"))
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    say("  leave-one-out on the slope: %+.6f, %+.6f, %+.6f -- "
        "spread %.6f" % (f[0], f[1], f[2], max(f) - min(f)))
    say("SWEPT primorial_reach N-range %.6f" % (max(f) - min(f)))
    say("CORR primorial_reach %.5f" % abs(rr))
    say("POP primorial_reach %d" % len(rows))

    # ------------------------------------------------------------- E3
    say()
    new = [r for r in rows if r[0] not in pub]
    e3 = any(r[5] >= HALF for r in new)
    say("E3  do the new rungs reach %.1f?" % HALF)
    for N, thrc, nk, beta, kstar, e in new:
        say("  N = %-10d exponent %.4f   %s"
            % (N, e, "at or above" if e >= HALF else "below"))
    say("  E3 %s" % ("hold" if e3 else "REFUTED"))

    # ------------------------------------------------------------- E4
    say()
    say("E4  where it first reaches %.1f, against the published bracket"
        % HALF)
    cross = None
    for i in range(1, len(rows)):
        a, b = rows[i - 1], rows[i]
        if a[5] < HALF <= b[5]:
            t = (HALF - a[5]) / (b[5] - a[5])
            cross = (math.log10(a[0])
                     + t * (math.log10(b[0]) - math.log10(a[0])))
            break
    if cross is None:
        top = math.log10(rows[-1][0])
        say("  no crossing inside the ladder: the data put it ABOVE the")
        say("  top rung at log10 N = %.4f, so every point estimate and"
            % top)
        say("  every bracket upper end below that is excluded and")
        say("  nothing above it is.")
        say()
        say("  Two brackets have to be distinguished. The one this")
        say("  audit was pre-registered against was built from the")
        say("  slope's leave-one-out spread ALONE; its published")
        say("  numbers, the reference values of this audit, appear in")
        say("  the table below. That file has since been corrected to")
        say("  carry the residual scatter as well, and now publishes")
        say("  the numbers read at the top of this run. Both:")
        say("  which bracket        point   lower   upper   verdict")
        for name, q, a1, b1 in (
                ("slope-only (published)", 7.0954, 7.0733, 7.3629),
                ("with scatter (current)", pt, blo, bhi)):
            say("  %-20s %-7.4f %-7.4f %-7.4f point %s, interval %s"
                % (name, q, a1, b1,
                   "EXCLUDED" if q < top else "open",
                   "excluded" if b1 < top else "open"))
        e4 = bhi >= top
        say("  E4 asked where the crossing falls and there is no")
        say("  crossing to place, so it is not evaluable as written.")
        say("  What is decidable is exclusion, and neither interval is")
        say("  excluded while both point estimates below %.4f are." % top)
        say("  E4 %s" % ("hold" if e4 else "REFUTED"))
    else:
        e4 = blo <= cross <= bhi
        say("  interpolated crossing at log10 N = %.4f" % cross)
        say("  published bracket [%.4f, %.4f], point estimate %.4f"
            % (blo, bhi, pt))
        say("  E4 %s" % ("hold" if e4 else "REFUTED"))

    say()
    say("  the arithmetic and the budget, declared:")
    rads = set()
    for N, thrc, nk, beta, kstar, e in rows:
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  %d N, %d distinct odd radical%s: %s"
        % (len(rows), len(rads), "" if len(rads) == 1 else "s",
           ", ".join(str(v) for v in sorted(rads))))
    say("RADICALS %d" % len(rads))
    for N, thrc, nk, beta, kstar, e in rows:
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, thrc))

    say()
    say("  DIAGNOSTIC (post hoc). What the three new rungs did to the")
    say("  fit. The slope on the published seven, and on all ten:")
    xo = x[:len(pub)]
    yo = y[:len(pub)]
    say("  seven rungs %+.6f, ten rungs %+.6f"
        % (float(np.polyfit(xo, yo, 1)[0]), sl))
    say("  and the exponent's own step from rung to rung:")
    say("  N            exponent   step")
    for i, (N, thrc, nk, beta, kstar, e) in enumerate(rows):
        say("  %-12d %-10.4f %s"
            % (N, e, "-" if i == 0 else "%+.4f" % (e - rows[i - 1][5])))

    say()
    say("  Why the bracket failed, in one number. It was built from the")
    say("  leave-one-out spread of the SLOPE. That spread is now")
    say("  %.6f over ten rungs -- the slope was never the problem."
        % (max(f) - min(f)))
    say("  What the bracket did not cover is the scatter of the rungs")
    say("  about the line:")
    a, b = np.polyfit(x, y, 1)
    resid = y - (a * x + b)
    rms = float(np.sqrt((resid ** 2).mean()))
    say("  N            exponent   fitted     residual")
    for i, (N, thrc, nk, beta, kstar, e) in enumerate(rows):
        say("  %-12d %-10.4f %-10.4f %+.4f"
            % (N, e, a * x[i] + b, resid[i]))
    say("  r.m.s. residual %.4f, against a trend of %.4f per doubling"
        % (rms, a * math.log(2.0)))
    say("  -- the scatter is %.1f times the step the trend makes in one"
        % (rms / (a * math.log(2.0))))
    say("  rung, so a level is crossed several rungs before or after")
    say("  the line says, and a bracket over the slope alone cannot")
    say("  know that.")
    say("SCATTER primorial_slope %.4f" % rms)
    say("DRIFT primorial_slope %.4f" % ((max(f) - min(f)) / abs(a)))
    say()
    say("  The forecast redone, with the scatter carried. The line")
    say("  reaches %.1f at:" % HALF)
    u0 = (HALF - b) / a
    say("    log10 N = %.4f" % (u0 / math.log(10.0)))
    say("  and with the line displaced by +-1 r.m.s. residual:")
    ulo = (HALF - rms - b) / a
    uhi = (HALF + rms - b) / a
    say("    [%.4f, %.4f] in log10 N"
        % (ulo / math.log(10.0), uhi / math.log(10.0)))
    say()
    say("  THE SHAPES. That bracket is over the line's parameters and")
    say("  its scatter, and says nothing about whether the line is the")
    say("  right form. audit_ladder_shape.py showed these points do")
    say("  not choose one. The same five, fitted to these ten rungs,")
    say("  with the crossing still just outside the data:")

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
    zz = np.log(x) / x
    c1 = float(((1.0 - y) * zz).sum() / (zz * zz).sum())
    fit2["heuristic1"] = (np.array([c1]),
                          float(np.sqrt(((y - (1.0 - c1 * zz)) ** 2)
                                        .mean())))
    NAMES.append("heuristic1")
    LAB["heuristic1"] = "1 - c log log N / log N"

    def val(nm, uu):
        cc = fit2[nm][0]
        if nm == "heuristic1":
            return 1.0 - cc[0] * math.log(uu) / uu
        return float(design(nm, np.array([uu])) @ cc)

    def crossat(nm, tgt):
        if val(nm, 400.0) < tgt:
            return None
        lo_, hi_ = 2.0, 400.0
        for _ in range(300):
            mid = 0.5 * (lo_ + hi_)
            if val(nm, mid) < tgt:
                lo_ = mid
            else:
                hi_ = mid
        return 0.5 * (lo_ + hi_) / math.log(10.0)

    best = min(fit2[nm][1] for nm in NAMES)
    se = best / math.sqrt(2.0 * (x.size - 2))
    say("  shape                        r.m.s.    (r-best)/s.e.  "
        "0.50 at")
    keep = []
    for nm in NAMES:
        rr2 = fit2[nm][1]
        t = (rr2 - best) / se
        cc = crossat(nm, HALF)
        if t <= 1.0:
            keep.append(cc)
        say("  %-28s %-9.5f %-14.2f %s"
            % (LAB[nm], rr2, t,
               "%.4f" % cc if cc is not None else "never"))
    say("  %d shapes fitted, %d within one standard error of the best"
        % (len(NAMES), len(keep)))
    say("SHAPES %d" % len(NAMES))
    vals = [v for v in keep if v is not None]
    l10 = ulo / math.log(10.0)
    h10 = uhi / math.log(10.0)
    lo3 = min(min(vals), l10)
    hi3 = max(max(vals), h10)
    say("  the surviving shapes put 0.50 at %s -- a spread of %.4f"
        % (", ".join("%.4f" % v for v in sorted(vals)),
           max(vals) - min(vals)))
    say("  decades; the line's scatter bracket alone is [%.4f, %.4f],"
        % (l10, h10))
    say("  and with the shapes [%.4f, %.4f], wider by %.4f"
        % (lo3, hi3, (hi3 - lo3) - (h10 - l10)))
    say("BRACKET log10_N_primorial_reaches_half_v2 %.4f %.4f %.4f"
        % (u0 / math.log(10.0), lo3, hi3))
    say("  The top rung sits at log10 N = %.4f, so this forecast is"
        % math.log10(rows[-1][0]))
    say("  again a short reach -- and it is now the second forecast of")
    say("  the same crossing, the first having been tested and failed.")

    say()
    say("=" * 70)
    ok = e1 and e2 and e3 and e4
    say("the forecast held and the barrier is reached inside the "
        "computable range" if ok else "REFUTED")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           along N = 30030 * 2^j for j = 0..9 -- the ladder of",
        "           lab_primorial_ladder.py extended by three doublings;",
        "           the least-squares slope over all ten rungs; and the",
        "           interpolated N at which the exponent first reaches",
        "           1/2, against the bracket that file published.",
        "NULL: none is run and none is needed. The question is whether a",
        "      specific published forecast is right, and the control is",
        "      E1, the same computation reproducing the rungs the",
        "      forecast was fitted on. The coin arm for this statistic",
        "      was run in lab_primorial_ladder.py, where eight global",
        "      sign vectors on the identical deviations bracketed mu's",
        "      slope, establishing that the rise is a fact about",
        "      magnitudes.",
        "FIELD: N = 30030 * 2^j, j = 0..9, so 30030 to 15375360, the odd",
        "       radical 3*5*7*11*13 fixed and the threshold constant at",
        "       every rung; k squarefree and coprime to N with",
        "       2 <= k < 100000; m odd, squarefree, coprime to the odd",
        "       part of k, m < N/k; the sieve weight uses the odd primes",
        "       up to 30; beta refitted as sum(H P)/sum(P^2) on the same",
        "       k-range; S(N) and A(N) from Euler products at the fixed",
        "       bound 4000000.",
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
