# -*- coding: utf-8 -*-
r"""
Is the failure at primorial-like N structural, or is it finite-N?

WHAT IS AT STAKE

Remark {#rem:residuearithmetic} found the conditional reduction of
{#rem:provablehalf} failing exactly where Proposition [prop:onesided]
said it would: at N = 1531530 and 1621620 the residue-only level sits
at 0.4808 and 0.4747, below the square-root barrier, while the
2^a 5^b family clears it by 0.06 to 0.08. Seven N of one size cannot
say whether that is a fact about the arithmetic or about the size.

The arithmetic is not obviously fatal. At N primorial to y the
one-sided margin collapses -- 1 - A(N) = sum_{p>y} 1/(p(p-1)), about
1/(y log y) -- while S(N) grows like log y, so the budget is of order
N/log N rather than N. Balancing sum_{k<K}(log k)c_R sqrt(N/k) against
N/log N still gives K of order N/log^4 N, whose exponent tends to one.
If that is right the exponent must RISE with N along a fixed
primorial-like radical, and the failure at 1.6e6 is the constants
winning at log N = 14.

So this sweeps a second family: N = 30030 * 2^j, radical 3*5*7*11*13
held fixed, over a factor 64 in N -- four times the lever the main
family has.

BACKS: Remark {#rem:primorialladder} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  The control: recomputing N = 1621620, which is in the arithmetic
      test set and shares this radical, reproduces the 0.4747 of
      results/audit_residue_arithmetic.txt to within 0.001.
  R2  The failure is not permanent: along the ladder the exponent
      rises, its least-squares slope against log N being positive.
  R3  And it rises at least as fast as the 2^a 5^b family's slope,
      read from results/audit_residue_level.txt, because the budget
      that has to be caught is growing more slowly here.
  R4  But it does not catch up within reach: no N in the ladder
      reaches 0.5. A refutation would be the good outcome -- it would
      put the barrier inside the computable range at a primorial-like
      radical.

REFUTATION RULE (fixed before the run)

  R1  REFUTED at 0.001, which would mean this is not the same
      measurement.
  R2  REFUTED if the slope is not positive. That is the one that
      matters: a flat or falling exponent over a factor 64 would say
      the primorial failure is structural and the conditional
      reduction cannot reach theta' > 1/2 at those N at all.
  R3  REFUTED if the slope is below the family's.
  R4  REFUTED if any N in the ladder reaches 0.5.

  All four gate.

  NULL: a coin arm on the identical deviations. R(N;k) =
  sum_m mu(m) delta(m,k) with delta = Lambda(N-mk) - beta w(m,k);
  replacing mu by 8 global sign vectors on the odd squarefree m, each
  held across all k as mu is, and carrying the result through the same
  crossing and the same fit, gives the slope a coin produces along the
  same ladder against the same fixed budget. That is what R2 and R3
  have to be read against: lab_residue_cancellation.py established
  that R's SIZE is exactly a coin's, so if its slope were also exactly
  a coin's the rise would be a fact about magnitudes and not about mu,
  and if it differs that is the thing to report.
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
OUT = os.path.join(ROOT, "results", "lab_primorial_ladder.txt")

BASE = 30030                       # 2*3*5*7*11*13
LADDER = [BASE * (1 << j) for j in range(7)]
CONTROL = 1_621_620
KCAP = 100_000
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


def read_control():
    p = os.path.join(ROOT, "results", "audit_residue_arithmetic.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("N            odd part               threshold")
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            continue
        if int(f[0]) == CONTROL:
            return float(f[4])
    return None


def read_family_slope():
    p = os.path.join(ROOT, "results", "audit_residue_level.txt")
    src = io.open(p, encoding="utf-8").read()
    return float(re.search(r"least-squares slope against log N = "
                           r"([+\-][\d.]+)", src).group(1))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubctl = read_control()
    famslope = read_family_slope()
    say("read the control exponent %.4f at N = %d and the family slope"
        % (pubctl, CONTROL))
    say("  %+.6f, from results/audit_residue_arithmetic.txt and "
        "results/audit_residue_level.txt" % famslope)

    NS = LADDER + [CONTROL]
    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = sieves(NMAX)
    sqf = mu != 0
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]

    oddsqf = np.zeros(NMAX + 1, dtype=bool)
    oddsqf[1::2] = True
    oddsqf &= sqf
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
    for N in NS:
        PN = factor_set(N)
        A_, S_ = artin, twin
        for q in sorted(PN):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))
        thrc = S_ * (1.0 - A_)

        ks, Hs, Ps, Cs = [], [], [], []
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
            Cs.append((ms, lam[vals], w))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        P = np.array(Ps)
        lw = np.log(ks.astype(float))
        beta = float((H * P).sum() / (P * P).sum())
        aR = np.abs(H - beta * P)
        aC = np.empty((COINS, ks.size))
        for i, (ms, L, w) in enumerate(Cs):
            d = L - beta * w
            aC[:, i] = np.abs(coinmat[:, ms].astype(np.float64) @ d)
        del Cs
        cum = np.cumsum(lw * aR)
        j = int(np.searchsorted(cum, thrc * N))
        kstar = int(ks[j]) if j < ks.size else None
        e = math.log(kstar) / math.log(N) if kstar else float("nan")
        ec = []
        for c in range(COINS):
            cc = np.cumsum(lw * aC[c])
            jc = int(np.searchsorted(cc, thrc * N))
            kc = int(ks[jc]) if jc < ks.size else None
            ec.append(math.log(kc) / math.log(N) if kc
                      else float("nan"))
        rows.append((N, thrc, S_, A_, ks.size, beta, kstar, e, ec))
        say("  N = %-9d thr %.6f  S %.4f  1-A %.6f  #k %-6d K*_R "
            "%-6s exp %.4f" % (N, thrc, S_, 1.0 - A_, ks.size,
                               str(kstar), e))

    lad = [r for r in rows if r[0] != CONTROL]
    ctl = [r for r in rows if r[0] == CONTROL][0]

    # ------------------------------------------------------------- R1
    say()
    say("R1  the control, a different N of the same radical")
    d = abs(ctl[7] - pubctl)
    r1 = d < 0.001
    say("  N = %d: exponent %.4f against the published %.4f, diff %.5f"
        % (CONTROL, ctl[7], pubctl, d))
    say("  R1 %s" % ("hold" if r1 else "REFUTED"))

    # ---------------------------------------------------------- R2/R3
    say()
    say("R2/R3  the ladder: N = %d * 2^j, radical held fixed" % BASE)
    say("  N            threshold  K*_R    exponent   log2(N/%d)" % BASE)
    for N, thrc, S_, A_, nk, beta, kstar, e, ec in lad:
        say("  %-12d %-10.6f %-7s %-10.4f %d"
            % (N, thrc, str(kstar), e, int(round(math.log2(N / BASE)))))
    x = np.log(np.array([r[0] for r in lad], dtype=float))
    y = np.array([r[7] for r in lad])
    sl = float(np.polyfit(x, y, 1)[0])
    rr = float(np.corrcoef(x, y)[0, 1])
    r2 = sl > 0.0
    r3 = sl >= famslope
    say("  slope against log N %+.6f, correlation %.5f" % (sl, rr))
    say("  R2 slope positive   %s" % ("hold" if r2 else "REFUTED"))
    say("  R3 slope at least the family's %+.6f   %s"
        % (famslope, "hold" if r3 else "REFUTED"))
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    say("  leave-one-out on the slope: %+.6f, %+.6f, %+.6f -- "
        "spread %.6f" % (f[0], f[1], f[2], max(f) - min(f)))

    say()
    say("  THE NULL. The same ladder with mu replaced by a coin on the")
    say("  identical deviations, %d global sign vectors:" % COINS)
    say("  who        exponents along the ladder                slope")
    say("  mu         %s %+.6f"
        % (" ".join("%.4f" % v for v in y), sl))
    cs = []
    for c in range(COINS):
        yc = np.array([r[8][c] for r in lad])
        sc = float(np.polyfit(x, yc, 1)[0])
        cs.append(sc)
        say("  coin %-6d %s %+.6f"
            % (c, " ".join("%.4f" % v for v in yc), sc))
    cs.sort()
    inband = cs[0] <= sl <= cs[-1]
    say("  coin slopes span %+.6f to %+.6f, median %+.6f; mu is %s"
        % (cs[0], cs[-1], float(np.median(cs)),
           "inside" if inband else "OUTSIDE"))
    say("  A coin rises too, because the budget is fixed while the")
    say("  magnitudes grow: the rise is a fact about sizes, and")
    say("  lab_residue_cancellation.py already established that R's")
    say("  sizes are exactly a coin's. What R3 compares is this slope")
    say("  with the same slope at the OTHER radical, and there the")
    say("  budget differs by five, which is the whole content.")

    # ------------------------------------------------------------- R4
    say()
    r4 = all(r[7] < 0.5 for r in lad)
    say("R4  does the ladder reach the barrier? max exponent %.4f   %s"
        % (max(r[7] for r in lad), "hold" if r4 else "REFUTED"))

    say()
    say("  the arithmetic and the budget, declared:")
    rads = set()
    for N, thrc, S_, A_, nk, beta, kstar, e, ec in rows:
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  %d N, %d distinct odd radical%s: %s"
        % (len(rows), len(rads), "" if len(rads) == 1 else "s",
           ", ".join(str(v) for v in sorted(rads))))
    say("RADICALS %d" % len(rads))
    for N, thrc, S_, A_, nk, beta, kstar, e, ec in rows:
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, thrc))

    say()
    say("  DIAGNOSTIC (post hoc). Why the budget behaves differently")
    say("  here. Along this ladder the radical is fixed, so S(N) and")
    say("  A(N) are CONSTANT and the threshold does not move at all:")
    say("  N            threshold  1 - A(N)")
    for N, thrc, S_, A_, nk, beta, kstar, e, ec in lad:
        say("  %-12d %-10.6f %.6f" % (N, thrc, 1.0 - A_))
    say("  So the rise in the exponent along the ladder is the level")
    say("  alone, with the budget held fixed -- which is the cleanest")
    say("  form of the question, and not the same experiment as the")
    say("  family sweep, where S(N)(1-A(N)) is also constant but at a")
    say("  five times larger value.")
    say()
    say("  What it would take to reach the barrier. The drift of the")
    say("  slope, and -- the thing an earlier version of this forecast")
    say("  left out and audit_primorial_reach.py then refuted it for --")
    say("  the scatter of the rungs about their own line:")
    aa, bb = np.polyfit(x, y, 1)
    resid = y - (aa * x + bb)
    rms = float(np.sqrt((resid ** 2).mean()))
    dr = (max(f) - min(f)) / abs(sl) if sl else float("nan")
    say("  N            exponent   fitted     residual")
    for i, r in enumerate(lad):
        say("  %-12d %-10.4f %-10.4f %+.4f"
            % (r[0], y[i], aa * x[i] + bb, resid[i]))
    say("  r.m.s. residual %.4f against a trend of %.4f per doubling,"
        % (rms, aa * math.log(2.0)))
    say("  i.e. %.1f times what the trend gains in one rung"
        % (rms / (aa * math.log(2.0))))
    say("DRIFT primorial_slope %.4f" % dr)
    say("SCATTER primorial_slope %.4f" % rms)
    say()
    say("  So the forecast carries both. The line reaches 0.5 at:")
    u0 = (0.5 - bb) / aa
    say("    log10 N = %.4f" % (u0 / math.log(10.0)))
    say("  displaced by one r.m.s. residual either way:")
    ulo = (0.5 - rms - bb) / aa
    uhi = (0.5 + rms - bb) / aa
    say("    [%.4f, %.4f]" % (ulo / math.log(10.0),
                             uhi / math.log(10.0)))
    say("  and with the slope at its leave-one-out extremes, for")
    say("  comparison -- this is the part that was published alone and")
    say("  turned out not to be the binding uncertainty:")
    for ff in (min(f), max(f)):
        if ff > 0:
            say("    slope %+.6f -> log10 N = %.4f"
                % (ff, ((0.5 - (y[-1] - ff * x[-1])) / ff)
                   / math.log(10.0)))
        else:
            say("    slope %+.6f -> never" % ff)
    say()
    say("  THE SHAPES. The bracket so far is over the line's own")
    say("  parameters and scatter; it says nothing about whether the")
    say("  line is the right form. audit_ladder_shape.py showed that")
    say("  eleven rungs do not choose the form, and here there are")
    say("  seven, with the crossing still outside the data -- the case")
    say("  where the choice should cost most. The same five shapes,")
    say("  fitted to these seven:")

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
    dof = x.size - 2
    se = best / math.sqrt(2.0 * dof)
    say("  shape                        r.m.s.    (r-best)/s.e.  "
        "0.50 at")
    keep = []
    for nm in NAMES:
        rr2 = fit2[nm][1]
        t = (rr2 - best) / se
        cc = crossat(nm, 0.5)
        if t <= 1.0:
            keep.append(cc)
        say("  %-28s %-9.5f %-14.2f %s"
            % (LAB[nm], rr2, t,
               "%.4f" % cc if cc is not None else "never"))
    say("  %d shapes fitted, %d within one standard error of the best"
        % (len(NAMES), len(keep)))
    say("SHAPES %d" % len(NAMES))
    vals = [v for v in keep if v is not None]
    lo3 = min(min(vals), ulo / math.log(10.0))
    hi3 = max(max(vals), uhi / math.log(10.0))
    say("  the surviving shapes put 0.50 at %s -- a spread of %.4f"
        % (", ".join("%.4f" % v for v in sorted(vals)),
           max(vals) - min(vals)))
    say("  decades")
    say("  the line's scatter bracket alone is [%.4f, %.4f]; adding"
        % (ulo / math.log(10.0), uhi / math.log(10.0)))
    say("  the shapes gives [%.4f, %.4f], wider by %.4f decades"
        % (lo3, hi3, (hi3 - lo3)
           - (uhi - ulo) / math.log(10.0)))
    say("BRACKET log10_N_primorial_reaches_half %.4f %.4f %.4f"
        % (u0 / math.log(10.0), lo3, hi3))

    say()
    say("=" * 70)
    ok = r1 and r2 and r3 and r4
    say("the primorial failure is finite-N: the exponent rises, and "
        "faster than the easy family's" if ok else "REFUTED")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           along N = 30030 * 2^j for j = 0..6 -- a ladder with",
        "           the odd radical 3*5*7*11*13 held fixed over a factor",
        "           64 in N -- together with N = 1621620 as a control of",
        "           the same radical; the least-squares slope of that",
        "           exponent against log N and its leave-one-out spread.",
        "NULL: a coin arm on the identical deviations. R = sum_m mu(m)",
        "      delta(m,k) with delta = Lambda(N-mk) - beta w(m,k);",
        "      replacing mu by 8 global sign vectors on the odd",
        "      squarefree m, each held across all k as mu is, and",
        "      carrying the result through the same crossing and the",
        "      same fit, gives the slope a coin produces along the same",
        "      ladder against the same fixed budget. That is the scale",
        "      R2 and R3 are read against.",
        "FIELD: N = 30030 * 2^j, j = 0..6, and N = 1621620; k squarefree",
        "       and coprime to N with 2 <= k < 100000; m odd, squarefree,",
        "       coprime to the odd part of k, m < N/k; the sieve weight",
        "       uses the odd primes up to 30; beta refitted as",
        "       sum(H P)/sum(P^2) on the same k-range; S(N) and A(N)",
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
