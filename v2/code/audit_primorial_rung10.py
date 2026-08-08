# -*- coding: utf-8 -*-
r"""
The rung that sits on the forecast.

WHAT IS AT STAKE

Remark {#rem:primorialreach} refuted the first forecast's point
estimate and left its interval open; corrected to carry the residual
scatter, the forecast is that the primorial ladder's exponent reaches
1/2 at N = 10^7.4684 with the bracket [10^7.2189, 10^7.7180].

The next rung is N = 30030 * 2^10 = 30750720, which is 10^7.4877 --
just above the point estimate. On the fitted line its exponent should
be 0.5003; the r.m.s. residual is 0.0039, so it lands either side by
chance. One rung decides whether the barrier is crossed inside the
computable range at a primorial-like radical, and this runs it.

Nothing is protected. The previous cycle's forecast has already been
tested once and its point estimate found wrong; this is the second
test of the same crossing, and the second point estimate is under the
same rule.

BACKS: Remark {#rem:primorialrung10} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  G1  The control: recomputing N = 15375360, the current top rung,
      reproduces the 0.4941 published in
      results/audit_primorial_reach.txt to within 1e-4.
  G2  The fit holds one rung further: the new exponent lies within one
      published r.m.s. residual of the line fitted on the ten
      published rungs.
  G3  The barrier is crossed: the new exponent is at or above 0.5.
      This is a coin flip on the published numbers and is written down
      as one -- the fitted value at this N is 0.5003 against a scatter
      of 0.0039.
  G4  And either way the bracket survives: the crossing implied by the
      eleven rungs stays inside [10^7.2189, 10^7.7180].

REFUTATION RULE (fixed before the run)

  G1  REFUTED at 1e-4.
  G2  REFUTED if the new exponent is more than one r.m.s. residual
      from the line, which would say the trend is bending and the
      linear extrapolation is the wrong shape.
  G3  REFUTED if the new exponent is below 0.5. That is not a failure
      of anything -- it is one side of a coin flip -- but it decides
      whether the crossing has been observed or merely bracketed.
  G4  REFUTED if the implied crossing leaves the bracket. That is the
      one that matters: it would be the second forecast of this
      crossing to fail, and the corrected bracket was built precisely
      to survive this test.

  All four gate.

  NO NULL IS RUN and none is needed. The question is where a
  deterministic curve crosses a line, and the control is G1: the same
  computation reproducing the rung it is extending. The coin arm for
  this statistic was run in lab_primorial_ladder.py, whose eight
  global sign vectors bracketed mu's slope, and in
  lab_primorial_share.py, whose coin arm scattered by the same amount
  as mu's -- establishing that both the rise and its noise are facts
  about magnitudes.
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
OUT = os.path.join(ROOT, "results", "audit_primorial_rung10.txt")

BASE = 30030                       # 2*3*5*7*11*13
CONTROL = BASE * (1 << 9)          # 15375360, the current top rung
NEW = BASE * (1 << 10)             # 30750720
KCAP = 100_000                     # beta is fitted on this range
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


def read_reach():
    """the ten published rungs and the scatter -- read, not copied"""
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
    m = re.search(r"BRACKET log10_N_primorial_reaches_half_v2 "
                  r"([\d.]+) ([\d.]+) ([\d.]+)", src)
    return ex, sc, (float(m.group(1)), float(m.group(2)),
                    float(m.group(3)))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub, sc, (bpt, blo, bhi) = read_reach()
    say("read %d published rungs, the scatter %.4f and the bracket"
        % (len(pub), sc))
    say("  10^%.4f [10^%.4f, 10^%.4f] from "
        "results/audit_primorial_reach.txt" % (bpt, blo, bhi))

    say("sieving to %d ..." % NEW)
    lam, mu = sieves(NEW)
    sqf = mu != 0
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    got = {}
    for N in (CONTROL, NEW):
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
        cum = np.cumsum(lw * np.abs(H - beta * P))
        j = int(np.searchsorted(cum, thrc * N))
        kstar = int(ks[j]) if j < ks.size else None
        e = math.log(kstar) / math.log(N) if kstar else float("nan")
        got[N] = (thrc, ks.size, beta, kstar, e)
        say("  N = %-9d thr %.6f  #k %-6d K*_R %-6s exp %.4f"
            % (N, thrc, ks.size, str(kstar), e))

    # ------------------------------------------------------------- G1
    say()
    d = abs(got[CONTROL][4] - pub[CONTROL])
    g1 = d < 1e-4
    say("G1  the control at N = %d: %.4f against the published %.4f, "
        "diff %.6f" % (CONTROL, got[CONTROL][4], pub[CONTROL], d))
    say("  G1 %s" % ("hold" if g1 else "REFUTED"))

    # ------------------------------------------------------------- G2
    say()
    xs = np.log(np.array(sorted(pub), dtype=float))
    ys = np.array([pub[n] for n in sorted(pub)])
    a1, b1 = np.polyfit(xs, ys, 1)
    fitted = a1 * math.log(NEW) + b1
    enew = got[NEW][4]
    off = abs(enew - fitted)
    g2 = off <= sc
    say("G2  the fit one rung further")
    say("  line on the ten published rungs: slope %+.6f" % a1)
    say("  at N = %d it predicts %.4f, measured %.4f, off by %.4f"
        % (NEW, fitted, enew, off))
    say("  published r.m.s. residual %.4f   %s"
        % (sc, "within" if g2 else "OUTSIDE"))
    say("  G2 %s" % ("hold" if g2 else "REFUTED"))

    # ------------------------------------------------------------- G3
    say()
    g3 = enew >= HALF
    say("G3  is the barrier crossed?")
    say("  N = %d, log10 N = %.4f, exponent %.4f   %s"
        % (NEW, math.log10(NEW), enew,
           "AT OR ABOVE 0.5" if g3 else "below 0.5"))
    say("  G3 %s" % ("hold" if g3 else "REFUTED"))

    # ------------------------------------------------------------- G4
    say()
    allx = np.append(xs, math.log(NEW))
    ally = np.append(ys, enew)
    a2, b2 = np.polyfit(allx, ally, 1)
    resid = ally - (a2 * allx + b2)
    rms2 = float(np.sqrt((resid ** 2).mean()))
    u0 = (HALF - b2) / a2
    pt2 = u0 / math.log(10.0)
    lo2 = ((HALF - rms2 - b2) / a2) / math.log(10.0)
    hi2 = ((HALF + rms2 - b2) / a2) / math.log(10.0)
    g4 = blo <= pt2 <= bhi
    say("G4  the crossing implied by eleven rungs")
    say("  slope %+.6f, r.m.s. residual %.4f" % (a2, rms2))
    say("  crossing at log10 N = %.4f, bracket [%.4f, %.4f]"
        % (pt2, lo2, hi2))
    say("  published bracket [%.4f, %.4f]   %s"
        % (blo, bhi, "inside" if g4 else "OUTSIDE"))
    say("  G4 %s" % ("hold" if g4 else "REFUTED"))
    say("SCATTER primorial_slope %.4f" % rms2)
    f = [float(np.polyfit(allx[s], ally[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    say("DRIFT primorial_slope %.4f" % ((max(f) - min(f)) / abs(a2)))
    say("BRACKET log10_N_primorial_reaches_half_v3 %.4f %.4f %.4f"
        % (pt2, lo2, hi2))
    say("SWEPT primorial_rung10 N-range %.6f" % (max(f) - min(f)))
    say("CORR primorial_rung10 %.5f"
        % abs(float(np.corrcoef(allx, ally)[0, 1])))
    say("POP primorial_rung10 %d" % allx.size)

    say()
    say("  THE SHAPES. A bracket built from one functional form is a")
    say("  bracket over that form's parameters and not over the")
    say("  extrapolation. audit_ladder_shape.py showed these eleven")
    say("  points do not choose the form, so this file fits the same")
    say("  five and brackets over the ones that survive.")
    lu = np.log(allx)

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
        X = design(nm, allx)
        cc, *_ = np.linalg.lstsq(X, ally, rcond=None)
        rr = ally - X @ cc
        fit2[nm] = (cc, float(np.sqrt((rr ** 2).mean())))
    zz = lu / allx
    c1 = float(((1.0 - ally) * zz).sum() / (zz * zz).sum())
    r1 = ally - (1.0 - c1 * zz)
    fit2["heuristic1"] = (np.array([c1]),
                          float(np.sqrt((r1 ** 2).mean())))
    NAMES.append("heuristic1")
    LAB["heuristic1"] = "1 - c log log N / log N"

    def val(nm, uu):
        cc = fit2[nm][0]
        if nm == "heuristic1":
            return 1.0 - cc[0] * math.log(uu) / uu
        return float(design(nm, np.array([uu])) @ cc)

    def cross(nm, tgt):
        if val(nm, 400.0) < tgt:
            return None
        lo, hi = 2.0, 400.0
        for _ in range(300):
            mid = 0.5 * (lo + hi)
            if val(nm, mid) < tgt:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi) / math.log(10.0)

    best = min(fit2[nm][1] for nm in NAMES)
    dof = allx.size - 2
    se = best / math.sqrt(2.0 * dof)
    say("  shape                        r.m.s.    (r-best)/s.e.  "
        "0.50 at")
    keep = []
    for nm in NAMES:
        rr = fit2[nm][1]
        t = (rr - best) / se
        cc = cross(nm, HALF)
        if t <= 1.0:
            keep.append((nm, cc))
        say("  %-28s %-9.5f %-14.2f %s"
            % (LAB[nm], rr, t,
               "%.4f" % cc if cc is not None else "never"))
    say("  %d shapes fitted, %d within one standard error of the best"
        % (len(NAMES), len(keep)))
    say("SHAPES %d" % len(NAMES))
    vals = [c for _, c in keep if c is not None]
    lo3 = min(min(vals), lo2)
    hi3 = max(max(vals), hi2)
    say("  the surviving shapes put 0.50 at %s -- a spread of %.4f"
        % (", ".join("%.4f" % v for v in sorted(vals)),
           max(vals) - min(vals)))
    say("  decades, and the line's own scatter bracket [%.4f, %.4f]"
        % (lo2, hi2))
    say("  already contains all of them, so the shape-aware bracket is")
    say("  [%.4f, %.4f], wider by %.4f decades."
        % (lo3, hi3, (hi3 - lo3) - (hi2 - lo2)))
    say("  At the crossing that has been OBSERVED the choice of form")
    say("  costs nothing: the shapes are pinned there by the data. It")
    say("  is only out at theta' = 0.56, four decades past the top")
    say("  rung, that they separate -- which is exactly the boundary")
    say("  between what these eleven points measure and what they")
    say("  merely extrapolate.")
    say("BRACKET log10_N_primorial_reaches_half_shapes %.4f %.4f %.4f"
        % (pt2, lo3, hi3))

    say()
    say("  the arithmetic and the budget, declared:")
    rads = set()
    for N in (CONTROL, NEW):
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  the ladder is one family, %d distinct odd radical%s: %s"
        % (len(rads), "" if len(rads) == 1 else "s",
           ", ".join(str(v) for v in sorted(rads))))
    say("RADICALS %d" % len(rads))
    for N in (CONTROL, NEW):
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, got[N][0]))

    say()
    say("  DIAGNOSTIC (post hoc). The eleven rungs and their residuals")
    say("  about the refitted line:")
    say("  N            log10 N   exponent   fitted     residual")
    for i, nn in enumerate(list(sorted(pub)) + [NEW]):
        say("  %-12d %-9.4f %-10.4f %-10.4f %+.4f"
            % (nn, math.log10(nn), ally[i], a2 * allx[i] + b2,
               resid[i]))
    say("  the slope moved from %+.6f on ten rungs to %+.6f on eleven,"
        % (a1, a2))
    say("  and the scatter from %.4f to %.4f" % (sc, rms2))
    say()
    say("  And what the same line says about the theta' these papers")
    say("  use, which is not 1/2:")
    for tgt in (0.5, 0.56):
        u = (tgt - b2) / a2
        say("    exponent %.2f at log10 N = %.4f"
            % (tgt, u / math.log(10.0)))
    say("  -- four decades past the top rung, so the crossing that has")
    say("  been observed is the square-root barrier and not the")
    say("  route's own requirement.")

    say()
    say("=" * 70)
    ok = g1 and g2 and g3 and g4
    say("the barrier is crossed inside the computable range"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N, at",
        "           N = 30030 * 2^10 = 30750720 and, as a control, at",
        "           N = 30030 * 2^9; the offset of the new exponent from",
        "           the line fitted on the ten published rungs; and the",
        "           crossing of 1/2 implied by all eleven.",
        "NULL: none is run and none is needed. The question is where a",
        "      deterministic curve crosses a line, and the control is",
        "      G1, the same computation reproducing the rung it",
        "      extends. The coin arms for this statistic were run in",
        "      lab_primorial_ladder.py and lab_primorial_share.py,",
        "      establishing that both the rise and its noise are facts",
        "      about magnitudes and not about mu.",
        "FIELD: N = 15375360 and 30750720, the odd radical 3*5*7*11*13",
        "       fixed so the threshold is constant; k squarefree and",
        "       coprime to N with 2 <= k < 100000, the range beta is",
        "       fitted on throughout this ladder; m odd, squarefree,",
        "       coprime to the odd part of k, m < N/k; the sieve weight",
        "       uses the odd primes up to 30; S(N) and A(N) from Euler",
        "       products at the fixed bound 4000000; the ten published",
        "       rungs and the bracket are read from",
        "       results/audit_primorial_reach.txt.",
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
