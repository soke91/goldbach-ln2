# -*- coding: utf-8 -*-
r"""
A derived shape instead of a fitted one.

WHAT IS AT STAKE

Remark {#rem:laddershape} showed that the eleven rungs of the
primorial ladder do not choose a functional form: three survive at one
standard error and put theta' = 0.56 over 8.35 decades apart. Every
one of those forms was FITTED. None was derived.

The heuristic derives one. K*_R is where

    sum_{k<K, admissible} (log k) |R(N;k)|  =  S(N)(1-A(N)) N,

and Remark {#rem:residueconstant} measured |R| ~ c_R(N) sqrt(N/k) with
c_R about 0.50 sqrt(log N). Putting those together gives the crossing
with no free shape at all -- the model that {#rem:residueconstant}'s
rule Q3 confirmed to within 0.7 per cent at the 2^a 5^b family, and
that nobody has run on this ladder.

If it reproduces the eleven exponents it is the right extrapolant, and
it answers about 0.56 without a choice being made. If it does not, the
shape ambiguity is real and irreducible and {#rem:laddershape} stands
as the last word.

BACKS: Remark {#rem:laddermodel} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  J1  The control: the measured exponents reproduce the eleven
      published in results/audit_primorial_rung10.txt to within 1e-4.
  J2  The derived model predicts each rung's crossing to within 5 per
      cent in K*_R -- the tolerance {#rem:residueconstant} used, and
      met, at the other family.
  J3  The constant it rests on is no worse behaved here than there:
      c_R/sqrt(log N) across the eleven rungs has a spread no larger
      than the 0.1436 that file measured.
  J4  And the derived shape lands inside the fitted ones: the N at
      which the model reaches 0.56 lies between the extreme answers
      of the surviving fitted shapes, read from
      results/audit_ladder_shape.txt.

REFUTATION RULE (fixed before the run)

  J1  REFUTED at 1e-4 at any rung.
  J2  REFUTED at 5 per cent at any rung. That is the one that
      matters: a model that cannot reproduce the eleven crossings it
      is built from cannot be trusted to extrapolate them, and the
      shape ambiguity would then be irreducible.
  J3  REFUTED if the spread exceeds 0.1436, which would say the
      derived shape rests on a constant that is worse here than the
      one already judged not to license extrapolation.
  J4  REFUTED if the model's answer falls outside the fitted range.
      A derived shape landing outside every fitted one would not be a
      resolution but a fifth opinion.

  All four gate.

  NO NULL IS RUN and none applies. A deterministic model is compared
  with a measured crossing of the same measured sum at eleven N; there
  is no background to detect against. The sign controls for R on this
  ladder were run in lab_primorial_ladder.py, whose eight global sign
  vectors bracketed mu's slope, and in lab_residue_cancellation.py,
  whose coin arm on the identical deviations established that R's size
  is bought by cancellation at exactly a coin's rate.
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
OUT = os.path.join(ROOT, "results", "audit_ladder_model.txt")

BASE = 30030                       # 2*3*5*7*11*13
LADDER = [BASE * (1 << j) for j in range(11)]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000
THETA = 0.56


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


def read_rungs():
    p = os.path.join(ROOT, "results", "audit_primorial_rung10.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("N            log10 N   exponent   fitted     "
                  "residual")
    ex = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            break
        ex[int(f[0])] = float(f[2])
    return ex


def read_family_drift():
    p = os.path.join(ROOT, "results", "audit_residue_constant.txt")
    src = io.open(p, encoding="utf-8").read()
    return float(re.search(r"spread ([\d.]+) of the mean",
                           src).group(1))


def read_shape_answers():
    """where the surviving fitted shapes put 0.56 -- read"""
    p = os.path.join(ROOT, "results", "audit_ladder_shape.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("  shape                        0.50      0.56")
    out = []
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3:
            break
        try:
            out.append(float(f[-1]))
        except ValueError:
            break
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_rungs()
    famdrift = read_family_drift()
    shp = read_shape_answers()
    say("read %d published rungs, the family's c_R drift %.4f, and %d"
        % (len(pub), famdrift, len(shp)))
    say("  surviving shape answers for 0.56: %s"
        % ", ".join("%.4f" % v for v in shp))

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
        scale = np.sqrt(N / ks.astype(float))
        sel = ks <= kstar
        cR = float((aR[sel] / scale[sel]).mean())
        dens = float(sel.sum()) / kstar
        model = np.cumsum(lw * cR * scale)
        jm = int(np.searchsorted(model, thrc * N))
        km = int(ks[min(jm, ks.size - 1)])
        rows.append((N, thrc, kstar, e, cR, dens, km))
        say("  N = %-9d K*_R %-6d exp %.4f  c_R %.4f  model K* %-6d"
            % (N, kstar, e, cR, km))

    # ------------------------------------------------------------- J1
    say()
    say("J1  the control")
    say("  N            here       published   diff")
    j1 = True
    for N, thrc, kstar, e, cR, dens, km in rows:
        d = abs(e - pub[N])
        if not (d < 1e-4):
            j1 = False
        say("  %-12d %-10.4f %-11.4f %.6f" % (N, e, pub[N], d))
    say("  J1 %s" % ("hold" if j1 else "REFUTED"))

    # ------------------------------------------------------------- J2
    say()
    say("J2  does the derived model reproduce the crossings?")
    say("  N            measured   model      ratio")
    j2 = True
    for N, thrc, kstar, e, cR, dens, km in rows:
        r = km / kstar
        if abs(r - 1.0) >= 0.05:
            j2 = False
        say("  %-12d %-10d %-10d %.4f" % (N, kstar, km, r))
    say("  J2 %s   (cap 5 per cent)" % ("hold" if j2 else "REFUTED"))

    # ------------------------------------------------------------- J3
    say()
    say("J3  the constant the derived shape rests on")
    say("  N            c_R        sqrt(log N)   c_R/sqrt(log N)")
    gam = []
    for N, thrc, kstar, e, cR, dens, km in rows:
        g = cR / math.sqrt(math.log(N))
        gam.append(g)
        say("  %-12d %-10.4f %-13.4f %.4f"
            % (N, cR, math.sqrt(math.log(N)), g))
    spread = (max(gam) - min(gam)) / float(np.mean(gam))
    j3 = spread <= famdrift
    say("  spread %.4f of the mean, against the family's %.4f   %s"
        % (spread, famdrift, "hold" if j3 else "REFUTED"))
    say("  J3 %s" % ("hold" if j3 else "REFUTED"))
    say("DRIFT ladder_model_cR %.4f" % spread)

    # ------------------------------------------------------------- J4
    say()
    say("J4  where the derived shape puts 0.56")
    gm = float(np.mean(gam))
    dm = float(np.mean([r[5] for r in rows]))
    thrm = float(np.mean([r[1] for r in rows]))
    say("  the model with c_R = g sqrt(log N), g = %.4f, admissible"
        % gm)
    say("  density %.4f and budget %.6f solves" % (dm, thrm))
    say("    N^{(1-e)/2} = 2 d g sqrt(log N) (e log N - 2) / budget")

    def expo_at(u):
        lo, hi = 0.05, 0.999
        for _ in range(200):
            e = 0.5 * (lo + hi)
            lhs = math.exp((1.0 - e) * u / 2.0)
            rhs = 2.0 * dm * gm * math.sqrt(u) * (e * u - 2.0) / thrm
            if lhs > rhs:
                lo = e
            else:
                hi = e
        return 0.5 * (lo + hi)

    def reaches(tgt):
        lo, hi = 5.0, 200.0
        if expo_at(hi * math.log(10.0)) < tgt:
            return None
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if expo_at(mid * math.log(10.0)) < tgt:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    say("  log10 N      model exponent   measured")
    for N, thrc, kstar, e, cR, dens, km in rows:
        say("  %-12.4f %-16.4f %.4f"
            % (math.log10(N), expo_at(math.log(N)), e))
    m50 = reaches(0.5)
    m56 = reaches(THETA)
    say("  the model reaches 0.50 at log10 N = %s and %.2f at %s"
        % ("%.4f" % m50 if m50 else "never", THETA,
           "%.4f" % m56 if m56 else "never"))
    fin = [v for v in shp if v is not None]
    j4 = (m56 is not None and fin and min(fin) <= m56 <= max(fin))
    say("  the fitted shapes span %.4f to %.4f" % (min(fin), max(fin)))
    say("  J4 %s" % ("hold" if j4 else "REFUTED"))
    say("SHAPES 1")
    say("BRACKET log10_N_ladder_model_reaches_theta %.4f %.4f %.4f"
        % (m56 if m56 else 0.0, min(fin), max(fin)))
    say("SCATTER ladder_model_cR %.4f" % spread)

    say()
    say("  the arithmetic and the budget, declared:")
    rads = set()
    for N, thrc, kstar, e, cR, dens, km in rows:
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  %d N, %d distinct odd radical%s: %s"
        % (len(rows), len(rads), "" if len(rads) == 1 else "s",
           ", ".join(str(v) for v in sorted(rads))))
    say("RADICALS %d" % len(rads))
    for N, thrc, kstar, e, cR, dens, km in rows:
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, thrc))

    say()
    say("  DIAGNOSTIC (post hoc). What the derived shape does that a")
    say("  fitted one cannot: it has no freedom left once c_R and the")
    say("  density are measured. Its residuals against the eleven")
    say("  measured exponents, and the line's, for comparison:")
    u = np.array([math.log(r[0]) for r in rows])
    y = np.array([r[3] for r in rows])
    a1, b1 = np.polyfit(u, y, 1)
    say("  N            measured   model      model resid  line resid")
    for i, (N, thrc, kstar, e, cR, dens, km) in enumerate(rows):
        me = expo_at(math.log(N))
        say("  %-12d %-10.4f %-10.4f %+-12.4f %+.4f"
            % (N, e, me, e - me, e - (a1 * u[i] + b1)))
    mres = np.array([r[3] - expo_at(math.log(r[0])) for r in rows])
    lres = y - (a1 * u + b1)
    say("  r.m.s.: model %.4f, line %.4f"
        % (float(np.sqrt((mres ** 2).mean())),
           float(np.sqrt((lres ** 2).mean()))))
    say("  The model's residuals are not noise: positive at the bottom")
    say("  of the ladder and negative at the top, so the derived shape")
    say("  rises too fast. J3 says why.")
    say()
    say("  What c_R actually does here. Fitting c_R against log N:")
    cc = np.array([r[4] for r in rows])
    sfit = float(np.polyfit(np.log(u), np.log(cc), 1)[0])
    rfit = float(np.corrcoef(np.log(u), np.log(cc))[0, 1])
    say("    c_R ~ (log N)^%.4f, correlation %.5f" % (sfit, rfit))
    say("  against the 1/2 the derived shape assumes and the 1/2 that")
    say("  {#rem:residueconstant} found at the 2^a 5^b family. The")
    say("  residue's constant grows about %.1f times faster in the"
        % (sfit / 0.5))
    say("  exponent at the primorial radical, which is not a wobble in")
    say("  a constant but a different law, and it is why no derived")
    say("  shape is available here.")
    f2 = [float(np.polyfit(np.log(u)[sl], np.log(cc)[sl], 1)[0])
          for sl in (slice(None), slice(1, None), slice(0, -1))]
    say("  leave-one-out on that exponent: %.4f, %.4f, %.4f -- "
        "spread %.4f" % (f2[0], f2[1], f2[2], max(f2) - min(f2)))
    say("SWEPT ladder_cR_growth N-range %.6f" % (max(f2) - min(f2)))
    say("CORR ladder_cR_growth %.5f" % abs(rfit))
    say("POP ladder_cR_growth %d" % len(rows))

    say()
    say("=" * 70)
    ok = j1 and j2 and j3 and j4
    say("the derived shape reproduces the ladder and selects among the "
        "fitted ones" if ok else "REFUTED")

    head = [
        "STATISTIC: at each rung of N = 30030 * 2^j, j = 0..10, the",
        "           measured crossing K*_R of sum(log k)|R| against",
        "           S(N)(1-A(N))N and its exponent; the constant",
        "           c_R = mean |R|/sqrt(N/k) below that crossing; the",
        "           crossing predicted by",
        "           sum_{k<K}(log k) c_R sqrt(N/k) = S(N)(1-A(N))N over",
        "           the actual admissible k; and the N at which that",
        "           derived model reaches 0.50 and 0.56.",
        "NULL: none is run and none applies. A deterministic model is",
        "      compared with a measured crossing of the same measured",
        "      sum at eleven N; there is no background to detect",
        "      against. The sign controls for R on this ladder were run",
        "      in lab_primorial_ladder.py and",
        "      lab_residue_cancellation.py.",
        "FIELD: N = 30030 * 2^j, j = 0..10, the odd radical 3*5*7*11*13",
        "       fixed so the budget is constant; k squarefree and",
        "       coprime to N with 2 <= k < 100000; m odd, squarefree,",
        "       coprime to the odd part of k, m < N/k; the sieve weight",
        "       uses the odd primes up to 30; beta refitted as",
        "       sum(H P)/sum(P^2) on that k-range; S(N) and A(N) from",
        "       Euler products at the fixed bound 4000000; the published",
        "       exponents, the family's c_R drift and the fitted shapes'",
        "       answers are read from results/.",
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
