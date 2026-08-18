# -*- coding: utf-8 -*-
r"""
The bracket on the lean forecast, and an independent measurement of
the lean itself.

WHAT IS AT STAKE

lab_lean_decay.py reports that the mass-weighted fraction f of k with
H(N;k) > 0 climbs towards one half, fits |0.5 - f| ~ N^{-0.1673}, and
quotes the N at which the fit reaches 0.01: 1.021e14. Its own
diagnostic already says the power law is not separable from
(log N)^{-c} over the accessible factor 32, and that the log law puts
the same event at 4.3e22.

Gate check G28 asks every forecast beyond the computed range to carry
a bracket. This supplies it, and does not take the measurement on
trust: f is recomputed by a different route. lab_lean_decay.py sums
H(N;k) = sum_{m < N/k, (m,k)=1} Lambda(N - mk) mu(m) directly; here it
is obtained from the dilation identity [eq:dilate], H(N;k) =
mu(k) A(N;k) with A(N;k) = sum_{n = N mod k} Lambda(n) mu(N-n), which
is a strided sum over a different array. If the two disagree, the
published lean is wrong and the forecast is moot.

The bracket has two independent sources and both are computed:
the fit's own parameter uncertainty, taken as the leave-one-out spread
the repository already sweeps, and the choice of law.

BACKS: Remark {#rem:leanbracket} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The dilation route reproduces the published mass-weighted f at
      every N to within 0.005 -- two routes, one number.
  V2  And so the refitted power-law exponent matches the published
      one to within 0.01.
  V3  The parameter bracket alone -- the forecast N refitted on each
      leave-one-out subset -- spans less than 4 decades. Estimated at
      roughly three and a half from the published exponent spread of
      0.0208 before running, so this threshold was not chosen blind.
  V4  The total bracket, parameter and law together, exceeds 3
      decades: the single number 1.021e14 cannot be quoted.

REFUTATION RULE (fixed before the run)

  V1  REFUTED if any N differs by 0.005 or more, which would mean one
      of the two routes to H is wrong.
  V2  REFUTED at 0.01.
  V3  REFUTED if the parameter bracket reaches 4 decades.
  V4  REFUTED if the total bracket is under 3 decades, in which case
      the published forecast is quotable as it stands.

  All four gate.

  NO NULL IS RUN here and none is needed: lab_lean_decay.py carries
  the coin arm as its reference level, and this script recomputes the
  same statistic on the same field by a second route rather than
  detecting anything against a background. What is added is arithmetic
  on a published fit.
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
OUT = os.path.join(ROOT, "results", "audit_lean_bracket.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
GOAL = 0.01


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


def read_published():
    """the published f and exponent and forecast -- read, not copied"""
    p = os.path.join(ROOT, "results", "lab_lean_decay.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("mass frac +")
    f = {}
    for ln in src[i:].splitlines()[2:]:
        g = ln.split()
        if len(g) < 4 or not g[0].isdigit():
            break
        f[int(g[0])] = float(g[2])
    b = float(re.search(r"\|0\.5 - f\| ~ N\^\{?-([\d.]+)", src).group(1))
    fc = float(re.search(r"\|0\.5 - f\| < 0\.01 first at N = ([\d.e+]+)",
                         src).group(1))
    lg = float(re.search(r"N = e\^([\d.]+) = ", src).group(1))
    return f, b, fc, lg


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubf, pubb, pubfc, publg = read_published()
    say("read from results/lab_lean_decay.txt: %d values of f, "
        "exponent %.4f," % (len(pubf), pubb))
    say("  forecast %.3e, and the log-law crossing at e^%.1f"
        % (pubfc, publg))

    NMAX = max(NS)
    say()
    say("sieving to %d ..." % NMAX)
    lam, mu = sieves(NMAX)
    sqf = mu != 0

    say()
    say("V1  the lean by the dilation route, H = mu(k) A(N;k)")
    say("  N            #k      f (dilation)   f (published)   diff")
    v1 = True
    fs = []
    for N in NS:
        PN = factor_set(N)
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if sqf[k] and all(k % q for q in PN)],
                      dtype=np.int64)
        f0 = np.zeros(N, dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        f0[1:] = lam[1:N] * mu[N - idx].astype(np.float64)
        del idx
        H = np.empty(ks.size)
        for i, k in enumerate(ks):
            kk = int(k)
            r = N % kk
            A = f0[r::kk].sum() if r else f0[kk::kk].sum()
            H[i] = mu[kk] * A
        del f0
        w = np.log(ks.astype(np.float64)) * np.abs(H)
        frac = float(w[H > 0].sum() / w.sum())
        fs.append(frac)
        d = abs(frac - pubf[N])
        if d >= 0.005:
            v1 = False
        say("  %-12d %-7d %-14.4f %-15.4f %.4f"
            % (N, ks.size, frac, pubf[N], d))
    say("  V1 %s" % ("hold" if v1 else "REFUTED"))

    # ------------------------------------------------------------- V2
    x = np.log(np.array(NS, dtype=float))
    y = np.log(np.array([abs(0.5 - v) for v in fs]))
    b_full = float(np.polyfit(x, y, 1)[0])
    say()
    say("V2  the power-law exponent refitted on the dilation numbers")
    v2 = abs(abs(b_full) - pubb) < 0.01
    say("  refitted %.4f against the published -%.4f   diff %.4f   %s"
        % (b_full, pubb, abs(abs(b_full) - pubb),
           "hold" if v2 else "REFUTED"))

    # both forecasts are returned as log10 N: the log law overflows a
    # double before it reaches 0.01 on some subsets
    L10 = math.log(10.0)

    def forecast_pow(xs, ys):
        s, c = np.polyfit(xs, ys, 1)
        return float((math.log(GOAL) - c) / s) / L10

    def forecast_log(xs, ys):
        """xs is already log N; the law is |0.5-f| ~ (log N)^{-c}"""
        s, c = np.polyfit(np.log(xs), ys, 1)
        return math.exp(float((math.log(GOAL) - c) / s)) / L10

    # ------------------------------------------------------------- V3
    say()
    say("V3  the parameter bracket: the same law refitted on each")
    say("  leave-one-out subset")
    say("  subset                   exponent    forecast N     log10")
    subs = [("all six N", slice(None)),
            ("without the smallest", slice(1, None)),
            ("without the largest", slice(0, -1))]
    pw = []
    for name, s in subs:
        e = float(np.polyfit(x[s], y[s], 1)[0])
        g = forecast_pow(x[s], y[s])
        pw.append(g)
        say("  %-24s %-11.4f 10^%-11.2f %.2f" % (name, e, g, g))
    pspan = max(pw) - min(pw)
    v3 = pspan < 4.0
    say("  parameter bracket %.2f decades   (cap 4)   %s"
        % (pspan, "hold" if v3 else "REFUTED"))

    # ------------------------------------------------------------- V4
    say()
    say("V4  the law bracket, and the total")
    lg = []
    for name, s in subs:
        g = forecast_log(x[s], y[s])
        lg.append(g)
        say("  (log N)^-c %-22s forecast N = 10^%.2f" % (name, g))
    allf = pw + lg
    tot = max(allf) - min(allf)
    v4 = tot > 3.0
    say("  total bracket %.2f decades   (floor 3)   %s"
        % (tot, "hold" if v4 else "REFUTED"))

    say()
    say("  Bracket lines, in the form the gate reads, in log10 N -- the")
    say("  log law overflows a double on some subsets. This file")
    say("  supplies them for results/lab_lean_decay.txt, whose forecast")
    say("  is quoted without one:")
    say("BRACKET log10_lean_reaches_0.01_powerlaw %.4f %.4f %.4f"
        % (pw[0], min(pw), max(pw)))
    say("BRACKET log10_lean_reaches_0.01_anylaw %.4f %.4f %.4f"
        % (pw[0], min(allf), max(allf)))
    say()
    say("  And the drift of the constant those brackets extrapolate,")
    say("  which gate check G33 reads. Within a law the constant is the")
    say("  decay exponent b, and the bracket is built from its own")
    say("  leave-one-out spread rather than from an assumed wobble, so")
    say("  the two cannot disagree:")
    bs = [abs(float(np.polyfit(x[sl], y[sl], 1)[0])) for _, sl in subs]
    dr = (max(bs) - min(bs)) / (sum(bs) / len(bs))
    say("  b runs %.4f to %.4f, a relative spread of %.4f"
        % (min(bs), max(bs), dr))
    say("DRIFT lean_exponent_b %.4f" % dr)
    say("  The wider of the two brackets above is not drift at all but")
    say("  model ambiguity: the power law and the log law are not")
    say("  separable over the accessible range, and that is a choice")
    say("  between shapes, not a wobble in a constant. It is reported")
    say("  as a bracket because the reader needs one number, but it")
    say("  should not be read as an error bar on a fit.")

    say()
    say("  DIAGNOSTIC (post hoc). Why a twelve per cent uncertainty in")
    say("  the exponent costs whole decades. The forecast solves")
    say("  |0.5-f| = c N^{-b} = %g, so log N = log(c/%g)/b: the exponent"
        % (GOAL, GOAL))
    say("  sits in the DENOMINATOR of a logarithm, and a relative error")
    say("  d in b multiplies log N by 1/(1-d), i.e. scales N by a power")
    say("  of itself.")
    e_all = [float(np.polyfit(x[s], y[s], 1)[0]) for _, s in subs]
    espan = (max(e_all) - min(e_all)) / abs(b_full)
    say("  relative spread in b            : %.4f" % espan)
    say("  so log10 N scales by up to      : %.4f"
        % (1.0 / (1.0 - espan)))
    say("  which on log10 N = %.2f is      : %.2f decades"
        % (pw[0], pw[0] * (1.0 / (1.0 - espan) - 1.0)))
    say("  measured one-sided spread       : %.2f decades"
        % (pw[0] - min(pw)))
    say("  The measured spread is the smaller of the two because the")
    say("  leave-one-out refits the intercept alongside the exponent,")
    say("  and a steeper decay comes with a smaller constant, which")
    say("  pulls the crossing back. The law choice has no such brake.")

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("the lean is confirmed by a second route and its forecast is "
        "a range, not a number" if ok else "REFUTED")

    head = [
        "STATISTIC: the mass-weighted fraction f of k < N^0.56 with",
        "           H(N;k) > 0, computed from the dilation identity",
        "           H(N;k) = mu(k) A(N;k) rather than from the inner sum",
        "           over m; the power-law exponent refitted on it; and",
        "           the N at which |0.5 - f| reaches 0.01 under each",
        "           leave-one-out subset and under each of the two laws.",
        "NULL: none is run and none is needed. lab_lean_decay.py carries",
        "      the coin arm as the reference level for this statistic;",
        "      this script recomputes the same quantity on the same",
        "      field by a second route and then does arithmetic on a",
        "      published fit, which is not a detection against a",
        "      background.",
        "FIELD: N = 2e5 to 6.4e6 by doubling, theta' = 0.56, k squarefree",
        "       and coprime to N with 2 <= k < N^0.56; Lambda and mu from",
        "       an integer sieve to 6.4e6; the published f, exponent and",
        "       forecasts are read from results/lab_lean_decay.txt.",
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
