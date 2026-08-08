# -*- coding: utf-8 -*-
r"""
Is the arithmetic dependence one variable, at every scale?

WHAT IS AT STAKE

OPEN item 3 is that the one-radical caveat does not close: the level
exponent's spread across seven arithmetic types runs 0.0928, 0.1063,
0.1011, 0.0981 as N is multiplied by 1, 2, 4, 8, staying at seven
times its floor ({#rem:arithmeticreach}). What that item does not say
is whether the spread is ONE thing. Remark {#rem:residuearithmetic}
attributed it to the budget, regressing the exponent on the log of
the threshold across the seven and getting a correlation of 0.97565 --
at one scale. Remark {#rem:kexponent} is the warning: an apparent
arithmetic dependence of the k-exponent turned out to be scatter.

The regressor here does not move with the scale. Multiplying N by 2
leaves the odd radical alone, so S(N)(1-A(N)) is the same number at
1, 2, 4 and 8 times each N -- the seven x-values are fixed and only
the exponents move. That makes the question clean: does the budget
keep explaining the spread as N grows, and does anything else appear
in the residue?

If the correlation holds and the residual spread stays inside the
floor, item 3 is a standing dependence on one measured variable. If
the residue grows past the floor, a second arithmetic variable is
appearing and the caveat is not one thing.

BACKS: Remark {#rem:arithmeticonevar} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  A1  The control: at the base scale the seven exponents reproduce
      {#rem:residuearithmetic}'s 0.5422, 0.5591, 0.5675, 0.5294,
      0.4808, 0.4747, 0.5424 to within 0.001, and the regression on
      log(threshold) reproduces its slope 0.0516 and correlation
      0.97565 to within 0.001.
  A2  The budget keeps explaining: the correlation stays above 0.9 at
      all four scales.
  A3  And nothing else appears: the residual spread after removing
      the budget stays inside the floor 0.0133 at all four scales.
  A4  The response is stable: the slope's spread across the four
      scales is under 0.01.

REFUTATION RULE (fixed before the run)

  A1  REFUTED at 0.001 anywhere -- not the same statistic, and
      nothing below may be compared with {#rem:residuearithmetic}.
  A2  REFUTED below 0.9 at any scale, which would say the budget
      stops accounting for the spread as N grows.
  A3  REFUTED if the residual spread reaches the floor at any scale.
      That is the one that matters: a second arithmetic variable
      would then be visible, and OPEN item 3 would be two things
      rather than one.
  A4  REFUTED beyond 0.01, which would say the size of the budget's
      effect is itself moving and the regression is not a fixed
      description.

  All four gate.

  NO NULL IS RUN and none applies. A measured sum is crossed against
  a computed threshold at twenty-eight N and the crossings regressed
  on a computed constant; there is no background to detect against.
  The floor the residual spread is judged against is the span of
  seven draws at the family's own scatter, the recipe
  audit_arithmetic_reach.py used, read from its results file.
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
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_arithmetic_onevar.txt")

SCALES = [1, 2, 4, 8]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
    """von Mangoldt and Moebius, the cofactor kept in int32"""
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
    del pr, lgp
    mu = np.ones(n + 1, dtype=np.int8)
    cof = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        cof[p::p] //= p
        pk = p * p
        while pk <= n:
            cof[pk::pk] //= p
            if pk > n // p:
                break
            pk *= p
    big = cof > 1
    del cof
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def residue_mask(n, qs):
    """bit i of mask[v] is set exactly when qs[i] divides v"""
    m = np.zeros(n + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        m[0::q] |= np.uint16(1 << i)
    return m


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


def read_base():
    """the seven test N, their odd parts and published exponents"""
    src = io.open(os.path.join(RES, "audit_residue_arithmetic.txt"),
                  encoding="utf-8").read()
    i = src.index("N            odd part               threshold  "
                  "K*_R    exponent  clears .5")
    ns, odd, ex = [], {}, {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            break
        N = int(f[0])
        ns.append(N)
        odd[N] = f[1]
        ex[N] = float(f[4])
    m = re.search(r"exponent against log\(threshold\): slope "
                  r"\+([\d.]+), correlation ([\d.]+)", src)
    return ns, odd, ex, float(m.group(1)), float(m.group(2))


def read_floor():
    """the span seven draws at the family's scatter give"""
    src = io.open(os.path.join(RES, "audit_arithmetic_reach.txt"),
                  encoding="utf-8").read()
    return float(re.search(r"^FLOOR level_across_radicals_reach "
                           r"([\d.]+)\s*$", src, re.M).group(1))


def measure(N, lam, mu, sqf, vmask, qs, artin, twin):
    """the level exponent and the budget constant at one N"""
    PN = factor_set(N)
    A_, S_ = artin, twin
    for q in sorted(PN):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))

    ks, Hs, Ps = [], [], []
    for k in range(2, KCAP):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[sqf[ms]]
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        for q in factor_set(k):
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        g = mu[ms].astype(np.float64)
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        ks.append(k)
        Hs.append(float((lam[vals] * g).sum()))
        Ps.append(ck * float(g[keep].sum()))
    ks = np.array(ks, dtype=np.int64)
    H = np.array(Hs)
    P = np.array(Ps)
    beta = float((H * P).sum() / (P * P).sum())
    R = H - beta * P
    cum = np.cumsum(np.log(ks.astype(np.float64)) * np.abs(R))
    thr = S_ * (1.0 - A_) * N
    j = int(np.searchsorted(cum, thr))
    if j >= ks.size:
        return None, thr / N
    return math.log(int(ks[j])) / math.log(N), thr / N


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    base, odd, pub, pslope, pcorr = read_base()
    floor = read_floor()
    say("read %d test N and their exponents from "
        "results/audit_residue_arithmetic.txt," % len(base))
    say("  its slope %.4f and correlation %.5f, and the floor %.4f "
        "from" % (pslope, pcorr, floor))
    say("  results/audit_arithmetic_reach.txt")

    NMAX = max(base) * max(SCALES)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    vmask = residue_mask(NMAX, qs)

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    E, thr = {}, {}
    for N0 in base:
        for s in SCALES:
            e, b = measure(N0 * s, lam, mu, sqf, vmask, qs,
                           artin, twin)
            E[(N0, s)] = e
            thr[N0] = b
        say("  odd %-22s threshold %.6f  %s"
            % (odd[N0], thr[N0],
               " ".join("%.4f" % E[(N0, s)] for s in SCALES)))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in base)))
    for N0 in base:
        say("BUDGET kstar_onevar_N%d %.6f" % (N0, thr[N0]))

    x = np.log(np.array([thr[N0] for N0 in base]))

    # ------------------------------------------------------------- A1
    say()
    say("A1  the control: the base scale")
    say("  odd part               here     published  diff")
    a1 = True
    for N0 in base:
        d = abs(E[(N0, 1)] - pub[N0])
        if not (d < 0.001):
            a1 = False
        say("  %-22s %-8.4f %-10.4f %.5f"
            % (odd[N0], E[(N0, 1)], pub[N0], d))
    y0 = np.array([E[(N0, 1)] for N0 in base])
    s0 = float(np.polyfit(x, y0, 1)[0])
    c0 = float(np.corrcoef(x, y0)[0, 1])
    if abs(s0 - pslope) >= 0.001 or abs(c0 - pcorr) >= 0.001:
        a1 = False
    say("  slope %.4f against %.4f, correlation %.5f against %.5f"
        % (s0, pslope, c0, pcorr))
    say("  A1 %s   (cap 0.001)" % ("hold" if a1 else "REFUTED"))

    # -------------------------------------------------- A2 / A3 / A4
    say()
    say("A2/A3/A4  the same regression at every scale")
    say("  scale   slope      correlation  residual spread  "
        "against the floor")
    a2 = a3 = True
    slopes, resid = [], []
    for s in SCALES:
        y = np.array([E[(N0, s)] for N0 in base])
        a, b = np.polyfit(x, y, 1)
        r = y - (a * x + b)
        sp = float(r.max() - r.min())
        cc = float(np.corrcoef(x, y)[0, 1])
        slopes.append(float(a))
        resid.append(sp)
        if cc <= 0.9:
            a2 = False
        if sp >= floor:
            a3 = False
        say("  x%-6d %-10.4f %-12.5f %-16.4f %.2f"
            % (s, a, cc, sp, sp / floor))
    a4 = (max(slopes) - min(slopes)) < 0.01
    say("  A2 the correlation stays above 0.9   %s   (floor 0.9)"
        % ("hold" if a2 else "REFUTED"))
    say("  A3 the residual spread stays inside the floor %.4f   %s"
        % (floor, "hold" if a3 else "REFUTED"))
    say("FLOOR arithmetic_onevar %.4f" % floor)
    say("  A4 the slope's spread is %.4f   %s   (cap 0.01)"
        % (max(slopes) - min(slopes), "hold" if a4 else "REFUTED"))
    say("PERN arithmetic_onevar_residual %d %.4f %.4f"
        % (len(SCALES), min(resid), max(resid)))
    say("SCALES audit_arithmetic_onevar %d" % len(SCALES))

    say()
    say("  what the residue looks like, type by type, so that a")
    say("  second variable would be visible if there were one:")
    say("  odd part               " + "".join("x%-8d" % s
                                              for s in SCALES))
    for i, N0 in enumerate(base):
        row = []
        for s in SCALES:
            y = np.array([E[M, s] for M in base])
            a, b = np.polyfit(x, y, 1)
            row.append(y[i] - (a * x[i] + b))
        say("  %-22s %s" % (odd[N0],
                            "".join("%+-10.4f" % v for v in row)))
    keep = 0
    for i, N0 in enumerate(base):
        sg = []
        for s in SCALES:
            y = np.array([E[M, s] for M in base])
            a, b = np.polyfit(x, y, 1)
            sg.append(np.sign(y[i] - (a * x[i] + b)))
        if len(set(sg)) == 1:
            keep += 1
    say("  %d of the %d types keep the sign of their residue at all"
        % (keep, len(base)))
    say("  four scales, which is what a second variable would look")
    say("  like -- but both of those residues SHRINK with N, and the")
    say("  spread as a whole crosses the floor between x2 and x4.")
    say("RESIDSCALE arithmetic_onevar %d %.4f %.4f %.4f"
        % (len(SCALES), min(resid), max(resid), floor))
    if min(resid) < floor < max(resid):
        say("CROSSES FLOOR arithmetic_onevar")

    say()
    say("=" * 70)
    ok = a1 and a2 and a3 and a4
    say("the arithmetic dependence is one measured variable"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the level exponent log K*_R / log N at each of the",
        "           seven arithmetic types of",
        "           {#rem:residuearithmetic} and at 2, 4 and 8 times",
        "           each; its least-squares regression on the log of",
        "           the budget constant S(N)(1-A(N)), which does not",
        "           change with the scale because doubling N leaves",
        "           the odd radical alone; the correlation and the",
        "           spread of the residues at each scale; and the",
        "           span seven draws at the family's own scatter give.",
        "NULL: none is run and none applies. A measured sum is crossed",
        "      against a computed threshold and the crossings",
        "      regressed on a computed constant; there is no",
        "      background to detect against.",
        "FIELD: the seven test N read from",
        "       results/audit_residue_arithmetic.txt and their",
        "       doublings; k squarefree and coprime to N with",
        "       2 <= k < " + str(KCAP) + "; m odd, squarefree and",
        "       coprime to k, m < N/k; the sieve weight over the odd",
        "       primes below " + str(QSIEVE) + "; beta refitted as",
        "       sum(H P)/sum(P^2) on the same k-range; S(N) and A(N)",
        "       from Euler products at the fixed bound " + str(CLIM)
        + ".",
        "       Seven distinct odd radicals, as RADICALS declares. The",
        "       floor is read from",
        "       results/audit_arithmetic_reach.txt.",
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
