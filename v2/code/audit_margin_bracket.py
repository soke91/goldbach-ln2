# -*- coding: utf-8 -*-
r"""
The bracket on the wall's extrapolation to 10^8, and an independent
check of the base it stands on.

WHAT IS AT STAKE

audit_margin.py extrapolates the wall's reciprocal one factor of 6.25
past the computed range: 99.325 * 6.25^0.43 = 218.42, against the
paper's "near 220", and records M3 as holding.

Two things are wrong with quoting that. The exponent 0.43 is the
PUBLISHED one, and the same audit refutes it as M2 -- refitting the
octave maxima gives 0.477527 on the grid anchored at 3e4 and 0.440619
on the grid anchored at 1.6e7, and the leave-one-out spread on the
first is 0.0173. So M3 mixes a recomputed base with a refuted
exponent. And gate check G28 asks every forecast beyond the computed
range to carry a bracket; this one has none.

The base is checked here without repeating the FFT that produced it.
audit_margin.py convolves mu * Lambda by a length-2^26 real FFT; that
transform is where a silent error would live, and rerunning it would
not find one. Instead C(N) is summed directly, term by term, at each
octave's argmax N: C(N) = sum_{n<N} Lambda(n) mu(N-n), one pass over
the sieve, sharing no arithmetic with the transform.

BACKS: Remark {#rem:marginbracket} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  Direct summation reproduces every published octave maximum at
      its published argmax N to within 1e-6 in max|C|/N. The FFT is
      not lying.
  U2  Refitting the two published octave tables reproduces the two
      published exponents to within 1e-4, so this script is reading
      them correctly.
  U3  The extrapolation bracket -- over both grid anchors and every
      leave-one-out subset of each -- spans less than a factor 1.5.
      Unlike the two forecasts already bracketed, this one reaches
      only 6.25 times past its data.
  U4  And the published "near 220" lies inside that bracket.

REFUTATION RULE (fixed before the run)

  U1  REFUTED at 1e-6 at any octave, which would mean the FFT and
      the direct sum disagree and one of them is wrong.
  U2  REFUTED at 1e-4 at either grid.
  U3  REFUTED if the bracket reaches a factor 1.5.
  U4  REFUTED if 220 falls outside the bracket. That is the outcome
      worth having: it would mean the published figure rests on the
      exponent 0.43 that M2 refutes, and not on any exponent the
      data support.

  All four gate.

  NO NULL IS RUN and none applies. A published deterministic
  extrapolation is refitted and perturbed, and one measured statistic
  is recomputed by a second route; there is no detection against a
  background. The control for max|C|/N itself is the Cauchy-Schwarz
  closed form already carried as M5 in audit_margin.py.
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
OUT = os.path.join(ROOT, "results", "audit_margin_bracket.txt")

XMAX = 16_000_000
STRETCH = 6.25            # 1.6e7 -> 1e8, the factor audit_margin uses
PUBLISHED = 220.0         # the paper's "near 220" -- a reference value


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


def read_tables():
    """the two octave tables and the two exponents -- read, not copied"""
    p = os.path.join(ROOT, "results", "audit_margin.txt")
    src = io.open(p, encoding="utf-8").read()

    def grab(anchor):
        i = src.index(anchor)
        rows = []
        for ln in src[i:].splitlines()[1:]:
            f = ln.split()
            if len(f) != 4 or not f[0].isdigit():
                if rows:
                    break
                continue
            rows.append((int(f[0]), float(f[2]), int(f[3])))
        return rows

    up = grab("octave top     count      max |C|/N     argmax N")
    j = src.index("On that grid:")
    dn = []
    for ln in src[j:].splitlines()[2:]:
        f = ln.split()
        if len(f) != 4 or not f[0].isdigit():
            if dn:
                break
            continue
        dn.append((int(f[0]), float(f[2]), int(f[3])))
    bup = float(re.search(r"max\|C\|/N ~ N\^\{-([\d.]+)\}", src).group(1))
    bdn = float(re.search(r"fitted b = ([\d.]+)", src).group(1))
    base = float(re.search(r"1/top = ([\d.]+)", src).group(1))
    bpub = float(re.search(r"published ([\d.]+)\s+REFUTED", src).group(1))
    return up, dn, bup, bdn, base, bpub


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    up, dn, bup, bdn, base, bpub = read_tables()
    say("read from results/audit_margin.txt: %d octaves anchored low, "
        "%d anchored high," % (len(up), len(dn)))
    say("  refitted exponents %.6f and %.6f, the published one %.4f,"
        % (bup, bdn, bpub))
    say("  and the base 1/max = %.3f" % base)

    say()
    say("sieving to %d ..." % XMAX)
    lam, mu = sieves(XMAX)

    # ------------------------------------------------------------- U1
    say()
    say("U1  the base, by direct summation instead of the FFT")
    say("  C(N) = sum_{n<N} Lambda(n) mu(N-n), one pass per argmax N")
    say("  octave top     argmax N     direct       from the FFT   diff")
    u1 = True
    seen = {}
    for top, val, arg in up + dn:
        if arg in seen:
            direct = seen[arg]
        else:
            c = float(np.dot(lam[1:arg], mu[arg - 1:0:-1]
                             .astype(np.float64)))
            direct = abs(c) / arg
            seen[arg] = direct
        d = abs(direct - val)
        if d >= 1e-6:
            u1 = False
        say("  %-14d %-12d %-12.8f %-14.6f %.2e"
            % (top, arg, direct, val, d))
    say("  U1 %s" % ("hold" if u1 else "REFUTED"))
    del lam, mu

    # ------------------------------------------------------------- U2
    def fit(rows, sel=slice(None)):
        t = np.array([r[0] for r in rows], dtype=float)[sel]
        v = np.array([r[1] for r in rows], dtype=float)[sel]
        return -float(np.polyfit(np.log(t), np.log(v), 1)[0])

    say()
    say("U2  the two published exponents, refitted here")
    say("  grid                     refitted    published   diff")
    u2 = True
    for name, rows, pub in (("anchored at 3e4, doubling", up, bup),
                            ("anchored at 1.6e7, halving", dn, bdn)):
        b = fit(rows)
        d = abs(b - pub)
        if d >= 1e-4:
            u2 = False
        say("  %-24s %-11.6f %-11.6f %.2e" % (name, b, pub, d))
    say("  U2 %s" % ("hold" if u2 else "REFUTED"))

    # ------------------------------------------------------------- U3
    say()
    say("U3  the bracket: every grid anchor, every leave-one-out subset")
    say("  grid                     subset                 b        "
        "%.2f^b * base" % STRETCH)
    subs = [("all octaves", slice(None)),
            ("without the smallest", slice(1, None)),
            ("without the largest", slice(0, -1))]
    vals = []
    for name, rows in (("anchored at 3e4", up),
                       ("anchored at 1.6e7", dn)):
        for sname, s in subs:
            b = fit(rows, s)
            e = base * STRETCH ** b
            vals.append(e)
            say("  %-24s %-22s %-8.6f %.2f" % (name, sname, b, e))
    lo, hi = min(vals), max(vals)
    u3 = hi / lo < 1.5
    say("  bracket [%.2f, %.2f], span %.4f   (cap 1.5)   %s"
        % (lo, hi, hi / lo, "hold" if u3 else "REFUTED"))

    # ------------------------------------------------------------- U4
    say()
    say("U4  does the published figure lie inside it?")
    u4 = lo <= PUBLISHED <= hi
    say("  published 'near %.0f' against [%.2f, %.2f]   %s"
        % (PUBLISHED, lo, hi, "hold" if u4 else "REFUTED"))
    say("  the value audit_margin.py prints, %.2f, is built on the"
        % (base * STRETCH ** bpub))
    say("  published exponent %.4f, which its own M2 refutes; on the"
        % bpub)
    say("  exponents the data support the extrapolation is %.2f to %.2f"
        % (lo, hi))

    say()
    say("  Bracket lines, in the form the gate reads. This file supplies")
    say("  them for results/audit_margin.txt, whose extrapolation is")
    say("  quoted without one:")
    say("BRACKET wall_reciprocal_at_1e8 %.4f %.4f %.4f"
        % (base * STRETCH ** fit(dn), lo, hi))
    say()
    say("  And the drift of the constant the bracket was built from,")
    say("  which gate check G33 reads. This bracket is not made by")
    say("  assuming a wobble: it IS the measured spread of the exponent")
    say("  across both grid anchors and every leave-one-out subset, so")
    say("  the wobble and the drift are the same number by")
    say("  construction and the bracket cannot understate itself.")
    ball2 = [fit(r, sel) for r in (up, dn) for _, sel in subs]
    dr = (max(ball2) - min(ball2)) / (sum(ball2) / len(ball2))
    say("  b runs %.6f to %.6f, a relative spread of %.4f"
        % (min(ball2), max(ball2), dr))
    say("DRIFT wall_exponent_b %.4f" % dr)

    say()
    say("  DIAGNOSTIC (post hoc). Why this bracket is narrow where the")
    say("  other two the repository has bracketed were not. The")
    say("  extrapolation multiplies the base by (reach)^b, so an")
    say("  absolute error d in the exponent costs exactly (reach)^d:")
    say("  the reach enters only through its logarithm, and this reach")
    say("  is a factor %.2f rather than millions." % STRETCH)
    ball = [fit(r, s) for r in (up, dn) for _, s in subs]
    espan = max(ball) - min(ball)
    say("  spread in b across grids and subsets : %.6f" % espan)
    say("  the factor that buys, %.2f^spread    : %.4f"
        % (STRETCH, STRETCH ** espan))
    say("  measured span                        : %.4f" % (hi / lo))
    say("  Reach is what sets a bracket, not the quality of the fit:")
    say("  the fit here is the WORST of the three -- its own exponent")
    say("  is refuted against the published one -- and its bracket is")
    say("  the tightest by orders of magnitude.")

    say()
    say("=" * 70)
    ok = u1 and u2 and u3 and u4
    say("the extrapolation is tight and the published figure sits in it"
        if ok else "REFUTED")

    head = [
        "STATISTIC: max |C(N)|/N at each octave's argmax N, recomputed by",
        "           direct summation rather than by FFT convolution; the",
        "           decay exponent refitted on each of the two published",
        "           octave grids and on each leave-one-out subset of",
        "           each; and the extrapolation base * 6.25^b to N = 1e8",
        "           under every one of them.",
        "NULL: none is run and none applies. A published deterministic",
        "      extrapolation is refitted and perturbed, and one measured",
        "      statistic is recomputed by a second route; there is no",
        "      background to detect against. The control for max|C|/N is",
        "      the Cauchy-Schwarz closed form carried as M5 in",
        "      audit_margin.py.",
        "FIELD: even N; Lambda and mu from an integer sieve to 1.6e7; the",
        "       octave tables, the two exponents and the base 1/max are",
        "       read from results/audit_margin.txt; the reach 6.25 and",
        "       the reference value 220 are the paper's.",
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
