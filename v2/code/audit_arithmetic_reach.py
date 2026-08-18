# -*- coding: utf-8 -*-
r"""
Does the arithmetic spread close with N, or is it standing?

WHAT IS AT STAKE

Every level result in these papers carries a caveat: the sweep has one
odd radical, and Remark {#rem:residuearithmetic} is why. It measures
the level exponent at seven N of comparable size and different odd
radicals and finds a spread of 0.0928, against a floor of 0.0134 --
6.9 times, so a real dependence -- with the two primorial-like N below
1/2 at 0.4808 and 0.4747.

That is one scale. Remarks {#rem:primorialrung11} and
{#rem:slopereach} have since measured how fast two of these radicals
rise: the primorial ladder at +0.007013 per unit log N and the 2^a5^b
family at +0.005112. If the hard radicals really rise faster, the
spread must be closing, and the caveat that qualifies every level
result in these papers would be a finite-N artefact rather than a
property of the method. If the spread is standing, it is not.

Multiplying N by 2 leaves its odd radical alone, so the same seven
types can be followed up the scale exactly. Four scales -- N, 2N, 4N,
8N -- give each radical its own rise and the spread its own trend.

The implementation is the bitmask one of audit_level_slope_reach.py,
independent of audit_residue_arithmetic.py's, so Q1 is a cross-check
of the published seven and not a rerun.

BACKS: Remark {#rem:arithmeticreach} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  The control: at the seven base N this implementation reproduces
      the published exponents to within 0.001.
  Q2  Every radical rises: all seven least-squares slopes of the
      exponent against log N, over the four scales, are positive.
  Q3  The spread across radicals closes: it is smaller at 8N than at
      N.
  Q4  But it has not closed: the spread at 8N still exceeds the floor
      that seven draws at the family's own scatter would give.
  Q5  And the identification is stable: the two primorial-like N are
      the bottom two of the seven at every one of the four scales.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED at 0.001 at any of the seven -- not the same statistic,
      and nothing below may be compared with
      {#rem:residuearithmetic}.
  Q2  REFUTED if any slope is negative. A radical whose level falls
      with N would be a harder obstruction than anything measured so
      far, and the "finite-N" reading of {#rem:primorialladder} would
      not extend to it.
  Q3  REFUTED if the spread is not smaller at 8N. The arithmetic
      dependence would then be standing, not finite-N, and the one-
      radical caveat could not be worked off by computing further.
  Q4  REFUTED if the spread at 8N is inside the floor. The dependence
      would already have closed at accessible N, and
      {#rem:residuearithmetic}'s reading would hold only at the scale
      it was measured.
  Q5  REFUTED if the bottom two change at any scale, which would say
      "primorial-like" is not what picks out the hard arithmetic.

  All five gate.

  NO NULL IS RUN and none applies. A measured sum is crossed against a
  computed threshold at twenty-eight N and the crossings compared;
  there is no background to detect against. The sign controls for this
  field were run in lab_residue_cancellation.py and
  lab_split_budget.py. The floor Q4 is judged against is the span of
  seven draws at the family's own scatter, the recipe
  audit_residue_arithmetic.py used, with the scatter read from
  results/audit_slope_significance.txt.
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
OUT = os.path.join(RES, "audit_arithmetic_reach.txt")

SCALES = [1, 2, 4, 8]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000
SIMS = 20_000
SEED = 20260808


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
    """the seven test N and their published exponents"""
    src = io.open(os.path.join(RES, "audit_residue_arithmetic.txt"),
                  encoding="utf-8").read()
    i = src.index("N            odd part               threshold  "
                  "K*_R    exponent  clears .5")
    ns, ex, odd = [], {}, {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            break
        N = int(f[0])
        ns.append(N)
        odd[N] = f[1]
        ex[N] = float(f[4])
    return ns, odd, ex


def count_scales(fname, header):
    """how many octaves of N a published table actually covers"""
    src = io.open(os.path.join(RES, fname), encoding="utf-8").read()
    i = src.index(header)
    oct_ = set()
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if not f or not f[0].isdigit():
            break
        oct_.add(int(math.log(int(f[0]), 2)))
    return len(oct_), sorted(oct_)


def read_family_scatter():
    """the family's own scatter about its trend -- read, not copied"""
    src = io.open(os.path.join(RES, "audit_slope_significance.txt"),
                  encoding="utf-8").read()
    return float(re.search(r"^SCATTER slope_audit_residue_level "
                           r"([\d.]+)\s*$", src, re.M).group(1))


def measure(N, lam, mu, sqf, vmask, qs, artin, twin):
    """the level exponent log K*_R / log N at one N"""
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

    base, odd, pub = read_base()
    sigma = read_family_scatter()
    say("read %d test N and their exponents from "
        "results/audit_residue_arithmetic.txt," % len(base))
    say("  and the family scatter %.4f from "
        "results/audit_slope_significance.txt" % sigma)

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

    E = {}
    for N0 in base:
        for s in SCALES:
            N = N0 * s
            e, bpn = measure(N, lam, mu, sqf, vmask, qs, artin, twin)
            E[(N0, s)] = e
            say("  N0 %-9d x%-2d = %-9d odd %-22s exponent %s"
                % (N0, s, N, odd[N0],
                   "none" if e is None else "%.4f" % e))
            say("BUDGET kstar_R_reach_N%d %.6f" % (N, bpn))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N0) if q > 2))
                  for N0 in base)))

    # ------------------------------------------------------------- Q1
    say()
    say("Q1  the control: the seven base N")
    say("  N            odd part               here     published  diff")
    q1 = True
    for N0 in base:
        e = E[(N0, 1)]
        d = abs(e - pub[N0])
        if not (d < 0.001):
            q1 = False
        say("  %-12d %-22s %-8.4f %-10.4f %.5f"
            % (N0, odd[N0], e, pub[N0], d))
    say("  Q1 %s   (cap 0.001)" % ("hold" if q1 else "REFUTED"))

    # ------------------------------------------------------------- Q2
    say()
    say("Q2  each radical followed up four scales")
    say("  four points each, so the standard errors are wide and the")
    say("  t column is printed to keep them from being over-read")
    say("  odd part               " + "".join("x%-8d" % s for s in SCALES)
        + " slope        t")
    q2 = True
    slopes = {}
    for N0 in base:
        ys = np.array([E[(N0, s)] for s in SCALES])
        xs = np.array([math.log(N0 * s) for s in SCALES])
        a, b = np.polyfit(xs, ys, 1)
        r = ys - (a * xs + b)
        se = math.sqrt(float((r ** 2).sum() / (xs.size - 2))
                       / float(((xs - xs.mean()) ** 2).sum()))
        slopes[N0] = float(a)
        if a <= 0.0:
            q2 = False
        say("  %-22s %s %+.6f    %.2f"
            % (odd[N0], "".join("%-9.4f" % y for y in ys), a,
               abs(a) / se))
    say("  Q2 every slope positive (min %+.6f)   %s"
        % (min(slopes.values()), "hold" if q2 else "REFUTED"))
    lo = min(base, key=lambda n: slopes[n])
    hi = max(base, key=lambda n: slopes[n])
    say("  slowest: odd %s at %+.6f; fastest: odd %s at %+.6f"
        % (odd[lo], slopes[lo], odd[hi], slopes[hi]))
    say("  the number of odd prime factors does not order these: the")
    say("  slowest has %d and the fastest %d."
        % (len([q for q in factor_set(lo) if q > 2]),
           len([q for q in factor_set(hi) if q > 2])))

    # ---------------------------------------------------------- Q3/Q4
    say()
    say("Q3/Q4  the spread across the seven radicals, by scale")
    rng = np.random.default_rng(SEED)
    draws = rng.normal(0.0, sigma, size=(SIMS, len(base)))
    floor = float((draws.max(axis=1) - draws.min(axis=1)).mean())
    say("  floor: expected span of %d draws at the family scatter "
        "%.4f" % (len(base), sigma))
    say("  over %d simulations = %.4f" % (SIMS, floor))
    say("FLOOR level_across_radicals_reach %.4f" % floor)
    say("  and how many octaves of N each span-against-floor in this")
    say("  repository actually covers, counted from its own N column:")
    for stem, hdr in (
            ("audit_residue_arithmetic",
             "N            odd part               threshold  "
             "K*_R    exponent  clears .5"),
            ("audit_residue_kexponent",
             "N            odd part               #odd  k-exp     "
             "level")):
        c, octs = count_scales(stem + ".txt", hdr)
        say("    %-26s %d octave(s): %s"
            % (stem, c, ", ".join("2^%d" % o for o in octs)))
        say("SCALES %s %d" % (stem, c))
        if c < 2:
            say("ONE SCALE %s" % stem)
    mine = set(int(math.log(N0 * sc, 2))
               for N0 in base for sc in SCALES)
    say("    %-26s %d octave(s)" % ("audit_arithmetic_reach", len(mine)))
    say("SCALES audit_arithmetic_reach %d" % len(mine))
    say("  scale   N about       spread    over floor")
    spreads = {}
    for s in SCALES:
        ys = [E[(N0, s)] for N0 in base]
        sp = max(ys) - min(ys)
        spreads[s] = sp
        say("  x%-6d %-12d %-9.4f %.2f"
            % (s, int(np.mean([N0 * s for N0 in base])), sp, sp / floor))
    q3 = spreads[SCALES[-1]] < spreads[SCALES[0]]
    q4 = spreads[SCALES[-1]] > floor
    say("MARGIN audit_arithmetic_reach %.4f %.4f"
        % (spreads[SCALES[-1]] - floor, floor))
    say("  Q3 the spread is smaller at x%d than at x%d "
        "(%.4f vs %.4f)   %s"
        % (SCALES[-1], SCALES[0], spreads[SCALES[-1]],
           spreads[SCALES[0]], "hold" if q3 else "REFUTED"))
    say("  Q4 and still above the floor (%.4f vs %.4f)   %s"
        % (spreads[SCALES[-1]], floor, "hold" if q4 else "REFUTED"))

    # ------------------------------------------------------------- Q5
    say()
    say("Q5  is it the same two at the bottom at every scale?")
    prim = sorted(base, key=lambda N0: len(
        [q for q in factor_set(N0) if q > 2]), reverse=True)[:2]
    say("  the two with the most odd prime factors: %s"
        % ", ".join("%d (odd %s)" % (N0, odd[N0]) for N0 in prim))
    q5 = True
    for s in SCALES:
        order = sorted(base, key=lambda N0: E[(N0, s)])
        bottom = set(order[:2])
        if bottom != set(prim):
            q5 = False
        say("  x%-3d bottom two: %s   %s"
            % (s, ", ".join(odd[N0] for N0 in order[:2]),
               "same" if bottom == set(prim) else "CHANGED"))
    say("  Q5 %s" % ("hold" if q5 else "REFUTED"))

    say()
    say("  what the spread is doing, for the record. Per unit log N")
    xs = np.log(np.array([np.mean([N0 * s for N0 in base])
                          for s in SCALES]))
    ys = np.array([spreads[s] for s in SCALES])
    a, b = np.polyfit(xs, ys, 1)
    r = ys - (a * xs + b)
    se = math.sqrt(float((r ** 2).sum() / (xs.size - 2))
                   / float(((xs - xs.mean()) ** 2).sum()))
    say("  the spread's least-squares slope against log N = %+.6f"
        % a)
    say("  r.m.s. residual %.4f, standard error %.6f, t = %.2f"
        % (float(np.sqrt((r ** 2).mean())), se, abs(a) / se))
    say("SCATTER slope_audit_arithmetic_reach %.4f"
        % float(np.sqrt((r ** 2).mean())))
    say("TSTAT slope_audit_arithmetic_reach %.2f" % (abs(a) / se))
    say("SPREAD slope_audit_arithmetic_reach %.4f"
        % float(xs.max() - xs.min()))
    if abs(a) / se < 2.0:
        say("UNRESOLVED SIGN slope_audit_arithmetic_reach")
    say("  Four scales is not a sweep and no N at which the spread")
    say("  would close is quoted; {#rem:forecastbracket} is the")
    say("  standing reason. What is measured is the sign at the two")
    say("  ends and whether the floor is still cleared.")

    say()
    say("=" * 70)
    ok = q1 and q2 and q3 and q4 and q5
    say("the arithmetic dependence is closing and has not closed"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R/log N, at",
        "           each of the seven arithmetic types of",
        "           {#rem:residuearithmetic} and at 2, 4 and 8 times",
        "           each -- multiplying by 2 leaves the odd radical",
        "           alone, so a type can be followed up the scale;",
        "           each type's least-squares slope against log N; the",
        "           spread across the seven at each scale; and the",
        "           span seven draws at the family's own scatter give.",
        "NULL: none is run and none applies. A measured sum is crossed",
        "      against a computed threshold at twenty-eight N and the",
        "      crossings compared; there is no background to detect",
        "      against. The sign controls for this field were run in",
        "      lab_residue_cancellation.py and lab_split_budget.py.",
        "FIELD: the seven test N read from",
        "       results/audit_residue_arithmetic.txt and their",
        "       doublings; k squarefree and coprime to N with",
        "       2 <= k < " + str(KCAP) + "; m odd, squarefree and",
        "       coprime to k, m < N/k; the sieve weight over the odd",
        "       primes below " + str(QSIEVE) + "; beta refitted as",
        "       sum(H P)/sum(P^2) on the same k-range; S(N) and A(N)",
        "       from Euler products at the fixed bound " + str(CLIM)
        + ";",
        "       the floor from " + str(SIMS) + " normal draws at seed "
        + str(SEED) + ",",
        "       the scatter read from",
        "       results/audit_slope_significance.txt.",
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
