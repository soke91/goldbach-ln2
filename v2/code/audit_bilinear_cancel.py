# -*- coding: utf-8 -*-
r"""
How much of sum a's smallness is cancellation over the short variable

WHAT IS AT STAKE

rem:bilinear put sum a's covered part in Type II form,

    S = - sum_{d <= D} mu(d) I(d),
    I(d) = sum_{K <= m < N/d} mu^2(m) Lambda(N-dm) log m,

with D = N^(1-theta) = N^0.44 and K = N^theta = N^0.56, the ranges
disjoint.  That names the object; it supplies no bound.  What the
classical tools give can be read off the dimensions before anything is
computed, and the arithmetic is worth writing down because it is not
encouraging:

  trivial                             exponent 1
  Bombieri-Vinogradov                 covers d <= N^(1-theta)
                                      unconditionally exactly when
                                      1-theta < 1/2, i.e. **when
                                      theta > 1/2** -- the regime the
                                      program needs -- but its error
                                      is N (log N)^-A, exponent 1
  GRH applied to each modulus d       sqrt(N) log^2 N per d, times
                                      N^(1-theta) moduli:
                                      exponent 1.5 - theta = 0.94
  measured (rem:denominator)          alpha = +0.717916
  the demand (rem:denominator)        e(l1) - theta/2 = +0.587483

**The demand sits below what GRH gives when applied to each modulus
separately.**  So meeting item 5 cannot come from per-modulus input of
any strength; it requires cancellation in the sum over d.  The measured
0.717916 is already well below 0.94, so that cancellation is happening.
This run measures how much.

Note what rem:thetalaw settled and this does not reopen: alpha =
1 - theta'/2 is a coincidence at theta' = 0.56 and was withdrawn as a
law, alpha rising with theta' at +0.978885 where that model wants
-1/2.  Nothing here fits alpha against theta'.

BACKS: Remark {#rem:bilinearcancel} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  THE GATE.  - sum_d mu(d) I(d) reproduces rem:bilinear's covered
      part, read from its POINT covered marker, to a relative 1e-12
      at every N here.

      *Disclosed.*  Y1 as first written read the covered part from
      the *table* of results/audit_bilinear.txt, which prints six
      decimals, and judged it at a relative 1e-12.  It came out
      3.53e-11, 8.31e-12 and 2.09e-12 -- refuted by the print bound
      and not by any disagreement, since 5e-7 on 13650 is 3.7e-11.
      That is TOL BELOW PRINT, which G75 exists for, and it is the
      **second** instance of exactly this defect after
      audit_which_floor's M1.  The script exits at Y1, so no verdict
      on Y2 to Y4 existed when it was found.  audit_bilinear.py now
      emits POINT covered_<N> at full double precision and the
      gate reads that -- ten digits was still a print bound of
      5e-7 absolute, the same wall one decimal further out;
      the tolerance is unchanged.
  Y2  **The d-sum cancels.**  The ratio sum_d |I(d)| / |sum_d mu(d)
      I(d)| exceeds 3 at the largest N.
  Y3  **And by about a square root.**  That ratio is within a factor
      of 2 of sqrt(#d), the number of contributing d.
  Y4  **And no single modulus is small enough.**  The largest |I(d)|
      exceeds l2 at every N, so the target is below what one term of
      the d-sum contributes and no per-modulus bound can reach it.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED outside 1e-12; nothing below is reported.
  Y2  REFUTED at 3 or below.  Then the d-sum is not where the
      smallness comes from, and since the m-sums are what per-modulus
      bounds control, the smallness would have to be inside them --
      which would make item 5 reachable by level-of-distribution
      input after all, and that would be the most consequential
      outcome this run could produce.
  Y3  REFUTED outside a factor of 2.  A ratio far **above** sqrt(#d)
      would mean the d-sum cancels better than random signs, which no
      bound would explain and which would need its own measurement;
      far below would mean the cancellation is weaker than
      square-root and the 0.94 is nearer the truth than 0.717916
      suggests.
  Y4  REFUTED if any N has max|I(d)| at or below l2.  Then a single
      modulus can be small enough and the per-modulus route is not
      excluded by size alone.

  **WHAT AN EXPONENT COMPARISON IS NOT.**  The four numbers above are
  arithmetic on the ranges, not theorems proved here; BV and GRH are
  cited for what they give, and the reading is that the demand lies
  below the per-modulus route, not that the demand is impossible.
  **No bound is established by this run and none may be claimed.**

  WHAT THIS CANNOT DO.  Three N, all small, so #d is 27 to 72 and a
  ratio compared against sqrt(#d) is compared against 5 to 8.  The
  cancellation measured here is the cancellation at these N; the
  exponents quoted are asymptotic statements from elsewhere and the
  two are not the same kind of thing.
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
OUT = os.path.join(ROOT, "results", "audit_bilinear_cancel.txt")
SRCB = os.path.join(ROOT, "results", "audit_bilinear.txt")
SRCD = os.path.join(ROOT, "results", "audit_denominator.txt")

THETA = 0.56
NS = [20_000, 50_000, 200_000]
RELID = 1e-12
CANCELCAP = 3.0
SQRTFACTOR = 2.0
DEC = 6


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
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
    rem = np.arange(n + 1, dtype=np.int64)
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


def pieces(N, lam, mu, sqf):
    """I(d) for every contributing d, and l2 over the k-range"""
    PN = factor_set(N)
    K = int(N ** THETA)
    out = []
    for d in range(1, (N - 1) // K + 1):
        md = int(mu[d])
        if md == 0 or any(d % q == 0 for q in PN):
            continue
        ms = np.arange(K, (N - 1) // d + 1, dtype=np.int64)
        if ms.size == 0:
            continue
        keep = sqf[ms]
        for q in factor_set(d) | PN:
            keep &= (ms % int(q)) != 0
        if d == 1:
            keep &= lam[ms] == 0.0
        ms = ms[keep]
        if ms.size == 0:
            continue
        v = lam[N - d * ms]
        out.append((d, md,
                    float((v * np.log(ms.astype(np.float64))).sum())))
        del ms, keep, v
    l2sq = 0.0
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        l2sq += (math.log(k)
                 * float((lam[N - ms * k]
                          * mu[ms].astype(np.float64)).sum())) ** 2
        del ms
    return out, math.sqrt(l2sq), K


def read_pub():
    src = io.open(SRCB, encoding="utf-8").read()
    cov = {}
    for N in NS:
        m = re.search(r"^POINT covered_%d ([-+]?[\d.eE+-]+)\s*$" % N,
                      src, re.M)
        if not m:
            raise SystemExit("no covered marker for N = %d" % N)
        cov[N] = float(m.group(1))
    return cov


HEAD = [
    "STATISTIC: the pieces I(d) of sum a's Type II form, the ratio of",
    "           their absolute sum to their signed sum, that ratio",
    "           against the square root of their number, and the",
    "           largest piece against l2.",
    "FIELD: N = %s; d over the squarefree d <= N/K coprime to N and m"
    % NS,
    "       over the squarefree m in [K, N/d) coprime to dN, with",
    "       K = floor(N^%.2f), the ranges of rem:bilinear. The covered"
    % THETA,
    "       part is READ from results/audit_bilinear.txt.",
    "DERIVED: BV covers d <= N^(1-theta) unconditionally exactly when",
    "         theta > 1/2; GRH applied per modulus gives exponent",
    "         1.5 - theta; neither is a bound established here.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    cov = read_pub()
    for N in NS:
        say("READ audit_bilinear.txt %d %.6f" % (N, cov[N]))
    say("  the covered part rem:bilinear published, "
        "from its ten-digit marker")
    say("  NOTE, disclosed: Y1 first read the six-decimal "
        "table and judged at")
    say("  a relative 1e-12 -- TOL BELOW PRINT, the second "
        "instance of the")
    say("  defect audit_which_floor's M1 had. The tolerance "
        "is unchanged; only")
    say("  the source of the digits is.")
    say("PRINTBOUND audit_bilinear_cancel %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, cancellation cap %.1f, square-root factor %.1f"
        % (THETA, CANCELCAP, SQRTFACTOR))
    say("  the exponents this run reads the dimensions for:")
    say("    trivial                  %.2f" % 1.0)
    say("    GRH per modulus          %.2f" % (1.5 - THETA))
    say("    BV covers the d-range    yes, since 1 - theta = %.2f "
        "< 0.5" % (1.0 - THETA))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    rows = []
    for N in NS:
        ps, l2, K = pieces(N, lam, mu, sqf)
        signed = -sum(md * v for _, md, v in ps)
        absum = sum(abs(v) for _, _, v in ps)
        biggest = max(abs(v) for _, _, v in ps)
        rows.append((N, ps, l2, K, signed, absum, biggest))
        say("  N = %-8d #d = %-4d signed %+.6f  sum|I| %.6f  "
            "max|I| %.6f  l2 %.6f"
            % (N, len(ps), signed, absum, biggest, l2))
    say("SCALES %d" % len(rows))

    # -------------------------------------------------------------- Y1
    say()
    say("Y1  the gate")
    y1 = True
    for N, ps, l2, K, signed, absum, big in rows:
        rel = abs(signed - cov[N]) / max(abs(cov[N]), 1.0)
        y1 &= rel <= RELID
        say("  N = %-8d here %+.6f against its %+.6f  relative %.2e"
            % (N, signed, cov[N], rel))
    say("  Y1 %s   (cap: %.0e relative)"
        % ("hold" if y1 else "REFUTED", RELID))
    if not y1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # ---------------------------------------------------------- Y2, Y3
    say()
    say("Y2, Y3  how much does the d-sum cancel?")
    say("    N          #d    sum|I|/|S|   sqrt(#d)   ratio")
    y2 = y3 = True
    for N, ps, l2, K, signed, absum, big in rows:
        r = absum / abs(signed)
        sq = math.sqrt(len(ps))
        f = r / sq
        if N == max(NS):
            y2 = r > CANCELCAP
        y3 &= (1.0 / SQRTFACTOR) <= f <= SQRTFACTOR
        say("  %-10d %-5d %10.4f  %9.4f  %.4f" % (N, len(ps), r, sq, f))
        say("POINT cancelratio_%d %.6f" % (N, r))
        say("POINT cancelsqrt_%d %.6f" % (N, f))
    say("  Y2 %s   (cap: above %.1f at the largest N)"
        % ("hold" if y2 else "REFUTED", CANCELCAP))
    say("  Y3 %s   (cap: within a factor %.1f of sqrt(#d))"
        % ("hold" if y3 else "REFUTED", SQRTFACTOR))
    refs = []
    for N, ps, l2, K, signed, absum, big in rows:
        vv = np.array([v for _, _, v in ps])
        refs.append(float(np.abs(vv).sum()
                          / math.sqrt(float((vv ** 2).sum()))
                          / math.sqrt(len(ps))))
    say("REFERENCE audit_bilinear_cancel %d %.4f %.4f"
        % (len(refs), min(refs), max(refs)))
    say("  sqrt(#d) is the count-based reference this run uses; the "
        "magnitude-based")
    say("  one for the same pieces is their own l1/l2, and the ratio "
        "of the two")
    say("  is printed above -- bounded by one, as Cauchy-Schwarz "
        "forces")

    # -------------------------------------------------------------- Y4
    say()
    say("Y4  is any single modulus small enough?")
    y4 = True
    for N, ps, l2, K, signed, absum, big in rows:
        ok = big > l2
        y4 &= ok
        d0 = max(ps, key=lambda t: abs(t[2]))[0]
        say("  N = %-8d max|I| %.4f at d = %-4d against l2 %.4f  %s"
            % (N, big, d0, l2, "above" if ok else "AT OR BELOW"))
        say("POINT maxpiece_%d %.6f" % (N, big / l2))
    say("  Y4 %s   (cap: max|I| above l2 at every N)"
        % ("hold" if y4 else "REFUTED"))

    say()
    say("=" * 70)
    say("Y1 %s  Y2 %s  Y3 %s  Y4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (y1, y2, y3, y4)))
    say()
    if y2 and y4:
        say("the smallness of sum a is cancellation over d, and no "
            "single modulus")
        say("is small enough for a per-modulus bound to reach the "
            "target. item 5")
        say("therefore needs cancellation between moduli, which is "
            "not what")
        say("level-of-distribution input supplies -- BV covers this "
            "d-range only")
        say("because theta > 1/2, and covering it is not the same as "
            "bounding it")
        say("below the target.")
    elif not y2:
        say("the d-sum is not where the smallness comes from, so it "
            "sits inside")
        say("the m-sums, which is what per-modulus bounds control. "
            "that would put")
        say("item 5 within reach of level-of-distribution input and "
            "is the most")
        say("consequential thing this run could have found.")
    else:
        say("a single modulus is already at or below l2, so size "
            "alone does not")
        say("exclude the per-modulus route and the reading above is "
            "not available.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
