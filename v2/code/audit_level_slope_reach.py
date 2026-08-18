# -*- coding: utf-8 -*-
r"""
Two more octaves: does the family's level slope resolve, or not?

WHAT IS AT STAKE

Remark {#rem:slopes} withdrew the reading of rule U4. The level
exponent of Remark {#rem:residuelevel} rises across the five N of
audit_residue_level.py at a least-squares slope of about +0.0047, but
that is 1.60 standard errors and its two-sigma interval contains zero.
Whether the margin over 1/2 opens or closes with N is, at present,
undetermined -- and it is the one thing about the conditional
reduction that a finite sweep could in principle settle.

The standard error of a slope falls like the spread of the abscissae.
Five octaves of N give a spread in log N of log 16; seven give log 64.
If the scatter about the trend stays what it is, that alone divides
the standard error by about 2.6 -- enough to carry +0.0047 past two
standard errors, if the slope is real. If it is not real, the extra
range will shrink the slope instead and the question stays open with
the answer "not by this method".

This script is an INDEPENDENT reimplementation, not a rerun. The
sieve weight is applied by a precomputed residue bitmask over
N - mk rather than by ten passes of a modulus test, and the constant
C_k = prod q/(q-1) is factored out of the sum instead of being carried
in the weights. Prediction N1 is the control that the two agree.

BACKS: Remark {#rem:slopereach} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  N1  The control: at the five N that audit_residue_level.py already
      covers, this independent implementation reproduces the published
      exponent log K*_R / log N to within 0.002.
  N2  The clearance is not an accessible-N artefact: at both new N,
      6400000 and 12800000, the exponent is still above 1/2.
  N3  Seven points settle the sign that {#rem:slopes} left open: the
      least-squares slope over all seven reaches two standard errors.
  N4  And it is positive -- the margin over 1/2 is opening.

REFUTATION RULE (fixed before the run)

  N1  REFUTED at 0.002 at any of the five, which would mean this is
      not the same statistic and nothing below may be compared with
      the published sweep.
  N2  REFUTED if the exponent reaches 0.5 from above at either new N.
      That would be the strongest possible result here and against the
      program: the residue-only clearance of {#rem:residuelevel} would
      be a feature of the range that was accessible when it was
      measured.
  N3  REFUTED if the slope over seven points stays below two standard
      errors, which would say the range needed to settle U4's question
      is out of reach of this computation and the question stays open.
  N4  REFUTED if the resolved slope is negative -- the margin closes,
      and the clearance at accessible N flatters the reduction exactly
      as U4 feared.

  All four gate.

  NO NULL IS RUN and none applies. A measured sum is crossed against a
  computed threshold and the crossing is located; there is no
  background to detect against. The sign controls for this field are
  lab_split_budget.py's size permutation and
  lab_residue_cancellation.py's coin arm, and the noise floor the
  slope is judged against is the scatter of this sweep's own
  residuals, exactly as in audit_slope_significance.py.
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
OUT = os.path.join(RES, "audit_level_slope_reach.txt")

OLD = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
NEW = [6_400_000, 12_800_000]
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
    """von Mangoldt and Moebius by one linear-ish pass each"""
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
    # mu by sign flips, with the unsieved cofactor picked up at the end
    mu = np.ones(n + 1, dtype=np.int8)
    cof = np.arange(n + 1, dtype=np.int64)
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


def read_published():
    """the five published exponents -- read from the results file"""
    src = io.open(os.path.join(RES, "audit_residue_level.txt"),
                  encoding="utf-8").read()
    i = src.index("N            budget factor   K*_R      "
                  "log K*_R/log N  clears .56")
    ex = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        ex[int(f[0])] = float(f[3])
    return ex


def measure(N, lam, mu, sqf, vmask, qs, artin, twin):
    """the level exponent of {#rem:residuelevel} at one N"""
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
        kb = 0
        ck = 1.0
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
        # the sieve weight is C_k on the survivors and 0 elsewhere
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
        return None
    kstar = int(ks[j])
    return (kstar, math.log(kstar) / math.log(N), S_, A_, beta,
            ks.size, thr / N)


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    rms = float(np.sqrt((r ** 2).mean()))
    s2 = float((r ** 2).sum() / (n - 2))
    se = math.sqrt(s2 / float(((x - x.mean()) ** 2).sum()))
    return float(a), rms, se, abs(float(a)) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d published exponents from "
        "results/audit_residue_level.txt" % len(pub))

    NS = OLD + NEW
    NMAX = max(NS)
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

    got = []
    for N in NS:
        r = measure(N, lam, mu, sqf, vmask, qs, artin, twin)
        if r is None:
            say("  N = %-12d no crossing below k = %d" % (N, KCAP))
            continue
        kstar, e, S_, A_, beta, nk, bpn = r
        got.append((N, kstar, e, beta, nk, bpn))
        say("  N = %-12d #k = %-7d beta = %.6f  K*_R = %-8d "
            "exponent = %.4f" % (N, nk, beta, kstar, e))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))

    rads = set()
    for N, kstar, e, beta, nk, bpn in got:
        rads.add(tuple(sorted(q for q in factor_set(N) if q > 2)))
    say("RADICALS %d" % len(rads))

    # ------------------------------------------------------------- N1
    say()
    say("N1  the control: the five N the published sweep covers")
    say("  N            here       published    diff")
    n1 = True
    for N, kstar, e, beta, nk, bpn in got:
        if N not in pub:
            continue
        d = abs(e - pub[N])
        if not (d < 0.002):
            n1 = False
        say("  %-12d %-10.4f %-12.4f %.4f" % (N, e, pub[N], d))
    say("  N1 %s   (cap 0.002)" % ("hold" if n1 else "REFUTED"))

    # ------------------------------------------------------------- N2
    say()
    say("N2  the two new octaves, against the operative budget")
    say("  S(N)(1-A(N))N :")
    say("  N            K*_R      log K*_R/log N  above .5")
    n2 = True
    for N, kstar, e, beta, nk, bpn in got:
        if N not in NEW:
            continue
        if e <= 0.5:
            n2 = False
        say("  %-12d %-9d %-15.4f %s"
            % (N, kstar, e, "yes" if e > 0.5 else "NO"))
    say("  N2 %s" % ("hold" if n2 else "REFUTED"))

    # ---------------------------------------------------------- N3/N4
    say()
    say("N3/N4  the slope {#rem:slopes} could not resolve, refitted")
    say("  over seven points instead of five")
    x7 = np.log(np.array([g[0] for g in got], dtype=np.float64))
    y7 = np.array([g[2] for g in got])
    a7, rms7, se7, t7 = fit(x7, y7)
    five = [g for g in got if g[0] in pub]
    x5 = np.log(np.array([g[0] for g in five], dtype=np.float64))
    y5 = np.array([g[2] for g in five])
    a5, rms5, se5, t5 = fit(x5, y5)
    say("  points   log N spread   slope        r.m.s.   s.e.       t")
    for lab, x, a, rms, se, t in (("five", x5, a5, rms5, se5, t5),
                                  ("seven", x7, a7, rms7, se7, t7)):
        say("  %-8s %-14.4f %-12.6f %-8.4f %-10.6f %.2f"
            % (lab, float(x.max() - x.min()), a, rms, se, t))
    say("SCATTER slope_audit_level_slope_reach %.4f" % rms7)
    say("TSTAT slope_audit_level_slope_reach %.2f" % t7)
    say("SPREAD slope_audit_level_slope_reach %.4f"
        % float(x7.max() - x7.min()))
    if t7 < 2.0:
        say("UNRESOLVED SIGN slope_audit_level_slope_reach")
    n3 = t7 >= 2.0
    n4 = a7 > 0.0
    say("  N3 seven points reach two standard errors (%.2f)   %s"
        % (t7, "hold" if n3 else "REFUTED"))
    say("  N4 and the slope is positive (%+.6f)   %s"
        % (a7, "hold" if n4 else "REFUTED"))
    say("  two-sigma interval over seven: [%+.6f, %+.6f]"
        % (a7 - 2 * se7, a7 + 2 * se7))

    say()
    say("  what the extra range actually bought. The spread in log N")
    say("  grew by the factor %.4f and the standard error moved by"
        % (float((x7.max() - x7.min()) / (x5.max() - x5.min()))))
    say("  the factor %.4f; had the scatter been unchanged it would"
        % (se7 / se5))
    say("  have moved by %.4f. The difference is the scatter, which"
        % float(np.sqrt(((x5 - x5.mean()) ** 2).sum()
                        / ((x7 - x7.mean()) ** 2).sum())))
    say("  went from %.4f to %.4f." % (rms5, rms7))

    say()
    say("=" * 70)
    ok = n1 and n2 and n3 and n4
    say("the sign of the level slope is settled over seven octaves"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the level exponent log K*_R / log N of Remark",
        "           {#rem:residuelevel}, recomputed by an independent",
        "           implementation at the five N it publishes and at",
        "           two further octaves, 6400000 and 12800000; then",
        "           the least-squares slope of that exponent against",
        "           log N over five points and over seven, each with",
        "           its r.m.s. residual and standard error",
        "           sqrt(sum r^2/(n-2) / sum (x-xbar)^2).",
        "NULL: none is run and none applies. A measured sum is crossed",
        "      against a computed threshold and the crossing located;",
        "      there is no background to detect against. The sign",
        "      controls for this field are lab_split_budget.py's size",
        "      permutation and lab_residue_cancellation.py's coin arm.",
        "FIELD: k squarefree and coprime to N, from 2 to "
        + str(KCAP) + "; m odd, squarefree and coprime to k, up to",
        "       (N-1)/k; the sieve weight over the odd primes below "
        + str(QSIEVE) + ";",
        "       the Euler products to " + str(CLIM) + ". Every N is",
        "       2^a 5^b, so the sweep has one odd radical, as the",
        "       RADICALS line above declares.",
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
