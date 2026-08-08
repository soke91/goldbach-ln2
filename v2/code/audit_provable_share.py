# -*- coding: utf-8 -*-
r"""
Is the elementary bound's overspend really getting worse?

WHAT IS AT STAKE

OPEN item 2 rests on Remark {#rem:provablehalf}'s price. Its rule W3
reports that at accessible N the classical bound alone would spend
13.98, 15.38, 16.82, 18.29, 19.83 times the budget -- "not merely
useless but getting worse across the sweep, since exp(-c sqrt(log x))
falls more slowly than (log K)^2 rises". That is a five-point trend
read off a table, and by the standard {#rem:slopes} set it has never
been put against its own standard error. It is also read with the
implied constant A set to one, while {#rem:provablearithmetic} showed
A is a per-N maximum that falls with N.

Three things follow that can be measured. Whether the rise is
resolved; whether it survives letting A be each N's own constant
rather than one; and whether the five points can tell a power of N
from a power of log N at all -- which is what the quoted mechanism
claims.

BACKS: Remark {#rem:provableshare} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  D1  The control: the five shares reproduce {#rem:provablehalf}'s
      13.98, 15.38, 16.82, 18.29, 19.83 to within 0.01.
  D2  The rise is real: the least-squares slope of log(share)
      against log N is positive and reaches two standard errors.
  D3  And it is not the constant: with A taken as each N's own
      maximum ratio instead of one, the slope is still positive at
      two standard errors.
  D4  But five points cannot tell the mechanism: fitting log(share)
      against log N and against log log N, the two r.m.s. residuals
      differ by less than the r.m.s.'s own standard error.

REFUTATION RULE (fixed before the run)

  D1  REFUTED at 0.01 -- not the same statistic, and nothing below
      may be compared with {#rem:provablehalf}.
  D2  REFUTED below two standard errors, which would mean "getting
      worse" is a reading of five points with no error bar and has to
      be withdrawn as rule U4 was.
  D3  REFUTED if the slope loses two standard errors, or turns
      negative, once A is per N. That is the one that matters: the
      overspend would then be growing only because the constant was
      frozen at its worst value, and the price of the elementary half
      would be improving with N rather than worsening.
  D4  REFUTED if one shape beats the other by more than that standard
      error, which would make the quoted mechanism testable here and
      the sweep long enough to decide it.

  All four gate.

  NO NULL IS RUN and none applies. A deterministic bound is evaluated
  on a measured summation set and divided by a computed budget; there
  is no background to detect against. The coin arm for this field is
  lab_elementary_provable.py's sixteen coins on the identical sifted
  set, which established that the looseness is the shape and not mu.
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
OUT = os.path.join(RES, "audit_provable_share.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
CZERO = 0.2098
THETA = 0.56
CLIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def moebius(n):
    """Moebius, the cofactor kept in int32"""
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
    return mu


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
    """the five shares of rule W3"""
    src = io.open(os.path.join(RES, "lab_elementary_provable.txt"),
                  encoding="utf-8").read()
    i = src.index("N            K = N^0.56   share of the budget")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        out[int(f[0])] = float(f[2])
    return out


def measure(N, mu, oddsqf, vmask, qs, artin, twin):
    """the bound's share of the budget, and A at this N"""
    PN = factor_set(N)
    A_, S_ = artin, twin
    for q in sorted(PN):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))

    ks, Pv, Lv, inner = [], [], [], []
    for k in range(2, KCAP):
        if not oddsqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        ms = np.flatnonzero(oddsqf[1:M + 1]) + 1
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        fk = factor_set(k)
        for q in fk:
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        ks.append(k)
        Pv.append(ck * abs(int(mu[ms[keep]].sum(dtype=np.int64))))
        Lv.append(math.prod(1.0 / (1.0 - 1.0 / p) for p in fk))
        inner.append(N // k)
    ks = np.array(ks, dtype=np.int64)
    Pv = np.array(Pv)
    Lv = np.array(Lv)
    inn = np.array(inner, dtype=np.float64)
    b = inn * np.exp(-CZERO * np.sqrt(np.log(inn))) * Lv
    amax = float((Pv / b).max())
    K = N ** THETA
    sel = ks < K
    spend = float((np.log(ks[sel].astype(np.float64)) * b[sel]).sum())
    thr = S_ * (1.0 - A_) * N
    return spend / thr, amax, K > KCAP


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se, abs(a) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d published shares from "
        "results/lab_elementary_provable.txt" % len(pub))

    NMAX = max(NS)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    mu = moebius(NMAX)
    oddsqf = (mu != 0)
    oddsqf[::2] = False
    vmask = residue_mask(NMAX, qs)

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    rows = []
    for N in NS:
        sh, am, trunc = measure(N, mu, oddsqf, vmask, qs, artin, twin)
        rows.append((N, sh, am, trunc))
        say("  N = %-10d share %-9.2f A here %-9.4f  %s"
            % (N, sh, am, "truncated at the cap" if trunc else ""))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    # ------------------------------------------------------------- D1
    say()
    say("D1  the control: the published shares")
    say("  N            here       published  diff")
    d1 = True
    for N, sh, am, trunc in rows:
        d = abs(sh - pub[N])
        if not (d < 0.01):
            d1 = False
        say("  %-12d %-10.2f %-10.2f %.4f" % (N, sh, pub[N], d))
    say("  D1 %s   (cap 0.01)" % ("hold" if d1 else "REFUTED"))

    # ---------------------------------------------------------- D2/D3
    say()
    say("D2/D3  is the rise resolved, and is it the constant?")
    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))
    y1 = np.log(np.array([r[1] for r in rows]))
    y2 = np.log(np.array([r[1] * r[2] for r in rows]))
    a1, r1, s1, t1 = fit(x, y1)
    a2, r2, s2, t2 = fit(x, y2)
    d2 = (a1 > 0.0) and (t1 >= 2.0)
    d3 = (a2 > 0.0) and (t2 >= 2.0)
    say("  arm            slope        s.e.       t")
    say("  A = 1          %+-12.6f %-10.6f %.2f" % (a1, s1, t1))
    say("  A per N        %+-12.6f %-10.6f %.2f" % (a2, s2, t2))
    say("SCATTER slope_provable_share %.4f" % r1)
    say("TSTAT slope_provable_share %.2f" % t1)
    say("SPREAD slope_provable_share %.4f" % float(x.max() - x.min()))
    if t1 < 2.0:
        say("UNRESOLVED SIGN slope_provable_share")
    say("SCATTER slope_provable_share_aperN %.4f" % r2)
    say("TSTAT slope_provable_share_aperN %.2f" % t2)
    say("SPREAD slope_provable_share_aperN %.4f"
        % float(x.max() - x.min()))
    if t2 < 2.0:
        say("UNRESOLVED SIGN slope_provable_share_aperN")
    say("  A itself runs %.4f down to %.4f over the sweep"
        % (max(r[2] for r in rows), min(r[2] for r in rows)))
    say("  D2 the rise is resolved with A = 1   %s"
        % ("hold" if d2 else "REFUTED"))
    say("  D3 and with A per N   %s" % ("hold" if d3 else "REFUTED"))
    say("FROZEN provableshare %+.6f %+.6f" % (a1, a2))
    if (a1 > 0) != (a2 > 0):
        say("TREND CONVENTION provableshare")
        say("  DIAGNOSTIC on D3 (post hoc). The two arms disagree in")
        say("  SIGN, so the direction of W3's trend is a fact about")
        say("  the convention and not about the object. Neither")
        say("  convention is wrong: {#rem:provablehalf} freezes A at")
        say("  the sweep maximum because it wants an upper bound, and")
        say("  says so -- A is a function of where one looks and the")
        say("  largest value is the conservative choice. What may not")
        say("  be said is that the overspend worsens, full stop; with")
        say("  each N's own constant it improves at %.2f standard"
            % t2)
        say("  errors. The forecast is untouched either way, since it")
        say("  uses the frozen constant deliberately.")
    say("PERN provableshare_A %d %.4f %.4f"
        % (len(rows), min(r[2] for r in rows),
           max(r[2] for r in rows)))
    say("  the two arms' shares at the same N differ by exactly A,")
    say("  and both extremes fall at the same N -- the smallest")
    say("  share with A frozen and the largest with A per N are both")
    say("  at N = %d -- so a paired ratio adds nothing here and only"
        % rows[0][0])
    say("  the constant's own range is declared:")

    # ------------------------------------------------------------- D4
    say()
    say("D4  can five points tell the mechanism?")
    lx = np.log(x)
    ap, rp, sp, tp = fit(x, y1)
    al, rl, sl, tl = fit(lx, y1)
    dof = len(rows) - 2
    ser = min(rp, rl) / math.sqrt(2.0 * dof)
    d4 = abs(rp - rl) <= ser
    say("  shape                 r.m.s.     exponent")
    say("  share ~ N^e           %-10.6f %+.6f" % (rp, ap))
    say("  share ~ (log N)^e     %-10.6f %+.6f" % (rl, al))
    say("  the gap is %.6f and the r.m.s.'s own standard error is"
        % abs(rp - rl))
    say("  %.6f on %d degrees of freedom, so the gap is %.2f of it"
        % (ser, dof, abs(rp - rl) / ser))
    say("SHAPEGAP provableshare %.6f %.6f" % (abs(rp - rl), ser))
    if abs(rp - rl) <= ser:
        say("SHAPES TIED provableshare")
    say("SHAPESURVIVE provableshare %d %d %.4f"
        % (len(rows), 2 if d4 else 1, abs(rp - rl)))
    say("SHAPECURRENT provableshare %d" % len(rows))
    top = math.log10(max(NS))
    lo2, hi2 = top, 60.0
    for _ in range(200):
        mid = 0.5 * (lo2 + hi2)
        uu = mid * math.log(10.0)
        gap = abs(math.exp(ap * uu + (y1.mean() - ap * x.mean()))
                  - math.exp(al * math.log(uu)
                             + (y1.mean() - al * lx.mean())))
        if gap > min(rp, rl) * math.exp(y1.mean()):
            hi2 = mid
        else:
            lo2 = mid
    say("  the two shapes part, in the share itself, at log10 N =")
    say("  %.4f, and no forecast is made from this fit -- the only"
        % (0.5 * (lo2 + hi2)))
    say("  value read off it is at the top of the data, %.4f." % top)
    say("TRUST provableshare %.4f %.4f" % (0.5 * (lo2 + hi2), top))
    say("  D4 %s" % ("hold" if d4 else "REFUTED"))

    say()
    say("=" * 70)
    ok = d1 and d2 and d3 and d4
    say("the overspend does worsen, and not because of the constant"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the share of the budget S(N)(1-A(N))N that the",
        "           classical bound",
        "           A x exp(-c sqrt(log x)) L(k) with x = N/k would",
        "           spend over the admissible k < N^" + str(THETA)
        + ", at A = 1",
        "           and at A each N's own maximum ratio |P|/bound;",
        "           each arm's least-squares slope against log N with",
        "           its standard error; and the r.m.s. of fitting the",
        "           share against log N and against log log N.",
        "NULL: none is run and none applies. A deterministic bound is",
        "      evaluated on a measured summation set and divided by a",
        "      computed budget. The coin arm for this field is",
        "      lab_elementary_provable.py's sixteen coins on the",
        "      identical sifted set.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k odd, squarefree",
        "       and coprime to N with 2 <= k < " + str(KCAP) + ", the cap",
        "       lab_elementary_provable.py uses; m odd, squarefree and",
        "       coprime to k, m <= (N-1)/k; the sieve weight over the",
        "       odd primes below " + str(QSIEVE) + "; c = " + str(CZERO),
        "       and the Euler products at the fixed bound "
        + str(CLIM) + ".",
        "       Every N is 2^a 5^b, one odd radical, as RADICALS says.",
        "       The published shares are read from",
        "       results/lab_elementary_provable.txt.",
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
