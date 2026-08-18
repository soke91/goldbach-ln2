# -*- coding: utf-8 -*-
r"""
Is 10^4785 a forecast, or a forecast about a frozen constant?

WHAT IS AT STAKE

OPEN item 2 is Remark {#rem:provablehalf}'s price: the classical bound
pays the budget at N = 10^5474.8, bracketed over the analytic exponent
and, by {#rem:provablearithmetic}, over the arithmetic. Remark
{#rem:provableshare} has since shown that the implied constant A is
not a constant -- it is each N's own maximum ratio and it falls from
1.2119 to 0.3487 across the accessible sweep, which reverses the sign
of the overspend's trend.

The forecast solves A dL I(u) = S(1-A) with A frozen at the sweep
maximum. If A keeps falling, the left side is smaller at large u than
the forecast assumes and the crossing moves. How far it moves depends
on A's shape over four thousand decades -- which is exactly the kind
of extrapolation this repository has spent cycles refusing to make.

So the question is not "what is the corrected forecast" but whether
one exists: if the two shapes A's own data cannot separate give
crossings that differ by more than the bracket already published, then
the forecast is undefined once A is allowed to move, and the published
figure is a statement about a frozen constant rather than about the
bound.

BACKS: Remark {#rem:provableforecast} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  The control: dL, the frozen A and the forecast reproduce
      {#rem:provablehalf}'s 0.3994, 1.2119 and 5474.8 to within
      0.001, 0.001 and one decade.
  E2  A's decay is real: the least-squares slope of log A against
      log N is negative and reaches two standard errors.
  E3  But its shape is not fixed: fitting log A against log N and
      against log log N, the two r.m.s. residuals differ by less than
      the r.m.s.'s own standard error.
  E4  And the forecast does not survive it: extrapolating A by the
      two shapes gives crossings that differ by more than one decade,
      or one of them gives none at all.

REFUTATION RULE (fixed before the run)

  E1  REFUTED at any cap -- not the same statistic, and nothing below
      may be compared with {#rem:provablehalf}.
  E2  REFUTED below two standard errors. A would then be constant to
      the precision this sweep can see and freezing it would be
      right.
  E3  REFUTED if one shape wins by more than that standard error,
      which would mean A's decay is pinned and a corrected forecast
      could be computed rather than merely bracketed.
  E4  REFUTED if the two crossings agree within a decade. That is the
      one that matters: the forecast would then be robust to letting
      A move, and {#rem:provablehalf}'s figure would stand as a
      statement about the bound and not about a convention.

  All four gate.

  NO NULL IS RUN and none applies. A deterministic equation is solved
  with measured inputs; there is no background to detect against. The
  coin arm for this field is lab_elementary_provable.py's sixteen
  coins on the identical sifted set.
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
OUT = os.path.join(RES, "audit_provable_forecast.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
CZERO = 0.2098
THETA = 0.56
CLIM = 4_000_000
UMAX = 40_000.0


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
    """dL, the frozen A and the forecast"""
    src = io.open(os.path.join(RES, "lab_elementary_provable.txt"),
                  encoding="utf-8").read()
    dl = float(re.search(r"density-times-L factor dL = ([\d.]+)",
                         src).group(1))
    a = float(re.search(r"constant A forced by the data = ([\d.]+)",
                        src).group(1))
    m = re.search(r"^BRACKET log10_N_elementary_provable "
                  r"([\d.]+) ([\d.]+) ([\d.]+)\s*$", src, re.M)
    return dl, a, float(m.group(1)), float(m.group(2)), \
        float(m.group(3))


def inputs(N, mu, oddsqf, vmask, qs, artin, twin):
    """A at this N, the density factor dL and the budget constant"""
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
    w = np.log(ks.astype(np.float64)) / ks.astype(np.float64)
    dl = float((w * Lv).sum()
               / (math.log(float(ks.max())) ** 2 / 2.0))
    return float((Pv / b).max()), dl, S_ * (1.0 - A_)


def integral(u, c):
    v = np.linspace(0.0, THETA * u, 20001)
    return float(np.trapezoid(v * np.exp(-c * np.sqrt(u - v)), v))


def solve(afun, dl, thrc):
    """where A(u) dL I(u) falls back through the budget constant"""
    lo, hi = 10.0, UMAX
    if afun(hi) * dl * integral(hi, CZERO) > thrc:
        return None
    if afun(lo) * dl * integral(lo, CZERO) <= thrc:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if afun(mid) * dl * integral(mid, CZERO) > thrc:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi) / math.log(10.0)


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(b), float(np.sqrt((r ** 2).mean())), se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pdl, pa, ppoint, pclo, pchi = read_published()
    say("read dL %.4f, A %.4f and the bracket %.1f [%.1f, %.1f] from"
        % (pdl, pa, ppoint, pclo, pchi))
    say("  results/lab_elementary_provable.txt")

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
        a, dl, thrc = inputs(N, mu, oddsqf, vmask, qs, artin, twin)
        rows.append((N, a, dl, thrc))
        say("  N = %-10d A %-9.4f dL %-9.4f S(1-A) %.6f"
            % (N, a, dl, thrc))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    dl = rows[-1][2]
    thrc = rows[-1][3]
    afrozen = max(r[1] for r in rows)

    # ------------------------------------------------------------- E1
    say()
    say("E1  the control")
    base = solve(lambda u: afrozen, dl, thrc)
    e1 = (abs(dl - pdl) < 0.001 and abs(afrozen - pa) < 0.001
          and base is not None and abs(base - ppoint) < 1.0)
    say("  dL %.4f against %.4f; A frozen %.4f against %.4f"
        % (dl, pdl, afrozen, pa))
    say("  forecast %.1f against %.1f"
        % (-1.0 if base is None else base, ppoint))
    say("  E1 %s   (cap 0.001 on dL, cap 0.001 on A, cap one decade)"
        % ("hold" if e1 else "REFUTED"))

    # ---------------------------------------------------------- E2/E3
    say()
    say("E2/E3  what A itself does")
    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))
    y = np.log(np.array([r[1] for r in rows]))
    ap, bp, rp, sp = fit(x, y)
    lx = np.log(x)
    al, bl, rl, sl = fit(lx, y)
    e2 = (ap < 0.0) and (abs(ap) / sp >= 2.0)
    dof = len(rows) - 2
    ser = min(rp, rl) / math.sqrt(2.0 * dof)
    e3 = abs(rp - rl) <= ser
    say("  shape                 slope        r.m.s.     s.e.")
    say("  log A ~ log N         %+-12.6f %-10.6f %.6f" % (ap, rp, sp))
    say("  log A ~ log log N     %+-12.6f %-10.6f %.6f" % (al, rl, sl))
    say("  A runs %.4f down to %.4f, at %.2f standard errors"
        % (max(r[1] for r in rows), min(r[1] for r in rows),
           abs(ap) / sp))
    say("  the two r.m.s. differ by %.6f against the r.m.s.'s own"
        % abs(rp - rl))
    say("  standard error %.6f on %d degrees of freedom, i.e. %.2f"
        % (ser, dof, abs(rp - rl) / ser))
    say("SHAPEGAP provableforecast %.6f %.6f" % (abs(rp - rl), ser))
    if abs(rp - rl) <= ser:
        say("SHAPES TIED provableforecast")
    say("  E2 the decay is resolved   %s" % ("hold" if e2 else
                                             "REFUTED"))
    say("  E3 and its shape is not   %s" % ("hold" if e3 else
                                            "REFUTED"))

    # ------------------------------------------------------------- E4
    say()
    say("E4  the forecast under each extrapolation of A")
    cp = solve(lambda u: math.exp(bp + ap * u), dl, thrc)
    cl = solve(lambda u: math.exp(bl + al * math.log(u)), dl, thrc)
    say("  A frozen at the sweep maximum   log10 N = %.1f"
        % (-1.0 if base is None else base))
    say("  A ~ N^%.6f                     log10 N = %s"
        % (ap, "none" if cp is None else "%.1f" % cp))
    say("  A ~ (log N)^%.6f               log10 N = %s"
        % (al, "none" if cl is None else "%.1f" % cl))
    if cp is None or cl is None:
        e4 = True
        gap = float("inf")
    else:
        gap = abs(cp - cl)
        e4 = gap > 1.0
    say("  the two differ by %s decades"
        % ("infinitely many -- one gives none" if gap == float("inf")
           else "%.1f" % gap))
    say("  E4 %s   (floor one decade)" % ("hold" if e4 else "REFUTED"))
    if cp is not None and cl is not None:
        say("SHAPESURVIVE provableforecast %d %d %.4f"
            % (len(rows), 2 if e3 else 1, gap))
        say("SHAPECURRENT provableforecast %d" % len(rows))
        say("TRUST provableforecast %.4f %.4f"
            % (math.log10(max(NS)), max(cp, cl)))
        if max(cp, cl) > math.log10(max(NS)):
            say("FORECAST OUTSIDE provableforecast")
    if cp is not None and cl is not None and base is not None:
        say("FORECAST BOTH provableforecast %.4f %.4f %.4f"
            % (base, min(cp, cl), max(cp, cl)))
        if not (min(cp, cl) <= base <= max(cp, cl)):
            say("FORECAST CONVENTION SPLIT provableforecast")
    say()
    say("  and what the two answers are answers TO, which is not the")
    say("  same question. A is a measured maximum ratio, not a")
    say("  constant any theorem supplies: freezing it at the largest")
    say("  observed value is the conservative choice for an UPPER")
    say("  BOUND on the spend, and letting it fall is a description")
    say("  of what the data do, which would eventually put the")
    say("  'bound' below |P| itself and stop being a bound at all.")
    say("  So %.1f is where the bound provably pays and %.1f to %.1f"
        % (base, min(cp, cl), max(cp, cl)))
    say("  is where the measured ratio would; the two differ by")
    say("  %.0f decades and neither corrects the other."
        % (base - max(cp, cl)))

    say()
    say("  what that leaves. The published %.1f is a forecast about"
        % ppoint)
    say("  a frozen constant, and the constant is not one: it falls")
    say("  at %.2f standard errors and its shape is tied at %.2f of"
        % (abs(ap) / sp, abs(rp - rl) / ser))
    say("  the r.m.s.'s own error. Nothing here corrects the figure,")
    say("  because correcting it would need A's shape over four")
    say("  thousand decades -- an extrapolation this repository")
    say("  refuses on its own standard, {#rem:forecastbracket}.")

    say()
    say("=" * 70)
    ok = e1 and e2 and e3 and e4
    say("the forecast is undefined once the constant is allowed to "
        "move" if ok else "REFUTED")

    head = [
        "STATISTIC: the implied constant",
        "           A = max_k |P| / [(N/k) exp(-c sqrt(log(N/k))) L(k)]",
        "           at each N; its least-squares decay against log N",
        "           and against log log N with their r.m.s. residuals;",
        "           and the N solving A dL I(u) = S(N)(1-A(N)) with A",
        "           frozen at the sweep maximum and with A",
        "           extrapolated by each shape.",
        "NULL: none is run and none applies. A deterministic equation",
        "      is solved with measured inputs. The coin arm for this",
        "      field is lab_elementary_provable.py's sixteen coins on",
        "      the identical sifted set.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k odd, squarefree",
        "       and coprime to N with 2 <= k < " + str(KCAP) + "; m odd,",
        "       squarefree and coprime to k, m <= (N-1)/k; the sieve",
        "       weight over the odd primes below " + str(QSIEVE) + ";",
        "       c = " + str(CZERO) + " and theta' = " + str(THETA)
        + "; the Euler",
        "       products at the fixed bound " + str(CLIM) + "; the",
        "       crossing searched over log N up to " + str(UMAX) + ".",
        "       Every N is 2^a 5^b, one odd radical, as RADICALS says.",
        "       The published dL, A and bracket are read from",
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
