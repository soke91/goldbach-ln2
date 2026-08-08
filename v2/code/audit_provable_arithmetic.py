# -*- coding: utf-8 -*-
r"""
The one exposure {#rem:provablehalf} names and does not test.

WHAT IS AT STAKE

Remark {#rem:provablehalf} prices the conditional reduction: if the
classical bound is used for the elementary half, it pays the budget
only at N = 10^5474.8. The forecast rests on three measured inputs --
the density-times-L factor dL, the implied constant A forced by the
data, and the budget constant S(N)(1-A(N)) -- and the remark brackets
only the analytic exponent c, over [10^2092.7, 10^13093.3], plus a
sensitivity sweep of A dL by a factor 2.1071 giving
[10^4838.5, 10^6139.9].

It also says, of dL, "the first does not drift at all here -- every N
in the sweep has the same odd radical, so the admissible k-set and dL
are identical across it, **a fact about this sweep and not a general
one**". That is an exposure named and left untested, and Remark
{#rem:arithmeticreach} has since shown the arithmetic dependence of
the level to be standing rather than closing, so it is worth testing.

Across arithmetic types both dL and the budget move, and they move
against each other: a primorial-like N admits fewer k, which lowers
A dL and brings the crossing in, and carries a thinner budget, which
pushes it out. Which wins is not readable off either number.

The implementation is the bitmask one of audit_level_slope_reach.py,
independent of lab_elementary_provable.py's residue-class masking, so
Z1 is a cross-check of the published dL, A and forecast.

BACKS: Remark {#rem:provablearithmetic} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  The control: at the family member N = 1600000 this
      implementation reproduces the published dL to within 0.001 and,
      with the published A, the forecast to within one decade.
  Z2  dL is not arithmetic-independent: its relative spread across the
      seven types exceeds 0.05, against the 0.0000 drift the family
      sweep reports for it.
  Z3  The arithmetic matters less than the analytic exponent: every
      type's forecast lies inside the published c-sweep bracket
      [10^2092.7, 10^13093.3].
  Z4  But more than the sensitivity the remark declares: the seven do
      not all lie inside the A dL sweep's [10^4838.5, 10^6139.9].

REFUTATION RULE (fixed before the run)

  Z1  REFUTED at 0.001 in dL or one decade in the forecast -- not the
      same statistic, and nothing below may be compared with
      {#rem:provablehalf}.
  Z2  REFUTED below 0.05. dL would then be a property of the cap and
      not of N's arithmetic, and the caveat the remark attaches to it
      would be empty.
  Z3  REFUTED if any type falls outside the c-sweep bracket. The
      arithmetic would then matter more than the classical constant,
      and no bracket in these papers would cover the forecast.
  Z4  REFUTED if all seven lie inside the A dL bracket. The declared
      sensitivity sweep would already cover the arithmetic and nothing
      would need to change.

  All four gate.

  NO NULL IS RUN and none applies to Z1-Z4. A deterministic bound is
  evaluated on a measured summation set and a crossing located; there
  is no background to detect against. The coin arm for this field is
  lab_elementary_provable.py's own sixteen coins on the identical
  sifted set, which established that the looseness is the shape and
  not mu; it is not repeated here.
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
OUT = os.path.join(RES, "audit_provable_arithmetic.txt")

FAMILY = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
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


def read_types():
    """the seven arithmetic test N and their odd parts"""
    src = io.open(os.path.join(RES, "audit_residue_arithmetic.txt"),
                  encoding="utf-8").read()
    i = src.index("N            odd part               threshold  "
                  "K*_R    exponent  clears .5")
    ns, odd = [], {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            break
        ns.append(int(f[0]))
        odd[int(f[0])] = f[1]
    return ns, odd


def read_published():
    """dL, A and the two brackets, read from the results file"""
    src = io.open(os.path.join(RES, "lab_elementary_provable.txt"),
                  encoding="utf-8").read()
    dl = float(re.search(r"density-times-L factor dL = ([\d.]+)",
                         src).group(1))
    a = float(re.search(r"constant A forced by the data = ([\d.]+)",
                        src).group(1))
    m = re.search(r"^BRACKET log10_N_elementary_provable "
                  r"([\d.]+) ([\d.]+) ([\d.]+)\s*$", src, re.M)
    point, clo, chi = (float(m.group(1)), float(m.group(2)),
                       float(m.group(3)))
    i = src.index("A dL scaled by        log10 N")
    vals = []
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) != 2:
            break
        try:
            vals.append((float(f[0]), float(f[1])))
        except ValueError:
            break
    alo = min(v for _s, v in vals)
    ahi = max(v for _s, v in vals)
    return dl, a, point, clo, chi, alo, ahi


def inputs(N, mu, oddsqf, vmask, qs, artin, twin):
    """dL, A and the budget constant at one N"""
    PN = factor_set(N)
    A_, S_ = artin, twin
    for q in sorted(PN):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))

    ks, Pv, Lv = [], [], []
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
    ks = np.array(ks, dtype=np.int64)
    Pv = np.array(Pv)
    Lv = np.array(Lv)
    inner = (N // ks).astype(np.float64)
    b = inner * np.exp(-CZERO * np.sqrt(np.log(inner))) * Lv
    w = np.log(ks.astype(np.float64)) / ks.astype(np.float64)
    dl = float((w * Lv).sum()
               / (math.log(float(ks.max())) ** 2 / 2.0))
    amax = float((Pv / b).max())
    return dl, amax, S_ * (1.0 - A_), ks.size


def integral(u, c):
    v = np.linspace(0.0, THETA * u, 20001)
    return float(np.trapezoid(v * np.exp(-c * np.sqrt(u - v)), v))


def solve(prod, thrc, c):
    """where A dL I(u) falls back through the budget constant"""
    lo, hi = 10.0, 1e7
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if prod * integral(mid, c) > thrc:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi) / math.log(10.0)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    base, odd = read_types()
    pdl, pa, ppoint, pclo, pchi, palo, pahi = read_published()
    say("read %d test N from results/audit_residue_arithmetic.txt"
        % len(base))
    say("  and dL %.4f, A %.4f, forecast %.1f, c-bracket [%.1f, %.1f],"
        % (pdl, pa, ppoint, pclo, pchi))
    say("  A dL bracket [%.1f, %.1f] from "
        "results/lab_elementary_provable.txt" % (palo, pahi))

    NMAX = max(max(base), max(FAMILY))
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

    got = []
    for N in base:
        dl, amax, thrc, nk = inputs(
            N, mu, oddsqf, vmask, qs, artin, twin)
        got.append((N, dl, amax, thrc, nk))
        say("  N = %-10d odd %-22s #k %-6d dL %.4f  A %.4f  "
            "S(1-A) %.6f" % (N, odd[N], nk, dl, amax, thrc))
        say("BUDGET provable_S1A_N%d %.6f" % (N, thrc))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(g[0]) if q > 2))
                  for g in got)))

    # ------------------------------------------------------------- Z1
    say()
    say("Z1  the control: the family member N = 1600000")
    fam = [g for g in got if g[0] == 1_600_000][0]
    ddl = abs(fam[1] - pdl)
    fc = solve(pa * fam[1], fam[3], CZERO)
    dfc = abs(fc - ppoint)
    z1 = (ddl < 0.001) and (dfc < 1.0)
    say("  dL %.4f against the published %.4f, diff %.5f"
        % (fam[1], pdl, ddl))
    say("  forecast with the published A: log10 N = %.1f against "
        "%.1f, diff %.2f" % (fc, ppoint, dfc))
    say("  Z1 %s   (cap 0.001 in dL, cap 1 decade in the forecast)"
        % ("hold" if z1 else "REFUTED"))

    # ------------------------------------------------------------- Z2
    say()
    say("Z2  does dL depend on the arithmetic?")
    dls = [g[1] for g in got]
    rel = (max(dls) - min(dls)) / float(np.mean(dls))
    z2 = rel > 0.05
    say("  dL by type: " + ", ".join("%.4f" % d for d in dls))
    say("  and the same statistic with the arithmetic held fixed --")
    say("  dL at the five N of the family sweep, computed here:")
    fdls = []
    for N in FAMILY:
        d1, _a1, _t1, _n1 = inputs(
            N, mu, oddsqf, vmask, qs, artin, twin)
        fdls.append(d1)
        say("    N = %-10d dL %.4f" % (N, d1))
    frel = (max(fdls) - min(fdls)) / float(np.mean(fdls))
    say("  spread %.4f   (floor %.4f)" % (rel, frel))
    say("FLOOR provable_dL_across_radicals %.4f" % frel)
    say("SCALES audit_provable_arithmetic 1")
    say("ONE SCALE audit_provable_arithmetic")
    say("  the k-set is fixed by N's odd radical and the cap, and the")
    say("  cap binds at every N here, so with the radical held the")
    say("  admissible k are literally the same set and dL cannot")
    say("  move: the floor is exact and not a simulation.")
    say("DRIFT provable_dL_across_radicals %.4f" % rel)
    say("  Z2 %s   (cap 0.05)" % ("hold" if z2 else "REFUTED"))

    # ---------------------------------------------------------- Z3/Z4
    say()
    say("Z3/Z4  the forecast, with every input measured per type")
    say("  odd part               dL       A        S(1-A)     "
        "log10 N")
    fcs = []
    for N, dl, amax, thrc, nk in got:
        v = solve(amax * dl, thrc, CZERO)
        fcs.append(v)
        say("  %-22s %-8.4f %-8.4f %-10.6f %.1f"
            % (odd[N], dl, amax, thrc, v))
    z3 = all(pclo <= v <= pchi for v in fcs)
    z4 = not all(palo <= v <= pahi for v in fcs)
    say("  the seven run %.1f to %.1f, a span of %.1f decades"
        % (min(fcs), max(fcs), max(fcs) - min(fcs)))
    say("  Z3 all inside the c-sweep bracket [%.1f, %.1f]   %s"
        % (pclo, pchi, "hold" if z3 else "REFUTED"))
    say("  Z4 not all inside the A dL bracket [%.1f, %.1f]   %s"
        % (palo, pahi, "hold" if z4 else "REFUTED"))
    out = [odd[g[0]] for g, v in zip(got, fcs)
           if not (palo <= v <= pahi)]
    say("  outside it: %s" % (", ".join(out) if out else "none"))
    say("BRACKET log10_N_provable_across_radicals %.4f %.4f %.4f"
        % (fcs[[g[0] for g in got].index(1_600_000)],
           min(fcs), max(fcs)))
    say("DRIFT log10_N_provable_across_radicals %.4f"
        % ((max(fcs) - min(fcs)) / float(np.mean(fcs))))
    say()
    say("  A held at the published A, which is the max over the whole")
    say("  family sweep and not a per-N quantity, so this row isolates")
    say("  the two channels the remark names -- dL and the budget:")
    say("  odd part               log10 N (A fixed)")
    fixed = []
    for N, dl, amax, thrc, nk in got:
        v = solve(pa * dl, thrc, CZERO)
        fixed.append(v)
        say("  %-22s %.1f" % (odd[N], v))
    say("  those run %.1f to %.1f, a span of %.1f decades"
        % (min(fixed), max(fixed), max(fixed) - min(fixed)))
    lo = min(min(fcs), min(fixed))
    hi = max(max(fcs), max(fixed))
    say("  so over both conventions the forecast runs %.1f to %.1f"
        % (lo, hi))
    say("ACROSS log10_N_elementary_provable %.4f %.4f" % (lo, hi))
    say("SENSITIVITY log10_N_elementary_provable %.4f %.4f"
        % (palo, pahi))
    if not (palo <= lo and hi <= pahi):
        say("SENSITIVITY UNDERSTATED log10_N_elementary_provable")

    say()
    say("  DIAGNOSTIC (post hoc). The two channels pull opposite")
    say("  ways, so the span is smaller than either. Against the")
    say("  family member, per type:")
    say("  odd part               A dL / family   S(1-A) / family   "
        "decades")
    fdl, fa, fth = fam[1], fam[2], fam[3]
    ffc = solve(fa * fdl, fth, CZERO)
    for (N, dl, amax, thrc, nk), v in zip(got, fcs):
        say("  %-22s %-15.4f %-17.4f %+.1f"
            % (odd[N], (amax * dl) / (fa * fdl), thrc / fth, v - ffc))
    say("  A thinner budget pushes the crossing out and a thinner")
    say("  k-set brings it in; at the primorial-like types both are")
    say("  thin at once.")

    say()
    say("=" * 70)
    ok = z1 and z2 and z3 and z4
    say("the arithmetic moves the forecast past its declared "
        "sensitivity" if ok else "REFUTED")

    head = [
        "STATISTIC: the three measured inputs to",
        "           {#rem:provablehalf}'s forecast -- the",
        "           density-times-L factor dL = sum L(k) log k / k",
        "           over the admissible k against log(kmax)^2/2, the",
        "           constant A = max_k |P|/[(N/k) exp(-c sqrt(log(N/k)))",
        "           L(k)], and the budget constant S(N)(1-A(N)) -- at",
        "           each of the seven arithmetic types of",
        "           {#rem:residuearithmetic}, and the N at which",
        "           A dL int_0^{theta' u} v exp(-c sqrt(u-v)) dv falls",
        "           back through that budget constant.",
        "NULL: none is run and none applies. A deterministic bound is",
        "      evaluated on a measured summation set and a crossing",
        "      located; there is no background to detect against. The",
        "      coin arm for this field is lab_elementary_provable.py's",
        "      sixteen coins on the identical sifted set.",
        "FIELD: the seven test N read from",
        "       results/audit_residue_arithmetic.txt; k squarefree,",
        "       odd and coprime to N with 2 <= k < " + str(KCAP) + ",",
        "       the cap lab_elementary_provable.py uses; m odd,",
        "       squarefree and coprime to k, m <= (N-1)/k; the sieve",
        "       weight over the odd primes below " + str(QSIEVE) + ";",
        "       c = " + str(CZERO) + " and theta' = " + str(THETA)
        + "; the Euler",
        "       products at the fixed bound " + str(CLIM) + ". Seven",
        "       distinct odd radicals, as RADICALS declares. The",
        "       published dL, A and brackets are read from",
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
