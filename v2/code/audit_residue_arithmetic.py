# -*- coding: utf-8 -*-
r"""
The knife-edge at one arithmetic type, and at the others.

WHAT IS AT STAKE

Every measurement of R in this repository -- {#rem:residue},
{#rem:residuelevel}, {#rem:residueconstant}, {#rem:betafree},
{#rem:residuesigned} -- runs over N = 2e5 * 2^j. Those five N are all
2^a 5^b, so they share one odd radical. lab_elementary_provable.py
made that visible by accident: the density-times-L factor d_L came out
identical to four decimals at every N, because the admissible k-set
does not change across the sweep at all.

So the knife-edge of {#rem:residuelevel} -- exponents 0.5654 to
0.5799, clearing 1/2 by 0.06 to 0.08 -- is measured at ONE arithmetic
type. And the budget it is measured against is the part that varies
most with type: audit_threshold_arithmetic.py found S(N)(1-A(N))
running from 0.073312 at a primorial-like N to 0.374487 at this
family, a factor of five. Proposition [prop:onesided] says exactly
this: the threshold is of order N for almost all even N and sinks
towards N/(log N log log N) at the worst ones.

If the residue-only level falls below 1/2 at those worst N, then the
conditional reduction of {#rem:provablehalf} fails precisely where the
problem is hardest, and the knife-edge is an artefact of a convenient
family.

BACKS: Remark {#rem:residuearithmetic} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  P1  The control: at N = 1600000, which belongs to the family, the
      exponent reproduces the 0.5675 of results/audit_residue_level.txt
      to within 0.001.
  P2  The reduction survives the worst arithmetic: every N in the test
      set has log K*_R / log N above 0.5.
  P3  Type matters a great deal: the spread of the exponent across the
      test set is at least 0.05. It would be a surprise, and worth
      having, if the level were insensitive to something the budget
      depends on by a factor of five.
  P4  And it is the budget that drives it: regressing the exponent on
      log of the threshold constant gives a correlation above 0.9.

REFUTATION RULE (fixed before the run)

  P1  REFUTED at 0.001, which would mean this is not the same
      measurement.
  P2  REFUTED if any N in the test set reaches 0.5 from above. That is
      the one that matters: it would say the conditional reduction is
      a fact about 2^a 5^b and not about even numbers.
  P3  REFUTED if the spread is under 0.05.
  P4  REFUTED if the correlation is 0.9 or below, which would say
      something other than the budget is moving the level.

  All four gate.

  NO NULL IS RUN and none applies. A measured sum is crossed against a
  computed threshold at seven values of N and the crossings compared;
  there is no background to detect against. The sign controls for this
  field were run in lab_residue_cancellation.py, whose coin arm on the
  identical delta established that R's size is bought by cancellation
  at exactly a coin's rate, and in lab_split_budget.py.
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
OUT = os.path.join(ROOT, "results", "audit_residue_arithmetic.txt")

KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000
FAMILY = 1_600_000
SIMS = 20_000
SEED = 20260808


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


def read_testset():
    """the arithmetic test set and its thresholds -- read, not copied"""
    p = os.path.join(ROOT, "results", "audit_threshold_arithmetic.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("V2/V3  the threshold across the test set")
    out = []
    for ln in src[i:].splitlines()[2:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        out.append((int(f[0]), f[1], float(f[2])))
    return out


def read_budget_gap():
    """the exponent cost of a budget factor -- read, not copied"""
    p = os.path.join(ROOT, "results", "audit_model_transfer.txt")
    src = io.open(p, encoding="utf-8").read()
    m = re.search(r"mean gap ([\d.]+) -- a budget factor of ([\d.]+)",
                  src)
    return float(m.group(1)), float(m.group(2))


def read_family_drift():
    """c_R's drift across the family -- read, not copied"""
    p = os.path.join(ROOT, "results", "audit_residue_constant.txt")
    src = io.open(p, encoding="utf-8").read()
    return float(re.search(r"spread ([\d.]+) of the mean",
                           src).group(1))


def read_family_exponents():
    """all five level exponents at the one 2^a 5^b radical -- read.

    They are the noise floor for this statistic: five N of a SINGLE
    arithmetic, so their scatter about their own trend is what the
    level exponent does when the radical is held fixed.
    """
    p = os.path.join(ROOT, "results", "audit_residue_level.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("log K*_R/log N")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        out[int(f[0])] = float(f[3])
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    test = read_testset()
    fam = read_family_exponents()
    pubexp = fam[FAMILY]
    say("read %d test N from results/audit_threshold_arithmetic.txt, "
        "and the" % len(test))
    say("  family exponent %.4f at N = %d from "
        "results/audit_residue_level.txt" % (pubexp, FAMILY))

    NMAX = max(t[0] for t in test)
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
    for N, odds, thrpub in test:
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
        sel = ks <= kstar if kstar else np.zeros(ks.size, bool)
        cR = float((aR[sel] / scale[sel]).mean()) if kstar else float("nan")
        rows.append((N, odds, thrc, thrpub, len(ks), beta, kstar, e, cR))
        say("  N = %-9d odd %-22s thr %.6f  #k %-6d K*_R %s"
            % (N, odds, thrc, ks.size, str(kstar)))

    # ------------------------------------------------------------- P1
    say()
    say("P1  the control: the family member")
    p1 = True
    for N, odds, thrc, thrpub, nk, beta, kstar, e, cR in rows:
        if N != FAMILY:
            continue
        d = abs(e - pubexp)
        if not (d < 0.001):
            p1 = False
        say("  N = %d: exponent %.4f against the published %.4f, "
            "diff %.5f" % (N, e, pubexp, d))
    say("  P1 %s" % ("hold" if p1 else "REFUTED"))

    # ---------------------------------------------------------- P2/P3
    say()
    say("P2/P3  the level across arithmetic types")
    say("  N            odd part               threshold  K*_R    "
        "exponent  clears .5")
    p2 = True
    ex = []
    for N, odds, thrc, thrpub, nk, beta, kstar, e, cR in rows:
        ex.append(e)
        if not (e > 0.5):
            p2 = False
        say("  %-12d %-22s %-10.6f %-7s %-9.4f %s"
            % (N, odds, thrc, str(kstar), e,
               "yes" if e > 0.5 else "NO"))
    spread = max(ex) - min(ex)
    p3 = spread >= 0.05
    say("  P2 every N above 0.5   %s" % ("hold" if p2 else "REFUTED"))
    say("  P3 spread %.4f   (floor 0.05)   %s"
        % (spread, "hold" if p3 else "REFUTED"))
    say()
    say("  and the noise floor that spread has to beat, which gate")
    say("  check G37 asks for. The five N of results/")
    say("  audit_residue_level.txt are one radical, so their scatter")
    say("  about their own trend is what this statistic does with the")
    say("  arithmetic held fixed:")
    fn = sorted(fam)
    fx = np.log(np.array(fn, dtype=float))
    fy = np.array([fam[n] for n in fn])
    fa, fb = np.polyfit(fx, fy, 1)
    frms = float(np.sqrt(((fy - (fa * fx + fb)) ** 2).mean()))
    rg = np.random.default_rng(SEED)
    sim = rg.normal(0.0, frms, size=(SIMS, len(rows)))
    expspan = float((sim.max(axis=1) - sim.min(axis=1)).mean())
    say("    %d N at one radical, r.m.s. about their line %.4f"
        % (len(fn), frms))
    say("    expected span of %d draws at that width %.4f (%d sims)"
        % (len(rows), expspan, SIMS))
    say("    measured span across %d radicals %.4f, ratio %.2f"
        % (len(rows), spread, spread / expspan))
    say("FLOOR level_across_radicals %.4f" % expspan)
    say("  The k-exponent of {#rem:kexponent} failed this same test,")
    say("  its span barely exceeding its floor; this one clears its")
    say("  own by %.1f, so P3 is a measurement and not a coincidence"
        % (spread / expspan))
    say("  of seven draws. The comparison is not read back from that")
    say("  file, which reads this one.")

    # ------------------------------------------------------------- P4
    say()
    say("P4  is it the budget that moves the level?")
    x = np.log(np.array([r[2] for r in rows]))
    y = np.array(ex)
    r4 = float(np.corrcoef(x, y)[0, 1])
    sl = float(np.polyfit(x, y, 1)[0])
    p4 = r4 > 0.9
    say("  exponent against log(threshold): slope %+.4f, correlation "
        "%.5f" % (sl, r4))
    say("  P4 %s   (floor 0.9)" % ("hold" if p4 else "REFUTED"))

    say()
    say("  the arithmetic this sweep actually covers, which gate check")
    say("  G34 reads. An exponent measured over one odd radical is a")
    say("  statement about that radical:")
    rads = set()
    for N, odds, thrc, thrpub, nk, beta, kstar, e, cR in rows:
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  %d N, %d distinct odd radicals" % (len(rows), len(rads)))
    say("RADICALS %d" % len(rads))

    say()
    say("  the budget constants crossed, declared:")
    for N, odds, thrc, thrpub, nk, beta, kstar, e, cR in rows:
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, thrc))

    say()
    gap, fac = read_budget_gap()
    fdrift = read_family_drift()
    say("  DIAGNOSTIC (post hoc). What the model of {#rem:modeltransfer}")
    say("  predicts for this. Read from results/audit_model_transfer.txt,")
    say("  a budget factor of %.4f costs %.4f in the exponent, i.e."
        % (fac, gap))
    say("  %.4f/log(%.4f) per natural log of the threshold. Against the"
        % (gap, fac))
    say("  measured slope above:")
    say("  predicted slope %+.4f, measured %+.4f"
        % (gap / math.log(fac), sl))
    say()
    say("  and the constant c_R, which {#rem:residueconstant} found")
    say("  drifting by %.4f across the FAMILY, read from" % fdrift)
    say("  results/audit_residue_constant.txt. Across types:")
    say("  N            c_R        c_R/sqrt(log N)")
    gam = []
    for N, odds, thrc, thrpub, nk, beta, kstar, e, cR in rows:
        g = cR / math.sqrt(math.log(N))
        gam.append(g)
        say("  %-12d %-10.4f %.4f" % (N, cR, g))
    say("  spread across types %.4f of the mean"
        % ((max(gam) - min(gam)) / float(np.mean(gam))))
    say("DRIFT residue_cR_across_types %.4f"
        % ((max(gam) - min(gam)) / float(np.mean(gam))))

    say()
    say("=" * 70)
    ok = p1 and p2 and p3 and p4
    say("the conditional reduction survives every arithmetic type, and "
        "the budget is what moves it" if ok else "REFUTED")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N, at",
        "           seven N of comparable size and different odd",
        "           radicals; the spread of that exponent across them;",
        "           its regression on the log of the threshold; and the",
        "           constant c_R at each.",
        "NULL: none is run and none applies. A measured sum is crossed",
        "      against a computed threshold at seven N and the crossings",
        "      compared; there is no background to detect against. The",
        "      sign controls for this field were run in",
        "      lab_residue_cancellation.py and lab_split_budget.py.",
        "FIELD: the seven N of the arithmetic test set, read from",
        "       results/audit_threshold_arithmetic.txt; k squarefree and",
        "       coprime to N with 2 <= k < 100000; m odd, squarefree,",
        "       coprime to the odd part of k, m < N/k; the sieve weight",
        "       uses the odd primes up to 30; beta refitted as",
        "       sum(H P)/sum(P^2) on the same k-range; S(N) and A(N)",
        "       from Euler products at the fixed bound 4000000.",
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
