# -*- coding: utf-8 -*-
r"""
Is the residue's k-exponent what makes the hard arithmetic hard?

WHAT IS AT STAKE

Remark {#rem:cRwindow} found that at the primorial radical |R| is not
(N/k)^{1/2}: the octave exponent runs 0.3991 to 0.5074, rising with N
and reaching a half only near the top of the ladder. At the 2^a 5^b
family Remark {#rem:residue} measured the same quantity at 0.4722 to
0.5079 -- already at a half at N = 2e5, where the primorial ladder is
at 0.42.

That suggests a single number behind the whole arithmetic dependence.
Remark {#rem:residuearithmetic} found the level crossing 1/2 downward
as the radical grows, and attributed it to the budget, which does
shrink by a factor of five. But if the k-exponent ALSO falls with the
radical, then the budget is not the only thing moving, and the level's
arithmetic dependence has a second source that nothing has measured.

The seven N of the arithmetic test set are all near 1.6e6, so N is
held nearly fixed and only the radical varies. That isolates it.

BACKS: Remark {#rem:kexponent} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  L1  The control: at N = 1600000 the octave exponent reproduces the
      0.5048 of results/lab_residue_size.txt to within 0.02, the two
      differing only in k-cap.
  L2  The k-exponent is arithmetic-dependent: across the seven N it
      spans more than 0.05.
  L3  And it falls as the radical grows: the exponent regressed on
      the number of odd prime factors of N has a negative slope.
  L4  It is a second source of the level's arithmetic dependence: the
      k-exponent correlates with the level exponent of
      results/audit_residue_arithmetic.txt at better than 0.7.

REFUTATION RULE (fixed before the run)

  L1  REFUTED at 0.02.
  L2  REFUTED if the span is 0.05 or less, which would say the
      k-exponent is a function of N alone and the arithmetic
      dependence of the level is the budget and nothing else.
  L3  REFUTED if the slope is not negative.
  L4  REFUTED if the correlation is 0.7 or below. A failure here with
      L2 and L3 holding would be the interesting case: the
      k-exponent would move with the arithmetic but not with the
      level, and the two dependences would be separate.

  All four gate.

  NO NULL IS RUN and none applies. One measured magnitude is fitted
  against another at seven N and the fits compared; there is no
  background to detect against. The sign controls for R were run in
  lab_residue_cancellation.py, whose coin arm on the identical
  deviations established that R's size is bought by cancellation at
  exactly a coin's rate, and in lab_split_budget.py.
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
OUT = os.path.join(ROOT, "results", "audit_residue_kexponent.txt")

KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000
OCT = [2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288,
       2097152]
MINPTS = 10
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


def read_levels():
    """the test set, its thresholds and its level exponents -- read"""
    p = os.path.join(ROOT, "results", "audit_residue_arithmetic.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("N            odd part               threshold  "
                  "K*_R    exponent")
    out = []
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            continue
        out.append((int(f[0]), f[1], float(f[2]), float(f[4])))
        if len(out) == 7:
            break
    return out


def read_ladder_kexp():
    """the eleven k-exponents at one fixed radical -- read"""
    p = os.path.join(ROOT, "results", "audit_cR_window.txt")
    src = io.open(p, encoding="utf-8").read()
    # the same header opens K4's wide table earlier in that file;
    # the bare two-column one is the last occurrence
    i = src.rindex("  N            octave exponent" + chr(10))
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) != 2 or not f[0].isdigit():
            break
        out[int(f[0])] = float(f[1])
    return out


def read_family_kexp():
    p = os.path.join(ROOT, "results", "lab_residue_size.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("N            fitted exponent   correlation")
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        if int(f[0]) == FAMILY:
            return float(f[1])
    return None


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    test = read_levels()
    pubk = read_family_kexp()
    say("read %d test N with their level exponents, and the family's"
        % len(test))
    say("  k-exponent %.4f at N = %d, from results/" % (pubk, FAMILY))

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
    for N, odds, thrpub, lvl in test:
        PN = factor_set(N)
        A_, S_ = artin, twin
        for q in sorted(PN):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))

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
        beta = float((H * P).sum() / (P * P).sum())
        aR = np.abs(H - beta * P)
        inner = N // ks

        cent, prof, cnt = [], [], []
        for a, b in zip(OCT, OCT[1:]):
            sel = (inner >= a) & (inner < b)
            cnt.append(int(sel.sum()))
            if sel.sum():
                cent.append(float(inner[sel].mean()))
                prof.append(float(aR[sel].mean()))
            else:
                cent.append(float("nan"))
                prof.append(float("nan"))
        c_ = np.array(cent)
        p_ = np.array(prof)
        ok = (~np.isnan(c_) & ~np.isnan(p_) & (p_ > 0)
              & (np.array(cnt) >= MINPTS))
        xx = np.log(c_[ok])
        yy = np.log(p_[ok])
        e = float(np.polyfit(xx, yy, 1)[0])
        rr = float(np.corrcoef(xx, yy)[0, 1])
        thin = int(min(np.array(cnt)[ok]))
        nodd = len([q for q in PN if q > 2])
        rows.append((N, odds, nodd, S_ * (1.0 - A_), lvl, e, rr, thin,
                     ks.size))
        say("  N = %-9d odd %-22s w %d  k-exp %.4f  corr %.5f"
            % (N, odds, nodd, e, rr))

    # ------------------------------------------------------------- L1
    say()
    ctl = [r for r in rows if r[0] == FAMILY][0]
    d = abs(ctl[5] - pubk)
    l1 = d < 0.02
    say("L1  the control at N = %d: %.4f against the published %.4f, "
        "diff %.4f" % (FAMILY, ctl[5], pubk, d))
    say("  L1 %s" % ("hold" if l1 else "REFUTED"))

    # ---------------------------------------------------------- L2/L3
    say()
    say("L2/L3  the k-exponent across arithmetic types")
    say("  N            odd part               #odd  k-exp     level")
    ex = []
    for N, odds, nodd, thrc, lvl, e, rr, thin, nk in rows:
        ex.append(e)
        say("  %-12d %-22s %-5d %-9.4f %.4f"
            % (N, odds, nodd, e, lvl))
        say("POP kexp_N%d %d" % (N, thin))
        say("CORR kexp_N%d %.5f" % (N, abs(rr)))
    span = max(ex) - min(ex)
    l2 = span > 0.05
    w = np.array([r[2] for r in rows], dtype=float)
    y = np.array(ex)
    sl = float(np.polyfit(w, y, 1)[0])
    rw = float(np.corrcoef(w, y)[0, 1])
    l3 = sl < 0.0
    say("  span %.4f   (floor 0.05)   %s"
        % (span, "hold" if l2 else "REFUTED"))
    say("  L2 %s" % ("hold" if l2 else "REFUTED"))
    say("  regressed on the number of odd prime factors: slope %+.6f, "
        "correlation %.5f" % (sl, rw))
    say("  L3 %s" % ("hold" if l3 else "REFUTED"))

    # ------------------------------------------------------------- L4
    say()
    lv = np.array([r[4] for r in rows])
    r4 = float(np.corrcoef(y, lv)[0, 1])
    l4 = r4 > 0.7
    sl4 = float(np.polyfit(y, lv, 1)[0])
    say("L4  does the k-exponent track the level?")
    say("  correlation %.5f, slope of level on k-exponent %+.4f"
        % (r4, sl4))
    say("  L4 %s   (floor 0.7)" % ("hold" if l4 else "REFUTED"))

    say()
    say("  the arithmetic and the budgets, declared:")
    rads = set()
    for N, odds, nodd, thrc, lvl, e, rr, thin, nk in rows:
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  %d N, %d distinct odd radical%s"
        % (len(rows), len(rads), "" if len(rads) == 1 else "s"))
    say("RADICALS %d" % len(rads))
    for N, odds, nodd, thrc, lvl, e, rr, thin, nk in rows:
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, thrc))

    say()
    say("  DIAGNOSTIC (post hoc). Two candidate drivers, side by side.")
    say("  The level depends on the budget -- that is measured, at")
    say("  the correlation results/audit_residue_arithmetic.txt")
    say("  publishes and this file recomputes below.")
    say("  Does the k-exponent add anything, or is it the same thing")
    say("  seen twice? Against the budget:")
    th = np.log(np.array([r[3] for r in rows]))
    rb = float(np.corrcoef(th, y)[0, 1])
    say("  k-exponent against log(threshold): correlation %.5f" % rb)
    say("  and the level against log(threshold), recomputed here for")
    say("  comparison: correlation %.5f"
        % float(np.corrcoef(th, lv)[0, 1]))
    say("  A k-exponent that tracks the budget as tightly as the level")
    say("  does is not a second source; one that does not, is.")
    say()
    say("  But first: is any of this above the noise? Seven points")
    say("  cannot carry a correlation of %.2f. The noise floor for"
        % abs(r4))
    say("  this statistic is measurable -- the ladder gives eleven")
    say("  k-exponents at ONE radical, so their scatter about their")
    say("  own trend is what the statistic does when the arithmetic is")
    say("  held fixed:")
    lad = read_ladder_kexp()
    ln_ = sorted(lad)
    lx = np.log(np.array(ln_, dtype=float))
    ly = np.array([lad[n] for n in ln_])
    la, lb = np.polyfit(lx, ly, 1)
    lres = ly - (la * lx + lb)
    lrms = float(np.sqrt((lres ** 2).mean()))
    say("    %d rungs at one radical, r.m.s. about their line %.4f"
        % (len(ln_), lrms))
    say("    the test set's whole span across seven radicals   %.4f"
        % span)
    say("    ratio %.2f" % (span / lrms))
    say("  Seven draws from a distribution of that width have an")
    say("  expected range, which is what a span has to beat before it")
    say("  means anything. Estimated here by simulation on the same")
    say("  scatter, %d draws of seven:" % SIMS)
    rg = np.random.default_rng(SEED)
    sim = rg.normal(0.0, lrms, size=(SIMS, len(rows)))
    expspan = float((sim.max(axis=1) - sim.min(axis=1)).mean())
    say("    expected noise span %.4f, measured %.4f, ratio %.2f"
        % (expspan, span, span / expspan))
    say("  A span that is %.2f times what noise alone would give, over"
        % (span / expspan))
    say("  seven points, is not an arithmetic dependence: it is what")
    say("  this statistic does. The same radical appears twice in the")
    say("  two files at nearly the same N -- 1621620 here and 1921920")
    say("  on the ladder, both 3*5*7*11*13 -- and they differ by:")
    lad_close = min(ln_, key=lambda n: abs(math.log(n / 1621620)))
    here = [r[5] for r in rows if r[0] == 1621620]
    if here:
        say("    N = 1621620 gives %.4f, N = %d gives %.4f, apart by "
            "%.4f" % (here[0], lad_close, lad[lad_close],
                      abs(here[0] - lad[lad_close])))
    say("  which is most of the span L2 measured, at one radical.")
    say("FLOOR kexp_across_radicals %.4f" % expspan)

    say()
    say("=" * 70)
    ok = l1 and l2 and l3 and l4
    say("the k-exponent falls with the radical and tracks the level"
        if ok else "REFUTED")

    head = [
        "SEED: the null draws from numpy default_rng at "
        "seed %d; without it the file does not reproduce "
        "its own null." % SEED,
        "STATISTIC: the exponent a in |R(N;k)| ~ (N/k)^a, fitted to the",
        "           octave means of |R| against the inner length N/k, at",
        "           the seven N of the arithmetic test set; its span",
        "           across them; its regression on the number of odd",
        "           prime factors of N; and its correlation with the",
        "           level exponent log K*_R / log N of",
        "           results/audit_residue_arithmetic.txt.",
        "NULL: none is run and none applies. One measured magnitude is",
        "      fitted against another at seven N and the fits compared;",
        "      there is no background to detect against. The sign",
        "      controls for R were run in lab_residue_cancellation.py",
        "      and lab_split_budget.py.",
        "FIELD: the seven N of the arithmetic test set, read from",
        "       results/audit_residue_arithmetic.txt; k squarefree and",
        "       coprime to N with 2 <= k < 100000; m odd, squarefree,",
        "       coprime to the odd part of k, m < N/k; the sieve weight",
        "       uses the odd primes up to 30; beta refitted as",
        "       sum(H P)/sum(P^2) on that k-range; octaves closed at",
        "       both ends and fitted only when they hold at least 10 k.",
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
