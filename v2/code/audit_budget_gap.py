# -*- coding: utf-8 -*-
r"""
Is the budget gap a constant, and is it the same for both halves?

WHAT IS AT STAKE

Remark {#rem:modeltransfer} fixes what the operative budget costs:
"a budget factor of 4.7009 costs 0.1677 in the exponent". That number
is read out of results/audit_model_transfer.txt by
audit_residue_level.py and quoted in the papers as a constant of the
method. It is a MEAN. The series behind it is

    0.1824, 0.1688, 0.1706, 0.1635, 0.1532,

which falls across the sweep by more than a tenth of itself, and
nothing has tested whether that fall is above its own noise.

Two things follow if it is. First, a constant quoted from a drifting
series is a statement about the middle of the accessible range, and
{#rem:residuearithmetic} has already used it as a model prediction
(+0.1084 against a measured +0.0516). Second, the gap was measured on
H, and every use of it is about R. Whether the transfer is a property
of the budget ratio -- in which case both halves pay the same -- or of
the half, has never been asked; the two crossings needed are already
published in different files and have never been differenced.

The implementation is the bitmask one of audit_level_slope_reach.py,
independent of audit_model_transfer.py's, so Y1 is a cross-check of
the published gap and not a rerun.

BACKS: Remark {#rem:budgetgap} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  The control: the H gap reproduces the published
      0.1824, 0.1688, 0.1706, 0.1635, 0.1532 to within 0.001.
  Y2  It is not a constant: the gap's least-squares slope against
      log N is negative and reaches two standard errors.
  Y3  The transfer belongs to the budget and not to the half: the
      same gap computed for R agrees with H's to within 0.01 at
      every N.
  Y4  And the two decline together: the difference of the two slopes
      is inside two standard errors of zero.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED at 0.001 at any N -- not the same statistic, and
      nothing below may be compared with {#rem:modeltransfer}.
  Y2  REFUTED if the slope fails to reach two standard errors. The
      gap would then be a constant to the precision this sweep can
      see, and quoting the mean would be right.
  Y3  REFUTED if the two gaps differ by more than 0.01 at any N.
      Then what the budget factor costs depends on which half is
      being truncated, and 0.1677 -- measured on H -- may not be
      used for R at all, which is what every use of it does.
  Y4  REFUTED if the slopes differ by more than two standard errors
      of their difference. The two halves would be losing the gap at
      different rates and no single gap could be quoted for both,
      even at one N.

  All four gate.

  NO NULL IS RUN and none applies. Two crossings of the same measured
  sum against two computed thresholds are located and differenced;
  there is no background to detect against. The sign controls for this
  field were run in lab_direct_level.py, whose mu-squared reference
  established that the level is bought by cancellation at all, and in
  lab_split_budget.py.
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
OUT = os.path.join(RES, "audit_budget_gap.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
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


def read_published_gap():
    """the H gap at each N, read from the results file"""
    src = io.open(os.path.join(RES, "audit_model_transfer.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log K*(S N)/log N   "
                  "log K*(S(1-A)N)/log N   gap")
    g = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        g[int(f[0])] = float(f[3])
    m = re.search(r"^\s*mean gap ([\d.]+)", src, re.M)
    return g, float(m.group(1))


def gaps(N, lam, mu, sqf, vmask, qs, artin, twin):
    """the exponent gap between the two budgets, for H and for R"""
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
    w = np.log(ks.astype(np.float64))
    big, small = S_ * N, S_ * (1.0 - A_) * N

    def exps(v):
        c = np.cumsum(w * np.abs(v))
        out = []
        for thr in (big, small):
            j = int(np.searchsorted(c, thr))
            out.append(None if j >= ks.size
                       else math.log(int(ks[j])) / math.log(N))
        return out

    eh = exps(H)
    er = exps(R)
    return eh, er, 1.0 / (1.0 - A_), small / N


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    rms = float(np.sqrt((r ** 2).mean()))
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), rms, se, abs(float(a)) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub, pubmean = read_published_gap()
    say("read %d published gaps and the mean %.4f from "
        "results/audit_model_transfer.txt" % (len(pub), pubmean))

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
        eh, er, fac, bpn = gaps(
            N, lam, mu, sqf, vmask, qs, artin, twin)
        got.append((N, eh, er, fac))
        say("  N = %-10d budget factor %.4f   H %.4f -> %.4f   "
            "R %.4f -> %.4f" % (N, fac, eh[0], eh[1], er[0], er[1]))
        say("BUDGET kstar_budgetgap_S1AN_N%d %.6f" % (N, bpn))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(g[0]) if q > 2))
                  for g in got)))

    # ------------------------------------------------------------- Y1
    say()
    say("Y1  the control: the H gap")
    say("  N            here     published  diff")
    y1 = True
    gh = []
    for N, eh, er, fac in got:
        v = eh[0] - eh[1]
        gh.append(v)
        d = abs(v - pub[N])
        if not (d < 0.001):
            y1 = False
        say("  %-12d %-8.4f %-10.4f %.5f" % (N, v, pub[N], d))
    say("  Y1 %s   (cap 0.001)" % ("hold" if y1 else "REFUTED"))

    # ------------------------------------------------------------- Y2
    say()
    say("Y2  is the gap a constant?")
    x = np.log(np.array([g[0] for g in got], dtype=np.float64))
    yh = np.array(gh)
    ah, rmsh, seh, th = fit(x, yh)
    say("  the H gap's least-squares slope against log N = %+.6f" % ah)
    say("  r.m.s. residual %.4f, standard error %.6f, t = %.2f"
        % (rmsh, seh, th))
    say("SCATTER slope_audit_budget_gap %.4f" % rmsh)
    say("TSTAT slope_audit_budget_gap %.2f" % th)
    say("SPREAD slope_audit_budget_gap %.4f" % float(x.max() - x.min()))
    if th < 2.0:
        say("UNRESOLVED SIGN slope_audit_budget_gap")
    y2 = (ah < 0.0) and (th >= 2.0)
    say("  the series and the mean that is quoted from it:")
    say("SERIES gap %d %.4f %.4f" % (len(gh), min(gh), max(gh)))
    if th >= 2.0:
        say("DRIFTS gap %+.6f %.2f" % (ah, th))
    else:
        say("FLAT gap %.2f" % th)
    say("  quoted mean %.4f, series %.4f to %.4f, so the mean is "
        "%.4f" % (pubmean, min(gh), max(gh), pubmean))
    say("  above the top of the range" if pubmean > max(gh) else
        ("  below the bottom of the range" if pubmean < min(gh)
         else "  inside the range but attained at no listed N"
              if pubmean not in gh else "  attained"))
    say("  Y2 the gap falls and clears two s.e.   %s"
        % ("hold" if y2 else "REFUTED"))

    # ---------------------------------------------------------- Y3/Y4
    say()
    say("Y3/Y4  the same gap for the residue")
    say("  N            H gap    R gap    |difference|")
    gr = []
    y3 = True
    for (N, eh, er, fac), v in zip(got, gh):
        u = er[0] - er[1]
        gr.append(u)
        d = abs(u - v)
        if not (d < 0.01):
            y3 = False
        say("  %-12d %-8.4f %-8.4f %.4f" % (N, v, u, d))
    say("  Y3 the two gaps agree to 0.01 at every N   %s"
        % ("hold" if y3 else "REFUTED"))
    yr = np.array(gr)
    ar, rmsr, ser, tr = fit(x, yr)
    say("  the R gap's slope %+.6f, standard error %.6f, t = %.2f"
        % (ar, ser, tr))
    sed = math.sqrt(seh ** 2 + ser ** 2)
    y4 = abs(ah - ar) <= 2.0 * sed
    say("  slopes %+.6f and %+.6f, difference %+.6f against "
        "%.6f = 2 s.e." % (ah, ar, ah - ar, 2.0 * sed))
    say("  Y4 the two decline together   %s"
        % ("hold" if y4 else "REFUTED"))
    say("SERIES gapR %d %.4f %.4f" % (len(gr), min(gr), max(gr)))
    if tr >= 2.0:
        say("DRIFTS gapR %+.6f %.2f" % (ar, tr))
    else:
        say("FLAT gapR %.2f" % tr)

    say()
    say("  DIAGNOSTIC (post hoc). What this does to")
    say("  {#rem:residuearithmetic}'s model comparison. It divided")
    say("  the quoted mean by the log of the budget factor to get a")
    say("  predicted response per natural log of the threshold, and")
    say("  measured about half that across seven radicals. With the")
    say("  gap taken at each N instead of averaged, the prediction")
    say("  runs:")
    say("  N            gap      /log(factor)")
    for (N, eh, er, fac), v in zip(got, gh):
        say("  %-12d %-8.4f %.4f" % (N, v, v / math.log(fac)))
    say("  so the prediction the arithmetic sweep should have been")
    say("  compared with, at its own N = 1600000, is %.4f and not"
        % (gh[3] / math.log(got[3][3])))
    say("  the %.4f the mean gives. The disagreement with the"
        % (pubmean / math.log(got[3][3])))
    say("  measured cross-arithmetic response survives either way.")

    say()
    say("=" * 70)
    ok = y1 and y2 and y3 and y4
    say("the gap is a declining quantity and both halves pay it alike"
        if ok else "REFUTED")

    head = [
        "STATISTIC: for each of H and R, the crossings of",
        "           sum_{k<K}(log k)|.| against S(N)N and against",
        "           S(N)(1-A(N))N, their exponents log K*/log N, and",
        "           the difference between them -- what the operative",
        "           budget costs in the exponent -- with the",
        "           least-squares slope of that difference against",
        "           log N, its r.m.s. residual and standard error.",
        "NULL: none is run and none applies. Two crossings of the same",
        "      measured sum against two computed thresholds are located",
        "      and differenced; there is no background to detect",
        "      against. The sign controls for this field were run in",
        "      lab_direct_level.py and lab_split_budget.py.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < " + str(KCAP) + "; m odd,",
        "       squarefree and coprime to k, m <= (N-1)/k; the sieve",
        "       weight over the odd primes below " + str(QSIEVE) + ";",
        "       beta refitted as sum(H P)/sum(P^2) on the same range;",
        "       S(N) and A(N) from Euler products at the fixed bound "
        + str(CLIM) + ".",
        "       Every N is 2^a 5^b, one odd radical, as RADICALS says.",
        "       The published gaps are read from",
        "       results/audit_model_transfer.txt.",
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
