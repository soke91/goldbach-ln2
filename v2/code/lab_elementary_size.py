# -*- coding: utf-8 -*-
r"""
The size law of the half that binds.

WHAT IS AT STAKE

Remark {#rem:splitbudget} found that the elementary half beta P, not
the residue R, is what exhausts the Goldbach budget: removing R
entirely moves theta' by about 0.02, removing beta P by three times
that.  So the operative question is what bounds
sum_{k<K}(log k)|P(N;k)|, and the first thing to know is the size of
|P(N;k)| itself.

P is a Mobius sum with no primes in it: with w(m,k) the sieve weight
for N - mk to be prime,

    P(N;k) = sum_{m <= N/k, m odd squarefree, (m,k)=1} mu(m) w(m,k),

and w is a deterministic function of (m,k) once N is fixed.  If |P|
obeys (N/k)^{1/2} then the split has bought a constant and not an
exponent -- every part of H would be square-root and the difficulty
would be the same shape in each.  If |P| is larger, the elementary
half is a harder object than the residue and the program's target
moves again.

BACKS: Remark {#rem:elemsize} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  |P| is square-root: fitting the octave means of |P| against N/k
      gives an exponent in [0.40, 0.60] at every N.
  Y2  mu achieves exactly what a random sign pattern does, no better:
      replacing mu(m) by a coin on the odd squarefree m with w
      untouched, the coin's exponent is also in [0.40, 0.60] and mu's
      is within 0.05 of the coins' median at every N.
  Y3  The sieve weight inflates rather than reshapes: the ratio of
      mean|P| to mean|F|, with F the unweighted sum_m mu(m) over the
      same range, has spread (max/min) under 2 across the octaves.
  Y4  So the split buys a constant and not an exponent: |a - e| is
      under 0.10, with e the residue's exponent read from
      results/lab_residue_size.txt.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED if the exponent leaves [0.40, 0.60] at any N. A
      refutation would move the program's target again, so it is the
      one that matters.
  Y2  REFUTED if the coin leaves that band, or if mu differs from the
      coins' median by 0.05 or more, at any N. The second would mean
      mu does better or worse than random signs in this sum.
  Y3  REFUTED if the spread reaches 2.
  Y4  REFUTED if |a - e| reaches 0.10.

  Y5 was added after Y1-Y4 had first run and is disclosed as such.
  The original binning ended in an open [32768, inf), which gate check
  G29 now forbids and lab_elementary_reach.py explains: an open tail
  folds the thinnest, furthest part of the range into one point with
  no abscissa. Closing it exposed what the tail had been hiding --
  the new top bin holds a single k at the largest N, and fitting
  through it gave 0.3249 with correlation 0.82468. A population floor
  is therefore applied, and a floor is a threshold, so it is swept.

  Y5  The exponent does not depend on where the population floor is
      put: sweeping it over 5, 10 and 20 k per octave moves the
      exponent by less than 0.05 at every N.

REFUTATION RULE for Y5 (fixed before the sweep was run)

  Y5  REFUTED if the exponent moves by 0.05 or more under the sweep,
      which would mean the fit is decided by the floor and not by the
      data.

  All five gate.

  NULL: the coin, eps(m) = +-1 on the odd squarefree m, with the sieve
  weight w(m,k), the summation range and the k-range untouched, so the
  sign pattern on the inner variable is the only difference. Eight
  draws. |P| is a magnitude whose mean is bounded away from zero, so
  the statistic stays well conditioned under it -- the criterion of
  {#rem:weightgapnull}.
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
OUT = os.path.join(ROOT, "results", "lab_elementary_size.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
# every bin is closed. The top bin used to be [32768, inf), and
# lab_elementary_reach.py showed what an open tail does to a fit of
# this kind: it folds the thinnest, furthest-out part of the range
# into one point with no abscissa, and that point moved the same
# measurement from 0.5178 to 0.3674 there. The largest N/k reached
# here is under 2097152, so the last edge closes the range.
OCT = [2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152]
COINS = 8
SEED = 20260808
MINPTS = 10                   # k per octave needed to fit it


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
    return pr, lam, mu


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


def loo(x, y, name, say):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    f = [float(np.polyfit(x[t], y[t], 1)[0])
         for t in (slice(None), slice(1, None), slice(0, -1))]
    sp = max(f) - min(f)
    say("  leave-one-out on %s: full %.4f, without the shortest octave "
        "%.4f," % (name, f[0], f[1]))
    say("  without the longest %.4f -- spread %.4f" % (f[2], sp))
    say("SWEPT %s octave-range %.4f" % (name, sp))
    return sp


def fit(cent, prof, cnt=None, floor=0):
    """the octave fit, over bins that hold enough k to mean anything.

    Closing the top bin is not enough on its own: at the largest N the
    bin [524288, 2097152) holds one k, and fitting through it dragged
    the exponent to 0.3249 with correlation 0.82468. lab_elementary_
    reach.py measured the same failure on a longer lever and fixed a
    population floor for it, so the same floor is applied here.
    """
    c = np.array(cent, dtype=float)
    y = np.array(prof, dtype=float)
    ok = ~np.isnan(c) & ~np.isnan(y) & (y > 0)
    if cnt is not None:
        ok &= np.array(cnt) >= floor
    x = np.log(c[ok])
    yy = np.log(y[ok])
    return x, yy, float(np.polyfit(x, yy, 1)[0])


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]

    rng = np.random.default_rng(SEED)
    coins = []
    for _ in range(COINS):
        c = np.zeros(NMAX + 1, dtype=np.float64)
        sel = sqf.copy()
        sel[::2] = False
        c[sel] = rng.choice([-1.0, 1.0], size=int(sel.sum()))
        coins.append(c)
    say("  %d coins on the odd squarefree m, sieve weight kept" % COINS)

    res = []
    for N in NS:
        PN = factor_set(N)
        ks, Ps, Fs = [], [], []
        Cs = [[] for _ in range(COINS)]
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
            Ps.append(float((g * w).sum()))
            Fs.append(float(g.sum()))
            for j in range(COINS):
                Cs[j].append(float((coins[j][ms] * w).sum()))
        ks = np.array(ks, dtype=np.int64)
        res.append((N, ks, np.array(Ps), np.array(Fs),
                    [np.array(c) for c in Cs], N // ks))
        say("  N = %-10d  #k = %d" % (N, ks.size))

    def profile(Minner, v):
        cent, prof, cnt = [], [], []
        for a, b in zip(OCT, OCT[1:]):
            sel = (Minner >= a) & (Minner < b)
            cnt.append(int(sel.sum()))
            if sel.sum():
                cent.append(float(Minner[sel].mean()))
                prof.append(float(np.abs(v[sel]).mean()))
            else:
                cent.append(float("nan"))
                prof.append(float("nan"))
        return cent, prof, cnt

    say()
    say("Y1  octave means of |P| against N/k, with the count of k in")
    say("  each bin; bins holding fewer than %d are not fitted" % MINPTS)
    say("  N/k octave   " + "  ".join("N=%-13d" % N for N in NS))
    tab = []
    for i, (a, b) in enumerate(zip(OCT, OCT[1:])):
        row, cells = [], []
        for N, ks, Ps, Fs, Cs, Minner in res:
            sel = (Minner >= a) & (Minner < b)
            n = int(sel.sum())
            v = float(np.abs(Ps[sel]).mean()) if n else float("nan")
            row.append(v)
            cells.append("%9.3f/%-5d" % (v, n) if n else
                         "%9s/%-5d" % ("-", 0))
        tab.append(row)
        say("  [%-6d,%-7s) %s" % (a, str(b), " ".join(cells)))
    y1 = True
    exps = []
    say("  N            exponent   correlation   bins fitted")
    for j, (N, ks, Ps, Fs, Cs, Minner) in enumerate(res):
        cent, prof, cnt = profile(Minner, Ps)
        x, yy, e = fit(cent, prof, cnt, MINPTS)
        r = float(np.corrcoef(x, yy)[0, 1])
        exps.append((x, yy, e, cnt, r))
        if not (0.40 <= e <= 0.60):
            y1 = False
        say("  %-12d %-10.4f %-13.5f %d" % (N, e, r, x.size))
    say("  Y1 %s" % ("hold" if y1 else "REFUTED"))

    say()
    say("Y2  the control: mu against coins on the same sum")
    say("  N            mu        coin min   coin median   coin max"
        "   |mu - median|")
    y2 = True
    for j, (N, ks, Ps, Fs, Cs, Minner) in enumerate(res):
        ce = []
        for c in Cs:
            cent, prof, cnt = profile(Minner, c)
            ce.append(fit(cent, prof, cnt, MINPTS)[2])
        med = float(np.median(ce))
        if not (0.40 <= min(ce) and max(ce) <= 0.60):
            y2 = False
        if abs(exps[j][2] - med) >= 0.05:
            y2 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %-10.4f %.4f"
            % (N, exps[j][2], min(ce), med, max(ce),
               abs(exps[j][2] - med)))
    say("  Y2 %s" % ("hold" if y2 else "REFUTED"))

    say()
    say("Y3  what the sieve weight does to the size")
    say("  N            mean|P|/mean|F| by octave        spread")
    y3 = True
    for N, ks, Ps, Fs, Cs, Minner in res:
        rr = []
        for a, b in zip(OCT, OCT[1:]):
            sel = (Minner >= a) & (Minner < b)
            if sel.sum() >= MINPTS and np.abs(Fs[sel]).mean() > 0:
                rr.append(float(np.abs(Ps[sel]).mean()
                                / np.abs(Fs[sel]).mean()))
        sp = max(rr) / min(rr)
        if sp >= 2.0:
            y3 = False
        say("  %-12d %-32s %.4f"
            % (N, " ".join("%.2f" % v for v in rr), sp))
    say("  Y3 %s" % ("hold" if y3 else "REFUTED"))

    say()
    rp = os.path.join(ROOT, "results", "lab_residue_size.txt")
    e_res = None
    if os.path.exists(rp):
        # exponent, correlation, and since the population floor was
        # added there is a count of fitted bins after them
        m = re.search(r"(?m)^\s*3200000\s+([\d.]+)\s+[\d.]+\s+\d+\s*$",
                      io.open(rp, encoding="utf-8").read())
        if m:
            e_res = float(m.group(1))
    if e_res is None:
        say("Y4  results/lab_residue_size.txt gives no exponent to read")
        y4 = False
    else:
        d = abs(exps[-1][2] - e_res)
        y4 = d < 0.10
        say("Y4  |P| exponent %.4f against the residue's %.4f, read from"
            % (exps[-1][2], e_res))
        say("    results/lab_residue_size.txt: gap %.4f  (cap 0.10)   %s"
            % (d, "hold" if y4 else "REFUTED"))

    say()
    say("Y5  the population floor is a threshold, so it is swept")
    say("  N            floor 5     floor 10    floor 20    spread")
    y5 = True
    for j, (N, ks, Ps, Fs, Cs, Minner) in enumerate(res):
        cent, prof, cnt = profile(Minner, Ps)
        ee = [fit(cent, prof, cnt, f)[2] for f in (5, 10, 20)]
        sp = max(ee) - min(ee)
        if sp >= 0.05:
            y5 = False
        say("  %-12d %-11.4f %-11.4f %-11.4f %.4f"
            % (N, ee[0], ee[1], ee[2], sp))
    say("  Y5 %s" % ("hold" if y5 else "REFUTED"))

    say()
    say("  The thinnest bin each fit actually stands on, which gate")
    say("  check G30 reads. A fit is only as good as its emptiest")
    say("  octave, and that number was invisible until it had gone")
    say("  wrong twice.")
    for j, N in enumerate(NS):
        x, yy, _, cnt, r = exps[j]
        loo(x, yy, "elementary_size_N%d" % N, say)
        say("POP elementary_size_N%d %d"
            % (N, min(c for c in cnt if c >= MINPTS)))
        say("CORR elementary_size_N%d %.5f" % (N, abs(r)))

    say()
    say("  DIAGNOSTIC (post hoc). The three parts of H side by side at")
    say("  the largest N, as octave means:")
    say("  N/k octave   |P|         |F|         |P|/|F|")
    N, ks, Ps, Fs, Cs, Minner = res[-1]
    for a, b in zip(OCT, OCT[1:]):
        sel = (Minner >= a) & (Minner < b)
        if not sel.sum():
            continue
        p = float(np.abs(Ps[sel]).mean())
        f = float(np.abs(Fs[sel]).mean())
        say("  [%-6d,%-7s) %-11.3f %-11.3f %.4f"
            % (a, str(b), p, f,
               p / f if f else float("nan")))

    say()
    say("=" * 70)
    say("Y1 %s  Y2 %s  Y3 %s  Y4 %s  Y5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (y1, y2, y3, y4, y5)))
    ok = y1 and y2 and y3 and y4 and y5
    say("the elementary half is square-root too, so the split buys a "
        "constant and not an exponent" if ok else "REFUTED")

    head = [
        "STATISTIC: octave means of |P(N;k)| against the length N/k of",
        "           the inner sum and the exponent fitted from them, for",
        "           mu and for 8 coins; the ratio of mean|P| to mean|F|",
        "           with F the unweighted Mobius sum over the same range;",
        "           and the gap to the residue's exponent read from",
        "           results/lab_residue_size.txt.",
        "NULL: the coin, eps(m) = +-1 on the odd squarefree m, with the",
        "      sieve weight w(m,k), the summation range and the k-range",
        "      untouched, eight draws. |P| is a magnitude whose mean is",
        "      bounded away from zero, so it stays well conditioned under",
        "      the control -- the criterion of [rem:weightgapnull].",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree, coprime",
        "       to N, 2 <= k < 30000; m odd squarefree, coprime to k,",
        "       m <= (N-1)/k; the sieve weight uses the odd primes up to",
        "       30; seed 20260808.",
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
