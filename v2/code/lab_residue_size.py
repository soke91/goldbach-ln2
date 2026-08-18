# -*- coding: utf-8 -*-
r"""
The size law of what is left after the elementary part is removed.

WHAT IS AT STAKE

Remark {#rem:predictable} split the dilated wall: with beta the
least-squares scale, H(N;k) = beta P(N;k) + R(N;k), where
P = sum_m mu(m) w(m,k) is an elementary sieve-weighted Mobius sum and
sum(log k)|beta P| is 0.86 to 0.96 of sum(log k)|H|.  The residue R
carries about half the demand and is the only part that is a
Mobius-prime correlation in any essential sense.

Everything the program still needs turns on its size.  Remark
{#rem:directlevel} showed that the measured level K*_H is exactly what
square-root cancellation in H would predict, and that the whole of
[eq:directcond] reduces to having it.  If R obeys |R| ~ (N/k)^{1/2}
then the square-root behaviour lives in the residue and the elementary
part is a separate, computable problem.  If R is larger, the split
does not help; if smaller, the elementary part is the whole
difficulty.

BACKS: Remark {#rem:residue} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  T1  The residue is centred: its mass-weighted fraction f+ of k with
      R > 0 lies inside the range spanned by 16 draws that keep |R|
      and randomise its signs, at every N.
  T2  The residue has square-root size: fitting the octave means of
      |R| against N/k gives an exponent in [0.40, 0.60] at every N.
  T3  The lean is in the elementary part: the mass-weighted f+ of P
      lies below the minimum of those same 16 draws at every N.
  T4  The exponent is robust: dropping the shortest octave, and
      dropping the longest, moves it by under 0.10.

REFUTATION RULE (fixed before the run)

  T1  REFUTED if f+ falls outside the draws' range at any N, which
      would mean the elementary part has not absorbed the lean.
  T2  REFUTED if the exponent leaves [0.40, 0.60] at any N. This is
      the one that decides whether the square-root behaviour is in the
      residue.
  T3  REFUTED if P's f+ is at or above the draws' minimum at any N.
  T4  REFUTED if either leave-one-out moves the exponent by 0.10 or
      more.

  T5 was added after T1-T4 had first run and is disclosed as such.
  The original binning ended in an open [32768, inf), which gate check
  G29 now forbids. Closing it exposes what the tail was hiding: on
  this same k-range the new top bins hold a single k at the largest N.
  lab_elementary_reach.py measured what that does -- an exponent of
  0.3674 where the well-populated bins give 0.5178 -- and fixed a
  population floor for it, so the same floor is applied here. A floor
  is a threshold, so it is swept.

  T5  The exponent does not depend on where the floor is put:
      sweeping it over 5, 10 and 20 k per octave moves the exponent by
      less than 0.05 at every N.

REFUTATION RULE for T5 (fixed before the sweep was run)

  T5  REFUTED if the exponent moves by 0.05 or more under the sweep,
      which would mean the fit is decided by the floor and not by the
      data.

  All five gate.

  THE CONTROL is the sign randomisation used by T1 and T3: every
  |R(N;k)| is held fixed and only the signs are redrawn, so the band
  it gives is the sampling range for f+ under no sign structure at
  all, computed rather than assumed.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_residue_size.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
# every bin is closed and every fitted bin is populated. The top
# bin used to be [32768, inf); gate check G29 forbids that now,
# and lab_elementary_reach.py measured why: an open tail folds
# the thinnest, furthest part of the range into one point with
# no abscissa, and on the same k-range it moved a fit of exactly
# this kind from 0.5178 to 0.3674. The largest N/k reached here
# is under 2097152, so the last edge closes the range.
OCT = [2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288,
       2097152]
MINPTS = 10                   # k per octave needed to fit it
DRAWS = 16
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
    """Refit dropping each end in turn, and report the spread."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    sp = max(f) - min(f)
    say("  leave-one-out on %s: full %.4f, without the shortest octave "
        "%.4f," % (name, f[0], f[1]))
    say("  without the longest %.4f -- spread %.4f" % (f[2], sp))
    say("SWEPT %s octave-range %.4f" % (name, sp))
    return f, sp


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

    res = []
    for N in NS:
        PN = factor_set(N)
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
            H = float((lam[vals] * g).sum())
            w = np.ones(ms.size, dtype=np.float64)
            for q in QS:
                if k % q == 0:
                    continue
                w *= np.where(vals % q == 0, 0.0, q / (q - 1.0))
            ks.append(k)
            Hs.append(H)
            Ps.append(float((g * w).sum()))
        ks = np.array(ks, dtype=np.int64)
        Hs = np.array(Hs)
        Ps = np.array(Ps)
        beta = float((Hs * Ps).sum() / (Ps * Ps).sum())
        R = Hs - beta * Ps
        res.append((N, ks, Hs, Ps, R, N // ks, beta))
        say("  N = %-10d  #k = %-7d beta = %.4f" % (N, ks.size, beta))

    say()
    say("T1/T3  is the lean in the elementary part or in the residue?")
    say("  N            f+ of R     f+ of P     draws min   draws max")
    t1 = t3 = True
    for i, (N, ks, Hs, Ps, R, Minner, beta) in enumerate(res):
        lw = np.log(ks.astype(np.float64))

        def frac(v):
            w = lw * np.abs(v)
            return float(w[v > 0].sum() / w.sum())

        fr, fp = frac(R), frac(Ps)
        rng = np.random.default_rng(SEED + i)
        aw = lw * np.abs(R)
        got = []
        for _ in range(DRAWS):
            e = rng.choice([-1.0, 1.0], size=R.size)
            got.append(float(aw[e > 0].sum() / aw.sum()))
        lo, hi = float(np.min(got)), float(np.max(got))
        if not (lo <= fr <= hi):
            t1 = False
        if fp >= lo:
            t3 = False
        say("  %-12d %-11.4f %-11.4f %-11.4f %.4f" % (N, fr, fp, lo, hi))
    say("  T1 the residue is centred      %s"
        % ("hold" if t1 else "REFUTED"))
    say("  T3 the elementary part leans   %s"
        % ("hold" if t3 else "REFUTED"))

    say()
    say("T2  the size law, octave means of |R| against N/k")
    say("  N/k octave   " + "  ".join("N=%-9d" % N for N in NS))
    prof = {N: [] for N in NS}
    cent = {N: [] for N in NS}
    cnts = {N: [] for N in NS}
    for a, b in zip(OCT, OCT[1:]):
        cells = []
        for N, ks, Hs, Ps, R, Minner, beta in res:
            sel = (Minner >= a) & (Minner < b)
            n = int(sel.sum())
            cnts[N].append(n)
            if n:
                prof[N].append(float(np.abs(R[sel]).mean()))
                # the bin centre is the MEAN of N/k inside it, not a
                # nominal geometric midpoint: the bins are wide and a
                # nominal centre sits where the k are not.
                cent[N].append(float(Minner[sel].mean()))
                cells.append("%9.2f/%-5d" % (prof[N][-1], n))
            else:
                prof[N].append(float("nan"))
                cent[N].append(float("nan"))
                cells.append("%9s/%-5d" % ("-", 0))
        say("  [%-6d,%-7s) %s" % (a, str(b), " ".join(cells)))

    def fitN(N, floor):
        y = np.array(prof[N])
        c = np.array(cent[N])
        ok = (~np.isnan(y) & (y > 0) & ~np.isnan(c)
              & (np.array(cnts[N]) >= floor))
        x = np.log(c[ok])
        yy = np.log(y[ok])
        return x, yy, float(np.polyfit(x, yy, 1)[0])

    t2 = True
    exps = []
    say("  N            fitted exponent   correlation   bins fitted")
    for j, N in enumerate(NS):
        x, yy, e = fitN(N, MINPTS)
        r = float(np.corrcoef(x, yy)[0, 1])
        exps.append((x, yy, e, r))
        if not (0.40 <= e <= 0.60):
            t2 = False
        say("  %-12d %-17.4f %-13.5f %d" % (N, e, r, x.size))
    say("  T2 %s" % ("hold" if t2 else "REFUTED"))

    say()
    say("T5  the population floor is a threshold, so it is swept")
    say("  N            floor 5     floor 10    floor 20    spread")
    t5 = True
    for N in NS:
        ee = [fitN(N, f)[2] for f in (5, 10, 20)]
        sp = max(ee) - min(ee)
        if sp >= 0.05:
            t5 = False
        say("  %-12d %-11.4f %-11.4f %-11.4f %.4f"
            % (N, ee[0], ee[1], ee[2], sp))
    say("  T5 %s" % ("hold" if t5 else "REFUTED"))

    say()
    say("T4  robustness of the exponent, and the thinnest bin each fit")
    say("  actually stands on, which gate check G30 reads. A fit is")
    say("  only as good as its emptiest octave, and that number was")
    say("  invisible here until it had gone wrong twice.")
    t4 = True
    for j, N in enumerate(NS):
        x, yy, _, r = exps[j]
        _, sp = loo(x, yy, "residue_size_N%d" % N, say)
        say("POP residue_size_N%d %d"
            % (N, min(c for c in cnts[N] if c >= MINPTS)))
        say("CORR residue_size_N%d %.5f" % (N, abs(r)))
        if sp >= 0.10:
            t4 = False
    say("  T4 %s" % ("hold" if t4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). The three parts side by side, octave")
    say("  means at the largest N:")
    say("  N/k octave   |H|         |beta P|    |R|         |R|/|H|")
    N, ks, Hs, Ps, R, Minner, beta = res[-1]
    for a, b in zip(OCT, OCT[1:]):
        sel = (Minner >= a) & (Minner < b)
        if not sel.sum():
            continue
        h = float(np.abs(Hs[sel]).mean())
        p = float(np.abs(beta * Ps[sel]).mean())
        r = float(np.abs(R[sel]).mean())
        say("  [%-6d,%-7s) %-11.2f %-11.2f %-11.2f %.4f"
            % (a, str(b), h, p, r, r / h))

    say()
    say("  Cross-check lines. lab_predictable_part.py fits the same beta")
    say("  on the same k-range.")
    for N, ks, Hs, Ps, R, Minner, beta in res:
        say("AGREE beta_HP N=%d %.6f 0.01" % (N, beta))

    say()
    say("=" * 70)
    say("T1 %s  T2 %s  T3 %s  T4 %s  T5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (t1, t2, t3, t4, t5)))
    ok = t1 and t2 and t3 and t4 and t5
    say("the square-root behaviour is in the residue and the lean is in "
        "the elementary part" if ok else "REFUTED")

    head = [
        "STATISTIC: with beta the least-squares scale through the origin,",
        "           R = H - beta P; the mass-weighted fraction f+ of k",
        "           with R > 0 and with P > 0, against 16 sign draws that",
        "           hold |R| fixed; the octave means of |R| against N/k",
        "           and the exponent fitted from them, with its",
        "           leave-one-out spread; and |H|, |beta P|, |R| side by",
        "           side.",
        "NULL: the sign randomisation of T1 and T3 -- every |R(N;k)| held",
        "      fixed, its sign redrawn, 16 draws. The band it gives is",
        "      the sampling range of f+ under no sign structure at all,",
        "      computed rather than assumed, which is what M2 asks for.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree, coprime",
        "       to N, 2 <= k < 30000, which covers K* at every N; m odd",
        "       squarefree, coprime to k, m <= (N-1)/k; the sieve weight",
        "       uses the odd primes up to 30; seed 20260808.",
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
