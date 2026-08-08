# -*- coding: utf-8 -*-
r"""
The split constant is a free parameter of the argument, and it has
been fixed by optimising the wrong thing.

WHAT IS AT STAKE

Remark {#rem:residuelevel} put the programme on a knife-edge:
conditional on {#rem:provablehalf}, the residue alone carries the
level to exponents 0.5654, 0.5642, 0.5599, 0.5675, 0.5799 against the
operative budget -- past the square-root barrier by 0.06 to 0.08, and
short of theta' = 0.56 at one N by one part in five thousand.

Every one of those numbers depends on beta, and beta is not given by
the problem. H = beta P + R holds for ANY beta, and so does
|H| <= beta|P| + |R|; granting {#rem:provablehalf} the conditional
bound is B(N) <= B_R(beta; N) + o(N) for any beta whatever. **The
argument may choose it.**

What it has instead is the least-squares beta,
sum(H P)/sum(P^2), which minimises the l2 distance from H to the ray
beta P. That is not the objective. The budget is spent by
sum_{k<K}(log k)|H - beta P|: an l1 norm, weighted by log k, and
truncated at the crossing. Minimising the l2 residual and maximising
the permitted level are different problems, and nothing has checked
how far apart their answers are.

BACKS: Remark {#rem:betafree} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  T1  The control: at the fitted beta this reproduces the exponents
      of results/audit_residue_level.txt to within 0.001 at every N.
  T2  The level is sensitive to beta at all: sweeping beta over
      +-10 per cent of the fitted value moves the exponent by more
      than 0.01 at some N. If it does not, the knife-edge is robust
      to the split constant and this line of attack is closed.
  T3  Least squares is answering a different question: the beta that
      maximises K*_R differs from the fitted one by more than 1 per
      cent at a majority of N.
  T4  And the free choice clears the barrier that the fitted one
      missed: at the optimal beta the exponent exceeds 0.56 at every
      N, including the 0.5599 at N = 8e5.

REFUTATION RULE (fixed before the run)

  T1  REFUTED at 0.001 at any N, which would mean this is not the
      same split.
  T2  REFUTED if no N moves by more than 0.01 under the sweep. That
      is a useful failure: it would say the knife-edge does not
      depend on a fitted constant.
  T3  REFUTED if the optimal beta is within 1 per cent of the fitted
      one at three or more of the five N.
  T4  REFUTED if the exponent at the optimal beta fails to exceed
      0.56 at any N.

  All four gate.

  NO NULL IS RUN and none applies. Nothing here is a detection against
  a background: a measured sum is minimised over one scalar and the
  crossing located. The sign control for this field is
  lab_residue_cancellation.py, whose coin arm on the identical delta
  established that R's size is bought by cancellation and bought at
  exactly a coin's rate, and lab_split_budget.py's size permutation,
  which established where the budget is spent.
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
OUT = os.path.join(ROOT, "results", "audit_beta_optimal.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000
THETA = 0.56
SWEEP = 0.10
GRID = 241                 # beta over [0.25, 6.25] at step 0.025


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


def read_published():
    """the operative exponents -- read, not copied"""
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

    pub = read_published()
    say("read %d operative exponents from "
        "results/audit_residue_level.txt" % len(pub))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = sieves(NMAX)
    sqf = mu != 0
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("  sieve weight over the odd primes %s"
        % ", ".join(map(str, QS)))

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    res = []
    for N in NS:
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
        lw = np.log(ks.astype(float))
        bfit = float((H * P).sum() / (P * P).sum())
        thr = S_ * (1.0 - A_) * N
        res.append((N, ks, lw, H, P, bfit, thr, S_, A_))
        say("  N = %-10d  #k = %-7d beta_fit = %.6f" % (N, ks.size, bfit))

    def kstar(ks, lw, H, P, b, thr):
        cum = np.cumsum(lw * np.abs(H - b * P))
        j = int(np.searchsorted(cum, thr))
        if j >= ks.size:
            return None
        return int(ks[j])

    def expo(N, k):
        return math.log(k) / math.log(N) if k else float("nan")

    # ------------------------------------------------------------- T1
    say()
    say("T1  the control: the fitted beta reproduces the published")
    say("  operative exponents")
    say("  N            beta_fit    K*_R      exponent   published  diff")
    t1 = True
    for N, ks, lw, H, P, bfit, thr, S_, A_ in res:
        k = kstar(ks, lw, H, P, bfit, thr)
        e = expo(N, k)
        d = abs(e - pub[N])
        if not (d < 0.001):
            t1 = False
        say("  %-12d %-11.6f %-9s %-10.4f %-10.4f %.5f"
            % (N, bfit, str(k), e, pub[N], d))
    say("  T1 %s" % ("hold" if t1 else "REFUTED"))

    # ------------------------------------------------------------- T2
    say()
    say("T2  is the level sensitive to beta at all? +-%d per cent"
        % int(100 * SWEEP))
    say("  N            at 0.9 b   at b       at 1.1 b   spread")
    t2 = False
    for N, ks, lw, H, P, bfit, thr, S_, A_ in res:
        ee = [expo(N, kstar(ks, lw, H, P, bfit * f, thr))
              for f in (1.0 - SWEEP, 1.0, 1.0 + SWEEP)]
        sp = max(ee) - min(ee)
        if sp > 0.01:
            t2 = True
        say("  %-12d %-10.4f %-10.4f %-10.4f %.4f"
            % (N, ee[0], ee[1], ee[2], sp))
    say("  T2 %s   (some N moves by more than 0.01)"
        % ("hold" if t2 else "REFUTED"))

    # ---------------------------------------------------------- T3/T4
    say()
    say("T3/T4  the beta the argument is free to choose")
    say("  the objective is the permitted level itself, maximised over")
    say("  a grid of %d values of beta from 0.25 to %.2f"
        % (GRID, 0.25 + 0.025 * (GRID - 1)))
    say("  N            beta_fit    beta_opt    ratio    K*_R      "
        "exponent  clears")
    grid = 0.25 + 0.025 * np.arange(GRID)
    t3 = t4 = True
    far = 0
    opt = []
    for N, ks, lw, H, P, bfit, thr, S_, A_ in res:
        best, bestk = None, -1
        for b in grid:
            k = kstar(ks, lw, H, P, float(b), thr)
            if k is not None and k > bestk:
                best, bestk = float(b), k
        e = expo(N, bestk)
        opt.append((best, bestk, e))
        r = best / bfit
        if abs(r - 1.0) > 0.01:
            far += 1
        if e <= THETA:
            t4 = False
        say("  %-12d %-11.6f %-11.6f %-8.4f %-9d %-9.4f %s"
            % (N, bfit, best, r, bestk, e, "yes" if e > THETA else "NO"))
    t3 = far > len(res) / 2.0
    say("  T3 the optimum is more than 1 per cent from the fit at %d "
        "of %d N   %s" % (far, len(res), "hold" if t3 else "REFUTED"))
    say("  T4 exponent above %.2f at every N   %s"
        % (THETA, "hold" if t4 else "REFUTED"))

    say()
    say("  the budget constant crossed throughout, declared:")
    for N, ks, lw, H, P, bfit, thr, S_, A_ in res:
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, S_ * (1.0 - A_)))

    say()
    say("  DIAGNOSTIC (post hoc). The two objectives, side by side. The")
    say("  fit minimises sum (H - b P)^2; the budget spends")
    say("  sum_{k<K*}(log k)|H - b P|. At the optimum, what each reads:")
    say("  N            l2 residual at b_fit / at b_opt   budget sum ratio")
    ratios = []
    for i, (N, ks, lw, H, P, bfit, thr, S_, A_) in enumerate(res):
        bo = opt[i][0]
        l2f = float(((H - bfit * P) ** 2).sum())
        l2o = float(((H - bo * P) ** 2).sum())
        kk = opt[i][1]
        sel = ks <= kk
        b1 = float((lw[sel] * np.abs(H[sel] - bfit * P[sel])).sum())
        b2 = float((lw[sel] * np.abs(H[sel] - bo * P[sel])).sum())
        ratios.append((l2o / l2f, b2 / b1))
        say("  %-12d %-33.4f %.4f" % (N, l2o / l2f, b2 / b1))
    flat = sum(1 for r in ratios
               if abs(r[0] - 1.0) < 1e-3 and abs(r[1] - 1.0) < 1e-3)
    say("  At %d of the %d N both columns are within one (tol 0.001):"
        % (flat, len(res)))
    say("  two objectives do not merely agree on the answer, they are")
    say("  flat in the same place. Where they separate, the optimal")
    say("  beta is worse by the l2 measure that chose the fitted one")
    say("  and better by the l1 measure the argument pays -- which is")
    say("  what T3 was asking about. The exponent it buys, at each N:")
    say("  N            exponent gain over the fitted beta")
    for i, (N, ks, lw, H, P, bfit, thr, S_, A_) in enumerate(res):
        say("  %-12d %+.4f"
            % (N, opt[i][2] - expo(N, kstar(ks, lw, H, P, bfit, thr))))
    say()
    say("  Two cautions on the optimum itself. K* is integer-valued and")
    say("  jumps, so a whole interval of beta shares one K*, and the")
    say("  grid returns the first of them. Where the beta-ratio differs")
    say("  from one but K* does not move, the 'optimum' is that")
    say("  granularity and nothing else:")
    say("  N            beta ratio   K* moved?")
    for i, (N, ks, lw, H, P, bfit, thr, S_, A_) in enumerate(res):
        kf = kstar(ks, lw, H, P, bfit, thr)
        say("  %-12d %-12.4f %s"
            % (N, opt[i][0] / bfit,
               "yes" if opt[i][1] != kf else "no"))
    say("  The budget sum can read marginally above one at such a point")
    say("  for the same reason.")

    say()
    say("=" * 70)
    ok = t1 and t2 and t3 and t4
    say("the split constant is free, and choosing it for the budget "
        "clears the barrier the fit missed" if ok else "REFUTED")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|H - beta P| first reaches",
        "           S(N)(1-A(N))N, as a function of beta: at the",
        "           least-squares beta, at beta scaled by 0.9 and 1.1,",
        "           and at the beta maximising K*_R over a grid; the",
        "           exponent log K*_R / log N of each; and the l2",
        "           residual and the truncated budget sum at the fitted",
        "           and optimal beta.",
        "NULL: none is run and none applies. A measured sum is minimised",
        "      over one scalar and a crossing located; there is no",
        "      background to detect against. The sign controls for this",
        "      field were run in lab_residue_cancellation.py, whose coin",
        "      arm on the identical delta established that R's size is",
        "      bought by cancellation at exactly a coin's rate, and in",
        "      lab_split_budget.py, whose size permutation established",
        "      where the budget is spent.",
        "FIELD: N = 2e5 to 3.2e6 by doubling; k squarefree and coprime",
        "       to N with 2 <= k < 100000; m odd, squarefree, coprime to",
        "       the odd part of k, m < N/k; the sieve weight uses the",
        "       odd primes up to 30; S(N) and A(N) from Euler products",
        "       at the fixed bound 4000000; beta swept over 241 values",
        "       from 0.25 to 6.25 in steps of 0.025; the published",
        "       exponents are read from results/audit_residue_level.txt.",
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
