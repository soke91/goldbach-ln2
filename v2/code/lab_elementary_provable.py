# -*- coding: utf-8 -*-
r"""
Which half of the wall has an unconditional bound in the right shape,
and where that bound starts to bite.

WHAT IS AT STAKE

Remark {#rem:splitbudget} put the target on the elementary half: at
the operative truncation beta B_P takes 0.81 to 0.87 of the Goldbach
budget and B_R takes about half. That is a statement about SIZE at
accessible N. It is not a statement about what can be PROVED, and the
two point in opposite directions.

P has no primes in it. Written out,

    P(N;k) = C_k * sum_{m in S} mu(m),
    S = { m < N/k : m odd, squarefree, (m,k) = 1,
                    m != N k^{-1} (mod q) for every odd prime q <= 29 }

every condition defining S is either multiplicative -- odd, squarefree,
coprime to k -- or a residue condition to a BOUNDED modulus. Coprimality
to k is not a residue condition to modulus k: in the Dirichlet series it
deletes Euler factors, contributing prod_{p|k}(1-1/p)^{-1} and nothing
worse. So the classical unconditional estimate

    |sum_{m <= x} mu(m)| <= A x exp(-c sqrt(log x))

is of the right shape for P uniformly in k, where for R -- an honest
Mobius-prime correlation of length N/k at level k -- no such bound is
available at all past k = N^{1/2}.

If that is right the target moves. sum_{k<K}(log k)|P| would be o(N)
for every fixed theta' < 1, the elementary half would be asymptotically
free, and the whole obstruction would sit in the residue that
{#rem:splitbudget} measured as the SMALLER half.

This script cannot prove the uniformity. It can do the two things that
decide whether the idea is worth pursuing: falsify the bound's shape
against the measurement if it is wrong, and, if it survives, compute
where the bound is actually strong enough to pay the budget.

BACKS: Remark {#rem:provablehalf} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The shape is not violated: with c at the classical 0.2098, the
      ratio |P(N;k)| / [ (N/k) exp(-c sqrt(log(N/k))) L(k) ] stays
      below 1 at every (N,k) measured, so no implied constant above 1
      is forced. L(k) = prod_{p|k} (1-1/p)^{-1}.
  W2  And the bound is loose by a growing factor, because the truth is
      square-root and the bound is nearly linear: that ratio falls
      monotonically across the octaves of N/k.
  W3  At accessible N the bound is useless: the bound's own budget
      share sum_{k<K}(log k)*bound / [S(N)(1-A(N))N] exceeds 10 at
      every N measured, with K = N^0.56.
  W4  It does eventually pay, and not anywhere near a computation:
      solving for the N at which the bound equals the budget gives
      log10 N above 1000.

REFUTATION RULE (fixed before the run)

  W1  REFUTED if the ratio reaches 1 at any (N,k). That is the
      interesting failure: it would mean the classical shape needs an
      implied constant above 1 here, and the constant would then be
      measured rather than assumed.
  W2  REFUTED if the octave profile of that ratio is not monotone
      decreasing at a majority of N.
  W3  REFUTED if the share falls to 10 or below at any N, which would
      mean the unconditional bound already has numerical content.
  W4  REFUTED if log10 N is at or below 1000, which would put the
      reduction within reach of an explicit computation.

  All four gate.

  NULL: a coin arm on the identical sifted set -- 16 global sign
  vectors eps(m) = +-1 on the odd squarefree m, each held fixed across
  all k exactly as mu is, same mask, same octaves -- carried through
  the same comparison against the same bound. It answers the question
  the ratio alone cannot: how much of the bound's looseness is the
  shape, (N/k) against (N/k)^{1/2}, and how much is mu. A coin is
  square-root too, so if the looseness is shape the coin sits the same
  distance from the bound as mu up to a constant; if mu were doing
  something the bound cannot see, it would not. The global vector is
  the convention lab_mu_vs_coin_size.py fixed, and using it here makes
  the two measurements of that constant comparable.
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
OUT = os.path.join(ROOT, "results", "lab_elementary_provable.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
OCT = [2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152]
MINPTS = 10
CZERO = 0.2098          # the classical exponent in exp(-c sqrt(log x))
CSWEEP = (0.15, 0.2098, 0.30)
COINS = 16
SEED = 20260808
THETA = 0.56
CLIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def mobius(n):
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
    return mu


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


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    mu = mobius(NMAX)
    sqf = mu != 0
    oddsqf = np.zeros(NMAX + 1, dtype=bool)
    oddsqf[1::2] = True
    oddsqf &= sqf
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("  sieve weight over the odd primes %s"
        % ", ".join(map(str, QS)))

    # the null uses GLOBAL sign vectors, one eps per m held across all
    # k exactly as mu is, which is the convention lab_mu_vs_coin_size.py
    # fixed and argued for: a fresh sign per (N,k) destroys the across-k
    # correlation mu has and is not the same control.
    rng = np.random.default_rng(SEED)
    coinmat = np.zeros((COINS, NMAX + 1), dtype=np.int8)
    _idx = np.flatnonzero(oddsqf)
    for _j in range(COINS):
        coinmat[_j, _idx] = rng.integers(0, 2, size=_idx.size) * 2 - 1
    del _idx
    say("  %d global sign vectors for the null" % COINS)
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
        thr = S_ * (1.0 - A_) * N

        ks, Pv, Lv, Cv = [], [], [], []
        for k in range(2, KCAP):
            if not sqf[k] or any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
            if M < 2:
                continue
            mask = oddsqf[1:M + 1].copy()
            fk = factor_set(k)
            for p in fk:
                if p > 2:
                    mask[p - 1::p] = False
            if not mask.any():
                continue
            C = 1.0
            for q in QS:
                if k % q == 0:
                    continue
                C *= q / (q - 1.0)
                r = (N * pow(k, -1, q)) % q
                mask[(r - 1) % q::q] = False
            ks.append(k)
            Pv.append(C * abs(int(mu[1:M + 1][mask]
                                  .sum(dtype=np.int64))))
            Lv.append(math.prod(1.0 / (1.0 - 1.0 / p) for p in fk))
            Cv.append(C * np.abs(coinmat[:, 1:M + 1][:, mask]
                                 .sum(axis=1, dtype=np.int64)))
        ks = np.array(ks, dtype=np.int64)
        res.append((N, ks, np.array(Pv), np.array(Lv), N // ks,
                    S_, A_, thr, np.array(Cv)))
        say("  N = %-10d  #k = %-7d threshold S(1-A)N = %.4e"
            % (N, ks.size, thr))

    def bound(inner, L, c):
        """A x exp(-c sqrt(log x)) L, with A = 1"""
        x = np.asarray(inner, dtype=float)
        return x * np.exp(-c * np.sqrt(np.log(x))) * L

    # ------------------------------------------------------------- W1
    say()
    say("W1  is the classical shape violated anywhere?")
    say("  ratio = |P| / [ (N/k) exp(-%.4f sqrt(log(N/k))) L(k) ]"
        % CZERO)
    say("  N            max ratio    at N/k       #k")
    w1 = True
    for N, ks, Pv, Lv, inner, S_, A_, thr, Cv in res:
        b = bound(inner, Lv, CZERO)
        r = Pv / b
        j = int(np.argmax(r))
        if r[j] >= 1.0:
            w1 = False
        say("  %-12d %-12.4f %-12d %d"
            % (N, float(r[j]), int(inner[j]), ks.size))
    say("  W1 %s" % ("hold" if w1 else "REFUTED"))
    say()
    say("  Where the violation sits, disclosed after the run rather")
    say("  than folded into W1. The maximum is attained at the SHORTEST")
    say("  inner sum available at each N -- 7, 13, 30, 63, 108, which is")
    say("  N over the k-cap -- and an estimate of the form")
    say("  A x exp(-c sqrt(log x)) has no asymptotic content at x = 7.")
    say("  Restricting the comparison to longer sums:")
    say("  cut-off x0   max ratio over every (N,k) with N/k >= x0")
    AFULL = max(float((Pv / bound(inner, Lv, CZERO)).max())
                for N, ks, Pv, Lv, inner, S_, A_, thr, Cv in res)
    for x0 in (2, 8, 32, 128):
        mx = 0.0
        for N, ks, Pv, Lv, inner, S_, A_, thr, Cv in res:
            sel = inner >= x0
            if not sel.any():
                continue
            mx = max(mx, float((Pv[sel] / bound(inner[sel], Lv[sel],
                                                CZERO)).max()))
        say("  %-12d %.4f" % (x0, mx))
    say("  So the shape survives everywhere an asymptotic estimate")
    say("  could speak, and the constant W4 uses below is the")
    say("  unrestricted one, %.4f, which is the conservative choice."
        % AFULL)

    # ------------------------------------------------------------- W2
    say()
    say("W2  how loose, by octave of the inner length")
    say("  N/k octave        " + " ".join("N=%-9d" % N for N in NS))
    prof = {}
    pops = {}
    for a, b_ in zip(OCT, OCT[1:]):
        row = []
        for i, (N, ks, Pv, Lv, inner, S_, A_, thr, Cv) in enumerate(res):
            sel = (inner >= a) & (inner < b_)
            n = int(sel.sum())
            if n < MINPTS:
                row.append("%11s" % "-")
                continue
            v = float((Pv[sel] / bound(inner[sel], Lv[sel],
                                       CZERO)).mean())
            prof.setdefault(i, []).append((a, v))
            pops.setdefault(i, []).append(n)
            row.append("%11.5f" % v)
        say("  %-17d %s" % (a, " ".join(row)))
    mono = 0
    for i in prof:
        v = [t[1] for t in prof[i]]
        mono += int(all(v[j] > v[j + 1] for j in range(len(v) - 1)))
    w2 = mono > len(res) / 2.0
    say("  monotone decreasing at %d of %d N   %s"
        % (mono, len(res), "hold" if w2 else "REFUTED"))
    say("  W2 %s" % ("hold" if w2 else "REFUTED"))
    for i, N in enumerate(NS):
        say("POP provable_looseness_N%d %d" % (N, min(pops[i])))

    # ------------------------------------------------------------- W3
    say()
    say("W3  what the bound alone would spend at accessible N")
    say("  the sum runs over the admissible k < N^%.2f that this sweep"
        % THETA)
    say("  reaches; where N^%.2f exceeds the k-cap the sum is truncated"
        % THETA)
    say("  there and the share is a LOWER bound on what it would spend.")
    say("  N            K = N^%.2f   share of the budget   truncated"
        % THETA)
    w3 = True
    for N, ks, Pv, Lv, inner, S_, A_, thr, Cv in res:
        K = N ** THETA
        sel = ks < K
        s = float((np.log(ks[sel].astype(float))
                   * bound(inner[sel], Lv[sel], CZERO)).sum())
        share = s / thr
        if share <= 10.0:
            w3 = False
        say("  %-12d %-12.0f %-21.2f %s"
            % (N, K, share, "yes" if K > KCAP else "no"))
    say("  W3 %s   (floor 10)" % ("hold" if w3 else "REFUTED"))

    # ------------------------------------------------------------- W4
    say()
    say("W4  where the bound would pay the budget")
    say("  Asymptotically the sum over admissible k is")
    say("    B_bound(N) = A dL N int_0^{theta' u} v exp(-c sqrt(u-v)) dv")
    say("  with u = log N and dL the mean of L(k)(log k)/k over the")
    say("  admissible k against the same integral without L and without")
    say("  the exponential, both measured here rather than assumed.")
    N, ks, Pv, Lv, inner, S_, A_, thr, Cv = res[-1]
    w = np.log(ks.astype(float)) / ks.astype(float)
    dL = float((w * Lv).sum() / (math.log(float(ks.max())) ** 2 / 2.0))
    Amax = AFULL
    thrc = S_ * (1.0 - A_)
    say("  measured density-times-L factor dL = %.4f" % dL)
    say("  implied constant A forced by the data = %.4f" % Amax)
    say("  budget constant S(N)(1-A(N)) = %.6f" % thrc)

    def integral(u, c):
        v = np.linspace(0.0, THETA * u, 20001)
        return float(np.trapezoid(v * np.exp(-c * np.sqrt(u - v)), v))

    def solve(c):
        lo, hi = 10.0, 1e7
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if Amax * dL * integral(mid, c) > thrc:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi) / math.log(10.0)

    say("  c          log10 N at which the bound pays")
    got = []
    for c in CSWEEP:
        g = solve(c)
        got.append(g)
        say("  %-10.4f %.1f" % (c, g))
    base = solve(CZERO)
    w4 = base > 1000.0
    say("  W4 log10 N = %.1f at the classical c   (floor 1000)   %s"
        % (base, "hold" if w4 else "REFUTED"))
    say("BRACKET log10_N_elementary_provable %.4f %.4f %.4f"
        % (base, min(got), max(got)))
    say()
    say("  And the drift of the constants that bracket does NOT sweep,")
    say("  which gate check G33 reads. The sweep above moves c, the")
    say("  analytic exponent, which is a literature constant and not a")
    say("  measurement. The forecast also rests on two things this")
    say("  script measures -- the density-times-L factor dL and the")
    say("  implied constant A -- and those have drifts of their own:")
    say("  N            dL         A forced by the data")
    dLs, Amaxs = [], []
    for NN, kk, PP, LL, inn, SS_, AA_, tt, CC in res:
        ww = np.log(kk.astype(float)) / kk.astype(float)
        d1 = float((ww * LL).sum()
                   / (math.log(float(kk.max())) ** 2 / 2.0))
        a1 = float((PP / bound(inn, LL, CZERO)).max())
        dLs.append(d1)
        Amaxs.append(a1)
        say("  %-12d %-10.4f %.4f" % (NN, d1, a1))
    ddL = (max(dLs) - min(dLs)) / float(np.mean(dLs))
    dA = (max(Amaxs) - min(Amaxs)) / float(np.mean(Amaxs))
    say("  relative spread: dL %.4f, A %.4f" % (ddL, dA))
    say("DRIFT elementary_provable_dL %.4f" % ddL)
    say("  dL does not drift at all: every N in this sweep has the same")
    say("  odd radical, so the admissible k-set and therefore dL are")
    say("  identical across it. That is a fact about the sweep, not a")
    say("  general one, and a sweep over different radicals would have")
    say("  to remeasure it.")
    say()
    say("  A is a different animal and is NOT declared as a drift. It")
    say("  falls monotonically with N because the maximum of the ratio")
    say("  is always attained at the shortest inner sum available,")
    say("  N over the k-cap, and that grows with N. It is a function of")
    say("  where one looks, not a wobbling constant, and the forecast")
    say("  uses the LARGEST value over all N, which is the conservative")
    say("  choice for an upper bound. What a factor of %.4f in A dL"
        % (1.0 + dA))
    say("  would move, either way:")
    say("  A dL scaled by        log10 N")
    for f2 in (1.0 / (1.0 + dA), 1.0, 1.0 + dA):
        lo2, hi2 = 10.0, 1e7
        for _ in range(200):
            mid = 0.5 * (lo2 + hi2)
            if f2 * Amax * dL * integral(mid, CZERO) > thrc:
                lo2 = mid
            else:
                hi2 = mid
        say("  %-21.4f %.1f" % (f2, 0.5 * (lo2 + hi2) / math.log(10.0)))
    say("  Compare the c-sweep's %.1f to %.1f. Both matter, and the"
        % (min(got), max(got)))
    say("  published bracket sweeps only c, so it is a bracket over the")
    say("  analytic exponent and not over everything the forecast rests")
    say("  on. Nothing here is close enough to a computable range for")
    say("  the difference to change a conclusion.")

    say()
    say("  THE NULL. Coins on the identical sifted set, measured")
    say("  against the identical bound. If the looseness is the shape")
    say("  -- (N/k) against (N/k)^{1/2} -- a coin sits the same")
    say("  distance from the bound as mu up to a constant, and that")
    say("  constant is the only thing mu contributes.")
    say("  The last column is the raw ratio of octave means, which is")
    say("  the statistic lab_mu_vs_coin_size.py reports; the two before")
    say("  it normalise by the bound first, which re-weights the k")
    say("  inside an octave by a factor of four in N/k and by L(k).")
    say("  N/k octave        mu / bound   coin / bound   normalised   raw")
    for a, v in prof[len(NS) - 1]:
        NN, kk, PP, LL, inn, SS_, AA_, tt, CC = res[-1]
        sel = (inn >= a) & (inn < a * 4)
        if sel.sum() < MINPTS:
            continue
        bb = bound(inn[sel], LL[sel], CZERO)
        cm = float(np.median((CC[sel] / bb[:, None]).mean(axis=0)))
        raw = float(PP[sel].mean()
                    / np.median(CC[sel].mean(axis=0)))
        say("  %-17d %-12.5f %-14.5f %-12.4f %.4f"
            % (a, v, cm, v / cm, raw))
    say("  The coin is square-root too, so both fall away from the")
    say("  bound at the same rate; what separates them is a constant.")

    say()
    say("  DIAGNOSTIC (post hoc). What this does and does not move.")
    say("  The bound is loose by the factor W2 prints, and that factor")
    say("  is the whole story: the truth is (N/k)^{1/2} and the bound")
    say("  is (N/k) damped by a sub-logarithmic factor, so the gap")
    say("  between them is itself a power of N/k. Numbers, at the")
    say("  largest N, as the ratio of the bound to the measurement:")
    say("  N/k octave        bound / |P|")
    for a, v in prof[len(NS) - 1]:
        say("  %-17d %.1f" % (a, 1.0 / v))
    say("  So the reduction is real in shape and empty in size: the")
    say("  elementary half is the half a classical estimate can reach,")
    say("  and the estimate reaches it %.0f orders of magnitude past"
        % base)
    say("  anything anyone will compute.")

    say()
    say("=" * 70)
    ok = w1 and w2 and w3 and w4
    say("the elementary half is the provable half, and the proof bites "
        "nowhere near here" if ok else "REFUTED")

    head = [
        "STATISTIC: the ratio of the measured |P(N;k)| to the classical",
        "           shape (N/k) exp(-c sqrt(log(N/k))) L(k) with",
        "           L(k) = prod_{p|k}(1-1/p)^{-1}, its octave profile,",
        "           the share of the Goldbach budget that shape alone",
        "           would spend at accessible N, and the N at which it",
        "           would pay that budget, swept over c.",
        "NULL: a coin arm on the identical sifted set -- 16 global",
        "      sign vectors on the odd squarefree m, each held across",
        "      all k as mu is, same mask, same octaves -- carried",
        "      through the same comparison against the same bound, so",
        "      that the share of the looseness owed to the shape can be",
        "      separated from the share owed to mu. Same convention as",
        "      lab_mu_vs_coin_size.py, so the two measurements of that",
        "      constant are comparable.",
        "FIELD: N = 2e5 to 3.2e6 by doubling; k squarefree and coprime",
        "       to N with 2 <= k < 30000; m odd, squarefree, coprime to",
        "       the odd part of k, m < N/k, avoiding N k^{-1} mod q for",
        "       every odd prime q <= 30 not dividing k; theta' = 0.56;",
        "       S(N) and A(N) from Euler products at the fixed bound",
        "       4000000; c swept over 0.15, 0.2098, 0.30.",
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
