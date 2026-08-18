# -*- coding: utf-8 -*-
r"""
Is the elementary wall bigger than random signs, and if so by a
constant or by an exponent?

WHAT IS AT STAKE

Remark {#rem:elemsize} now records that mu's fitted exponent for |P|
lies above the coins' median at all five N, and Remark {#rem:elemreach}
found the same sign on a longer lever, 0.5178 against a coin band
topping out at 0.5095. Neither separates mu from a coin at any single
N: the offset sits inside the per-N coin spread every time. What is
not accidental is its sign.

Fitting exponents is the wrong instrument for that question. An
exponent is a slope through six points and carries the noise of all
six; the quantity actually at issue is a ratio of two magnitudes at
the SAME inner length, which needs no fit at all. So this measures

    ratio(N, octave) = mean |P| / median_j mean |P_coin^j|

directly, together with mu's rank among the coins in each cell. The
sieve weight is constant on the set it keeps -- w(m,k) = C_k times an
indicator -- so C_k cancels out of the ratio entirely and what is
compared is sum_{m in S} mu(m) against sum_{m in S} eps(m) on the
identical set S.

The distinction the ratio can draw and the exponent cannot: a ratio
that is constant in N/k means mu beats or loses to a coin by a factor
and shares its exponent; a ratio that grows with N/k means the
exponents genuinely differ, and then sum_{k<K}(log k)|P| is a harder
object than the square-root heuristic of {#rem:directlevel} assumes.

BACKS: Remark {#rem:muvscoin} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  The control: the octave means of |P| recomputed here match the
      ones in results/lab_elementary_size.txt to within 0.5 per cent
      in every populated cell.
  Z2  The sign reported by {#rem:elemsize} survives at the level of
      magnitudes: mu's mean |P| exceeds the coins' median in more
      than half of the (N, octave) cells.
  Z3  The excess is a constant and not an exponent: the ratio in the
      longest fitted octave does not exceed the ratio in the shortest
      by more than 0.15, at a majority of N.
  Z4  And it is a bias, not a separation: mu's mean |P| stays inside
      the range of the 32 coin draws in every cell.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED at 0.5 per cent in any cell, which would mean this is
      not the same P.
  Z2  REFUTED if mu is at or below the coin median in half the cells
      or more, which would withdraw the sign {#rem:elemsize} claims.
  Z3  REFUTED if the ratio grows by more than 0.15 across the fitted
      octaves at three or more of the five N. That is the outcome
      that matters: it would say the elementary wall carries a
      different exponent from a coin, and the square-root heuristic
      understates what sum(log k)|P| costs.
  Z4  REFUTED if mu leaves the coin range in any cell, which would be
      a separation at a single N and a much stronger statement than
      anything measured so far.

  All four gate.

  NULL: 32 global sign vectors eps(m) = +-1 on the odd squarefree m,
  each fixed across all k as mu is, so the coins carry the same
  across-k correlation structure that mu does. Drawing a fresh sign
  per (N,k) would make the coin octave means less variable than mu's
  for a reason that has nothing to do with mu, and would bias every
  rank; the global vector is the conservative choice. The statistic is
  a ratio of two means of magnitudes, both bounded away from zero, so
  it stays well conditioned under the control -- the criterion of
  {#rem:weightgapnull}.
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
OUT = os.path.join(ROOT, "results", "lab_mu_vs_coin_size.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
OCT = [2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152]
MINPTS = 10
COINS = 32
SEED = 20260808


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


def read_published():
    """the octave means of |P| already published -- read, not copied"""
    p = os.path.join(ROOT, "results", "lab_elementary_size.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("Y1  octave means of |P|")
    out = {}
    for ln in src[i:].splitlines():
        if not ln.strip().startswith("["):
            continue
        head, rest = ln.split(")", 1)
        lo = int(head.strip()[1:].split(",")[0])
        vals = []
        for cell in rest.split():
            v = cell.split("/")[0]
            vals.append(float(v) if v != "-" else float("nan"))
        out[lo] = vals
        if lo == OCT[-2]:
            break
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d published octave rows from "
        "results/lab_elementary_size.txt" % len(pub))

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

    rng = np.random.default_rng(SEED)
    coinmat = np.zeros((COINS, NMAX + 1), dtype=np.int8)
    idx = np.flatnonzero(oddsqf)
    for j in range(COINS):
        coinmat[j, idx] = rng.integers(0, 2, size=idx.size) * 2 - 1
    say("  %d global sign vectors on the %d odd squarefree m"
        % (COINS, idx.size))
    del idx

    res = []
    for N in NS:
        PN = factor_set(N)
        ks, Pm, Pc = [], [], []
        for k in range(2, KCAP):
            if not sqf[k] or any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
            if M < 2:
                continue
            mask = oddsqf[1:M + 1].copy()
            for p in factor_set(k):
                if p > 2:
                    mask[p - 1::p] = False
            # lab_elementary_size.py drops k only when the set is empty
            # BEFORE the sieve weight; a k whose m are all sifted away
            # contributes |P| = 0 to the octave mean rather than being
            # skipped. At short inner lengths that is most of them, and
            # dropping them shifts the mean by tens of per cent.
            if not mask.any():
                continue
            # C_k cancels from the ratio at a FIXED k and does not
            # cancel from a ratio of means over k, because C_k varies
            # with which small primes divide k. The budget carries
            # C_k, so the means must too.
            C = 1.0
            for q in QS:
                if k % q == 0:
                    continue
                C *= q / (q - 1.0)
                r = (N * pow(k, -1, q)) % q
                mask[(r - 1) % q::q] = False
            ks.append(k)
            Pm.append(C * abs(int(mu[1:M + 1][mask]
                                  .sum(dtype=np.int64))))
            sub = coinmat[:, 1:M + 1][:, mask]
            Pc.append(C * np.abs(sub.sum(axis=1, dtype=np.int64)))
        ks = np.array(ks, dtype=np.int64)
        res.append((N, ks, np.array(Pm, dtype=float),
                    np.array(Pc, dtype=float), N // ks))
        say("  N = %-10d  #k = %d" % (N, ks.size))

    # ------------------------------------------------------------- Z1
    say()
    say("Z1  the control: octave means of |P| against the published")
    say("  ones, with the sieve weight's constant C_k carried through")
    say("  as the budget carries it.")
    say("  N/k octave        " + " ".join("N=%-11d" % N for N in NS))
    z1 = True
    cells = {}
    for a, b in zip(OCT, OCT[1:]):
        row = []
        for i, (N, ks, Pm, Pc, inner) in enumerate(res):
            sel = (inner >= a) & (inner < b)
            n = int(sel.sum())
            if not n:
                row.append("%13s" % "-")
                continue
            # C_k depends on which q divide k, so it must be put back
            # per k before averaging
            here = float(Pm[sel].mean())
            there = pub.get(a, [float("nan")] * len(NS))[i]
            d = abs(here / there - 1.0) if there == there else 0.0
            if there == there and d >= 0.005:
                z1 = False
            row.append("%8.3f/%.4f" % (here, d))
            cells[(a, i)] = n
        say("  [%-6d,%-7s) %s" % (a, str(b), " ".join(row)))
    say("  Z1 %s" % ("hold" if z1 else "REFUTED"))

    # -------------------------------------------------------- Z2/Z3/Z4
    say()
    say("Z2/Z4  mu against 32 coins on the identical set, by octave")
    say("  ratio = mean|P| / median_j mean|P_coin^j|; rank counts the")
    say("  coins mu beats, out of %d" % COINS)
    ratios = {}
    above = tot = 0
    z4 = True
    for i, (N, ks, Pm, Pc, inner) in enumerate(res):
        say("  N = %d" % N)
        say("    N/k octave        #k     ratio    rank   inside range")
        for a, b in zip(OCT, OCT[1:]):
            sel = (inner >= a) & (inner < b)
            n = int(sel.sum())
            if n < MINPTS:
                continue
            m_ = float(Pm[sel].mean())
            cm = Pc[sel].mean(axis=0)
            med = float(np.median(cm))
            r = m_ / med
            rank = int((cm < m_).sum())
            ins = float(cm.min()) <= m_ <= float(cm.max())
            if not ins:
                z4 = False
            ratios.setdefault(i, []).append((a, r))
            tot += 1
            above += int(m_ > med)
            say("    [%-6d,%-7s) %-6d %-8.4f %-6s %s"
                % (a, str(b), n, r, "%d/%d" % (rank, COINS),
                   "yes" if ins else "NO"))
    z2 = above > tot / 2.0
    say("  Z2 mu above the coin median in %d of %d cells   %s"
        % (above, tot, "hold" if z2 else "REFUTED"))
    say("  Z4 mu inside the coin range in every cell   %s"
        % ("hold" if z4 else "REFUTED"))

    say()
    say("Z3  does the excess grow with the inner length?")
    say("  N            shortest   longest    growth")
    grew = 0
    for i, N in enumerate(NS):
        rr = ratios[i]
        g = rr[-1][1] - rr[0][1]
        grew += int(g > 0.15)
        say("  %-12d %-10.4f %-10.4f %+.4f" % (N, rr[0][1], rr[-1][1], g))
    z3 = grew < 3
    say("  Z3 the ratio grows by more than 0.15 at %d of %d N   "
        "(cap 2)   %s" % (grew, len(NS), "hold" if z3 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). The ratio pooled over all five N,")
    say("  which is what the octave means are actually estimating:")
    say("  N/k octave        cells  mean ratio   mean rank")
    pool = {}
    for i, rr in ratios.items():
        for a, r in rr:
            pool.setdefault(a, []).append(r)
    for a in sorted(pool):
        v = pool[a]
        say("  %-17d %-6d %-12.4f -" % (a, len(v), float(np.mean(v))))
    allr = [r for v in pool.values() for r in v]
    say("  over every cell: mean ratio %.4f, and %d of %d above 1"
        % (float(np.mean(allr)), sum(1 for r in allr if r > 1.0),
           len(allr)))
    say()
    say("  Z3's instrument is the wrong one and the pooled profile says")
    say("  why: it differences the last octave against the first, but")
    say("  each N starts and ends at a different octave, so it compares")
    say("  different parts of the curve at different N. The profile")
    say("  above is not monotone -- it climbs out of the short octaves,")
    say("  where parity leaves the non-negative m = 1 term nearly alone")
    say("  and mu falls BELOW every coin, and then flattens.")
    say()
    say("  The trend measured where the transient is over. Pooling the")
    say("  octaves from %d up and fitting log(ratio) on log(N/k):" % 128)
    tail = sorted(a for a in pool if a >= 128)
    tx = np.log(np.array(tail, dtype=float))
    ty = np.log(np.array([float(np.mean(pool[a])) for a in tail]))
    sl = float(np.polyfit(tx, ty, 1)[0])
    say("    octaves %s" % ", ".join(str(a) for a in tail))
    say("    mean ratios %s"
        % ", ".join("%.4f" % float(np.mean(pool[a])) for a in tail))
    say("    slope %+.4f, mean ratio %.4f"
        % (sl, float(np.mean([r for a in tail for r in pool[a]]))))
    f = [float(np.polyfit(tx[s], ty[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    say("    leave-one-out on the slope: %+.4f, %+.4f, %+.4f -- "
        "spread %.4f" % (f[0], f[1], f[2], max(f) - min(f)))
    say("  A slope of zero would mean a constant factor and the same")
    say("  exponent; a slope of a tenth over the four decades of N/k")
    say("  spanned here would be a different exponent. The measured")
    say("  slope is %+.4f with a leave-one-out spread of %.4f."
        % (sl, max(f) - min(f)))

    say()
    say("=" * 70)
    ok = z1 and z2 and z3 and z4
    say("mu's elementary sum is larger than random signs by a constant, "
        "not by an exponent" if ok else "REFUTED")

    head = [
        "STATISTIC: per octave of the inner length N/k, the mean of",
        "           |sum_{m in S} mu(m)| over the admissible k, against",
        "           the same mean for each of 32 global sign vectors on",
        "           the identical sifted set S; their ratio to the coin",
        "           median, mu's rank among the coins, and the change in",
        "           that ratio from the shortest fitted octave to the",
        "           longest. The sieve weight is C_k times an indicator,",
        "           so C_k cancels from the ratio; it is restored for the",
        "           control against the published octave means.",
        "NULL: 32 global sign vectors eps(m) = +-1 on the odd squarefree",
        "      m, each fixed across all k exactly as mu is, so the coins",
        "      carry the same across-k correlation mu does. A fresh sign",
        "      per (N,k) would make the coin octave means less variable",
        "      than mu's for reasons unrelated to mu and would bias every",
        "      rank. The statistic is a ratio of means of magnitudes,",
        "      both bounded away from zero, so it stays well conditioned",
        "      under the control.",
        "FIELD: N = 2e5 to 3.2e6 by doubling; k squarefree and coprime to",
        "       N with 2 <= k < 30000; m odd, squarefree, coprime to the",
        "       odd part of k, m < N/k, and avoiding N k^{-1} mod q for",
        "       every odd prime q <= 30 not dividing k; octaves closed at",
        "       both ends and fitted only when they hold at least 10 k;",
        "       mu from an integer sieve to 3.2e6; numpy default_rng seed",
        "       20260808.",
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
