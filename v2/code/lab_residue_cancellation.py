# -*- coding: utf-8 -*-
r"""
Does the residue buy cancellation, or is it small because its terms
are?

WHAT IS AT STAKE

Remark {#rem:provablehalf} moved the target: the elementary half has a
classical unconditional bound in the right shape and the residue has
none, so R is what the route turns on. Remark {#rem:residue} measured
R's size, |R| ~ (N/k)^{1/2}, and its lean. Neither says whether that
size is bought or given.

Write the residue out. With w(m,k) = C_k on the sifted set S and 0 off
it,

    R(N;k) = H - beta P = sum_m mu(m) delta(m,k),
    delta(m,k) = Lambda(N - mk) - beta w(m,k),

so R is a Mobius sum against the deviation of the von Mangoldt weight
from its own sieve prediction. Two very different worlds give
|R| ~ (N/k)^{1/2}:

  (i)  delta is already small -- sum_m |delta| is itself of order
       (N/k)^{1/2} -- and no cancellation happens at all;
  (ii) sum_m |delta| is of order N/k and mu cancels it down to the
       square root, which is the hard thing nobody can prove.

The distinction is the whole question of whether R is hard, and it is
settled by comparing |R| with the l1 and l2 norms of its own summands.
l2 is the scale a random sign pattern would reach; l1 is the scale of
no cancellation at all.

BACKS: Remark {#rem:residuecancel} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The control: beta refitted here reproduces the published
      AGREE beta_HP values to within 1 per cent, so this is the same
      split.
  V2  World (ii), not world (i): the octave means of the l1 norm
      sum_m |delta| are of order N/k, not (N/k)^{1/2} -- fitting them
      against N/k gives an exponent above 0.85 at every N.
  V3  And mu cancels that down to exactly what random signs would
      give: the octave means of |R| divided by the l2 norm
      (sum_m delta^2)^{1/2} lie in [0.5, 1.5] at every populated
      octave.
  V4  So the gain is square-root in the length: fitting the octave
      means of (l1 norm)/|R| against N/k gives an exponent in
      [0.35, 0.65] at every N.

REFUTATION RULE (fixed before the run)

  V1  REFUTED at 1 per cent at any N.
  V2  REFUTED if the l1 exponent is at or below 0.85 at any N. That
      would be world (i): the residue would be small because its terms
      are, and its square-root size would be no evidence of
      cancellation at all.
  V3  REFUTED if the ratio leaves [0.5, 1.5] anywhere. Below would
      mean mu does BETTER than random signs on delta -- extra
      cancellation, and R easier than feared. Above would mean it does
      worse.
  V4  REFUTED if the exponent leaves [0.35, 0.65] at any N.

  All four gate.

  NULL: 16 global sign vectors eps(m) = +-1 on the odd squarefree m,
  each held fixed across all k exactly as mu is, summed against the
  IDENTICAL delta. That is the control V3 is read against: it fixes
  what "what random signs would give" means on this particular
  sequence of deviations, rather than assuming the l2 norm times
  sqrt(2/pi). Same convention as lab_mu_vs_coin_size.py.
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
OUT = os.path.join(ROOT, "results", "lab_residue_cancellation.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
OCT = [2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152]
MINPTS = 10
COINS = 16
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


def read_beta():
    p = os.path.join(ROOT, "results", "lab_predictable_part.txt")
    src = io.open(p, encoding="utf-8").read()
    return {int(m.group(1)): float(m.group(2))
            for m in re.finditer(r"AGREE beta_HP N=(\d+) ([\d.]+)", src)}


def fit(cent, prof, cnt, floor):
    c = np.array(cent, dtype=float)
    y = np.array(prof, dtype=float)
    ok = (~np.isnan(c) & ~np.isnan(y) & (y > 0)
          & (np.array(cnt) >= floor))
    x = np.log(c[ok])
    yy = np.log(y[ok])
    return x, yy, float(np.polyfit(x, yy, 1)[0])


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_beta()
    say("read %d published beta values from "
        "results/lab_predictable_part.txt" % len(pub))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = sieves(NMAX)
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
        ks, Hs, Ps = [], [], []
        store = []
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
            L = lam[vals]
            w = np.ones(ms.size, dtype=np.float64)
            for q in QS:
                if k % q == 0:
                    continue
                w *= np.where(vals % q == 0, 0.0, q / (q - 1.0))
            ks.append(k)
            Hs.append(float((L * g).sum()))
            Ps.append(float((g * w).sum()))
            store.append((ms, g, L, w))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        P = np.array(Ps)
        beta = float((H * P).sum() / (P * P).sum())

        Rv, T1, T2, Cv = [], [], [], []
        for i, (ms, g, L, w) in enumerate(store):
            d = L - beta * w
            Rv.append(abs(float((g * d).sum())))
            T1.append(float(np.abs(d).sum()))
            T2.append(float(math.sqrt(float((d * d).sum()))))
            Cv.append(np.abs(coinmat[:, ms].astype(np.float64) @ d))
        res.append((N, ks, np.array(Rv), np.array(T1), np.array(T2),
                    np.array(Cv), N // ks, beta, H, P))
        say("  N = %-10d  #k = %-7d beta = %.6f" % (N, ks.size, beta))
        del store

    # ------------------------------------------------------------- V1
    say()
    say("V1  the control: is this the published split?")
    say("  N            beta (here)   beta (published)   ratio")
    v1 = True
    for N, ks, Rv, T1, T2, Cv, inner, beta, H, P in res:
        r = beta / pub[N]
        if abs(r - 1.0) >= 0.01:
            v1 = False
        say("  %-12d %-13.6f %-18.6f %.4f" % (N, beta, pub[N], r))
    say("  V1 %s" % ("hold" if v1 else "REFUTED"))

    def octaves(inner, v):
        cent, prof, cnt = [], [], []
        for a, b in zip(OCT, OCT[1:]):
            sel = (inner >= a) & (inner < b)
            cnt.append(int(sel.sum()))
            if sel.sum():
                cent.append(float(inner[sel].mean()))
                prof.append(float(np.asarray(v)[sel].mean()))
            else:
                cent.append(float("nan"))
                prof.append(float("nan"))
        return cent, prof, cnt

    # ------------------------------------------------------------- V2
    say()
    say("V2  how big is the sum of |delta| -- the no-cancellation size?")
    say("  N            l1 exponent   correlation   thinnest bin")
    v2 = True
    for N, ks, Rv, T1, T2, Cv, inner, beta, H, P in res:
        cent, prof, cnt = octaves(inner, T1)
        x, y, e = fit(cent, prof, cnt, MINPTS)
        r = float(np.corrcoef(x, y)[0, 1])
        if e <= 0.85:
            v2 = False
        say("  %-12d %-13.4f %-13.5f %d"
            % (N, e, r, min(c for c in cnt if c >= MINPTS)))
        say("POP residue_l1_N%d %d"
            % (N, min(c for c in cnt if c >= MINPTS)))
    say("  V2 %s   (floor 0.85)" % ("hold" if v2 else "REFUTED"))

    # ------------------------------------------------------------- V3
    say()
    say("V3  |R| against the l2 norm of its own summands, by octave")
    say("  and against 16 coins summed on the IDENTICAL delta")
    v3 = True
    for N, ks, Rv, T1, T2, Cv, inner, beta, H, P in res:
        say("  N = %d" % N)
        say("    N/k octave        #k     |R|/l2    coin/l2   |R|/coin")
        for a, b in zip(OCT, OCT[1:]):
            # a k whose whole m-set is sifted away has delta == 0 and
            # no ratio; it carries no information either way
            sel = (inner >= a) & (inner < b) & (T2 > 0)
            n = int(sel.sum())
            if n < MINPTS:
                continue
            q = float((Rv[sel] / T2[sel]).mean())
            cq = float(np.median((Cv[sel] / T2[sel][:, None])
                                 .mean(axis=0)))
            if not (0.5 <= q <= 1.5):
                v3 = False
            say("    [%-6d,%-7s) %-6d %-9.4f %-9.4f %.4f"
                % (a, str(b), n, q, cq, q / cq))
    say("  the scale both are near is E|Z|/sigma for a random sign")
    say("  sum, sqrt(2/pi) = %.4f" % math.sqrt(2.0 / math.pi))
    say("  V3 %s   (band [0.5, 1.5])" % ("hold" if v3 else "REFUTED"))

    # ------------------------------------------------------------- V4
    say()
    say("V4  the cancellation gain l1/|R| against the inner length")
    say("  N            gain exponent   correlation   gain at the top")
    v4 = True
    for N, ks, Rv, T1, T2, Cv, inner, beta, H, P in res:
        cent, prof, cnt = octaves(inner, T1 / np.maximum(Rv, 1e-12))
        x, y, e = fit(cent, prof, cnt, MINPTS)
        r = float(np.corrcoef(x, y)[0, 1])
        if not (0.35 <= e <= 0.65):
            v4 = False
        say("  %-12d %-15.4f %-13.5f %.1f"
            % (N, e, r, math.exp(float(y[-1]))))
        say("POP residue_gain_N%d %d"
            % (N, min(c for c in cnt if c >= MINPTS)))
        f = [float(np.polyfit(x[s], y[s], 1)[0])
             for s in (slice(None), slice(1, None), slice(0, -1))]
        say("SWEPT residue_gain_N%d octave-range %.4f"
            % (N, max(f) - min(f)))
        say("CORR residue_gain_N%d %.5f" % (N, abs(r)))
    say("  V4 %s   (band [0.35, 0.65])" % ("hold" if v4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). V4's instrument is the wrong one and")
    say("  its correlations say so: -0.90, -0.82, -0.76, -0.56, +0.98,")
    say("  wandering in sign. The gain l1/|R| is a ratio whose")
    say("  DENOMINATOR passes near zero -- a k at which R happens to")
    say("  cancel almost completely contributes an enormous gain -- so")
    say("  the octave mean of the ratio is set by outliers. The")
    say("  aggregate gain, mean(l1) over mean(|R|), has no such")
    say("  denominator problem and is what V4 meant to ask:")
    say("  N            aggregate gain exponent   correlation")
    for N, ks, Rv, T1, T2, Cv, inner, beta, H, P in res:
        c1, p1, n1 = octaves(inner, T1)
        c2, p2, n2 = octaves(inner, Rv)
        cent = c1
        prof = [a / b if b > 0 else float("nan")
                for a, b in zip(p1, p2)]
        x, y, e = fit(cent, prof, n1, MINPTS)
        r = float(np.corrcoef(x, y)[0, 1])
        say("  %-12d %-25.4f %.5f" % (N, e, r))
    say("  That is the square root the pre-registration was after, and")
    say("  it is what V3 already says in another form: the aggregate")
    say("  gain is l1/l2 times l2/|R|, and the second factor is the")
    say("  constant near 1/0.8 that a coin also gives.")

    say()
    say("  DIAGNOSTIC (post hoc). What delta actually looks like. It is")
    say("  Lambda(N-mk), which is log of a prime or zero, minus beta C_k")
    say("  on the sifted set and zero off it, so it takes essentially")
    say("  two values. At the largest N, by octave:")
    N, ks, Rv, T1, T2, Cv, inner, beta, H, P = res[-1]
    say("  N/k octave        l1          l2          l1/l2     sqrt(N/k)")
    for a, b in zip(OCT, OCT[1:]):
        sel = (inner >= a) & (inner < b)
        if sel.sum() < MINPTS:
            continue
        t1 = float(T1[sel].mean())
        t2 = float(T2[sel].mean())
        say("  %-17d %-11.1f %-11.1f %-9.2f %.1f"
            % (a, t1, t2, t1 / t2, math.sqrt(float(inner[sel].mean()))))
    say("  l1/l2 is the square root of the effective number of terms;")
    say("  compare it with sqrt(N/k) in the last column to see how much")
    say("  of the range actually carries a non-zero deviation.")

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("the residue buys square-root cancellation and buys exactly "
        "what random signs would" if ok else "REFUTED")

    head = [
        "STATISTIC: for each admissible k, the residue",
        "           R = sum_m mu(m) delta(m,k) with",
        "           delta = Lambda(N-mk) - beta w(m,k); the l1 norm",
        "           sum_m |delta| and the l2 norm (sum delta^2)^{1/2} of",
        "           its own summands; the octave means of |R|/l2 and of",
        "           the gain l1/|R|; and the exponents of l1 and of the",
        "           gain against the inner length N/k.",
        "NULL: 16 global sign vectors on the odd squarefree m, each held",
        "      across all k as mu is, summed against the IDENTICAL",
        "      delta, so that 'what random signs give' is measured on",
        "      this sequence of deviations rather than assumed from the",
        "      l2 norm. Same convention as lab_mu_vs_coin_size.py.",
        "FIELD: N = 2e5 to 3.2e6 by doubling; k squarefree and coprime",
        "       to N with 2 <= k < 30000; m odd, squarefree, coprime to",
        "       the odd part of k, m < N/k; the sieve weight uses the",
        "       odd primes up to 30; beta refitted here as",
        "       sum(H P)/sum(P^2) on the same k-range and checked",
        "       against results/lab_predictable_part.txt; octaves closed",
        "       at both ends and fitted only when they hold at least 10",
        "       k; numpy default_rng seed 20260808.",
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
