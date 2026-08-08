# -*- coding: utf-8 -*-
r"""
The control the split never had.

WHAT IS AT STAKE

Remark {#rem:predictable} reported that beta P carries 0.86 to 0.96 of
sum(log k)|H|, and Remark {#rem:residue} built the whole remaining
picture on the resulting split H = beta P + R.  Neither ran a null.
lab_survivor_selection.py permuted signs to test whether P predicts
the SIGN of H, which it does; nothing tested whether the SIZE claim
means anything.

It might not.  P is a smooth-ish sum of about N/(2k) terms with
bounded weights, and |H| is a sum of the same length.  Any predictor
of that shape, fitted by least squares, will absorb some fixed share
of the mass simply because both grow like the square root of the
length.  The question is whether mu is doing the work.

Two controls answer it, and they are the two ways to break P without
changing its shape:

  * the coin, mu(m) replaced by a fixed +-1 on the odd squarefree m,
    with the sieve weights w(m,k) untouched -- same length, same
    weights, no arithmetic;
  * the flat predictor, w(m,k) replaced by 1 with mu kept -- the plain
    Mobius sum over the coprime odd range, which isolates what the
    sieve weight itself is worth.

BACKS: Remark {#rem:splitnull} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  The coin cannot do it: its residual share
      sum(log k)|H - beta_c P_c| / sum(log k)|H| exceeds mu's by at
      least 0.20 at every N.
  U2  And it carries no correlation: r(H, P_coin) lies in [-0.1, 0.1]
      at every N, against about 0.85 for r(H, P).
  U3  The sieve weight earns its keep on the sign: r(H, P) exceeds
      r(H, P_flat) at every N.
  U4  And on the mass: P's residual share is below P_flat's at every N.

REFUTATION RULE (fixed before the run)

  U1  REFUTED if the gap is under 0.20 at any N. This is the one that
      decides whether the split of [rem:predictable] means anything:
      if a coin absorbs as much mass, the 0.86 to 0.96 is a property
      of the shape and not of mu.
  U2  REFUTED if the coin correlation leaves that band at any N.
  U3  REFUTED if the flat predictor matches or beats P at any N.
  U4  REFUTED likewise on the residual share.

  All four gate.  U3 and U4 failing would not void the split -- it
  would say the sieve weight is decoration and the plain Mobius sum is
  the predictor, which is a cleaner statement, not a weaker one.
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
OUT = os.path.join(ROOT, "results", "lab_split_null.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 30_000
QSIEVE = 30
COINS = 8
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


def shares(lw, H, P):
    """least-squares scale through the origin, and the two mass shares"""
    d = float((P * P).sum())
    b = float((H * P).sum() / d) if d > 0 else 0.0
    tot = float((lw * np.abs(H)).sum())
    return (b,
            float((lw * np.abs(H - b * P)).sum()) / tot,
            float((lw * np.abs(b * P)).sum()) / tot)


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
    say("  %d fixed coins on the odd squarefree m, sieve weights kept"
        % COINS)

    res = []
    for N in NS:
        PN = factor_set(N)
        ks, Hs, Ps, Fs = [], [], [], []
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
            Hs.append(float((lam[vals] * g).sum()))
            Ps.append(float((g * w).sum()))
            Fs.append(float(g.sum()))
            for j in range(COINS):
                Cs[j].append(float((coins[j][ms] * w).sum()))
        ks = np.array(ks, dtype=np.int64)
        lw = np.log(ks.astype(np.float64))
        Hs = np.array(Hs)
        Ps = np.array(Ps)
        Fs = np.array(Fs)
        Cs = [np.array(c) for c in Cs]
        res.append((N, ks, lw, Hs, Ps, Fs, Cs))
        say("  N = %-10d  #k = %d" % (N, ks.size))

    say()
    say("U1  what each predictor leaves behind")
    say("  N            mu residual   coin residual (median)   gap")
    u1 = True
    for N, ks, lw, Hs, Ps, Fs, Cs in res:
        _, rm, _ = shares(lw, Hs, Ps)
        rc = [shares(lw, Hs, c)[1] for c in Cs]
        med = float(np.median(rc))
        if med - rm < 0.20:
            u1 = False
        say("  %-12d %-13.4f %-24.4f %.4f" % (N, rm, med, med - rm))
    say("  U1 %s" % ("hold" if u1 else "REFUTED"))

    say()
    say("U2  correlation with H")
    say("  N            mu        coin min   coin max   flat")
    u2 = True
    for N, ks, lw, Hs, Ps, Fs, Cs in res:
        rmu = float(np.corrcoef(Hs, Ps)[0, 1])
        rc = [float(np.corrcoef(Hs, c)[0, 1]) for c in Cs]
        rf = float(np.corrcoef(Hs, Fs)[0, 1])
        if not (-0.1 <= min(rc) and max(rc) <= 0.1):
            u2 = False
        say("  %-12d %-9.4f %-10.4f %-10.4f %.4f"
            % (N, rmu, min(rc), max(rc), rf))
    say("  U2 %s" % ("hold" if u2 else "REFUTED"))

    say()
    say("U3/U4  is the sieve weight worth anything?")
    say("  N            r(H,P)    r(H,flat)   P residual   flat residual")
    u3 = u4 = True
    for N, ks, lw, Hs, Ps, Fs, Cs in res:
        rmu = float(np.corrcoef(Hs, Ps)[0, 1])
        rf = float(np.corrcoef(Hs, Fs)[0, 1])
        _, rp, _ = shares(lw, Hs, Ps)
        _, rr, _ = shares(lw, Hs, Fs)
        if rmu <= rf:
            u3 = False
        if rp >= rr:
            u4 = False
        say("  %-12d %-9.4f %-11.4f %-12.4f %.4f" % (N, rmu, rf, rp, rr))
    say("  U3 sieve weight wins on correlation   %s"
        % ("hold" if u3 else "REFUTED"))
    say("  U4 and on the residual share          %s"
        % ("hold" if u4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). Why U2's band was wrong. The naive")
    say("  sampling error of a correlation over n points is 1/sqrt(n),")
    say("  which at these n would be about %.4f -- but the inner sums"
        % (1.0 / math.sqrt(res[-1][1].size)))
    say("  are NESTED: consecutive k share almost all their m, so the")
    say("  effective degrees of freedom are the number of octaves, not")
    say("  the number of k. The coin's own spread is the empirical band,")
    say("  and mu sits far outside it:")
    say("  N            mu        coin band          mu - band edge")
    for N, ks, lw, Hs, Ps, Fs, Cs in res:
        rmu = float(np.corrcoef(Hs, Ps)[0, 1])
        rc = [float(np.corrcoef(Hs, c)[0, 1]) for c in Cs]
        say("  %-12d %-9.4f [%+.4f, %+.4f]   %+.4f"
            % (N, rmu, min(rc), max(rc), rmu - max(rc)))
    say("  So U2 is refuted on its stated band and not on its content:")
    say("  the band should have been computed from the coins, which is")
    say("  what M2 asks and what U1 does correctly.")

    say()
    say("  DIAGNOSTIC 2 (post hoc). The fitted scales, which say how much")
    say("  of |H| each predictor is being asked to carry:")
    say("  N            beta(mu)  beta(flat)  |beta P| share  "
        "|beta_f F| share")
    for N, ks, lw, Hs, Ps, Fs, Cs in res:
        bm, _, sm = shares(lw, Hs, Ps)
        bf, _, sf = shares(lw, Hs, Fs)
        say("  %-12d %-9.4f %-11.4f %-15.4f %.4f" % (N, bm, bf, sm, sf))

    say()
    say("  Cross-check lines. lab_predictable_part.py and")
    say("  lab_residue_size.py fit the same beta on the same k-range.")
    for N, ks, lw, Hs, Ps, Fs, Cs in res:
        say("AGREE beta_HP N=%d %.6f 0.01" % (N, shares(lw, Hs, Ps)[0]))

    say()
    say("=" * 70)
    ok = u1 and u2 and u3 and u4
    say("the split is mu's doing: a coin of the same shape absorbs "
        "nothing" if ok else "REFUTED")

    head = [
        "STATISTIC: for each of three predictors -- P = sum mu(m)w(m,k),",
        "           the same with mu replaced by a fixed coin on the odd",
        "           squarefree m, and the flat sum_m mu(m) with w = 1 --",
        "           the least-squares scale through the origin, the",
        "           residual share sum(log k)|H - beta P|/sum(log k)|H|,",
        "           and the correlation with H.",
        "NULL: the coin. mu(m) is replaced by a fixed +-1 on the odd",
        "      squarefree m, eight draws, with the sieve weights w(m,k)",
        "      and the summation range untouched, so the predictor keeps",
        "      its length, its weights and its shape and loses only the",
        "      arithmetic. This is the control [rem:predictable] and",
        "      [rem:residue] were built without.",
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
