# -*- coding: utf-8 -*-
r"""
Where the sign lean of the dilated walls comes from.

WHAT IS AT STAKE

Remark {#rem:signmass} located the failure of cancellation across
dilations: counting moduli the signs of H(N;k) are balanced, but
weighting by contribution, 68 to 78 percent of the mass is negative,
and splitting at the median of |H| the large dilated walls are the
ones that lean.  Remark {#rem:leandecay} then showed the lean shrinks
as N grows and called it a finite-N effect.  Neither says WHY the
lean is there, and without a mechanism "finite-N effect" is a label.

[eq:dilate] supplies the mechanism to test.  H(N;k) sums m over
1 <= m < N/k, so large k means a SHORT inner sum, and over a short
range sum_{m<=M} mu(m) is the Mertens function M(M), which is negative
for most small M.  If Lambda(N-mk) were flat in m the sign of H(N;k)
would simply be the sign of M(floor(N/k)).  The prediction is
therefore concrete: the lean should live at large k, weaken as N/k
grows, and the sign of H(N;k) should track the sign of the Mertens
function at N/k.

BACKS: Remark {#rem:leanmertens} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  P1  In the shortest inner sums, 2 <= N/k < 4, the mass-weighted
      fraction f+ of k with H(N;k) > 0 is under 0.35 at every N.
  P2  It rises with the length of the inner sum: f+ is monotone
      increasing across the octaves of N/k at every N.
  P3  The sign of H tracks the Mertens function: over 2 <= N/k <= 1000
      the fraction of k with sign H(N;k) = sign M(floor(N/k)) exceeds
      the largest of 16 permutation draws at every N.
  P4  The lean is concentrated where the inner sum is short:
      |1/2 - f+| restricted to N/k >= 1000 is smaller than the same
      quantity restricted to N/k < 1000, at every N.

REFUTATION RULE (fixed before the run)

  P1  REFUTED if f+ reaches 0.35 in that octave at any N.
  P2  REFUTED by a single fall across the octaves.
  P3  REFUTED if the agreement fails to beat every draw at any N. This
      is the one that decides the mechanism: without it the lean is
      not the Mertens function and the explanation is wrong.
  P4  REFUTED if the long-inner-sum lean is at least as large.

  All four gate.

  THE CONTROL is P3's permutation. The Mertens signs attached to the
  distinct values of floor(N/k) are shuffled among themselves, 16
  times, which destroys the claimed correspondence while preserving
  both marginal sign distributions exactly -- so the baseline it
  gives already accounts for the fact that most H are negative AND
  most M(x) are negative, which a naive 1/2 would not.
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
OUT = os.path.join(ROOT, "results", "lab_lean_mechanism.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
OCT = [2, 4, 8, 16, 32, 64, 128, 256, 1024, 4096]
MCAP = 1000
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


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    mert = np.cumsum(mu.astype(np.int64))          # M(x) = sum_{m<=x} mu(m)

    rows = []
    for N in NS:
        PN = factor_set(N)
        ks = np.array([k for k in range(2, N // 2)
                       if sqf[k] and all(k % q for q in PN)],
                      dtype=np.int64)
        f0 = np.zeros(N, dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        f0[1:] = lam[1:N] * mu[N - idx].astype(np.float64)
        A = np.empty(ks.size, dtype=np.float64)
        for i, k in enumerate(ks):
            k = int(k)
            r = N % k
            A[i] = f0[r::k].sum() if r else f0[k::k].sum()
        del f0
        H = mu[ks].astype(np.float64) * A          # [eq:dilate]
        w = np.log(ks.astype(np.float64)) * np.abs(H)
        Minner = N // ks                            # length of the m-sum
        rows.append((N, ks, H, w, Minner))
        say("  N = %-10d  #k = %d" % (N, ks.size))

    say()
    say("P1/P2  mass-weighted f+ by the length of the inner sum N/k")
    say("  range of N/k   " + "  ".join("N=%-9d" % N for N in NS))
    tab = []
    for a, b in zip(OCT, OCT[1:] + [1 << 62]):
        row = []
        for N, ks, H, w, Minner in rows:
            sel = (Minner >= a) & (Minner < b)
            tot = float(w[sel].sum())
            row.append(float(w[sel & (H > 0)].sum()) / tot
                       if tot > 0 else float("nan"))
        tab.append((a, b, row))
        say("  [%-6d,%-7s) %s"
            % (a, "inf" if b > NMAX else str(b),
               "  ".join("%-11.4f" % v for v in row)))
    p1 = all(v < 0.35 for v in tab[0][2])
    say("  P1  f+ under 0.35 in the shortest octave   %s"
        % ("hold" if p1 else "REFUTED"))
    p2 = True
    for j in range(len(NS)):
        col = [t[2][j] for t in tab if not math.isnan(t[2][j])]
        if not all(col[i] < col[i + 1] for i in range(len(col) - 1)):
            p2 = False
    say("  P2  monotone across the octaves   %s"
        % ("hold" if p2 else "REFUTED"))

    say()
    say("P3  does the sign of H track the Mertens function at N/k?")
    say("  N            #k in range   agreement   draws max   draws median")
    p3 = True
    rng = np.random.default_rng(SEED)
    for N, ks, H, w, Minner in rows:
        sel = (Minner >= 2) & (Minner <= MCAP) & (H != 0)
        mm = Minner[sel]
        sh = np.sign(H[sel])
        sm = np.sign(mert[mm]).astype(np.float64)
        agree = float((sh == sm).mean())
        uniq = np.unique(mm)
        base = np.searchsorted(uniq, mm)
        vals = np.sign(mert[uniq]).astype(np.float64)
        got = []
        for _ in range(DRAWS):
            got.append(float((sh == vals[rng.permutation(uniq.size)][base])
                             .mean()))
        got = np.array(got)
        if agree <= float(got.max()):
            p3 = False
        say("  %-12d %-13d %-11.4f %-11.4f %.4f"
            % (N, int(sel.sum()), agree, float(got.max()),
               float(np.median(got))))
    say("  P3 %s" % ("hold" if p3 else "REFUTED"))

    say()
    say("P4  where the lean lives")
    say("  N            |1/2-f+| at N/k<1000   at N/k>=1000")
    p4 = True
    for N, ks, H, w, Minner in rows:
        out = []
        for sel in ((Minner < MCAP), (Minner >= MCAP)):
            tot = float(w[sel].sum())
            f = float(w[sel & (H > 0)].sum()) / tot if tot > 0 else 0.5
            out.append(abs(0.5 - f))
        if out[1] >= out[0]:
            p4 = False
        say("  %-12d %-21.4f %.4f" % (N, out[0], out[1]))
    say("  P4 %s" % ("hold" if p4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). What the profile actually is.")
    say("  The premise of the Mertens reading is only half true --")
    say("  M(x) is negative at most small x but not at all of them:")
    say("  range of x   fraction with M(x) < 0")
    for lo, hi in ((2, 30), (2, 100), (2, 1000), (2, 10000)):
        neg = float((mert[lo:hi + 1] < 0).mean())
        say("  [%-5d,%-6d] %.4f" % (lo, hi, neg))
    say("  and P3 shows the correspondence is absent anyway: the")
    say("  agreement sits near the permutation median, and the best")
    say("  draw beats it at every N. The lean is not the Mertens")
    say("  function seen through a short window.")
    say()
    say("  What the octave table shows instead is a definite profile in")
    say("  N/k with a POSITIVE lean at the shortest sums, a minimum in")
    say("  the middle, and a slow return toward 1/2:")
    say("  N            argmin octave   f+ there   f+ at [2,4)")
    for j, (N, ks, H, w, Minner) in enumerate(rows):
        col = [t[2][j] for t in tab]
        i0 = int(np.nanargmin(np.array(col)))
        say("  %-12d [%-6d,%-7s) %-10.4f %.4f"
            % (N, tab[i0][0],
               "inf" if tab[i0][1] > NMAX else str(tab[i0][1]),
               col[i0], col[0]))
    say("  The positive lean at N/k < 4 has an elementary cause: N is")
    say("  even, so N - mk is even whenever mk is, and Lambda vanishes")
    say("  on even numbers above 2. For odd k only odd m survive, so")
    say("  the m = 1 term stands almost alone there and it is")
    say("  nonnegative. That is parity, not cancellation, and it says")
    say("  the short-sum end of the profile carries no information")
    say("  about the wall.")
    say()
    say("  Cross-check lines. lab_lean_oddmertens.py recomputes the")
    say("  same octave profile while testing the odd-Mertens")
    say("  predictor against it.")
    for j, N in enumerate(NS):
        say("AGREE lean_f_16_32 N=%d %.6f 0.02" % (N, tab[3][2][j]))

    say()
    say("=" * 70)
    ok = p1 and p2 and p3 and p4
    say("the sign lean of the dilated walls is the Mertens function at "
        "the length of the inner sum" if ok else "REFUTED")

    head = [
        "STATISTIC: the mass-weighted fraction f+ of k with H(N;k) > 0,",
        "           binned by the length N/k of the inner sum; the",
        "           fraction of k at which sign H(N;k) equals",
        "           sign M(floor(N/k)) with M the Mertens function, over",
        "           2 <= N/k <= 1000; the same under 16 permutations of",
        "           the Mertens signs; and |1/2 - f+| split at N/k = 1000.",
        "NULL: the permutation of P3. The Mertens signs attached to the",
        "      distinct values of floor(N/k) are shuffled among",
        "      themselves, which destroys the claimed correspondence",
        "      while preserving both marginal sign distributions exactly.",
        "      A naive 1/2 would be the wrong baseline here, because most",
        "      H are negative and most M(x) are negative, so agreement is",
        "      high by chance; the permutation accounts for that.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k over the squarefree",
        "       k < N/2 coprime to N; H(N;k) = mu(k)A(N;k) by [eq:dilate];",
        "       weights (log k)|H(N;k)|; seed 20260808.",
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
