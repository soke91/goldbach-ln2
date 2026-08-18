# -*- coding: utf-8 -*-
r"""
paper/wall_v3.md, Remark {#rem:supp} and the E1 ratios in
Section {#sec:supply} -- the type-II dilate field at N = 10^8.

WHAT IS UNDER TEST

    D(k)    = sum_{sqrt N < m <= N/k} mu(m) mu(N - mk),
    supp(k) = #{ m in that range : mu(m) != 0 and mu(N-mk) != 0 },
    M_k     = floor(N/k) - floor(sqrt N),

at N = 10^8 = 2^8 5^8, over three dyadic bands of k running up to
N^{1/3} = 464.  The remark says:

  (a) prod_p (1 - 2/p^2) = 0.32263 is NOT the right density, because at
      p = 2 and p = 5 the condition p^2 | N collapses the second
      congruence onto the first;
  (b) the true density supp(k)/M_k ranges from 0 to 0.594, and lands
      within 0.01 of 0.32263 for only 1, 1 and 3 of the 58, 116 and 232
      values of k;
  (c) D(k) vanishes identically exactly when 4 | k or 25 | k -- 2799 of
      the 9999 values k < sqrt N, or 28.0%;
  (d) aggregated, sum_k supp(k) / sum_k M_k is 0.3303, 0.3298, 0.3320 on
      the three bands, high by 2.2 to 2.9 percent.

Section {#sec:supply} adds that sum_k |D(k)|^2 / sum_k M_k is
0.34, 0.39, 0.32 over the three bands, and that normalised by the
surviving support instead it is 1.04, 1.17, 0.96 -- "exact square-root
cancellation on the terms that survive".  No script for any of it
exists here.

THE BANDS ARE NOT NAMED IN THE TEXT AND ARE INFERRED HERE

The text says "three dyadic bands" with 58, 116 and 232 values of k,
and E1 is stated for K <= N^{1/3}.  Since a dyadic band K <= k < 2K
holds K values, the bands are [58,116), [116,232), [232,464), and
464 = floor(N^{1/3}).  That inference is checked, not assumed: Z0
fails if the three bands do not have exactly 58, 116, 232 members.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z0  The bands [58,116), [116,232), [232,464) hold 58, 116, 232 values
      of k, and 464 = floor(N^{1/3}).
  Z1  prod_p (1 - 2/p^2) = 0.32263.
  Z2  max_k supp(k)/M_k = 0.594 over the three bands, and the minimum
      is exactly 0.
  Z3  |supp(k)/M_k - 0.32263| < 0.01 for exactly 1, 1 and 3 of the
      58, 116 and 232 values.
  Z4  #{k < sqrt N : 4|k or 25|k} = 2799 = 28.0% of 9999; and D(k) = 0
      exactly for every k in the bands with 4|k or 25|k, and for no
      other k in the bands.
  Z5  sum_k supp(k) / sum_k M_k = 0.3303, 0.3298, 0.3320, each between
      2.2% and 2.9% above 0.32263.
  Z6  sum_k |D(k)|^2 / sum_k M_k = 0.34, 0.39, 0.32, and
      sum_k |D(k)|^2 / sum_k supp(k) = 1.04, 1.17, 0.96.

REFUTATION RULE (fixed before the run)

  Z0  REFUTED if any band count differs, or floor(N^{1/3}) != 464.
  Z1  REFUTED if the product differs from 0.32263 by more than 0.000005.
  Z2  REFUTED if the maximum differs from 0.594 by more than 0.0005, or
      if the minimum is not 0.
  Z3  REFUTED if any of the three counts differs.
  Z4  REFUTED if the count is not 2799, or if the vanishing set is not
      exactly {4|k or 25|k} inside the bands.
  Z5  REFUTED if any entry differs by more than 0.00005, or if any
      excess falls outside [2.2%, 2.9%].
  Z6  REFUTED if any of the six entries differs by more than 0.005.

  All seven gate.

BACKS: the E1 arm of Conjecture {#conj:L} in paper/wall_v3.md, and
only that arm: the support density and the vanishing set of D(k).
The pair statistics, the (v_2,v_3) cells, the blind mask prediction
and the Wishart spectrum are NOT covered here.

CITED BY: {#rem:e1row} in paper/.
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
OUT = os.path.join(ROOT, "results", "audit_support_density.txt")

N = 100_000_000
BANDS = [(58, 116), (116, 232), (232, 464)]
PUB_MAXDENS = 0.594
PUB_NEAR = [1, 1, 3]
PUB_AGG = [0.3303, 0.3298, 0.3320]
PUB_R_M = [0.34, 0.39, 0.32]
PUB_R_S = [1.04, 1.17, 0.96]


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def mobius_big(n):
    """mu on [0,n] using the sqrt-sieve plus a cofactor array."""
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
    big = rem > 1                       # one surviving prime factor > sqrt n
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return mu


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    root = int(math.isqrt(N))
    cube = int(round(N ** (1.0 / 3.0)))
    while (cube + 1) ** 3 <= N:
        cube += 1
    while cube ** 3 > N:
        cube -= 1

    say("N = %d   sqrt N = %d   floor(N^{1/3}) = %d" % (N, root, cube))
    z0 = (cube == 464 and all(hi - lo == n for (lo, hi), n in
                              zip(BANDS, (58, 116, 232))))
    say("Z0  bands %s hold %s values   %s"
        % (BANDS, [hi - lo for lo, hi in BANDS],
           "hold" if z0 else "REFUTED"))

    prod = 1.0
    for p in primes_upto(2_000_000):
        prod *= 1.0 - 2.0 / (int(p) ** 2)
    z1 = abs(prod - 0.32263) <= 5e-6
    say("Z1  prod_p (1 - 2/p^2) = %.7f   published 0.32263   %s"
        % (prod, "hold" if z1 else "REFUTED"))

    say()
    say("sieving mu to %d ..." % N)
    mu = mobius_big(N)
    say("done; nonzero fraction = %.6f (6/pi^2 = %.6f)"
        % (float(np.count_nonzero(mu)) / N, 6.0 / math.pi ** 2))

    say()
    say("  band          k-range      sum M_k       sum supp     agg dens"
        "   pub      excess   sum|D|^2/sumM  pub    /sum supp  pub")
    say("  " + "-" * 110)
    aggs, rM, rS = [], [], []
    dens_all, zero_ok = [], True
    near = []
    percell = []
    for bi, (lo, hi) in enumerate(BANDS):
        sM = sS = sD2 = 0
        cnt = 0
        cell = []
        for k in range(lo, hi):
            top = N // k
            m = np.arange(root + 1, top + 1, dtype=np.int64)
            Mk = m.size
            a = mu[m]
            b = mu[N - m * k]
            surv = np.count_nonzero(a * b)
            d = int((a.astype(np.int64) * b.astype(np.int64)).sum())
            sM += Mk
            sS += surv
            sD2 += d * d
            dens = surv / Mk
            dens_all.append(dens)
            if abs(dens - 0.32263) < 0.01:
                cnt += 1
            dead = (k % 4 == 0) or (k % 25 == 0)
            if dead != (d == 0 and surv == 0):
                zero_ok = False
            cell.append((k, Mk, surv, d))
        percell.append(cell)
        near.append(cnt)
        agg = sS / sM
        aggs.append(agg)
        rM.append(sD2 / sM)
        rS.append(sD2 / sS)
        say("  %-13d [%d,%d)%s %-13d %-12d %-11.4f %-8.4f %-8.2f%% "
            "%-14.4f %-6.2f %-10.4f %.2f"
            % (bi + 1, lo, hi, " " * (6 - len(str(lo)) - len(str(hi))),
               sM, sS, agg, PUB_AGG[bi],
               100.0 * (agg / 0.32263 - 1.0), sD2 / sM, PUB_R_M[bi],
               sD2 / sS, PUB_R_S[bi]))

    mx, mn = max(dens_all), min(dens_all)
    z2 = abs(mx - PUB_MAXDENS) <= 5e-4 and mn == 0.0
    say()
    say("Z2  density range over the three bands: [%.6f, %.6f]   "
        "published [0, 0.594]   %s" % (mn, mx, "hold" if z2 else "REFUTED"))

    z3 = near == PUB_NEAR
    say("Z3  within 0.01 of 0.32263: %s   published %s   %s"
        % (near, PUB_NEAR, "hold" if z3 else "REFUTED"))

    ks = np.arange(1, root, dtype=np.int64)
    ndead = int(np.count_nonzero((ks % 4 == 0) | (ks % 25 == 0)))
    z4 = ndead == 2799 and zero_ok
    say("Z4  #{k < sqrt N : 4|k or 25|k} = %d of %d = %.1f%%; "
        "vanishing set exact inside the bands: %s   %s"
        % (ndead, root - 1, 100.0 * ndead / (root - 1), zero_ok,
           "hold" if z4 else "REFUTED"))

    e5 = max(abs(a - b) for a, b in zip(aggs, PUB_AGG))
    exc = [100.0 * (a / 0.32263 - 1.0) for a in aggs]
    z5 = e5 <= 5e-5 and all(2.2 <= e <= 2.9 for e in exc)
    say("Z5  aggregate density max |deviation| = %.6f (tol 0.00005); "
        "excess %s%%   %s"
        % (e5, ", ".join("%.2f" % e for e in exc),
           "hold" if z5 else "REFUTED"))

    e6 = max(max(abs(a - b) for a, b in zip(rM, PUB_R_M)),
             max(abs(a - b) for a, b in zip(rS, PUB_R_S)))
    z6 = e6 <= 5e-3
    say("Z6  E1 ratios max |deviation| = %.5f (tol 0.005)   %s"
        % (e6, "hold" if z6 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc, not a pre-registered test). Z5 and Z6 are")
    say("  aggregates, and Z0-Z4 -- which are per-k -- all hold, so the")
    say("  per-k field agrees and only the aggregation can differ. The")
    say("  text writes Z5 as a ratio of sums but \\S{#sec:R4} states that")
    say("  moduli with supp(k)=0 are excluded, and also distinguishes a")
    say("  ratio-of-sums from an unweighted mean of per-modulus ratios.")
    say("  All four combinations:")
    say("  band  ratio-of-sums   same, live k   mean of ratios   "
        "same, live k    published")
    for bi, cell in enumerate(percell):
        live = [c for c in cell if c[2] > 0]
        a1 = sum(c[2] for c in cell) / sum(c[1] for c in cell)
        a2 = sum(c[2] for c in live) / sum(c[1] for c in live)
        a3 = sum(c[2] / c[1] for c in cell) / len(cell)
        a4 = sum(c[2] / c[1] for c in live) / len(live)
        say("  %-5d %-15.4f %-14.4f %-16.4f %-15.4f %.4f"
            % (bi + 1, a1, a2, a3, a4, PUB_AGG[bi]))
    say("  and for the E1 ratio sum|D|^2 normalised four ways:")
    say("  band  /sum M_k       /sum M_k live  /sum supp        "
        "mean of per-k    published /M  /supp")
    for bi, cell in enumerate(percell):
        live = [c for c in cell if c[2] > 0]
        s2 = sum(c[3] * c[3] for c in cell)
        b1 = s2 / sum(c[1] for c in cell)
        b2 = s2 / sum(c[1] for c in live)
        b3 = s2 / sum(c[2] for c in cell)
        b4 = sum(c[3] * c[3] / c[2] for c in live) / len(live)
        say("  %-5d %-14.4f %-14.4f %-16.4f %-16.4f %-13.2f %.2f"
            % (bi + 1, b1, b2, b3, b4, PUB_R_M[bi], PUB_R_S[bi]))

    say()
    say("  DIAGNOSTIC 2 (post hoc). No aggregation convention fits, and the")
    say("  per-k field agrees, so the remaining free choice is where a")
    say("  dyadic band starts. 'K <= k < 2K' and 'K < k <= 2K' both hold K")
    say("  values, so the counts 58/116/232 do not pin the endpoints.")
    say("  Scanning the offset:")
    cache = {}
    for k in range(50, 470):
        top = N // k
        m = np.arange(root + 1, top + 1, dtype=np.int64)
        a = mu[m]
        b = mu[N - m * k]
        cache[k] = (m.size, int(np.count_nonzero(a * b)),
                    int((a.astype(np.int64) * b.astype(np.int64)).sum()))
    say("  offset  aggregate density              sum|D|^2/sum M_k")
    say("          band1   band2   band3          band1   band2   band3")
    best = None
    for off in range(-4, 5):
        rows, rows2 = [], []
        for lo, hi in BANDS:
            ks = range(lo + off, hi + off)
            sM = sum(cache[k][0] for k in ks)
            sS = sum(cache[k][1] for k in ks)
            s2 = sum(cache[k][2] ** 2 for k in ks)
            rows.append(sS / sM)
            rows2.append(s2 / sM)
        err = max(abs(a - b) for a, b in zip(rows, PUB_AGG))
        say("  %+3d     %.4f  %.4f  %.4f         %.4f  %.4f  %.4f  %s"
            % (off, rows[0], rows[1], rows[2], rows2[0], rows2[1],
               rows2[2], "<- matches" if err <= 5e-5 else ""))
        if best is None or err < best[0]:
            best = (err, off, rows, rows2)
    say("  closest offset %+d, max |deviation| %.5f" % (best[1], best[0]))
    say("  at that offset, sum|D|^2 normalised by the support instead:")
    say("    %s   against published %s"
        % (", ".join("%.4f" % (b / a) for a, b in zip(best[2], best[3])),
           ", ".join("%.2f" % v for v in PUB_R_S)))
    say("  and the range of sum|D|^2/sum M_k over all nine offsets:")
    for bi in range(3):
        vals = []
        for off in range(-4, 5):
            lo, hi = BANDS[bi]
            ks = range(lo + off, hi + off)
            vals.append(sum(cache[k][2] ** 2 for k in ks)
                        / sum(cache[k][0] for k in ks))
        say("    band %d: [%.4f, %.4f]   published %.2f"
            % (bi + 1, min(vals), max(vals), PUB_R_M[bi]))

    say()
    say("=" * 74)
    ok = z0 and z1 and z2 and z3 and z4 and z5 and z6
    say("Z0 %s  Z1 %s  Z2 %s  Z3 %s  Z4 %s  Z5 %s  Z6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (z0, z1, z2, z3, z4, z5, z6)))
    say("Remark {#rem:supp} and the E1 ratios reproduce" if ok
        else "REFUTED")

    head = [
        "STATISTIC: for the type-II dilate field at N = 10^8, per modulus",
        "           k: the surviving support supp(k) = #{m : mu(m) != 0 and",
        "           mu(N-mk) != 0} over sqrt N < m <= N/k, its ratio to the",
        "           range length M_k, the aggregate sum supp / sum M_k per",
        "           dyadic band, the count of k whose density lies within",
        "           0.01 of prod_p(1-2/p^2), the set of k with D(k) = 0,",
        "           and sum_k |D(k)|^2 normalised by sum M_k and by",
        "           sum supp(k).",
        "FIELD: N = 10^8 = 2^8 5^8; three dyadic bands k in [58,116),",
        "       [116,232), [232,464) with 464 = floor(N^{1/3}); m ranges",
        "       over (10^4, N/k]; mu on [0,10^8] by a sqrt-sieve with a",
        "       cofactor array, so the single prime factor above 10^4 is",
        "       picked up exactly; prod_p(1-2/p^2) over p < 2*10^6.",
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
