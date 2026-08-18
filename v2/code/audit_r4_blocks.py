# -*- coding: utf-8 -*-
r"""
paper/wall_v3.md, Section {#sec:R4} -- "the divisor switch does not
localize".

WHAT IS UNDER TEST

With D(k) = sum_{sqrt N < m <= N/k} mu(m) mu(N-mk) and supp(k) the
number of surviving terms, over the whole band k < sqrt N on which the
type-II field is non-empty, the section states:

  (a) the exact identity
        sum_{k>=1} sum_{m: mk<=N-1} mu(m) mu(N-mk)
          = sum_{u<N} mu(N-u) sum_{m|u} mu(m) = mu(N-1),
      "verified by brute force at six values of N";
  (b) moduli with supp(k) = 0 are excluded from the blocks and from the
      autocorrelation alike, and there are 2799 such k of 9999;
  (c) r(B) = < S_B(j)^2 / sum_{k in block j} supp(k) >_j, the unweighted
      mean of per-block ratios, reads 1.024, 1.017, 1.016, 1.024 at
      B = 1, 2, 4, 8;
  (d) the ratio-of-sums normalisation instead gives a B=1 baseline of
      1.78, "because the five smallest k carry a third of
      sum_k supp(k)";
  (e) the lag-1 autocorrelation of D(k)/sqrt(supp(k)) over the surviving
      moduli reads +0.0055 against a 400-draw permutation null of
      standard deviation 0.0121, i.e. +0.5 standard errors;
  (f) keeping the dead moduli with D(k)/sqrt(supp(k)) set to 0, the same
      statistic reads -0.0108 against its own null of 0.0094, or -1.2
      standard errors;
  (g) consecutive surviving moduli are 1, 2 or 3 apart, so "lag-1" names
      position in the surviving sequence and not a fixed gap in k.

No script for any of it exists here.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  D1  The identity of (a) holds exactly at six values of N.
  D2  The surviving moduli number 9999 - 2799 = 7200.
  D3  r(B) = 1.024, 1.017, 1.016, 1.024 at B = 1, 2, 4, 8.
  D4  The ratio-of-sums baseline at B = 1 is 1.78, and the five
      smallest surviving k carry between 30% and 36% of sum_k supp(k).
  D5  The lag-1 autocorrelation over surviving moduli is +0.0055, its
      400-draw permutation null has standard deviation 0.0121, and the
      ratio is +0.5.
  D6  With dead moduli kept at 0 the statistic is -0.0108, its null
      0.0094, and the ratio -1.2.
  D7  Consecutive surviving moduli are 1, 2 or 3 apart -- the maximum
      gap is exactly 3.

REFUTATION RULE (fixed before the run)

  D1  REFUTED by a single N where the two sides differ.
  D2  REFUTED if the count is not 7200.
  D3  REFUTED if any of the four differs by more than 0.0005.
  D4  REFUTED if the baseline differs from 1.78 by more than 0.005, or
      if the five-smallest share falls outside [30%, 36%].
  D5  REFUTED if the point estimate differs by more than 0.00005, the
      null standard deviation by more than 0.001, or the ratio by more
      than 0.1.
  D6  Same tolerances as D5 against -0.0108, 0.0094, -1.2.
  D7  REFUTED if the maximum gap is not 3.

  All seven gate.  The permutation nulls use a fixed seed, so D5 and D6
  are reproducible here; they are not expected to match the source draw
  for draw, only in standard deviation, which is why the tolerance on
  the null is looser than on the point estimate.

CITED BY: {#rem:r4null} in paper/.
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
OUT = os.path.join(ROOT, "results", "audit_r4_blocks.txt")

N = 100_000_000
PUB_R = {1: 1.024, 2: 1.017, 4: 1.016, 8: 1.024}
PUB_ROS = 1.78
PUB_AC = (0.0055, 0.0121, 0.5)
PUB_AC0 = (-0.0108, 0.0094, -1.2)
DRAWS = 400
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def mobius_big(n):
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


def lag1(x):
    x = x - x.mean()
    d = float((x * x).sum())
    return float((x[:-1] * x[1:]).sum()) / d if d else 0.0


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    # ------------------------------------------------------------- D1
    say("D1  the exact identity, brute force")
    say("=" * 70)
    d1 = True
    small = mobius_big(200_000)
    for Nv in (1000, 4096, 10_000, 30_030, 65_536, 100_000):
        tot = 0
        for k in range(1, Nv):
            mm = (Nv - 1) // k
            if mm < 1:
                break
            a = small[1:mm + 1].astype(np.int64)
            b = small[Nv - k:Nv - k * mm - 1:-k].astype(np.int64)
            tot += int((a * b[:a.size]).sum())
        want = int(small[Nv - 1])
        ok = tot == want
        d1 = d1 and ok
        say("  N = %-8d  sum = %-4d  mu(N-1) = %-4d  %s"
            % (Nv, tot, want, "ok" if ok else "MISMATCH"))
    say("  D1 %s" % ("hold" if d1 else "REFUTED"))
    del small

    say()
    say("sieving mu to %d ..." % N)
    mu = mobius_big(N)
    root = int(math.isqrt(N))

    say("computing D(k) and supp(k) for k = 1 .. %d ..." % (root - 1))
    D = np.zeros(root, dtype=np.int64)
    SU = np.zeros(root, dtype=np.int64)
    for k in range(1, root):
        m1 = N // k
        if m1 <= root:
            continue
        a = mu[root + 1:m1 + 1]
        stop = N - k * m1 - 1
        b = mu[N - k * (root + 1)::-k][:a.size] if stop < 0 else \
            mu[N - k * (root + 1):stop:-k]
        b = b[:a.size]
        pr = a.astype(np.int8) * b
        SU[k] = int(np.count_nonzero(pr))
        D[k] = int(pr.astype(np.int64).sum())

    ks = np.arange(1, root, dtype=np.int64)
    su, dd = SU[1:root], D[1:root]
    live = su > 0
    nlive = int(live.sum())
    d2 = nlive == 7200
    say()
    say("D2  surviving moduli = %d of %d (dead %d)   %s"
        % (nlive, root - 1, root - 1 - nlive, "hold" if d2 else "REFUTED"))

    gaps = np.diff(ks[live])
    d7 = int(gaps.max()) == 3
    say("D7  gaps between consecutive surviving moduli: max %d, "
        "values %s   %s"
        % (int(gaps.max()), sorted(set(int(g) for g in gaps)),
           "hold" if d7 else "REFUTED"))

    # ------------------------------------------------------------- D3/D4
    dl, sl = dd[live], su[live]
    say()
    say("D3  block ratios over the surviving moduli")
    say("  B    blocks   r(B) = <S_B^2 / sum supp>_j    published")
    got = {}
    for B in (1, 2, 4, 8):
        nb = nlive // B
        s = dl[:nb * B].reshape(nb, B).sum(axis=1).astype(np.float64)
        w = sl[:nb * B].reshape(nb, B).sum(axis=1).astype(np.float64)
        got[B] = float((s * s / w).mean())
        say("  %-4d %-8d %-30.6f %.3f" % (B, nb, got[B], PUB_R[B]))
    e3 = max(abs(got[B] - PUB_R[B]) for B in PUB_R)
    d3 = e3 <= 5e-4
    say("  D3  max |deviation| = %.6f  (tol 0.0005)   %s"
        % (e3, "hold" if d3 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). B=1 reproduces exactly, so D(k) and")
    say("  supp(k) agree per modulus and only the blocking can differ.")
    say("  'Block' is not defined in the text. Two constructions:")
    say("  B    consecutive survivors   k-ranges of width B   published")
    for B in (1, 2, 4, 8):
        nb = nlive // B
        s = dl[:nb * B].reshape(nb, B).sum(axis=1).astype(np.float64)
        w = sl[:nb * B].reshape(nb, B).sum(axis=1).astype(np.float64)
        v1 = float((s * s / w).mean())
        edges = np.arange(1, root + B, B)
        v2s = []
        for lo in edges[:-1]:
            sel = (ks >= lo) & (ks < lo + B) & live
            if not sel.any():
                continue
            ww = float(su[sel].sum())
            if ww <= 0:
                continue
            ss = float(dd[sel].sum())
            v2s.append(ss * ss / ww)
        v2 = float(np.mean(v2s))
        say("  %-4d %-23.6f %-21.6f %.3f" % (B, v1, v2, PUB_R[B]))

    ros = float((dl.astype(np.float64) ** 2).sum() / sl.sum())
    order = np.argsort(-sl)
    share = float(sl[order[:5]].sum()) / float(sl.sum())
    d4 = abs(ros - PUB_ROS) <= 5e-3 and 0.30 <= share <= 0.36
    say()
    say("D4  ratio-of-sums baseline at B=1 = %.6f   published 1.78" % ros)
    say("    five largest-supp moduli carry %.4f of sum supp   %s"
        % (share, "hold" if d4 else "REFUTED"))

    # ------------------------------------------------------------- D5/D6
    rng = np.random.default_rng(SEED)
    say()
    say("D5/D6  lag-1 autocorrelation of D(k)/sqrt(supp(k))")
    say("  convention                     statistic   null sd    ratio"
        "    published")
    out = []
    for lab, x in (("surviving moduli only",
                    dl.astype(np.float64) / np.sqrt(sl)),
                   ("dead moduli kept at 0",
                    np.where(su > 0,
                             dd.astype(np.float64)
                             / np.sqrt(np.where(su > 0, su, 1)), 0.0))):
        r = lag1(x)
        nulls = np.empty(DRAWS)
        for i in range(DRAWS):
            nulls[i] = lag1(rng.permutation(x))
        sd = float(nulls.std(ddof=1))
        out.append((r, sd, r / sd))
        pub = PUB_AC if "surviving" in lab else PUB_AC0
        say("  %-30s %+.6f   %.6f   %+.3f   %+.4f / %.4f / %+.1f"
            % (lab, r, sd, r / sd, pub[0], pub[1], pub[2]))
    d5 = (abs(out[0][0] - PUB_AC[0]) <= 5e-5
          and abs(out[0][1] - PUB_AC[1]) <= 1e-3
          and abs(out[0][2] - PUB_AC[2]) <= 0.1)
    d6 = (abs(out[1][0] - PUB_AC0[0]) <= 5e-5
          and abs(out[1][1] - PUB_AC0[1]) <= 1e-3
          and abs(out[1][2] - PUB_AC0[2]) <= 0.1)
    say("  D5 %s   D6 %s" % ("hold" if d5 else "REFUTED",
                             "hold" if d6 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). Both point estimates reproduce exactly;")
    say("  only D6's null does not. 'Permutation null' does not say what")
    say("  is permuted when the dead moduli are kept: the whole vector,")
    say("  or only the surviving entries with the zeros held in place.")
    x0 = np.where(su > 0,
                  dd.astype(np.float64)
                  / np.sqrt(np.where(su > 0, su, 1)), 0.0)
    pos = np.flatnonzero(live)
    nulls2 = np.empty(DRAWS)
    for i in range(DRAWS):
        y = x0.copy()
        y[pos] = rng.permutation(x0[pos])
        nulls2[i] = lag1(y)
    say("    permuting the whole 9999-vector : sd %.6f" % out[1][1])
    say("    permuting survivors, zeros fixed: sd %.6f   published 0.0094"
        % float(nulls2.std(ddof=1)))
    say("    the point estimate is the same either way: %+.6f" % out[1][0])

    say()
    say("=" * 70)
    ok = d1 and d2 and d3 and d4 and d5 and d6 and d7
    say("D1 %s  D2 %s  D3 %s  D4 %s  D5 %s  D6 %s  D7 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (d1, d2, d3, d4, d5, d6, d7)))
    say("Section {#sec:R4} reproduces" if ok else "REFUTED")

    head = [
        "STATISTIC: (a) the complete-divisor-sum identity sum_k sum_m",
        "           mu(m)mu(N-mk) against mu(N-1); (b) the count of",
        "           surviving moduli and the gaps between them; (c) the",
        "           block ratio r(B) as the unweighted mean of",
        "           S_B(j)^2 / sum_{k in block} supp(k); (d) the",
        "           ratio-of-sums baseline and the share of sum supp",
        "           carried by the five largest moduli; (e),(f) the lag-1",
        "           autocorrelation of D(k)/sqrt(supp(k)) with and without",
        "           the dead moduli, each against a 400-draw permutation",
        "           null.",
        "FIELD: N = 10^8; k = 1 .. 9999, i.e. the whole band k < sqrt N on",
        "       which the type-II field is non-empty; m over",
        "       (10^4, N/k]; mu on [0,10^8] by a sqrt-sieve with a",
        "       cofactor array; permutation nulls with numpy default_rng",
        "       seed 20260808; (a) at N = 1000, 4096, 10^4, 30030, 65536,",
        "       10^5 with mu sieved to 2*10^5.",
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
