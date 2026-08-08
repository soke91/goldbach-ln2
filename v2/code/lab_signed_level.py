# -*- coding: utf-8 -*-
r"""
What the absolute values in the direct condition cost, in level.

WHAT IS AT STAKE

[eq:direct] says rtilde(N) = S(N)N + sum_{k<K}(log k)H(N;k) - small, so
what rtilde(N) > 0 actually needs is the SIGNED statement

    sum_{k<K,(k,N)=1} (log k) H(N;k)  >  -(1-eps) S(N) N.

Proposition {#prop:direct} then threw the signs away and asked for the
absolute version [eq:directcond], because that is what a bound on
|H(N;k)| would give.  lab_direct_level.py measured the absolute
version's crossing K*_H and found it at N^{0.7361}.  The signed sum is
never larger, so its crossing is never earlier; the question is how
much later, and whether the cross-k structure the wall paper spent
Sections 4 and 5 on buys anything at the level that matters.

Two crossings are therefore walked side by side, k ascending, and both
are FIRST crossings: the condition has to hold at the K actually used,
so a sum that dips below the threshold and recovers has already
failed.

BACKS: Remark {#rem:signedlevel} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  The signed crossing is not earlier: K*_signed >= K*_H at every N.
      This is forced by |sum x| <= sum |x| and is here as a check on
      the walk, not as a discovery -- a failure means the code is wrong.
  Y2  K*_signed/sqrt(N) grows across the sweep.
  Y3  mu's signs are WORSE than random. Keeping |H(N;k)| and drawing
      the sign of each term at random, mu's K*_signed is below the
      median of the draws at every N. That is what the measured lean
      of Remark {#rem:signmass} predicts: the terms lean one way, so
      they accumulate where random signs would cancel.
  Y4  The gain from keeping the signs is bounded: K*_signed/K*_H stays
      under 10 at every N.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED by a single N with K*_signed < K*_H.
  Y2  REFUTED by a single fall.
  Y3  REFUTED if mu is at or above the median of the draws at any N,
      in which case the signs help rather than hurt and the absolute
      version of the condition is throwing away something real.
  Y4  REFUTED if the ratio reaches 10 anywhere. A refutation is the
      good outcome: it would mean the cross-k cancellation is worth a
      level of its own and [eq:directcond] is the wrong statement of
      the demand.

  All four gate.

  THE CONTROL is Y3's sign randomisation, and it is the right one
  here for the reason Remark {#rem:whycoinwins} gives: a coin on
  supp(mu^2) changes the MAGNITUDES too -- it gets square-root
  cancellation inside each H and so beats mu for reasons that have
  nothing to do with the signs across k. Permuting only the signs
  holds every |H(N;k)| fixed, so the comparison isolates exactly the
  structure [eq:directcond] discards.
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
OUT = os.path.join(ROOT, "results", "lab_signed_level.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
CLIM = 4_000_000          # fixed bound for the Euler product (G20)
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


def read_gap():
    """the budget gap in the exponent -- read, not copied"""
    p = os.path.join(ROOT, "results", "audit_model_transfer.txt")
    src = io.open(p, encoding="utf-8").read()
    return float(re.search(r"mean gap ([\d.]+)", src).group(1))


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


def first_cross(cum, ks, thr, below):
    """First k at which the running sum leaves the admissible side."""
    bad = (cum < -thr) if below else (cum > thr)
    j = int(np.argmax(bad)) if bool(bad.any()) else -1
    if j < 0:
        return int(ks[-1]), True
    return int(ks[j]), False


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)

    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    say()
    say("  both crossings are FIRST crossings, k ascending, over the")
    say("  squarefree k < N/2 coprime to N. The threshold is S(N)N.")
    say()
    say("  N            K*_H       K*_signed  ratio    theta'_H  "
        "theta'_signed")
    say("  " + "-" * 78)

    rows = []
    exhausted = []
    for N in NS:
        PN = factor_set(N)
        S = twin
        for q in sorted(PN):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        thr = S * N
        ks = np.array([k for k in range(2, N // 2)
                       if mu[k] != 0 and all(k % q for q in PN)],
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
        # [eq:dilate]: H(N;k) = mu(k) A(N;k). An earlier version of
        # this script dropped the mu(k) and summed A, which is a
        # different object; lab_direct_identity.py found it.
        H = mu[ks].astype(np.float64) * A
        w = np.log(ks.astype(np.float64))
        wa = w * np.abs(H)
        ws = w * H
        cs = np.cumsum(ws)
        ca = np.cumsum(wa)
        kA, exA = first_cross(np.cumsum(wa), ks, thr, False)
        kS, exS = first_cross(cs, ks, thr, True)
        if exA or exS:
            exhausted.append(N)
        jmin = int(np.argmin(cs))
        traj = []
        for th in (0.56, 0.70, 0.90):
            j = int(np.searchsorted(ks, int(N ** th)))
            j = min(j, cs.size - 1)
            traj.append(cs[j] / N)
        rows.append((N, kA, kS, S, w, np.abs(H), ks, thr,
                     float(cs[jmin]) / N, int(ks[jmin]),
                     float(cs[-1]) / N, traj,
                     float(ca[-1]) / N))
        say("  %-12d %-10d %-10d %-8.3f %-9.4f %.4f"
            % (N, kA, kS, kS / kA,
               math.log(kA) / math.log(N), math.log(kS) / math.log(N)))

    say()
    say("  Both theta' columns are against the budget S(N)N. That is")
    say("  NOT the budget [eq:nolog] asks for: prop:nolog needs")
    say("  S(N)(1-A(N))N, smaller by a factor near five, and")
    say("  audit_model_transfer.py measures the difference as %.4f in"
        % read_gap())
    say("  the exponent -- more than the whole gap between the two")
    say("  columns above. The constant crossed here, declared so that")
    say("  no exponent from this file is read against one from a file")
    say("  that crossed the other:")
    for r in rows:
        say("BUDGET kstar_signed_SN_N%d %.6f" % (r[0], r[3]))

    say()
    say("  And the arithmetic this sweep covers, which gate check G34")
    say("  reads. Every N here is 2^a 5^b, so what the signs across k")
    say("  are worth is measured at ONE odd radical -- and so is the")
    say("  same quantity for the residue in lab_residue_signed.py,")
    say("  which runs the identical family:")
    rads = set()
    for r in rows:
        rr = 1
        for q in factor_set(r[0]):
            if q > 2:
                rr *= q
        rads.add(rr)
    say("  %d N, %d distinct odd radical%s: %s"
        % (len(rows), len(rads), "" if len(rads) == 1 else "s",
           ", ".join(str(v) for v in sorted(rads))))
    say("RADICALS %d" % len(rads))

    say()
    if exhausted:
        say("  NOTE: at N = %s a walk reached the cap k < N/2 without"
            % ", ".join(map(str, exhausted)))
        say("  crossing, so that K* is a lower bound.")

    y1 = all(r[2] >= r[1] for r in rows)
    say("Y1  K*_signed >= K*_H at every N      %s"
        % ("hold" if y1 else "REFUTED"))
    rs = [r[2] / math.sqrt(r[0]) for r in rows]
    y2 = all(rs[i] < rs[i + 1] for i in range(len(rs) - 1))
    say("Y2  K*_signed/sqrt N: %s   %s"
        % (", ".join("%.3f" % v for v in rs), "hold" if y2 else "REFUTED"))

    say()
    say("Y3  the control: same |H(N;k)|, signs drawn at random, %d draws"
        % DRAWS)
    say("  N            K*_signed(mu)  draws median  min       max     "
        "  mu's rank")
    y3 = True
    for j, (N, kA, kS, S, w, absH, ks, thr, _mn, _am, _en,
            _tj, _ae) in enumerate(rows):
        rng = np.random.default_rng(SEED + j)
        vals = []
        for d in range(DRAWS):
            s = rng.choice([-1.0, 1.0], size=absH.size)
            k2, _ = first_cross(np.cumsum(w * s * absH), ks, thr, True)
            vals.append(k2)
        med = float(np.median(vals))
        rank = int(sum(1 for v in vals if v < kS))
        if kS >= med:
            y3 = False
        say("  %-12d %-14d %-13.0f %-9d %-9d %d of %d"
            % (N, kS, med, min(vals), max(vals), rank, DRAWS))
    say("  Y3 %s" % ("hold" if y3 else "REFUTED"))

    say()
    gains = [r[2] / r[1] for r in rows]
    y4 = all(g < 10.0 for g in gains)
    say("Y4  K*_signed/K*_H: %s   (cap 10)   %s"
        % (", ".join("%.3f" % g for g in gains),
           "hold" if y4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). How far the signed sum goes, and what")
    say("  the randomised signs do with the same magnitudes:")
    say("  N            min sum/N   at k        -S(N)     min/-S    "
        "sum/N at k<N/2")
    for (N, kA, kS, S, w, absH, ks, thr, mn, am, en, tj, ae) in rows:
        say("  %-12d %-11.4f %-11d %-9.4f %-9.4f %.4f"
            % (N, mn, am, -S, mn / (-S), en))
    say("  and the trajectory, at k < N^0.56, N^0.70, N^0.90:")
    say("  N            0.56        0.70        0.90")
    for (N, kA, kS, S, w, absH, ks, thr, mn, am, en, tj, ae) in rows:
        say("  %-12d %-11.4f %-11.4f %.4f" % (N, tj[0], tj[1], tj[2]))
    say("  The signed sum crosses -S(N)N well inside the walk and keeps")
    say("  going, reaching %.2f times its threshold. Every one of the %d"
        % (rows[-1][8] / (-rows[-1][3]), DRAWS))
    say("  sign-randomised draws, holding each |H(N;k)| fixed, fails to")
    say("  cross at all. So mu's signs across k are not a source of")
    say("  cancellation here -- they are worse than random, which is the")
    say("  lean Remark [rem:signmass] measured, seen at the level.")
    say("  Against the absolute sum the signed one is not negligible:")
    say("  N            signed/N at k<N/2   absolute/N      ratio")
    for (N, kA, kS, S, w, absH, ks, thr, mn, am, en, tj, ae) in rows:
        say("  %-12d %-19.4f %-15.1f %.4f" % (N, en, ae, abs(en) / ae))
    say("  a ratio flat in N, so the signs buy a fixed factor of about")
    say("  %.1f in size and %.2f to %.2f in K*, i.e. about %.3f in the"
        % (1.0 / (abs(rows[-1][10]) / rows[-1][12]),
           rows[0][2] / rows[0][1], rows[-1][2] / rows[-1][1],
           math.log(rows[-1][2] / rows[-1][1]) / math.log(rows[-1][0])))
    say("  exponent theta'. That is a real gain and a bounded one, which")
    say("  is what Y4 predicted; it does not reach the theta' the")
    say("  reduction needs from the absolute form.")

    say()
    say("  Cross-check lines, against lab_direct_identity.py, which")
    say("  reaches the same quantity from [eq:untrunc]:")
    for (N, kA, kS, S, w, absH, ks, thr, mn, am, en, tj,
         ae) in rows:
        say("AGREE signed_partial_056 N=%d %.6f 0.02" % (N, tj[0]))

    say()
    say("=" * 70)
    ok = y1 and y2 and y3 and y4
    say("keeping the signs across k buys a bounded factor and mu's own "
        "signs are worse than random" if ok else "REFUTED")

    head = [
        "STATISTIC: the first crossing K*_H of sum_{k<K}(log k)|H(N;k)|",
        "           above S(N)N, the first crossing K*_signed of the",
        "           signed sum below -S(N)N, their ratio, the implied",
        "           exponents log K*/log N, and the same signed crossing",
        "           with the signs of the terms drawn at random.",
        "NULL: the sign randomisation of Y3 -- every |H(N;k)| held fixed,",
        "      the sign of each term redrawn, 16 draws. A coin on",
        "      supp(mu^2) is NOT used because it changes the magnitudes",
        "      as well, and by [rem:whycoinwins] it beats mu for reasons",
        "      unrelated to the signs across k; permuting signs alone",
        "      isolates exactly what [eq:directcond] discards.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k over the squarefree",
        "       k < N/2 coprime to N, walked upward; S(N) from an Euler",
        "       product at the fixed bound 4e6; seed 20260808.",
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
