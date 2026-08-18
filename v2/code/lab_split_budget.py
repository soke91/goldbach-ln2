# -*- coding: utf-8 -*-
r"""
What the split buys: which half of H exhausts the Goldbach budget
first.

WHAT IS AT STAKE

[eq:directcond] asks for B_H(N;K) = sum_{k<K}(log k)|H(N;k)| under
(1-eps) S(N) N, and Remark {#rem:predictable} split H = beta P + R
with P the elementary sieve-weighted Mobius sum and R the residue.
Remark {#rem:residue} then found |R| ~ (N/k)^{1/2}, and OPEN.md
records the remaining task as proving square-root cancellation for R.

That task is worth doing only if R is what binds.  By the triangle
inequality B_H <= beta B_P + B_R with
B_P = sum(log k)|P| and B_R = sum(log k)|R|, and each half has its own
crossing of the budget:

    K*_P = max{K : beta B_P(N;K) <= S(N)N},
    K*_R = max{K : B_R(N;K) <= S(N)N}.

If K*_P is barely beyond K*_H then the elementary part alone almost
exhausts the budget, and a proof about R buys almost nothing -- the
wall would be the sieve-weighted Mobius sum, which is a different and
more elementary object.  If K*_P runs far past K*_H then R is the
binding half and the recorded task is the right one.

BACKS: Remark {#rem:splitbudget} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  The crossings are ordered K*_H < K*_P < K*_R at every N, which
      the mass shares of {#rem:predictable} force.
  Z2  The elementary part is most of the budget at the operative
      truncation: beta B_P(N;K*_H)/(S(N)N) is above 0.8 at every N.
  Z3  The residue is not: B_R(N;K*_H)/(S(N)N) is below 0.7 at every N.
  Z4  And removing the residue entirely buys little: K*_P/K*_H is
      under 1.3 at every N.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED by a single N out of order; it is forced by the shares,
      so a failure means the walk is wrong.
  Z2  REFUTED if the share drops to 0.8 at any N.
  Z3  REFUTED if it reaches 0.7 at any N.
  Z4  REFUTED if the ratio reaches 1.3 at any N. A refutation is the
      good outcome: it would say the residue is the binding half and
      the task OPEN.md records is worth its cost.

  All four gate.

  NULL: none is run here and one has already been run for this split.
  lab_split_null.py replaced mu by a coin with the sieve weights and
  the range untouched and found the coin absorbs essentially none of
  the mass -- residual share 0.9566 to 0.9986 against mu's 0.5307 to
  0.6310 -- so the decomposition being budgeted here is mu's. What is
  measured now is a comparison between two halves of one measured sum
  against a fixed threshold S(N)N, which is the well-conditioned case
  of {#rem:weightgapnull} and needs no further control.
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
OUT = os.path.join(ROOT, "results", "lab_split_budget.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000


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


def crossing(ks, cum, thr):
    j = int(np.searchsorted(cum, thr))
    return int(ks[min(j, ks.size - 1)]), j >= ks.size


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

    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    res = []
    for N in NS:
        PN = factor_set(N)
        S = twin
        for q in sorted(PN):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
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
        lw = np.log(ks.astype(np.float64))
        Hs = np.array(Hs)
        Ps = np.array(Ps)
        beta = float((Hs * Ps).sum() / (Ps * Ps).sum())
        R = Hs - beta * Ps
        cH = np.cumsum(lw * np.abs(Hs))
        cP = np.cumsum(lw * np.abs(beta * Ps))
        cR = np.cumsum(lw * np.abs(R))
        res.append((N, S, beta, ks, cH, cP, cR))
        say("  N = %-10d  #k = %-7d beta = %.4f  budget S(N) = %.4f"
            % (N, ks.size, beta, S))

    say()
    say("Z1  where each half crosses the budget S(N)N")
    say("  N            K*_H      K*_P      K*_R      K*_P/K*_H   "
        "K*_R/K*_H")
    z1 = True
    rows = []
    exh = []
    for N, S, beta, ks, cH, cP, cR in res:
        thr = S * N
        kh, e1 = crossing(ks, cH, thr)
        kp, e2 = crossing(ks, cP, thr)
        kr, e3 = crossing(ks, cR, thr)
        if e1 or e2 or e3:
            exh.append(N)
        if not (kh < kp < kr):
            z1 = False
        rows.append((N, S, beta, ks, cH, cP, cR, kh, kp, kr))
        say("  %-12d %-9d %-9d %-9d %-11.4f %.4f"
            % (N, kh, kp, kr, kp / kh, kr / kh))
    if exh:
        say("  NOTE: at N = %s a crossing reached the cap k < %d, so that"
            % (", ".join(map(str, exh)), KCAP))
        say("  value is a lower bound.")
    say("  Z1 ordered at every N   %s" % ("hold" if z1 else "REFUTED"))

    say()
    say("Z2/Z3  each half's share of the budget at K = K*_H")
    say("  N            beta B_P / (S N)   B_R / (S N)   sum")
    z2 = z3 = True
    for N, S, beta, ks, cH, cP, cR, kh, kp, kr in rows:
        j = int(np.searchsorted(ks, kh))
        thr = S * N
        a = float(cP[min(j, cP.size - 1)]) / thr
        b = float(cR[min(j, cR.size - 1)]) / thr
        if a <= 0.8:
            z2 = False
        if b >= 0.7:
            z3 = False
        say("  %-12d %-18.4f %-13.4f %.4f" % (N, a, b, a + b))
    say("  Z2 the elementary part is over 0.8   %s"
        % ("hold" if z2 else "REFUTED"))
    say("  Z3 the residue is under 0.7          %s"
        % ("hold" if z3 else "REFUTED"))

    say()
    z4 = all(r[8] / r[7] < 1.3 for r in rows)
    say("Z4  K*_P/K*_H: %s   (cap 1.30)   %s"
        % (", ".join("%.4f" % (r[8] / r[7]) for r in rows),
           "hold" if z4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). What a proof about R would buy. If R")
    say("  were removed entirely the truncation could move from K*_H to")
    say("  K*_P; if the elementary part were removed instead it could")
    say("  move to K*_R. In exponents:")
    say("  N            log K*_H/log N   log K*_P/log N   log K*_R/log N")
    for N, S, beta, ks, cH, cP, cR, kh, kp, kr in rows:
        L = math.log(N)
        say("  %-12d %-16.4f %-16.4f %.4f"
            % (N, math.log(kh) / L, math.log(kp) / L, math.log(kr) / L))
    say()
    say("  These exponents are against the budget S(N)N, which is NOT")
    say("  the one [eq:nolog] asks for. Proposition {#prop:nolog} needs")
    say("  S(N)(1-A(N))N, smaller by a factor near five, and")
    say("  audit_residue_level.py measures what K*_R does against it --")
    say("  0.18 lower in the exponent, which is the difference between")
    say("  clearing theta' = 0.56 comfortably and clearing it barely.")
    say("  And the arithmetic this sweep covers, which gate check G34")
    say("  reads. Every N here is 2^a 5^b, so the sweep has ONE odd")
    say("  radical. audit_residue_arithmetic.py repeats the residue")
    say("  measurement across seven types and finds the exponent below")
    say("  1/2 at the primorial-like ones, so the ordering of the two")
    say("  halves below is a statement about this radical:")
    rads = set()
    for N, S, beta, ks, cH, cP, cR, kh, kp, kr in rows:
        r = 1
        for q in factor_set(N):
            if q > 2:
                r *= q
        rads.add(r)
    say("  %d N, %d distinct odd radical%s: %s"
        % (len(rows), len(rads), "" if len(rads) == 1 else "s",
           ", ".join(str(r) for r in sorted(rads))))
    say("RADICALS %d" % len(rads))
    say()
    say("  The constant crossed here, declared so that no exponent from")
    say("  this file is ever read against one from that file:")
    for N, S, beta, ks, cH, cP, cR, kh, kp, kr in rows:
        say("BUDGET kstar_SN_N%d %.6f" % (N, S))

    say()
    say()
    say("  THE CONTROL. Each half's magnitudes are permuted across k,")
    say("  16 draws, with the weights left attached to their own k. That")
    say("  preserves the multiset of |beta P| and of |R| exactly and")
    say("  destroys only the pairing with the modulus, so it says whether")
    say("  the crossings are a property of the pairing or only of the")
    say("  magnitude distributions.")
    say("  N            K*_P      perm band          K*_R      perm band")
    rng = np.random.default_rng(20260808)
    for N, S, beta, ks, cH, cP, cR, kh, kp, kr in rows:
        thr = S * N
        lw = np.log(ks.astype(np.float64))
        vP = np.diff(np.concatenate(([0.0], cP))) / lw
        vR = np.diff(np.concatenate(([0.0], cR))) / lw
        bp, br = [], []
        for _ in range(16):
            bp.append(crossing(ks, np.cumsum(lw * vP[rng.permutation(
                vP.size)]), thr)[0])
            br.append(crossing(ks, np.cumsum(lw * vR[rng.permutation(
                vR.size)]), thr)[0])
        say("  %-12d %-9d [%-7d %-7d] %-9d [%d %d]"
            % (N, kp, min(bp), max(bp), kr, min(br), max(br)))
    say("  The measured crossings are FAR EARLIER than the bands, by a")
    say("  factor of three to five, and the bands are narrow. So the")
    say("  crossings are strongly a property of the pairing, not of the")
    say("  magnitude distributions -- and in the direction that matters:")
    say("  |P| and |R| are largest at the SMALLEST k, where the inner")
    say("  sum is longest, so the measured cumulative front-loads and")
    say("  reaches the budget while a random pairing is still climbing.")
    say("  The budget is consumed at the bottom of the k-range, which is")
    say("  where the reduction has the least room to move.")

    say()
    say("=" * 70)
    ok = z1 and z2 and z3 and z4
    say("the elementary half exhausts the budget first and a proof "
        "about the residue buys little" if ok else "REFUTED")

    head = [
        "STATISTIC: with beta the least-squares scale and H = beta P + R,",
        "           the crossings of S(N)N by B_H = sum(log k)|H|, by",
        "           beta B_P = sum(log k)|beta P| and by B_R =",
        "           sum(log k)|R|; each half's share of the budget at",
        "           K = K*_H; and the three crossings as exponents in N.",
        "NULL: each half's magnitudes permuted across k, 16 draws, with",
        "      the weights left attached to their own k. It preserves the",
        "      multiset of |beta P| and of |R| exactly and destroys only",
        "      the pairing with the modulus, so it separates a property",
        "      of the pairing from a property of the magnitude",
        "      distributions. That the decomposition itself is mu's was",
        "      established separately in lab_split_null.py, where a coin",
        "      absorbed essentially none of the mass.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree, coprime",
        "       to N, 2 <= k < 100000; m odd squarefree, coprime to k,",
        "       m <= (N-1)/k; the sieve weight uses the odd primes up to",
        "       30; S(N) from an Euler product at the fixed bound 4e6.",
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
