# -*- coding: utf-8 -*-
r"""
Which m survive in H(N;k), and why the lean is deeper than counting
says.

WHAT IS AT STAKE

Remark {#rem:leanodd} established that the sign of H(N;k) is the sign
of the odd Mertens function at the length of the inner sum, and left
one thing open.  The model behind it -- a contributing k has one
surviving m, drawn uniformly from the admissible ones, so
f+ = (1 + Modd(x)/count(x))/2 -- tracks the shape but UNDER-predicts
the depth, by -0.1686 at N/k in [16,32) and -0.2699 at [64,128).  The
sign of each term is settled; how the surviving terms are selected is
not.

There is an obvious candidate.  A term survives when Lambda(N - mk) is
nonzero, and N - mk shrinks as m grows: at the top of the range it is
of order k, where primes are far denser than near N.  So the survivors
should be biased towards LARGE m -- and over a short range the large
odd squarefree m are mostly primes, with mu = -1.  That would deepen
the lean exactly as observed, and it would make m = 1, the one
guaranteed mu = +1, the least likely survivor of all, since N - k is
the largest value in the range.

The sharper form of the same idea replaces counting by the density
each term actually has: the standard sieve weight for N - mk to be
prime, w(m,k) = prod_{q odd, q | k not} [N != mk mod q] * q/(q-1),
giving the predictor P(N;k) = sum_m mu(m) w(m,k).

BACKS: Remark {#rem:survivors} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  m = 1 is under-represented: among sampled k with exactly one
      surviving m, the fraction whose survivor is m = 1 is below the
      uniform expectation, the mean of 1/count(N/k) over those k, at
      every N.
  R2  The survivors lean large: the mean of m/(N/k) over all survivors
      exceeds 1/2 at every N.
  R3  The sieve-weighted predictor beats counting: the agreement of
      sign H with sign P exceeds the agreement with sign Modd at
      every N.
  R4  And it is not an artefact of the marginals: that agreement also
      exceeds the largest of 16 permutation draws at every N.

REFUTATION RULE (fixed before the run)

  R1  REFUTED if the fraction is at or above the uniform expectation
      at any N. This is the direct test of the proposed mechanism.
  R2  REFUTED if the mean is at or below 1/2 at any N.
  R3  REFUTED if the sieve weight fails to improve on counting at any
      N, in which case density is not what selects the survivors.
  R4  REFUTED if it fails to beat every draw at any N.

  All four gate.

  THE CONTROL is R4's permutation, the design used in
  {#rem:leanmertens} and {#rem:leanodd}: the predictor's signs are
  shuffled among the distinct k, preserving the marginal sign
  distributions exactly.
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
OUT = os.path.join(ROOT, "results", "lab_survivor_selection.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
MLO, MHI = 2, 1000            # range of N/k sampled
SAMPLE = 60_000
QSIEVE = 30                   # sieve weight uses odd primes up to this
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

    MAJS = []          # the predictor's own majority share

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    odd = (np.arange(NMAX + 1) & 1).astype(np.int64)
    modd = np.cumsum(mu.astype(np.int64) * odd)
    del odd
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]

    say("  sieve weight over the odd primes %s"
        % ", ".join(map(str, QS)))

    res = []
    for N in NS:
        PN = factor_set(N)
        rng = np.random.default_rng(SEED + NS.index(N))
        lo = max(2, N // MHI)
        hi = N // MLO
        cand = np.array([k for k in range(lo, hi)
                         if sqf[k] and all(k % q for q in PN)],
                        dtype=np.int64)
        if cand.size > SAMPLE:
            cand = np.sort(rng.choice(cand, SAMPLE, replace=False))

        nsurv, one_is_1, one_tot = 0, 0, 0
        relpos, oneexp = [], []
        sh, sp, sm = [], [], []
        for k in cand:
            k = int(k)
            M = N // k
            if M < MLO:
                continue
            ms = np.arange(1, M + 1, 2, dtype=np.int64)
            ms = ms[sqf[ms]]
            for q in factor_set(k):
                if q > 2:
                    ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            lv = lam[N - ms * k]
            surv = ms[lv > 0]
            H = float((lv * mu[ms]).sum())
            # sieve-weighted predictor
            w = np.ones(ms.size, dtype=np.float64)
            for q in QS:
                if k % q == 0:
                    continue
                w *= np.where((N - ms * k) % q == 0, 0.0, q / (q - 1.0))
            P = float((mu[ms].astype(np.float64) * w).sum())
            if H != 0 and modd[M] != 0 and P != 0:
                sh.append(math.copysign(1.0, H))
                sp.append(math.copysign(1.0, P))
                sm.append(float(np.sign(modd[M])))
            if surv.size:
                nsurv += surv.size
                relpos.extend((surv / float(M)).tolist())
            if surv.size == 1:
                one_tot += 1
                one_is_1 += int(surv[0] == 1)
                # the uniform expectation must be taken over the SAME
                # k that the numerator counts, not over a separate
                # sweep: k with N/k near 2 have one admissible m and
                # would push the expectation to 1.
                oneexp.append(1.0 / ms.size)
        sh = np.array(sh)
        sp = np.array(sp)
        sm = np.array(sm)
        for _s in (sp, sm):
            MAJS.append(max(float((_s > 0).mean()),
                            float((_s < 0).mean())))
        res.append((N, cand.size, one_is_1, one_tot,
                    float(np.mean(relpos)) if relpos else float("nan"),
                    sh, sp, sm,
                    float(np.mean(oneexp)) if oneexp else float("nan")))
        say("  N = %-10d sampled k = %-8d single-survivor k = %d"
            % (N, cand.size, one_tot))

    say()
    say("R1  is m = 1 under-represented among single survivors?")
    say("  N            m=1 fraction   uniform expectation   ratio")
    r1 = True
    for N, nk, o1, ot, rp, sh, sp, sm, e in res:
        frac = o1 / ot if ot else float("nan")
        if frac >= e:
            r1 = False
        say("  %-12d %-14.4f %-21.4f %.4f" % (N, frac, e, frac / e))
    say("  R1 %s" % ("hold" if r1 else "REFUTED"))

    say()
    say("R2  do the survivors lean towards large m?")
    say("  N            mean m/(N/k)   uniform would give")
    r2 = all(r[4] > 0.5 for r in res)
    for N, nk, o1, ot, rp, sh, sp, sm, e in res:
        say("  %-12d %-14.4f %.4f" % (N, rp, 0.5))
    say("  R2 %s" % ("hold" if r2 else "REFUTED"))

    say()
    say("R3/R4  the sieve-weighted predictor")
    say("  N            sign P    sign Modd   draws max   #k")
    r3 = r4 = True
    rng = np.random.default_rng(SEED)
    for N, nk, o1, ot, rp, sh, sp, sm, e in res:
        ap = float((sh == sp).mean())
        am = float((sh == sm).mean())
        got = [float((sh == sp[rng.permutation(sp.size)]).mean())
               for _ in range(DRAWS)]
        mx = float(np.max(got))
        if ap <= am:
            r3 = False
        if ap <= mx:
            r4 = False
        say("  %-12d %-9.4f %-11.4f %-11.4f %d"
            % (N, ap, am, mx, sh.size))
    say("  R3 beats counting            %s" % ("hold" if r3 else "REFUTED"))
    say("  R4 beats every permutation   %s" % ("hold" if r4 else "REFUTED"))

    say()
    say("=" * 70)
    ok = r1 and r2 and r3 and r4
    say("the survivors are selected by prime density, which biases them "
        "to large m and deepens the lean" if ok else "REFUTED")

    mj = max(MAJS)
    say()
    say("  the predictor's own majority sign share, at its worst over "
        "everything")
    say("  reported above: %.4f. An agreement is only a measurement "
        "where the" % mj)
    say("  predictor has variance; where it takes one sign almost "
        "everywhere,")
    say("  the agreement is the other side's marginal rate read back.")
    say("MARGINAL lab_survivor_selection %.4f" % mj)
    if mj >= 0.9:
        say("DEGENERATE lab_survivor_selection")

    head = [
        "STATISTIC: among sampled k with 2 <= N/k <= 1000, the fraction of",
        "           single-survivor k whose survivor is m = 1 against the",
        "           uniform expectation; the mean relative position",
        "           m/(N/k) of survivors; and the agreement of sign H(N;k)",
        "           with the sign of the sieve-weighted predictor",
        "           P = sum_m mu(m) w(m,k), against the agreement with the",
        "           odd Mertens function and against 16 permutations.",
        "NULL: the permutation of R4, the design of [rem:leanodd]: the",
        "      predictor's signs are shuffled across the sampled k,",
        "      preserving the marginal sign distributions exactly, so the",
        "      baseline accounts for both sides being predominantly",
        "      negative.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree, coprime",
        "       to N, with 2 <= N/k <= 1000, subsampled to 60000 per N;",
        "       m odd squarefree, coprime to k, m <= N/k; the sieve weight",
        "       uses the odd primes up to 30; seed 20260808.",
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
