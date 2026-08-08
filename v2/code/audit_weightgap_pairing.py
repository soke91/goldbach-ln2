# -*- coding: utf-8 -*-
r"""
The control that IS well conditioned for the weight-gap ratio.

WHAT IS AT STAKE

Remark {#rem:weightgapnull} ran a coin on lab_weight_gap.py's
statistics and found it unusable for two of them: a coin drives both
|sum H| and |sum (log k)H| to square-root size, so their ratio is a
quotient of two near-zero quantities and takes any value.  Those two
claims -- the ratio itself, and the effective modulus read off the
profile in j -- were left standing uncontrolled, and {#rem:weightgap}
still carries them.

There is a control that keeps them well conditioned.  Permute the
values H(N;k) across k, leaving the weights (log k)^j attached to
their own k.  A permutation does not change a plain sum, so sum H is
EXACTLY invariant: the numerator of the ratio is pinned, and only the
weighted denominator moves.  What the permutation destroys is the
pairing between the weight and the value -- which is the only thing
the claims can be about.

The suspicion it tests is concrete.  If H and log k are unpaired then
sum (log k)^j H is about mean((log k)^j) times sum H, so the ratio is
about 1/mean(log k) and the profile in j is automatically near
geometric with an "effective modulus" equal to the geometric mean of
the k-range.  If that is what the numbers are, the reading in
{#rem:weightgap} says nothing about concentration.

BACKS: Remark {#rem:pairingnull} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The control is well conditioned: sum H is invariant under the
      permutation to machine precision, and the permuted ratio's
      spread over 16 draws is under 0.05 of its median.
  W2  The ratio is not about the pairing: mu's |sum H|/|sum (log k)H|
      lies inside the range spanned by the 16 permutations, at
      every N.
  W3  Nor is the profile: mu's spread of the four consecutive ratios
      in j lies inside the permutations' range, at every N.
  W4  And the effective modulus is just the geometric mean of the
      k-range: |log k* - mean(log k)| / log K is under 0.05 at
      every N.

REFUTATION RULE (fixed before the run)

  W1  REFUTED if sum H moves by 1e-12 relative, or if the spread
      reaches 0.05 of the median. The first would mean the
      permutation is misapplied; the second that this control is no
      better conditioned than the coin.
  W2  REFUTED if mu falls outside the permutations' range at any N --
      the good outcome, since it would mean the ratio does carry
      pairing information.
  W3  Likewise for the profile.
  W4  REFUTED if the gap reaches 0.05 of log K at any N, which would
      mean k* is not simply the geometric mean.

  All four gate.

  NULL: the permutation of H across k, weights left in place, 16
  draws. It preserves the multiset of values exactly -- and therefore
  sum H, sum |H| and every moment of |H| -- and destroys only the
  correspondence between a modulus and its dilated wall.
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
OUT = os.path.join(ROOT, "results", "audit_weightgap_pairing.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
JS = [0.0, 0.25, 0.5, 0.75, 1.0]
THETA = 0.56
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


def stats(lg, H):
    """the ratio and the profile spread, as lab_weight_gap.py forms them"""
    row = [abs(float(((lg ** j) * H).sum())) for j in JS]
    ratio = row[0] / max(row[-1], 1e-300)
    rr = [row[i + 1] / row[i] for i in range(len(row) - 1) if row[i] > 0]
    if len(rr) != len(JS) - 1:
        return ratio, float("inf"), float("nan")
    spread = max(rr) / min(rr)
    logk = float(np.mean(rr)) ** 4.0
    return ratio, spread, logk


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0

    res = []
    for N in NS:
        PN = factor_set(N)
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if sqf[k] and all(k % q for q in PN)])
        lg = np.log(ks.astype(float))
        f0 = np.zeros(N, dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        f0[1:] = lam[1:N] * mu[N - idx].astype(np.float64)
        A = np.empty(ks.size)
        for i, k in enumerate(ks):
            r = N % int(k)
            A[i] = f0[r::int(k)].sum() if r else f0[int(k)::int(k)].sum()
        del f0
        H = mu[ks].astype(np.float64) * A
        res.append((N, K, ks, lg, H))
        say("  N = %-10d  K = %-7d #k = %d" % (N, K, ks.size))

    say()
    say("W1  the control is well conditioned")
    say("  N            sum H drift   permuted ratio spread / median")
    w1 = True
    perm = {}
    for N, K, ks, lg, H in res:
        rng = np.random.default_rng(SEED + NS.index(N))
        base = float(H.sum())
        got, drift = [], 0.0
        for _ in range(DRAWS):
            Hp = H[rng.permutation(H.size)]
            drift = max(drift, abs(float(Hp.sum()) - base)
                        / max(abs(base), 1e-300))
            got.append(stats(lg, Hp))
        perm[N] = got
        r = [g[0] for g in got]
        rel = (max(r) - min(r)) / float(np.median(r))
        if drift >= 1e-12 or rel >= 0.05:
            w1 = False
        say("  %-12d %-13.3e %.4f" % (N, drift, rel))
    say("  W1 %s" % ("hold" if w1 else "REFUTED"))

    say()
    say("W2  the ratio |sum H| / |sum (log k)H|")
    say("  N            mu        perm min   perm median   perm max")
    w2 = True
    for N, K, ks, lg, H in res:
        m = stats(lg, H)[0]
        p = [g[0] for g in perm[N]]
        if not (min(p) <= m <= max(p)):
            w2 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, m, min(p), float(np.median(p)), max(p)))
    say("  W2 mu inside the permutations' range   %s"
        % ("hold" if w2 else "REFUTED"))

    say()
    say("W3  spread of the four consecutive ratios in j")
    say("  N            mu        perm min   perm median   perm max")
    w3 = True
    for N, K, ks, lg, H in res:
        m = stats(lg, H)[1]
        p = [g[1] for g in perm[N]]
        if not (min(p) <= m <= max(p)):
            w3 = False
        say("  %-12d %-9.4f %-10.4f %-13.4f %.4f"
            % (N, m, min(p), float(np.median(p)), max(p)))
    say("  W3 %s" % ("hold" if w3 else "REFUTED"))

    say()
    say("W4  what the effective modulus actually is")
    say("  N            log k*    mean log k   |gap|/log K   k*/K")
    w4 = True
    for N, K, ks, lg, H in res:
        lk = stats(lg, H)[2]
        ml = float(lg.mean())
        g = abs(lk - ml) / math.log(K)
        if g >= 0.05:
            w4 = False
        say("  %-12d %-9.4f %-12.4f %-13.4f %.4f"
            % (N, lk, ml, g, math.exp(lk) / K))
    say("  W4 k* is the geometric mean of the k-range   %s"
        % ("hold" if w4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). What an unpaired sum would give. If")
    say("  H and log k are independent then sum (log k)^j H is about")
    say("  mean((log k)^j) times sum H, so the ratio is 1/mean(log k)")
    say("  and the profile is geometric by construction:")
    say("  N            1/mean(log k)   mu's ratio   permuted median")
    for N, K, ks, lg, H in res:
        p = [g[0] for g in perm[N]]
        say("  %-12d %-15.4f %-12.4f %.4f"
            % (N, 1.0 / float(lg.mean()), stats(lg, H)[0],
               float(np.median(p))))

    say()
    say("=" * 70)
    ok = w1 and w2 and w3 and w4
    say("the ratio and the effective modulus carry no pairing "
        "information" if ok else "REFUTED")

    head = [
        "STATISTIC: |sum H|/|sum (log k)H|, the spread of the four",
        "           consecutive ratios of |sum (log k)^j H| in j, and the",
        "           effective modulus log k* read off that profile, for mu",
        "           and for 16 permutations of H across k; together with",
        "           the mean of log k over the same range.",
        "NULL: the permutation of H(N;k) across k with the weights left",
        "      attached to their own k, 16 draws. It preserves the",
        "      multiset of values exactly -- hence sum H, sum |H| and",
        "      every moment of |H| -- and destroys only the pairing",
        "      between a modulus and its dilated wall. Unlike the coin of",
        "      [rem:weightgapnull] it leaves the ratio's numerator pinned,",
        "      which is what makes it usable here.",
        "FIELD: N = 2e5 through 3.2e6 by doubling, theta' = 0.56; k over",
        "       the squarefree k < N^0.56 coprime to N; seed 20260808.",
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
