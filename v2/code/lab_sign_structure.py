# -*- coding: utf-8 -*-
r"""
Where the cross-k correlation of the dilated walls actually sits.

WHAT IS AT STAKE

Remark {#rem:nocrossk} measured the cross-k gain
G = sum(log k)|H| / |sum(log k)H| at 1.8 to 2.8 where independent
signs would give sqrt(#k) = 17.7 to 38.5, and Remark {#rem:mertens}
closed the obvious mechanism (a common Mertens factor).  Two
explanations remain and they are distinguishable:

  (a) HETEROGENEITY. If a few k carried almost all of sum(log k)|H|,
      the sum would fail to cancel for arithmetic reasons alone, with
      no correlation at all. Rule T4 tested a crude form of this and
      failed -- the top decile carries only about 35% -- but the right
      test is not a decile share, it is the exact random-sign null:
      keep the magnitudes |H(N;k)| and re-sign them at random.

  (b) CORRELATION. If the signs of the LARGE terms agree while the
      small ones are balanced, the count-weighted positive fraction can
      sit near 1/2 -- as it does, 0.38 to 0.48 -- while the
      MASS-weighted fraction is far from it. Algebraically, if a
      fraction f of sum(log k)|H| carries one sign then
      |sum(log k)H| = |2f-1| sum(log k)|H|, so G = 1/|2f-1| and
      G = 2.8 forces f near 0.68 or 0.32.

The random-sign null separates them exactly, and it is nearly free:
the magnitudes are already computed.

BACKS: Remark {#rem:signmass} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The random-sign null gives a gain G_null > 3 G_mu at every N --
      heterogeneity alone does not explain the missing cancellation.
  V2  The mass-weighted positive fraction f satisfies |2f - 1| >= 0.30
      at every N, i.e. f is at least 0.65 or at most 0.35.
  V3  The count-weighted positive fraction stays in [0.35, 0.55] at
      every N -- near balanced, so the asymmetry is in the weighting
      and not in the raw count.
  V4  The association is between magnitude and sign: splitting k at the
      median of |H|, the positive fraction differs between the two
      halves by more than 0.10.

REFUTATION RULE (fixed before the run)

  V1  REFUTED if G_null <= 3 G_mu at any N.
  V2  REFUTED if |2f - 1| < 0.30 at any N.
  V3  REFUTED if the count-weighted fraction leaves [0.35, 0.55].
  V4  REFUTED if the two halves differ by 0.10 or less at any N.

  All four gate.  V1 is the one that matters: it is the exact null for
  "is there any correlation at all".
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
OUT = os.path.join(ROOT, "results", "lab_sign_structure.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
THETA = 0.56
DRAWS = 200
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


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    rng = np.random.default_rng(SEED)

    say()
    say("  N          #k     G_mu     G_null (mean of %d)   ratio    "
        "mass frac +   count frac +   low/high halves" % DRAWS)
    say("  " + "-" * 100)
    v1 = v2 = v3 = v4 = True
    coin_diag = []
    for N in NS:
        v, PN, d = N, set(), 2
        while d * d <= v:
            if v % d == 0:
                PN.add(d)
                while v % d == 0:
                    v //= d
            d += 1
        if v > 1:
            PN.add(v)
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if mu[k] != 0 and all(k % q for q in PN)])
        lg = np.log(ks.astype(float))

        H = np.empty(ks.size)
        for i, k in enumerate(ks):
            k = int(k)
            m = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
            m = m[np.gcd(m, k) == 1]
            H[i] = float((lam[N - m * k] * mu[m]).sum())

        w = lg * np.abs(H)
        tot = float(w.sum())
        S = float((lg * H).sum())
        G_mu = tot / max(abs(S), 1e-300)

        gs = []
        for t in range(DRAWS):
            e = rng.integers(0, 2, size=H.size) * 2.0 - 1.0
            gs.append(tot / max(abs(float((e * w).sum())), 1e-300))
        G_null = float(np.mean(gs))

        pos = H > 0
        fmass = float(w[pos].sum() / tot)
        fcount = float(pos.mean())
        med = float(np.median(np.abs(H)))
        lowhalf = float(pos[np.abs(H) <= med].mean())
        highhalf = float(pos[np.abs(H) > med].mean())

        if G_null <= 3.0 * G_mu:
            v1 = False
        if abs(2 * fmass - 1) < 0.30:
            v2 = False
        if not (0.35 <= fcount <= 0.55):
            v3 = False
        if abs(highhalf - lowhalf) <= 0.10:
            v4 = False

        say("  %-10d %-6d %-8.3f %-21.2f %-8.1f %-13.4f %-14.4f "
            "%.4f / %.4f"
            % (N, ks.size, G_mu, G_null, G_null / G_mu, fmass, fcount,
               lowhalf, highhalf))
        arr = []
        for t in range(2):
            sig = np.zeros(N + 1, dtype=np.float64)
            supp = np.flatnonzero(mu[:N + 1] != 0)
            sig[supp] = rng.integers(0, 2, size=supp.size) * 2.0 - 1.0
            Ae = np.empty(ks.size)
            for i, k in enumerate(ks):
                k = int(k)
                m = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
                Ae[i] = float((lam[N - m * k] * sig[m * k]).sum())
            we = lg * np.abs(Ae)
            pe = Ae > 0
            mede = float(np.median(np.abs(Ae)))
            arr.append((float(we[pe].sum() / we.sum()), float(pe.mean()),
                        float(pe[np.abs(Ae) <= mede].mean()),
                        float(pe[np.abs(Ae) > mede].mean()),
                        float(we.sum() / max(abs(float((lg * Ae).sum())),
                                             1e-300))))
        coin_diag.append((N, arr))

    say()
    say("V1  G_null > 3 G_mu at every N                        %s"
        % ("hold" if v1 else "REFUTED"))
    say("V2  |2f_mass - 1| >= 0.30 at every N                  %s"
        % ("hold" if v2 else "REFUTED"))
    say("V3  count fraction in [0.35, 0.55] at every N         %s"
        % ("hold" if v3 else "REFUTED"))
    say("V4  halves differ by more than 0.10 at every N        %s"
        % ("hold" if v4 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc, and it is the control this finding")
    say("  needs). V1-V4 establish that the association exists; they do")
    say("  not say it is about mu. The random re-signing above holds the")
    say("  MAGNITUDES fixed, so it cannot answer that. Here the whole")
    say("  field is rebuilt from coin signs on the support of mu^2, and")
    say("  the same three statistics are read off it:")
    say("  N          mass frac +   count frac +   low/high halves   G")
    for N, arr in coin_diag:
        for (fm, fc, lo, hi, g) in arr:
            say("  %-10d %-13.4f %-14.4f %.4f / %-8.4f %.3f"
                % (N, fm, fc, lo, hi, g))

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("the missing cancellation is a correlation between magnitude and "
        "sign, not heterogeneity" if ok else "REFUTED")

    head = [
        "STATISTIC: the cross-k gain G = sum(log k)|H| / |sum(log k)H|",
        "           for mu; the same gain when the magnitudes |H(N;k)| are",
        "           kept and the signs are drawn at random; the fraction",
        "           of sum(log k)|H| carried by the k with H > 0 and the",
        "           fraction of k with H > 0; and that fraction split at",
        "           the median of |H|.",
        "NULL: the random re-signing IS the null, and it is exact for the",
        "      question asked -- it holds the magnitudes, and hence all",
        "      heterogeneity, fixed and destroys only the sign pattern.",
        "      200 draws per N.",
        "FIELD: N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6 with theta' = 0.56, so k",
        "       runs over the squarefree k < N^0.56 coprime to N; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an integer",
        "       sieve to 3.2e6; numpy default_rng seed 20260808.",
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
