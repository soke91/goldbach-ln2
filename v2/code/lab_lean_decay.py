# -*- coding: utf-8 -*-
r"""
Is the negative lean of the large dilated walls a finite-N effect?

WHAT IS AT STAKE

Remark {#rem:signmass} localised the missing cross-k cancellation: the
mass-weighted fraction of k with H(N;k) > 0 is far below one half, so
the large dilated walls lean negative, and that is why E_3 is negative
at every N measured.  It left "why" open.  But its own table already
carries a hint it did not read: the fraction runs

    0.2273, 0.2228, 0.2735, 0.3068, 0.3207

at N = 2e5 to 3.2e6 -- rising towards one half.  If that is the whole
story then the lean is not a structural fact about mu at all but the
finite-N error of Theorem {#thm:C}, which Remark {#rem:thetasweep}
measured at 0.17 to 0.46 of N and which decays; and the "no cross-k
cancellation" of Remark {#rem:nocrossk} is a statement about the
accessible range and not an asymptotic one, since

    G = 1 / |2f - 1|  ->  infinity  as  f -> 1/2.

That would change what the finding means, so it is tested here with a
sixth N added to lengthen the lever.

BACKS: Remark {#rem:leandecay} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The mass fraction rises: f(6.4e6) > f(2e5) + 0.05.
  W2  |0.5 - f| fitted as N^{-b} gives b in (0.05, 0.35).
  W3  G(N) = 1/|2f-1| rises with N: G(6.4e6) > G(2e5).
  W4  The decay is slow: extrapolating the fit, the N at which
      |0.5 - f| falls below 0.01 exceeds 10^10.

REFUTATION RULE (fixed before the run)

  W1  REFUTED if the rise is 0.05 or less.  If the fraction is flat,
      the lean is structural and Remark {#rem:nocrossk} stands as an
      asymptotic claim.
  W2  REFUTED if b leaves (0.05, 0.35).
  W3  REFUTED if G does not rise.
  W4  REFUTED if the extrapolated N is 10^10 or below.

  All four gate.
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
OUT = os.path.join(ROOT, "results", "lab_lean_decay.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
DRAWS = 2
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


def loo(x, y, name, say):
    """Refit dropping each end in turn, and report the spread.

    Every exponent this repository quotes is a slope over four to eight
    values of N, and audit_truncation_exponent.py showed what such a
    slope is worth when nobody varies the free parameter that defines
    it. For a direct fit the free parameter is the N-range, so the
    cheapest honest check is to refit without the smallest N and
    without the largest and print how far the answer moves.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    sp = max(f) - min(f)
    say("  leave-one-out on %s: full %.4f, without the smallest N "
        "%.4f," % (name, f[0], f[1]))
    say("  without the largest %.4f -- spread %.4f" % (f[2], sp))
    say("SWEPT %s N-range %.4f" % (name, sp))
    return sp


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
    say("  N          #k     mass frac +   |0.5 - f|   G = 1/|2f-1|   "
        "coin frac + (mean of %d)" % DRAWS)
    say("  " + "-" * 88)
    fs, gs = [], []
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
        f = float(w[H > 0].sum() / w.sum())
        G = 1.0 / max(abs(2 * f - 1), 1e-300)
        fs.append(f)
        gs.append(G)

        cf = []
        for t in range(DRAWS):
            sig = np.zeros(N + 1, dtype=np.float64)
            supp = np.flatnonzero(mu[:N + 1] != 0)
            sig[supp] = rng.integers(0, 2, size=supp.size) * 2.0 - 1.0
            Ae = np.empty(ks.size)
            for i, k in enumerate(ks):
                k = int(k)
                m = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
                Ae[i] = float((lam[N - m * k] * sig[m * k]).sum())
            we = lg * np.abs(Ae)
            cf.append(float(we[Ae > 0].sum() / we.sum()))
        say("  %-10d %-6d %-13.4f %-11.4f %-15.3f %.4f"
            % (N, ks.size, f, abs(0.5 - f), G, float(np.mean(cf))))

    say()
    w1 = fs[-1] > fs[0] + 0.05
    say("W1  f rises by %.4f from %d to %d   (floor 0.05)   %s"
        % (fs[-1] - fs[0], NS[0], NS[-1], "hold" if w1 else "REFUTED"))

    dev = np.array([abs(0.5 - f) for f in fs])
    b = -float(np.polyfit(np.log(np.array(NS, dtype=float)),
                          np.log(dev), 1)[0])
    w2 = 0.05 < b < 0.35
    say("W2  |0.5 - f| ~ N^{-%.4f}   (band (0.05, 0.35))   %s"
        % (b, "hold" if w2 else "REFUTED"))
    loo(np.log(np.array(NS, dtype=float)), np.log(dev),
        "lean_decay", say)

    w3 = gs[-1] > gs[0]
    say("W3  G rises from %.3f to %.3f   %s"
        % (gs[0], gs[-1], "hold" if w3 else "REFUTED"))

    A = math.exp(float(np.polyfit(np.log(np.array(NS, dtype=float)),
                                  np.log(dev), 1)[1]))
    Ncrit = (A / 0.01) ** (1.0 / b) if b > 0 else float("inf")
    w4 = Ncrit > 1e10
    say("W4  |0.5 - f| < 0.01 first at N = %.3e   (floor 1e10)   %s"
        % (Ncrit, "hold" if w4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). W2's power law is not the only fit, and")
    say("  the distinction decides the extrapolation. The lean is the")
    say("  finite-N error of Theorem {#thm:C}, whose error term is a LOG")
    say("  power, so |0.5 - f| ~ (log N)^{-c} is at least as natural:")
    L = np.log(np.array(NS, dtype=float))
    c = -float(np.polyfit(np.log(L), np.log(dev), 1)[0])
    Ac = math.exp(float(np.polyfit(np.log(L), np.log(dev), 1)[1]))
    r_pow = float(np.corrcoef(np.log(np.array(NS, dtype=float)),
                              np.log(dev))[0, 1])
    r_log = float(np.corrcoef(np.log(L), np.log(dev))[0, 1])
    say("    power law  |0.5-f| ~ N^{-%.4f}        r = %.5f" % (b, r_pow))
    say("    log law    |0.5-f| ~ (log N)^{-%.4f}   r = %.5f" % (c, r_log))
    Lc = (Ac / 0.01) ** (1.0 / c)
    say("    reaching 0.01: N = %.3e under the power law, N = e^%.1f = "
        "%.1e under the log law" % (Ncrit, Lc, math.exp(min(Lc, 700))))
    say("  Over a factor 32 in N the two are not separable, so the RATE")
    say("  is not determined here -- only that the lean decays.")

    say()
    say("  what this means for Remark {#rem:nocrossk}: G = 1/|2f-1|, so")
    say("  the cross-k cancellation the accessible range does not show")
    say("  arrives as f -> 1/2. Projected G at a few N:")
    for e in (7, 9, 12, 20):
        Nv = 10.0 ** e
        dv = A * Nv ** (-b)
        say("    N = 10^%-3d  |0.5 - f| = %.4f   G = %.1f"
            % (e, dv, 1.0 / (2 * dv)))

    say()
    say("=" * 70)
    ok = w1 and w2 and w3 and w4
    say("the lean is a finite-N effect and decays, slowly" if ok
        else "REFUTED")

    head = [
        "STATISTIC: the mass-weighted fraction f of k with H(N;k) > 0,",
        "           where the mass is (log k)|H(N;k)|; its distance from",
        "           one half; the implied cross-k gain G = 1/|2f-1|; the",
        "           exponent b in |0.5 - f| ~ N^{-b}; and the N at which",
        "           the fit reaches |0.5 - f| = 0.01.",
        "NULL: the coin arm is carried alongside as a reference level --",
        "      a field rebuilt from coin signs on the support of mu^2 has",
        "      no reason to lean, so its f is the scale against which mu's",
        "      departure and its decay are read. Two draws per N.",
        "FIELD: N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6, 6.4e6 with theta' = 0.56,",
        "       so k runs over the squarefree k < N^0.56 coprime to N; m",
        "       over 1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to 6.4e6; numpy default_rng seed 20260808.",
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
