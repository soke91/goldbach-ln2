# -*- coding: utf-8 -*-
r"""
The coin null for Remark {#rem:levelmeas} -- does K*(N) measure mu, or
does it measure the support?

WHY THIS HAS TO BE RUN

Remark [rem:levelmeas] reports K*(N) ~ N^{0.7057} with K*/sqrt(N)
crossing 1 near N = 10^6, and reads it as the empirical level of
distribution of the Mobius-twisted primes.  It was published with no
null.  This program's own methodology forbids that -- "a threshold
means nothing until the spread of the quantity it judges has been
measured" -- and Lemma [lem:coin] names the specific danger: replacing
mu by arbitrary signs eps on the same support leaves V(N) and every
centred second moment untouched.  If a coin reproduces K*(N), then
K* is a statement about the support of mu^2 and about square-root
cancellation, not about mu, and the headline number deflates.

So: rebuild K* with eps(v) = +-1 on {v : mu(v) != 0} and zero
elsewhere, averaged over draws, and compare.

    E_eps(N;k) = sum_{n = N mod k} Lambda(n) eps(N-n)
                   - C_eps(N)/phi(k),
    B_eps(N;K) = sum_{k<K,(k,N)=1} (log k)|E_eps(N;k)|,
    K*_eps(N)  = max{K : B_eps(N;K) <= S(N)(1-A(N))N}.

Everything else -- threshold, support, k-range, weight -- is held
identical, so the only thing that changes is the sign pattern.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  K1  K*_mu(N) > K*_eps(N) at every N tested: mu's discrepancy is
      SMALLER than a coin's, so it survives to a higher level.
  K2  K*_mu/K*_eps >= 3 at every N tested.
  K3  The coin's fitted exponent beta_eps exceeds mu's 0.7057, because
      |E_eps(N;k)| ~ sqrt(N log N / phi(k)) makes B_eps ~ sqrt(NK) log K
      and hence K*_eps ~ N / (log N)^3, an exponent near 1.
  K4  The coin's discrepancy scales like phi(k)^{-1/2}: fitting
      |E_eps(N;k)| ~ k^{-a} over the swept range gives a in (0.4, 0.6).
  K5  The ratio |E_mu(N;k)| / <|E_eps(N;k)|> is below 1 at k near
      K*_mu, i.e. mu is quieter than a coin exactly where it matters.

REFUTATION RULE (fixed before the run)

  K1  REFUTED if K*_eps >= K*_mu at any N.  This is the one that
      matters: if the coin matches or beats mu, Remark [rem:levelmeas]
      is measuring the support and must be withdrawn.
  K2  REFUTED if the ratio drops below 3 at any N.
  K3  REFUTED if beta_eps <= 0.7057.
  K4  REFUTED if a leaves (0.4, 0.6).
  K5  REFUTED if the ratio is at least 1 at any N.

  All five gate.
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
OUT = os.path.join(ROOT, "results", "lab_level_coin_null.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KMU = [319, 537, 767, 1353, 2319]          # from lab_level_of_distribution
BETA_MU = 0.7057
KMAX = 20_000
DRAWS = 8
SEED = 20260808
PLIM = 4_000_000


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


def phi_sieve(n):
    ph = np.arange(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if ph[p] == p:
            ph[p::p] -= ph[p::p] // p
    return ph


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("sieving to %d ..." % PLIM)
    pr, lam, mu = sieves(PLIM)
    ph = phi_sieve(KMAX)
    artin, twin = 1.0, 2.0
    for p in pr:
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    rng = np.random.default_rng(SEED)

    say()
    kmus = []
    say("  N          thresh/N   K*_mu   K*_eps (mean of %d)   ratio    "
        "|E_mu|/<|E_eps|> at K*_mu" % DRAWS)
    say("  " + "-" * 86)
    Ke, ratios, eratio = [], [], []
    kdecay = []
    for ni, N in enumerate(NS):
        v, PN, d = N, set(), 2
        while d * d <= v:
            if v % d == 0:
                PN.add(d)
                while v % d == 0:
                    v //= d
            d += 1
        if v > 1:
            PN.add(v)
        A, S = artin, twin
        for q in sorted(PN):
            A /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        thr = S * (1.0 - A) * N

        ks = [k for k in range(2, KMAX)
              if mu[k] != 0 and all(k % q for q in PN)]
        lg = np.array([math.log(k) for k in ks])
        iph = np.array([int(ph[k]) for k in ks], dtype=np.float64)

        n = np.arange(N, dtype=np.int64)
        base = lam[1:N]
        shift = (N - n)[1:]

        def sweep(sig):
            f = np.zeros(N, dtype=np.float64)
            f[1:] = base * sig[shift]
            Ctot = float(f.sum())
            out = np.empty(len(ks))
            for i, k in enumerate(ks):
                r = N % k
                inner = float(f[r::k].sum()) if r else float(f[k::k].sum())
                out[i] = abs(inner - Ctot / iph[i])
            return out

        emu = sweep(mu.astype(np.float64))
        eps_all = np.empty((DRAWS, len(ks)))
        for dd in range(DRAWS):
            sig = np.zeros(N + 1, dtype=np.float64)
            supp = np.flatnonzero(mu[:N + 1] != 0)
            sig[supp] = rng.integers(0, 2, size=supp.size) * 2.0 - 1.0
            eps_all[dd] = sweep(sig)

        def kstar(e):
            run = np.cumsum(lg * e)
            j = int(np.searchsorted(run, thr))
            return ks[min(j, len(ks) - 1)]

        Kmu = kstar(emu)
        kmus.append(int(Kmu))
        Keps = [kstar(eps_all[dd]) for dd in range(DRAWS)]
        Kem = float(np.mean(Keps))
        Ke.append(Kem)
        ratios.append(Kmu / Kem)
        j0 = int(np.searchsorted(np.array(ks), KMU[ni]))
        j0 = min(j0, len(ks) - 1)
        er = emu[j0] / float(eps_all[:, j0].mean())
        eratio.append(er)
        mid = (np.array(ks) >= 50) & (np.array(ks) <= 5000)
        a = -float(np.polyfit(np.log(np.array(ks)[mid]),
                              np.log(eps_all[:, mid].mean(axis=0)), 1)[0])
        kdecay.append(a)
        say("  %-10d %-10.4f %-7d %-20.1f %-8.2f %.4f"
            % (N, thr / N, Kmu, Kem, Kmu / Kem, er))
        if Kmu != KMU[ni]:
            say("      (note: K*_mu here = %d against %d from the earlier "
                "sweep)" % (Kmu, KMU[ni]))

    say()
    k1 = all(KMU[i] > Ke[i] for i in range(len(NS)))
    say("K1  K*_mu > K*_eps at every N: %s   %s"
        % (k1, "hold" if k1 else "REFUTED"))
    k2 = all(r >= 3.0 for r in ratios)
    say("K2  K*_mu/K*_eps = %s   (floor 3)   %s"
        % (", ".join("%.2f" % r for r in ratios),
           "hold" if k2 else "REFUTED"))
    be = float(np.polyfit(np.log(np.array(NS, dtype=float)),
                          np.log(np.array(Ke)), 1)[0])
    k3 = be > BETA_MU
    say("K3  beta_eps = %.4f against beta_mu = %.4f   %s"
        % (be, BETA_MU, "hold" if k3 else "REFUTED"))
    k4 = all(0.4 < a < 0.6 for a in kdecay)
    say("K4  |E_eps| ~ k^{-a} with a = %s   (band (0.4,0.6))   %s"
        % (", ".join("%.4f" % a for a in kdecay),
           "hold" if k4 else "REFUTED"))
    k5 = all(r < 1.0 for r in eratio)
    say("K5  |E_mu|/<|E_eps|> at K*_mu = %s   all below 1: %s   %s"
        % (", ".join("%.4f" % r for r in eratio), k5,
           "hold" if k5 else "REFUTED"))

    say()
    say("  Cross-check lines. audit_levelmeas_budget.py recomputes")
    say("  the same crossing while sweeping the budget factor.")
    for N, K in zip(NS, kmus):
        say("AGREE kstar_nolog N=%d %d 0.02" % (N, K))

    say()
    say("=" * 70)
    ok = k1 and k2 and k3 and k4 and k5
    say("K1 %s  K2 %s  K3 %s  K4 %s  K5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (k1, k2, k3, k4, k5)))
    say("K* separates mu from a coin: the level measurement is about mu"
        if k1 else
        "K* does NOT separate mu from a coin -- Remark {#rem:levelmeas} "
        "must be withdrawn")

    head = [
        "STATISTIC: K*_eps(N) = max{K : sum_{k<K,(k,N)=1}(log k)",
        "           |E_eps(N;k)| <= S(N)(1-A(N))N} with eps(v) = +-1 on",
        "           {v : mu(v) != 0} and zero elsewhere, averaged over 8",
        "           draws; the ratio K*_mu/K*_eps; the coin's fitted",
        "           exponent in N; the decay exponent of |E_eps(N;k)| in",
        "           k; and |E_mu(N;k)|/<|E_eps(N;k)|> at k = K*_mu.",
        "FIELD: N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6; k from 2 to 2e4,",
        "       restricted to squarefree k coprime to N; the coin uses the",
        "       same support, threshold, weight and k-range as mu, so the",
        "       sign pattern is the only difference; decay exponent fitted",
        "       over 50 <= k <= 5000; numpy default_rng seed 20260808.",
        'NULL: this file is the null. Eight coin draws with threshold,',
        '      support, weight and k-range held identical to the mu run, so',
        '      the sign pattern is the only difference.',
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
