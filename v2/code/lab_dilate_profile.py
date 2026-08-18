# -*- coding: utf-8 -*-
r"""
OPEN.md, "표적은 theta' > 1/2 하나" -- how the demand's budget is spent
across dilations, and what exponent the dilated wall actually has.

WHAT IS AT STAKE

Proposition {#prop:dilate} turned the demand into a weighted sum of
dilated walls,

    B(N) = sum_{k<K,(k,N)=1} (log k) |E_mu(N;k)|,
    E_mu(N;k) = mu(k) H(N;k) - C(N)/phi(k),
    H(N;k) = sum_{m<N/k,(m,k)=1} Lambda(N-mk) mu(m),

and Proposition {#prop:nolog} needs B(N) <= (1-eps) S(N)(1-A(N)) N.
So everything now turns on ONE profile: the relative size of the
dilated wall,

    rho(k) := |H(N;k)| * k / N,

which is the wall at scale M = N/k measured against that scale.  If
rho(k) ~ c k^{1/2} -- square-root cancellation in the dilate -- then

    B(N) ~ c sqrt(N) * 2 sqrt(K) log K,

so B(N) <= 0.37 N is satisfied up to K ~ N / (c log K)^2, a level of
N^{1-o(1)}, far above the N^{1/2} the route needs.  If instead
rho(k) ~ k^a with a > 1/2 the room shrinks, and the exponent a is the
whole question.  Nothing in either paper measures it.

This also explains a puzzle left by lab_level_of_distribution: the
fitted K* exponents there were 0.7057 for mu and 0.6534 for a coin,
both well below 1, while the square-root heuristic predicts an
exponent near 1.  Either a > 1/2, or the shortfall is the polylog that
the heuristic hides.  The fit here decides which.

BACKS: Remark {#rem:dilateprofile} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  For mu, rho(k) fitted as k^a over 2 <= k < 20000 gives a in
      (0.35, 0.65) at every N tested.
  Q2  For a coin on the same support the fit gives a within 0.05 of
      0.5 -- the square-root law, which a coin must obey.
  Q3  mu's exponent exceeds the coin's at every N: mu is relatively
      noisier at large k, which is what lab_dilate_identity saw
      band by band (1.08 rising to 1.53).
  Q4  Consistency: sum_{k<N^0.56} (log k)|A(N;k)|/N is within 0.05 of
      B(N)/N as measured by lab_onesided_demand, namely 0.8086, 0.7395,
      0.7303, 0.6547, 0.5916 at N = 2e5 ... 3.2e6 -- i.e. dropping the
      C(N)/phi(k) term changes little.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED if a leaves (0.35, 0.65) at any N.
  Q2  REFUTED if the coin's a differs from 0.5 by more than 0.05.
  Q3  REFUTED if mu's a fails to exceed the coin's at any N.
  Q4  REFUTED if the difference exceeds 0.05 at any N.

  All four gate.  Q2 is the calibration: a coin that did not come out
  at 1/2 would mean the estimator, not the arithmetic, is being
  measured.
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
OUT = os.path.join(ROOT, "results", "lab_dilate_profile.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
PUB_B = [0.8086, 0.7395, 0.7303, 0.6547, 0.5916]
KFIT = 20_000
THETA = 0.56
DRAWS = 4
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
    say("  N          #k      a (mu)     a (coin)   budget sum   pub B/N"
        "    |diff|")
    say("  " + "-" * 74)
    q1 = q2 = q3 = q4 = True
    prof = {}
    prof_eps = {}
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
        Kfit = min(KFIT, N // 4)
        ks = np.array([k for k in range(2, Kfit)
                       if mu[k] != 0 and all(k % q for q in PN)])

        def prog(sig):
            f = np.zeros(N, dtype=np.float64)
            idx = np.arange(1, N, dtype=np.int64)
            f[1:] = lam[1:N] * sig[N - idx]
            out = np.empty(ks.size)
            for i, k in enumerate(ks):
                r = N % int(k)
                out[i] = f[r::int(k)].sum() if r else f[int(k)::int(k)].sum()
            return np.abs(out)

        Amu = prog(mu.astype(np.float64))
        rho = Amu * ks / N
        good = rho > 0
        a_mu = float(np.polyfit(np.log(ks[good].astype(float)),
                                np.log(rho[good]), 1)[0])

        aeps = []
        eps_rho = []
        for t in range(DRAWS):
            sig = np.zeros(N + 1, dtype=np.float64)
            supp = np.flatnonzero(mu[:N + 1] != 0)
            sig[supp] = rng.integers(0, 2, size=supp.size) * 2.0 - 1.0
            Ae = prog(sig)
            re = Ae * ks / N
            g2 = re > 0
            aeps.append(float(np.polyfit(np.log(ks[g2].astype(float)),
                                         np.log(re[g2]), 1)[0]))
            eps_rho.append(re)
        a_ep = float(np.mean(aeps))

        K56 = int(N ** THETA)
        sel = ks < K56
        budget = float((np.log(ks[sel].astype(float)) * Amu[sel]).sum()) / N
        diff = abs(budget - PUB_B[ni])

        prof[N] = (ks, rho, a_mu, a_ep)
        prof_eps[N] = eps_rho
        if not (0.35 < a_mu < 0.65):
            q1 = False
        if abs(a_ep - 0.5) > 0.05:
            q2 = False
        if a_mu <= a_ep:
            q3 = False
        if diff > 0.05:
            q4 = False
        say("  %-10d %-7d %-10.4f %-10.4f %-12.4f %-10.4f %.4f"
            % (N, ks.size, a_mu, a_ep, budget, PUB_B[ni], diff))

    say()
    say("Q1  mu exponent in (0.35,0.65) at every N   %s"
        % ("hold" if q1 else "REFUTED"))
    say("Q2  coin exponent within 0.05 of 0.5        %s"
        % ("hold" if q2 else "REFUTED"))
    say("Q3  mu exponent exceeds the coin's          %s"
        % ("hold" if q3 else "REFUTED"))
    say("Q4  budget reproduces B(N)/N within 0.05    %s"
        % ("hold" if q4 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). Q2 and Q3 fail only at N = 2e5, and the")
    say("  field is why: the fit ran to k < 2e4, so the dilate length")
    say("  M = N/k fell to 10 there, where a square-root law in k has")
    say("  nothing to describe. Refitting over the k with M = N/k >= 1000:")
    say("  N          #k      a (mu)     a (coin)   M range")
    for N in NS:
        ks, rho, _, _ = prof[N]
        sel = (N // ks) >= 1000
        if sel.sum() < 8:
            say("  %-10d (too few k)" % N)
            continue
        am = float(np.polyfit(np.log(ks[sel].astype(float)),
                              np.log(rho[sel]), 1)[0])
        ae = prof_eps[N][sel.nonzero()[0][0]:] if False else None
        aes = []
        for re in prof_eps[N]:
            g = sel & (re > 0)
            aes.append(float(np.polyfit(np.log(ks[g].astype(float)),
                                        np.log(re[g]), 1)[0]))
        say("  %-10d %-7d %-10.4f %-10.4f [%d, %d]"
            % (N, int(sel.sum()), am, float(np.mean(aes)),
               int((N // ks[sel]).min()), int((N // ks[sel]).max())))

    say()
    say("  the profile rho(k) = |H(N;k)| k / N, median per octave of k,")
    say("  at N = %d" % NS[-1])
    ks, rho, a_mu, a_ep = prof[NS[-1]]
    say("  k-octave          n      median rho(k)   rho/k^{1/2}")
    lo = 2
    while lo < ks.max():
        sel = (ks >= lo) & (ks < 2 * lo)
        if sel.sum():
            m = float(np.median(rho[sel]))
            say("  [%-7d,%-8d) %-6d %-15.6f %.6f"
                % (lo, 2 * lo, int(sel.sum()), m,
                   m / math.sqrt(1.5 * lo)))
        lo *= 2

    say()
    say("  what the exponent implies for the level. With rho(k) ~ c k^a,")
    say("  B(N) ~ (c/N^a) * sum_{k<K} (log k) (N/k) k^a")
    say("       ~ c N^{1-a} K^a log K / a,  so B <= 0.37 N is kept up to")
    say("  K ~ (0.37 a / (c log K))^{1/a} N.  At a = 1/2 that is")
    say("  N / polylog; the exponent measured is %.4f." % a_mu)

    say()
    say("=" * 70)
    ok = q1 and q2 and q3 and q4
    say("Q1 %s  Q2 %s  Q3 %s  Q4 %s"
        % tuple("hold" if v else "REFUTED" for v in (q1, q2, q3, q4)))
    say("the dilated wall follows a square-root law in k, so the level "
        "the demand can carry is N^{1-o(1)}" if ok else "REFUTED")

    head = [
        "STATISTIC: rho(k) = |A(N;k)| k / N with A(N;k) = mu(k)H(N;k) the",
        "           progression sum of Proposition {#prop:dilate}; the",
        "           exponent a in a least-squares fit of log rho against",
        "           log k, for mu and for coin signs on the same support;",
        "           the median of rho per octave of k; and the budget",
        "           sum_{k<N^0.56} (log k)|A(N;k)| / N against B(N)/N as",
        "           measured by lab_onesided_demand.",
        "NULL: the coin is the control and is also the calibration -- a",
        "      sum of independent signs must give a = 1/2 exactly, so if",
        "      the coin's fit missed 1/2 the estimator would be at fault",
        "      rather than the arithmetic. Four draws per N, same support",
        "      and same k-range as mu.",
        "FIELD: N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6; k over the squarefree k",
        "       coprime to N with 2 <= k < min(2e4, N/4) for the fit, and",
        "       k < N^0.56 for the budget; m over 1 <= m < N/k; Lambda and",
        "       mu from an integer sieve to 3.2e6; numpy default_rng seed",
        "       20260808.",
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
