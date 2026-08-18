# -*- coding: utf-8 -*-
r"""
OPEN.md, "표적은 theta' > 1/2 하나" -- the empirical level of
distribution of the Mobius-twisted primes, in the exact norm the
Goldbach route consumes.

WHAT IS AT STAKE

Proposition {#prop:nolog} reduced the demand to a constant-factor bound,

    B(N; K) := sum_{k<K, (k,N)=1} (log k) |E_mu(N;k)|
             <= (1-eps) S(N) (1 - A(N)) N,

and Remark {#rem:relocate} observed that this leaves the whole
difficulty on the LEVEL axis: Huang-Li need K = N^{theta'} with a
single theta' > 1/2, and Bombieri-Vinogradov gives nothing past
N^{1/2}.  What nobody in either paper has measured is where the truth
actually sits on that axis.

So define the empirical level

    K*(N) := max { K : B(N; K) <= S(N)(1 - A(N)) N },

the largest truncation at which the consumed quantity is still small
enough for Goldbach, and ask how it grows.  If K*(N) >> N^{1/2} the
route's demand is true and only a proof is missing; if K*(N) << N^{1/2}
the demand is not even empirically supported at the level it needs,
which would be a far stronger statement than anything in either paper.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  J1  K*(N) exists at every N tested, i.e. B(N; K_max) exceeds the
      threshold before the sweep ends.
  J2  Fitting log K*(N) against log N gives a slope beta in (0.4, 0.7).
  J3  K*(N) > sqrt(N) at every N tested: the truth supports a level
      above the Bombieri-Vinogradov barrier in this norm.
  J4  K*(N)/sqrt(N) is increasing in N.
  J5  At K = floor(N^0.56), B(N)/N reproduces lab_onesided_demand:
      0.8086, 0.7395, 0.7303, 0.6547, 0.5916 at
      N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6.

REFUTATION RULE (fixed before the run)

  J1  REFUTED if the threshold is never crossed at some N.
  J2  REFUTED if beta leaves (0.4, 0.7).
  J3  REFUTED if K*(N) <= sqrt(N) at any N tested.
  J4  REFUTED if K*/sqrt(N) fails to increase between consecutive N.
  J5  REFUTED if any of the five differs by more than 0.001.

  All five gate.  J5 is the tie to the existing measurement.  J3 is the
  question, and it is written the way I expect it to go; a refutation
  there is the more interesting outcome and is reported as such.
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
OUT = os.path.join(ROOT, "results", "lab_level_of_distribution.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
PUB56 = [0.8086, 0.7395, 0.7303, 0.6547, 0.5916]
PLIM = 4_000_000
KCAP = 300_000


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
    say("phi to %d ..." % KCAP)
    ph = phi_sieve(KCAP)

    artin, twin = 1.0, 2.0
    for p in pr:
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    say()
    say("  N          sqrt N    thresh/N   K*        K*/sqrt N   "
        "B at N^0.56 /N   pub")
    say("  " + "-" * 78)
    Ks, got56 = [], []
    j1 = True
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
        A, S = artin, twin
        for q in sorted(PN):
            A /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        thr = S * (1.0 - A) * N

        n = np.arange(N, dtype=np.int64)
        f = np.zeros(N, dtype=np.float64)
        f[1:] = lam[1:N] * mu[(N - n)[1:]]
        C = float(f.sum())

        Kmax = min(KCAP, N // 2)
        K56 = int(N ** 0.56)
        run = 0.0
        Kstar = None
        b56 = None
        for k in range(2, Kmax):
            if k == K56:
                b56 = run
            if mu[k] != 0 and all(k % q for q in PN):
                r = N % k
                inner = float(f[r::k].sum()) if r else float(f[k::k].sum())
                run += math.log(k) * abs(inner - C / int(ph[k]))
            if Kstar is None and run > thr:
                Kstar = k
        if b56 is None:
            b56 = run
        if Kstar is None:
            j1 = False
            Kstar = Kmax
        Ks.append(Kstar)
        got56.append(b56 / N)
        say("  %-10d %-9.1f %-10.4f %-9d %-11.4f %-16.4f %.4f"
            % (N, math.sqrt(N), thr / N, Kstar,
               Kstar / math.sqrt(N), b56 / N, PUB56[NS.index(N)]))

    say()
    xs = np.log(np.array(NS, dtype=float))
    ys = np.log(np.array(Ks, dtype=float))
    beta = float(np.polyfit(xs, ys, 1)[0])
    j2 = 0.4 < beta < 0.7
    say("J1  threshold crossed at every N: %s   %s"
        % (j1, "hold" if j1 else "REFUTED"))
    say("J2  K*(N) ~ N^%.4f   (band (0.4, 0.7))   %s"
        % (beta, "hold" if j2 else "REFUTED"))
    rat = [Ks[i] / math.sqrt(NS[i]) for i in range(len(NS))]
    j3 = all(r > 1.0 for r in rat)
    say("J3  K*/sqrt N = %s   all above 1: %s   %s"
        % (", ".join("%.4f" % r for r in rat), j3,
           "hold" if j3 else "REFUTED"))
    j4 = all(rat[i] < rat[i + 1] for i in range(len(rat) - 1))
    say("J4  K*/sqrt N increasing: %s   %s"
        % (j4, "hold" if j4 else "REFUTED"))
    e5 = max(abs(a - b) for a, b in zip(got56, PUB56))
    j5 = e5 <= 1e-3
    say("J5  B at N^0.56 max deviation %.6f  (tol 0.001)   %s"
        % (e5, "hold" if j5 else "REFUTED"))

    say()
    say("  What the numbers say about the level axis:")
    say("  the route needs K = N^{theta'} with theta' > 1/2; measured,")
    say("  the consumed quantity stays under the Goldbach threshold only")
    say("  up to K* ~ N^%.4f." % beta)
    if beta < 0.5:
        say("  That is BELOW the barrier, so at these N the demand is not")
        say("  merely unproved at the level it needs -- it is false there.")
        say("  Whether that survives N -> infinity is the question, and")
        say("  K*/sqrt N %s over the range measured."
            % ("rises" if j4 else "does not rise"))
    else:
        say("  That is at or above the barrier, so the demand is")
        say("  empirically supported and only a proof is missing.")

    say()
    say("=" * 70)
    ok = j1 and j2 and j3 and j4 and j5
    say("J1 %s  J2 %s  J3 %s  J4 %s  J5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (j1, j2, j3, j4, j5)))
    say("the empirical level clears the barrier" if ok else "REFUTED")

    head = [
        "STATISTIC: K*(N) = max{K : sum_{k<K,(k,N)=1} (log k)|E_mu(N;k)|",
        "           <= S(N)(1-A(N))N}, the largest truncation at which the",
        "           quantity Proposition {#prop:nolog} bounds is still",
        "           under the Goldbach threshold; its exponent in N;",
        "           K*/sqrt(N); and B(N; N^0.56)/N as a tie to",
        "           lab_onesided_demand.",
        "FIELD: N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6; k swept from 2 up to",
        "       min(3e5, N/2); E_mu(N;k) by direct enumeration of the",
        "       progression n = N mod k; Lambda and mu from an integer",
        "       sieve to 4e6, phi from a sieve to 3e5; A(N) and S(N) as",
        "       Euler products over p < 4e6.",
        'NULL: run afterwards, in lab_level_coin_null.py, and it REFUTES the',
        '      reading given here: a coin on the same support reaches a',
        '      HIGHER K* at every N. K* measures the support and square-root',
        '      cancellation, not mu. See Remark {#rem:levelmeas}.',
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
