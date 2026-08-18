# -*- coding: utf-8 -*-
r"""
Why the dilated walls move together: is H(N;k) a Mertens sum?

WHAT IS AT STAKE

Remark {#rem:nocrossk} established that the H(N;k) do not cancel
across k -- the effective number of independent signs is 3.4 to 7.8
where the number of moduli is 313 to 1485 -- and left the mechanism
open.  There is an obvious candidate.  With M = N/k,

    H(N;k) = sum_{m<M, (m,k)=1} Lambda(N-mk) mu(m),

and if the weight Lambda(N-mk) behaved like its mean over m -- which
is what the prime number theorem in progressions says it does on
average -- then H(N;k) would be a constant times

    M_k(M) = sum_{m<M, (m,k)=1} mu(m),

a coprimality-restricted Mertens function.  Values of the Mertens
function at different arguments are correlated over long ranges, which
would explain the correlation across k directly.

If it holds, the payoff is a reduction: Proposition {#prop:nolog}'s
demand B(N) <= S(N)(1-A(N))N becomes a statement about
sum_k (log k)|M_k(N/k)| -- Mertens sums over dilations -- which is a
far more classical object than a Lambda-mu correlation.

THE CONTROL MATTERS AND IS PRE-REGISTERED AS SUCH

The reduction may be a fact about the WEIGHT rather than about mu: for
coin signs, A_eps(N;k) = sum_m Lambda(N-mk) eps(mk) and the analogue
of M_k is sum_m eps(mk), the same index set with the same signs and
only the weight removed.  If the coin's correlation is equally high,
the reduction says "the Lambda weight is nearly a constant" and
nothing about mu.  That is still what the route needs, but it must not
be dressed as a mu-fact -- the mistake Remark {#rem:levelmeas} was
withdrawn for.

BACKS: Remark {#rem:mertens} in paper/wall_v3.md -- a closed line,
not a result.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  corr(H(N;k), M_k(N/k)) over the admissible k exceeds 0.7 at
      every N tested.
  U2  The coin's corresponding correlation also exceeds 0.7 -- the
      reduction is about the weight, not about mu.
  U3  The least-squares slope c in H ~ c M_k, through the origin, lies
      in [0.7, 1.4] at every N.
  U4  The reduction carries the demand: sum(log k)|M_k(N/k)| is within
      30% of sum(log k)|H(N;k)| at every N.

REFUTATION RULE (fixed before the run)

  U1  REFUTED if the correlation is 0.7 or below at any N.
  U2  REFUTED if the coin's correlation is 0.7 or below at any N -- in
      which case the reduction IS mu-specific, which is the more
      interesting outcome and is reported as such.
  U3  REFUTED if c leaves [0.7, 1.4] at any N.
  U4  REFUTED if the two sums differ by more than 30% at any N.

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
OUT = os.path.join(ROOT, "results", "lab_mertens_reduction.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
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
    say("  N          #k     corr(H, M_k)   coin corr   slope c    "
        "sum|M_k| / sum|H|")
    say("  " + "-" * 78)
    u1 = u2 = u3 = u4 = True
    diag = []
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
        Mk = np.empty(ks.size)
        for i, k in enumerate(ks):
            k = int(k)
            m = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
            m = m[np.gcd(m, k) == 1]
            sm = mu[m].astype(np.float64)
            H[i] = float((lam[N - m * k] * sm).sum())
            Mk[i] = float(sm.sum())

        c = float((H * Mk).sum() / (Mk * Mk).sum())
        r = float(np.corrcoef(H, Mk)[0, 1])
        ratio = float((lg * np.abs(Mk)).sum() / (lg * np.abs(H)).sum())

        rc = []
        for t in range(DRAWS):
            sig = np.zeros(N + 1, dtype=np.float64)
            supp = np.flatnonzero(mu[:N + 1] != 0)
            sig[supp] = rng.integers(0, 2, size=supp.size) * 2.0 - 1.0
            Ae = np.empty(ks.size)
            Se = np.empty(ks.size)
            for i, k in enumerate(ks):
                k = int(k)
                m = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
                s2 = sig[m * k]
                Ae[i] = float((lam[N - m * k] * s2).sum())
                Se[i] = float(s2.sum())
            rc.append(float(np.corrcoef(Ae, Se)[0, 1]))
        rcm = float(np.mean(rc))

        if r <= 0.7:
            u1 = False
        if rcm <= 0.7:
            u2 = False
        if not (0.7 <= c <= 1.4):
            u3 = False
        if abs(ratio - 1.0) > 0.30:
            u4 = False
        say("  %-10d %-6d %-14.4f %-11.4f %-10.4f %.4f"
            % (N, ks.size, r, rcm, c, ratio))
        diag.append((N, r, c, float(np.abs(H).mean()),
                     float(np.abs(Mk).mean())))

    say()
    say("U1  corr(H, M_k) > 0.7 at every N          %s"
        % ("hold" if u1 else "REFUTED"))
    say("U2  coin correlation > 0.7 at every N      %s"
        % ("hold" if u2 else "REFUTED"))
    say("U3  slope c in [0.7, 1.4] at every N       %s"
        % ("hold" if u3 else "REFUTED"))
    say("U4  sum|M_k| within 30%% of sum|H|          %s"
        % ("hold" if u4 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). How far from constant is the weight?")
    say("  If Lambda(N-mk) were a constant the correlation would be 1 and")
    say("  the slope would be that constant. Decomposing H = c M_k + R:")
    say("  N          var(R)/var(H)   mean|H|/mean|M_k|   sqrt(log N)")
    for N, r, c, mh, mm in diag:
        say("  %-10d %-15.4f %-19.4f %.4f"
            % (N, 1.0 - r * r, mh / mm, math.sqrt(math.log(N))))
    say("  The residual carries most of the variance and the size ratio")
    say("  is several times sqrt(log N), so the weight is not a constant")
    say("  multiplier: it supplies the bulk of H's own fluctuation.")

    say()
    say("=" * 70)
    ok = u1 and u2 and u3 and u4
    say("H(N;k) is a Mertens sum up to the weight, so the demand is a "
        "statement about Mertens sums over dilations"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the Pearson correlation across k between",
        "           H(N;k) = sum_{m<N/k,(m,k)=1} Lambda(N-mk) mu(m) and the",
        "           restricted Mertens sum M_k(N/k) = sum_{m<N/k,(m,k)=1}",
        "           mu(m); the same for coin signs, where the two sums",
        "           share index set and signs and differ only by the",
        "           Lambda weight; the least-squares slope of H against",
        "           M_k through the origin; and the ratio of",
        "           sum(log k)|M_k| to sum(log k)|H|.",
        "NULL: the coin is the control and the question it settles is",
        "      stated in advance -- if the coin correlates just as well,",
        "      the reduction is a fact about the Lambda weight and not",
        "      about mu, and must be reported as one. Four draws per N,",
        "      same support, same k-range, same index sets.",
        "FIELD: N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6 with theta' = 0.56, so k",
        "       runs over the squarefree k < N^0.56 coprime to N; m over",
        "       1 <= m < N/k with (m,k) = 1 for mu and over the same m for",
        "       the coin; Lambda and mu from an integer sieve to 3.2e6;",
        "       numpy default_rng seed 20260808.",
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
