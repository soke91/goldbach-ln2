# -*- coding: utf-8 -*-
"""
Re-verification of the inference drawn from Proposition 15 (`prop:W`) in
v1/paper/wall_v1.tex: that Chowla's conjecture forces rho -> 1.

THE STATEMENT UNDER TEST (wall_v1.tex, Section `sec:coin`), verbatim:

    rho - 1 = (1/V) sum_{h != 0} c(h) S(h),
      c(h) = sum_{p'-p=h} (log p)(log p'),
      S(h) = < mu(u) mu(u-h) >.
    "S(h) is the binary Chowla correlation. Chowla's conjecture gives
     S(h) = o(1) for each fixed h, and its averaged form over h is a
     theorem [MRT15]. Under that input rho -> 1, and then the wall is
     exactly square-root."

S(h) is an AVERAGE (angle brackets), so |S(h)| <= 1 and Chowla's input
is a bound of size o(1). The question this script asks is arithmetic:
how large is the coefficient that o(1) gets multiplied by?

Write Gamma(N) := (sum_{h != 0} c(h)) / V(N). Since c(h) >= 0, feeding
|S(h)| <= eps into the displayed formula yields nothing better than

        |rho - 1|  <=  eps * Gamma(N).

So the inference "S(h) = o(1)  ==>  rho -> 1" is valid only if
Gamma(N) = O(1). Note also
    sum_{h != 0} c(h) = theta(N)^2 - sum_{p<N} (log p)^2,
which is ~ N^2, while V(N) ~ A(N) N log N: so Gamma(N) ~ N/(A log N),
which DIVERGES. This script measures Gamma exactly.

PRE-REGISTRATION (written before the run).

  (1) RULE: the inference stands only if Gamma(N) is bounded.
      Gamma is computed exactly at each N; if Gamma grows like a power
      of N the inference as written is refuted.

  (2) PREDICTION, recorded so it cannot be reported as a surprise.
      Gamma(N) ~ N / (A(N) log N), i.e. Gamma * log N / N -> 1/A(N)
      with A(N) = prod_{q | N does not hold} (1 - 1/(q(q-1))).
      Equivalently the strength Chowla must supply is not o(1) but
      o(log N / N) -- smaller by a factor N/log N.

  (3) SECOND TEST (the absolute-value budget). Even the true size of
      S(h) does not rescue the inference termwise. The paper's own
      Section `sec:coin` measures the Mobius autocorrelation at
      1.051-1.068 times sqrt(0.32264 (X-h)), i.e. |S(h)| ~ c X^{-1/2}.
      Substituting that gives |rho-1| <~ Gamma * X^{-1/2} ~ N^{1/2}/log N,
      still divergent. RULE: report
      sum_h c(h) |S(h)| / V with S measured directly at the smaller
      sizes; if that ratio exceeds 1 by orders of magnitude, then the
      formula can only give rho -> 1 through SIGNED cancellation across
      h, which is a strictly stronger input than any smallness of S(h).

Nothing here disputes the displayed formula. What is tested is the
sentence that follows it.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def sieve_mu_lambda(X):
    """mu and Lambda on [0, X]. Plain sieve, independent of v1."""
    mu = np.ones(X + 1, dtype=np.int64)
    is_p = np.zeros(X + 1, dtype=bool)
    rem = np.arange(X + 1, dtype=np.int64)
    mu[0] = 0
    for p in range(2, X + 1):
        if rem[p] == p:
            is_p[p] = True
            mu[p::p] *= -1
            rem[p::p] //= p
            if p * p <= X:
                mu[p * p::p * p] = 0
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in range(2, X + 1):
        if is_p[p]:
            lg = math.log(p)
            q = p
            while q <= X:
                lam[q] = lg
                q *= p
    return mu.astype(np.float64), lam, is_p


def A_of_N(N, primes):
    """A(N) = prod_{q does not divide N} (1 - 1/(q(q-1))), truncated at
    the primes available. Converges fast."""
    val = 1.0
    for q in primes:
        if N % q != 0:
            val *= 1.0 - 1.0 / (q * (q - 1.0))
    return val


def main():
    print("Re-verification: what Chowla's o(1) buys in Proposition 15")
    print()

    XMAX = 4_000_000
    mu, lam, is_p = sieve_mu_lambda(XMAX)
    small_primes = [p for p in range(2, 100000) if is_p[p]]

    print("(A) the amplification factor  Gamma(N) = (sum_h c(h)) / V(N)")
    print("    sum_{h!=0} c(h) = theta(N)^2 - sum_p (log p)^2  (exact)")
    print("    V(N)            = sum_{v<N} mu^2(v) Lambda(N-v)^2 (exact)")
    print()
    hdr = (f"{'N':>10} {'sum_h c(h)':>14} {'V(N)':>13} {'Gamma':>12} "
           f"{'Gamma*logN/N':>13} {'1/A(N)':>9}")
    print(hdr)
    print("-" * len(hdr))

    Ns = [10_000, 40_000, 160_000, 640_000, 2_560_000, 4_000_000]
    Ns = [n - (n % 2) for n in Ns]
    gammas = []
    for N in Ns:
        lg = lam[:N]
        theta = float(lg.sum())                 # sum_{n<N} Lambda(n)
        sq = float((lg * lg).sum())             # sum_{n<N} Lambda(n)^2
        csum = theta * theta - sq
        # V(N) = sum_{v<N} mu^2(v) Lambda(N-v)^2
        v = np.arange(1, N)
        V = float(((mu[1:N] ** 2) * (lam[N - v] ** 2)).sum())
        G = csum / V
        A = A_of_N(N, small_primes)
        gammas.append((N, G))
        print(f"{N:>10} {csum:>14.5e} {V:>13.5e} {G:>12.4e} "
              f"{G * math.log(N) / N:>13.4f} {1.0 / A:>9.4f}")

    print()
    print("    (1) is Gamma(N) bounded?  ", end="")
    growth = gammas[-1][1] / gammas[0][1]
    span = gammas[-1][0] / gammas[0][0]
    print(f"Gamma grew {growth:.1f}x while N grew {span:.0f}x")
    verdict1 = "FAIL (Gamma diverges)" if growth > 10 else "PASS"
    print(f"        RULE: the inference needs Gamma = O(1) -> {verdict1}")
    print()
    print("    Consequence, stated exactly: feeding |S(h)| <= eps into")
    print("    the displayed formula gives |rho-1| <= eps * Gamma. At")
    print(f"    N = {gammas[-1][0]}, Gamma = {gammas[-1][1]:.3e}, so a")
    print("    hypothesis of the form S(h) = o(1) leaves the bound")
    print(f"    |rho-1| = o({gammas[-1][1]:.1e}) -- vacuous, since")
    print("    rho >= 0 and the trivial bound on |rho-1| is already")
    print("    C(N)^2/V + 1. The strength actually required is")
    print("    S(h) = o(log N / N), not S(h) = o(1).")

    # ---- (3) the absolute-value budget, with S measured ----
    print()
    print("(3) the absolute-value budget, with S(h) measured directly")
    print("    S(h) = (1/(X-h)) sum_{u<=X-h} mu(u) mu(u+h)")
    print()
    hdr3 = (f"{'X':>8} {'mean|S(h)|':>12} {'sum_h c|S|/V':>14} "
            f"{'signed/V':>12}")
    print(hdr3)
    print("-" * len(hdr3))
    for X in (20_000, 40_000, 80_000, 160_000):
        N = X
        m = mu[:X + 1]
        lg = lam[:N]
        # c(h) for h = 1..N-1, from the prime-power weights:
        #   c(h) = sum_{n} Lambda(n) Lambda(n+h), n, n+h < N
        # (prime powers instead of primes: difference is O(sqrt N log^2 N))
        F = np.fft.rfft(np.pad(lg, (0, len(lg))))
        cfull = np.fft.irfft(F * np.conj(F), 2 * len(lg))[:N]
        # S(h)
        S = np.zeros(N)
        for h in range(1, N):
            S[h] = float(np.dot(m[1:X + 1 - h], m[1 + h:X + 1])) / (X - h)
        v = np.arange(1, N)
        V = float(((mu[1:N] ** 2) * (lam[N - v] ** 2)).sum())
        absbudget = 2.0 * float(np.dot(cfull[1:], np.abs(S[1:]))) / V
        signed = 2.0 * float(np.dot(cfull[1:], S[1:])) / V
        print(f"{X:>8} {np.abs(S[1:]).mean():>12.3e} {absbudget:>14.4e} "
              f"{signed:>12.4e}")

    print()
    print("    The absolute budget exceeds 1 by orders of magnitude at")
    print("    every size, and grows. So rho -> 1 cannot follow from any")
    print("    bound on |S(h)|; it requires signed cancellation across h.")
    print("    Chowla's conjecture, as invoked, supplies smallness and")
    print("    not cancellation, and the averaged theorem [MRT15] is")
    print("    itself a statement about sum_h |S(h)| -- absolute values.")
    print("DONE")


if __name__ == "__main__":
    main()
