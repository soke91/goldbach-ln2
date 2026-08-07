# -*- coding: utf-8 -*-
"""
Direct test of the conclusion of Theorem 1 (`thm:A`) of
v1/paper/wall_v1.tex and of `theorem_A.tex`.

THE STATEMENT UNDER TEST, verbatim:

    Fix theta' in (1/2,1) and put K = N^{theta'}. Then for every A>0,
      sup_{1<=t<N} | sum_{k<K, (k,N)=1} mu(k) E_mu(t;k) |  <<_{A,theta'}
      N (log N)^{-A}.

with
    E_mu(t;k) = sum_{n<=t, n=N mod k} Lambda(n) mu(N-n)
                - (1/phi(k)) sum_{n<=t} Lambda(n) mu(N-n).

METHOD HERE. Written from the statement. The whole object is computed
EXACTLY, for every t at once, with no appeal to the proof's
decomposition. Since n = N mod k is k | N-n, put u = N-n; then

    T_1(t) = sum_{u >= N-t} Lambda(N-u) mu(u) sigma_K(u)
             - C(t) * B(K),
    sigma_K(u) = sum_{k|u, k<K, (k,N)=1} mu(k),
    B(K)       = sum_{k<K, (k,N)=1} mu(k)/phi(k),
    C(t)       = sum_{n<=t} Lambda(n) mu(N-n).

sigma_K is built by sieving: for each admissible k, add mu(k) to every
multiple of k. That costs sum_{k<K} N/k and is exact in integer
arithmetic. Taking a cumulative sum in n then gives T_1(t) for EVERY t
simultaneously, so the supremum in the statement is attained, not
sampled.

v1's own numerics (theorem_A.tex, "Numerical verification") never
compute T_1 itself: they compute the residual R and its predicted main
term MT and report that the two agree to 1-4%. That is a check of
Steps 3-6, not of the theorem's conclusion. This script closes that
gap.

PRE-REGISTRATION (written before the run).

  (1) RULE. sup_t |T_1(t)| / N must fall as N grows. If it is flat or
      rising over the tested decades, the theorem's conclusion is not
      visible at accessible N and the measurement says so; if it falls
      faster than any power of log N (i.e. like a power of N) the
      theorem is true but is not what the data are showing either.
      Both outcomes are reported, neither is a refutation on its own:
      the theorem is asymptotic.

  (2) PREDICTION, recorded so it cannot be reported as a surprise.
      theorem_A.tex measures the RESIDUAL to decay like N^{1-theta'/2},
      i.e. R/N ~ N^{-theta'/2} = N^{-0.28} at theta'=0.56, and reports
      R and MT agreeing to 1-4%. T_1 is what is left after that
      cancellation, so I predict sup_t|T_1|/N is SMALLER than R/N by
      roughly that 1-4%, i.e. of order 1e-3, and decaying at least as
      fast as N^{-0.28}.

  (3) SECOND TEST, the one that can embarrass either side. The
      dominant term of the proof is the main term MT, killed by the
      cancellation in sum_m f(m)/m (Lemma 15 + Proposition 17). The
      script also reports the same quantity with the sign of mu(k)
      DISCARDED,
          Tabs(t) = sum_{k<K,(k,N)=1} |E_mu(t;k)|,
      which is what EH_mu(N^{theta'}) asserts is small. Theorem 1 is
      explicitly NOT an instance of EH_mu, so Tabs/N must stay large
      while T_1/N falls. If both fall together, the theorem is weaker
      than advertised -- it would be measuring nothing that the
      absolute-value form does not already give.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

THETA = 0.56


def sieve_mu_lambda_phi(X):
    mu = np.ones(X + 1, dtype=np.int64)
    is_p = np.zeros(X + 1, dtype=bool)
    rem = np.arange(X + 1, dtype=np.int64)
    phi = np.arange(X + 1, dtype=np.float64)
    mu[0] = 0
    for p in range(2, X + 1):
        if rem[p] == p:
            is_p[p] = True
            mu[p::p] *= -1
            rem[p::p] //= p
            phi[p::p] -= phi[p::p] / p
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
    return mu.astype(np.int64), lam, phi


def run(N, mu, lam, phi):
    K = int(N ** THETA)
    # admissible k: k < K, (k,N) = 1, mu(k) != 0
    ks = np.arange(1, K, dtype=np.int64)
    ks = ks[mu[1:K] != 0]
    ks = ks[np.gcd(ks, N) == 1]

    # sigma_K(u) for u = 1..N-1
    sigma = np.zeros(N, dtype=np.int64)
    for k in ks:
        sigma[k::k] += int(mu[k])

    u = np.arange(1, N)
    # terms indexed by n = N-u, i.e. n = 1..N-1
    term_D = lam[N - u] * mu[u] * sigma[1:]          # index j -> u=j+1
    term_C = lam[N - u] * mu[u]
    # order by n ascending: n = N-u, so u descending
    D_by_n = term_D[::-1]
    C_by_n = term_C[::-1]
    Dcum = np.cumsum(D_by_n)          # D(t) for t = 1..N-1
    Ccum = np.cumsum(C_by_n)          # C(t)

    B = float((mu[ks] / phi[ks]).sum())
    T1 = Dcum - Ccum * B
    supT1 = float(np.abs(T1).max())

    # the absolute-value form, at t = N-1 only (that is what EH_mu asks)
    # E_mu(N;k) = sum_{u=0 mod k} Lambda(N-u) mu(u) - C(N)/phi(k)
    CN = float(Ccum[-1])
    Tabs = 0.0
    for k in ks:
        s = float((lam[N - np.arange(k, N, k)]
                   * mu[np.arange(k, N, k)]).sum())
        Tabs += abs(s - CN / phi[k])
    return K, len(ks), supT1, abs(CN), Tabs, B


def main():
    print("Direct test of Theorem 1: sup_t |T_1(t)| computed exactly")
    print(f"theta' = {THETA}")
    print()

    NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
    XMAX = max(NS)
    mu, lam, phi = sieve_mu_lambda_phi(XMAX)

    hdr = (f"{'N':>9} {'K':>7} {'#k':>7} {'sup|T_1|':>12} "
           f"{'sup|T_1|/N':>12} {'|C(N)|/N':>10} {'B(K)':>10} "
           f"{'Tabs/N':>9}")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for N in NS:
        N -= N % 2
        K, nk, s, cn, ta, B = run(N, mu, lam, phi)
        rows.append((N, s / N, ta / N))
        print(f"{N:>9} {K:>7} {nk:>7} {s:>12.4e} {s/N:>12.4e} "
              f"{cn/N:>10.4e} {B:>10.3e} {ta/N:>9.3f}")

    print()
    print("(1) does sup|T_1|/N fall?  fitted exponent b in N^{-b}:")
    ln = np.log([r[0] for r in rows])
    lt = np.log([r[1] for r in rows])
    b = -np.polyfit(ln, lt, 1)[0]
    print(f"      sup|T_1|/N  ~  N^(-{b:.4f})")
    la = np.log([r[2] for r in rows])
    ba = -np.polyfit(ln, la, 1)[0]
    print(f"      Tabs/N      ~  N^(-{ba:.4f})")
    print()
    print("    for reference, over this range:")
    for A in (1, 2, 3):
        lo = math.log(rows[0][0]) ** -A
        hi = math.log(rows[-1][0]) ** -A
        print(f"      (log N)^-{A}: {lo:.4e} -> {hi:.4e}"
              f"   (a factor {lo/hi:.2f})")
    print(f"      measured    : {rows[0][1]:.4e} -> {rows[-1][1]:.4e}"
          f"   (a factor {rows[0][1]/rows[-1][1]:.2f})")
    print()
    print("(3) the signed sum must fall while the absolute-value form")
    print("    (what EH_mu asserts) does not. Compare the two exponents")
    print("    and the Tabs/N column above.")
    print("DONE")


if __name__ == "__main__":
    main()
