# -*- coding: utf-8 -*-
"""
Direct test of Theorem C, via its alpha-independence prediction
(increment 209).

Theorem C states, unconditionally,

    E_3(alpha) = r(N) - S(N)(N - C(N)) + O_A(N (log N)^{-A}),

where E_3(alpha) = Sum_{k<K,(k,N)=1} mu(k) log k * E_mu(N;k),
K = (N-1)/alpha, and

    E_mu(N;k) = Sum_{n<N, n = N mod k} Lambda(n) mu(N-n)
                - (1/phi(k)) Sum_{n<N} Lambda(n) mu(N-n).

The left side depends on the truncation K; the right side

    R(N) := r(N) - S(N)(N - C(N))

does NOT -- every term in it is alpha-free. So the theorem predicts
that **E_3(K) is independent of K up to the error term**, with the
spread across K bounding that error from below. That is a sharp,
falsifiable consequence of a result this program claims, and it has
never been tested directly: increment 191's check compared the
constants B(K) -> -S(N), not the identity itself.

PRE-REGISTERED (fixed before the run):
  CONFIRMED  iff for each N the values E_3(K)/N over the swept K agree
             with R(N)/N to within 0.05, and their spread does not
             grow with K.
  PROBLEM    iff E_3(K)/N drifts systematically with K by more than
             0.05, or sits far from R(N)/N: then either the error term
             is large at accessible N, or the identity needs review.

Scope: Theorem C is derived in the Corollary-1 regime theta' > 1/2.
Values theta < 0.5 are swept anyway, and reported separately, since
the theorem does not cover them.
"""
import numpy as np
import math


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]
            sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8)
    mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1)
    for p in primes:
        q = int(p); lp = math.log(int(p))
        while q <= X:
            lam[q] = lp
            q *= int(p)
    phi = np.arange(X + 1, dtype=np.int64)
    for p in primes:
        phi[p::p] -= phi[p::p] // p
    return mu, lam, phi, primes, spf


def singular(N, spf):
    C2 = 0.6601618158468696
    S = 2 * C2
    n = N
    while n > 1:
        p = int(spf[n])
        if p > 2:
            S *= (p - 1) / (p - 2)
        while n % p == 0:
            n //= p
    return S


def main():
    X = 1_200_000
    print(f"sieving to {X} ...", flush=True)
    mu, lam, phi, primes, spf = sieve(X)

    thetas = [0.30, 0.40, 0.50, 0.56, 0.62, 0.70]
    for N in (399_998, 799_998, 1_199_998):
        idx = np.arange(1, N)
        f = lam[1:N] * mu[N - idx].astype(np.float64)   # f[n-1]
        C = float(f.sum())
        r = float(np.dot(lam[1:N], lam[N - idx]))
        S = singular(N, spf)
        R = r - S * (N - C)
        print(f"\nN = {N}   C/N = {C/N:+.5f}   r/(S N) = {r/(S*N):.5f}"
              f"   R/N = {R/N:+.5f}")
        print(f"  {'theta':>6} {'K':>8} {'E3/N':>10} {'(E3-R)/N':>10}")
        vals = []
        for th in thetas:
            K = int(N ** th)
            tot = 0.0
            Bsum = 0.0
            for k in range(2, K):
                if mu[k] == 0 or math.gcd(k, N) != 1:
                    continue
                n0 = N % k
                if n0 == 0:
                    n0 = k
                s = float(f[n0 - 1::k].sum())
                w = mu[k] * math.log(k)
                tot += w * s
                Bsum += w / phi[k]
            E3 = tot - C * Bsum
            vals.append(E3 / N)
            print(f"  {th:>6.2f} {K:>8} {E3/N:>10.5f} "
                  f"{(E3-R)/N:>10.5f}")
        v = np.array(vals)
        cov = v[2:]              # theta >= 0.50, the theorem's regime
        print(f"  spread over theta>=0.5: {cov.max()-cov.min():.5f}"
              f"   max |E3-R|/N there: "
              f"{max(abs(x - R/N) for x in cov):.5f}")

        # WHY it drifts: for squarefree k, mu(k)*mu(mk) = mu(m), so the
        # k-weight in E_3's first part is log k with NO sign -- the same
        # mu(k)^2 = 1 mechanism as Theorem A. Check that identity, then
        # compare the drift against the size it predicts.
        K = int(N ** 0.62)
        lhs = rhs = 0.0
        for k in range(2, K):
            if mu[k] == 0 or math.gcd(k, N) != 1:
                continue
            n0 = N % k
            if n0 == 0:
                n0 = k
            lhs += mu[k] * float(f[n0 - 1::k].sum())
            ms = np.arange(1, (N - 1) // k + 1)
            ms = ms[np.gcd(ms, k) == 1]
            rhs += float(np.dot(mu[ms].astype(np.float64),
                                lam[N - k * ms]))
        print(f"  identity check at K={K}: "
              f"sum_k mu(k)*AP = {lhs:.4f}, "
              f"sum_k sum_m mu(m)Lam(N-mk) = {rhs:.4f}, "
              f"diff = {abs(lhs-rhs):.3e}")
        pred = 2 * math.sqrt(N * K) * math.log(K)
        print(f"  predicted |E3| ~ 2 sqrt(NK) log K = {pred:.3e}"
              f"   measured |E3| = {abs(vals[4])*N:.3e}"
              f"   ratio = {abs(vals[4])*N/pred:.3f}")

    # verdict-ok: criterion: states the rule in advance
    print("\nPre-registered: CONFIRMED if |E3-R|/N <= 0.05 with no")
    print("systematic drift in K; PROBLEM otherwise.")
    print("DONE")


if __name__ == "__main__":
    main()
