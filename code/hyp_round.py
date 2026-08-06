# -*- coding: utf-8 -*-
"""
Hypothesis round (increment 210): three new hypotheses about the wall,
each with a pre-registered kill-test.

H1 -- WEIGHTED SIGNED SUMS. The chain consumes the SIGNED sum
Sum_{k~K} b_k D(k), not the L^2 norm; Cauchy-Schwarz costs sqrt(K) and
throws away b's structure. R4 tested only b = 1. The Vaughan weights
that actually occur are log k, mu(k), Lambda(k). One of them is
special: over the FULL ranges, Sum_k mu(k) D(k) = Sum_u mu(N-u)
(mu*mu)(u), a linear mu-sum against a bounded multiplicative function
-- tractable in principle. Does any weight show extra cancellation in
a restricted band?
  ALIVE iff some b gives |S_b| / (||b||_2 * rms(D)) <= 0.3 at two
        bands and two N; DEAD if every weight sits at the random-sign
        value ~1.

H2 -- A MASK ON C(N). Conjecture L says every mu-field is
(computable mask) x (Gaussian). C(N) ~ N^{0.503} was measured pooled.
If C(N)/sqrt(N) has a mask, its spread should depend on N's small
prime structure.
  ALIVE iff sd(C/sqrt N) differs by a factor >= 2 between bins of N
        (by v_2(N), by 3|N); DEAD if all bins agree to within 20%.

H3 -- CONCENTRATION. If Sum_k |D(k)|^2 were carried by few k, those k
could be isolated and handled separately.
  Participation ratio P = (Sum|D|^2)^2 / (K * Sum|D|^4); for Gaussian
  P = 1/3.
  ALIVE iff P <= 0.15 (far more concentrated than Gaussian); DEAD if
  P is near 1/3.
"""
import numpy as np
import math

from e1_forge_r4 import mobius_upto


def field(mu, N, ks):
    SQ = int(N ** 0.5)
    D = np.zeros(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        hi = N // k
        if hi <= SQ:
            continue
        ms = np.arange(SQ + 1, hi + 1, dtype=np.int64)
        D[i] = float((mu[ms].astype(np.int64)
                      * mu[N - k * ms].astype(np.int64)).sum())
    return D


def h1_h3(mu, N, bands):
    print(f"\n--- N = {N} ---")
    print(f"{'K':>7} {'#k':>5} {'rms D':>9} "
          f"{'b=1':>8} {'b=log':>8} {'b=mu':>8} {'b=mulog':>8} "
          f"{'b=Lam':>8} {'P':>7}")
    for K, nk in bands:
        ks = np.arange(K, K + nk, dtype=np.int64)
        D = field(mu, N, ks)
        rms = math.sqrt(float(np.dot(D, D)) / len(D))
        if rms == 0:
            continue
        lam = np.zeros(len(ks))
        for i, k in enumerate(ks):
            k = int(k)
            p = 2
            q = k
            fac = []
            while p * p <= q:
                if q % p == 0:
                    fac.append(p)
                    while q % p == 0:
                        q //= p
                p += 1
            if q > 1:
                fac.append(q)
            lam[i] = math.log(fac[0]) if len(fac) == 1 else 0.0
        weights = {
            "b=1": np.ones(len(ks)),
            "b=log": np.log(ks.astype(np.float64)),
            "b=mu": mu[ks].astype(np.float64),
            "b=mulog": mu[ks].astype(np.float64)
                       * np.log(ks.astype(np.float64)),
            "b=Lam": lam,
        }
        out = []
        for name, b in weights.items():
            nb = math.sqrt(float(np.dot(b, b)))
            if nb == 0:
                out.append(float('nan')); continue
            out.append(abs(float(np.dot(b, D))) / (nb * rms))
        m2 = float(np.dot(D, D)) / len(D)
        m4 = float(np.dot(D ** 2, D ** 2)) / len(D)
        P = m2 * m2 / m4
        print(f"{K:>7} {len(ks):>5} {rms:>9.1f} "
              + " ".join(f"{v:>8.3f}" for v in out) + f" {P:>7.3f}")


def h2(mu, lam, spf, Xmax):
    print("\n--- H2: is there a mask on C(N)? ---")
    bins = {"v2=1": [], "v2=2": [], "v2>=3": [], "3|N": [], "3!|N": []}
    for N in range(400_000, 400_000 + 2 * 700, 2):
        idx = np.arange(1, N)
        C = float(np.dot(lam[1:N], mu[N - idx].astype(np.float64)))
        c = C / math.sqrt(N)
        v2 = 0
        t = N
        while t % 2 == 0:
            v2 += 1; t //= 2
        bins["v2=1" if v2 == 1 else ("v2=2" if v2 == 2 else "v2>=3")
             ].append(c)
        bins["3|N" if N % 3 == 0 else "3!|N"].append(c)
    print(f"{'bin':>8} {'#N':>5} {'mean':>9} {'sd':>9}")
    for k, v in bins.items():
        a = np.array(v)
        print(f"{k:>8} {a.size:>5} {a.mean():>9.3f} {a.std():>9.3f}")
    r = (max(np.std(bins['v2=1']), np.std(bins['v2=2']))
         / max(min(np.std(bins['v2=1']), np.std(bins['v2=2'])), 1e-9))
    print(f"  v2 bin sd ratio = {r:.3f}   "
          f"(ALIVE if >= 2, DEAD if <= 1.2)")


def main():
    N1, N2 = 19_999_998, 9_999_998
    mu = mobius_upto(N1)
    h1_h3(mu, N1, [(300, 400), (1000, 400), (3000, 400)])
    h1_h3(mu, N2, [(300, 400), (1000, 400)])
    del mu

    X = 900_000
    from thmC_alpha_scan import sieve
    mu2, lam2, phi2, primes2, spf2 = sieve(X)
    h2(mu2, lam2, spf2, X)

    # verdict-ok: criterion: states the rule in advance, must be fixed text
    print("\nPre-registered: H1 ALIVE if some weight <= 0.3; H2 ALIVE "
          "if bin sd ratio >= 2; H3 ALIVE if P <= 0.15.")
    print("DONE")


if __name__ == "__main__":
    main()
