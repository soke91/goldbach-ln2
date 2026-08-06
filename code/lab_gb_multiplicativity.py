# -*- coding: utf-8 -*-
"""
Transform Lab, session 2 (increment 219): how far is G_b from
multiplicative, and is the truncation what breaks it?

G_b(u) = Sum_{k | u, k ~ K} b_k mu(u/k) is the companion factor in the
switched form Sum_k b_k D(k) = Sum_u mu(N-u) G_b(u). Tao's reduction,
the engine of the Helfgott-Radziwill route, needs multiplicativity of
the function being shifted. This measures whether G_b has any.

The untruncated version is exactly (b * mu)(u), a Dirichlet
convolution, hence multiplicative whenever b is. For b = 1 it is
[u = 1]. So any departure from multiplicativity is caused by the
truncation k ~ K and by nothing else -- which is the point.

Statistic: over coprime pairs (u, v) with uv in range,
    dev = | G(uv) - G(u)G(v) | / (|G(u)G(v)| + 1)
NULLS on the same line as the criteria:
  - a multiplicative function has dev = 0 identically;
  - a function with no multiplicative structure has G(uv) essentially
    independent of G(u)G(v), so corr(G(uv), G(u)G(v)) = 0 with
    SE 1/sqrt(n).
  MULTIPLICATIVE   iff dev = 0 for all pairs.
  SOME STRUCTURE   iff corr >= 4/sqrt(n).
  NONE             iff corr sits at the null.
"""
import numpy as np
import math


def mobius_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool); pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p * p::p] = False
            mu[p::p] *= -1
            mu[p * p::p * p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
    mu[val > 1] *= -1
    mu[0] = 0
    return mu


def divisors(u):
    ds = []
    i = 1
    while i * i <= u:
        if u % i == 0:
            ds.append(i)
            if i != u // i:
                ds.append(u // i)
        i += 1
    return ds


def main():
    X = 400_000
    mu = mobius_upto(X)
    rng = np.random.default_rng(20260909)

    for K, name in ((60, "K=60"), (300, "K=300")):
        lo, hi = K, 2 * K

        def G(u):
            return sum(mu[u // k] for k in divisors(u)
                       if lo <= k < hi)

        A, B = [], []
        devs = []
        tried = 0
        while len(A) < 4000 and tried < 400000:
            tried += 1
            u = int(rng.integers(2, 600))
            v = int(rng.integers(2, 600))
            if math.gcd(u, v) != 1 or u * v > X:
                continue
            guv, gu, gv = G(u * v), G(u), G(v)
            A.append(guv); B.append(gu * gv)
            devs.append(abs(guv - gu * gv) / (abs(gu * gv) + 1))
        A = np.array(A, dtype=float); B = np.array(B, dtype=float)
        n = A.size
        exact = int(np.count_nonzero(np.array(devs) == 0))
        if A.std() > 0 and B.std() > 0:
            r = float(np.corrcoef(A, B)[0, 1])
        else:
            r = float('nan')
        z = r * math.sqrt(n) if r == r else float('nan')
        print(f"{name}: n={n}  exact matches G(uv)=G(u)G(v): "
              f"{exact}/{n} ({exact/n:.3f})")
        print(f"        mean dev = {np.mean(devs):.4f}   "
              f"corr(G(uv), G(u)G(v)) = {r:+.4f}  z = {z:+.2f} "
              f"(null 0, SE {1/math.sqrt(n):.4f})")
        print(f"        verdict: "
              + ("MULTIPLICATIVE" if exact == n else
                 ("SOME STRUCTURE" if abs(z) >= 4 else
                  "NO MULTIPLICATIVE STRUCTURE")))
    print("\nThe untruncated G is (b*mu)(u), multiplicative by")
    print("construction; whatever is measured above is caused by the")
    print("truncation k ~ K alone.")
    print("DONE")


if __name__ == "__main__":
    main()
