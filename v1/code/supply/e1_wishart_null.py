# -*- coding: utf-8 -*-
"""
Increment 154: the CORRECT null -- Gram/Wishart, not iid entries.

Inc. 153 killed the mean hypothesis (E[C] = 0.006 +- 0.005; centering
does not move lambda_max). Diagnosis of inc. 152's design: the real
matrix is a GRAM matrix C = R R^T whose entries share rows; an
iid-entry null underestimates lambda_max even for perfectly random
rows (Wishart vs Wigner). Correct factorization-law null: random
row VECTORS with the real support pattern (same zero positions in
R), iid signs elsewhere, then the same pipeline
(C = RR^T, normalize by sqrt(NZ), zero diagonal). If this null
reproduces 36.76, the z = 9 'gap' was a null-design category error
(correction #26) and the law closes at matrix level after all.
"""
import numpy as np, time

def mobius_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool); pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p*p::p] = False
            mu[p::p] *= -1
            mu[p*p::p*p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
    mu[val > 1] *= -1
    return mu

def primes_in(lo, hi):
    sieve = np.ones(hi - lo + 1, dtype=bool)
    if lo <= 1:
        sieve[:max(0, 2 - lo)] = False
    for p in range(2, int(hi ** 0.5) + 1):
        start = max(p * p, ((lo + p - 1) // p) * p)
        sieve[start - lo::p] = False
    return np.nonzero(sieve)[0] + lo

def lam_of(R):
    Rf = R.astype(np.float32)
    C = Rf @ Rf.T
    NZ = np.abs(Rf) @ np.abs(Rf.T)
    np.fill_diagonal(C, 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        Ch = np.where(NZ > 50, C / np.sqrt(np.maximum(NZ, 1)), 0.0)
    np.fill_diagonal(Ch, 0.0)
    Ch = (Ch + Ch.T) / 2
    return np.abs(np.linalg.eigvalsh(Ch)).max()

def main():
    rng = np.random.default_rng(20260815)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    K0, K1 = 3000, 3400
    ks = np.arange(K0, K1); K = len(ks)
    P0 = N // (2 * K1); P1 = 2 * P0
    ps = primes_in(P0, P1)
    pp = ps[ps <= (N - 2) // (K1 - 1)]

    R = np.empty((K, len(pp)), dtype=np.int8)
    for i, k in enumerate(ks):
        R[i] = mu[N - pp * int(k)]
    lam_real = lam_of(R)
    print(f"real lambda_max = {lam_real:.2f}", flush=True)

    support = (R != 0)
    lams = []
    for d in range(16):
        signs = rng.choice(np.array([-1, 1], dtype=np.int8),
                           size=R.shape)
        Rn = np.where(support, signs, 0).astype(np.int8)
        lams.append(lam_of(Rn))
        if d % 4 == 3:
            print(f"  wishart null {d+1}/16  t={time.time()-t0:.0f}s",
                  flush=True)
    lams = np.array(lams)
    print(f"wishart null lambda_max = {lams.mean():.2f} +- "
          f"{lams.std():.2f}  range [{lams.min():.2f}, {lams.max():.2f}]",
          flush=True)
    z = (lam_real - lams.mean()) / max(lams.std(), 1e-9)
    print(f"real vs wishart null: z = {z:.2f}", flush=True)
    print("(iid-entry null was 32.23 +- 0.50 -- category error check)",
          flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
