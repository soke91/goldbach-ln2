# -*- coding: utf-8 -*-
"""
Increment 152: is the 15% spectral excess pure mask geometry?

Inc. 151: lambda_max(real) / lambda_max(sign-shuffled control) = 1.152.
Sign-shuffling preserves magnitudes but destroys any support-pattern
correlation. The factorization law's EXACT null is: same support
pattern (which pairs live, with what support), entries = fresh
Gaussians at half-normal scale. If lambda_max under that null matches
the real 36.76, the excess is mask geometry and the law closes at
matrix level. Otherwise the gap is genuine matrix-level structure.
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

def main():
    rng = np.random.default_rng(20260813)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    K0, K1 = 3000, 3400
    ks = np.arange(K0, K1); K = len(ks)
    P0 = N // (2 * K1); P1 = 2 * P0
    ps = primes_in(P0, P1)
    pp = ps[ps <= (N - 2) // (K1 - 1)]

    rows = np.empty((K, len(pp)), dtype=np.int8)
    for i, k in enumerate(ks):
        rows[i] = mu[N - pp * int(k)]
    R = rows.astype(np.float32)
    C = R @ R.T
    NZ = np.abs(R) @ np.abs(R.T)
    np.fill_diagonal(C, 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        Chat = np.where(NZ > 50, C / np.sqrt(np.maximum(NZ, 1)), 0.0)
    np.fill_diagonal(Chat, 0.0)
    Chat = (Chat + Chat.T) / 2
    live = (Chat != 0)
    w = np.linalg.eigvalsh(Chat)
    lam_real = np.abs(w).max()
    print(f"real lambda_max = {lam_real:.2f}", flush=True)

    # exact factorization-law null: same live pattern, fresh N(0,1)
    lams = []
    iu = np.triu_indices(K, 1)
    live_u = live[iu]
    for d in range(24):
        g = np.zeros(len(live_u), dtype=np.float64)
        g[live_u] = rng.standard_normal(live_u.sum())
        M = np.zeros((K, K))
        M[iu] = g
        M = M + M.T
        lams.append(np.abs(np.linalg.eigvalsh(M)).max())
        if d % 6 == 5:
            print(f"  null draws {d+1}/24  t={time.time()-t0:.0f}s",
                  flush=True)
    lams = np.array(lams)
    print(f"null lambda_max = {lams.mean():.2f} +- {lams.std():.2f}  "
          f"range [{lams.min():.2f}, {lams.max():.2f}]", flush=True)
    zsc = (lam_real - lams.mean()) / max(lams.std(), 1e-9)
    print(f"real vs null: z = {zsc:.2f}", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
