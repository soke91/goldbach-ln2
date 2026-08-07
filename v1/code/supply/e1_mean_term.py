# -*- coding: utf-8 -*-
"""
Increment 153: is the z=9 spectral excess a MEAN (main term)?

Inc. 152: real lambda_max exceeds the exact factorization null by
z = 9. Pairwise stats (|C|, C^2, kurtosis) were Gaussian-exact, but
E[C] was never measured -- and the classical dispersion computation
subtracts a main term precisely because such sums have a
singular-series-type nonzero mean. A tiny mean m ~ 0.03 on live
entries adds a near-rank-one component lifting lambda_max by
~ K * live_frac * m without touching |C|-statistics.

Test: (1) measure the mean of live normalized entries (SE ~ 0.005);
(2) recompute lambda_max after global centering and after
sign-class centering (entries grouped by the (-1)^{...} local sign
pattern of the pair); (3) compare with the null 32.23 +- 0.50.
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
    rng = np.random.default_rng(20260814)
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

    iu = np.triu_indices(K, 1)
    lu = live[iu]
    e = Chat[iu][lu]
    n = len(e)
    m = e.mean(); se = e.std() / np.sqrt(n)
    print(f"live entries n = {n}", flush=True)
    print(f"(1) mean of live entries = {m:.4f} +- {se:.4f}  "
          f"(z vs 0: {m/se:.2f})", flush=True)

    w0 = np.linalg.eigvalsh(Chat)
    print(f"raw lambda_max = {np.abs(w0).max():.2f}  "
          f"(null 32.23 +- 0.50)", flush=True)

    # (2a) global centering on live entries
    Cg = Chat.copy()
    Cg[live] -= m
    np.fill_diagonal(Cg, 0.0)
    wg = np.linalg.eigvalsh((Cg + Cg.T) / 2)
    print(f"(2a) globally centered lambda_max = {np.abs(wg).max():.2f}",
          flush=True)

    # (2b) class centering by k-parity-ish local signature:
    # class of pair = (k mod 2, k' mod 2, k mod 3, k' mod 3) unordered
    sig2 = ks % 2; sig3 = ks % 3
    cls = {}
    I, J = np.nonzero(np.triu(live, 1))
    keys = [tuple(sorted(((sig2[i], sig3[i]), (sig2[j], sig3[j]))))
            for i, j in zip(I, J)]
    vals = Chat[I, J]
    Cc = Chat.copy()
    from collections import defaultdict
    acc = defaultdict(list)
    for t, v in zip(keys, vals):
        acc[t].append(v)
    means = {t: np.mean(v) for t, v in acc.items()}
    print("(2b) class means:", flush=True)
    for t, v in sorted(acc.items(), key=lambda x: -len(x[1])):
        mm = np.mean(v); s = np.std(v)/np.sqrt(len(v))
        print(f"    {str(t):40s} n={len(v):6d}  mean={mm:+.4f}"
              f" (z={mm/max(s,1e-9):+.1f})", flush=True)
    adj = np.array([means[t] for t in keys])
    Cc[I, J] = vals - adj
    Cc[J, I] = Cc[I, J]
    np.fill_diagonal(Cc, 0.0)
    wc = np.linalg.eigvalsh(Cc)
    print(f"(2b) class-centered lambda_max = {np.abs(wc).max():.2f}",
          flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
