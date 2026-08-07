# -*- coding: utf-8 -*-
"""
Increment 151: matrix-level structure of the Gaussian part.

The chain consumes |Sum_{k != k'} xi_k xi_k' C_{k,k'}| for signs xi
from the EH_mu realization; the worst case over signs is governed by
the spectral norm of the C-matrix. If the (mask-normalized) C-matrix
behaves like a GOE random matrix, lambda_max ~ 2 sqrt(K v) and the
worst-sign bilinear form is ~ K lambda_max ~ 2 K^{3/2} sqrt(v) --
far below the K^2-scale trivial bound, with a sqrt(K) margin over
even the needed bound. Test:
  - build the full C matrix on a k-band (all pairs),
  - normalized entries Chat = C / sqrt(support) where support > 0,
  - compare lambda_max(Chat) with the GOE prediction 2 sqrt(K v),
    v = off-diagonal entry variance of Chat,
  - eigenvector delocalization: participation ratio of top vectors
    (GOE: PR ~ K * 2/3... report PR/K),
  - control: same matrix with entries sign-shuffled (exact GOE ref).
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
    rng = np.random.default_rng(20260812)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    K0, K1 = 3000, 3400  # 400 consecutive k
    ks = np.arange(K0, K1)
    K = len(ks)
    P0 = N // (2 * K1); P1 = 2 * P0
    ps = primes_in(P0, P1)
    pmax = (N - 2) // (K1 - 1)
    pp = ps[ps <= pmax]
    print(f"K = {K}, n_p = {len(pp)}", flush=True)

    # precompute mu(N - p k) rows
    rows = np.empty((K, len(pp)), dtype=np.int8)
    for i, k in enumerate(ks):
        rows[i] = mu[N - pp * int(k)]
        if i % 100 == 0:
            print(f"rows {i}/{K}  t={time.time()-t0:.0f}s", flush=True)
    R = rows.astype(np.float32)
    C = R @ R.T                      # C[i,j] = sum_p mu()mu()
    NZ = (np.abs(R) @ np.abs(R.T))   # support counts
    np.fill_diagonal(C, 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        Chat = np.where(NZ > 50, C / np.sqrt(np.maximum(NZ, 1)), 0.0)
    np.fill_diagonal(Chat, 0.0)
    Chat = (Chat + Chat.T) / 2

    off = Chat[np.triu_indices(K, 1)]
    live = off[off != 0]
    v = live.var()
    frac_live = len(live) / len(off)
    print(f"live entries {frac_live:.3f}, entry var v = {v:.3f}",
          flush=True)

    w, V = np.linalg.eigh(Chat)
    lam = np.abs(w).max()
    goe = 2 * np.sqrt(K * v * frac_live)  # variance diluted by mask
    print(f"lambda_max = {lam:.2f}  GOE pred = {goe:.2f}  "
          f"ratio = {lam/goe:.3f}", flush=True)

    # participation ratio of top-3 |eigenvalue| vectors
    idx = np.argsort(-np.abs(w))[:3]
    for j in idx:
        vec = V[:, j]
        pr = 1.0 / np.sum(vec**4) / K
        print(f"  eig {w[j]:8.2f}  PR/K = {pr:.3f} (GOE ~ 0.33, "
              f"localized ~ 1/K)", flush=True)

    # control: sign-shuffled matrix (same magnitudes, random signs)
    mag = np.abs(Chat[np.triu_indices(K, 1)])
    sgn = rng.choice([-1.0, 1.0], size=mag.shape)
    Ctrl = np.zeros_like(Chat)
    Ctrl[np.triu_indices(K, 1)] = mag * sgn
    Ctrl = Ctrl + Ctrl.T
    wc = np.linalg.eigvalsh(Ctrl)
    lamc = np.abs(wc).max()
    print(f"control lambda_max = {lamc:.2f}  real/control = "
          f"{lam/lamc:.3f}", flush=True)

    # the actual bilinear worst over random sign draws (sanity)
    best = 0.0
    for _ in range(2000):
        xi = rng.choice([-1.0, 1.0], size=K)
        best = max(best, abs(xi @ Chat @ xi))
    print(f"random-sign bilinear max (2000 draws) = {best:.1f}  "
          f"vs K*lambda_max = {K*lam:.1f}", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
