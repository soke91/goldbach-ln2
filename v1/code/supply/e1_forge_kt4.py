# -*- coding: utf-8 -*-
"""
Forge kill-test K4 (increment 176): can the k-average replace the
N-average?

The N-averaged theorem is provable because summing over N linearizes
the pair constraint into shift-correlations depending only on the
shift h = (m - m')k -- the N-average supplies decorrelation across
(m, m') pairs. A descent to fixed N via the k-average requires the
dual field
    T(m, delta) = Sum_{k ~ K} mu(N - mk) mu(N - (m+delta)k)
to retain, within a single N, the shift-structure the N-average
exploits. Pre-registered tests (600 m x 4 delta, K-average over 2000
k):
  (a) coherence: |mean_m T_norm(m, delta)| vs 1/sqrt(M) random
      prediction -- z >= 4 at two deltas = alive coordinate;
  (b) shift-invariance: autocorrelation of m -> T_norm(m, delta) at
      lag 1 -- |rho| >= 0.15 at two deltas = alive coordinate;
  (c) baseline: T_norm unit-Gaussian (Conjecture L check).
No flags => DEAD: the N-average is load-bearing and the k-average
supplies no substitute decorrelation structure.
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

def main():
    rng = np.random.default_rng(20260901)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    ks = np.arange(2000, 4000, dtype=np.int64)  # k-average, 2000 terms
    M0, MM = 40000, 600                          # m in [40000, 40600)
    DELTAS = [1, 2, 6, 30]

    flags = 0
    for delta in DELTAS:
        Tn = np.full(MM, np.nan)
        for i in range(MM):
            m = M0 + i
            a = mu[N - m * ks].astype(np.int64)
            b = mu[N - (m + delta) * ks].astype(np.int64)
            t = a * b
            nz = int(np.count_nonzero(t))
            if nz > 100:
                Tn[i] = t.sum() / np.sqrt(nz)
        ok = ~np.isnan(Tn)
        v = Tn[ok]; n = len(v)
        mean = v.mean(); se = v.std() / np.sqrt(n)
        zc = mean / max(se, 1e-12)
        m2 = float(np.mean(v ** 2))
        # lag-1 autocorr on consecutive m (both viable)
        pairs = [(Tn[i], Tn[i+1]) for i in range(MM - 1)
                 if ok[i] and ok[i+1]]
        if len(pairs) > 50:
            aa = np.array([p[0] for p in pairs])
            bb = np.array([p[1] for p in pairs])
            rho = float(np.corrcoef(aa, bb)[0, 1])
        else:
            rho = float('nan')
        f = (abs(zc) >= 4) or (abs(rho) >= 0.15)
        flags += int(f)
        print(f"delta={delta:3d}  n={n:4d}  coher_z={zc:+.2f}  "
              f"m2={m2:.3f}  rho1={rho:+.3f}"
              f"{'  <-- FLAG' if f else ''}", flush=True)
        print(f"  t={time.time()-t0:.0f}s", flush=True)

    print(f"=== K4 KILL-TEST: {flags} flags ===", flush=True)
    print("verdict:",
          "ALIVE coordinate found" if flags >= 2 else
          ("MARGINAL (repeat)" if flags == 1 else
           "DEAD (the N-average is load-bearing; the k-average "
           "supplies no substitute decorrelation)"), flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
