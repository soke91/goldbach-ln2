# -*- coding: utf-8 -*-
"""
Increment 169: settle the suspect cell -- N = 899_999_998, K [3000,4000).

Two independent samples (150 k and 300 k) both read 0.83. A ratio
below unit cannot come from mixing unit cells across k; it requires
negative correlation WITHIN each D(k) sum. Settle with 800 k
(non-overlapping seed) + bootstrap; if z < -3, structure candidate 5
(the first on-support over-cancellation); else multiple-comparison
tail, grid closes.
"""
import numpy as np, time

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i*i::i] = False
    return np.nonzero(s)[0]

def mobius_full(B, seg=20_000_000):
    bp = primes_upto(int(B ** 0.5) + 1)
    out = np.empty(B, dtype=np.int8)
    for lo in range(0, B, seg):
        hi = min(lo + seg, B)
        mu = np.ones(hi - lo, dtype=np.int8)
        val = np.arange(lo, hi, dtype=np.int64)
        if lo == 0:
            val[0] = 1
        for p in bp:
            p = int(p)
            start = ((lo + p - 1) // p) * p - lo
            mu[start::p] *= -1
            val[start::p] //= p
            pp = p * p
            if pp < hi:
                start2 = ((lo + pp - 1) // pp) * pp - lo
                mu[start2::pp] = 0
        mu[val > 1] *= -1
        out[lo:hi] = mu
    out[0] = 0
    if B > 1:
        out[1] = 1
    return out

def main():
    rng = np.random.default_rng(20260829)
    N = 899_999_998
    t0 = time.time()
    mu = mobius_full(N + 1)
    print(f"mu ready t={time.time()-t0:.0f}s", flush=True)
    SQ = int(N ** 0.5)
    K0, K1 = 3000, 4000
    ks = rng.permutation(np.arange(K0, K1))[:800]
    D = np.zeros(len(ks)); sup = np.zeros(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
        t = mu[ms].astype(np.int64) * mu[N - k * ms]
        D[i] = t.sum(); sup[i] = np.count_nonzero(t)
        if i % 100 == 99:
            print(f"k {i+1}/800  t={time.time()-t0:.0f}s", flush=True)
    ratio = D.dot(D) / sup.sum()
    boots = []
    idx0 = np.arange(len(ks))
    for _ in range(3000):
        idx = rng.choice(idx0, size=len(ks), replace=True)
        boots.append(D[idx].dot(D[idx]) / sup[idx].sum())
    bs = float(np.std(boots))
    z = (ratio - 1) / bs
    print(f"SETTLE: ratio = {ratio:.3f} +- {bs:.3f}  z = {z:+.2f}",
          flush=True)
    # per-k r1 distribution summary for extra diagnostics
    r1 = np.abs(D) / np.sqrt(np.maximum(sup, 1))
    print(f"mean r1 = {r1.mean():.3f}  (0.798)  frac r1>2: "
          f"{np.mean(r1 > 2):.3f}", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
