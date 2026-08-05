# -*- coding: utf-8 -*-
"""
Increment 168: the definitive E1 grid -- 4 N x 2 bands x 300 k with
bootstrap SE. Replaces the noisy 120-150-k readings of increments
163/166 (whose challenges 3 and 4 both settled as estimator noise).
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
    rng = np.random.default_rng(20260828)
    t0 = time.time()
    table = []
    for N in (999_999_998, 899_999_998, 799_999_996, 1_999_999_998):
        print(f"=== N = {N} ===", flush=True)
        mu = mobius_full(N + 1)
        print(f"  mu ready t={time.time()-t0:.0f}s", flush=True)
        SQ = int(N ** 0.5)
        for K0, K1 in ((3000, 4000), (10000, 20000)):
            ks = rng.choice(np.arange(K0, K1), size=300, replace=False)
            D = np.zeros(len(ks)); sup = np.zeros(len(ks))
            for i, k in enumerate(ks):
                k = int(k)
                ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
                t = mu[ms].astype(np.int64) * mu[N - k * ms]
                D[i] = t.sum(); sup[i] = np.count_nonzero(t)
            ratio = D.dot(D) / sup.sum()
            boots = []
            idx0 = np.arange(len(ks))
            for _ in range(2000):
                idx = rng.choice(idx0, size=len(ks), replace=True)
                boots.append(D[idx].dot(D[idx]) / sup[idx].sum())
            bs = float(np.std(boots))
            z = (ratio - 1) / bs
            line = (f"N={N}  K[{K0},{K1}): ratio={ratio:.3f} "
                    f"+- {bs:.3f}  z={z:+.2f}")
            table.append(line)
            print("  " + line, flush=True)
        del mu
    print("=== FINAL GRID ===", flush=True)
    for line in table:
        print("  " + line, flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
