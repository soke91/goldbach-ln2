# -*- coding: utf-8 -*-
"""
Increment 165: settle the elevated-N question with power.

Inc. 164 refuted the parity-aliasing mechanism (corr +0.061) and the
elevation regressed on resample (1.424 -> 1.113). Remaining doubt:
three bands moved together at inc. 163. Settle: 600 k in the mid band
at the elevated N and at the control N -- ratio SE ~ 0.06 (empirical,
bootstrap) decides between "heavy-tail noise" and "real elevation".
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
    rng = np.random.default_rng(20260825)
    t0 = time.time()
    for N in (799_999_996, 699_999_994):
        print(f"=== N = {N} ===", flush=True)
        mu = mobius_full(N + 1)
        print(f"  mu ready t={time.time()-t0:.0f}s", flush=True)
        SQ = int(N ** 0.5)
        K0, K1 = 3000, 4000
        ks = rng.choice(np.arange(K0, K1), size=600, replace=False)
        D = np.zeros(len(ks)); sup = np.zeros(len(ks))
        for i, k in enumerate(ks):
            k = int(k)
            ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
            t = mu[ms].astype(np.int64) * mu[N - k * ms]
            D[i] = t.sum(); sup[i] = np.count_nonzero(t)
            if i % 100 == 99:
                print(f"  k {i+1}/600  t={time.time()-t0:.0f}s",
                      flush=True)
        ratio = D.dot(D) / sup.sum()
        # bootstrap SE
        boots = []
        idx0 = np.arange(len(ks))
        for _ in range(2000):
            idx = rng.choice(idx0, size=len(ks), replace=True)
            boots.append(D[idx].dot(D[idx]) / sup[idx].sum())
        boots = np.array(boots)
        print(f"  RATIO = {ratio:.3f}  bootstrap SE = {boots.std():.3f}"
              f"  z vs 1.0 = {(ratio-1)/boots.std():+.2f}", flush=True)
        del mu
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
