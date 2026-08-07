# -*- coding: utf-8 -*-
"""
Increment 167: anatomy of the low band at v3(N) = 4.

Inc. 166: N = 1_999_999_998 = 2 x 3^4 x 37 x 333667 reads
E1_ratio 0.729 / r1 0.581 in K [1e4, 2e4] -- downward. Cell mixing
cannot push a mean of unit cells below unit, so if this is not noise
it is a genuine negative cross-correlation between m-residue-class
sub-sums (multi-level 3-adic sign locking) -- the first candidate
structure OUTSIDE the mask. Control: N = 999_999_998 (v3 = 0) read
0.922 in the same band.

Measure at the low band (200 k):
  (1) per-v3(k) class ratios;
  (2) for each k, split D(k) into 9 sub-sums by m mod 9; variance
      budget Var(D) vs Sum_c Var-contribution, plus the mean pairwise
      correlation of sub-sums (computed per k over the class vector,
      then averaged);
  (3) bootstrap SE of the band ratio to size the deviation.
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

def v3(n, cap=5):
    v = 0
    while n % 3 == 0 and v < cap:
        n //= 3; v += 1
    return v

def main():
    rng = np.random.default_rng(20260827)
    N = 1_999_999_998
    t0 = time.time()
    print("building mu to 2e9 (2GB)...", flush=True)
    mu = mobius_full(N + 1)
    print(f"mu ready t={time.time()-t0:.0f}s", flush=True)
    SQ = int(N ** 0.5)
    K0, K1 = 10000, 20000
    ks = rng.choice(np.arange(K0, K1), size=200, replace=False)

    D = np.zeros(len(ks)); sup = np.zeros(len(ks))
    subs = np.zeros((len(ks), 9)); subsup = np.zeros((len(ks), 9))
    v3k = np.zeros(len(ks), dtype=int)
    for i, k in enumerate(ks):
        k = int(k)
        ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
        t = mu[ms].astype(np.int64) * mu[N - k * ms]
        D[i] = t.sum(); sup[i] = np.count_nonzero(t)
        c = (ms % 9).astype(int)
        for cc in range(9):
            m = c == cc
            subs[i, cc] = t[m].sum()
            subsup[i, cc] = np.count_nonzero(t[m])
        v3k[i] = v3(k)
        if i % 40 == 39:
            print(f"  k {i+1}/200  t={time.time()-t0:.0f}s", flush=True)

    ratio = D.dot(D) / sup.sum()
    boots = []
    idx0 = np.arange(len(ks))
    for _ in range(2000):
        idx = rng.choice(idx0, size=len(ks), replace=True)
        boots.append(D[idx].dot(D[idx]) / sup[idx].sum())
    bs = np.std(boots)
    print(f"band ratio = {ratio:.3f}  bootstrap SE = {bs:.3f}  "
          f"z vs 1.0 = {(ratio-1)/bs:+.2f}", flush=True)

    print("--- per v3(k) ---", flush=True)
    for c in range(4):
        m = v3k == c
        if m.sum() < 8:
            continue
        s = sup[m].sum()
        if s < 1:
            print(f"  v3(k)={c}: n={m.sum():3d}  annihilated", flush=True)
            continue
        r = D[m].dot(D[m]) / s
        r1m = np.mean(np.abs(D[m]) / np.sqrt(np.maximum(sup[m], 1)))
        print(f"  v3(k)={c}: n={m.sum():3d}  ratio={r:.3f}  "
              f"r1={r1m:.3f}  supfrac_k={np.mean(sup[m]):.0f}", flush=True)

    print("--- m mod 9 sub-sum budget ---", flush=True)
    # per-k: D = sum_c subs; unit prediction per class:
    # E[subs_c^2] = subsup_c. cross-term per k:
    cross = (D**2 - (subs**2).sum(axis=1))
    pred_diag = subsup.sum(axis=1)
    diag = (subs**2).sum(axis=1)
    print(f"  mean diag/pred = {np.mean(diag / np.maximum(pred_diag,1)):.3f}"
          f"  (unit if classes individually Gaussian)", flush=True)
    print(f"  mean cross/diag = {np.mean(cross / np.maximum(diag,1)):+.3f}"
          f"  (0 if classes independent; negative = anti-corr)",
          flush=True)
    # which classes are live?
    print("  class support fractions (mean over k): " +
          " ".join(f"{c}:{np.mean(subsup[:,c]/np.maximum(sup,1)):.2f}"
                   for c in range(9)), flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
