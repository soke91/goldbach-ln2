# -*- coding: utf-8 -*-
"""
Increment 164: anatomy of the elevated N = 799_999_996 (= 2^2 x odd).

Inc. 163: E1_ratio = 1.202 / 1.424 / 1.153 (mid band +3.7 sigma) with
depressed r1 -- heavy tail across k. Candidate mechanism: at v2(N)=2
the m-parity sub-sums of D(k) = Sum_m mu(m) mu(N-mk) are aliases of a
common underlying field, giving POSITIVE cross-correlation and
variance elevation (the integer-field mask was never blind-verified
at v2(N)=2; only the prime-indexed pair field was).

Tests:
  (1) per-class ratios: group k by v2(k) in {0, 1, >=2}; per class
      Sum D^2 / Sum sup, and the top-5 outlier k's contribution.
  (2) parity sub-sums: D = D_even + D_odd (m even / m odd);
      measure corr(D_even, D_odd) across k and the variance budget
      Var(D) vs Var(D_even) + Var(D_odd). Positive correlation
      confirms the aliasing mechanism.
  (3) control: same decomposition at the unit-consistent
      N = 699_999_994 (v2 = 1).
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

def vq(n, q, cap=2):
    v = 0
    while n % q == 0 and v < cap:
        n //= q; v += 1
    return v

def analyze(N, mu, rng, t0):
    SQ = int(N ** 0.5)
    K0, K1 = 3000, 4000
    ks = rng.choice(np.arange(K0, K1), size=250, replace=False)
    D = np.zeros(len(ks)); sup = np.zeros(len(ks))
    De = np.zeros(len(ks)); Do = np.zeros(len(ks))
    v2k = np.zeros(len(ks), dtype=int)
    for i, k in enumerate(ks):
        k = int(k)
        ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
        t = mu[ms].astype(np.int64) * mu[N - k * ms]
        D[i] = t.sum(); sup[i] = np.count_nonzero(t)
        even = (ms % 2 == 0)
        De[i] = t[even].sum(); Do[i] = t[~even].sum()
        v2k[i] = vq(k, 2)
        if i % 50 == 49:
            print(f"  k {i+1}/250  t={time.time()-t0:.0f}s", flush=True)

    print(f"  overall ratio = {D.dot(D)/sup.sum():.3f}", flush=True)
    for c in (0, 1, 2):
        m = v2k == c
        if m.sum() < 10:
            continue
        s = sup[m].sum()
        r = D[m].dot(D[m]) / max(s, 1)
        print(f"  v2(k)={c}: n={m.sum():3d}  ratio={r:.3f}", flush=True)
    # outliers
    contrib = D**2 / np.maximum(sup, 1)
    idx = np.argsort(-contrib)[:5]
    print("  top-5 k by D^2/sup: " +
          ", ".join(f"k={int(ks[j])}(v2={v2k[j]}):{contrib[j]:.1f}"
                    for j in idx), flush=True)
    # parity budget
    ve = De.var(); vo = Do.var(); vd = D.var()
    cc = np.corrcoef(De, Do)[0, 1]
    print(f"  parity split: Var(D)={vd:.0f} vs Var(De)+Var(Do)="
          f"{ve+vo:.0f}  corr(De,Do)={cc:+.3f}", flush=True)

def main():
    rng = np.random.default_rng(20260824)
    t0 = time.time()
    for N in (799_999_996, 699_999_994):
        print(f"=== N = {N} (v2={vq(N,2)}) ===", flush=True)
        mu = mobius_full(N + 1)
        print(f"  mu ready t={time.time()-t0:.0f}s", flush=True)
        analyze(N, mu, rng, t0)
        del mu
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
